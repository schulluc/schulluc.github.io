import { add, sub, mul, dot, cross, norm, dist } from "./vec.js";

export const SHAPE_DEFAULT = {
  hW: 1.40,      // α-helix curled ribbon half-width (Å)
  sW: 1.26,      // β-strand shaft half-width (Å)
  sT: 0.24,      // β-strand 3D architectural slab thickness (Å)
  head: 2.15,    // β-strand arrowhead half-width (Å)
  coilR: 0.32,   // loop/coil rope half-width (Å)
  flat: 0.78,    // β-strand pleat smoothing (0..1)
  arrows: true,
};

function computeHelixAxis(CA, s, e) {
  const len = e - s + 1;
  const ax = new Array(len).fill(null);
  for (let i = s + 2; i <= e - 2; i++) {
    ax[i - s] = mul(
      add(
        add(add(mul(CA[i - 2], 0.5), CA[i - 1]), add(CA[i], CA[i + 1])),
        mul(CA[i + 2], 0.5)
      ),
      0.25
    );
  }
  const def = ax.map((v, k) => (v ? k : -1)).filter(k => k >= 0);
  if (def.length < 2) {
    const c = mul(
      CA.slice(s, e + 1).reduce((a, p) => add(a, p), [0, 0, 0]),
      1 / len
    );
    const d = norm(sub(CA[e], CA[s]));
    return CA.slice(s, e + 1).map(p => add(c, mul(d, dot(sub(p, c), d))));
  }
  const f = def[0];
  const l = def[def.length - 1];
  const stepF = sub(ax[f], ax[f + 1]);
  const stepL = sub(ax[l], ax[l - 1]);
  for (let k = f - 1; k >= 0; k--) ax[k] = add(ax[k + 1], stepF);
  for (let k = l + 1; k < len; k++) ax[k] = add(ax[k - 1], stepL);
  // Smooth helix axis so cylindrical ribbon normals are silky and kink-free
  for (let pass = 0; pass < 3; pass++) {
    const sm = ax.map((v, k) =>
      k === 0 || k === len - 1 ? v : mul(add(add(ax[k - 1], mul(v, 2)), ax[k + 1]), 0.25)
    );
    for (let k = 0; k < len; k++) ax[k] = sm[k];
  }
  return ax;
}

export function buildGeometry(model, shapeOverrides = {}) {
  const SH = { ...SHAPE_DEFAULT, ...shapeOverrides };
  const totalRes = model.chains.reduce((acc, ch) => acc + ch.ss.length, 0);
  if (totalRes < 5) throw new Error("Too few residues to build ribbon geometry.");

  const SUB = totalRes > 2200 ? 3 : Math.max(4, Math.min(9, Math.round(2400 / totalRes)));
  const allPts = [];
  let hCount = 0;
  let eCount = 0;

  for (let ci = 0; ci < model.chains.length; ci++) {
    const ch = model.chains[ci];
    const SS = ch.ss;
    const N = SS.length;
    if (N < 5) continue;

    for (let i = 0; i < N; i++) {
      if (SS[i] === "H") hCount++;
      else if (SS[i] === "E") eCount++;
    }

    const CA0 = new Array(N);
    const O0 = new Array(N);
    for (let i = 0; i < N; i++) {
      CA0[i] = [ch.ca[i * 3], ch.ca[i * 3 + 1], ch.ca[i * 3 + 2]];
      O0[i] = [ch.o[i * 3], ch.o[i * 3 + 1], ch.o[i * 3 + 2]];
    }

    // Smooth pleated Cα zig-zag inside β-strands for clean Richardson flat arrows
    let CA = CA0.map(p => p.slice());
    for (let pass = 0; pass < 4 && SH.flat > 0; pass++) {
      const next = CA.map((p, i) => {
        if (SS[i] !== "E" || i === 0 || i === N - 1) return p;
        const avg = mul(add(add(CA[i - 1], mul(p, 2)), CA[i + 1]), 0.25);
        return add(mul(p, 1 - SH.flat), mul(avg, SH.flat));
      });
      CA = next;
    }

    // Compute smooth α-helix axes for cylindrical twist-free curled ribbons
    const AX = new Array(N).fill(null);
    for (let i = 0; i < N; ) {
      if (SS[i] !== "H") {
        i++;
        continue;
      }
      let j = i;
      while (j + 1 < N && SS[j + 1] === "H") j++;
      const axRun = computeHelixAxis(CA0, i, j);
      for (let k = 0; k < axRun.length; k++) AX[i + k] = axRun[k];
      i = j + 1;
    }

    // Build twist-free per-residue SIDE vectors:
    // - On α-helices: purely cylindrical (cross(radialOut, tangent)) so helices NEVER Möbius-flip!
    // - On β-strands & coils: peptide C=O projected perpendicular to local backbone tangent, sign-aligned.
    const SIDE = new Array(N);
    for (let i = 0; i < N; i++) {
      const tLocal = norm(sub(CA[Math.min(N - 1, i + 1)], CA[Math.max(0, i - 1)]));
      let sv;
      if (SS[i] === "H" && AX[i]) {
        const radOut = norm(sub(CA0[i], AX[i]));
        sv = norm(cross(tLocal, radOut));
      } else {
        const co = sub(O0[i], CA0[i]);
        const perp = sub(co, mul(tLocal, dot(co, tLocal)));
        sv = norm(perp);
      }
      if (i > 0 && SS[i] !== "H" && dot(sv, SIDE[i - 1]) < 0) {
        sv = mul(sv, -1);
      }
      SIDE[i] = sv;
    }

    // Smooth SIDE strongly inside β-strands so sheets are flat, sleek, and untwisted
    for (let pass = 0; pass < 7; pass++) {
      const sm = SIDE.map((v, i) => {
        if (i === 0 || i === N - 1 || SS[i] !== "E") return v;
        const prev = SS[i - 1] === "E" ? SIDE[i - 1] : v;
        const next = SS[i + 1] === "E" ? SIDE[i + 1] : v;
        return norm(add(add(mul(prev, 0.25), mul(v, 0.5)), mul(next, 0.25)));
      });
      for (let i = 0; i < N; i++) SIDE[i] = sm[i];
    }

    const pts = [];
    for (let i = 0; i < N - 1; i++) {
      const p0 = CA[Math.max(0, i - 1)];
      const p1 = CA[i];
      const p2 = CA[i + 1];
      const p3 = CA[Math.min(N - 1, i + 2)];
      const gap = dist(CA0[i], CA0[i + 1]) > 4.4;

      // Ensure SIDE interpolation across H<->C boundary doesn't flip 180°
      let sA = SIDE[i];
      let sB = SIDE[i + 1];
      if (dot(sA, sB) < 0) sB = mul(sB, -1);

      for (let s = 0; s < SUB; s++) {
        const t = s / SUB;
        const t2 = t * t;
        const t3 = t2 * t;
        const c = k =>
          0.5 *
          (2 * p1[k] +
            (-p0[k] + p2[k]) * t +
            (2 * p0[k] - 5 * p1[k] + 4 * p2[k] - p3[k]) * t2 +
            (-p0[k] + 3 * p1[k] - 3 * p2[k] + p3[k]) * t3);
        const pos = [c(0), c(1), c(2)];
        const axPt =
          AX[i] && AX[i + 1]
            ? add(mul(AX[i], 1 - t), mul(AX[i + 1], t))
            : AX[i] || AX[i + 1] || null;
        pts.push({
          p: pos,
          ax: axPt,
          side: norm(add(mul(sA, 1 - t), mul(sB, t))),
          res: i,
          subIdx: s,
          t,
          gap,
          chain: ci,
        });
      }
    }
    pts.push({
      p: CA[N - 1],
      ax: AX[N - 1],
      side: SIDE[N - 1],
      res: N - 1,
      subIdx: 0,
      t: 0,
      gap: true,
      chain: ci,
    });

    for (let k = 0; k < pts.length; k++) {
      const q = pts[k];
      const prevP = k > 0 && !pts[k - 1].gap ? pts[k - 1].p : q.p;
      const nextP = k < pts.length - 1 && !q.gap ? pts[k + 1].p : q.p;
      const T = norm(sub(nextP, prevP));

      const i = q.res;
      const t = q.t;
      const prevSS = i > 0 ? SS[i - 1] : "C";
      const ownSS = SS[i];
      const nextSS = i < N - 1 ? SS[i + 1] : "C";

      let S, Nrm;
      if (ownSS === "H" && q.ax) {
        // Pure cylindrical spiral frame: Nrm points radially OUTWARD from helix axis!
        const rad = sub(q.p, q.ax);
        const radPerp = sub(rad, mul(T, dot(rad, T)));
        Nrm = norm(radPerp);
        S = norm(cross(Nrm, T));
        q.helixSign = 1; // Nrm is guaranteed radially outward (+1)
      } else {
        S = norm(sub(q.side, mul(T, dot(q.side, T))));
        Nrm = norm(cross(T, S));
        q.helixSign = 1;
      }
      q.tan = T;
      q.side = S;
      q.nrm = Nrm;

      let kind = "C";
      let wIn = 0;
      let wOut = 0;
      let th = 0;

      if (ownSS === "H") {
        kind = "H";
        const enterW =
          prevSS !== "H" ? SH.coilR + (SH.hW - SH.coilR) * Math.min(1, t / 0.85) : SH.hW;
        const exitW =
          nextSS !== "H" ? SH.hW - (SH.hW - SH.coilR) * Math.max(0, (t - 0.1) / 0.9) : SH.hW;
        wIn = wOut = prevSS !== "H" ? enterW : exitW;
        th = 0.05;
      } else if (ownSS === "E") {
        kind = "E";
        th = SH.sT;
        if (nextSS === "E") {
          wIn = wOut = SH.sW;
        } else if (SH.arrows) {
          const wCurr = SH.head * (1 - t) + SH.coilR * t;
          wIn = wCurr;
          wOut = wCurr;
          if (q.subIdx === 0) {
            wIn = SH.sW;
            wOut = SH.head;
            q.arrowShoulder = true;
            q.shaftW = SH.sW;
          }
        } else {
          wIn = wOut = SH.sW * (1 - 0.5 * t) + SH.coilR * 0.5 * t;
        }
        if (prevSS !== "E" && q.subIdx === 0) {
          q.strandStart = true;
        }
      }

      q.kind = kind;
      q.wIn = wIn;
      q.wOut = wOut;
      q.th = th;
    }

    for (const p of pts) allPts.push(p);
  }

  const radii = allPts.map(q => Math.hypot(q.p[0], q.p[1], q.p[2])).sort((a, b) => a - b);
  const rMax = radii[radii.length - 1] || 20;
  const r96 = radii[Math.floor(radii.length * 0.96)] || rMax;
  const R = Math.min(rMax, r96 * 1.12) + 2.6;

  return {
    pts: allPts,
    R,
    SUB,
    shape: SH,
    residues: totalRes,
    chains: model.chains.length,
    ssCounts: { H: hCount, E: eCount },
  };
}
