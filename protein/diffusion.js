import { SHAPE_DEFAULT } from "./geometry.js";

const smoothstep = t => {
  const x = Math.max(0, Math.min(1, t));
  return x * x * (3 - 2 * x);
};

const posMod = (x, m) => ((x % m) + m) % m;

/**
 * Deterministic confined 3D Brownian polypeptide random walk (Cα–Cα virtual bond = 3.72 Å < 4.4 Å)
 * so the SE(3) noise state remains a continuous disordered chain in both pen-and-ink and ASCII.
 */
function buildPolymerCloud(maxN) {
  const pts = new Array(maxN);
  let x = 0, y = 0, z = 0;
  let dx = 1, dy = 0, dz = 0;
  const bond = 3.72;

  const hash = s => {
    const h = Math.sin(s * 12.9898 + 78.233) * 43758.5453;
    return h - Math.floor(h);
  };

  for (let i = 0; i < maxN; i++) {
    pts[i] = [x, y, z];
    const u1 = hash(i * 3.1 + 1.7);
    const u2 = hash(i * 5.7 + 4.3);
    const theta = u1 * Math.PI * 2;
    const rz = u2 * 2 - 1;
    const rxy = Math.sqrt(Math.max(0, 1 - rz * rz));
    let rx = rxy * Math.cos(theta);
    let ry = rz;
    let rzz = rxy * Math.sin(theta);

    // Soft harmonic confinement toward origin so radius stays ~18–21 Å
    const distOrigin = Math.hypot(x, y, z) || 1;
    const pull = Math.max(0, (distOrigin - 11.5) / 10.0) * 0.68;
    rx -= (x / distOrigin) * pull;
    ry -= (y / distOrigin) * pull;
    rzz -= (z / distOrigin) * pull;

    // Slight chain persistence so coils look organic rather than sawtooth
    dx = dx * 0.34 + rx * 0.66;
    dy = dy * 0.34 + ry * 0.66;
    dz = dz * 0.34 + rzz * 0.66;
    const dl = Math.hypot(dx, dy, dz) || 1;
    dx /= dl;
    dy /= dl;
    dz /= dl;

    x += dx * bond;
    y += dy * bond;
    z += dz * bond;
  }

  let cx = 0, cy = 0, cz = 0;
  for (let i = 0; i < maxN; i++) {
    cx += pts[i][0];
    cy += pts[i][1];
    cz += pts[i][2];
  }
  cx /= maxN;
  cy /= maxN;
  cz /= maxN;
  for (let i = 0; i < maxN; i++) {
    pts[i][0] -= cx;
    pts[i][1] -= cy;
    pts[i][2] -= cz;
  }
  return pts;
}

function extractContacts(model) {
  const flat = [];
  for (const ch of model.chains) {
    for (let i = 0; i < ch.ss.length; i++) {
      flat.push([ch.ca[3 * i], ch.ca[3 * i + 1], ch.ca[3 * i + 2]]);
    }
  }
  const contacts = [];
  for (let i = 0; i < flat.length; i += 2) {
    for (let j = i + 7; j < flat.length; j += 3) {
      const d = Math.hypot(
        flat[i][0] - flat[j][0],
        flat[i][1] - flat[j][1],
        flat[i][2] - flat[j][2]
      );
      if (d < 7.8) contacts.push([i, j]);
    }
  }
  return { count: flat.length, contacts: contacts.slice(0, 85) };
}

/**
 * Attach scroll-driven de novo protein diffusion & refolding trajectory:
 * Stop 0: 2PAB (prealbumin β-sandwich dimer) — 0% ASCII on mobile at rest
 * Stop 1: SE(3) disordered polypeptide diffusion prior + pairwise restraint graph (sweeps up to 50% ASCII)
 * Stop 2: 1TIM (de novo-style (β/α)8 TIM barrel enzyme)
 * Stop 3: SE(3) disordered polypeptide diffusion prior + pairwise restraint graph
 */
export function attachDiffusionScroll(view, { baseModel, altModel, initialStop = 0 }) {
  const NUM_STOPS = 4;
  let dialTarget = initialStop;
  let dialSmooth = initialStop;
  let lastInputAt = -1e9;
  let isTouching = false;

  // Measure natural bounding radii of both folded structures for seamless zero-jump scaling
  view.setModelFast(altModel, null, null);
  const radiusB = view.getInfo()?.R || 26.5;
  view.setModelFast(baseModel, null, null);
  const radiusA = view.getInfo()?.R || 26.0;
  const radiusCloud = 26.4;

  const infoA = extractContacts(baseModel);
  const infoB = extractContacts(altModel);
  const cloud = buildPolymerCloud(Math.max(infoA.count, infoB.count));

  let restState = "2pab";
  let latestFrame = { cloudAmt: 0, useB: false, currentPts: null };

  function updateGeometry(u) {
    if (u < 0.008 || u > 3.992) {
      if (restState !== "2pab") {
        restState = "2pab";
        view.setModelFast(baseModel, null, radiusA);
      }
      return { cloudAmt: 0, useB: false, currentPts: null };
    }
    if (Math.abs(u - 2.0) < 0.008) {
      if (restState !== "1tim") {
        restState = "1tim";
        view.setModelFast(altModel, null, radiusB);
      }
      return { cloudAmt: 0, useB: true, currentPts: null };
    }

    restState = "morph";
    let useB = false;
    let cloudAmt = 0;
    if (u < 1.0) {
      useB = false;
      cloudAmt = smoothstep(u);
    } else if (u < 2.0) {
      useB = true;
      cloudAmt = smoothstep(2.0 - u);
    } else if (u < 3.0) {
      useB = true;
      cloudAmt = smoothstep(u - 2.0);
    } else {
      useB = false;
      cloudAmt = smoothstep(4.0 - u);
    }

    const src = useB ? altModel : baseModel;
    const baseR = useB ? radiusB : radiusA;
    const interpR = baseR + (radiusCloud - baseR) * cloudAmt;

    const keepSS = cloudAmt < 0.56;
    const foldFactor = Math.max(0, 1 - cloudAmt * 1.35);
    const shapeOverride = {
      hW: SHAPE_DEFAULT.coilR * 0.72 + (SHAPE_DEFAULT.hW - SHAPE_DEFAULT.coilR * 0.72) * foldFactor,
      sW: SHAPE_DEFAULT.coilR * 0.72 + (SHAPE_DEFAULT.sW - SHAPE_DEFAULT.coilR * 0.72) * foldFactor,
      head: SHAPE_DEFAULT.coilR * 0.80 + (SHAPE_DEFAULT.head - SHAPE_DEFAULT.coilR * 0.80) * foldFactor,
      sT: SHAPE_DEFAULT.sT * foldFactor,
      coilR: SHAPE_DEFAULT.coilR * (1 - 0.16 * cloudAmt),
    };

    const currentPts = [];
    const outChains = [];
    let globalIdx = 0;

    for (let ci = 0; ci < src.chains.length; ci++) {
      const ch = src.chains[ci];
      const n = ch.ss.length;
      const caArr = new Array(n * 3);
      const oArr = new Array(n * 3);
      const sepSign = !useB ? (ci === 0 ? -1 : 1) : 0;
      const sepShift = sepSign * Math.sin(cloudAmt * Math.PI) * 3.8;

      for (let i = 0; i < n; i++) {
        const cp = cloud[globalIdx];
        const fx = ch.ca[3 * i] + sepShift * 0.5;
        const fy = ch.ca[3 * i + 1];
        const fz = ch.ca[3 * i + 2] - sepShift * 0.7;

        const x = fx + (cp[0] - fx) * cloudAmt;
        const y = fy + (cp[1] - fy) * cloudAmt;
        const z = fz + (cp[2] - fz) * cloudAmt;

        const odx = ch.o[3 * i] - ch.ca[3 * i];
        const ody = ch.o[3 * i + 1] - ch.ca[3 * i + 1];
        const odz = ch.o[3 * i + 2] - ch.ca[3 * i + 2];

        caArr[3 * i] = x;
        caArr[3 * i + 1] = y;
        caArr[3 * i + 2] = z;
        oArr[3 * i] = x + odx;
        oArr[3 * i + 1] = y + ody;
        oArr[3 * i + 2] = z + odz;

        currentPts.push([x, y, z]);
        globalIdx++;
      }

      outChains.push({
        id: ch.id,
        ss: keepSS ? ch.ss : "C".repeat(n),
        ca: caArr,
        o: oArr,
      });
    }

    view.setModelFast(
      { id: src.id, title: src.title, residues: globalIdx, chains: outChains },
      shapeOverride,
      interpR
    );

    return { cloudAmt, useB, currentPts };
  }

  function onTick(dt, now) {
    const idleMs = now - lastInputAt;
    if (!isTouching && idleMs > 320) {
      const snapped = Math.round(dialTarget);
      dialTarget += (snapped - dialTarget) * (1 - Math.exp(-dt * 6.5));
    }

    dialSmooth += (dialTarget - dialSmooth) * (1 - Math.exp(-dt * 11.0));
    const u = posMod(dialSmooth, NUM_STOPS);
    latestFrame = updateGeometry(u);

    // Smoothly sweep ASCII into the protein (0% -> up to scrollMax=50%) upon scrolling
    const asciiProgress = smoothstep(Math.min(1, Math.max(0, Math.abs(dialSmooth))));
    view.setScrollFast(asciiProgress);

    return (
      isTouching ||
      Math.abs(dialTarget - dialSmooth) > 5e-4 ||
      (idleMs <= 320 && Math.abs(dialTarget - Math.round(dialTarget)) > 1e-3)
    );
  }

  function onBeforeMask({ rctx, dpr, R, rot, pixX, pixY, inkPalette }) {
    const { cloudAmt, useB, currentPts } = latestFrame;
    if (!currentPts || cloudAmt < 0.03) return;

    const vis = smoothstep((cloudAmt - 0.05) / 0.65);
    const contacts = useB ? infoB.contacts : infoA.contacts;
    const proj = currentPts.map(p => {
      const rp = rot(p);
      return { x: pixX(rp), y: pixY(rp), z: rp[2] };
    });

    // Dashed pairwise SE(3) distance-restraint lines
    rctx.save();
    rctx.setLineDash([2.2 * dpr, 3.4 * dpr]);
    rctx.lineWidth = 0.85 * dpr;
    for (const [i, j] of contacts) {
      if (i >= proj.length || j >= proj.length) continue;
      const a = proj[i];
      const b = proj[j];
      const cue = Math.max(0, Math.min(1, ((a.z + b.z) * 0.5 + R) / (2 * R)));
      rctx.strokeStyle = inkPalette[Math.round((0.18 + 0.45 * cue) * 31)];
      rctx.globalAlpha = vis * (0.32 + 0.45 * cue);
      rctx.beginPath();
      rctx.moveTo(a.x, a.y);
      rctx.lineTo(b.x, b.y);
      rctx.stroke();
    }
    rctx.setLineDash([]);

    // Crisp 3D Cα coordinate crosshairs in the diffusion cloud
    const crossLen = 2.6 * dpr;
    rctx.lineWidth = 1.05 * dpr;
    for (let i = 0; i < proj.length; i += 2) {
      const pt = proj[i];
      const cue = Math.max(0, Math.min(1, (pt.z + R) / (2 * R)));
      rctx.strokeStyle = inkPalette[Math.round((0.35 + 0.65 * cue) * 31)];
      rctx.globalAlpha = vis * (0.45 + 0.55 * cue);
      rctx.beginPath();
      rctx.moveTo(pt.x - crossLen, pt.y);
      rctx.lineTo(pt.x + crossLen, pt.y);
      rctx.moveTo(pt.x, pt.y - crossLen);
      rctx.lineTo(pt.x, pt.y + crossLen);
      rctx.stroke();
    }
    rctx.restore();
  }

  view.set({ onTick, onBeforeMask });
  onTick(0.016, performance.now());
  if (initialStop !== 0) {
    view.render();
  }

  function onWheel(e) {
    e.preventDefault();
    const rawDy = e.deltaMode === 1 ? e.deltaY * 28 : e.deltaMode === 2 ? e.deltaY * 320 : e.deltaY;
    const rawDx = e.deltaMode === 1 ? e.deltaX * 28 : e.deltaX;
    const delta = Math.abs(rawDy) >= Math.abs(rawDx) ? rawDy : rawDx;
    const dStops = Math.max(-0.55, Math.min(0.55, delta / 280));

    if (view.getOptions().mode === "scroll") {
      dialTarget = Math.max(0, dialTarget + dStops);
    } else {
      dialTarget += dStops;
    }
    lastInputAt = performance.now();
    view.nudgeSpin(dStops * 0.42, 0);
    view.kick();
  }

  function onKeyDown(e) {
    if (e.key === "ArrowDown" || e.key === "ArrowRight") {
      e.preventDefault();
      dialTarget = Math.round(dialTarget) + 1;
      lastInputAt = performance.now();
      view.nudgeSpin(0.35, 0);
      view.kick();
    } else if (e.key === "ArrowUp" || e.key === "ArrowLeft") {
      e.preventDefault();
      const next = Math.round(dialTarget) - 1;
      dialTarget = view.getOptions().mode === "scroll" ? Math.max(0, next) : next;
      lastInputAt = performance.now();
      view.nudgeSpin(-0.35, 0);
      view.kick();
    }
  }

  // Full-page inertial touch scrolling across the entire viewport (including the protein stage)
  let touchY = null;
  let touchLastT = 0;
  let touchVelStops = 0;

  function onTouchStart(e) {
    if (e.touches.length !== 1) {
      touchY = null;
      isTouching = false;
      return;
    }
    touchY = e.touches[0].clientY;
    touchLastT = performance.now();
    touchVelStops = 0;
    isTouching = true;
  }

  function onTouchMove(e) {
    if (touchY === null || e.touches.length !== 1) return;
    e.preventDefault();
    const now = performance.now();
    const dt = Math.max(0.004, (now - touchLastT) / 1000);
    const curY = e.touches[0].clientY;
    const dy = touchY - curY; // swipe up (scroll down) -> dy > 0
    touchY = curY;
    touchLastT = now;

    if (Math.abs(dy) < 0.2) return;
    const dStops = Math.max(-0.45, Math.min(0.45, dy / 220));
    const instVel = dStops / dt;
    touchVelStops = touchVelStops * 0.45 + instVel * 0.55;

    dialTarget = Math.max(0, dialTarget + dStops);
    lastInputAt = now;
    view.nudgeSpin(dStops * 0.42, 0);
    view.kick();
  }

  function onTouchEnd() {
    if (!isTouching) return;
    isTouching = false;
    touchY = null;
    const now = performance.now();
    // Apply smooth flick coasting if released shortly after moving
    if (now - touchLastT < 90 && Math.abs(touchVelStops) > 0.25) {
      const coast = Math.max(-1.15, Math.min(1.15, touchVelStops * 0.14));
      dialTarget = Math.max(0, dialTarget + coast);
      lastInputAt = now;
      view.kick();
    }
  }

  window.addEventListener("wheel", onWheel, { passive: false });
  window.addEventListener("keydown", onKeyDown);
  window.addEventListener("touchstart", onTouchStart, { passive: true });
  window.addEventListener("touchmove", onTouchMove, { passive: false });
  window.addEventListener("touchend", onTouchEnd, { passive: true });
  window.addEventListener("touchcancel", onTouchEnd, { passive: true });

  return {
    setStop(stopVal) {
      dialTarget = stopVal;
      dialSmooth = stopVal;
      lastInputAt = performance.now();
      onTick(0.016, performance.now());
      view.render();
      view.kick();
    },
    dispose() {
      window.removeEventListener("wheel", onWheel);
      window.removeEventListener("keydown", onKeyDown);
      window.removeEventListener("touchstart", onTouchStart);
      window.removeEventListener("touchmove", onTouchMove);
      window.removeEventListener("touchend", onTouchEnd);
      window.removeEventListener("touchcancel", onTouchEnd);
    },
  };
}
