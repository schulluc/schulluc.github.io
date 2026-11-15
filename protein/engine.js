import { add, mul, dot, norm } from "./vec.js";
import { buildGeometry, SHAPE_DEFAULT } from "./geometry.js";
import { parseStructureText, fetchRCSB, fetchAlphaFold } from "./parse.js";

// 8x8 Bayer ordered halftone matrix (normalized 0..1)
const BAYER8 = [
   0, 48, 12, 60,  3, 51, 15, 63,
  32, 16, 44, 28, 35, 19, 47, 31,
   8, 56,  4, 52, 11, 59,  7, 55,
  40, 24, 36, 20, 43, 27, 39, 23,
   2, 50, 14, 62,  1, 49, 13, 61,
  34, 18, 46, 30, 33, 17, 45, 29,
  10, 58,  6, 54,  9, 57,  5, 53,
  42, 26, 38, 22, 41, 25, 37, 21,
].map(v => (v + 0.5) / 64);

const LIGHT = norm([-0.42, 0.56, 0.71]);
const HALFTONE_RAMP = " .,:;-=+*x#%@";
const STIPPLE_RAMP = " ·.:;+x#";

export const DEFAULTS = {
  paper: "#ffffff",
  ink: "#1a1a1a",
  transparent: false,

  // Motion & camera
  period: 16,           // seconds per full 360° turntable revolution (0 = static)
  tiltX: -0.56,          // +43.5° tilts the (β/α)8 barrel mouth toward camera like Richardson Fig. 1b
  rotY: 0.42,
  zoom: 1.0,
  interactive: true,
  dragOrbit: true,
  touchOrbit: true,
  externalScroll: false,

  // Option 2C interactive blend modes:
  // "lens" | "scroll" | "both" | "curtain" | "split" | "ink" | "ascii"
  mode: "lens",
  invert: false,
  lensRadius: 0.30,
  splitPos: 0.50,
  scrollRange: 220,
  scrollMax: 0.40,          // fraction of protein width drawn as ASCII left-to-right at full scroll

  // Aesthetic controls
  asciiStyle: "richardson", // "richardson" | "hybrid" | "halftone" | "stipple"
  detail: 70,
  lineWeight: 1.32,
  haloWidth: 4.2,
  depthFade: 0.0,           // 0..0.85 atmospheric ink fade toward paper for huge complexes
  hatchHelix: true,
  slabStrand: true,
  asciiHalo: true,

  font: '"JetBrains Mono", "SF Mono", "Menlo", "IBM Plex Mono", ui-monospace, monospace',
};

function parseHexRgb(hex, fallback) {
  const s = String(hex || "").trim().replace("#", "");
  if (s.length === 3) {
    return [parseInt(s[0] + s[0], 16), parseInt(s[1] + s[1], 16), parseInt(s[2] + s[2], 16)];
  }
  if (s.length === 6) {
    return [parseInt(s.slice(0, 2), 16), parseInt(s.slice(2, 4), 16), parseInt(s.slice(4, 6), 16)];
  }
  return fallback;
}

function buildInkPalette(inkHex, paperHex) {
  const ink = parseHexRgb(inkHex, [26, 26, 26]);
  const pap = parseHexRgb(paperHex, [255, 255, 255]);
  const pal = new Array(32);
  for (let i = 0; i < 32; i++) {
    const t = i / 31; // 0 = paper, 1 = full ink
    const r = Math.round(pap[0] + (ink[0] - pap[0]) * t);
    const g = Math.round(pap[1] + (ink[1] - pap[1]) * t);
    const b = Math.round(pap[2] + (ink[2] - pap[2]) * t);
    pal[i] = `rgb(${r},${g},${b})`;
  }
  return pal;
}

function distToSegSq(px, py, ax, ay, bx, by) {
  const vx = bx - ax, vy = by - ay;
  const wx = px - ax, wy = py - ay;
  const l2 = vx * vx + vy * vy;
  if (l2 < 1e-8) return wx * wx + wy * wy;
  let t = (wx * vx + wy * vy) / l2;
  if (t < 0) t = 0;
  else if (t > 1) t = 1;
  const dx = px - (ax + t * vx), dy = py - (ay + t * vy);
  return dx * dx + dy * dy;
}

function pointInQuad(px, py, x0, y0, x1, y1, x2, y2, x3, y3) {
  const c0 = (x1 - x0) * (py - y0) - (y1 - y0) * (px - x0);
  const c1 = (x2 - x1) * (py - y1) - (y2 - y1) * (px - x1);
  const c2 = (x3 - x2) * (py - y2) - (y3 - y2) * (px - x2);
  const c3 = (x0 - x3) * (py - y3) - (y0 - y3) * (px - x3);
  return (c0 >= 0 && c1 >= 0 && c2 >= 0 && c3 >= 0) || (c0 <= 0 && c1 <= 0 && c2 <= 0 && c3 <= 0);
}

function dirGlyph(angleRad, cue) {
  let a = angleRad % Math.PI;
  if (a < 0) a += Math.PI;
  if (cue < 0.22) return ":";
  if (a < 0.39 || a >= 2.75) return "-";
  if (a < 1.18) return "\\";
  if (a < 1.96) return "|";
  return "/";
}

const lerp2 = (a, b, t) => [a[0] + (b[0] - a[0]) * t, a[1] + (b[1] - a[1]) * t];

export function createProteinFold(container, options = {}) {
  const reduceMotion =
    typeof matchMedia === "function" && matchMedia("(prefers-reduced-motion: reduce)").matches;
  const o = {
    ...DEFAULTS,
    ...(reduceMotion && options.period === undefined ? { period: 0 } : {}),
    ...options,
  };
  let shape = { ...SHAPE_DEFAULT, ...(options.shape || {}) };
  const listeners = {};
  const emit = (ev, data) => (listeners[ev] || []).forEach(fn => fn(data));

  if (getComputedStyle(container).position === "static") {
    container.style.position = "relative";
  }

  const root = document.createElement("div");
  root.className = "protein-fold-stage";
  root.style.cssText =
    "position:absolute;inset:0;overflow:hidden;display:flex;align-items:center;justify-content:center;user-select:none;-webkit-user-select:none;touch-action:pinch-zoom;";
  root.setAttribute("role", "img");

  const inner = document.createElement("div");
  inner.style.cssText = "position:relative;width:max-content;line-height:1;";

  const rib = document.createElement("canvas");
  rib.style.cssText = "position:absolute;left:0;top:0;pointer-events:none;";
  rib.setAttribute("aria-hidden", "true");

  const out = document.createElement("pre");
  out.style.cssText =
    "margin:0;padding:0;white-space:pre;line-height:1;font-weight:600;position:relative;letter-spacing:0;";
  out.setAttribute("aria-hidden", "true");

  const probe = document.createElement("span");
  probe.textContent = "MMMMMMMMMM";
  probe.style.cssText =
    "position:absolute;visibility:hidden;white-space:pre;line-height:1;font-weight:600;left:-9999px;top:0;";

  inner.append(rib, out);
  root.append(inner, probe);
  container.append(root);

  const rctx = rib.getContext("2d");
  const maskCanvas = document.createElement("canvas");
  const mctx = maskCanvas.getContext("2d");

  let model = null;
  let G = null;
  let COLS = 108, ROWS = 70, cw = 7.2, lh = 12, aspect = 0.6;
  let rotY = o.rotY;
  let rotX = o.tiltX;
  let zoom = o.zoom;
  let panX = 0, panY = 0;

  let pointerInside = false;
  let targetMX = 0.5, targetMY = 0.5;
  let smoothMX = 0.5, smoothMY = 0.5;
  let hoverAmt = 0;
  let scrollAmt = 0;
  let smoothSpanX0 = -1, smoothSpanX1 = -1;

  let mimg = null;
  let ditherBuf = null;
  let zbuf = null;
  let kindBuf = null;   // 1:coil, 2:H-front, 3:H-back, 4:E-front, 5:E-slab
  let edgeBuf = null;
  let angBuf = null;
  let shadeBuf = null;
  let cueBuf = null;
  let occBuf = null;
  let inkPalette = buildInkPalette(o.ink, o.paper);
  let lastText = "";

  function alloc() {
    maskCanvas.width = COLS;
    maskCanvas.height = ROWS;
    mimg = mctx.createImageData(COLS, ROWS);
    const n = COLS * ROWS;
    ditherBuf = new Float32Array(n);
    for (let j = 0; j < ROWS; j++) {
      for (let i = 0; i < COLS; i++) {
        const b = BAYER8[(j & 7) * 8 + (i & 7)];
        const jitter = ((Math.sin(i * 12.9898 + j * 78.233) * 43758.5453) % 1 + 1) % 1;
        ditherBuf[j * COLS + i] = b * 0.84 + jitter * 0.16;
      }
    }
    zbuf = new Float32Array(n);
    kindBuf = new Uint8Array(n);
    edgeBuf = new Uint8Array(n);
    angBuf = new Float32Array(n);
    shadeBuf = new Float32Array(n);
    cueBuf = new Float32Array(n);
    occBuf = new Uint8Array(n);
  }

  function applyStyle() {
    out.style.color = o.ink;
    out.style.fontFamily = probe.style.fontFamily = o.font;
    root.style.background = o.transparent ? "transparent" : o.paper;
    root.style.cursor = o.interactive && o.dragOrbit ? "grab" : "default";
    inkPalette = buildInkPalette(o.ink, o.paper);
  }

  function fit() {
    const W = root.clientWidth;
    const H = root.clientHeight;
    if (!W || !H) return;
    probe.style.fontSize = "100px";
    const ratio = probe.getBoundingClientRect().width / 10 / 100 || 0.6;
    const targetDetail = Math.max(24, Math.min(130, Math.round(o.detail)));
    const effH = Math.min(H, W * 1.08);
    let fs = effH / targetDetail;
    let rows = Math.max(24, Math.floor(H / fs));
    while (rows > 24 && rows * Math.floor(W / (fs * ratio)) > 65000) {
      fs *= 1.04;
      rows = Math.max(24, Math.floor(H / fs));
    }
    probe.style.fontSize = fs + "px";
    cw = probe.getBoundingClientRect().width / 10 || fs * ratio;
    lh = fs;
    aspect = cw / lh;
    const cols = Math.max(20, Math.floor(W / cw));
    if (cols !== COLS || rows !== ROWS || !zbuf) {
      COLS = cols;
      ROWS = rows;
      smoothSpanX0 = -1;
      smoothSpanX1 = -1;
      alloc();
    }
    out.style.fontSize = fs + "px";
    const dpr = Math.min(2.5, window.devicePixelRatio || 1);
    rib.style.width = COLS * cw + "px";
    rib.style.height = ROWS * lh + "px";
    rib.width = Math.round(COLS * cw * dpr);
    rib.height = Math.round(ROWS * lh * dpr);
  }

  const baseScale = () =>
    G ? Math.min((COLS * 0.46) / G.R, (ROWS * 0.46) / (G.R * aspect)) : 1;

  function lensFactorAt(i, j) {
    if (hoverAmt <= 0.002) return 0;
    const dx = (i + 0.5) * cw - smoothMX * (COLS * cw);
    const dy = (j + 0.5) * lh - smoothMY * (ROWS * lh);
    const distPx = Math.hypot(dx, dy);
    const minDim = Math.min(COLS * cw, ROWS * lh);
    const rad = minDim * o.lensRadius * (0.35 + 0.65 * hoverAmt);
    const feather = Math.max(20, rad * 0.46);
    const raw = 1 - (distPx - (rad - feather * 0.5)) / feather;
    return Math.max(0, Math.min(1, raw)) * hoverAmt;
  }

  function scrollLRFactorAt(i, box) {
    if (scrollAmt <= 0.001) return 0;
    const x0 = smoothSpanX0 >= 0 ? smoothSpanX0 : box.x0;
    const x1 = smoothSpanX1 > x0 ? smoothSpanX1 : box.x1;
    const u = (i + 0.5 - x0) / Math.max(1, x1 - x0);
    const maxFrac = o.scrollMax !== undefined ? o.scrollMax : 0.40;
    const feather = 0.09;
    const frontier = -feather * 0.5 + scrollAmt * (maxFrac + feather * 0.5);
    const raw = (frontier + feather * 0.5 - u) / feather;
    return Math.max(0, Math.min(1, raw));
  }

  function asciiTargetAt(i, j, box) {
    let f = 0;
    if (o.mode === "ascii") {
      f = 1;
    } else if (o.mode === "ink") {
      f = 0;
    } else if (o.mode === "split") {
      const x0 = smoothSpanX0 >= 0 ? smoothSpanX0 : box.x0;
      const x1 = smoothSpanX1 > x0 ? smoothSpanX1 : box.x1;
      const u = (i + 0.5 - x0) / Math.max(1, x1 - x0);
      f = Math.max(0, Math.min(1, (o.splitPos + 0.08 - u) / 0.16));
    } else if (o.mode === "lens") {
      f = lensFactorAt(i, j);
    } else if (o.mode === "scroll") {
      f = scrollLRFactorAt(i, box);
    } else if (o.mode === "both") {
      f = Math.max(scrollLRFactorAt(i, box), lensFactorAt(i, j));
    } else if (o.mode === "curtain") {
      if (hoverAmt > 0.002) {
        const u = (i + 0.5 - box.x0) / Math.max(1, box.x1 - box.x0);
        const v = (j + 0.5 - box.y0) / Math.max(1, box.y1 - box.y0);
        const diag = u * 0.82 + v * 0.18;
        const sweep = hoverAmt * 1.35 - 0.18;
        f = Math.max(0, Math.min(1, (sweep + 0.18 - diag) / 0.36));
      }
    }
    if (f > 0 && f < 1) f = f * f * (3 - 2 * f);
    return o.invert ? 1 - f : f;
  }

  function rasterizeQuadToGrid(
    cx0, cy0, cx1, cy1, cx2, cy2, cx3, cy3,
    zVal, cue, shade, kindCode, railAngle, capStart, capShoulder
  ) {
    const pad = 0.52;
    const i0 = Math.max(0, Math.floor(Math.min(cx0, cx1, cx2, cx3) - pad));
    const i1 = Math.min(COLS - 1, Math.ceil(Math.max(cx0, cx1, cx2, cx3) + pad));
    const j0 = Math.max(0, Math.floor(Math.min(cy0, cy1, cy2, cy3) - pad));
    const j1 = Math.min(ROWS - 1, Math.ceil(Math.max(cy0, cy1, cy2, cy3) + pad));

    const edgeThreshSq = kindCode === 1 ? 0.25 : 0.34;

    for (let j = j0; j <= j1; j++) {
      const py = j + 0.5;
      for (let i = i0; i <= i1; i++) {
        const px = i + 0.5;
        const dLeftSq = distToSegSq(px, py, cx0, cy0, cx3, cy3);
        const dRightSq = distToSegSq(px, py, cx1, cy1, cx2, cy2);
        const minRailSq = dLeftSq < dRightSq ? dLeftSq : dRightSq;

        const inside = kindCode !== 1 && pointInQuad(px, py, cx0, cy0, cx1, cy1, cx2, cy2, cx3, cy3);
        if (!inside && minRailSq > (kindCode === 1 ? 0.23 : 0.19)) continue;

        const idx = j * COLS + i;
        const isEdge =
          kindCode === 1 ||
          minRailSq <= edgeThreshSq ||
          (capStart && distToSegSq(px, py, cx0, cy0, cx1, cy1) <= 0.28) ||
          (capShoulder && distToSegSq(px, py, cx0, cy0, cx1, cy1) <= 0.30);

        const zEff = zVal + (isEdge ? 0.08 : 0);
        if (zEff < zbuf[idx]) continue;

        zbuf[idx] = zEff;
        kindBuf[idx] = kindCode;
        edgeBuf[idx] = isEdge ? 1 : 0;
        angBuf[idx] = railAngle;
        shadeBuf[idx] = shade;
        cueBuf[idx] = cue;
      }
    }
  }

  function pickAsciiChar(idx) {
    const kind = kindBuf[idx];
    if (!kind || occBuf[idx]) return " ";
    const isEdge = edgeBuf[idx] === 1;
    const rawCue = cueBuf[idx];
    const cue = o.depthFade > 0 ? rawCue * (1 - o.depthFade * 0.45 * (1 - rawCue)) : rawCue;
    const shade = shadeBuf[idx];

    if (o.asciiStyle === "richardson") {
      if (isEdge || kind === 5) {
        if (o.depthFade > 0.25 && rawCue < 0.28) return "·";
        return dirGlyph(angBuf[idx], cue);
      }
      if (kind === 3) {
        // Concave inner curl of α-helix -> parallel Rapidograph line-shading glyphs
        if (o.depthFade > 0.25 && rawCue < 0.34) return "·";
        return cue > 0.58 && shade < 0.48 ? "=" : ":";
      }
      // Clean airy white paper fronts on strands and outer helices
      return " ";
    }

    if (o.asciiStyle === "hybrid") {
      if (isEdge || kind === 5) {
        return dirGlyph(angBuf[idx], cue);
      }
      if (kind === 3) {
        return shade < 0.45 ? "=" : ":";
      }
      const d = (1 - shade) * (0.35 + 0.65 * cue);
      if (d > 0.52) return ":";
      if (d > 0.28) return "·";
      return " ";
    }

    if (o.asciiStyle === "stipple") {
      const ramp = STIPPLE_RAMP;
      const n = ramp.length - 1;
      let darkness = isEdge ? 0.88 : kind === 3 || kind === 5 ? 0.58 : 0.16 + (1 - shade) * 0.26;
      darkness *= 0.38 + 0.62 * cue;
      return ramp[Math.max(1, Math.min(n, Math.round(darkness * n)))];
    }

    // "halftone"
    const ramp = HALFTONE_RAMP;
    const n = ramp.length - 1;
    let darkness =
      kind === 5
        ? 0.84
        : kind === 3
        ? 0.56 + 0.28 * (1 - shade)
        : isEdge
        ? 0.76
        : 0.18 + 0.54 * (1 - shade);
    darkness *= 0.32 + 0.68 * cue;
    return ramp[Math.max(1, Math.min(n, Math.round(darkness * n)))];
  }

  function strokeMainRails(p) {
    const q = p.q;
    rctx.moveTo(q[0][0], q[0][1]);
    rctx.lineTo(q[3][0], q[3][1]);
    rctx.moveTo(q[1][0], q[1][1]);
    rctx.lineTo(q[2][0], q[2][1]);
    if (p.strandStart) {
      rctx.moveTo(q[0][0], q[0][1]);
      rctx.lineTo(q[1][0], q[1][1]);
    }
    if (p.shoulderWings) {
      for (const w of p.shoulderWings) {
        rctx.moveTo(w[0], w[1]);
        rctx.lineTo(w[2], w[3]);
      }
    }
  }

  function render() {
    if (!G || !zbuf) return;
    const { pts, R } = G;
    const nPts = pts.length;

    const cy = Math.cos(rotY), sy = Math.sin(rotY);
    const cx = Math.cos(rotX), sx = Math.sin(rotX);
    const rot = v => {
      const x = v[0] * cy + v[2] * sy;
      const z0 = -v[0] * sy + v[2] * cy;
      return [x, v[1] * cx - z0 * sx, v[1] * sx + z0 * cx];
    };

    const S = baseScale() * zoom;
    const ox = COLS * 0.5;
    const oy = ROWS * 0.5;

    const dpr = rib.width / (COLS * cw);
    const kpx = S * cw * dpr;
    const OX = ox * cw * dpr;
    const OY = oy * lh * dpr;

    const cellX = v => ox + (v[0] + panX) * S;
    const cellY = v => oy - (v[1] + panY) * S * aspect;
    const pixX = v => OX + (v[0] + panX) * kpx;
    const pixY = v => OY - (v[1] + panY) * kpx;

    const box = {
      x0: Math.max(0, ox + (panX - R) * S),
      x1: Math.min(COLS, ox + (panX + R) * S),
      y0: Math.max(0, oy - (panY + R) * S * aspect),
      y1: Math.min(ROWS, oy - (panY - R) * S * aspect),
    };

    const P = new Array(nPts);
    const SN = new Array(nPts);
    const NN = new Array(nPts);
    for (let k = 0; k < nPts; k++) {
      P[k] = rot(pts[k].p);
      SN[k] = rot(pts[k].side);
      NN[k] = rot(pts[k].nrm);
    }

    const ropePerpPx = new Array(nPts);
    const ropePerpCell = new Array(nPts);
    const rCoilPx = G.shape.coilR * kpx;
    const rCoilCellX = G.shape.coilR * S;
    const rCoilCellY = G.shape.coilR * S * aspect;
    for (let k = 0; k < nPts; k++) {
      const prevK = k > 0 && !pts[k - 1].gap ? k - 1 : k;
      const nextK = k < nPts - 1 && !pts[k].gap ? k + 1 : k;
      let dx = pixX(P[nextK]) - pixX(P[prevK]);
      let dy = pixY(P[nextK]) - pixY(P[prevK]);
      const dl = Math.hypot(dx, dy) || 1;
      dx /= dl; dy /= dl;
      const nx = -dy, ny = dx;
      ropePerpPx[k] = [nx * rCoilPx, ny * rCoilPx];
      ropePerpCell[k] = [nx * rCoilCellX, ny * rCoilCellY];
    }

    zbuf.fill(-1e9);
    kindBuf.fill(0);
    edgeBuf.fill(0);
    occBuf.fill(0);

    const polys = [];
    let prevMainPoly = null;
    const needInk = (o.mode !== "ascii" || o.invert) && !o.skipDefaultInk;
    const needCustomInk = typeof o.onBeforeMask === "function" && (o.mode !== "ascii" || o.invert);
    const needAscii = o.mode !== "ink" || o.invert;

    for (let k = 0; k < nPts - 1; k++) {
      const a = pts[k];
      if (a.gap) {
        prevMainPoly = null;
        continue;
      }
      const b = pts[k + 1];
      const A = P[k], B = P[k + 1];
      const zc = (A[2] + B[2]) * 0.5;
      const cue = Math.max(0, Math.min(1, (zc + R) / (2 * R)));
      const railAngle = Math.atan2(pixY(B) - pixY(A), pixX(B) - pixX(A));

      if (a.kind === "C" || a.wOut <= 0.01) {
        const ax = pixX(A), ay = pixY(A), bx = pixX(B), by = pixY(B);
        const pa = ropePerpPx[k], pb = ropePerpPx[k + 1];
        if (needInk) {
          const poly = {
            z: zc,
            type: "rope",
            cue,
            prev: prevMainPoly,
            next: null,
            q: [
              [ax + pa[0], ay + pa[1]],
              [ax - pa[0], ay - pa[1]],
              [bx - pb[0], by - pb[1]],
              [bx + pb[0], by + pb[1]],
            ],
          };
          if (prevMainPoly) prevMainPoly.next = poly;
          prevMainPoly = poly;
          polys.push(poly);
        }
        if (needAscii) {
          const cax = cellX(A), cay = cellY(A), cbx = cellX(B), cby = cellY(B);
          const cpa = ropePerpCell[k], cpb = ropePerpCell[k + 1];
          rasterizeQuadToGrid(
            cax + cpa[0], cay + cpa[1],
            cax - cpa[0], cay - cpa[1],
            cbx - cpb[0], cby - cpb[1],
            cbx + cpb[0], cby + cpb[1],
            zc, cue, 0.55, 1, railAngle, false, false
          );
        }
        continue;
      }

      const wa = a.wOut;
      const wb = b.res === a.res || b.kind === a.kind ? b.wIn : G.shape.coilR;
      const sa = SN[k], sb = SN[k + 1];
      const na = NN[k], nb = NN[k + 1];
      const nf = norm(add(na, nb));
      const shade = Math.abs(dot(nf, LIGHT));

      // Cylindrical radial normal has helixSign=+1 outward; < -0.05 means concave inner curl!
      const isHelixBack = a.kind === "H" && nf[2] * a.helixSign < -0.05;

      const hasSlab = o.slabStrand && a.kind === "E" && a.th > 0.02 && a.th * kpx >= 0.85 * dpr;
      const slabDir = nf[2] >= 0 ? 1 : -1;
      const topOffA = hasSlab ? mul(na, slabDir * a.th) : [0, 0, 0];
      const topOffB = hasSlab ? mul(nb, slabDir * (b.kind === "E" ? b.th : a.th)) : [0, 0, 0];

      const c1 = add(add(A, topOffA), mul(sa, wa));
      const c2 = add(add(A, topOffA), mul(sa, -wa));
      const c3 = add(add(B, topOffB), mul(sb, -wb));
      const c4 = add(add(B, topOffB), mul(sb, wb));

      let shoulderWings = null;
      if (a.arrowShoulder) {
        const sIn1 = add(add(A, topOffA), mul(sa, a.shaftW));
        const sIn2 = add(add(A, topOffA), mul(sa, -a.shaftW));
        shoulderWings = [
          [pixX(sIn1), pixY(sIn1), pixX(c1), pixY(c1)],
          [pixX(sIn2), pixY(sIn2), pixX(c2), pixY(c2)],
        ];
      }

      if (needInk) {
        if (hasSlab && Math.abs(nf[2]) < 0.94) {
          const botOffA = mul(na, -slabDir * a.th);
          const botOffB = mul(nb, -slabDir * (b.kind === "E" ? b.th : a.th));
          for (const sg of [1, -1]) {
            const sideViewZ = (sa[2] + sb[2]) * 0.5 * sg;
            if (sideViewZ < -0.10) continue;
            const tA = sg > 0 ? c1 : c2;
            const tB = sg > 0 ? c4 : c3;
            const bA = add(add(A, botOffA), mul(sa, sg * wa));
            const bB = add(add(B, botOffB), mul(sb, sg * wb));
            polys.push({
              z: zc - 0.08,
              type: "slab",
              cue,
              q: [
                [pixX(tA), pixY(tA)],
                [pixX(bA), pixY(bA)],
                [pixX(bB), pixY(bB)],
                [pixX(tB), pixY(tB)],
              ],
            });
          }
        }

        const mainPoly = {
          z: zc,
          type: "ribbon",
          kind: a.kind,
          cue,
          hatch: o.hatchHelix && isHelixBack,
          strandStart: !!a.strandStart,
          shoulderWings,
          prev: prevMainPoly,
          next: null,
          q: [
            [pixX(c1), pixY(c1)],
            [pixX(c2), pixY(c2)],
            [pixX(c3), pixY(c3)],
            [pixX(c4), pixY(c4)],
          ],
        };
        if (prevMainPoly) prevMainPoly.next = mainPoly;
        prevMainPoly = mainPoly;
        polys.push(mainPoly);
      }

      if (needAscii) {
        const ascKind = a.kind === "H" ? (isHelixBack ? 3 : 2) : 4;
        rasterizeQuadToGrid(
          cellX(c1), cellY(c1),
          cellX(c2), cellY(c2),
          cellX(c3), cellY(c3),
          cellX(c4), cellY(c4),
          zc, cue, shade, ascKind, railAngle, !!a.strandStart, !!a.arrowShoulder
        );
        if (hasSlab && Math.abs(nf[2]) < 0.90) {
          const botOffA = mul(na, -slabDir * a.th * 1.25);
          const botOffB = mul(nb, -slabDir * (b.kind === "E" ? b.th : a.th) * 1.25);
          const sg = sa[2] + sb[2] >= 0 ? 1 : -1;
          const bA = add(add(A, botOffA), mul(sa, sg * wa));
          const bB = add(add(B, botOffB), mul(sb, sg * wb));
          const tA = sg > 0 ? c1 : c2;
          const tB = sg > 0 ? c4 : c3;
          rasterizeQuadToGrid(
            cellX(tA), cellY(tA),
            cellX(bA), cellY(bA),
            cellX(bB), cellY(bB),
            cellX(tB), cellY(tB),
            zc - 0.05, cue, 0.25, 5, railAngle, false, false
          );
        }
      }
    }

    if (needAscii && o.asciiHalo) {
      const zGap = Math.max(3.4, R * 0.11);
      for (let j = 1; j < ROWS - 1; j++) {
        for (let i = 1; i < COLS - 1; i++) {
          const idx = j * COLS + i;
          if (!kindBuf[idx]) continue;
          const z0 = zbuf[idx] + zGap;
          if (
            zbuf[idx - 1] > z0 ||
            zbuf[idx + 1] > z0 ||
            zbuf[idx - COLS] > z0 ||
            zbuf[idx + COLS] > z0
          ) {
            occBuf[idx] = 1;
          }
        }
      }
    }

    if (needAscii) {
      let minC = COLS, maxC = -1;
      for (let j = 0; j < ROWS; j++) {
        const rowOff = j * COLS;
        for (let i = 0; i < COLS; i++) {
          if (kindBuf[rowOff + i]) {
            if (i < minC) minC = i;
            if (i > maxC) maxC = i;
          }
        }
      }
      if (maxC >= minC) {
        if (smoothSpanX0 < 0) {
          smoothSpanX0 = minC;
          smoothSpanX1 = maxC + 1;
        } else {
          smoothSpanX0 += (minC - smoothSpanX0) * 0.22;
          smoothSpanX1 += (maxC + 1 - smoothSpanX1) * 0.22;
        }
      }
    }

    const md = mimg.data;
    const lines = new Array(ROWS);
    for (let j = 0; j < ROWS; j++) {
      let rowStr = "";
      for (let i = 0; i < COLS; i++) {
        const idx = j * COLS + i;
        const f = asciiTargetAt(i, j, box);
        const showAscii = f >= 0.999 ? true : f <= 0.001 ? false : ditherBuf[idx] < f;
        md[idx * 4] = 0;
        md[idx * 4 + 1] = 0;
        md[idx * 4 + 2] = 0;
        md[idx * 4 + 3] = showAscii ? 0 : 255;
        rowStr += showAscii && kindBuf[idx] ? pickAsciiChar(idx) : " ";
      }
      lines[j] = rowStr;
    }
    const fullText = lines.join("\n");
    if (fullText !== lastText) {
      out.textContent = fullText;
      lastText = fullText;
    }

    // Draw Jane S. Richardson pen-and-ink layer on <canvas>
    rctx.setTransform(1, 0, 0, 1, 0, 0);
    rctx.clearRect(0, 0, rib.width, rib.height);

    if (needInk && polys.length > 0) {
      polys.sort((p, q) => p.z - q.z);
      const rScale = Math.min(1, Math.pow(27 / Math.max(18, R), 0.50));
      const baseLw = o.lineWeight * dpr * rScale;
      const haloLw = o.haloWidth * dpr * rScale;
      const df = Math.max(0, Math.min(0.95, o.depthFade || 0));

      for (let m = 0; m < polys.length; m++) {
        const p = polys[m];
        const q = p.q;
        const lw = baseLw * (0.72 + 0.45 * p.cue);
        const back = Math.max(0, (0.66 - p.cue) / 0.66);
        const inkStrength = df > 0 ? Math.max(0.08, 1 - df * Math.pow(back, 1.15)) : 1;
        const strokeInk = inkStrength >= 0.99 ? o.ink : inkPalette[Math.round(inkStrength * 31)];

        if (p.type === "slab") {
          // Crisp architectural white side-slab bounded by fine ink borders (Richardson Fig. 1b/1c)
          rctx.fillStyle = o.paper;
          rctx.strokeStyle = o.paper;
          rctx.lineWidth = 0.8 * dpr;
          rctx.lineJoin = "round";
          rctx.beginPath();
          rctx.moveTo(q[0][0], q[0][1]);
          rctx.lineTo(q[1][0], q[1][1]);
          rctx.lineTo(q[2][0], q[2][1]);
          rctx.lineTo(q[3][0], q[3][1]);
          rctx.closePath();
          rctx.fill();
          rctx.stroke();

          rctx.strokeStyle = strokeInk;
          rctx.lineWidth = lw * 0.78;
          rctx.lineCap = "round";
          rctx.beginPath();
          rctx.moveTo(q[1][0], q[1][1]);
          rctx.lineTo(q[2][0], q[2][1]);
          rctx.stroke();
          continue;
        }

        // 1. White crossing-gap halo along lateral rails (breaks background ink lines cleanly)
        if (haloLw > 0.2) {
          rctx.strokeStyle = o.paper;
          rctx.lineWidth = haloLw;
          rctx.lineCap = "butt";
          rctx.beginPath();
          rctx.moveTo(q[0][0], q[0][1]);
          rctx.lineTo(q[3][0], q[3][1]);
          rctx.moveTo(q[1][0], q[1][1]);
          rctx.lineTo(q[2][0], q[2][1]);
          rctx.stroke();
        }

        // 2. Crisp opaque paper body fill
        rctx.fillStyle = o.paper;
        rctx.strokeStyle = o.paper;
        rctx.lineWidth = 0.9 * dpr * rScale;
        rctx.lineJoin = "round";
        rctx.beginPath();
        rctx.moveTo(q[0][0], q[0][1]);
        rctx.lineTo(q[1][0], q[1][1]);
        rctx.lineTo(q[2][0], q[2][1]);
        rctx.lineTo(q[3][0], q[3][1]);
        rctx.closePath();
        rctx.fill();
        rctx.stroke();

        // 3. Longitudinal spiral engraving tracks on concave inside of α-helix curls
        if (p.hatch) {
          rctx.strokeStyle = strokeInk;
          rctx.globalAlpha = (0.52 + 0.26 * p.cue) * (0.55 + 0.45 * inkStrength);
          rctx.lineWidth = Math.max(0.55 * dpr, lw * 0.50);
          rctx.lineCap = "round";
          rctx.beginPath();
          for (const t of [0.27, 0.50, 0.73]) {
            const s0 = lerp2(q[0], q[1], t);
            const s1 = lerp2(q[3], q[2], t);
            rctx.moveTo(s0[0], s0[1]);
            rctx.lineTo(s1[0], s1[1]);
          }
          rctx.stroke();
          rctx.globalAlpha = 1;
        }

        // 4. Primary Rapidograph ink rail outlines (+ 1-slice neighbor bridge so own halo never nicks own curve!)
        rctx.strokeStyle = strokeInk;
        rctx.lineWidth = p.type === "rope" ? lw * 0.86 : lw;
        rctx.lineCap = "round";
        rctx.lineJoin = "round";
        rctx.beginPath();
        strokeMainRails(p);
        if (p.prev && p.z >= p.prev.z) strokeMainRails(p.prev);
        if (p.next && p.z >= p.next.z) strokeMainRails(p.next);
        rctx.stroke();
      }
    }

    if (needCustomInk) {
      o.onBeforeMask({
        rctx,
        width: rib.width,
        height: rib.height,
        dpr,
        rotY,
        rotX,
        panX,
        panY,
        zoom,
        S,
        kpx,
        OX,
        OY,
        R,
        rot,
        pixX,
        pixY,
        model,
        G,
        inkPalette,
        options: o,
        interacting: dragging || Math.abs(velY) > 0.6,
      });
    }

    if (needInk || needCustomInk) {
      mctx.putImageData(mimg, 0, 0);
      rctx.globalCompositeOperation = "destination-in";
      rctx.imageSmoothingEnabled = false;
      rctx.drawImage(maskCanvas, 0, 0, rib.width, rib.height);
      rctx.globalCompositeOperation = "source-over";
    }
  }

  let raf = 0;
  let prev = performance.now();
  let visible = true;
  let dragging = false;
  let destroyed = false;
  let velY = 0, velX = 0;

  function updateScroll() {
    if (typeof window === "undefined") return;
    const doc = document.documentElement;
    const maxPageScroll = Math.max(1, (doc ? doc.scrollHeight : 0) - window.innerHeight);
    const effRange = Math.max(40, Math.min(o.scrollRange, maxPageScroll * 0.85));
    scrollAmt = Math.max(0, Math.min(1, window.scrollY / effRange));
  }

  function loop(now) {
    raf = 0;
    if (destroyed) return;
    const dt = Math.min(0.05, Math.max(0.001, (now - prev) / 1000));
    prev = now;

    const targetHover = pointerInside ? 1 : 0;
    const rate = 1 - Math.exp(-dt * (pointerInside ? 9.5 : 5.5));
    hoverAmt += (targetHover - hoverAmt) * rate;
    const posRate = 1 - Math.exp(-dt * 14);
    smoothMX += (targetMX - smoothMX) * posRate;
    smoothMY += (targetMY - smoothMY) * posRate;

    if ((o.mode === "scroll" || o.mode === "both") && !o.externalScroll) updateScroll();

    if (!dragging) {
      const baseOmega = o.period > 0.1 ? (Math.PI * 2) / o.period : 0;
      const damp = Math.exp(-dt * 5.0);
      velY = velY * damp + baseOmega * (1 - damp);
      velX = velX * damp;
      rotY += velY * dt;
      rotX = Math.max(-1.45, Math.min(1.45, rotX + velX * dt));
    }

    let keepTicking = false;
    if (typeof o.onTick === "function") {
      keepTicking = !!o.onTick(dt, now);
    }

    render();

    const needsNextFrame =
      visible &&
      (o.period > 0.1 ||
        keepTicking ||
        dragging ||
        Math.abs(velY) > 1e-3 ||
        Math.abs(hoverAmt - targetHover) > 1e-3 ||
        o.mode === "scroll" ||
        o.mode === "both");
    if (needsNextFrame) raf = requestAnimationFrame(loop);
  }

  const kick = () => {
    if (!raf && visible && !destroyed) {
      prev = performance.now();
      raf = requestAnimationFrame(loop);
    }
  };

  const io =
    typeof IntersectionObserver === "function"
      ? new IntersectionObserver(entries => {
          visible = entries.some(e => e.isIntersecting);
          if (visible) kick();
        })
      : null;
  if (io) io.observe(root);

  const ro =
    typeof ResizeObserver === "function"
      ? new ResizeObserver(() => {
          fit();
          render();
          kick();
        })
      : null;
  if (ro) ro.observe(root);

  const onWinScroll = () => {
    if ((o.mode === "scroll" || o.mode === "both") && !o.externalScroll) {
      updateScroll();
      kick();
    }
  };
  window.addEventListener("scroll", onWinScroll, { passive: true });
  if (document.fonts && document.fonts.ready) {
    document.fonts.ready.then(() => {
      if (!destroyed) {
        fit();
        render();
      }
    });
  }

  let lastPX = 0, lastPY = 0;
  let lastTouchAt = -1e9;

  function updatePointerNorm(e) {
    const rect = inner.getBoundingClientRect();
    if (!rect.width || !rect.height) return;
    targetMX = Math.max(0, Math.min(1, (e.clientX - rect.left) / rect.width));
    targetMY = Math.max(0, Math.min(1, (e.clientY - rect.top) / rect.height));
  }

  root.addEventListener("pointerenter", e => {
    if (!o.interactive) return;
    if (e.pointerType !== "mouse" || performance.now() - lastTouchAt < 800) return;
    pointerInside = true;
    updatePointerNorm(e);
    smoothMX = targetMX;
    smoothMY = targetMY;
    kick();
  });

  root.addEventListener("pointermove", e => {
    if (!o.interactive) return;
    if (e.pointerType !== "mouse" || performance.now() - lastTouchAt < 800) return;
    pointerInside = true;
    updatePointerNorm(e);
    if (dragging) {
      const dx = e.clientX - lastPX;
      const dy = e.clientY - lastPY;
      lastPX = e.clientX;
      lastPY = e.clientY;
      const dRY = dx * 0.0075;
      const dRX = dy * 0.0075;
      rotY += dRY;
      rotX = Math.max(-1.45, Math.min(1.45, rotX + dRX));
      velY = dRY / 0.016;
      velX = dRX / 0.016;
    }
    kick();
  });

  root.addEventListener("pointerleave", () => {
    if (!dragging) pointerInside = false;
    kick();
  });

  root.addEventListener("pointerdown", e => {
    if (!o.interactive || !o.dragOrbit || (e.pointerType === "touch" && !o.touchOrbit) || (e.pointerType === "mouse" && e.button !== 0)) return;
    root.setPointerCapture(e.pointerId);
    dragging = true;
    pointerInside = true;
    updatePointerNorm(e);
    lastPX = e.clientX;
    lastPY = e.clientY;
    velY = 0;
    velX = 0;
    root.style.cursor = "grabbing";
    kick();
  });

  const endDrag = e => {
    if (!dragging) return;
    dragging = false;
    try {
      root.releasePointerCapture(e.pointerId);
    } catch (_) {}
    root.style.cursor = o.interactive ? "grab" : "default";
    const rect = root.getBoundingClientRect();
    if (
      e.pointerType === "touch" ||
      e.clientX < rect.left ||
      e.clientX > rect.right ||
      e.clientY < rect.top ||
      e.clientY > rect.bottom
    ) {
      pointerInside = false;
    }
    kick();
  };
  root.addEventListener("pointerup", endDrag);
  root.addEventListener("pointercancel", endDrag);

  root.addEventListener("dblclick", () => {
    if (!o.interactive) return;
    rotX = o.tiltX;
    zoom = o.zoom;
    panX = 0;
    panY = 0;
    kick();
  });

  function loadModel(bakedModel) {
    model = bakedModel;
    smoothSpanX0 = -1;
    smoothSpanX1 = -1;
    G = buildGeometry(model, shape);
    const label = `${model.id || "Protein"} — ${model.title || `${G.residues} residues`}`;
    root.setAttribute("aria-label", label);
    fit();
    render();
    kick();
    emit("load", { model, geometry: G });
    return G;
  }

  applyStyle();
  fit();

  return {
    loadModel,
    loadText(text, parseOpts = {}) {
      const parsed = parseStructureText(text, parseOpts);
      return loadModel(parsed);
    },
    async loadPDB(pdbId, parseOpts = {}) {
      const parsed = await fetchRCSB(pdbId, parseOpts);
      return loadModel(parsed);
    },
    async loadUniProt(acc, parseOpts = {}) {
      const parsed = await fetchAlphaFold(acc, parseOpts);
      return loadModel(parsed);
    },
    set(nextOpts = {}) {
      const prevDetail = o.detail;
      Object.assign(o, nextOpts);
      if (nextOpts.tiltX !== undefined && !dragging) rotX = o.tiltX;
      if (nextOpts.zoom !== undefined) zoom = o.zoom;
      if ((nextOpts.mode === "scroll" || nextOpts.mode === "both") && !o.externalScroll) updateScroll();
      applyStyle();
      if (nextOpts.detail !== undefined && nextOpts.detail !== prevDetail) fit();
      render();
      kick();
    },
    setShape(nextShape = {}) {
      shape = { ...shape, ...nextShape };
      if (model) G = buildGeometry(model, shape);
      render();
      kick();
    },
    setModelFast(bakedModel, shapeOverride = null, fixedR = null) {
      model = bakedModel;
      G = buildGeometry(model, shapeOverride ? { ...shape, ...shapeOverride } : shape);
      if (fixedR) G.R = fixedR;
    },
    setSkipDefaultInk(flag) {
      o.skipDefaultInk = !!flag;
    },
    render,
    nudgeSpin(dVelY = 0, dVelX = 0) {
      velY += dVelY;
      velX += dVelX;
      kick();
    },
    kick,
    setRotation(ry, rx = rotX) {
      rotY = ry;
      rotX = rx;
      render();
    },
    setPointer(nx, ny, active = true) {
      pointerInside = active;
      targetMX = smoothMX = nx;
      targetMY = smoothMY = ny;
      hoverAmt = active ? 1 : 0;
      render();
    },
    setTouchLens(clientX, clientY, active = true) {
      lastTouchAt = performance.now();
      if (!active) {
        pointerInside = false;
        kick();
        return;
      }
      const rect = inner.getBoundingClientRect();
      if (rect.width && rect.height) {
        const rx = (clientX - rect.left) / rect.width;
        const ry = (clientY - 24 - rect.top) / rect.height;
        targetMX = Math.max(0.18, Math.min(0.82, rx));
        targetMY = Math.max(0.18, Math.min(0.82, ry));
        if (!pointerInside && hoverAmt < 0.05) {
          smoothMX = targetMX;
          smoothMY = targetMY;
        }
      }
      pointerInside = true;
      kick();
    },
    setScroll(amt) {
      scrollAmt = Math.max(0, Math.min(1, amt));
      render();
    },
    setScrollFast(amt) {
      scrollAmt = Math.max(0, Math.min(1, amt));
    },
    getText() {
      return lastText;
    },
    getOptions() {
      return { ...o, shape: { ...shape } };
    },
    getInfo() {
      return model ? { id: model.id, title: model.title, ...G } : null;
    },
    on(ev, fn) {
      (listeners[ev] ||= []).push(fn);
    },
    destroy() {
      destroyed = true;
      if (raf) cancelAnimationFrame(raf);
      if (io) io.disconnect();
      if (ro) ro.disconnect();
      window.removeEventListener("scroll", onWinScroll);
      root.remove();
    },
  };
}
