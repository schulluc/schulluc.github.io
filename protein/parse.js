import { dist, sub, dot, norm } from "./vec.js";

// Secondary structure fallback from Cα distances alone (for structures without HELIX/SHEET).
export function assignSS(P) {
  const n = P.length;
  const ss = new Array(n).fill("C");
  if (n < 6) return ss.join("");
  const d = (i, j) => dist(P[i], P[j]);

  const hw = i =>
    i + 4 < n &&
    Math.abs(d(i, i + 2) - 5.45) < 0.65 &&
    Math.abs(d(i, i + 3) - 5.20) < 0.75 &&
    Math.abs(d(i, i + 4) - 6.20) < 0.85;

  for (let i = 0; i + 5 < n; i++) {
    if (hw(i) && hw(i + 1)) {
      for (let k = i; k <= i + 5; k++) ss[k] = "H";
    }
  }

  const ext = new Array(n).fill(false);
  for (let i = 0; i + 2 < n; i++) {
    if (d(i, i + 2) > 6.1 && (i + 3 >= n || d(i, i + 3) > 8.4)) {
      ext[i] = ext[i + 1] = ext[i + 2] = true;
    }
  }

  const cell = 5.5;
  const grid = new Map();
  const keyOf = p => `${Math.floor(p[0] / cell)},${Math.floor(p[1] / cell)},${Math.floor(p[2] / cell)}`;
  P.forEach((p, i) => {
    const k = keyOf(p);
    if (!grid.has(k)) grid.set(k, []);
    grid.get(k).push(i);
  });

  const ok = (a, b) => a >= 0 && b >= 0 && a < n && b < n && Math.abs(a - b) > 2 && d(a, b) < 5.5;
  for (let i = 0; i < n; i++) {
    if (!ext[i] || ss[i] === "H") continue;
    const gx = Math.floor(P[i][0] / cell);
    const gy = Math.floor(P[i][1] / cell);
    const gz = Math.floor(P[i][2] / cell);
    let paired = false;
    for (let ax = -1; ax <= 1 && !paired; ax++) {
      for (let ay = -1; ay <= 1 && !paired; ay++) {
        for (let az = -1; az <= 1 && !paired; az++) {
          const bucket = grid.get(`${gx + ax},${gy + ay},${gz + az}`);
          if (!bucket) continue;
          for (const j of bucket) {
            if (Math.abs(i - j) < 3 || !ext[j] || ss[j] === "H" || d(i, j) > 5.5) continue;
            if (ok(i + 1, j + 1) || ok(i - 1, j - 1) || ok(i + 1, j - 1) || ok(i - 1, j + 1)) {
              paired = true;
              break;
            }
          }
        }
      }
    }
    if (paired) ss[i] = "E";
  }

  for (let i = 1; i < n - 1; i++) {
    if (ss[i] === "C" && ss[i - 1] === "E" && ss[i + 1] === "E" && ext[i]) ss[i] = "E";
  }

  const clean = (type, minLen) => {
    for (let i = 0; i < n; ) {
      if (ss[i] !== type) {
        i++;
        continue;
      }
      let j = i;
      while (j < n && ss[j] === type) j++;
      if (j - i < minLen) {
        for (let k = i; k < j; k++) ss[k] = "C";
      }
      i = j;
    }
  };
  clean("H", 4);
  clean("E", 3);
  return ss.join("");
}

function alignStrandAxisToY(chains) {
  const vecs = [];
  for (const ch of chains) {
    const ss = ch.ss;
    const n = ss.length;
    for (let i = 0; i < n; ) {
      if (ss[i] !== "E") {
        i++;
        continue;
      }
      let j = i;
      while (j < n && ss[j] === "E") j++;
      if (j - i >= 3) {
        const i0 = i * 3;
        const i1 = (j - 1) * 3;
        const v = [ch.ca[i1] - ch.ca[i0], ch.ca[i1 + 1] - ch.ca[i0 + 1], ch.ca[i1 + 2] - ch.ca[i0 + 2]];
        const l = Math.hypot(v[0], v[1], v[2]);
        if (l > 1) vecs.push([v[0] / l, v[1] / l, v[2] / l]);
      }
      i = j;
    }
  }
  if (vecs.length < 3) return;
  const ref = vecs[0];
  let sx = 0, sy = 0, sz = 0;
  for (const v of vecs) {
    const s = dot(v, ref) < 0 ? -1 : 1;
    sx += v[0] * s;
    sy += v[1] * s;
    sz += v[2] * s;
  }
  const ln = Math.hypot(sx, sy, sz);
  if (ln < 0.5) return;
  const u = [sx / ln, sy / ln, sz / ln];
  // cross(u, [0, 1, 0]) = [-u[2], 0, u[0]]
  let vx = -u[2], vy = 0, vz = u[0];
  const s = Math.hypot(vx, vy, vz);
  const c = u[1];
  if (s < 1e-6) return;
  vx /= s; vz /= s;
  const k = 1 - c;
  const M = [
    [c + vx * vx * k, vx * vy * k - vz * s, vx * vz * k + vy * s],
    [vy * vx * k + vz * s, c + vy * vy * k, vy * vz * k - vx * s],
    [vz * vx * k - vy * s, vz * vy * k + vx * s, c + vz * vz * k],
  ];
  const rotArr = arr => {
    for (let i = 0; i < arr.length; i += 3) {
      const x = arr[i], y = arr[i + 1], z = arr[i + 2];
      arr[i] = M[0][0] * x + M[0][1] * y + M[0][2] * z;
      arr[i + 1] = M[1][0] * x + M[1][1] * y + M[1][2] * z;
      arr[i + 2] = M[2][0] * x + M[2][1] * y + M[2][2] * z;
    }
  };
  for (const ch of chains) {
    rotArr(ch.ca);
    rotArr(ch.o);
  }
}

function finalizeChains(id, title, rawChains, align = true) {
  let tot = 0, cx = 0, cy = 0, cz = 0;
  for (const ch of rawChains) {
    for (const r of ch.rows) {
      cx += r.CA[0]; cy += r.CA[1]; cz += r.CA[2];
      tot++;
    }
  }
  if (!tot) throw new Error("No Cα coordinates found in structure.");
  cx /= tot; cy /= tot; cz /= tot;

  const chains = rawChains.map(ch => {
    const ca = new Array(ch.rows.length * 3);
    const o = new Array(ch.rows.length * 3);
    ch.rows.forEach((r, idx) => {
      const x = r.CA[0] - cx, y = r.CA[1] - cy, z = r.CA[2] - cz;
      ca[idx * 3] = x; ca[idx * 3 + 1] = y; ca[idx * 3 + 2] = z;
      if (r.O) {
        o[idx * 3] = r.O[0] - cx; o[idx * 3 + 1] = r.O[1] - cy; o[idx * 3 + 2] = r.O[2] - cz;
      } else {
        o[idx * 3] = x + 1.2; o[idx * 3 + 1] = y; o[idx * 3 + 2] = z;
      }
    });
    return { id: ch.id, ss: ch.ss, ca, o };
  });

  if (align && chains.length <= 4) alignStrandAxisToY(chains);
  return { id: id || "PROTEIN", title: title || "", residues: tot, chains };
}

export function parsePDB(text, { chains: wantedChains = "all", align = true } = {}) {
  const chainsMap = new Map();
  const chainOrder = [];
  const helix = [];
  const sheet = [];
  const titleParts = [];
  let pdbId = null;

  const wantSet = Array.isArray(wantedChains) ? new Set(wantedChains) : null;
  const wantAll = wantedChains === "all";

  for (const l of text.split(/\r?\n/)) {
    const rec = l.slice(0, 6);
    if (rec === "ENDMDL") break;
    if (rec === "HEADER" && l.length >= 66) pdbId = l.slice(62, 66).trim() || pdbId;
    if (rec === "TITLE ") titleParts.push(l.slice(10, 70).trim());
    if (rec === "HELIX ") {
      helix.push([l[19], parseInt(l.slice(21, 25), 10), parseInt(l.slice(33, 37), 10)]);
      continue;
    }
    if (rec === "SHEET ") {
      sheet.push([l[21], parseInt(l.slice(22, 26), 10), parseInt(l.slice(33, 37), 10)]);
      continue;
    }
    if (rec !== "ATOM  ") continue;
    const name = l.slice(12, 16).trim();
    if (name !== "CA" && name !== "O") continue;
    const alt = l[16];
    if (alt !== " " && alt !== "A" && alt !== "1") continue;
    const ch = l[21];
    if (wantSet && !wantSet.has(ch)) continue;
    if (!chainsMap.has(ch)) {
      if (!wantSet && !wantAll && chainOrder.length >= 1) continue;
      chainsMap.set(ch, { order: [], map: new Map() });
      chainOrder.push(ch);
    }
    const bucket = chainsMap.get(ch);
    const key = l.slice(22, 27);
    let r = bucket.map.get(key);
    if (!r) {
      r = { num: parseInt(l.slice(22, 26), 10) };
      bucket.map.set(key, r);
      bucket.order.push(key);
    }
    r[name] = [
      parseFloat(l.slice(30, 38)),
      parseFloat(l.slice(38, 46)),
      parseFloat(l.slice(46, 54)),
    ];
  }

  const rawChains = [];
  for (const ch of chainOrder) {
    const bucket = chainsMap.get(ch);
    const rows = bucket.order.map(k => bucket.map.get(k)).filter(r => r.CA);
    if (rows.length < 5) continue;
    const ssMap = new Map();
    for (const [c, a, b] of helix) if (c === ch) for (let i = a; i <= b; i++) ssMap.set(i, "H");
    for (const [c, a, b] of sheet) if (c === ch) for (let i = a; i <= b; i++) ssMap.set(i, "E");
    const ss = ssMap.size
      ? rows.map(r => ssMap.get(r.num) || "C").join("")
      : assignSS(rows.map(r => r.CA));
    rawChains.push({ id: ch, rows, ss });
  }

  const title = titleParts.join(" ").replace(/\s+/g, " ").trim();
  return finalizeChains(pdbId, title, rawChains, align);
}

function cifTokens(line) {
  const out = [];
  const re = /'([^']*)'|"([^"]*)"|(\S+)/g;
  let m;
  while ((m = re.exec(line))) out.push(m[1] ?? m[2] ?? m[3]);
  return out;
}

function parseCifLoop(lines, prefix) {
  let idx = lines.findIndex(l => l.startsWith(prefix));
  if (idx < 0) return { cols: [], rows: [] };
  const cols = [];
  while (idx < lines.length && lines[idx].startsWith(prefix)) {
    cols.push(lines[idx].trim().slice(prefix.length));
    idx++;
  }
  const rows = [];
  for (let i = idx; i < lines.length; i++) {
    const s = lines[i].trim();
    if (!s) continue;
    if (s[0] === "#" || s[0] === "_" || s.startsWith("loop_")) break;
    const t = cifTokens(s);
    if (t.length >= cols.length) rows.push(t);
  }
  return { cols, rows };
}

export function parseCIF(text, { chains: wantedChains = "all", align = true } = {}) {
  const lines = text.split(/\r?\n/);
  let pdbId = "CIF";
  if (lines[0] && lines[0].startsWith("data_")) {
    pdbId = lines[0].slice(5).trim().split("-")[0] || "CIF";
  }

  const helix = [];
  const hLoop = parseCifLoop(lines, "_struct_conf.");
  if (hLoop.cols.length) {
    const icT = hLoop.cols.indexOf("conf_type_id");
    const icC = hLoop.cols.includes("beg_auth_asym_id") ? hLoop.cols.indexOf("beg_auth_asym_id") : hLoop.cols.indexOf("beg_label_asym_id");
    const icS = hLoop.cols.includes("beg_auth_seq_id") ? hLoop.cols.indexOf("beg_auth_seq_id") : hLoop.cols.indexOf("beg_label_seq_id");
    const icE = hLoop.cols.includes("end_auth_seq_id") ? hLoop.cols.indexOf("end_auth_seq_id") : hLoop.cols.indexOf("end_label_seq_id");
    for (const r of hLoop.rows) {
      if (icT >= 0 && !r[icT].startsWith("HELX")) continue;
      const a = parseInt(r[icS], 10), b = parseInt(r[icE], 10);
      if (!Number.isNaN(a) && !Number.isNaN(b)) helix.push([r[icC], a, b]);
    }
  }

  const sheet = [];
  const sLoop = parseCifLoop(lines, "_struct_sheet_range.");
  if (sLoop.cols.length) {
    const icC = sLoop.cols.includes("beg_auth_asym_id") ? sLoop.cols.indexOf("beg_auth_asym_id") : sLoop.cols.indexOf("beg_label_asym_id");
    const icS = sLoop.cols.includes("beg_auth_seq_id") ? sLoop.cols.indexOf("beg_auth_seq_id") : sLoop.cols.indexOf("beg_label_seq_id");
    const icE = sLoop.cols.includes("end_auth_seq_id") ? sLoop.cols.indexOf("end_auth_seq_id") : sLoop.cols.indexOf("end_label_seq_id");
    for (const r of sLoop.rows) {
      const a = parseInt(r[icS], 10), b = parseInt(r[icE], 10);
      if (!Number.isNaN(a) && !Number.isNaN(b)) sheet.push([r[icC], a, b]);
    }
  }

  let i = lines.findIndex(l => l.startsWith("_atom_site."));
  if (i < 0) throw new Error("No _atom_site block in mmCIF file.");
  const cols = [];
  while (i < lines.length && lines[i].startsWith("_atom_site.")) {
    cols.push(lines[i].trim().slice(11));
    i++;
  }
  const c = n => cols.indexOf(n);
  const pick = (a, b) => (c(a) >= 0 ? c(a) : c(b));
  const ix = {
    g: c("group_PDB"),
    atom: c("label_atom_id"),
    alt: c("label_alt_id"),
    chain: pick("auth_asym_id", "label_asym_id"),
    seq: pick("auth_seq_id", "label_seq_id"),
    ins: c("pdbx_PDB_ins_code"),
    x: c("Cartn_x"),
    y: c("Cartn_y"),
    z: c("Cartn_z"),
    model: c("pdbx_PDB_model_num"),
  };

  const wantSet = Array.isArray(wantedChains) ? new Set(wantedChains) : null;
  const wantAll = wantedChains === "all";
  const chainsMap = new Map();
  const chainOrder = [];
  let firstModel = null;

  for (; i < lines.length; i++) {
    const l = lines[i].trim();
    if (!l) continue;
    if (l[0] === "#" || l[0] === "_" || l.startsWith("loop_")) break;
    const t = cifTokens(l);
    if (t.length < cols.length) continue;
    if (ix.g >= 0 && t[ix.g] !== "ATOM") continue;
    if (ix.model >= 0) {
      if (firstModel === null) firstModel = t[ix.model];
      else if (t[ix.model] !== firstModel) break;
    }
    const name = t[ix.atom];
    if (name !== "CA" && name !== "O") continue;
    if (ix.alt >= 0 && !".?A1".includes(t[ix.alt])) continue;
    const ch = t[ix.chain];
    if (wantSet && !wantSet.has(ch)) continue;
    if (!chainsMap.has(ch)) {
      if (!wantSet && !wantAll && chainOrder.length >= 1) continue;
      chainsMap.set(ch, { order: [], map: new Map() });
      chainOrder.push(ch);
    }
    const bucket = chainsMap.get(ch);
    const ins = ix.ins >= 0 && !".?".includes(t[ix.ins]) ? t[ix.ins] : "";
    const key = t[ix.seq] + ins;
    let r = bucket.map.get(key);
    if (!r) {
      r = { num: parseInt(t[ix.seq], 10) };
      bucket.map.set(key, r);
      bucket.order.push(key);
    }
    r[name] = [parseFloat(t[ix.x]), parseFloat(t[ix.y]), parseFloat(t[ix.z])];
  }

  const rawChains = [];
  for (const ch of chainOrder) {
    const bucket = chainsMap.get(ch);
    const rows = bucket.order.map(k => bucket.map.get(k)).filter(r => r.CA);
    if (rows.length < 5) continue;
    const baseCh = ch.split("-")[0];
    const ssMap = new Map();
    for (const [cId, a, b] of helix) if (cId === ch || cId === baseCh) for (let k = a; k <= b; k++) ssMap.set(k, "H");
    for (const [cId, a, b] of sheet) if (cId === ch || cId === baseCh) for (let k = a; k <= b; k++) ssMap.set(k, "E");
    const ss = ssMap.size
      ? rows.map(r => ssMap.get(r.num) || "C").join("")
      : assignSS(rows.map(r => r.CA));
    rawChains.push({ id: ch, rows, ss });
  }
  return finalizeChains(pdbId, "", rawChains, align);
}

export function parseStructureText(text, opts = {}) {
  if (/^data_/m.test(text.slice(0, 4000)) || text.includes("_atom_site.")) {
    return parseCIF(text, opts);
  }
  return parsePDB(text, opts);
}

export async function fetchRCSB(pdbId, opts = {}) {
  const code = pdbId.trim().toUpperCase();
  for (const ext of ["pdb", "cif"]) {
    const res = await fetch(`https://files.rcsb.org/download/${code}.${ext}`);
    if (res.ok) {
      const text = await res.text();
      const parsed = parseStructureText(text, opts);
      parsed.id = code;
      return parsed;
    }
  }
  throw new Error(`Could not fetch ${code} from RCSB PDB.`);
}

export async function fetchAlphaFold(uniprotId, opts = {}) {
  const acc = uniprotId.trim().toUpperCase();
  const metaRes = await fetch(`https://alphafold.ebi.ac.uk/api/prediction/${acc}`);
  if (!metaRes.ok) throw new Error(`AlphaFold DB entry not found for ${acc}`);
  const meta = await metaRes.json();
  const entry = Array.isArray(meta) ? meta[0] : meta;
  const pdbUrl = entry?.pdbUrl || entry?.cifUrl;
  if (!pdbUrl) throw new Error(`No structure URL in AlphaFold DB for ${acc}`);
  const sRes = await fetch(pdbUrl);
  if (!sRes.ok) throw new Error(`Failed downloading AlphaFold structure for ${acc}`);
  const text = await sRes.text();
  const parsed = parseStructureText(text, opts);
  parsed.id = acc;
  parsed.title = entry?.uniprotDescription || acc;
  return parsed;
}
