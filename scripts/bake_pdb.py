#!/usr/bin/env python3
"""
bake_pdb.py — Convert any PDB/mmCIF file (or RCSB 4-letter PDB code / biological assembly)
into a compact, zero-dependency ES module for the Richardson Ink + ASCII protein renderer.

Examples:
  python3 scripts/bake_pdb.py 1TIM --chains A -o protein/structures/1tim-a.js
  python3 scripts/bake_pdb.py 7QSW --assembly 1 --chains all --decimals 1 -o protein/structures/7qsw-l8s8.js
  python3 scripts/bake_pdb.py path/to/my_design.pdb -o protein/structures/design1.js
"""

import argparse
import gzip
import json
import math
import os
import re
import sys
import urllib.request


def fetch_rcsb(pdb_id: str, assembly: int = 0) -> str:
    code = pdb_id.strip().upper()
    if assembly > 0:
        url = f"https://files.rcsb.org/download/{code}-assembly{assembly}.cif.gz"
        req = urllib.request.Request(url, headers={"User-Agent": "protein-baker/1.0"})
        with urllib.request.urlopen(req, timeout=30) as resp:
            return gzip.decompress(resp.read()).decode("utf-8", errors="replace")

    for ext in ("pdb", "cif"):
        url = f"https://files.rcsb.org/download/{code}.{ext}"
        try:
            req = urllib.request.Request(url, headers={"User-Agent": "protein-baker/1.0"})
            with urllib.request.urlopen(req, timeout=25) as resp:
                return resp.read().decode("utf-8", errors="replace")
        except Exception:
            continue
    raise RuntimeError(f"Could not download {code} (.pdb or .cif) from RCSB.")


def assign_ss_ca(ca_pts):
    """Fallback C-alpha secondary structure assignment (P-SEA-like) when HELIX/SHEET are absent."""
    n = len(ca_pts)
    ss = ["C"] * n
    if n < 6:
        return "".join(ss)

    def d(i, j):
        return math.dist(ca_pts[i], ca_pts[j])

    def hw(i):
        return (
            i + 4 < n
            and abs(d(i, i + 2) - 5.45) < 0.65
            and abs(d(i, i + 3) - 5.20) < 0.75
            and abs(d(i, i + 4) - 6.20) < 0.85
        )

    for i in range(n - 5):
        if hw(i) and hw(i + 1):
            for k in range(i, i + 6):
                ss[k] = "H"

    ext = [False] * n
    for i in range(n - 2):
        if d(i, i + 2) > 6.1 and (i + 3 >= n or d(i, i + 3) > 8.4):
            ext[i] = ext[i + 1] = ext[i + 2] = True

    def ok(a, b):
        return 0 <= a < n and 0 <= b < n and abs(a - b) > 2 and d(a, b) < 5.5

    for i in range(n):
        if not ext[i] or ss[i] == "H":
            continue
        for j in range(n):
            if abs(i - j) < 3 or not ext[j] or ss[j] == "H" or d(i, j) > 5.5:
                continue
            if ok(i + 1, j + 1) or ok(i - 1, j - 1) or ok(i + 1, j - 1) or ok(i - 1, j + 1):
                ss[i] = "E"
                break

    for i in range(1, n - 1):
        if ss[i] == "C" and ss[i - 1] == "E" and ss[i + 1] == "E" and ext[i]:
            ss[i] = "E"

    for kind, min_len in (("H", 4), ("E", 3)):
        i = 0
        while i < n:
            if ss[i] != kind:
                i += 1
                continue
            j = i
            while j < n and ss[j] == kind:
                j += 1
            if j - i < min_len:
                for k in range(i, j):
                    ss[k] = "C"
            i = j

    return "".join(ss)


def parse_pdb(text: str, wanted_chains=None, want_all=False):
    helix = []
    sheet = []
    chains_map = {}
    chain_order = []
    title_parts = []
    pdb_id = None

    for line in text.splitlines():
        rec = line[:6]
        if rec == "HEADER" and len(line) >= 66:
            pdb_id = line[62:66].strip() or None
        elif rec == "TITLE ":
            title_parts.append(line[10:70].strip())
        elif rec == "ENDMDL":
            break
        elif rec == "HELIX ":
            try:
                helix.append((line[19], int(line[21:25]), int(line[33:37])))
            except ValueError:
                pass
        elif rec == "SHEET ":
            try:
                sheet.append((line[21], int(line[22:26]), int(line[33:37])))
            except ValueError:
                pass
        elif rec == "ATOM  ":
            name = line[12:16].strip()
            if name not in ("CA", "O"):
                continue
            alt = line[16]
            if alt not in (" ", "A", "1"):
                continue
            ch = line[21]
            if wanted_chains and ch not in wanted_chains:
                continue
            if ch not in chains_map:
                if not want_all and wanted_chains is None and len(chain_order) >= 1:
                    continue
                chains_map[ch] = ([], {})
                chain_order.append(ch)
            order, res_dict = chains_map[ch]
            key = line[22:27]
            if key not in res_dict:
                res_dict[key] = {"num": int(line[22:26])}
                order.append(key)
            res_dict[key][name] = [
                float(line[30:38]),
                float(line[38:46]),
                float(line[46:54]),
            ]

    out_chains = []
    for ch in chain_order:
        order, res_dict = chains_map[ch]
        rows = [res_dict[k] for k in order if "CA" in res_dict[k]]
        if len(rows) < 5:
            continue
        ss_by_num = {}
        for c, a, b in helix:
            if c == ch:
                for r in range(a, b + 1):
                    ss_by_num[r] = "H"
        for c, a, b in sheet:
            if c == ch:
                for r in range(a, b + 1):
                    ss_by_num[r] = "E"
        ss = "".join(ss_by_num.get(r["num"], "C") for r in rows) if ss_by_num else assign_ss_ca([r["CA"] for r in rows])
        out_chains.append({"id": ch, "rows": rows, "ss": ss})

    title = re.sub(r"\s+", " ", " ".join(title_parts).strip()).capitalize() if title_parts else ""
    return {"id": pdb_id, "title": title, "chains": out_chains}


def cif_tokens(line: str):
    return [
        m.group(1) or m.group(2) or m.group(3)
        for m in re.finditer(r"'([^']*)'|\"([^\"]*)\"|(\S+)", line)
    ]


def parse_cif_loop(lines, prefix):
    idx = next((i for i, l in enumerate(lines) if l.startswith(prefix)), -1)
    if idx < 0:
        return [], []
    cols = []
    while idx < len(lines) and lines[idx].startswith(prefix):
        cols.append(lines[idx].strip()[len(prefix):])
        idx += 1
    rows = []
    for i in range(idx, len(lines)):
        s = lines[i].strip()
        if not s:
            continue
        if s[0] in "#_" or s.startswith("loop_"):
            break
        t = cif_tokens(s)
        if len(t) >= len(cols):
            rows.append(t)
    return cols, rows


def parse_cif(text: str, wanted_chains=None, want_all=False, include_all_models=False):
    lines = text.splitlines()
    pdb_id = None
    if lines and lines[0].startswith("data_"):
        pdb_id = lines[0][5:].strip().split("-")[0] or None

    # Secondary structure loops (_struct_conf and _struct_sheet_range)
    helix = []
    h_cols, h_rows = parse_cif_loop(lines, "_struct_conf.")
    if h_cols:
        ic_type = h_cols.index("conf_type_id") if "conf_type_id" in h_cols else -1
        ic_ch = h_cols.index("beg_auth_asym_id") if "beg_auth_asym_id" in h_cols else h_cols.index("beg_label_asym_id")
        ic_s = h_cols.index("beg_auth_seq_id") if "beg_auth_seq_id" in h_cols else h_cols.index("beg_label_seq_id")
        ic_e = h_cols.index("end_auth_seq_id") if "end_auth_seq_id" in h_cols else h_cols.index("end_label_seq_id")
        for r in h_rows:
            if ic_type >= 0 and not r[ic_type].startswith("HELX"):
                continue
            try:
                helix.append((r[ic_ch], int(r[ic_s]), int(r[ic_e])))
            except ValueError:
                pass

    sheet = []
    s_cols, s_rows = parse_cif_loop(lines, "_struct_sheet_range.")
    if s_cols:
        ic_ch = s_cols.index("beg_auth_asym_id") if "beg_auth_asym_id" in s_cols else s_cols.index("beg_label_asym_id")
        ic_s = s_cols.index("beg_auth_seq_id") if "beg_auth_seq_id" in s_cols else s_cols.index("beg_label_seq_id")
        ic_e = s_cols.index("end_auth_seq_id") if "end_auth_seq_id" in s_cols else s_cols.index("end_label_seq_id")
        for r in s_rows:
            try:
                sheet.append((r[ic_ch], int(r[ic_s]), int(r[ic_e])))
            except ValueError:
                pass

    cols, atom_rows = parse_cif_loop(lines, "_atom_site.")
    if not cols:
        raise RuntimeError("No _atom_site block in mmCIF file.")

    def col(a, b=None):
        if a in cols:
            return cols.index(a)
        if b and b in cols:
            return cols.index(b)
        return -1

    ix_g = col("group_PDB")
    ix_atom = col("label_atom_id")
    ix_alt = col("label_alt_id")
    ix_ch = col("auth_asym_id", "label_asym_id")
    ix_seq = col("auth_seq_id", "label_seq_id")
    ix_ins = col("pdbx_PDB_ins_code")
    ix_x, ix_y, ix_z = col("Cartn_x"), col("Cartn_y"), col("Cartn_z")
    ix_mod = col("pdbx_PDB_model_num")

    chains_map = {}
    chain_order = []
    first_model = None

    for t in atom_rows:
        if ix_g >= 0 and t[ix_g] != "ATOM":
            continue
        if ix_mod >= 0 and not include_all_models:
            if first_model is None:
                first_model = t[ix_mod]
            elif t[ix_mod] != first_model:
                break
        name = t[ix_atom]
        if name not in ("CA", "O"):
            continue
        if ix_alt >= 0 and t[ix_alt] not in (".", "?", "A", "1"):
            continue
        ch_raw = t[ix_ch]
        ch = f"{ch_raw}_m{t[ix_mod]}" if (include_all_models and ix_mod >= 0 and t[ix_mod] != "1") else ch_raw
        if wanted_chains and ch not in wanted_chains and ch_raw not in wanted_chains:
            continue
        if ch not in chains_map:
            if not want_all and wanted_chains is None and len(chain_order) >= 1:
                continue
            chains_map[ch] = ([], {}, ch_raw)
            chain_order.append(ch)
        order, res_dict, _ = chains_map[ch]
        ins = t[ix_ins] if (ix_ins >= 0 and t[ix_ins] not in (".", "?")) else ""
        key = f"{t[ix_seq]}{ins}"
        if key not in res_dict:
            res_dict[key] = {"num": int(t[ix_seq])}
            order.append(key)
        res_dict[key][name] = [float(t[ix_x]), float(t[ix_y]), float(t[ix_z])]

    out_chains = []
    for ch in chain_order:
        order, res_dict, ch_raw = chains_map[ch]
        rows = [res_dict[k] for k in order if "CA" in res_dict[k]]
        if len(rows) < 5:
            continue
        base_ch = ch_raw.split("-")[0]
        ss_by_num = {}
        for c, a, b in helix:
            if c in (ch_raw, base_ch):
                for r in range(a, b + 1):
                    ss_by_num[r] = "H"
        for c, a, b in sheet:
            if c in (ch_raw, base_ch):
                for r in range(a, b + 1):
                    ss_by_num[r] = "E"
        ss = "".join(ss_by_num.get(r["num"], "C") for r in rows) if ss_by_num else assign_ss_ca([r["CA"] for r in rows])
        out_chains.append({"id": ch, "rows": rows, "ss": ss})

    return {"id": pdb_id, "title": "", "chains": out_chains}


def rot_matrix_align_to_y(axis):
    ax, ay, az = axis
    ln = math.hypot(ax, ay, az) or 1.0
    u = [ax / ln, ay / ln, az / ln]
    vx, vy, vz = -u[2], 0.0, u[0]
    s = math.hypot(vx, vy, vz)
    c = u[1]
    if s < 1e-6:
        return [[1, 0, 0], [0, 1 if c > 0 else -1, 0], [0, 0, 1 if c > 0 else -1]]
    vx /= s
    vz /= s
    k = 1.0 - c
    return [
        [c + vx * vx * k, vx * vy * k - vz * s, vx * vz * k + vy * s],
        [vy * vx * k + vz * s, c + vy * vy * k, vy * vz * k - vx * s],
        [vz * vx * k - vy * s, vz * vy * k + vx * s, c + vz * vz * k],
    ]


def apply_rot(M, p):
    return [
        M[0][0] * p[0] + M[0][1] * p[1] + M[0][2] * p[2],
        M[1][0] * p[0] + M[1][1] * p[1] + M[1][2] * p[2],
        M[2][0] * p[0] + M[2][1] * p[1] + M[2][2] * p[2],
    ]


def compute_strand_axis(chains):
    vecs = []
    for ch in chains:
        ss = ch["ss"]
        rows = ch["rows"]
        n = len(ss)
        i = 0
        while i < n:
            if ss[i] != "E":
                i += 1
                continue
            j = i
            while j < n and ss[j] == "E":
                j += 1
            if j - i >= 3:
                p0 = rows[i]["CA"]
                p1 = rows[j - 1]["CA"]
                dx, dy, dz = p1[0] - p0[0], p1[1] - p0[1], p1[2] - p0[2]
                ln = math.hypot(dx, dy, dz)
                if ln > 1.0:
                    vecs.append([dx / ln, dy / ln, dz / ln])
            i = j
    if len(vecs) < 3:
        return None
    ref = vecs[0]
    sx = sy = sz = 0.0
    for v in vecs:
        sign = -1.0 if (v[0] * ref[0] + v[1] * ref[1] + v[2] * ref[2]) < 0 else 1.0
        sx += v[0] * sign
        sy += v[1] * sign
        sz += v[2] * sign
    ln = math.hypot(sx, sy, sz)
    return [sx / ln, sy / ln, sz / ln] if ln > 0.5 else None


def bake(parsed, align_mode="auto", label=None, decimals=2):
    chains = parsed["chains"]
    if not chains:
        raise RuntimeError("No valid protein chains found.")

    all_ca = [r["CA"] for ch in chains for r in ch["rows"]]
    n_tot = len(all_ca)
    cx = sum(p[0] for p in all_ca) / n_tot
    cy = sum(p[1] for p in all_ca) / n_tot
    cz = sum(p[2] for p in all_ca) / n_tot

    for ch in chains:
        for r in ch["rows"]:
            r["CA"] = [r["CA"][0] - cx, r["CA"][1] - cy, r["CA"][2] - cz]
            if "O" in r:
                r["O"] = [r["O"][0] - cx, r["O"][1] - cy, r["O"][2] - cz]

    axis = None
    if align_mode.startswith("axis:"):
        axis = [float(x) for x in align_mode[5:].split(",")]
    elif align_mode in ("auto", "strands") and len(chains) <= 4:
        axis = compute_strand_axis(chains)

    if axis is not None:
        M = rot_matrix_align_to_y(axis)
        for ch in chains:
            for r in ch["rows"]:
                r["CA"] = apply_rot(M, r["CA"])
                if "O" in r:
                    r["O"] = apply_rot(M, r["O"])

    compact_chains = []
    for ch in chains:
        ca_flat = []
        o_flat = []
        for r in ch["rows"]:
            ca = r["CA"]
            ca_flat.extend([round(ca[0], decimals), round(ca[1], decimals), round(ca[2], decimals)])
            if "O" in r:
                o = r["O"]
                o_flat.extend([round(o[0], decimals), round(o[1], decimals), round(o[2], decimals)])
            else:
                o_flat.extend([round(ca[0] + 1.2, decimals), round(ca[1], decimals), round(ca[2], decimals)])
        compact_chains.append({"id": ch["id"], "ss": ch["ss"], "ca": ca_flat, "o": o_flat})

    return {
        "id": label or parsed["id"] or "PROTEIN",
        "title": parsed["title"],
        "residues": n_tot,
        "chains": compact_chains,
    }


def main():
    ap = argparse.ArgumentParser(description="Bake PDB/mmCIF into compact ES module.")
    ap.add_argument("source", help="Path to .pdb/.cif(.gz) file OR 4-letter RCSB PDB ID")
    ap.add_argument("--assembly", type=int, default=0, help="RCSB biological assembly number (e.g. 1)")
    ap.add_argument("--chains", default=None, help="Comma-separated chain IDs or 'all'")
    ap.add_argument("--align", default="auto", help="'auto', 'strands', 'none', or 'axis:x,y,z'")
    ap.add_argument("--decimals", type=int, default=2, help="Coordinate decimal places (default 2, use 1 for huge complexes)")
    ap.add_argument("--title", default=None, help="Override structure title")
    ap.add_argument("-o", "--output", required=True, help="Output .js or .json path")
    args = ap.parse_args()

    if os.path.isfile(args.source):
        if args.source.endswith(".gz"):
            with gzip.open(args.source, "rt", encoding="utf-8", errors="replace") as f:
                text = f.read()
        else:
            with open(args.source, "r", encoding="utf-8", errors="replace") as f:
                text = f.read()
        default_id = os.path.basename(args.source).split(".")[0].split("-")[0].upper()
    elif re.fullmatch(r"[A-Za-z0-9]{4}", args.source.strip()):
        default_id = args.source.strip().upper()
        text = fetch_rcsb(default_id, assembly=args.assembly)
    else:
        sys.exit(f"Source '{args.source}' is neither an existing file nor a 4-char PDB code.")

    want_all = bool(args.chains and args.chains.lower() == "all")
    wanted = None if (args.chains is None or want_all) else {c.strip() for c in args.chains.split(",") if c.strip()}

    is_cif = text.lstrip().startswith("data_") or "_atom_site." in text
    parsed = (
        parse_cif(text, wanted, want_all=want_all, include_all_models=(args.assembly > 0))
        if is_cif
        else parse_pdb(text, wanted, want_all=want_all)
    )
    baked = bake(parsed, align_mode=args.align, label=default_id, decimals=args.decimals)
    if args.title:
        baked["title"] = args.title

    os.makedirs(os.path.dirname(os.path.abspath(args.output)), exist_ok=True)
    payload = json.dumps(baked, separators=(",", ":"))
    if args.output.endswith(".json"):
        out_text = payload + "\n"
    else:
        out_text = f"// Generated by scripts/bake_pdb.py ({baked['id']}, {baked['residues']} residues)\nexport default {payload};\n"

    with open(args.output, "w", encoding="utf-8") as f:
        f.write(out_text)

    size_kb = len(out_text.encode("utf-8")) / 1024
    ch_ids = ",".join(c["id"] for c in baked["chains"])
    print(f"Baked {baked['id']} (chains: {ch_ids}, {baked['residues']} res) -> {args.output} ({size_kb:.1f} KB)")


if __name__ == "__main__":
    main()
