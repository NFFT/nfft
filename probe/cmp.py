import re, sys, collections
def load(log):
    g = collections.defaultdict(list); sec = "?"
    for line in open(log):
        s = line.rstrip()
        if s.endswith(":") and " " not in s: sec = s[:-1]; continue
        m = re.match(r"^(\S+)\s.*?(init\S*(?: \([^)]*\))?)\s*, m =\s*(\d+),\s*(\S+)\s+-> (OK|FAIL)\s+(\S+) \(\s*(\S+)\)", s)
        if not m: continue
        name, init, mm, op, _, err, bd = m.groups()
        mod = name.split("_")[0].split(".")[0]
        if op.endswith(("_1d", "_2d", "_3d")): continue
        g[(mod, sec.replace("check_", ""), op, int(mm))].append((float(err), float(bd)))
    return g
old, new = load(sys.argv[1]), load(sys.argv[2])
mod_filter = sys.argv[3] if len(sys.argv) > 3 else None
print(f"{'mod':4} {'section':20} {'op':8} {'m':>2} {'n':>3} {'old e/b':>7} {'new e/b':>7} {'new err':>7} {'new bd':>7}")
rows = []
for k in new:
    if k not in old or (mod_filter and k[0] != mod_filter): continue
    o = max(e / b for e, b in old[k]); n = max(e / b for e, b in new[k])
    e, b = max(new[k], key=lambda eb: eb[0] / eb[1])
    rows.append((k, o, n, e, b, len(new[k])))
for k, o, n, e, b, cnt in sorted(rows, key=lambda r: (r[0][0], r[0][3], r[0][1])):
    sec = k[1].replace("_file", "").replace("adjoint_", "adj_")
    print(f"{k[0]:4} {sec:20} {k[2][:8]:8} {k[3]:2d} {cnt:3d} {o:7.1e} {n:7.1e} {e:7.1e} {b:7.1e}")
