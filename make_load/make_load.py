# QSMA CLOAD

MASS_FILE = "nodal_mass_two_parts.txt"
MODE_FILE = "modeshape.rpt"
LOAD_OUT  = "load_2490.inp"         # output load or force file

ALPHA_MAX = 2490.2      # macroslip load

# READ NODAL MASSES
nodal_mass = {}

with open(MASS_FILE, "r") as f:
    next(f)
    for line in f:
        line = line.strip()
        if not line:
            continue
        p = line.split(",")
        if len(p) == 3:
            nodal_mass[(p[0], int(p[1]))] = float(p[2])

print("nodal masses read:", len(nodal_mass))

# Read mode shape
phi = {}
part = None

with open(MODE_FILE, "r", errors="ignore") as f:
    for line in f:
        l = line.strip()
        if "part: PART-1-1" in l:
            part = "PART-1-1"
            continue
        if "part: PART-2-1" in l:
            part = "PART-2-1"
            continue
        if part is None:
            continue
        p = l.split()
        if len(p) == 4:
            try:
                phi[(part, int(p[0]))] = (float(p[2]), float(p[3]))
            except ValueError:
                pass

print("mode shape nodes read:", len(phi))


# Check  (phi^T)*M*phi = 1 (should be one or close to one)
gm = 0.0
missing = 0

for key, m in nodal_mass.items():
    if key not in phi:
        missing += 1
        continue
    u1, u2 = phi[key]
    gm += m * (u1 * u1 + u2 * u2)

print("nodes with no phi:", missing)
print("generalized mass :", gm)

if missing > 0:
    raise SystemExit("mass file and mode shape file do not match.")

if abs(gm - 1.0) > 1e-4:
    raise SystemExit("generalized mass is not 1.0")


# Load  ( f = Alpha_max*m*phi )
lines = 0
fsum1 = 0.0
fsum2 = 0.0

with open(LOAD_OUT, "w") as f:
    f.write("*load\n")
    for (inst, node) in sorted(nodal_mass.keys()):
        m = nodal_mass[(inst, node)]
        u1, u2 = phi[(inst, node)]

        f1 = ALPHA_MAX * m * u1
        f2 = ALPHA_MAX * m * u2

        fsum1 += f1
        fsum2 += f2

        f.write("{}.{}, 1, {:.9e}\n".format(inst, node, f1))
        f.write("{}.{}, 2, {:.9e}\n".format(inst, node, f2))
        lines += 2

print()
print("  Load        :", LOAD_OUT)
print("  alpha_max   :", ALPHA_MAX)
print("  lines       :", lines)
print("  sum Fx      : {:.6e} N".format(fsum1))
print("  sum Fy      : {:.6e} N".format(fsum2))