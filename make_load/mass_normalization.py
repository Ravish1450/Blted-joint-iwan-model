# GENERALIZED MASS CHECK
MASS_FILE = "nodal_mass_two_parts.txt"
MODE_FILE = "modeshape.rpt"

# read nodal masses
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
print("total mass:", sum(nodal_mass.values()), "kg")

# read mode shape
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

# generalized mass
gm = 0.0
missing = 0

for key, m in nodal_mass.items():
    if key not in phi:
        missing += 1
        continue
    u1, u2 = phi[key]
    gm += m * (u1 * u1 + u2 * u2)

print()
print("nodes matched:", len(nodal_mass) - missing)
print("nodes with no phi:", missing)
print("GENERALIZED MASS =", gm)