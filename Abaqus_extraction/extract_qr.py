# EXTRACT MODAL AMPLITUDE q_r FROM THE ODB

from __future__ import print_function
from odbAccess import openOdb
import numpy as np

# settings
ODB_FILE  = "Job-57.odb"          
STEP_NAME = "Step-3"
ALPHA_MAX = 2490.2                # Load that caused the onset of macroslip

MASS_FILE = "nodal_mass_two_parts.txt"
MODE_FILE = "modeshape.rpt"

OUT_FILE  = "qr_curve.txt"

# READ NODAL MASS AND MODE SHAPE
nodal_mass = {}
f = open(MASS_FILE, "r")
f.readline()
for line in f:
    line = line.strip()
    if not line:
        continue
    p = line.split(",")
    if len(p) == 3:
        nodal_mass[(p[0], int(p[1]))] = float(p[2])
f.close()

print("nodal masses read:", len(nodal_mass))

phi = {}
part = None

f = open(MODE_FILE, "r")
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
f.close()
print("mode shape nodes read:", len(phi))

# GLOBAL INDEX - one slot per node
keys = sorted(nodal_mass.keys())
index = {}
for i in range(len(keys)):
    index[keys[i]] = i

N = len(keys)
MP1 = np.zeros(N)
MP2 = np.zeros(N)

gm = 0.0
for i in range(N):
    k = keys[i]
    if k not in phi:
        continue
    m = nodal_mass[k]
    u1, u2 = phi[k]
    MP1[i] = m * u1
    MP2[i] = m * u2
    gm += m * (u1 * u1 + u2 * u2)

print("generalized mass check:", gm)
if abs(gm - 1.0) > 1e-4:
    raise SystemExit("STOP: generalized mass is not 1.0")

# tip node = largest |phi2|
tipkey = ("PART-1-1", 15)    # Selected by observing in Abaqus; node 15 of PART-1-1
tipidx = index[tipkey]
print("tip node:", tipkey, " phi2 =", phi[tipkey][1])


# OPEN ODB
odb = openOdb(ODB_FILE, readOnly=True)
step = odb.steps[STEP_NAME]
frames = step.frames
print("frames in", STEP_NAME, ":", len(frames))

# LOOP FRAMES
idxcache = {}
ref1 = None
ref2 = None
results = []
missing_reported = False

for n in range(len(frames)):
    fr = frames[n]
    blocks = fr.fieldOutputs["U"].bulkDataBlocks

    U1 = np.zeros(N + 1)
    U2 = np.zeros(N + 1)

    for blk in blocks:
        inst = blk.instance.name
        labels = np.asarray(blk.nodeLabels)
        data = np.asarray(blk.data)

        sig = (inst, len(labels), int(labels[0]), int(labels[-1]))

        if sig in idxcache:
            idx = idxcache[sig]
        else:
            idx = np.zeros(len(labels), dtype=int)
            miss = 0
            for i in range(len(labels)):
                k = (inst, int(labels[i]))
                if k in index:
                    idx[i] = index[k]
                else:
                    idx[i] = N
                    miss += 1
            if miss and not missing_reported:
                print("  warning: %d odb nodes not in mass file" % miss)
                missing_reported = True
            idxcache[sig] = idx

        U1[idx] = data[:, 0]
        U2[idx] = data[:, 1]

    U1 = U1[:N]
    U2 = U2[:N]

    if ref1 is None:
        ref1 = U1.copy()
        ref2 = U2.copy()

    d1 = U1 - ref1
    d2 = U2 - ref2

    q = float(np.dot(MP1, d1) + np.dot(MP2, d2))
    t = fr.frameValue
    alpha = ALPHA_MAX * t

    results.append((n, t, alpha, q, d2[tipidx]))

    if n % 10 == 0:
        print("  frame %3d  t=%.6f  alpha=%9.2f  q_r=%12.6e" % (n, t, alpha, q))
odb.close()

# WRITE OUT
f = open(OUT_FILE, "w")
f.write("frame,time,alpha,q_r,tip_u2\n")
for r in results:
    f.write("%d,%.8e,%.8e,%.8e,%.8e\n" % r)
f.close()

print()
print("written:", OUT_FILE, " rows:", len(results))

print()
print("low-amplitude check  (alpha / q_r should equal omega_0^2 = 1.452e6):")
for n in range(1, min(6, len(results))):
    a = results[n][2]
    q = results[n][3]
    if q != 0.0:
        print("  alpha=%8.2f   q_r=%12.6e   alpha/q_r=%12.6e" % (a, q, a / q))
