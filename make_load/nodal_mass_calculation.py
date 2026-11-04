from collections import defaultdict

# INPUT FILES

INPUT_FILE = "connectivity.inp"
EVOL_FILE = "evol.rpt"
OUTPUT_FILE = "nodal_mass_two_parts.txt"

DENSITY = 7850.0       # kg/m^3
NODES_PER_ELEMENT = 4  # CPE4

# STEP 1: READ ELEMENT CONNECTIVITY

connectivity = {}

current_instance = None
inside_element_block = False
with open(INPUT_FILE, "r") as f:
    for line in f:
        line = line.strip()
        if not line:
            continue
        # Detect instance
        if line.lower().startswith("*instance"):
            if "name=part-1-1" in line.lower():
                current_instance = "PART-1-1"
            elif "name=part-2-1" in line.lower():
                current_instance = "PART-2-1"
            else:
                current_instance = None
            inside_element_block = False
            continue

        # End of instance
        if line.lower().startswith("*end instance"):
            current_instance = None
            inside_element_block = False
            continue

        # Detect CPE4 element block
        if line.lower().startswith("*element"):
            if "type=cpe4" in line.lower() and current_instance is not None:
                inside_element_block = True
            else:
                inside_element_block = False
            continue

        # Another keyword ends element block
        if line.startswith("*"):
            inside_element_block = False
            continue

        # Read element connectivity
        if inside_element_block and current_instance is not None:
            values = [v.strip() for v in line.split(",")]
            if len(values) >= 5:
                try:
                    element = int(values[0])
                    nodes = [
                        int(values[1]),
                        int(values[2]),
                        int(values[3]),
                        int(values[4])
                    ]
                    connectivity[(current_instance, element)] = nodes

                except ValueError:
                    pass
print("Connectivity read:")
print("  PART-1-1 elements:",
      sum(1 for key in connectivity if key[0] == "PART-1-1"))

print("  PART-2-1 elements:",
      sum(1 for key in connectivity if key[0] == "PART-2-1"))


# STEP 2: READ EVOL
evol = {}
current_instance = None
with open(EVOL_FILE, "r") as f:
    for line in f:
        line = line.strip()
        if not line:
            continue
        upper = line.upper()
        # Detect EVOL part heading
        if "FIELD OUTPUT REPORTED FOR WHOLE ELEMENT FOR PART:" in upper:
            if "PART-1-1" in upper:
                current_instance = "PART-1-1"
            elif "PART-2-1" in upper:
                current_instance = "PART-2-1"
            else:
                current_instance = None
            continue

        # Read element + EVOL
        if current_instance is not None:
            parts = line.split()
            if len(parts) >= 2:
                try:
                    element = int(parts[0])
                    volume = float(parts[1])
                    evol[(current_instance, element)] = volume
                except ValueError:
                    pass

print("EVOL read:")
print("  PART-1-1 elements:",
      sum(1 for key in evol if key[0] == "PART-1-1"))

print("  PART-2-1 elements:",
      sum(1 for key in evol if key[0] == "PART-2-1"))


# CHECK TOTAL MASS DIRECTLY FROM EVOL

total_element_mass = 0.0
for key, volume in evol.items():
    total_element_mass += DENSITY * volume

print("\nDirect mass from EVOL =", total_element_mass, "kg")

# End of verification step.

# STEP 3: CALCULATE NODAL MASS

nodal_mass = defaultdict(float)
missing_evol = 0
missing_connectivity = 0

for key, nodes in connectivity.items():

    instance, element = key

    # Check EVOL exists
    if key not in evol:

        missing_evol += 1
        continue

    volume = evol[key]

    # Element mass
    element_mass = DENSITY * volume

    # Equal mass distribution to four CPE4 nodes
    mass_per_node = element_mass / NODES_PER_ELEMENT

    # Add contribution to each node
    for node in nodes:

        nodal_mass[(instance, node)] += mass_per_node

# STEP 4: OUTPUT
with open(OUTPUT_FILE, "w") as f:
    f.write("Part,Node,Nodal_Mass_kg\n")
    for (instance, node), mass in sorted(nodal_mass.items()):
        f.write(
            "{},{},{:.12e}\n".format(
                instance,
                node,
                mass
            )
        )

# STEP 5: CHECKS

part1_nodes = sum(
    1 for key in nodal_mass
    if key[0] == "PART-1-1"
)

part2_nodes = sum(
    1 for key in nodal_mass
    if key[0] == "PART-2-1"
)

total_mass = sum(nodal_mass.values())

print("NODAL MASS CALCULATION")
print("PART-1-1 nodes:", part1_nodes)
print("PART-2-1 nodes:", part2_nodes)
print("Total physical nodes:", len(nodal_mass))
print("Missing EVOL:", missing_evol)
print("Total nodal mass:", total_mass, "kg")
print("\nOutput:")
print(OUTPUT_FILE)