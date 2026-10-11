"""Default settings of the pipeline, defined once.

The command line (runner.py), the library functions and the $VS_* environment
fallbacks all read these values, so a default can only be changed here. The
container entry point (entry.sh) passes an option only when its environment
variable is set, so it does not repeat them either.

This module imports nothing from the pipeline, so every other module can use
it without import cycles.
"""

# Docking (AutoDock Vina)
SEED = 42
EXHAUSTIVENESS = 8              # Vina's own default
BINDING_MODES = 9               # Vina's own default
ENERGY_RANGE = 3.0              # kcal/mol; Vina's own default
PROCESSES = 1                   # ligands docked in parallel (CLI)
COMMAND_TIMEOUT = 900           # seconds per external command
# Accepted affinity window: scores >= MIN_VALID_AFFINITY (or below the floor
# in docking.AFFINITY_FLOOR) are reported as OUT_OF_WINDOW.
MIN_VALID_AFFINITY = -1.0       # kcal/mol
MAX_FAILURE_FRACTION = 0.5      # exit 3 above this fraction of failed ligands

# Receptor and docking box
PH = 7.4                        # receptor protonation pH
PADDING = 6.0                   # Angstrom around the selected pocket(s)
WATER_DISTANCE = 4.0            # Angstrom; --keep-waters radius around the box center
AUTO_FLEX_RADIUS = 8.0          # Angstrom; --auto-flexible candidates around the box center
AUTO_FLEX_MAX_RESIDUES = 10

# Ligand library
MINIMIZE_STEPS = 250            # MMFF94 steps
MAX_COMPOUNDS = 20              # PubChem compounds prepared per run
CUSTOM_MW_MIN = 200.0           # Da, --library custom
CUSTOM_MW_MAX = 500.0           # Da, --library custom (Lipinski's limit)
CUSTOM_LOGP_MAX = 5.0           # --library custom (Lipinski's limit)

# Analysis and reports
TOP_N = 20                      # compounds in Top_hits.txt
RMSD_THRESHOLD = 2.0            # Angstrom, pose clustering
