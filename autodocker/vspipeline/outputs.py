"""Names of the files and folders AutoDocker writes in the output directory.

Defined once so the modules that write a file and the ones that read it (and
the benchmark scripts, see scripts/_benchlib.py) cannot disagree.
"""

# Ranked results and summaries
RANKING_CSV = "ranking.csv"
TOP_HITS_TXT = "Top_hits.txt"
METRICS_TXT = "metrics.txt"
RESULTS_FULL_TXT = "Results_full.txt"
HTML_REPORT = "results_report.html"
POSE_CLUSTERS_JSON = "pose_clusters.json"
CONSENSUS_RANKING_CSV = "consensus_ranking.csv"
CONSENSUS_SUMMARY_JSON = "consensus_summary.json"
RUN_MANIFEST_JSON = "run_manifest.json"

# Ligand preparation
LIGANDS_DIR = "ligands"
LIGANDS_MINIMIZED_DIR = "ligands_minimized"
LIGAND_METADATA_JSON = "ligand_metadata.json"

# Receptor, box and frame
PROTEIN_CLEAN_PDB = "protein_clean.pdb"
RECEPTOR_PDBQT = "receptor.pdbqt"
RIGID_RECEPTOR_PDBQT = "rigid.pdbqt"   # receptor without the flexible side chains
FLEX_PDBQT = "flex.pdbqt"              # flexible side chains
GRID_CONF = "grid.conf"
GRID_BOX_PML = "grid_box.pml"
POCKET_SUMMARY_CSV = "pocket_summary.csv"
FRAME_OFFSET_TXT = "frame_offset.txt"

# Docking
DOCKED_DIR = "docked"   # <name>_out.pdbqt poses and <name>_vina.log per ligand

# ranking.csv columns, in file order
COL_LIGAND = "Ligand"
COL_AFFINITY = "Binding_Affinity"
COL_SIMSCORE = "SimScore"
COL_BINDING_MODES = "Binding_Modes"
COL_MODE_COVERAGE = "Mode_Coverage"
COL_MW = "MW"
COL_LOGP = "LogP"
COL_HBA = "HBA"
COL_HBD = "HBD"
COL_TPSA = "TPSA"
COL_ROTORS = "Rotors"
COL_SOURCE_URL = "ZINC_URL"   # historical name; holds the PubChem URL when known
COL_DOCKED_FILE = "Docked_File"
COL_STATUS = "Status"
COL_ENGINE = "Engine"
COL_MINIMIZED = "Minimized"
RANKING_COLUMNS = (
    COL_LIGAND, COL_AFFINITY, COL_SIMSCORE, COL_BINDING_MODES, COL_MODE_COVERAGE,
    COL_MW, COL_LOGP, COL_HBA, COL_HBD, COL_TPSA, COL_ROTORS, COL_SOURCE_URL,
    COL_DOCKED_FILE, COL_STATUS, COL_ENGINE, COL_MINIMIZED,
)
