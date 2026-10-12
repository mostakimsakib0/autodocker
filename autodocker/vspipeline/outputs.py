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
# --box-from-het: the input receptor without the reference HETATM group (the
# file receptor preparation reads) and that group alone (input frame).
RECEPTOR_NO_REFERENCE_PDB = "receptor_no_reference.pdb"
REFERENCE_LIGAND_PDB = "reference_ligand.pdb"

# fpocket's own output directory (next to PROTEIN_CLEAN_PDB, named by fpocket)
FPOCKET_OUT_DIR = PROTEIN_CLEAN_PDB[:-len(".pdb")] + "_out"

# Docking
DOCKED_DIR = "docked"   # <name>_out.pdbqt poses and <name>_vina.log per ligand

# --replicates N > 1: replicate k of a ligand is <name>_rep<k> in DOCKED_DIR
# (poses, log, checkpoint key), docked from a copy of the ligand in
# DOCKED_DIR/REPLICATE_INPUTS_DIR; per-replicate scores and the run summary.
REPLICATE_TAG = "_rep"
REPLICATE_INPUTS_DIR = "replicate_inputs"
REPLICATE_SCORES_CSV = "replicate_scores.csv"
REPLICATE_SUMMARY_JSON = "replicate_summary.json"

# Outputs derived from a run's options and results. Every run removes the ones
# left by an earlier run in the same output directory before it starts, so a
# file that this run does not write (--no-html-report, no --consensus, no
# --cluster-poses, no fpocket with an explicit box, rigid docking, a failed
# report, ...) can never be mistaken for one of its results. The resumable
# docking results (DOCKED_DIR and the checkpoint) are never in this list.
DERIVED_OUTPUTS = (
    RANKING_CSV, TOP_HITS_TXT, METRICS_TXT, RESULTS_FULL_TXT, HTML_REPORT,
    POSE_CLUSTERS_JSON, CONSENSUS_RANKING_CSV, CONSENSUS_SUMMARY_JSON,
    POCKET_SUMMARY_CSV, FPOCKET_OUT_DIR, GRID_CONF, GRID_BOX_PML, FRAME_OFFSET_TXT,
    RIGID_RECEPTOR_PDBQT, FLEX_PDBQT, RECEPTOR_NO_REFERENCE_PDB, REFERENCE_LIGAND_PDB,
    REPLICATE_SCORES_CSV, REPLICATE_SUMMARY_JSON,
)

# Per-ligand files in DOCKED_DIR: <ligand name> + suffix
POSE_SUFFIX = "_out.pdbqt"                  # Vina/QuickVina poses
VINA_LOG_SUFFIX = "_vina.log"               # exact command + engine output
SMINA_POSE_SUFFIX = "_smina.pdbqt"          # --smina-only poses
SMINA_LOG_SUFFIX = "_smina.log"
VINA_TOP_POSE_SUFFIX = "_vina_top.pdbqt"    # --consensus: Vina's top pose
VINA_TOP_FLEX_SUFFIX = "_vina_top_flex.pdbqt"  # its flexible side chains
SMINA_RESCORE_SUFFIX = "_smina_rescore.pdbqt"  # --consensus: SMINA rescoring output

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
# Appended to RANKING_COLUMNS only with --replicates N > 1 (Binding_Affinity
# is then the best replicate's score; these summarise the in-window replicates).
COL_N_REPLICATES = "n_replicates"
COL_SCORE_MEAN = "score_mean"
COL_SCORE_SD = "score_sd"
COL_SCORE_MIN = "score_min"
COL_SCORE_MAX = "score_max"
RANKING_REPLICATE_COLUMNS = (
    COL_N_REPLICATES, COL_SCORE_MEAN, COL_SCORE_SD, COL_SCORE_MIN, COL_SCORE_MAX,
)

# Files whose SHA-256 run_manifest.json records under "outputs" (those that
# exist when the manifest is written), plus every file below
# HASHED_OUTPUT_DIRS (hidden files excluded). Not hashed: the manifest
# itself, pocket_summary.csv (fpocket's Monte-Carlo volume differs run to
# run), fpocket's own FPOCKET_OUT_DIR and the checkpoint/lock files.
HASHED_OUTPUTS = (
    RANKING_CSV, TOP_HITS_TXT, METRICS_TXT, RESULTS_FULL_TXT, HTML_REPORT,
    POSE_CLUSTERS_JSON, CONSENSUS_RANKING_CSV, CONSENSUS_SUMMARY_JSON,
    REPLICATE_SCORES_CSV, REPLICATE_SUMMARY_JSON,
    RECEPTOR_PDBQT, PROTEIN_CLEAN_PDB, RIGID_RECEPTOR_PDBQT, FLEX_PDBQT,
    RECEPTOR_NO_REFERENCE_PDB, REFERENCE_LIGAND_PDB,
    GRID_CONF, GRID_BOX_PML, FRAME_OFFSET_TXT, LIGAND_METADATA_JSON,
)
HASHED_OUTPUT_DIRS = (DOCKED_DIR, LIGANDS_DIR, LIGANDS_MINIMIZED_DIR)
