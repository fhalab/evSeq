# Import 3rd party modules
import numpy as np
import pandas as pd
import random
import os

# Get location of the globals.py file
GLOBALS_DIR = os.path.dirname(os.path.abspath(__file__))

# Range of BP ind start
BP_IND_START_MIN = 0
BP_IND_START_MAX = 400000

# Range of AA ind start
AA_IND_START_MIN = 0
AA_IND_START_MAX = 10000

# Reference sequence bounds (IN NUMBER OF AAS!!)
MIN_REFSEQ_LEN = 50
MAX_REFSEQ_LEN = 251

# Bounds on readlength
MIN_READLENGTH = 150
MAX_READLENGTH = 301
assert (MIN_REFSEQ_LEN * 3) <= MIN_READLENGTH

# Bounds on quality cutoffs
MIN_BP_QUAL_CUTOFF = 20
MAX_BP_QUAL_CUTOFF = 36

MIN_GLOBAL_QUAL_CUTOFF = 15
MAX_GLOBAL_QUAL_CUTOFF = 36

MAX_QUAL_ALLOWED = 41

# Bounds on sequence length cutoffs
MIN_SEQLEN_CUTOFF = 0.2
MAX_SEQLEN_CUTOFF = 0.8

# Bounds on number of indels added
MIN_INDELS_ADDED = 1
MAX_INDELS_ADDED = 4

# Bounds on primer lengths
PRIMER_MIN_LEN = 10
PRIMER_MAX_LEN = 41

# Bounds on the frameshift for a reference seq
FRAMESHIFT_MIN = 0
FRAMESHIFT_MAX = 3

# Bounds on position identification
MIN_VARIABLE_THRESH = 0.01
MAX_VARIABLE_THRESH = 0.3
MIN_VARIABLE_COUNT = 2
MAX_VARIABLE_COUNT = 21
assert MIN_VARIABLE_COUNT > 1 # We cannot have it at 1 or 0 for data generation to work

# Number of variants in a well
MIN_N_VARIANTS = 1
MAX_N_VARIANTS = 4

# Number of reads in a well
MIN_N_READS = 0
MAX_N_READS = 101

# Bounds on number of mutations per sequence (as a fraction
# of the number of amino acids captured by the read-length)
MIN_PERC_MUTATED = 0.0
MAX_PERC_MUTATED = 0.1

# Number of amino acids that can have noise added to them (as
# a fraction of the number of amino acids that have been mutated)
MIN_NOISE_PERC = 0.0
MAX_NOISE_PERC = 0.5

# Number of amino acids to rescue in overlapping regions (as
# a fraction of the total number of rescue options)
RESCUE_FREQ = 0.5

# Number of noisy reads added to the fastq
MIN_DUD_READS = 0
MAX_DUD_READS = 11

# Base for quality score calculation
Q_SCORE_BASE = 33

# How many mutations can be next to each other at max?
MAX_BORDERING_MUTS = 2

# How frequently do we convert a good well into a dead well?
DEAD_WELL_FREQ = 0.1

# Allowed nucleotides
ALLOWED_NUCLEOTIDES = ("A", "T", "C", "G")

# Random number generators
RANDOM_SEED = sum(ord(char) for char in "PDawg")
NP_RNG = np.random.default_rng(RANDOM_SEED)
RANDOM_RNG = random.Random(RANDOM_SEED)

# Barcode file
INDEX_DF = pd.read_csv(os.path.join(GLOBALS_DIR, "../../evSeq/util/index_map.csv"))
N_INDICES = len(INDEX_DF)
N_PLATES = len(INDEX_DF.IndexPlate.unique())

# Names of the columns in the refseq file
REFSEQ_COL_NAMES = (
    "PlateName",
    "IndexPlate",
    "Well",
    "FPrimer",
    "RPrimer",
    "VariableRegion",
    "FrameDistance",
    "BpIndStart",
    "AaIndStart"   
)

# Names of the columns in the output files
DECOUPLED_AA_COL_NAMES = (
    "IndexPlate",
    "Plate",
    "Well",
    "AaPosition",
    "Aa",
    "AlignmentFrequency",
    "WellSeqDepth",
    "Flags"
)

COUPLED_AA_COL_NAMES = (
    "IndexPlate",
    "Plate",
    "Well",
    "VariantCombo",
    "SimpleCombo",
    "VariantsFound",
    "AlignmentFrequency",
    "WellSeqDepth",
    "VariantSequence",
    "Flags"
)

# Get columns we will use for sorting
DECOUPLED_SORT_COLS = ["IndexPlate", "Well", "AaPosition", "Aa"]
COUPLED_SORT_COLS = ["IndexPlate", "Well", "SimpleCombo", "AlignmentFrequency"]

# Save location for test runs
SAVELOC = os.path.join(GLOBALS_DIR, "..", "test_data")