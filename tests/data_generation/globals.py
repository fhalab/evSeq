# Import 3rd party modules
import numpy as np
import pandas as pd
import random

# Import codon table for support
from .data_generation_utils import CustomCodonTable

# Reference sequence bounds (IN NUMBER OF AAS!!)
MIN_REFSEQ_LEN = 50
MAX_REFSEQ_LEN = 251

# Bounds on readlength
MIN_READLENGTH = 150
MAX_READLENGTH = 301
assert (MIN_REFSEQ_LEN * 3) <= MIN_READLENGTH

# Bounds on quality cutoffs
MIN_BP_QUAL_CUTOFF = 15
MAX_BP_QUAL_CUTOFF = 36

MIN_GLOBAL_QUAL_CUTOFF = 15
MAX_GLOBAL_QUAL_CUTOFF = 36

MAX_QUAL_ALLOWED = 41

# Bounds on sequence length cutoffs
MIN_SEQLEN_CUTOFF = 0.0
MAX_SEQLEN_CUTOFF = 1.0

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
MIN_VARIABLE_THRESH = 0.0
MAX_VARIABLE_THRESH = 1.0
MIN_VARIABLE_COUNT = 0
MAX_VARIABLE_COUNT = 21

# Number of variants in a well
MIN_N_VARIANTS = 1
MAX_N_VARIANTS = 4

# Number of reads in a well
MIN_N_READS = 0
MAX_N_READS = 101

# Bounds on number of mutations per sequence (as a fraction
# of the number of amino acids captured by the read-length)
MIN_PERC_MUTATED = 0.0
MAX_PERC_MUTATED = 0.2

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

# Allowed nucleotides
ALLOWED_NUCLEOTIDES = ("A", "T", "C", "G")

# Codon table and allowed amino acid characters
CODON_TABLE = CustomCodonTable()
ALLOWED_AAS = tuple(sorted(list(CODON_TABLE.aa_to_codon.keys())))
INT_TO_AA = dict(enumerate(ALLOWED_AAS))
N_AAS = len(ALLOWED_AAS)

# Random number generators
RANDOM_SEED = sum(ord(char) for char in "PDawg")
NP_RNG = np.random.default_rng(RANDOM_SEED)
RANDOM_RNG = random.Random(RANDOM_SEED)

# Barcode file
INDEX_DF = pd.read_csv("../evSeq/util/index_map.csv")
N_INDICES = len(INDEX_DF)
N_PLATES = len(INDEX_DF.IndexPlate.unique())