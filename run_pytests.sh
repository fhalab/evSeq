#!/bin/bash
# Handles running pytests. Will activate conda environments and skip relevant
# folders that we don't want included in tests.

# Activate the testing conda environment
source ~/anaconda3/etc/profile.d/conda.sh
conda activate evSeq_exact

# Run pytests, skipping appropriate directories
pytest ./tests/*.py --ignore ./tests/test_data/ --ignore ./tests/data_generation/