# Load the parser utils file and other modules
import numpy as np
import pandas as pd
import pickle
from tqdm import tqdm
from multiprocessing import Pool
import gzip
import os
import os.path

# Import ssSeq modules
from ssSeq.GlobalSetup import *
import ssSeq.Classes as ss_utils

# Ignore divide by 0
np.seterr(invalid="ignore")

# First check if file paths are gzipped or not
gzipped_f = True if ForwardReads_Filepath[-2:] == 'gz' else False
gzipped_r = True if ReverseReads_Filepath[-2:] == 'gz' else False

# Assign proper parser
if gzipped_f:
    forward_parser = gzip.open(ForwardReads_Filepath, "rt")
else:
    forward_parser = open(ForwardReads_Filepath, "r")

if gzipped_r:
    reverse_parser = gzip.open(ReverseReads_Filepath, "rt")
else:
    reverse_parser = open(ReverseReads_Filepath, "r")

# Loop over the forward file and start building objects
with forward_parser as f:

    # Load the file to memory and count how many lines are in it
    file_list = list(f)
    n_lines_f = len(file_list)

    # Create a list to store objects
    seq_pairs = []

    # Create a dictionary so we can find seq_pairs later
    id_to_pair = {}

    # Loop over all lines, feeding blocks of 4 into a seq_pair object
    info_block = []
    print("Parsing forward read file...")
    for i, line in enumerate(tqdm(file_list), 1):

        # Strip whitespace from the line
        line = line.strip()

        # If this is a line divisble by 4, create a new seq_pair object and append to the growing list
        if i % 4 == 0:

            # Append the last needed line to the info block
            info_block.append(line)

            # Create the object and append
            seq_pair = ss_utils.seq_pair(info_block)
            seq_pairs.append(seq_pair)

            # Add information to a dictionary
            id_to_pair[seq_pair.id] = seq_pair

            # Reset the info block list
            info_block = []

        # If this is not divisible by 4, add to the info block list
        else:
            info_block.append(line)

# Identify the uids in seq_pairs
uids = set(id_to_pair.keys())

# Loop over the reverse file and append partners to each file
with reverse_parser as f:

    # Load the file to memory
    file_list = list(f)
    n_lines_r = len(file_list)

    # Initialize the orphan list. Allocate the total number of sequences.
    total_sequences = int(np.floor((n_lines_f + n_lines_r)/4))
    orphans = [None for _ in range(total_sequences)]
    orphan_ind = 0

    # Loop over each line in the file.
    info_block = [None, None, None, None]
    print("Parsing reverse reads...")
    for i, line in enumerate(tqdm(file_list), 1):

        # Strip whitespace from line
        line = line.strip()

        # Determine what line we're on
        line_check = i % 4

        # Generate the index for info_block
        block_ind = line_check - 1

        # Identify if this is a line divisible by 4. If so, append to the appropriate seq_pair
        if line_check == 0:

            # Complete the info block
            info_block[block_ind] = line

            # Identify the uid
            uid = ss_utils.create_id(id_line)

            # See if there is a match for the uid. If not, add to the orphan list
            if uid not in uids:
                orphans[orphan_ind] = ss_utils.seq_pair(info_block, start_f = False)
                orphan_ind += 1

            # If there is a match for the uid, append it to the appropriate partner
            else:
                id_to_pair[uid].attach_partner(info_block)

            # Reset the info block
            info_block = [None, None, None, None]

        # Identify if this is the first line. If so, this is the id line.
        elif line_check == 1:

            # Set the line as the id line
            id_line = line

            # Add to the info block
            info_block[block_ind] = line

        # If anything else, append to the info block
        else:
            info_block[block_ind] = line

# Remove the file from memory
del(file_list)

# Eliminate all "None" from the orphan list
orphans = [orphan for orphan in orphans if orphan is not None]

# Parse each pair and separate into plates and wells. Discard orphan pairs.
well_pairs = {}
for pair in seq_pairs:

    # Test to see if the well is an orphan. If it is, append to orphan list and move on
    if pair.is_orphan():
        orphans.append(pair)
        continue

    # Determine if the combination of forward and reverse barcode has been determined.
    # If the combo has not been seen yet, create a list in the dictionary
    well_ID = (pair.f_barcode, pair.r_barcode)
    if well_ID not in well_pairs.keys():
        well_pairs[well_ID] = [pair]

    # If the combo has been seen, add to the list that's growing
    else:
        well_pairs[well_ID].append(pair)


# Create wells for each list in the well_pairs dictionary
all_wells = [ss_utils.well(pair_list) for pair_list in well_pairs.values()]
failed_wells = [well for well in all_wells if well.real_well is not True]
true_wells = [well for well in all_wells if well.real_well is True]

# Parse each well and separate wells by plate
plate_wells = {}
for well in true_wells:

    # Determine if the plate has been seen yet. If not, create a list
    # in the dictionary
    if well.plate not in plate_wells.keys():
        plate_wells[well.plate] = [well]

    # If the plate has been seen, just append to the list.
    else:
        plate_wells[well.plate].append(well)

# Generate plates
plates = [ss_utils.plate(well_list) for well_list in plate_wells.values()]
n_plates = len(plates)

# Loop over each plate
n_plates = len(plates)
print("Processing plates...")
for i, plate in enumerate(plates):
    desc = "{} of {}".format(i + 1, n_plates)
    plate.process(i, ts_mode, n_jobs, desc)

# If in troubleshoot mode, save binaries of everything
if ts_mode:
    with open(output_location + "/Pickles/pickles.pkl", "wb") as f:
        pickle.dump([orphans, failed_wells, plates], f)

# Output a more concise final summary for each plate
for plate in plates:

    # Find output location of Variant Info
    file_output = os.path.join(output_location, "Summaries")

    # Get full Variant Info for current plate
    df = pd.read_csv(file_output+"/{}_VariantInfo.csv".format(plate.name))

    # Find the max value for Alignment Frequency for each well
    df_max = df.groupby('Well')['AlignmentFrequency'].max()

    # Merge with full DataFrame
    df_full = df.merge(df_max)

    # Round for easier reading
    df_full['AlignmentFrequency'] = np.round(df_full['AlignmentFrequency'].values, 2)

    # Save
    df_full.to_csv(file_output+"/{}_MaxInfo.csv".format(plate.name))