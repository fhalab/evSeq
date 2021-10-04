"""
Contains the wrapper for making a synthetic evSeq run for stress testing.
"""
# Import evSeq objects
from .config_generator import Config
from .well_generator import FakeWell
from .globals import (
    REFSEQ_COL_NAMES,
    INDEX_DF,
    SAVELOC,
    DECOUPLED_AA_COL_NAMES,
    COUPLED_AA_COL_NAMES
    )
from .stress_tests import calculate_parent_counts, check_well_is_parent

# Import 3rd party modules
import os
import subprocess
import numpy as np
import pandas as pd

# Class that holds information for a test run
class FakeRun():
    def __init__(self, detailed = True):
        """
        fakewells: A list of fully prepared FakeWell objects.
        detailed: Whether or not we are using a detailed refseq file. 
        """
        # Store instance variables
        self.detailed = detailed
        
        # Set a random configuration
        self.config = Config(detailed = detailed)
        
        # Build wells
        self._build_wells()
        
    def build_fastq(self):
        
        # Loop over all wells. Generate reads
        forward_reads = []
        reverse_reads = []
        for well in self.wells:
            
            # Produce reads if this is not a dud
            if not well.dud_well:
                f_well_reads, r_well_reads = well.build_all_reads()
                forward_reads.extend(f_well_reads)
                reverse_reads.extend(r_well_reads)
        
        # Return the fastq files ready for processing
        with open(self.r1_saveloc, "w") as f:
            f.write("\n".join(forward_reads))
        with open(self.r2_saveloc, "w") as f:
            f.write("\n".join(reverse_reads))
            
    def build_refseq(self, include_nnn):
        
        # Record the latest status of including NNN
        self.include_nnn = include_nnn
        
        # Loop over all wells. Only consider the first well of each
        # plate if not using detailed refseqs. Otherwise, consider
        # all wells
        n_refs = len(self.config.refseqs)
        observed_plates = set()
        refseq_output = [None] * n_refs
        counter = 0
        for well in self.wells:
                        
            # Create an entry for all wells if we are running in detailed mode
            # or if we have not observed the plate before
            if self.detailed or (well.platename not in observed_plates):
                
                # Add to the output
                refseq_output[counter] = well.build_refseq_entry(include_nnn)
                
                # Increment the counter and record the observed plate
                observed_plates.add(well.platename)
                counter += 1
        
        # Create a dataframe of the entries
        refseq_df = pd.DataFrame(refseq_output,
                                 columns = REFSEQ_COL_NAMES)
        
        # Sort the dataframe to be plate then row
        refseq_df.sort_values(by = [REFSEQ_COL_NAMES[1], 
                                    REFSEQ_COL_NAMES[2]],
                              inplace = True)
        
        # Drop the well column if this is not detailed
        if not self.detailed:
            refseq_df.drop(columns = ["Well"], inplace = True)
        
        # Save df 
        refseq_df.to_csv(self.refseq_saveloc, index = False)
        
    def _build_wells(self):
        
        # Build wells
        self.wells = [None] * len(INDEX_DF)
        unique_plates_found = set()
        refseq_counter = -1
        for i, row in enumerate(INDEX_DF.itertuples()):
            
            # Add only if detailed or we have not seen the plate before
            if self.detailed or (row.IndexPlate not in unique_plates_found):
                refseq_counter += 1
                unique_plates_found.add(row.IndexPlate)
        
            # Create a new well
            well = FakeWell(self.config,
                            self.config.refseqs[refseq_counter])
            
            # Assign information to the well
            well.platename = row.IndexPlate
            well.wellname = row.Well
            well.f_barcode = row.FBC
            well.r_barcode = row.RBC
            
            # Record
            self.wells[i] = well
            
    def build_output_counts(self):
        """
        Builds output files for the different `OutputCounts`
        """
        pass
    
    def run_evseq(self):
        """
        Wraps the evSeq command line interface to run the program with
        the provided fastq and refseq files
        """
        # Write the command for the subprocess. We save all outputs as this 
        # is for debugging
        command = [
            "evSeq",
            self.refseq_saveloc,
            SAVELOC,
            "--output", SAVELOC,
            "--return_alignments",
            "--keep_parsed_fastqs",
            "--average_q_cutoff", str(self.config.average_q_cutoff),
            "--bp_q_cutoff", str(self.config.bp_q_cutoff),
            "--length_cutoff", str(self.config.length_cutoff),
            "--variable_thresh", str(self.config.variable_thresh),
            "--variable_count", str(self.config.variable_count)
        ]
        
        # Add a flag for whether this is a detailed run or not
        if self.detailed:
            command.append("--detailed_refseq")
            
        # Run the command. Fail if errors are thrown.
        with open(self.test_stdout,"wb") as out, open(self.test_stderr,"wb") as err:      
            subprocess.run(command, 
                           check = True,
                           stdout = out,
                           stderr = err)
            
    def gather_expected_positions(self):

        # Extract positions we expect to be variable. This is given by the "NNN" 
        # sequence. If we are running a detailed test with NNN, then this is new for every
        # sequence; if not detailed, it is the refseq of the first sequence. If we are
        # not including NNN, then there is no expected variation
        if self.include_nnn:
            if self.detailed:
                expected_nnn_positions = {(well.platename, well.wellname): set(well.variants[0].mutated_positions)
                                          if not well.dud_well else set() for well in self.wells }
            else:
                expected_nnn_positions = {well.platename: set(well.variants[0].mutated_positions)
                                          if not well.dud_well and well.wellname == "A01"
                                          else set()
                                          for well in self.wells if well.wellname == "A01"}        
        else:
            if self.detailed:
                expected_nnn_positions = {(well.platename, well.wellname): set() for well in self.wells}
            else:
                expected_nnn_positions = {well.platename: set() for well in self.wells}

        # Extract all other positions where we added variability.
        expected_mut_positions = {}
        for well in self.wells:

            # Skip dud wells
            if well.dud_well:
                continue

            # Get all the mutated positions
            all_expected = set(np.concatenate([variant.mutated_positions for
                                               variant in well.variants]))

            # Get the nnn positions for the well
            platewell_key = (well.platename, well.wellname)
            nnn_key = platewell_key if self.detailed else well.platename
            well_nnn_pos = expected_nnn_positions[nnn_key] if self.include_nnn else set()

            # Record all positions that are not already captured by "nnn"
            expected_mut_positions[platewell_key] = all_expected - well_nnn_pos
            
        return expected_nnn_positions, expected_mut_positions
    
    # Build the expected decoupled results
    def _build_expected_decoupled_aa(self, well, nnn_positions, all_positions, is_parent):
        
        # If there are no positions OR the positions all give
        # the parent amino acid, this is a parent well
        if is_parent:

            # Record results
            return [[
                well.platename,
                well.platenickname,
                well.wellname,
                "#PARENT#",
                "#PARENT#",
                1.0, # We expect 100% alignment frequency given how we wrote the test code
                calculate_parent_counts(well),
                "#PARENT#",                    
            ]]

        # Loop over all positions and variants. For each mutated position, record the 
        # counts, the AA identity, the position, and any flags.
        all_pos_results = []
        for pos in all_positions:

            # Create output containers
            aa_to_count = {}

            # Determine position-uniform traits
            new_aa_pos = well.refseq.aa_ind_start + pos
            flag = "Unexpected Variation" if pos not in nnn_positions else np.nan

            # Loop over all variants
            for variant in well.variants:

                # Gather the AA id at the mutated position along with the expected 
                # counts
                mut_aa = variant.base_mut_aa_seq[pos]
                if mut_aa not in aa_to_count:
                    aa_to_count[mut_aa] = variant.expected_aa_counts[pos]
                else:
                    aa_to_count[mut_aa] += variant.expected_aa_counts[pos]

            # Get the frequency of all counts
            total_count = sum(aa_to_count.values())

            # If the only aa is the same as the parent and this is not a "NNN"
            # position, continue. There is no entry as evSeq will not find it.
            if (
                (len(aa_to_count) == 1) and
                (pos not in nnn_positions) and
                (mut_aa == well.refseq.aa_refseq[pos])
            ):
                continue

            # Record results
            all_pos_results.extend([
                [
                    well.platename,
                    well.platenickname,
                    well.wellname,
                    str(new_aa_pos),
                    aa,
                    spec_count / total_count,
                    total_count,
                    flag
                ]
                for aa, spec_count in aa_to_count.items()
            ])
            
        return all_pos_results

    # Build the expected coupled results
    def _build_expected_coupled_aa(self, well, nnn_positions, all_positions, is_parent):
        
        # If this is a parent well, format as appropriate
        if is_parent:
            
            return [[
                well.platename,
                well.platenickname,
                well.wellname,
                "#PARENT#",
                "#PARENT#",
                0,
                1.0,
                calculate_parent_counts(well),
                "".join(well.refseq.aa_refseq),
                "#PARENT#"
            ]]
                
        # Get a sorted list of all positions. 
        all_positions_sorted = sorted(list(all_positions))

        # Identify indices where the aa is (1) the same as the parent for all variants,
        # and (2) not an "NNN" position, continue. Remove these from the variant pool
        # as evSeq will not find it.
        no_use_positions = []
        for pos in all_positions_sorted:
            same_as_ref_check = all(variant.base_mut_aa_seq[pos] == well.refseq.aa_refseq[pos]
                                    for variant in well.variants)
            if (pos not in nnn_positions) and same_as_ref_check:
                no_use_positions.append(pos)
        no_use_positions = set(no_use_positions)
                
        adjusted_positions = [pos + well.refseq.aa_ind_start
                            for pos in all_positions_sorted
                            if pos not in no_use_positions]
            
            
        # Get variant counts and frequencies
        combo_counts = [variant.expected_combo_counts for variant in well.variants]
        total_counts = sum(combo_counts)
        frequencies = combo_counts / total_counts
            
        # Loop over all variants
        well_res = [None] * len(well.variants)
        for i, variant in enumerate(well.variants):
            
            # Grab the amino acid identities for both the variant and reference
            # at all positions. 
            all_ref_aas = [None] * len(all_positions_sorted)
            all_variant_aas = all_ref_aas.copy()
            for j, pos in enumerate(all_positions_sorted):
                all_ref_aas[j] = "?" if pos in nnn_positions else well.refseq.aa_refseq[pos]
                all_variant_aas[j] = variant.base_mut_aa_seq[pos]
            
            # Build the variant combo
            variant_combo = "_".join([f"{ref_aa}{pos}{mut_aa}" for ref_aa, pos, mut_aa in
                                    zip(all_ref_aas, adjusted_positions, all_variant_aas)])
                
            # Build the expected outputs
            well_res[i] = [
                well.platename,
                well.platenickname,
                well.wellname,
                variant_combo,
                "".join(all_variant_aas),
                len(adjusted_positions),
                frequencies[i],
                total_counts,
                "".join(variant.base_mut_aa_seq)
            ] 
            
        return well_res
    
    def build_expected_aa(self):
    
        # Get the expected mutated positions
        expected_nnn_positions, expected_mut_positions = self.gather_expected_positions()    

        # Create a list for storing expected results
        expected_decoupled_aa = []
        expected_coupled_aa = []

        # Loop over all wells
        for well in self.wells:

            # Skip dud wells
            if well.dud_well:
                continue

            # Gather expected mutation positions
            mut_key = (well.platename, well.wellname)
            nnn_key = mut_key if self.detailed else well.platename
            nnn_positions = expected_nnn_positions[nnn_key]
            mut_positions = expected_mut_positions[mut_key]
            all_positions = nnn_positions | mut_positions
            
            # If there are no positions OR the positions all give
            # the parent amino acid, this is a parent well
            is_parent = check_well_is_parent(well, all_positions, nnn_positions)
            
            # Build the expected decoupled results
            expected_decoupled_aa.extend(
                self._build_expected_decoupled_aa(well, nnn_positions, all_positions, is_parent)
            )
            
            # Build the expected coupled results
            expected_coupled_aa.extend(
                self._build_expected_coupled_aa(well, nnn_positions, all_positions, is_parent)
            )
                            
        # Format the output
        expected_decoupled_aa_df = pd.DataFrame(expected_decoupled_aa,
                                                columns = DECOUPLED_AA_COL_NAMES)
        expected_coupled_aa_df = pd.DataFrame(expected_coupled_aa,
                                              columns = COUPLED_AA_COL_NAMES)

        # Sort results
        expected_decoupled_aa_df.sort_values(by = ["IndexPlate", "Well", "AaPosition", "Aa"],
                                             inplace = True)
        expected_coupled_aa_df.sort_values(by = ["IndexPlate", "Well", "AlignmentFrequency"],
                                           inplace = True)
        
        return expected_decoupled_aa_df, expected_coupled_aa_df

    @property
    def refseq_saveloc(self):
        return os.path.join(SAVELOC, "testinput_refseq.csv")
    
    @property
    def r1_saveloc(self):
        return os.path.join(SAVELOC, "testinput_R1_allreads.fastq")
    
    @property
    def r2_saveloc(self):
        return os.path.join(SAVELOC, "testinput_R2_allreads.fastq")
    
    @property
    def test_stdout(self):
        return os.path.join(SAVELOC, "test_stdout.txt")
    
    @property
    def test_stderr(self):
        return os.path.join(SAVELOC, "test_stderr.txt")