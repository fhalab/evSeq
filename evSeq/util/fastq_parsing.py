import gzip
from pathlib import Path, PurePath
import tqdm
import glob
import random
from Bio import SeqIO

def trim_fastqs(folder, new_length=150):
    """A function to take a folder with fastqs that harbored returned reads 
    longer than desired and write new fastqs with the reads trimmed to 
    the desired length.

    Parameters:
    -----------
    folder: string, path to folder containing fastq files to trim. 
        fastq files to be trimmed must have the *_R*_*.fastq or 
        *_R*_*.fastq.gz naming convention. This will be the location of a new 
        folder called trimmed_reads harboring the trimmed fastqs, which 
        will be named trimmed_ORIGINAL_NAME.fastq.
    new_length: int, default 150, 
        length to trim reads to, will be used for both forward and reverse 
        read files. Reads shorter than the new_length will remain the same 
        length.
    """
    
    # get path for making folder within current folder for trimmed reads
    new_folder = Path(folder).joinpath('trimmed_reads')
    
    # make a directory to store the trimmed files
    try:
        new_folder.mkdir()
    except FileExistsError:
        print('Could not create directory because it already exists. New files will be written to the existing directory:', new_folder)
    
    # get the files with glob and pathlib
    files = [PurePath(file) for file in glob.glob(folder+'*_R*_*')] 
    
    # iterate through each file grabbed by glob --> will grab everything with *_R*_* in it
    for file in files:
        
        print('Loading reads from', file)
        trimmed_reads = []
        
        # check that the files grabbed are .fastq files
        if '.fastq' in file.suffixes:
        
            # open with gzip if still zipped
            if '.gz' in file.suffixes:
                with gzip.open(file, "rt") as handle:
                    for record in tqdm.tqdm(SeqIO.parse(handle, "fastq")):
                        
                        # trim read and append
                        new_record = record[:new_length]
                        trimmed_reads.append(new_record)
            
            # open regularly if not a .gz file
            else:
                with open(file, "rt") as handle:
                    for record in tqdm.tqdm(SeqIO.parse(handle, "fastq")):
                        
                        # trim read and append
                        new_record = record[:new_length]
                        trimmed_reads.append(new_record)
        
        # raise error if a .fastq or .fastq.gz is not the file type passed
        else:
            raise ValueError("You must pass a .fastq or a .fastq.gz file.")

        # make the new file path and write to file with Biopython SeqIO
        new_file_path = file.parent.joinpath('trimmed_reads')
        new_file = new_file_path.joinpath('trimmed_'+PurePath(file.stem).stem+'.fastq')
        SeqIO.write(trimmed_reads, new_file, "fastq")

def downsample_fastqs(folder, num_entries=100000):
    """A function to take a folder with fastqs that harbored more returned 
    reads than desired and write new fastqs with fewer reads that have been
    downsampled randomly. Downsampling the reads reduces the compute burden 
    of evSeq, but also reduces the data available.

    Parameters:
    -----------
    folder: string, path to folder containing fastq files to downsample. 
        fastq files to be downsampled must have the *_R*_*.fastq or 
        *_R*_*.fastq.gz naming convention. This folder will be the location 
        of a new folder called downsampled_reads harboring the downsampled 
        fastqs, which will be named downsampled_ORIGINAL_NAME.fastq.
    num_entries: int, default 100000, 
        number of reads to downsample to. This will be used for both forward 
        and reverse read files, which will have been downsampled in the same 
        way.
    """
    
    # get path for making folder within current folder for trimmed reads
    new_folder = Path(folder).joinpath('downsampled_reads')
    
    # make a directory to store the trimmed files
    try:
        new_folder.mkdir()
    except FileExistsError:
        print('Could not create directory because it already exists. New files will overwrite old ones and be written to the existing directory:', new_folder)
    
    # get the files with glob and pathlib
    if len(glob.glob(folder+'*_R1_*')) == 1 and len(glob.glob(folder+'*_R2_*')) == 1:
        fwd_file = PurePath(glob.glob(folder+'*_R1_*')[0])
        rev_file = PurePath(glob.glob(folder+'*_R2_*')[0])

    else:
        raise ValueError('Too many or too few .fastq files in directory.')     

    # blank dictionaries to hold reads      
    fwd_reads = {}
    rev_reads = {}
    ds_fwd_reads = []
    ds_rev_reads = []
    
    # check that the files grabbed are .fastq files
    if ('.fastq' not in fwd_file.suffixes) or ('.fastq' not in rev_file.suffixes):
        raise ValueError("You must pass a fwd and rev .fastq or a .fastq.gz file.")
    
    else:
    
        print('Loading forward reads...')
        # open with gzip if still zipped
        if '.gz' in fwd_file.suffixes:
            with gzip.open(fwd_file, "rt") as handle:
                for record in tqdm.tqdm(SeqIO.parse(handle, "fastq")):
                    fwd_reads[record.id] = record
        
        # open regularly if not a .gz file
        else:
            with open(fwd_file, "rt") as handle:
                for record in tqdm.tqdm(SeqIO.parse(handle, "fastq")):
                    fwd_reads[record.id] = record
        
        # downsample the record ids
        downsample_keys = random.sample(fwd_reads.keys(), num_entries)
        
        print('Downsampling forward reads...')
        for key in downsample_keys:
            ds_fwd_reads.append(fwd_reads[key])

        # make the new file path and write to file with Biopython SeqIO
        new_fwd_file_path = fwd_file.parent.joinpath('downsampled_reads')
        new_fwd_file = new_fwd_file_path.joinpath('downsampled_'+PurePath(fwd_file.stem).stem+'.fastq')
        SeqIO.write(ds_fwd_reads, new_fwd_file, "fastq")

        print('Loading reverse reads...')
        # open with gzip if still zipped
        if '.gz' in rev_file.suffixes:
            with gzip.open(rev_file, "rt") as handle:
                for record in tqdm.tqdm(SeqIO.parse(handle, "fastq")):
                    rev_reads[record.id] = record

        # open regularly if not a .gz file
        else:
            with open(rev_file, "rt") as handle:
                for record in tqdm.tqdm(SeqIO.parse(handle, "fastq")):
                    rev_reads[record.id] = record

        print('Downsampling forward reads...')
        # get the down-sampled forward reads
        ds_rev_reads = []
        for key in downsample_keys:
            ds_rev_reads.append(rev_reads[key])

        # make the new file path and write to file with Biopython SeqIO
        new_rev_file_path = rev_file.parent.joinpath('downsampled_reads')
        new_rev_file = new_rev_file_path.joinpath('downsampled_'+PurePath(rev_file.stem).stem+'.fastq')
        SeqIO.write(ds_rev_reads, new_rev_file, "fastq")