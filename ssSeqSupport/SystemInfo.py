# Import needed directories
import os.path
from multiprocessing import cpu_count

# Get the number of CPUs available on the computer
NCpus = cpu_count()
    
# Get the the running location of ssSeq and the logfile
Homedir = os.path.dirname(os.path.realpath(__file__))
Logfilename = os.path.join(Homedir, "ssSeqLog.log")