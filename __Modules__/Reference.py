import os
import shutil
import subprocess
from pathlib import Path
#from __Modules__.Modules import *

#---------------------#
def check_ext(filename, ext):
    """
    Check the extension (only the last extension !)
    check_ext(filename, "fasta")
    """
    path = Path(filename)
    if path.suffix == ".{0}".format(str(ext)):
        return True
    else : False


def make_dir(dirname):
    dir = dirname
    if Path(dir).is_dir():
        pass
    else : os.mkdir(dir)

#---------------------#
class Reference(File):

    def bwa_index(self):

        if check_ext(self.ledger["Path"], "fasta") :

            #Create a directory with bwa_index
            index_dir = os.getcwd() + "/Database/"
            make_dir(index_dir)
            #Prepare the index
            index_p = index_dir +  Path(self.ledger["Path"]).stem
            cmd = "bwa index -p{0} {1}".format(index_p, self.ledger["Path"])
            #Run BWA
            subprocess.run(cmd , shell = True)
        else :
            print("This is not a Fasta file")