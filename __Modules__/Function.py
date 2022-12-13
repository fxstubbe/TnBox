import os
import shutil
import numpy as np
import pandas as pd
import seaborn as sns
from pathlib import Path
from __Modules__.Library import *



##### -------- ---------------- ---------------- ---------------- -------- A few all purpose functions  -------- ---------------- ---------------- ---------------- -------- #####

def check_ext(filename, ext):
    """
    Check the extension (only the last extension !)
    check_ext(filename, "fasta")
    """
    path = Path(filename)
    if path.suffix == ".{0}".format(str(ext)):
        return True
    else : False


def is_fasta(filename):
 """
 Check if the extension is fasta or fna
 """
 path = Path(filename)
 try:
     ext = path.suffix
     if ext == ".fasta" or ext == ".fna":
         return True
     else : return False
 except IndexError:
     return False


def is_fastq(filename):
    """
    Check if the extension is fastq.gz
    """
    path = Path(filename)
    print(path.suffix)
    if path.suffix == ".fastq" : 
        print("File is fastq")
        return True
    else : 
        try:
            ext = "{0}{1}".format(path.stem.split('.')[1] , path.suffix)
            if ext == "fastq.gz":
                print("File is fastq.gz")
                return True
            else : 
                print("Wrong extension")
                return False
        except IndexError:
            return False

def check_ext(filename_p, extensions, comp = False) :

    filepath = Path(filename_p)
    print( f"suffix = {filepath.suffixes}")
    print(f"extensions checked = {extensions}")
    
    # increment for each extension that match
    catch = 0

    if any((match := item) in extensions for item in filepath.suffixes): 
            catch += 1
            if comp is True and ".gz" not in filepath.suffixes : 
                catch -= 1
                
    # Check if an extension was valid
    if catch > 0 :
        print("File has valid extension")
        return True
    else : 
        print("File of unvalid extension")
        return False

def make_dir(dirname):
    dir = dirname
    if Path(dir).is_dir():
        pass
    else : os.mkdir(dir)


def ConvertTuple(my_tuple):

    # Convert a tuple of string into a tuple of integers
    # If an element cannot be converted or is missing return a None Object

    Flag = 0 # Catch flags
    int_metrics = list()

    for i, item in enumerate(my_tuple):

        try :
            int_metrics.append(abs(int(item)))

        except ValueError :
            Flag += 1

    if Flag == 0 :
        return tuple(int_metrics)
    else : 
        return None


def switchButtonstate(button):
    if button['state'] == NORMAL :
        button['state'] = DISABLED
    else : 
        button['state'] = NORMAL





