
######
####
##
# Author : François-Xavier STUBBE 
# Institution : Université de Namur
#
##
####
######

##### -------- ---------------- ---------------- ---------------- -------- Libraries and packages needed  -------- ---------------- ---------------- ---------------- -------- #####

#Standard library packages

import os
import sys
import shutil 
import time 
import subprocess 
import threading 
from pathlib import Path
from tkinter import * 
from tkinter import filedialog as fd
import importlib.util

#Import my packages

from __Modules__.tkin2 import *
from __Modules__.Library import *
from __Modules__.Function import *
from __Modules__.Indexing import *
from __Modules__.TnSeq import *
from __Modules__.TnIF import *
from __Modules__.GFF_parser import *

##### -------- ---------------- ---------------- ---------------- -------- Libraries and packages needed  -------- ---------------- ---------------- ---------------- -------- #####



# -------- -------- BBAMP -------- -------- 

# BBMAP is tared into __Modules__
#Untar to make the tool available


bbamp_untar = f"tar -xf {os.getcwd()}/__Modules__/BBMap_39.01.tar.gz -C {os.getcwd()}/__Modules__/ "
subprocess.run(bbamp_untar , shell = True)

##### -------- ---------------- ---------------- ---------------- -------- Create the TnBox instance  -------- ---------------- ---------------- ---------------- -------- #####


root = Tk()
e = TnBox(root)
root.mainloop()



