

import os
import shutil
import subprocess
from pathlib import Path
from __Modules__.Function import *

##### -------- ---------------- ---------------- ---------------- -------- A File class  -------- ---------------- ---------------- ---------------- -------- #####

class File(object):
    """
    Create an object that contains a path and a description
    Few methods :
        * Create an index if fasta
        * Align if FastQ (needs to be provided with an index)
        * Compute coverage
    """
    def __init__(self, input, description = ""):

        if Path(input).is_file():

           self.temp_dir = os.getcwd() + "/Temp/"
           self.bam = os.getcwd() + "/bam/"
           self.data = os.getcwd() + "/data/"
           self.reference = ""
           self.ledger = {
           "Description" : description,
           "Path" : input,
           "Library" : Path(input).stem.split(".")[0]
           }
        else :
            print("File doesn't exist !")

    def clean_temp(self):
        if Path(self.temp_dir).is_dir():
            shutil.rmtree(self.temp_dir)
        else :
            pass


##### -------- ---------------- ---------------- ---------------- -------- A Reference class  -------- ---------------- ---------------- ---------------- -------- #####

class Reference(File):

    def bwa_index(self):

        if check_ext(self.ledger["Path"], [".fasta", ".fna"]) :

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

##### -------- ---------------- ---------------- ---------------- -------- A Libray class  -------- ---------------- ---------------- ---------------- -------- #####

class Library(File):

    #inhirits from File

    def trim_adapters(self): 

        #Needs to make sure the tool is installed and in a useful directory, maybe use github ? !!!
        bbduk = f"{os.getcwd()}/__Modules__/bbmap/bbduk.sh"
        adapters = f"{os.getcwd()}/__Modules__/bbmap/resources/adapters.fa"

        if check_ext(self.ledger["Path"], [".fastq", ".fastq.gz"]) :

            #Catching files in Temp
            make_dir(self.temp_dir)

            #Create file to handle
            Trimmed = "{0}{1}_trimmed.fastq.gz".format(self.temp_dir , self.ledger["Library"])

            #Create bbduk functions
            adapter_trimming = "{0}  -Xmx2g in={1} out={2} ref={3} minlen=25 qtrim=rl trimq=10 ktrim=r k=25 mink=11".format(bbduk, self.ledger["Path"], Trimmed ,adapters)

            #Run the command
            subprocess.run(adapter_trimming , shell = True)

    def quality_trim(self): 

        #Needs to make sure the tool is installed and in a useful directory, maybe use github ? !!!
        bbduk = f"{os.getcwd()}/__Modules__/bbmap/bbduk.sh"

        if is_fastq(self.ledger["Path"]) :

            #Create file to handle
            Trimmed = "{0}{1}_trimmed.fastq.gz".format(self.temp_dir , self.ledger["Library"])
            Cleaned = "{0}Cleaned_{1}.fastq.gz".format(self.temp_dir , self.ledger["Library"])

            #Create bbduk functions
            cmd_qc= "{0} in={1} out={2} trimq=20 qtrim=rl k=23 minlen=35".format(bbduk, Trimmed, Cleaned)

            #Run the command
            subprocess.run(cmd_qc , shell = True)


    def trim_transposon(self, transposon):

        #Needs to make sure the tool is installed and in a useful directory, maybe use github ? !!!
        bbduk = f"{os.getcwd()}/__Modules__/bbmap/bbduk.sh"

        # Needs to make sure it's a fastq.gz file !
        if is_fastq(self.ledger["Path"]) :

            #Let's construct the bbduk commands
            Trimmed = "{0}{1}_trimmed.fastq.gz".format(self.temp_dir , self.ledger["Library"])
            Matched = "{0}TnMatched_{1}.fastq.gz".format(self.temp_dir , self.ledger["Library"])
            Unmatched =  "{0}UnMatched_{1}.fastq.gz".format(self.temp_dir , self.ledger["Library"])
            Cleaned = "{0}Cleaned_{1}.fastq.gz".format(self.temp_dir , self.ledger["Library"])

            #Command to pass
            cmd_filter= "{0} in={1} out={2} outm={3} literal={4} k=23 mm=f".format(bbduk, Trimmed, Unmatched , Matched , transposon)
            cmd_qc= "{0} in={1} out={2} ftl=31 trimq=20 qtrim=rl k=23 minlen=35".format(bbduk, Matched, Cleaned)
            
            #Run bbduk
            subprocess.run(cmd_filter , shell = True)
            subprocess.run(cmd_qc , shell = True)
        else :
            print("Not a FASTQ file")


    def get_reference(self, reference):
        ref_path = "{0}/Database/{1}.amb".format(os.getcwd(), str(reference))

        if Path(ref_path).is_file():
            #print("Index is found")
            self.reference = "{0}/Database/{1}".format(os.getcwd(), str(reference))
            return True
        else :
            #print("No index found")
            return False

    def bwa_aln(self, reference ,thread=4):
        #Make sure the index exists
        Cleaned = "{0}Cleaned_{1}.fastq.gz".format( self.temp_dir , self.ledger["Library"])
        Aligned = "{}.bam".format(self.bam + self.ledger["Library"])

        if self.get_reference(reference) and Path(Cleaned).is_file() :
            make_dir(self.bam)
            Output = "{}".format(self.bam, )
            #Prepare the commande
            cmd_bwa = "bwa mem -t{0} {1} {2} | samtools sort -@{0} -o {3}".format(thread, self.reference, Cleaned, Aligned)
            #Run the command
            subprocess.run(cmd_bwa , shell = True)
        else :
            print("Something went wrong")

    def get_TA(self):
        make_dir(self.data)
        Aligned = "{}.bam".format(self.bam + self.ledger["Library"])
        Output = "{}_TA.txt".format(self.data + self.ledger["Library"])

        if Path(Aligned).is_file():
            cmd_TA = "bedtools genomecov -d -5 -ibam {0} > {1}".format(Aligned,Output)
            subprocess.run(cmd_TA , shell = True)
        else :
            print("Something went wrong")

    def get_TnIF(self):
        make_dir(self.data)
        Aligned = "{}.bam".format(self.bam + self.ledger["Library"])
        Output = "{}_TnIF.txt".format(self.data + self.ledger["Library"])

        if Path(Aligned).is_file():
            cmd_TnIF = "bedtools genomecov -d -ibam {0} > {1}".format(Aligned,Output)
            subprocess.run(cmd_TnIF , shell = True)
        else :
            print("Something went wrong")