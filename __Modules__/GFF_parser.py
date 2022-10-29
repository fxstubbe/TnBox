import os
import shutil
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from pathlib import Path
from __Modules__.Library import *


##### -------- ---------------- ---------------- ---------------- -------- A GFF parser  -------- ---------------- ---------------- ---------------- -------- #####


def get_attribute(data, atrfield):
    attributes = []
    my_attributes = data["attribute"].str.split(";")
    for index, row in my_attributes.iteritems():
        attribute = list(filter(lambda a: atrfield in a, row))
        if len(attribute) :
            attributes.append(attribute[0].replace(f"{atrfield}=", ""))
        else :
            attributes.append("NA")
    data[atrfield] = attributes

def parse_gff(my_gff):
    #Read the GFF
    gff = pd.read_csv(my_gff, sep="\t", header=None,comment="#", names=("seqname", "source", "feature", "start","end", "score", "strand", "frame", "attribute"),)
    #Select desirable attributes
    for atr in ["gene_biotype", "Name", "ID", "Parent", "locus_tag", "old_locus_tag", "product"]:
        get_attribute(gff, atr)
    #Claean up a little bit
    for index,row in gff.iterrows():
        if gff.iloc[index, gff.columns.get_loc("feature")] == "CDS" :
            index_loc = gff.index[gff["ID"] == gff.iloc[index, gff.columns.get_loc("Parent")]]
            gff.iloc[index_loc, gff.columns.get_loc("product") ] = gff.iloc[index,gff.columns.get_loc("product")]
        else : pass
    #Print the file
    genes = gff[gff["feature"].isin(["gene", "pseudogene"])]
    genes = genes.drop(["score", "frame", "attribute", "source", "ID", "Parent"],1)
    #Print the file
    index_dir = os.getcwd() + "/Annotation/"
    file_name = Path(my_gff).stem
    make_dir(index_dir)
    genes.to_csv(f"{index_dir}{file_name}.csv")
