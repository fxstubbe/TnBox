# TnBox 

## A Tn-seq Tool Box

Since its introduction in micriobiology labs, massively parallel sequencing coupled with transposson mutagenesis (Tn-seq) has become an important tool for molecular microbiologist to gain insight into bacterial fitness. Tools for the anaylsis of such data exist (e.g. [Tn-seq Explorer][explorer]) but they're not trivial to use for non bioinformatician. 

TnBox is a user friendly annotation based software for the analysis of high throughput mutant libraries (Tn-seq). It's aim is to provide an easy to use tool enabling wet lab researchers to gain insight into their Tn-seq data. It provides a set of tools allowing the user to parse and align sequencing reads on a prokaryotic genome. It then implements two algorithms (Rslide & TnIF) for the study of gene essentiality and bacterial fitness and converts the computed output into comprehensive excel spreadsheets. Furthermore, visualisation tools are provided to allow the user to quickly explore the computed data.


TnBox has been designed to analyse miniTn5 libraries provided by the FASTERIS [in-house transposon sequencing service][fasteris]. Even if, in principle, any Tn-seq should work with this pipeline, we have not tested other libraries type.

For support, questions or requests, please contact: francois-xavier.stubbe@unamur.be

## Installation & quick start

TnBox comes as a simple python program that works on both macOS and Linux systems. However, a few dependencies need to be installed and set in your PATH. Similarly, a couple python packages need to be downloaded. There are many ways to install these but the eaisest way might be to use the python package manager [Anaconda][conda]. Since Anaconda comes along a version of python, it is not necessary to have python pre-installed on your machine.


### Clone TnBox

For simplicity of use, we will download TnBox on the deskop. Firstly, open a terminal an place yourself on the desktop. This can be done using cd (change directory)  as follow

```sh
cd ./Desktop
git clone https://github.com/fxstubbe/TnBox
```

Once on the Desktop, download the repository with the following command

```sh
git clone https://github.com/fxstubbe/TnBox
```

Alternatively, you can also download TnBox using the green button at page's top.

### Install depedencies and packages
### #Creating a virtual environment (recommended)

It is good practice to create a separate virtual environment for each project. An environment installation file is provided with tnbox. Open a terminal, go into the tnbox repository and create the environment using conda. 

```sh
cd Desktop/TnBox/
```

The environment contains all the tools that TnBox needs to analyze your Tn-seq experiment.

```sh
conda env create -f environment_setup.yml
```

Before using TnBox, remember to activate the virtual environment. Simply paste the folllowing command in your terminal 

```sh
conda activate tnbox
```

To close the tnbox environment, simply run

```sh
conda deactivate 
```



#### Installation in main environment (not recommended)
##### Install dependencies with conda on linux/macOS systems

- [Burrow Wheeler Aligner][bwa] 
- [Samtools][samtools]
- [Bedtools][bedtools]
- [Git][git]
- [BioPython][biopython]

Paste the following commands in a terminal

```sh
conda config --add channels bioconda
conda config --add channels conda-forge
conda install -y git
conda install -y biopython
conda install -y bwa
conda install -y bedtools
conda install -y samtools=1.9
```
##### Install python packages with conda on linux/macOS systems

- [Numpy][python]
- [Pandas][bwa] 
- [Matplotlib][bwa] 
- [Searborn][samtools]

Paste the following commands (one by one) in a terminal

```sh
conda install -y -c anaconda numpy
conda install -y -c anaconda pandas
conda install -y -c anaconda matplotlib
conda install -y -c anaconda searborn
```


### Run TnBox 

##### Launch TnBox

TnBox is a program that needs to be launched from the terminal. Opeen a terminal and direct yourself in the TnBox folder. Assuming you've installed TnBox on your Desktop, use the command below. If you installed TnBox elsewhere, then you'll need to use the folder path.

```sh
cd Desktop/TnBox/
```

If you created a virtual environment using the environment_setup.yml file, then you need to activate the environment. Simply paste the following command in the terminal.

```sh
conda activate tnbox
```

You're now ready to start TnBox. Simply use python to launch the program.

```sh
python TnBox.py
```

If TnBox interface opened ... congratulations ! TnBox is now working on your machine !

## Analysis

![](https://github.com/fxstubbe/TnBox/blob/main/Images/Tnbox.png)

### 1. Process libraries

The first step int any high throughput sequencing project is to map the sequencing reads onto a genome of interest. To do so, reads must first be trimmed of adapter and quality trimmed. Then, the processed reads are aliged onto an indexed genome using an aligner. TnBox uses [BBMAP][bbmap] for the filtering and [bwa][bwa] for the alignment. 

#### 1.1 Transposon fishing

We recommend sequencing of miniTn5 libraries with the [in-house transposon sequencing service][fasteris] provided by fasteris.

In the past, Tn-seq libraries where often sequenced on Illumina HiSeq intruments. With those machines a nested strategy using a specific primer matching the 3'end of the transposon. Therefore, the first read base was the transposon insertion site. On newer machines (e.g. NextSeq, NovaSeq), sequencing reactions with specific primers aren't reliable. Libraires are now sequenced from the illumina adapters, which then re-sequence the end of the transposon. To separate specific from aspecific reads, it is necessary to separate reads containing the transposon from reads that do not (up to 35% of the library). By checking the tickbox "miniTn5", TnBox will filter transposon containing reads and trim out the transposon sequence prior to mapping (bottom panel, image below). If omitted, lots of non-specific reads will map and weakened downstrean analysis. 

Both Mariner and miniTn5 are available to parse. However, this methode has only been used and confirmed on miniTn5 libraires.

![](https://github.com/fxstubbe/TnBox/blob/main/Images/aspecific_mapping.png)

#### 1.2 Select or Add a reference

As previoulsy stated, sequencing reads need to be mapped onto a reference genome. Genomes (.fasta , .fna) can be dowaloaded from NCBI (e.g. [*Brucella abortus*][abortus]). A few references are provided with TnBox. Once added, your reference genome will appear in the list. Select the genome of interest (highlighted in blue) and proceed to the enxt step. 


#### 1.3 Add sequencing files 

Add your sequencing reads (.fastq.gz, .fastq) otherwise TnBox will fail. Tn-seq are usually sequence in single-end but, if your reads are paired-end, only provide the forward reads (R1) as this is where the transposon lies. Similarly, sometimes, sequencing reads will be split over mutiple file. You can simply concatenate them together as follow : 

```sh
cat file_1 file_2 > concatenated_file
```

#### Start the anaysis

Click on start. That's it. If nothing starts, make sure all the previous steps are correctly fulfilled.

This process is time consuming. It largely depends on the amount of reads, the genome size ... but also your machine performance. Be patient, make yourself a cofee and come back later. TnBox provides a visual cue of where it is in the process. When the library name is orange, it's under process. When it turns green, TnBox is done with this library.

*FYI : TnBox will produced both an aligmnet file (.bam) and two coverage file (more details below). Those are respectively stored in ./TnBox/bam and ./TnBox/data*

### 2. Get transposons insertion sites 


`Rslide or Sliding Window approach`

This approach implements an alternative to the R200 metric introduced by [Sternon et al.,][sternon]. In their approach, the coverage file is split into windows of size 200 that are slide by 5 nucleotide. For example a 3 278 307 bp genome is split into 655,662 windows. For each window, the coverage sum is computed. This method has the advantage of being annotation independent and could be used to re-annotate genomic features as well as identifying essential protein regions.


The TnBox Rslide algorithm has for objective to provide an insight in wether or not a given gene is  essential in a given condition (preferably growth on plate). To do so, it computes the R window strategy described above. It then assigns for each gene and *Essentiality score* equivalent to the number of windows having a sum of 0 (aka, no transposons jumped in that genomic region). In other words; a gene having an Essentiality score of 20 means that , in that gene, there are 20x 200nt wide windows without any transposons. We consider that having a R200 score > 0 indicates the gene as essential.

`TnIF or Insertion Density approach`

The `TnIF` algorithm which has for goal to detect essentiality variations across several condition. To do so, TnBox computes the *transposon insertion freqeuncy (TnIF)* defined by [Potemberg et al., ][potemberg]. Briefly, for each gene the average log10(r+1)/l (where r is the read count on a given nucleotide and l is the gene length) is computed. For efficient comparison, an indexed table can be generated (see index section). 

#### 2.1 Choose or Add a reference (.gff)
Now that the reads have been aligned to the genome of interest, it's time extract where the transposons are. The first step is to select the appropirate reference or the add your reference of choice in gff format. For example, you can dowload Brucella abortus reference [here][abortus]. 

#### 2.2 Define metrics for the algorithms

Both the `Rslide or Sliding Window approach` and `TnIF or Insertion Density approach` can be fine tuned with different parameters.


Parameter  | algorithm | Default | Description
------------- | ------------- | ------------- | ------------- 
`Trim End`  |  RSlide, TnIF |  10% | 5' and 3'end trimming (10% = computing over central 80%)
`R Windows` | Rslide | 200 nt | Size of the sliding window
`Slide` | Rlside | 5 nt | increment between each sliding window




#### 2.3 Select files and algorithm

When processing the libraries, TnBox generated 2 types of files :
-   **TA (anchor) files**  
When the coverage was computed, only the read  5'end was kept. Therefore, only the transposon insertion site (TA) is kept. The term *TA* is used in reference to the mariner transposon which insert in TA rich regions. This is the preferred type for the `Rslide` algorithm.

-  **TnIF (coverage) files** 
Instead of only mapping the 5'end, the whole read is counted in the coverage file. Even if the the method is slightly less sensistive for the detection of essential genes with tge `Rslide` algorithmm, it performs better when using the `TnIF`. This is the preferred type for the `TnIF` algorithm. Please do not use the `TnIF` if your mutational library isn't saturating.

First and foremost, make sure you have selected the right reference in the previous section. Then, select the files in the listbox of interest (TA or TnIF). Simply press on the algorithm of choice and choose where you wanna save the output table. In the table, each selected file will be added as a column (genes are rows). Once the button unlock, the file has been saved. Feel free to launch multiple algorithms simultanously (might significantly slow down your computer).

### 3. Indexing & Delta (TnIF only)

If you're using the `Rslide` algorithm, skip this panel.

#### 3.1 Indexing

TnIf reprenset a bimodal distribution where the second peak of the distribution correspond to non-essential genes. To allow correct library comparison, it is best to *index* your library(ies) of interest on a control. To do so, an easy method is to simply divide by the mode of the second peak. I could easily be done manually but to simplify  the process, TnBox implements an indexing algorithm. 

Simply load the file you desire to index, choose the library you wish to use as reference and start the indexing. 

It is to note that indexing will perform poorly on libraries with low diversity (e.g. bottleneck effect). 

#### 3.1 Delta

For each gene, a ΔTnIF (TnIFcdt−TnIFCTRL) value was calculated, where TnIF was computed for the tested condition (TnIFcdt) and the control condition (TnIFCTRL). The frequency distribution of ΔTnIF values was plotted for both chromosomes and for each tested condition (S6 Fig), to identify the main peak of unaffected ΔTnIF values and its standard deviation. 2% of ΔTnIF values at each extremity were removed to avoid an influence of extreme values, the standard deviation was calculated on this distribution. Depending on the conditions tested, the standard deviation ranged from 0.049 and 0.244. The ΔTnIF values larger than 0.5 were thus selected as significant, since they correspond to 2 to 5 standard deviations from the mode, designating genes for which the TnIF value was decreased compared to the control condition.

### 4. Explore

![](https://github.com/fxstubbe/TnBox/blob/main/Images/comp.png)

## References

There will be a couple of references

## License

MIT

**Free Software, Hell Yeah!**

[//]: # (These are reference links used in the body of this note and get stripped out when the markdown processor does its job. There is no need to format nicely because it shouldn't be seen. Thanks SO - http://stackoverflow.com/questions/4823468/store-comments-in-markdown-syntax)

[bwa]: https://sourceforge.net/projects/bio-bwa/files/
[samtools]: http://www.htslib.org/
[bedtools]: https://bedtools.readthedocs.io/en/latest/
[python]: https://www.python.org/
[conda]: https://www.anaconda.com/
[git]: https://git-scm.com/
[biopython]: https://biopython.org/
[pandas]: https://pandas.pydata.org/
[numpy]: https://numpy.org/
[plt]: https://matplotlib.org/
[sns]: https://seaborn.pydata.org/
[bbmap]: https://github.com/BioInfoTools/BBMap
[abortus]:https://www.ncbi.nlm.nih.gov/genome/520
[fasteris]:https://www.fasteris.com/en-us/NGS/DNA-sequencing/Tn-Seq
[sternon]:https://journals.asm.org/doi/10.1128/IAI.00312-18
[potemberg]: https://journals.plos.org/plospathogens/article?id=10.1371/journal.ppat.1010621
[explorer]: https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0126070
