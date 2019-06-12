---
title: PyPan
tags: [Github, Notebooks, Python]
created: '2019-01-28T09:54:16.863Z'
modified: '2019-06-12T13:36:23.416Z'
---

# PyPan

PyPan is a python library designed for pangenome construction, using local and global alignment, Augustus and SNAP for gene predictions, and a graph based method for genes clustering. PyPan takes assemblies and a reference sequence as inputs and automaticaly removed reference segments, predict genes from non-reference segments (NR segs) and compare new genes within a population. PyPan is designed to take advantage of multi-core system, which is a requirement for most pangenome construction.

PyPan is mainly composed of 4 steps :
* **Step 1** : Assemblies comparaison against the reference sequence to retrieve solely non-reference (species) segments
* **Step 2** : Genes prediction using both SNAP and Augustus
* **Step 3** : Genes cleaning based on sequence length and sequence complexity (Seg)
* **Step 4** : Genes comparison and graph clustering using Blast

## Installation

Download / clone the github project and follow instructions below.

## Requirements

### Softwares

* [Blast+](https://www.ncbi.nlm.nih.gov/guide/howto/run-blast-local/)
* [Seg](ftp://ftp.ncbi.nih.gov/pub/seg/seg/)
* [SNAP](https://github.com/KorfLab/SNAP)
* [Augustus](http://augustus.gobics.de/)

### Python libraries

For now, PyPan has only been tested using python 3.6. I highly recommend to use the same version.

* [Pandas](https://pandas.pydata.org/)
* [Biopython](https://biopython.org/)
* [Networkx](https://networkx.github.io/) 
* [PyBlast](https://github.com/jsgounot/PyBlast)

## Pre configuration and methods

### Input files

PyPan needs at least an assembly for each strain (fasta format) and a reference sequence (fasta format).

### Identity score (Step 1)

By default, PyPan takes a single identity score for all samples, corresponding to the identity threshold for which a segment is considered to come from the reference species based on a similarity value. Since this value can highly diverged within a single population, a sample specific threshold can also be given, corresponding to a 2 columns, no header, tabulated file (see example). Note that the sample threshold add up to the species threshold. For example, if we have a sample with a value of 2, and a species threshold of 5, a final treshold of 100 - (2 + 5) will be used for the sample. 

### Prediction software (Step 2)

Manual training of both Augustus and SNAP are requiered before running PyPan. Look at the specific documentation for these tools. For both of them, several pre-trained datasets are available for the most known species. 

### Samples similarity methods (Step 4)

For now, three methods are available to cluster genes with similar sequences : 

* **UNI** : Similar to Uniprot method, based on the HSP query coverage and the HSP identity percentage (blastp)
* **L2G** : A method combining global alignment (blastp), extracting best hits, and preforming global alignment afterward
* **HSPFuse** : A method which combine HSPs from the same query / subject pair to provide a fast similarity score (blastp)

UNI method is by far the fastest method and can even provide the best result when comparing results from both methods without any extra optimization. If you don't know which method to use for your first try, use UNI.

**About UNI**

UNI method is kind of similar to the UniRef method and select pair based on blastp results. For one pair to be selected, it has to have a coverage > `uni_cov` (default : 80) and an identity percentage > `uni_ide` (default : 90). These values can be changed in the configuration file for each project.

**About L2G and global alignment**

For global alignment, biopython [`pairwise2`](http://biopython.org/DIST/docs/api/Bio.pairwise2-module.html) method is used by default with the function `pairwise2.align.globalxx`. The default PyPan algorithm grab the best alignment for each pair with the highest query identity percentage (QIP : corresponding to the number of well aligned bases, divided by the query size). Only global alignment showing a global alignment similarity > `L2G_sim` are retained. I highly encourage to check and use other parameters, since default configuration can leads to sporious alignments in which most aligned bases are surrounded by gaps and therefore to high QIP, especialy when a small query is compared to a long subject sequence.

**About HSPFuse**

This method combines all HSPs from the same query / subject pair to provide a fast similarity score. While this method can provide score close to global alignment without predifined parameters, it can also lead to different results, especially in case of repetitive sequences. See *Merging HSP* in [PyBlast documentation](https://github.com/jsgounot/PyBlast).

**Clustering pairs**

Once pairs have been determined, clustering is based on a network analysis. For this, networks are constructed in which each node corresponds to a sequence and each bond correspond to a similarity found between two sequences. Once done, the center of each network (or graph) is determined with the [`networkx.center`](https://networkx.github.io/documentation/networkx-1.10/reference/generated/networkx.algorithms.distance_measures.center.html) function and this sequence correspond to the representive of the group. These networks correspond to the `groups.json` file in the result directory. Note that several sequences from the same sample can be found in the same group., statistic of the number of representatives and the number of groups in which each sample has been found are shown in the `countseq.tsv` result file. I recommand to look at this file before any analysis, since theses numbers can potentially reveal contaminations.

## PyPan configuration file

PyPan uses a json formated file as a configuration file. An example of this file is given in the example folder. Here is a description of the different fields :

- PROJECTS : Informations links to each species
  * genomes : Assemblies path (str, fastas)
  * reference : Reference path (str, fasta)
  * outdir : Outdir path (str)
  * exception : List of samples to not use, names correspond to new formated names (see below) (list)
  * pid_file : Tabulated file with sample specific identity threshold (str)
  * pid_value : Single threshold of identity threshold for the whole species (see above) (int, float)
  * augustus : Augustus species name to use (str)
  * snap : Path of the hmm file to use (str)
  * sm_method : Sample similarity method to use (str : L2G or UNI)
  * seg_threshold : Seg threshold corresponding to the minimum low complexity size (% protein size) for a protein to be retained (float, default : 50)
  * uni_cov : Coverage value to use for the UNI method (float, default : 80)
  * uni_ide : Identity percentage to use for the UNI method (float, default : 90)
  * L2G_sim : Global alignment similarity score to use for the L2G method (float, default : 90)
  * HSPF_sim : Semi global alignment similarity score to use for the HSPFuse method (float, default : 98)

- SOFTPATH : Softwares executable path
  * snap (str)
  * augustus (str)
  * seg (str)
 
Since blast is managed by BioPython, it is assumed that blast executables are in the PATH. See BioPython documentation is case of errors.

## Running PyPan

### Using the command line with an example

While PyPan pipeline can be modified using a simple python script, a basic command line implementation is proposed with the ```basic_launcher``` script. When all software are installed, you can test your PyPan installation using this example :

```bash

cd ./path/to/pypan
tar -zxvf example.tar.gz
python basic_launcher.py Example/configuration.json sace --ncore 4

```

Once it's over and if no error is printed, all results should be available at ```./Example/OutDir/PanGenome/```

### Output files 

In **bold**, files which could be particulary usefull for pangenome analysis.

* **countseq.tsv** : Number of sequences for each sample found at relative steps before and during species merging
* *cs.raw.tsv* : Result of the blast analysis, used for merging genes by the graph based method
* **groups.json** : Groups of genes with keys corresponding to representatives sequences and therefore found in the final pangenomes and values corresponding to similar sequences. This file give you the link if a sequence in your pangenome is found in one or multiple samples.
* *names.tsv* : A transition file which make the link between the name given for each sample in previous step and the name used beforme merging
* **sequences.clean.*.fa** : Representatives proteic and nucleic sequences

### Tweaking options using python

The default pipeline is composed of only 5 functions calls corresponding to each step :

```python3

import os
import pypan

# we setup the configuration file
pypan.utils.set_config(config_file_path)

def pipeline() :

    ncore = 10 # number of CPU
    species = "Sace"
    
    # sample name based on path
    funname = lambda fname : os.path.basename(fname).splitext()[0] 
    
    fnames = pypan.nrfinder.nrclean(species, funname=funname, ncore=ncore)
    pypan.prediction.launch_fnames(species, fnames, ncore=ncore)
    pypan.soft_seg.run(species, fnames, ncore=ncore)
    pypan.sample_merge.run(species, fnames, ncore=ncore)

```

The parameters used in each function can be retrieved easily using the docstring :

```python

import pypan
help(pypan.nrfinder.nrclean)

```

However, most of the time important parameters are the ones specified in the configuration file. Nevertheless, you can still tuned your pipeline with different blast parameters for example :

```python

import os
import pypan

# we setup the configuration file
pypan.utils.set_config(config_file_path)

def pipeline() :

    ncore = 10 # number of CPU
    species = "Sace"
    
    # sample name based on path
    funname = lambda fname : os.path.basename(fname).splitext()[0] 
    
    fnames = pypan.nrfinder.nrclean(species, funname=funname, ncore=ncore, word_size=10)
    pypan.prediction.launch_fnames(species, fnames, ncore=ncore)
    pypan.soft_seg.run(species, fnames, ncore=ncore)
    pypan.sample_merge.run(species, fnames, ncore=ncore, evalue=3)

```

## Default blast options

| Module       	| Blast kind | Function 	| evalue 	| gapextend 	| gapopen 	| penalty 	| reward 	| word_size 	|
|--------------	|----------- |----------	|--------	|-----------	|---------	|---------	|--------	|-----------	|
| nrfinder     	| Blastn     | nrclean  	| 2      	| 5         	| 5       	| -5      	| 1      	| 11        	|
| sample_sim   	| Blastn     | run      	| 2      	| 5         	| 5       	| -5      	| 1      	| 11        	|
| sample_merge 	| Blastp     | UNI      	| 2      	|           	|         	|         	|        	| 3         	|
| sample_merge 	| Blastp     | L2G      	| 2      	|           	|         	|         	|        	| 3         	|
| sample_merge 	| Blastp     | HSPFuse  	| 2      	|           	|         	|         	|        	| 3         	|

## Some notes about multiprocessing

PyPan is designed for cluster analysis. This way, multiprocessing can be used in almost all parts of the pipeline. Most of the multicore use is done when blast and global alignment are done and are carried by the PyBlast library. [Look at the specific github for more informations](https://github.com/jsgounot/PyBlast).
