# -*- coding: utf-8 -*-
# @Author: jsgounot
# @Date:   2019-09-07 16:48:21
# @Last modified by:   jsgounot
# @Last Modified time: 2019-11-13 12:04:05

import fire
import os, glob, json, random
import itertools, collections
import pandas as pd
import processing

def _json2tsv(json_path, outfile=None, split=True, samples=[], value=1) :

    if not os.path.isfile(json_path) :
        print ("ERROR - File not found : %s" %(json_path))
        exit()

    with open(json_path) as f :
        jdata = json.load(f)

    if split :
        jdata = {sequence : {sample.split("::")[0] for sample in samples}
            for sequence, samples in jdata.items()}

    else :
        jdata = {sequence : set(samples) for sequence, samples in jdata.items()}

    all_samples = set.union(* jdata.values())
    print ("Found %i samples" %(len(all_samples)))

    jdata = [{"name" : sequence, ** {sample : value if sample in samples else 0 for sample in all_samples}}
        for sequence, samples in jdata.items()]

    jdata = pd.DataFrame(jdata).set_index("name")
    outfile = outfile or json_path[:-4] + "tsv"
    jdata.to_csv(outfile, sep="\t")

# ----------------------------------------------------------------------------------------------------------------------------

class TSVFile() :

    def __init__(self, tsv_file) :
        print ("Read tsv file : %s" %(tsv_file))
        self.df = pd.read_csv(tsv_file, sep="\t", index_col=0)
        self.path = tsv_file

        print ("%i samples" %(len(self.samples)))
        print ("%i genes" %(len(self.genes)))
        print ("-" * 20)

    @property
    def samples(self) :
        return set(self.df.columns)

    @property
    def genes(self) :
        return set(self.df.index)

    @staticmethod
    def fuse_tsv(tsvfiles, missing=0) :
        return pd.concat((tsvfile.df for tsvfile in tsvfiles)).fillna(missing).astype(int)

def _make_distribution(data, nsamples, sample_size) :
    results = []

    print ("Done with sample size %i" %(nsamples))    

    for count in range(0, sample_size) :
        samples = random.sample(list(data), nsamples)
        genes = [data[sample] for sample in samples]

        core = len(set.intersection(* genes))
        pan = len(set.union(* genes))
        results.append({"nsample" : nsamples, "core" : core, "pan" : pan})

    print ("Done with sample size %i" %(nsamples))
    return pd.DataFrame(results)

def _saturation(tsv_files, outfile, sample_size, force, ncore) :
    tsv_files = [TSVFile(tsvfile) for tsvfile in glob.glob(tsv_files)]

    for pair in itertools.combinations(tsv_files, 2) :
        t1, t2 = pair
        
        samples_diff = t1.samples ^ t2.samples
        if samples_diff :
            print ("WARNING : %i samples are found in one file but not the other" %(len(samples_diff)))
            print ("Files : %s - %s" %(t1.path, t2.path))
            print ("Sample diff : %s" %(str(samples_diff)))
            if not force :
                print ("If you want to override this, use force=True")
                exit()

        genes_overlapp = t1.genes & t2.genes
        if genes_overlapp :
            print ("WARNING : %i genes were found in both files" %(len(genes_overlapp)))
            print ("Files : %s - %s" %(t1.path, t2.path))
            print ("Genes : %s" %(str(genes_overlapp)))
            if not force :
                print ("If you want to override this, use force=True")
                exit()

    df = TSVFile.fuse_tsv(tsv_files)
    all_samples = set(df.columns)
    total_count = len(df.index)
    df = df[df != 1].dropna(how='all', axis=0)
    coregenome_size = total_count - len(df.index)

    print ("TOTAL MATRIX")
    print ("%i samples" %(len(all_samples)))
    print ("%i genes (core : %i, variable : %i)" %(total_count, coregenome_size, len(df.index)))
    print ("Launch distribution (sample size = %i)" %(sample_size))

    genes = collections.defaultdict(set)
    for gene, row in df.iterrows() :
        samples = row[row != 0].index
        for sample in samples :
            genes[sample].add(gene)

    args = [processing.FunArgs(_make_distribution, genes, nsamples, sample_size)
        for nsamples in range(1, len(all_samples))]

    res = processing.mproc(args, ncore)
    res = pd.concat(res)

    res["core"] = res["core"] + coregenome_size
    res["pan"] = res["pan"] + coregenome_size

    res.to_csv(outfile, sep="\t")
    print ("Saved at %s" %(outfile))

# ----------------------------------------------------------------------------------------------------------------------------

class Manager(object) :

    def json2tsv(self, json_path, outfile=None, split=True, samples=[], value=1) :
        """     
        Transform json file produced by pypan to a tsv file
        
        Arguments:
            json_path {str} -- Path to your json file
        
        Keyword Arguments:
            outfile {[type]} -- Path to your outfile file, if None the file is created using the same name of the input (default={None})
            split {bool} -- Split the sample name of the json file, specific to pypan, set to false if you want to use this function with a json file not produced by pypan (default: {True})
            samples {list} -- Additionnal sample names to include within the tsv file (default: {[]})
            value {int} -- Default value for gene found (default : {1})
        """

        _json2tsv(json_path, outfile, split, samples, value)

    def saturation(self, tsv_files, outfile, sample_size=5000, force=False, ncore=1) :
        """        
        Make the saturation distribution for the core and the pangenome based on matrix files.
        This function can be slow depending of your matrix size and the sample size (corresponding to the number n of random genomes combinations)
        It can be usefull to first try with a slow number of sample size (like 10) before using a higher number

        Arguments:
            tsv_files {[str]} -- Path of your tsv files, use quote if more than one
            outfile {[type]} -- Path of your result file (a tsv file)
        
        Keyword Arguments:
            sample_size {number} -- Sample size to use, corresponding to the number n of random genomes combinations (default: {5000})
            force {bool} -- Ignore error with TSV files checking (default: {False})
            ncore {number} -- Number of core to use (default: {1})
        """

        _saturation(tsv_files, outfile, sample_size, force, ncore)

if __name__ == '__main__' :
    manager = Manager()
    fire.Fire(Manager)