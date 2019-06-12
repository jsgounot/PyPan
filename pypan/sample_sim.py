# -*- coding: utf-8 -*-
# @Author: jsgounot
# @Date:   2019-02-11 19:40:34
# @Last modified by:   jsgounot
# @Last Modified time: 2019-06-12 15:12:31

import os, glob
import pandas as pd
pd.set_option('display.width', 1000)

bname = os.path.basename
dname = os.path.dirname

from Bio import SeqIO
from pypan import utils, processing
from pyblast import BCLine6, TMPFasta

def extract_rm(df, pid) :
    df = df[df["qseqid"] != df["sseqid"]]
    df = df[df["pident"] >= pid]

    def get_min_seq(line) :
        if line["qlength"] >= line["slength"] :
            return pd.Series({"minid" : line["sseqid"], "minlen" : line["slen"]})
        else :
            return pd.Series({"minid" : line["qseqid"], "minlen" : line["qlen"]})

    # uggly but do the work
    df = df.merge(df.apply(get_min_seq, axis=1), left_index=True, right_index=True)
    df = df[df["length"] > (df["minlen"] - 50)]

    #print (df["minid"].unique())
    #print (df["minid"].nunique())

    return df["minid"].unique()

def run_sample(fname, pid, ** kwargs) :
    BCLine6.mbdb(fname, dbtype="nucl")

    blast_opt = {"gapopen" : "5",  "gapextend" : "5", "penalty" : "-5",
    "reward" : "1", "evalue" : "2", "word_size" : "11", ** kwargs}

    df = BCLine6("blastn", query=fname, subject=fname, ** blast_opt).run()  
    torm = [] if df.empty else extract_rm(df, pid)

    fdata = SeqIO.parse(fname, "fasta")
    fdata = [seq for seq in fdata if seq.id not in torm]

    sample = bname(dname(fname))
    outfile = os.path.join(dname(fname), sample + ".nr.nosim.fa")
    with open(outfile, "w") as handle : SeqIO.write(fdata, handle, "fasta")
    return outfile

def run(species, fnames, ncore=1, ** kwargs) :
    """Remove similar sequence within the same fasta file to make the process
    faster afterwards.
        
    Arguments:
        species {[str]} -- [Species to use corresponding to the name on the configuration file]
        fnames {[list]} -- [Fasta files to clean]
        ** kwargs -- [blast parameters to adjust]
    
    Keyword Arguments:
        ncore {int} -- [Number of core to use] (default: {1})

    Returns:
        [list] -- [Cleaned fasta files path generated]
    """

    try : pid = float(utils.get_field("PROJECT", species, "SS_PID"))
    except utils.PyPanUtilsException : pid = 98

    args = [processing.FunArgs(run_sample, fname, pid, ** kwargs) for fname in fnames]
    return processing.mproc(args, ncore=ncore)