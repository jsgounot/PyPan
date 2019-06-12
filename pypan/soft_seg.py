# -*- coding: utf-8 -*-
# @Author: jsgounot
# @Date:   2019-01-22 15:33:21
# @Last modified by:   jsgounot
# @Last Modified time: 2019-06-12 15:15:19

import os, glob, subprocess
import pandas as pd
pd.set_option('display.width', 1000)

bname = os.path.basename
dname = os.path.dirname

from Bio import SeqIO

from pypan import prediction, processing
from pypan.utils import get_soft, get_field, PyPanUtilsException

def get_result(fname) :
    fname = os.path.join(dname(fname), "Seg", "sequences.clean.fasta")
    return fname

def write_fdata(fdata, outfile) :
    # seg does a segmentation fault if sequence length is higher than ~ 10kb
    # which is of course bad protein to begin with
    # we remove them in the initial file, and they will be removed at the end of the process

    fdata = [sequence for sequence in fdata if len(sequence) < 10000]
    SeqIO.write(fdata, outfile, "fasta")

def run_seg(fname, outdir) :
    # launch seg for a fasta file

    seg = get_soft("seg")
    outfile = os.path.join(outdir, "seg_results.fasta")
    cmd = "%s %s -a -l -n > %s" %(seg, fname, outfile)
    subprocess.check_call(cmd, shell=True)

    segdata = SeqIO.parse(outfile, "fasta")
    seglist = []

    for sequence in segdata :
        sname, region = sequence.id, sequence.seq
        gname = sname[:sname.index("(")]
        start, end = sname[sname.index("(")+1:-1].split("-")
        seglist.append({"gname" : gname, "start" : start, "end" : end, "len" : len(region)})

    outfile = outfile = os.path.join(outdir, "seg_results.tsv")
    df = pd.DataFrame(seglist, columns=["gname", "start", "end", "len"])
    df.to_csv(outfile, sep="\t")

    return df

def process_segresult(fasta, df, outdir) :
    # parse the seg result files and look for low complexity regions
    # compare the segresults with the original fasta file

    fdata = SeqIO.parse(fasta, "fasta")
    fdata = pd.Series({sequence.id : len(sequence) for sequence in fdata}).rename("LenSeq")

    sum_size = df.groupby("gname")["len"].sum()
    count = df.groupby("gname").size().rename("count")

    df = pd.concat([fdata, sum_size, count], axis=1, join="outer")
    df = df.fillna(0)
    df["prc"] = (df["len"] / df["LenSeq"]) * 100
    df["remaining"] = df["LenSeq"] - df["len"]

    outfile = outfile = os.path.join(outdir, "seg_results.remaining.tsv")
    df.to_csv(outfile, sep="\t")

    return df

def clean_fasta(species, fasta, df, outdir) :
    # clean the original fasta file based on the threshold

    try : threshold = float(get_field("PROJECT", species, "seg_threshold"))
    except PyPanUtilsException : threshold = 50

    to_keep = set(df[df["remaining"] >= threshold].index)
    fdata = SeqIO.parse(fasta, "fasta")
    fdata = (sequence for sequence in fdata if sequence.id in to_keep)

    outfile = outfile = os.path.join(outdir, "sequences.clean.fasta")
    SeqIO.write(fdata, outfile, "fasta")

def run_fname(species, fname) :
    print ("seg", fname)
    fdata = prediction.get_fdata(fname) 

    outdir = os.path.join(dname(fname), "Seg")
    if not os.path.isdir(outdir) : os.mkdir(outdir)
    fname = os.path.join(outdir, "sequences.fasta")
    write_fdata(fdata, fname)

    df = run_seg(fname, outdir)
    df = process_segresult(fname, df, outdir)
    clean_fasta(species, fname, df, outdir)
    return fname, df

def run(species, fnames, ncore=1) :
    """Run the SEG tool
    
    Remove sequence showing low complexity regions
    
    Arguments:
        species {[str]} -- [Species to use corresponding to the name on the configuration file]
        fnames {[list]} -- [Fasta files prior to prediction, this script will automaticaly retrieve results from prediction tools]
    
    Keyword Arguments:
        ncore {number} -- [Number of core to use] (default: {1})
    """

    fnames = glob.glob(fnames) if isinstance(fnames, str) else fnames
    args = [processing.FunArgs(run_fname, species, fname) for fname in fnames]
    processing.mproc(args, ncore=ncore)
