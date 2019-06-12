# -*- coding: utf-8 -*-
# @Author: jsgounot
# @Date:   2019-01-22 16:20:12
# @Last modified by:   jsgounot
# @Last Modified time: 2019-06-12 15:18:39

import os, glob, json
from collections import Counter, defaultdict
import pandas as pd
pd.set_option('display.width', 1500)

bname = os.path.basename
dname = os.path.dirname

from Bio import SeqIO
from Bio.Seq import Seq

import networkx as nx

from pypan import soft_seg, prediction
from pypan.utils import get_outdir, get_field, PyPanUtilsException
from pyblast import Local2Global, BCLine6, BCLineHSPFuse, TMPFasta

def graph_groups(df) :
    # here we don't care about which on is the best, we just take the center of the graph
    # this function is the one which is the most similar to the 1002 one

    dseq = {}
    graphs = set()

    for smax, smin in zip(* (df["qseqid"], df["sseqid"])) :

        g1 = dseq.get(smax, None)
        g2 = dseq.get(smin, None)
        update_seq = [smin, smax]

        if g1 is None and g2 is None :
            graph = nx.Graph()
            graph.add_edge(smax, smin)
            graphs.add(graph)

        elif g1 is g2 :
            g1.add_edge(smax, smin)
            continue

        elif g1 is None :
            g2.add_edge(smax, smin)
            graph = g2

        elif g2 is None :
            g1.add_edge(smax, smin)         
            graph = g1

        else :
            g1.add_edge(smax, smin)
            graph = nx.compose(g1, g2)
            graphs.remove(g1)
            graphs.remove(g2)
            graphs.add(graph)

            update_seq.extend(list(g1.nodes()))
            update_seq.extend(list(g2.nodes()))

        for sname in update_seq :
            dseq[sname] = graph

    return {nx.center(graph)[0] : set(graph.nodes) for graph in graphs}

def localGlobal(species, fname, ncore, ** kwargs) :
    blast_opt = {"evalue" : "2", "word_size" : "3", ** kwargs.pop("bkwargs", {})}
    l2g = Local2Global(str(fname), str(fname))
    df = l2g.run("blastp", ncore=ncore, chunksize=30, bkwargs=blast_opt, match_fun=None, ** kwargs)

    df = df[["qseqid", "sseqid", "qpos", "spos", "gSCORE", "gQPID", "gSPID"]]
    df = df[df["qseqid"] != df["sseqid"]]

    try : threshold = float(get_field("PROJECT", species, "L2G_sim"))
    except PyPanUtilsException : threshold = 90

    df = df[(df["gQPID"] >= threshold) | (df["gSPID"] >= threshold)]
    df = df.drop_duplicates(["qseqid", "sseqid"])
    return df

def unirefStyle(species, fname, ncore, def_uni_ide=90, def_uni_cov=80, ** kwargs) :
    blast_opt = {"evalue" : "2", "word_size" : "3", "query" : str(fname), "db" : str(fname), ** kwargs}

    bc = BCLine6("blastp", ** blast_opt)
    df = bc.run(ncore=ncore, chunksize=50)

    # By default we use uniref90 values
    # At least 90 sim and 80 coverage

    try : uni_ide = float(get_field("PROJECT", species, "uni_ide"))
    except PyPanUtilsException : uni_ide = def_uni_ide

    try : uni_cov = float(get_field("PROJECT", species, "uni_cov"))
    except PyPanUtilsException : uni_cov = def_uni_cov

    df = df[df["pident"] >= uni_ide]
    df = df[(df["qcov"] >= uni_cov) | (df["scov"] >= uni_cov)]

    df = df.drop_duplicates(["qseqid", "sseqid"])
    return df

def HSPFuse(species, fname, ncore, simprc=98, ** kwargs) :
    blast_opt = {"evalue" : "2", "word_size" : "3", "query" : str(fname), "db" : str(fname), ** kwargs}

    bc = BCLineHSPFuse("blastp", ** blast_opt)
    df = bc.run(ncore=ncore, chunksize=50)

    try : simprc = float(get_field("PROJECT", species, "HSPF_SIM"))
    except PyPanUtilsException : pass

    df = df[(df["qpos_prc"] >= simprc) | (df["spos_prc"] >= simprc)]
    return df

def compare_sequences(species, fdata, method=None, outdir=None, ncore=1, ** kwargs) :
    outdir = outdir or get_outdir(species, "PanGenome")
    db_outdir = os.path.join(outdir, "DBSMerge")
    if not os.path.isdir(db_outdir) : os.mkdir(db_outdir)

    fname = TMPFasta(fdata, dir=db_outdir, delete=False)
    BCLine6.mbdb(str(fname), dbtype="prot")

    method = method or get_field("PROJECT", species, "sm_method")
    methods = {"L2G" : localGlobal, "UNI" : unirefStyle, "HSPFuse" : HSPFuse}
    if method not in methods : raise Exception("Method not found : %s" %(method))
    df = methods[method](species, fname, ncore, ** kwargs)
    
    outfile = os.path.join(outdir, "cs.raw.tsv")
    df.to_csv(outfile, sep="\t")

    # Make groups and fill it with non used sequence (sequence without a similarity)
    groups = graph_groups(df)
    notfound = {sequence.id for sequence in fdata} - set.union(* groups.values())
    groups.update({sequence : set((sequence, )) for sequence in notfound})

    # count sample
    cs_previous = pd.Series(Counter(sequence.id.split("::")[0] for sequence in fdata)).rename("Previous")
    cs_groups = pd.Series(Counter(sample for group in groups.values() for sample in {seq.split("::")[0] for seq in group})).rename("Groups")
    cs_rep = pd.Series(Counter(seq.split("::")[0] for seq in groups)).rename("Representative")
    
    df = pd.concat([cs_previous, cs_groups, cs_rep], axis=1).fillna(0).astype(int)
    df["%Representative"] = df["Representative"] * 100 / df["Representative"].sum()

    # Explanation of the table and why we dont find Groups + Representative = Previous
    # -> If 2 sequences or more from the same strains are found in one group, the sample will be counted only once

    # We save everything
    outfile = os.path.join(outdir, "groups.json")
    groups = {name : list(names) for name, names in groups.items()}
    with open(outfile, "w") as f :
        json.dump(groups, f, indent=4)

    outfile = os.path.join(outdir, "countseq.tsv")
    df.to_csv(outfile, sep="\t")

    fdata = [seqrecord for seqrecord in fdata if seqrecord.id in set(groups)]
    outfile = os.path.join(outdir, "sequences.clean.prot.fa")
    SeqIO.write(fdata, outfile, "fasta")

    return groups

class SCounter() :

    def __init__(self) :
        self.count = 0

    def __call__(self) :
        self.count += 1
        return self.count

def merge_fnames(sfnames, minsize, maxsize) :

    names = []
    fdata = []

    # we gather results
    for sample, fnames in sfnames.items() :
        sc = SCounter()
        fnames = [fnames] if isinstance(fnames, str) else fnames
        for fname in fnames :
            for sequence in SeqIO.parse(fname, "fasta") :

                if not (minsize <= len(sequence.seq) <= maxsize) : continue

                new_name = "Seq%i" %(sc())
                new_name = "%s::%s" %(sample, new_name)

                names.append({"sample" : sample, "name" : sequence.id, "new_name" : new_name})
                sequence.id = new_name
                sequence.name = sequence.id
                sequence.description = sequence.id
                sequence.seq = Seq(str(sequence.seq).replace("?", "N"))
                fdata.append(sequence)

    names = pd.DataFrame(names)
    return fdata, names

def write_nucl(species, fnames, names, used) :
    
    # we map new names
    newnames = defaultdict(dict)
    for idx, row in names.iterrows() :
        sample, name, newname = row["sample"], row["name"], row["new_name"]
        newnames[sample][name] = newname

    fdata = []
    for fname in fnames :
        sample = bname(dname(fname))
        for seqrecord in prediction.get_fdata(fname, prot=False) :
            newname = newnames[sample].get(seqrecord.description)
            if newname is None : continue # removed by seg
            if newname not in used : continue
            seqrecord.id = seqrecord.name = seqrecord.description = newname
            fdata.append(seqrecord)

    outfile = os.path.join(get_outdir(species, "PanGenome"), "sequences.clean.nucl.fa")
    SeqIO.write(fdata, outfile, "fasta")

def run(species, fnames, ncore=1, ** kwargs) :
    """[summary]
    
    [description]
    
    Arguments:
        species {[str]} -- [Species to use corresponding to the name on the configuration file]
        fnames {[list]} -- [Fasta files prior to prediction, this script will automaticaly retrieve results from seg results]
        ** kwargs -- [Kwargs transfered to each method accordingly. Most of the time it is blast parameters but can also be additionnal parameters in case of L2G. See Local2Global.run in this case (PyBlast lib)]
    
    Keyword Arguments:
        ncore {number} -- [Number of core to use] (default: {1})
    """

    sfnames = {bname(dname(fname)) : soft_seg.get_result(fname) for fname in fnames}
    print ("sm", sfnames)

    minsize = float(get_field("PROJECT", species, "prot_minsize"))
    maxsize = float(get_field("PROJECT", species, "prot_maxsize"))
    fdata, names = merge_fnames(sfnames, minsize, maxsize)

    # We save names
    outfile = os.path.join(get_outdir(species, "PanGenome"), "names.tsv")
    names.to_csv(outfile, sep="\t")

    # We launch the process
    groups = compare_sequences(species, fdata, ncore=ncore, ** kwargs)

    # We save nucleic data
    write_nucl(species, fnames, names, set(groups))
