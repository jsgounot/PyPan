# -*- coding: utf-8 -*-
# @Author: jsgounot
# @Date:   2019-01-22 10:44:06
# @Last modified by:   jsgounot
# @Last Modified time: 2019-06-12 15:09:35

import os, glob, shutil
import pandas as pd
pd.set_option('display.width', 1000)

bname = os.path.basename
dname = os.path.dirname

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from pypan import utils, processing
from pyblast import BCLine6, TMPFasta
from pypan.pysegs import Segment, SegList

class NRException(Exception) :
    pass

def get_pid(species, samples) :
    try : pidfile = utils.get_field("PROJECT", species, "pid_file")
    except utils.PyPanUtilsException : pidfile = None
    pidvalue = utils.get_field("PROJECT", species, "pid_value")

    if not os.path.isfile(pidfile) :
        error_message = "PID File not found : %s" %(pidfile)
        raise utils.PyPanUtilsException(error_message)

    if pidfile is None or not os.path.isfile(pidfile) : 
        return {sample : 100 - float(pidvalue) for sample in samples}

    df = pd.read_csv(pidfile, sep="\t", names=["sample", "pid"])
    df = df.set_index("sample")["pid"].to_dict()

    samples = {sample : None for sample in samples}

    for sample in samples :
        if sample not in df :
            raise NRException("Sample not found in the pid file : %s" %(sample))
        samples[sample] = 100 - (float(pidvalue) + float(df[sample]))

    return samples

def get_samples(species, funname=None) :
    funname = funname or (lambda fname : bname(fname))
    samples = utils.get_field("PROJECT", species, "genomes")
    exceptions = set(utils.get_field("PROJECT", species, "exceptions"))

    samples = glob.glob(samples)
    samples = {funname(fname) : fname for fname in samples}
    samples = {sample : fname for sample, fname in samples.items() if sample not in exceptions}

    return samples

def make_ref(species, ref) :
    outdir = utils.get_outdir(species, "Reference")
    newref = os.path.join(outdir, bname(ref))
    shutil.copyfile(ref, newref)
    return newref

def uniform_clean(df, fdata, minlength, mindist=100) :
    fsegs = {sname : Segment(0, len(seq)) for sname, seq in fdata.items()}

    for qseqid, sdf in df.groupby("qseqid") :
        seglist = SegList(zip(sdf["qstart"], sdf["qend"]))

        # we remove segments with a size < mindist
        # this way we join nr sequences which have a distance shorter than mindist 
        seglist = SegList([element for element in seglist if len(element) >= mindist])

        tokeep = fsegs[qseqid] - seglist

        for segment in tokeep :
            if len(segment) < minlength : continue
            start, end = segment.start, segment.stop
            
            idname = qseqid + "_%i_%i" %(start, end)
            sequence = fdata[qseqid][start:end]
            sequence.id = idname

            if set(sequence.seq) == set("N") : continue
            fdata[idname] = sequence

        fdata.pop(qseqid)

    return fdata

def side_step(fdata, ref, sidelist, cutsize, sidepid, sidecov, bkwargs) :

    nfdata = []
    
    for qseqid, segments in sidelist.items() :
        sequence = fdata[qseqid]
        for segment in segments :
            subseq = sequence[segment.start : segment.stop]
            for i in range(0, len(segment), cutsize) :
                cutseq = subseq[i:i+cutsize]
                if len(cutseq) == 0 or set(cutseq) == set("N") : continue

                start = segment.start + i
                name = qseqid + ":%i:%i" %(start, start + len(cutseq))

                cutseq.id = cutseq.name = cutseq.description = name
                nfdata.append(cutseq)

    fname = TMPFasta(nfdata)
    blast_opt = {"gapopen" : "5",  "gapextend" : "5", "penalty" : "-5",
    "reward" : "1", "evalue" : "2", "word_size" : "11", ** bkwargs}
    
    df = BCLine6("blastn", query=str(fname), db=ref, ** blast_opt).run()
    df = df[(df["pident"] >= sidepid) & (df["qcov"] >= sidecov)]

    sidelist = {}
    for name in df["qseqid"] :
        name = name.split(":")
        contig, start, end = ":".join(name[:-2]), name[-2], name[-1]
        sidelist.setdefault(contig, SegList()).append(Segment(start, end))

    return sidelist

def complexe_clean(df, fdata, minlength, ref, sidepid, sidesize=1000, cutsize=250, sidecov=75, mindist=100, bkwargs={}) :

    # Clean long HSP with a side step process (another blast where long HSP are cut based on cutsize and sidesize args)

    # minlength = Minimum size of a NR segment to be used
    # sidesize = Minimum HSP size for an HSP to be subject to another blast
    # cutsize = Slice size for selected HSP
    # sidecov = %Coverage of the sliced HSP to be selected as an reference segment
    # mindist = Minimum distance of the NR segment

    # Note for myself : 2 changes compared to V0
    # 1. We discard sequence with high N prc
    # 2. We cut merged HSP and not each HSP individualy

    mainlist = {}
    sidelist = {}

    for qseqid, sdf in df.groupby("qseqid") :
        
        # One difference here compared to V0 is that HSP are first merged
        # and after cut, while the sidesize was search for individual HSP previously

        slist = SegList(zip(sdf["qstart"], sdf["qend"]))

        for segment in slist :
            dic = sidelist if len(segment) >= sidesize else mainlist
            dic.setdefault(qseqid, []).append(segment)  

    if sidelist : sidelist = side_step(fdata, ref, sidelist, cutsize, sidepid, sidecov, bkwargs)
    mainlist = {qseqid : SegList(segments) for qseqid, segments in mainlist.items()}

    for contig, seglist in sidelist.items() :
        mainlist.setdefault(contig, SegList([])).extend(seglist)

    # we remove segments with a size < mindist
    # this way we join nr sequences which have a distance shorter than mindist 
    mainlist = {contig : SegList([element for element in elements if len(element) >= mindist])
        for contig, elements in mainlist.items()}

    fsegs = {sname : Segment(0, len(seq)) for sname, seq in fdata.items()}

    for qseqid, seglist in mainlist.items() :

        tokeep = fsegs[qseqid] - seglist

        for segment in tokeep :
            if len(segment) < minlength : continue
            start, end = segment.start, segment.stop
            
            idname = qseqid + "_%i_%i" %(start, end)
            sequence = fdata[qseqid][start:end]
            sequence.id = idname

            if set(sequence.seq) == set("N") : continue
            fdata[idname] = sequence

        fdata.pop(qseqid)

    return fdata

def check_round(fname, ref, species, sample, seglength=50, ** bkwargs) :

    fdata = SeqIO.to_dict(SeqIO.parse(str(fname), "fasta"))
    pid = get_pid(species, [sample])[sample]

    # Blast
    blast_opt = {"gapopen" : "5",  "gapextend" : "5", "penalty" : "-5",
    "reward" : "1", "evalue" : "2", "word_size" : "11", ** bkwargs}
    df = BCLine6("blastn", query=str(fname), db=ref, ** blast_opt).run()

    df = df[(df["pident"] >= pid) & (df["length"] >= seglength)]
    df["qstart"] = df["qstart"] - 1

    intervals = {contig : SegList(zip(sdf["qstart"], sdf["qend"])) for contig, sdf in df.groupby("qseqid")}

    return df, intervals

def nrclean_sample(species, sample, fname, ref, pid, minlength=200, seglength=50, 
        toclean=True, nround=None, complexe=True, outdir=None, ** bkwargs) :

    print ("Sample : %s - NRound : %i" %(sample, nround))
    outdir = outdir or utils.get_outdir(species, "Samples", sample)

    if nround == 0 :
        outfile = os.path.join(outdir, sample + ".nr.fa")
        shutil.copyfile(str(fname), outfile)
        return outfile

    if toclean :
        fdata = BCLine6.clean_fasta(str(fname), astmp=False)
        fdata = [sequence for sequence in fdata if len(sequence) >= minlength]
        fname = TMPFasta(fdata)
        fdata = {sequence.id : sequence for sequence in fdata}

    else :
        fdata = SeqIO.to_dict(SeqIO.parse(str(fname), "fasta"))

    # Blast
    blast_opt = {"gapopen" : "5",  "gapextend" : "5", "penalty" : "-5",
    "reward" : "1", "evalue" : "2", "word_size" : "11", ** bkwargs}
    df = BCLine6("blastn", query=str(fname), db=ref, ** blast_opt).run()

    df = df[(df["pident"] >= pid) & (df["length"] >= seglength)]
    df["qstart"] = df["qstart"] - 1

    # if we don't have DF, it means that the file is cleaned
    if df.empty :
        outfile = os.path.join(outdir, sample + ".nr.fa")
        shutil.copyfile(str(fname), outfile)
        return outfile

    # we remove reference segments
    args = (df, fdata, minlength)
    fdata = complexe_clean(* args, ref, pid, bkwargs=bkwargs) if complexe else uniform_clean(* args)

    # fname = TMPFasta(sorted(fdata.values(), key = lambda x : int(x.id.split("_")[0][3:])), wodesc=True)
    fname = TMPFasta(fdata.values(), wodesc=True)
    nround = nround - 1 if nround else None

    return nrclean_sample(species, sample, fname, ref, pid, minlength, seglength, 
        toclean=False, nround=nround, complexe=complexe, outdir=outdir, ** bkwargs)

def nrclean(species, funname=None, cpref=True, ncore=1, ** kwargs) :
    """Remove reference segments from fasta files compared ot the reference. Most of the options here are from the configuration file.
    
    Arguments:
        species {[str]} -- [Species to use corresponding to the name on the configuration file]
        ** kwargs -- [blast parameters to adjust]
    
    Keyword Arguments:
        funname -- [callable which will transform the samples path to a name, by default basenames are used] (default: {None})
        cpref {bool} -- [Copy the reference sequence to the outdir directory] (default: {True})
        ncore {int} -- [Number of core to use] (default: {1})

    Returns:
        [list] -- [Cleaned fasta files path generated]

    """


    # Create outdir
    utils.get_outdir(species, "Samples")

    # Copy reference
    ref = utils.get_field("PROJECT", species, "reference")
    if not os.path.isfile(ref) : raise NRException("Reference file not found : %s" %(ref))
    ref = make_ref(species, ref) if cpref else ref
    BCLine6.mbdb(ref, dbtype = "nucl")

    # grab samples and pid
    samples = get_samples(species, funname)
    pid = get_pid(species, samples)
    
    # Show samples
    for sample, fname in samples.items() :
        print (sample, fname, pid[sample])
    print ("Total : %i samples" %(len(samples)))

    # NRounds
    try : nround = int(utils.get_field("PROJECT", species, "nrounds"))
    except utils.PyPanUtilsException : nround = None

    # Multiprocess
    args = [processing.FunArgs(nrclean_sample, species, sample, fname, ref, pid[sample], nround=nround, ** kwargs) for sample, fname in samples.items()]
    args = sorted(args, key = lambda arg : arg.args[1])
    return processing.mproc(args, ncore=ncore)