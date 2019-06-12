# -*- coding: utf-8 -*-
# @Author: jsgounot
# @Date:   2019-01-22 13:47:40
# @Last modified by:   jsgounot
# @Last Modified time: 2019-06-12 15:12:07

import os, shutil
import pandas as pd

bname = os.path.basename
dname = os.path.dirname

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from pypan import utils, processing
from pypan.predictors.augustus import Augustus
from pypan.predictors.snap import SNAP

class PredException(Exception) :
    pass

def get_fdata(fname, prot=True) :
    kind = "prot" if prot else "nucl"
    
    augustus_fname = os.path.join(dname(fname), "CDSPrediction", bname(fname) + ".augustus.%s.fasta" %(kind))
    fdata = list(SeqIO.parse(augustus_fname, "fasta"))

    snap_fname = os.path.join(dname(fname), "CDSPrediction", bname(fname) + ".snap.%s.fasta" %(kind))
    fdata += list(SeqIO.parse(snap_fname, "fasta"))

    return fdata

def launch_augustus(fname, augspe, fdata) :

    # Launch Augustus
    outfile = fname + ".augustus.gb"
    aug = Augustus(fname, augspe, outfile, run_soft=True)

    proteins = [SeqRecord(id=gene.name, seq=Seq(gene.sequences["protein"]), description="") for gene in aug.get_features()]
    outfile_prot = fname + ".augustus.prot.fasta"
    SeqIO.write(proteins, outfile_prot, "fasta")

    nucleic = [SeqRecord(id=gene.name, seq=Seq(gene.childs[0].sequences["CDS"]), description="") for gene in aug.get_features()]
    outfile_nucl = fname + ".augustus.nucl.fasta"
    SeqIO.write(nucleic, outfile_nucl, "fasta")

    infos = [{"gene" : gene.name, "prot_length" : len(gene.sequences["protein"]), "origin" : gene.chromosome,
    "cds_length" : len(gene.childs[0].sequences["CDS"])} for gene in aug.get_features()]
    df = pd.DataFrame(infos)
    
    if not df.empty :
        df["lenOrigin"] = df["origin"].map({sequence.id : len(sequence.seq) for sequence in fdata})

    outfile_info = fname + ".augustus.prot.info"
    df.to_csv(outfile_info, sep="\t")

def launch_snap(fname, sdir) :
    snap = SNAP(fname, sdir)
    snap.run()

def launch_fname(species, fname) :
    
    # Directory
    outdir = os.path.join(dname(fname), "CDSPrediction")
    if not os.path.isdir(outdir) : os.mkdir(outdir)

    # Augustus : copy fasta and read it
    augustus_species = augspe = utils.get_field("PROJECT", species, "augustus")
    nfname = os.path.join(outdir, bname(fname))
    shutil.copyfile(fname, nfname)
    fdata = SeqIO.parse(nfname, "fasta")

    launch_augustus(nfname, augspe, fdata)
    launch_snap(nfname, species)

def launch_fnames(species, fnames, ncore=1) :
    """Launch the CDS prediction tools
       
    Arguments:
        species {[str]} -- [Species to use corresponding to the name on the configuration file]
        fnames {[list]} -- [Fasta files to clean]
    
    Keyword Arguments:
        ncore {number} -- [Number of core to use] (default: {1})
    """

    fnames = glob.glob(fnames) if isinstance(fnames, str) else fnames
    args = [processing.FunArgs(launch_fname, species, fname) for fname in fnames]
    processing.mproc(args, ncore=ncore)