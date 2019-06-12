# -*- coding: utf-8 -*-
# @Author: jsgounot
# @Date:   2019-01-22 13:42:03
# @Last modified by:   jsgounot
# @Last Modified time: 2019-01-22 15:29:56

import os, subprocess
from Bio import SeqIO
from pypan.utils import get_soft
from pypan.predictors.features import Feature
from pypan.pysegs import SegList

class Augustus() :

    def __init__(self, fasta, species_augustus, outfile, run_soft=True) :
        self.fasta = fasta
        self.outfile = outfile
        if run_soft : self.run(species_augustus)

    @property
    def path(self):
        return get_soft("augustus")
    
    def run(self, species_augustus) :
        cmd = "%s --species=%s %s > %s" %(self.path, species_augustus, self.fasta, self.outfile)
        subprocess.check_call(cmd, shell=True)

    @staticmethod
    def revcomp(sequence) :
        complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
        return "".join(complement.get(letter, "?") for letter in sequence[::-1])

    def get_features(self) :
        fname = self.outfile
        fasta = self.fasta
        genes = []
        fdata = SeqIO.to_dict(SeqIO.parse(fasta, "fasta"))

        with open(fname) as f :
            for line in f :
                if line.startswith("# start gene") :
                    gene = Augustus._read_gene_lines(f)
                    genes.append(gene)


        for gene in genes :
            for transcript in gene.childs :
                fpositions = transcript.childs_positions()
                fpositions = sorted(fpositions, key = lambda x : x[0])
                chromosome = transcript.chromosome

                sequence = "".join(str(fdata[chromosome][segment.start:segment.stop].seq)
                for segment in SegList(fpositions))
                if transcript.strand == "-" : sequence = Augustus.revcomp(sequence)

                transcript.sequences["CDS"] = sequence

        return genes

    @staticmethod
    def _read_gene_lines(handler) :
        feature = None
        protein = ""

        for line in handler :
            line = line.strip()

            if line.startswith("# end gene") :
                feature.sequences["protein"] = protein
                return feature

            if line.startswith("#") :
                # protein sequence
                line = line.split()[-1]
                line = line[1:] if line.startswith("[") else line
                line = line[:-1] if line.endswith("]") else line
                protein += line

            else :
                line = line.split("\t")
                chromosome, soft, kind, start, end, _, strand, _, info = line
                start = int(start) - 1
                name = info if "transcript_id" not in info else None
                f = Feature(name, kind, chromosome, strand, start, end)

                if name is None :
                    tid = info.split('"')[1]
                    transcript = feature.find_feature(tid)
                    transcript.add_child(f)

                else :
                    if kind == "gene" :
                        feature = f
                    else :
                        feature.add_child(f)
