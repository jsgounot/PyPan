# -*- coding: utf-8 -*-
# @Author: jsgounot
# @Date:   2019-01-22 13:42:09
# @Last modified by:   jsgounot
# @Last Modified time: 2019-01-23 12:59:56

import os, subprocess
from tempfile import NamedTemporaryFile
import pandas as pd
from pypan.utils import get_soft, get_field

class SNAP() :

    def __init__(self, fasta, species) :
        self.fasta = fasta
        self.hmm = get_field("PROJECT", species, "snap")


    @property
    def path(self):
        return get_soft("snap")

    def run(self, outfile=None) :
        if outfile is None : outfile = self.fasta

        tmpfile = NamedTemporaryFile(delete=False).name
        nucl = outfile + ".snap.nucl.fasta"
        prot = outfile + ".snap.prot.fasta"

        cmd = "%s %s %s -aa %s -tx %s > %s" %(self.path, self.hmm, self.fasta, prot, nucl, tmpfile)
        subprocess.check_call(cmd, shell=True)

        outfile = outfile + ".snap.info.tsv"
        self.parse_res(tmpfile, outfile, delete=True)

    def parse_res(self, resfile, outfile, delete=True) :
        header = ["Label", "Start", "End", "Strand", "Score", "5' overhang", "3' overhang", "Frame", "Group"]
        data = []
        with open(resfile) as f :
            for line in f :
                if line.startswith("scoring....") or line.startswith(">") :
                    continue
                else :
                    line = line.strip().split("\t")
                    line = {column : line[idx] for idx, column in enumerate(header)}
                    data.append(line)

        df = pd.DataFrame(data)
        df.to_csv(outfile, sep="\t")
        if delete : os.remove(resfile)