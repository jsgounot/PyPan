# -*- coding: utf-8 -*-
# @Author: jsgounot
# @Date:   2019-01-22 15:14:11
# @Last modified by:   jsgounot
# @Last Modified time: 2019-01-22 15:14:40

class Feature() :

    def __init__(self, name, kind, chromosome, strand, start, end) :
        self.name = name
        self.chromosome = chromosome
        self.kind = kind
        self.strand = strand
        self.position = (int(start), int(end))
        self.sequences = {}
        self.childs = []

    def add_child(self, Feature) :
        self.childs.append(Feature)

    def find_feature(self, name) :
        if name != self.name :
            for feature in self.childs :
                return feature.find_feature(name)
        else :
            return self

    def childs_positions(self) :
        used = {"CDS", "start_codon", "stop_codon"}
        return [child.position for child in self.childs if child.kind in used]

    @staticmethod
    def features_kind(features) :
        kinds = set()
        for feature in features :
            kinds.add(feature.kind)
            kinds = kinds | Feature.features_kind(feature.childs)
        return kinds

    def length(self) :
        return self.position[1] - self.position[0]

    def posi_to_interval(self) :
        return interval[self.position[0], self.position[1]]

    def posi_overlapping(self, feature, asprc=False) :
        if self.chromosome != feature.chromosome :
            return 0

        else :
            common = self.posi_to_interval() & feature.posi_to_interval()
            length = sum(c.sup - c.inf for c in common)

            if asprc : length = (length / self.length()) * 100
            return length

    def share_edge(self, feature) :
        if self.chromosome != feature.chromosome :
            return False, False

        else :
            sstart, send = self.position
            fstart, fend = feature.position
            return sstart == fstart, send == fend