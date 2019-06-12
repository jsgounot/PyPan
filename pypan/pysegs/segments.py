class Segment() :
    
    def __init__(self, start, stop) :
        self.start, self.stop = sorted((int(start), int(stop)))

    def __str__(self) :
        return "Segment(%i - %i)" %(self.start, self.stop)

    def __repr__(self) :
        return str(self)

    def __contains__(self, value) :
        return self.overlapp(other) if isinstance(other, Segment) else self.isin(other)

    def __hash__(self) :
        return hash((self.start, self.stop))

    def __eq__(self, other) :
        return hash(self) == hash(other) if isinstance(other, Segment) else False

    def __len__(self) :
        return self.stop - self.start

    def isin(self, coord) :
        return self.start <= coord <= self.stop 

    def segin(self, segment) :
        return segment.isin(self.start) or segment.isin(self.stop)

    def overlapp(self, segment) :
        return self.segin(segment) or segment.segin(self)

    def overlapp_count(self, other, prc=False) :
        if isinstance(other, Segment) :
            size = len(self) - sum(len(subseg) for subseg in self - other)

        elif isinstance(other, SegList) :
            size = other.overlapp_count(self, prc=False)

        else :
            raise ValueError("Can only calculate overlapp against Segment or SegList objects")

        size = size if size >= 0 else 0
        if prc : size = size * 100 / len(self)
        return size

    def __add__(self, other) :
        
        if isinstance(other, SegList) :
            return other + self

        elif not isinstance(other, Segment) :
            raise ValueError("Can only add another segment")

        if self.overlapp(other) :
            return SegList((Segment(min((self.start, other.start)), max(self.stop, other.stop)), ),
                skip_check=True)

        else :
            return SegList((self, other), skip_check=True)

    def __sub__(self, other) :
        
        if isinstance(other, SegList) :
            return SegList([self]) - other

        elif not isinstance(other, Segment) :
            raise ValueError("Can only substract another segment")

        if self.overlapp(other) :
            segments = []

            if self.start < other.start :
                segments.append(Segment(self.start, other.start))

            if self.stop > other.stop :
                start = other.stop
                segments.append(Segment(start, self.stop))

            return SegList(segments)

        else :
            return SegList((self, other))

class SegList() :
    def __init__(self, segments=[], skip_check=False) :
        self._segments = []

        if not segments :
            pass

        elif isinstance(segments, SegList) :
            self._segments = segments._segments.deepcopy()

        elif not skip_check :
            segments = [segment if isinstance(segment, Segment) else Segment(* segment)
            for segment in segments]
            self.extend(segments)

        else :
            self._segments = list(segments)

    def __bool__(self) :
        return bool(self._segments)

    def __hash__(self) :
        return hash(tuple((seg.start, seg.stop) for seg in self))

    def __eq__(self, other) :
        return hash(self) == hash(other) if isinstance(other, SegList) else False

    def __iter__(self) :
        yield from self._segments

    def __getitem__(self, * args, ** kwargs) :
        return self._segments.__getitem__(* args, ** kwargs)

    def __str__(self) :
        return "SegList(%s)" %(self._segments)

    def __repr__(self) :
        return str(self)

    def __len__(self) :
        return len(self._segments)

    def __add__(self, other) :

        if isinstance(other, Segment) :
            return self._add_segments(SegList((other, )))

        elif isinstance(other, SegList) :
            return self._add_segments(other)

        else :
            raise ValueError("Only segment or SegList can be added")

    def __sub__(self, other) :

        if isinstance(other, Segment) :
            return self._sub_segments(SegList((other, )))

        elif isinstance(other, SegList) :
            return self._sub_segments(other)

        else :
            raise ValueError("Only segment or SegList can be used")

    def __contains__(self, value) :
        return any(value in segment for segment in self)

    def append(self, segment) :
        if not isinstance(segment, Segment) :
            raise ValueError("Append only segment object")

        self.extend(SegList((segment, )))

    def min_coor(self, default=None) :
        if not self._segments : return default
        return min(segment.start for segment in self._segments)

    def max_coor(self, default=None) :
        if not self._segments : return default
        return max(segment.stop for segment in self._segments)        

    def cum_size(self) :
        return sum(len(segment) for segment in self)

    def _create_from_list(self, seglist) :
        seglist = sorted(seglist, key = lambda x : x.start)
        lenlist = len(seglist)
        idx = 0

        if len(seglist) == 1 : return seglist

        while idx + 1 != lenlist :
            res = seglist[idx] + seglist[idx + 1]
            
            if len(res) == 1 :
                seglist[idx:idx+2] = list(res)
                lenlist -= 1

            else :
                idx += 1

        return seglist

    def _add_segments(self, seglist) :
        copy_sl = SegList(self._segments)
        copy_sl.extend(seglist)
        return copy_sl

    def _sub_segments(self, seglist) :
        copy_sl = SegList(self._segments)
        copy_sl.remove(seglist)
        return copy_sl

    def pairwise_finder(self, seglist) :
        # find either overlapping segments
        # or the idx position of the segmentation list
        # if not found

        if not isinstance(seglist, SegList) :
            raise ValueError()

        add_idx = 0

        for other_segment in seglist :
            start_seg = other_segment.start
            found = False
            segments = []

            for idx, self_segment in enumerate(self[add_idx:], start=add_idx) :

                if found : break

                elif self_segment.overlapp(other_segment) :
                    add_idx = idx
                    found = True

                    for self_segment in self[idx:] :
                        if not self_segment.overlapp(other_segment) :
                            break

                        segments.append(self_segment)

                    break

                elif start_seg < self_segment.start :
                    found = True
                    add_idx = idx
                    break

            if not found : idx += 1                    
            yield idx, other_segment, segments

    def overlapp_count(self, other, prc=False) :

        if isinstance(other, Segment) or isinstance(other, SegList) :
            size = self.cum_size() - sum(len(subseg) for subseg in self - other)

        else :
            raise ValueError("Can only calculate overlapp against Segment or SegList objects")

        if prc : size = size * 100 / self.cum_size()
        return size

    def extend(self, seglist) :
       
        if not isinstance(seglist, SegList) :
            seglist = SegList(self._create_from_list(seglist), skip_check=True)

        if not self :
            self._segments = list(seglist)
            return

        for idx, other_segment, self_segments in self.pairwise_finder(seglist) :
            new_segment = sum(self_segments, other_segment)
            new_segment = [new_segment] if isinstance(new_segment, Segment) else list(new_segment)           
            self._segments[idx:idx+len(self_segments)] = new_segment
    
    def remove(self, seglist) :

        if not isinstance(seglist, SegList) :
            seglist = SegList(self._create_from_list(seglist), skip_check=True)

        if not self or not seglist:
            return

        for idx, other_segment, self_segments in self.pairwise_finder(seglist) :

            new_segments = []
            for segment in self_segments :
                new_seg = segment - other_segment
                if new_seg : new_segments.extend(list(new_seg))

            self._segments[idx:idx+len(self_segments)] = list(new_segments)
