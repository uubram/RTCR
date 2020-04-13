from collections import namedtuple
import string

Barcode = namedtuple("Barcode",["sample_id", "master", "slave"])

DNA_COMPLEMENT = string.maketrans("ATCG","TAGC")

IUPAC = {
    "A" : ["A"],
    "C" : ["C"],
    "G" : ["G"],
    "T" : ["T"],
    "U" : ["U"],
    "R" : ["A", "G"],
    "Y" : ["C", "T"],
    "S" : ["G", "C"],
    "W" : ["A", "T"],
    "K" : ["G", "T"],
    "M" : ["A", "C"],
    "B" : ["C", "G", "T"],
    "D" : ["A", "G", "T"],
    "H" : ["A", "C", "T"],
    "V" : ["A", "C", "G"],
    "N" : ["A", "C", "G", "T"]
}

# TODO: use seq module for taking reverse complement
def revcomp(seq):
    return seq.translate(DNA_COMPLEMENT)[::-1]

class Adapter(object):
    def __init__(self, adapter_str):
        seed_ranges = get_uppercase_ranges(adapter_str, skipchars = "N")
        if len(seed_ranges) == 0:
            raise Exception("No seed sequence provided for adapter.")
        if len(seed_ranges) > 1:
            raise Exception("Adapter contains multiple seeds.")

        self._seed_start = seed_ranges[0][0]
        self._seed_end = seed_ranges[0][1]
        self._seed = adapter_str[self._seed_start : self._seed_end]
        assert(self._seed.isupper())
        if any([not ch in "ACTG" for ch in self._seed]):
            raise Exception("Adapter seed contains ambiguous or non-DNA bases.") 

        self._UMI_ranges = get_ranges(adapter_str,'N')

        if len(self._UMI_ranges) > 0:
            self._UMI_length = sum([r[1] - r[0] for r in self._UMI_ranges])
            self._UMI_start = self._UMI_ranges[0][0]
            self._UMI_end = self._UMI_ranges[-1][1]
        else:
            self._UMI_length = 0
            self._UMI_start = None
            self._UMI_end = None

        self._adapter_str = adapter_str
        self._seq = adapter_str.upper()
        self._rc_seq = revcomp(self._seq)
        self._rc_seed = revcomp(self._seed)
        self._rc_seed_start = len(self._adapter_str) - self._seed_end
        self._rc_seed_end = self._rc_seed_start + len(self._seed)

        if any([not ch in IUPAC for ch in self._seq]):
            raise Exception("Adapter contains non-IUPAC codes.")

    def has_UMI(self):
        return len(self._UMI_ranges) > 0

    @property
    def seed(self):
        return self._seed

    @property
    def seq(self):
        return self._seq

    @property
    def UMI_ranges(self):
        return self._UMI_ranges

    @property
    def UMI_length(self):
        return self._UMI_length

    @property
    def UMI_start(self):
        return self._UMI_start

    @property
    def UMI_end(self):
        return self._UMI_end

    def get_UMI(self, seq, qual, offset):
        return ("".join([seq[offset + r[0] : offset + r[1]] \
                for r in self._UMI_ranges]),\
                "".join([qual[offset + r[0] : offset + r[1]] \
                for r in self._UMI_ranges]))

    # Locates adapter in sequence
    def locate_in(self, seq, max_mm, search_rc):
        matches = sorted(fuzzy_findall(query = self._seq, target = seq,
                seed_start = self._seed_start, seed_end = self._seed_end,
                max_mm = max_mm))
        matches_rc = None
        if search_rc:
            matches_rc = sorted(fuzzy_findall(query = self._rc_seq,
                target = seq, seed_start = self._rc_seed_start,
                seed_end = self._rc_seed_end, max_mm = max_mm))
        return matches, matches_rc

    def __str__(self):
        return self._adapter_str

# Finds uppercase ranges in a string. 
def get_uppercase_ranges(s, skipchars = ""):
    ranges = []
    start = 0
    end = 0
    while True:
        if end == len(s) or s[end] in skipchars or not s[end].isupper():
            if end > start:
                ranges.append((start,end))
            if end == len(s):
                return ranges
            start = end + 1
        end += 1

# Finds ranges of a given character in a string.
def get_ranges(s, ch):
    assert(len(ch) == 1)
    ranges = []
    start = 0
    end = 0
    while True:
        if end == len(s) or s[end] != ch:
            if end > start:
                ranges.append((start,end))
            if end == len(s):
                return ranges
            start = end + 1
        end += 1

# Finds longest match between two sequences with at most max_mm mismatches.
# Note: this method does not perform an alignment.
# query is allowed to contain IUPAC codes. Whereas target should only contain
# ACTGU
def fuzzy_match_len(query, target, max_mm, direction, qstart, tstart):
    assert(direction in (1,-1))
    match_len = 0
    mm = 0 # number of mismatches
    qpos = qstart
    tpos = tstart
    while True:
        if target[tpos] in IUPAC[query[qpos]]:
            match_len = (qpos - qstart) + 1
        else:
            if mm == max_mm:
                break
            mm += 1
        qpos += direction
        tpos += direction

        if qpos < 0 or tpos < 0 or qpos >= len(query) or tpos >= len(target):
            break
    return match_len, mm

# Finds all matches of query with target, allowing for overlapped matches.
def overlapped_findall(query, target):
    start = 0
    while True:
        start = target.find(query, start) + 1
        if start > 0:
            yield start - 1
        else:
            break

# Finds all matches between query's seed sequence and the target sequence,
# returning tuples with the number of mismatches between query and target
# for each seed match found, and the position at which the query matched the
# target.
def fuzzy_findall(query, target, seed_start, seed_end, max_mm):
    # Search the target for the seed sequence
    seed = query[seed_start: seed_end]
    for i in overlapped_findall(seed, target):
        match_start = i - seed_start
        match_len, mm = fuzzy_match_len(query, target, max_mm, 1, 0,
                match_start)
        if match_len == len(query):
            yield (mm, match_start)
