# Description: Input/output routines.

import logging
logger = logging.getLogger(__name__)
logger.addHandler(logging.NullHandler())

import os
import gzip
from collections import namedtuple
from itertools import groupby
from seq import Sequence, QSequence
from clone import Clone, CloneSet
from allele import Allele, AlleleContainer
from barcode import Adapter, Barcode
from math import log

def filesize(f):
    """Return size of file.

    :f: file-like object
    """
    pos = f.tell()
    f.seek(0, os.SEEK_END)
    size = f.tell()
    f.seek(pos, os.SEEK_SET)
    return size

def identity(x):
    return x

def solexa2phred(qual):
    """Convert vector of Sanger quality scores to Phred quality scores.
    
    The converted Phred score is rounded to the nearest integer and returned
    as an int."""
    return [ int(round(10 * log(10 ** (q_solexa / 10.0) + 1, 10))) \
            for q_solexa in qual]

def decode_ascii(qual_str, offset, min_q, max_q):
    """Convert string of ascii characters to corresponding quality scores."""
    qual = [None] * len(qual_str)
    for i,qch in enumerate(qual_str):
        q = ord(qch) - offset
        if q < min_q or q > max_q:
            raise Exception("Quality score (%s) not in [%s, %s]"%(
                    q, min_q, max_q))
        qual[i] = q
    return qual

MIN_PHRED_QUAL = 0
MAX_PHRED_QUAL = 41

FastqRecord = namedtuple("FastqRecord", [
    "id", # sequence identifier
    "desc", # description
    "seq", # sequence
    "qual_str", # quality (ascii encoded)
    ])
QualityEncoding = namedtuple("QualityEncoding", [
    "name",
    "qual2phred",
    "offset",
    "min_q",
    "max_q"])

SangerEncoding = QualityEncoding("Sanger", identity, 33, 0, 40)
SolexaEncoding = QualityEncoding("Solexa", solexa2phred, 64, -5, 40)
Illumina1_3_Encoding = QualityEncoding("Illumina1.3", identity, 64, 0, 40)
Illumina1_8_Encoding = QualityEncoding("Illumina1.8", identity, 33, 0, 41)

Encodings = [SangerEncoding, SolexaEncoding, Illumina1_3_Encoding,
        Illumina1_8_Encoding]

ASCII2PHRED = {
    encoding : lambda qual_str:
        qual2phred(decode_ascii(qual_str, offset, min_q, max_q))
    for encoding, qual2phred, offset, min_q, max_q in Encodings}

def zopen(fn, *args, **kwargs):
    """Open normal or zipped file."""
    if fn[-3:] == ".gz":
        return gzip.open(fn, *args, **kwargs)
    else:
        return open(fn, *args, **kwargs)

class FileFormat(object):
    """Interface for classes that implement reading and/or writing a particular
    file format.
    """
    def __init__(self):
        raise NotImplementedError("Cannot instantiate FileFormat.")

    @staticmethod
    def records_in(handle):
        """Generator function yielding records from handle. Type of records is
        implementation dependent.

        :handle: a file-like object, supporting the "with handle:" construct
        """
        raise NotImplementedError()

    @staticmethod
    def records_out(handle, records):
        """Writes records to handle.

        :handle: a file-like object, supporting the "with handle:" construct
        """
        raise NotImplementedError()

class FastaFormat(FileFormat):
    @staticmethod
    def records_in(handle):
        """Generator function to iterate over fasta records as Sequence
        objects.
        """
        with handle:
            header = None
            for (is_header, lines) in groupby(handle, \
                    lambda line : line[0] == '>'):
                if is_header:
                    header = lines.next().rstrip()[1:]
                else:
                    if header:
                        seq = "".join([line.rstrip() \
                                for line in lines]).upper()
                        yield Sequence(name = header, seq = seq)
                    else: # skip comment section
                        continue
    @staticmethod
    def records_out(handle, records):
        with handle:
            for record in records:
                handle.write(">%s\n%s\n"%(record.name, record.seq))

class FastqFormat(FileFormat):
    class FastqOutputStream(object):
        def __init__(self, handle):
            self._handle = handle

        def __add__(self, record):
            self._handle.write("\n".join(["@" + record.name, record.seq, "+",
                "".join([chr(q + 33) for q in record.qual])]) + "\n")
            return self

    @staticmethod
    def records_in(handle, encoding = "Illumina1.8"):
        """Generator function to iterate over fastq records as QSequence
        objects.

        :encoding: type of base quality scores. If None, FastqRecord objects
        will be returned.
        """
        with handle:
            i = 0 # define i (for in case handle is empty)
            for i, line in enumerate(handle):
                if not line:
                    break

                if i % 4 == 0:
                    if not line[0] == '@':
                        raise Exception("@ expected, got %s."%line[0])
                    #fields = line.rstrip()[1:].split(' ')
                    #record_id = fields[0]
                    #desc = ' '.join(fields[1:])
                    record_id = line.rstrip()[1:]
                elif i % 4 == 1:
                    seq = line.rstrip()
                elif i % 4 == 2:
                    if line[0] != '+':
                        raise Exception("+ expected, got %s."%line[0])
                elif i % 4 == 3:
                    qual_str = line.rstrip()
                    if encoding is None:
                        fields = record_id.split(" ")
                        yield FastqRecord(id = fields[0],
                                desc = "".join(fields[1:]),
                                seq = seq.upper(),
                                qual_str = qual_str)
                    else:
                        yield QSequence(name = record_id, seq = seq.upper(),
                            qual = ASCII2PHRED[encoding](qual_str),
                            force_N_to_zero = True)
                i += 1

            if i % 4 != 0:
                raise Exception("Unexpected end of file.")

    @staticmethod
    def records_out(handle, records):
        """Output QSequence-like objects as FastaQ records, using Illumina1.8
        encoding.

        :records: sequence-like object containing QSequence-like objects. If
        None, a FastqOutputStream object.
        """
        if records is None:
            return FastqFormat.FastqOutputStream(handle)
        with handle:
            for record in records:
                handle.write("\n".join(["@" + record.name, record.seq, "+",
                    "".join([chr(q + 33) for q in record.qual])]) + "\n")

SAMRecord = namedtuple("SAMRecord", [
    "QNAME",    # Query template NAME
    "FLAG",     # bitwise FLAG
    "RNAME",    # Reference sequence NAME
    "POS",      # 1-based leftmost mapping POSition of first (mis)matching base
    "MAPQ",     # MAPping Quality
    "CIGAR",    # CIGAR string
    "RNEXT",    # Ref. name of the mate/next read
    "PNEXT",    # Position of the mate/next read 
    "TLEN",     # observed Template LENgth
    "SEQ",      # segment SEQuence
    "QUAL",     # ASCII of Phred-scaled base QUALity+33
])

def string2SAMRecord(s):
    fields = s.rstrip().split('\t')
    if len(fields) < 11:
        return None

    return SAMRecord(
            QNAME   = fields[0],
            FLAG    = int(fields[1]),
            RNAME   = fields[2],
            POS     = int(fields[3]),
            MAPQ    = int(fields[4]),
            CIGAR   = fields[5],
            RNEXT   = fields[6],
            PNEXT   = int(fields[7]),
            TLEN    = int(fields[8]),
            SEQ     = fields[9],
            QUAL    = fields[10])

class SAMFormat(FileFormat):
    @staticmethod
    def records_in(handle):
        with handle:
            for line in handle:
                rec = string2SAMRecord(line)
                if not rec is None:
                    yield rec

class LegacyCloneSetFormat(FileFormat):
    """Reads and writes cloneset files of RTCR version 0.2.4 or older."""
    @staticmethod
    def records_in(handle):
        def get_fieldnr(names, fields):
            if isinstance(names, basestring):
                names = [names]
            for name in names:
                if name in fields:
                    return fields.index(name)

        cloneset = CloneSet()
        with handle:
            header = handle.readline().rstrip().split("\t")
            nt_jnc_col = get_fieldnr(["nt_jnc",
                "Junction nucleotide sequence"], header)
            orig_depth_col = get_fieldnr("orig_depth", header)
            depth_col = get_fieldnr("depth", header)
            qual_col = get_fieldnr("qvec", header)

            for line in handle:
                fields = line.rstrip().split("\t")
                nt_jnc = fields[nt_jnc_col]
                orig_depth = 0 if orig_depth_col is None else \
                        int(fields[orig_depth_col])
                depth = int(fields[depth_col])
                qual = map(int, fields[qual_col].split("|"))
                clone = Clone((nt_jnc, qual, depth), refpos = 0,
                        orig_count = orig_depth)
                cloneset.add(clone)
        return cloneset

    @staticmethod
    def records_out(handle, cloneset):
        """Warning: saving new CloneSet objects may result in information loss
        because the old format does not contain position frequency information.
        """
        with handle:
            handle.write("vid\tjid\tnt_jnc\torig_depth\tdepth\tve\tjs\tqvec\n")
            for clone in cloneset:
                if hasattr(clone, "v"):
                    vid = clone.v.allele.name
                    jid = clone.j.allele.name
                    ve = clone.v.end
                    js = clone.j.start
                else:
                    vid = None
                    jid = None
                    ve = 0
                    js = 0
                handle.write("\t".join(map(str,
                        [vid,
                         jid,
                         clone.seq,
                         clone.orig_count,
                         clone.count,
                         ve,
                         js,
                         "|".join(map(str, clone.qual))])) + "\n")

class LegacyAlleleContainerFormat(FileFormat):
    """Reads immune reference files of RTCR version 0.2.4 or older."""
    @staticmethod
    def records_in(handle):
        fa = FastaFormat.records_in(handle)
        species = "HomoSapiens"
        functionality = ""
        acc_nr = ""
        for rec in fa:
            fields = rec.name.split('|')
            name = fields[0]
            region = name[3]
            if region == 'V' or region == 'J':
                refpos = int(fields[1].split("refpos:")[1])
            else:
                refpos = -1
            yield Allele(species, name, functionality, refpos, acc_nr, rec.seq)

class IMGTAlleleContainerFormat(FileFormat):
    """Reads IMGT immune reference fasta."""
    @staticmethod
    def records_in(handle):
        fa = FastaFormat(handle)
        for rec in fa:
            fields = rec.name.split("|")
            if len(fields) <= 5:
                raise ValueError("In an IMGT reference fasta, \
record titles contain at least 5 fields delimited by \"|\".")
            acc_nr, name, species, functionality, region = fields[:5]
                
            # Turn sequence to upper case
            seq = rec.seq.upper()

            # Concatenate and capitalize species name
            # e.g. "Homo sapiens" becomes "HomoSapiens"
            species = species.split(" ")
            species[1] = species[1].capitalize() # e.g. Sapiens
            species = ''.join(species)

            # Create short name for region, e.g. "V-REGION"
            region = region.replace("-REGION","") 

            region_from_name = name[3]
            if region != region_from_name and region_from_name != 'C':
                raise ValueError("Region from allele name (%s) does \
not match region field (%s)."%(region, region_from_name))

            refpos = -1
                
            # Searching for reference position in V and J regions.
            if region == 'V':
                # For aligned alleles the conserved cysteine
                # is at position 104 (aa index).
                # So the CDR3 start (right after the cys),
                # is at 104 * 3 = 312
                offset = seq[:312].count('.') # Count gaps
                refpos = 312 - offset
            elif region == 'J':
                # Search for conserved Phe/Trp-Gly-X-Gly motif
                motif="((TT[TC])|(TGG))GG[TCAG]{4}GG[TCAG]{1}"
                m = re.search(motif,seq)
                if m == None: # motif not found
                    liberal_motif = \
                        "[TCAG]{3}GG[TCAG]{4}[TCAG]{1}"
                    m = re.search(liberal_motif,seq)
                    if m == None:
                        raise Exception("Unable to find X-Gly-X-Gly motif in \
J-REGION allele in %s"%name)
                    else:
                        logger.warning("Had to use liberal X-Gly-X-Gly motif \
to find J-REGION reference position in %s"%name)
                if seq.count('.'):
                    raise Exception("Gaps \".\" found in J-REGION allele %s"%\
                            name)
                refpos = m.start()

            seq = seq.replace('.','')
            yield Allele(species, name, functionality, refpos, acc_nr, seq)

class AlleleContainerFormat(FileFormat):
    @staticmethod
    def records_in(handle):
        """Yields immune receptor reference sequences from RTCR tsv file."""
        with handle:
            header = handle.readline()
            if header != "Species\tAllele name\tFunctionality\t\
Reference position\tAccession number\tSequence\n":
                raise ValueError("Wrong header for a RTCR reference file.")
            for line in handle:
                fields = line.rstrip().split("\t")
                if len(fields) != 6:
                    raise ValueError("Wrong number of fields")
                species, name, functionality, refpos, acc_nr, seq = fields
                yield Allele(species, name, functionality, int(refpos),
                    acc_nr, seq)

    @staticmethod
    def records_out(handle, records):
        """Saves AlleleContainer object to file.
        
        :records: AlleleContainer instance
        """
        with handle:
            handle.write("Species\tAllele name\tFunctionality\t\
Reference position\tAccession number\tSequence\n")
            for a in sorted(records.alleles.itervalues(), key = lambda a:
                    (a.species, a.name, a.functionality, a.refpos, a.acc_nr,
                        a.seq)):
                handle.write("%s\t%s\t%s\t%s\t%s\t%s\n"%(
                    a.species, a.name, a.functionality, a.refpos, a.acc_nr,
                    a.seq))

def read_reference(fn):
    """Reads immune reference sequences and creates an AlleleContainer object.
    """
    ref = AlleleContainer()
    get_allele = None

    with zopen(fn, 'r') as f:
        line = f.readline().rstrip()
        if line[0] == '>':
            if "refpos" in line:
                get_allele = LegacyAlleleContainerFormat.records_in 
            else:
                get_allele = IMGTAlleleContainerFormat.records_in
        elif "Species" in line:
            get_allele = AlleleContainerFormat.records_in
    if get_allele == None:
        raise ValueError("Unknown filetype")
    
    with zopen(fn,'r') as handle:
        for allele in get_allele(handle):
            ref.add(allele)
    return ref

def junctions2cloneset(handle, phred_encoding):
    """Converts a jnc.tsv file to a CloneSet object.
    
    Note, VJ annotation is lost.
    """
    cloneset = CloneSet()
    simple_clones = {} 
    with handle:
        header  = handle.readline().rstrip()
        assert header == "read_name\tvid\tjid\tnt_jnc\tq_jnc\tjnc_ve\tjnc_js\t\
nt_jnc_gr\tq_jnc_gr\tjnc_ve_gr\tjnc_js_gr\tv_score\tj_score"
        for line in handle:
            read_name, vid, jid, nt_jnc, q_jnc, jnc_ve, jnc_js, nt_jnc_gr,\
                    q_jnc_gr, jnc_ve_gr, jnc_js_gr,\
                    v_score,j_score = line.rstrip().split('\t')

            simple_clone = simple_clones.get(nt_jnc, [0, [0]*len(nt_jnc)])
            simple_clone[0] += 1 # increase clone's count
            q_vec = simple_clone[1]
            for i in xrange(len(nt_jnc)):
                q_vec[i] = max(q_vec[i], ord(q_jnc[i]) - phred_encoding)
            simple_clones[nt_jnc] = simple_clone
    for nt_jnc, (count, q_vec) in simple_clones.iteritems():
        cloneset.add(Clone((nt_jnc, q_vec, count)))
    return cloneset

class BarcodeFormat(FileFormat):
    @staticmethod
    def records_in(handle):
        """Return the fields in a barcode file as a 3-tuple containing
        sample ID (string), master adapter (Adapter), and slave adapter
        (Adapter).
        
        For example:
        ("S2-1-alpha", "AGACAcagtggtatcaacgcagagtNNNNtNNNNtNNNNtctt",
        "ACTGAgggtcagggttctggatat")
        
        If no slave adapter (3rd element in the tuple) is provided, it will be
        None.
        """
        with handle:
            while True:
                line = handle.readline().rstrip()
                if not line:
                    break
                fields = line.split("\t")
                sample_id = fields[0]
                master_adapter = Adapter(fields[1])
                
                if len(fields) == 2:
                    yield Barcode(sample_id, master_adapter, None)
                elif len(fields) == 3:
                    yield Barcode(sample_id, master_adapter,
                            Adapter(fields[2]))
                else:
                    raise Exception("Expected 2 or 3 fields, got %s"%len(fields))
