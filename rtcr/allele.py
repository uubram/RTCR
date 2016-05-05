import collections
import re

class Allele(object):
    """Holds sequence and other information about a segment allele.
    
    Attributes:
        :name: full name of the allele, e.g. TRBV1*01
        :species: species the allele is part of, e.g. HomoSapiens
        :gene: e.g. TRB (always 3 characters)
        :region: e.g. V (always a single character)
        :segment: e.g. 10-2 (string)
        :refpos: reference position, 0-based. For V region it corresponds to
            first base after the conserved Cys, i.e. start of the CDR3.
            For J region it corresponds to the first base of the
            F/W-G-X-G motif, i.e. index of the first base after CDR3.
    """
    def __init__(self, species, name, functionality, refpos, acc_nr, seq):
        gene, region, segment, allele = Allele.parse_allele_name(name)

        if not isinstance(name, basestring):
            raise ValueError("Allele name should be a string")

        if not isinstance(species, basestring):
            raise ValueError("Allele species should be a string")

        if not gene[:2] in ("TR","IG"):
            raise ValueError("Allele gene should start with TR \
or IG.")
        
        if not region in tuple("VDJC"):
            raise ValueError("Allele region should be one of V,D,\
J, or C.")

        if segment == "" and region != 'C':
            raise ValueError("Allele segment id missing.")

        if re.search("[^ACTGN]", seq) != None:
            raise ValueError("Allele sequence should consist only \
of ACTGN characters.")
        
        if not isinstance(refpos, int):
            raise ValueError("Allele reference position should \
be an integer.")

        if refpos < -1:
            raise ValueError("Illegal reference position. \
Allele: %s; Sequence length: %s; Reference position: %s"%(
    name, len(seq), refpos))

        if not isinstance(acc_nr,str):
            raise ValueError("Accession number should be a string.")
        
        if not functionality in ("F","P","[F]","(F)","[P]","(P)",\
                "ORF","[ORF]",""):
            raise ValueError("Unknown type of functionality \
\"%s\" for allele %s."%(functionality,name))
        self.name = name
        self.species = species
        self.gene = gene
        self.region = region
        self.segment = segment
        self.allele = allele
        self.seq = seq
        self.refpos = refpos # 0-based reference position
        self.acc_nr = acc_nr # accession number
        self.functionality = functionality

    def __eq__(self, other):
        if isinstance(other, Allele):
            return self.name == other.name and self.species == other.species \
                    and self.gene == other.gene \
                    and self.region == other.region \
                    and self.segment == other.segment \
                    and self.allele == other.allele \
                    and self.seq == other.seq \
                    and self.refpos == other.refpos \
                    and self.acc_nr == other.acc_nr \
                    and self.functionality == other.functionality
        else:
            raise NotImplementedError()

    def __str__(self):
        return "%s:%s"%(self.species,self.name)

    @staticmethod
    def parse_allele_name(name):
        gene = name[:3]
        region = name[3]
        segment = name.split("*")[0][4:]
        allele = name.split("*")[1]
        return (gene, region, segment, allele)

class AlleleContainer(object):
    """Container for Allele objects.
    
    To identify an allele, there is the following hierarchy:
    Species -> Gene -> Region -> Segment -> Allele
    e.g. HomoSapiens -> TRB -> V -> 1 -> 1 (this  would correspond to TRBV1*01)
    """

    def __init__(self, alleles = None):
        """alleles: dict, or list with Allele objects"""

        self._species = set() # species in this container 
        self._genes = set()
        self._regions = set()

        self._alleles = {}
        self._alleles_by_name = {}
        if not alleles is None:
            if isinstance(alleles, dict):
                alleles = alleles.itervalues()
            for allele in alleles:
                self.add(allele)
    
    @property
    def species(self):
        return self._species

    @property
    def genes(self):
        return self._genes

    @property
    def regions(self):
        return self._regions

    @property
    def alleles(self):
        return self._alleles

    def get_alleles(self, name = None,
            species = None,
            gene = None,
            region = None,
            segment = None,
            allele = None,
            functionality = None,
            acc_nr = None):
        """Returns a list of all alleles corresponding to the keywords that
        are not None. Keywords are allowed to be lists.
        """
        def matches(x, y):
            if isinstance(y, basestring):
                return x == y
            elif isinstance(y, collections.Iterable):
                    return x in y
            elif y is None:
                return False
            else:
                raise ValueError("not None or iterable")

        matches = lambda x, y: (not y is None and x in y) or y is None
        res = []
        for a in self.alleles.itervalues():
            if matches(a.name, name) and \
                    matches(a.species, species) and \
                    matches(a.gene, gene) and \
                    matches(a.region, region) and \
                    matches(a.segment, segment) and \
                    matches(a.allele, allele) and \
                    matches(a.functionality, functionality) and \
                    matches(a.acc_nr, acc_nr):
                res.append(a)
        return res

    def get_slice(self,**kwargs):
        """Same as get_alleles method, except it returns an AlleleContainer
        object with the selected alleles.
        """
        return AlleleContainer(self.get_alleles(**kwargs))

    def __getitem__(self,key):
        """Returns Allele with name corresponding to the key. If multiple
        alleles match, an exception is raised.
        """
        return self._alleles_by_name[key]

    def _add_allele(self, allele):
        if allele.acc_nr in self._alleles:
            raise ValueError("Allele %s already exists."%allele)
        self._alleles[str(allele)] = allele
        self._species.add(allele.species)
        self._genes.add(allele.gene)
        self._regions.add(allele.region)

        # alleles_by_name dict is used for faster access to alleles
        # by their name. However, names that refer to multiple alleles
        # are removed. This functionality is primarily for when
        # only a single species is in the AlleleContainer object.
        if allele.name in self._alleles_by_name:
            del self._alleles_by_name[allele.name]
        else:
            self._alleles_by_name[allele.name] = allele

    def _add_container(self, container):
        for allele in self.container:
            self.add(allele)

    def add(self, other):
        if isinstance(other, Allele):
            self._add_allele(other)
        elif isinstance(other, AlleleContainer):
            self._add_container(other)
        else:
            raise NotImplementedError()

class RestrictedAlleleContainer(AlleleContainer):
    """Container for Allele-like objects, but only alleles that are from a
    particular species/gene/region/segment are allowed.
    """

    def __init__(self, species = None, gene = None, region = None,
            segment = None, **kwds):
        super(RestrictedAlleleContainer, self).__init__(**kwds)
        self._r_species = species
        self._r_gene = gene
        self._r_region = region
        self._r_segment = segment

    def _add_allele(self, allele):
        if (allele.species == self._r_species or self._r_species is None) and \
                (allele.gene == self._r_gene or self._r_gene is None) and \
                (allele.region == self._r_region or self._r_region is None) \
                and \
                (allele.segment == self._r_segment or \
                self._r_segment is None):
            super(RestrictedAlleleContainer, self)._add_allele(allele)
        else:
            raise ValueError("Allele does not meet restrictions")

    def _add_container(self, container):
        if self._r_species == container._r_species and \
                self._r_gene == container._r_gene and \
                self._r_region == container._r_region and \
                self.r_segment == container._r_segment:
            for allele in container:
                super(RestrictedAlleleContainer, self)._add_allele(allele)
        else:
            for allele in container:
                self._add_allele(allele)
