import sys
import re

class BasicRowParser:
    
    """
    Baseclass to access Contig, position and reference
    """
    
    def __init__(self, row, startpos=None):
        
        """
        The optional startpos argument is used to get rid of any extra fields before the actual vcf part of the rows
        """
        
        row = row.rstrip("\n")
        
        self.row = row
        self.fields = row.split("\t")
        
        if startpos is not None:
            
            self.fields = self.fields[startpos:]
        
        if row[0] != "#":
            
            self.valid_row = True
            
            self.contig_name = self.fields[0]
            self.pos = int(self.fields[1])
            self.ref = self.fields[3]
            
            #################################################
            # deal with C.elegans and expectatus:
	    
	    if "scaffold" in self.contig_name:
		self.con = int(self.fields[0].strip("scaffold"))
                self.conpos = "Contig" + str(self.con) + "\t" + str(self.pos)
	    
	    elif "CHROMOSOME" not in self.contig_name:
                self.con = int(self.fields[0].lstrip("Contig"))
                self.conpos = "Contig" + str(self.con) + "\t" + str(self.pos)
	    
            else:
                self.con = "NA"
                self.conpos = self.contig_name + "\t" + str(self.pos)
                
        else:
            self.valid_row = False

    
    ##########################################################
    # type questions
    
    def is_valid_row(self):
        """ test if the row contains all fields """
        return self.valid_row
    
    def is_ref_known(self):
        """ test if the ref is actually a base and not just "N" """

        if "N" in self.ref:
            return False
        else:
            return True
        
    def get_conpos(self):
        
        return self.conpos
    
    def get_con(self):
        return self.con
    
    def get_pos(self):
        return self.pos
    
    def get_alt(self):
        return self.alt
    
    def get_ref(self):
        return self.ref
    
    def cutoff_test(tested_trait, direction, cutoff):
        """
        Tests if the tested trait is [above/below] the cutoff
        """
        
        if direction == "max":
            if tested_trait <= cutoff:
                return True
            
        elif direction == "min":
            if tested_trait >= cutoff:
                return True
        else:
            return False


##############################################################################################################################################################################
##############################################################################################################################################################################


class ShortVCFrow(BasicRowParser):
    
    """
    Each instance corresponds to one site from a shortened VCF as in the 104 strains "corrected" files
    Contig0 317     .       T       T       31.5
 
    """
    
    def __init__(self, row):
        
        BasicRowParser.__init__(self, row)
        
        self.alt = self.fields[4]
        self.qual = float(self.fields[5])
        self.fq = None
        self.dp = None
        
    ##########################################################
    # type questions
    
    def is_ref(self):
        """ test if the alt equals the ref base """
        
        if self.alt == self.ref:
            return True
        else:
            return False
    
    def is_indel(self):
        
        if len(self.ref) == 1 and len(self.alt) == 1:
            return False
        else:
            return True
    
    
##############################################################################################################################################################################
##############################################################################################################################################################################


class RawBCFrow(BasicRowParser):
    
    """
    Each instance corresponds to one site from BCFtools view without SNP calling
    Contig11        123     .       T       X       0       .       DP=12;I16=8,3,1,0,359,12051,16,256,590,33950,60,3600,167,3273,21,441    PL:DP   0,33,208:12
 
    """
    
    def __init__(self, row):
        
        BasicRowParser.__init__(self, row)
        
        self.alt = self.fields[4]
        self.qual = float(self.fields[5])
        self.fq = None
        self.dp = None
        
    ##########################################################
    # type questions
    
    def is_ref(self):
        """ test if the row contains only ref bases """
        
        if self.alt == "X":
            return True
        else:
            return False
    
    def get_altbase(self):
        """ returns the bases that dont correspond to the ref """
        
        allbases = self.alt.split(",")
        altbases = [x for x in allbases if x != "X"]
        
        return altbases
    
    def is_snp(self):
        """ test if the alt field is anything other than "X" """
        
        if self.alt == "X":
            return False
        else:
            return True
    
    

############################################################################################################################################################################################################################################################################################################################################################
##############################################################################################################################################################################

class AltBaseType():
    
    """
    This class collects informations about one type (A/C/T/G/ref) of base in a RawBaseRow
    """
    
    def __init__(self, letter, bases, quals):
        
        self.numquals = []
        
        #print "parsing", letter, bases, quals
        
        if letter == "ref":
            self.upper = "."
            self.lower = ","
        else:   
            self.upper = letter.lower()
            self.lower = letter.upper()
        
        for i in range(0,len(bases)):
            base = bases[i]
            if base == self.upper or base == self.lower:
                if len(bases) == len(quals):
                    qual = quals[i]
                    self.numquals.append(qual)
                
        self.worstqual = None
        self.meanqual = None
        len_qual = 1
        
        if len(self.numquals) != 0:
            len_qual = len(self.numquals)       
            self.meanqual = sum(self.numquals)/float(len_qual)
            self.worstqual = min(self.numquals)
        
        #print letter, "meanqual", self.meanqual, "worst", self.worstqual
        
        self.up_count = bases.count(self.upper)
        self.low_count =  bases.count(self.lower)        
        self.total_count = self.up_count + self.low_count
    

class RawBaseRow(BasicRowParser):
      
    """
    Each instance corresponds to one site from samtools mpileup without SNP calling or bcftools
    Contig0 31 T    11    G.........^#.     ,>HFGHHGBHH
    """
        
    def __init__(self, row):
        
        BasicRowParser.__init__(self, row)
        self.ref = self.fields[2]
        self.qual = int(self.fields[3])
        
        self.bases = self.fields[4]
        self.basequals = self.fields[5]
        self.allnumquals = [(ord(x) - 33) for x in self.basequals]
        
        valid_baseletters = ["A","a", "C","c", "T","t", "G","g", ".", ",", "*", "N", "n", "+", "-"]
        bases = [base for base in self.bases if base in valid_baseletters]
        
        ##############################################################
        
        self.A = AltBaseType("A", bases, self.allnumquals)
        self.C = AltBaseType("C", bases, self.allnumquals)
        self.T = AltBaseType("T", bases, self.allnumquals)
        self.G = AltBaseType("G", bases, self.allnumquals)
        self.N = AltBaseType("N", bases, self.allnumquals)
        self.ref = AltBaseType("ref", bases, self.allnumquals)
        
        ##############################################################

        
        Abase = bases.count("A") + bases.count("a") 
        Cbase = bases.count("C") + bases.count("c") 
        Tbase = bases.count("T") + bases.count("t") 
        Gbase = bases.count("G") + bases.count("g")
        Nbase = bases.count("N") + bases.count("n")
        refbase = bases.count(".") + bases.count(",")
        indels = bases.count("+") + bases.count("-") + bases.count("*")
        
	self.coverage = float(sum([Abase, Cbase,Tbase,Gbase,refbase, Nbase]))
        total_count = float(sum([Abase, Cbase,Tbase,Gbase,refbase])) # should this include Ns or not?
	
        if total_count < 1:
            total_count = 1
        
        self.total_count = total_count
        self.ref_count = refbase
        
        self.alt_count = sum([Abase, Cbase, Tbase, Gbase])   
        
        self.indels = indels
        self.indel_size = None
        if self.indels != 0:
	    indel_re = re.search(("([+-]+)([\d]+)"), self.bases)

            if indel_re != None:
                
                self.indel_size = int(indel_re.groups()[1])
                self.alt_count -= self.indel_size
		self.total_count -= self.indel_size
    
	self.alt_percent = float(self.alt_count) / total_count
	
	if self.total_count < 1:
            self.total_count = 1 
	self.ref_percent = self.ref_count / self.total_count
	self.coverage = total_count
	
    ##################################################################################################
    # Base type methods
    
    def is_both_strands(self, base):
        
        cand = getattr(self, base)
        if cand.up_count > 0 and cand.low_count > 0:
            return True     
        else:
            return False
        
    def is_present(self, base):
        
        if base.total_count != 0:
            return True
        else:
            return False
        
    def cutoff_worst_qual(self, base, cutoff):
        
        if base.worstqual >= cutoff:
            return True
        else:
            return False
    
    ##################################################################################################
    
    def is_valid_celegans_denver_snp(self, base):
    
        """
        tests if the snp is valid according to the denver 2011 paper
        """
        
        if len(base) == 1:       
            cand = getattr(self, base.upper())
            
            if self.coverage >= 3 and self.coverage <= 25 and cand.worstqual >= 25 and self.is_both_strands(base): 
                return True
            else:
                return False
            
        else:
            return False
    
    ##################################################################################################
    
    def is_ref(self):
        """
        Tests if the position is likely to be a reference base.
        """
        
        if self.indels <= 1 and self.ref_count >= 5 and self.alt_count <= 2 and self.ref_percent > 0.8:
            return True        
        else:
            return False

    
    ##################################################################################################
    
    def ref_percent_cutoff(direction, cutoff):
        """
        Tests if the percentage of ref bases is (> min) or (< max)
        """  
        return cutoff_test(self.ref_percent, direction, cutoff)
            
    def ref_count_cutoff(direction, cutoff):
        """
        Tests if the number of ref bases is (> min) or (< max)
        """  
        return cutoff_test(self.ref_count, direction, cutoff)
        
    def alt_count_cutoff(direction, cutoff):
        """
        Tests if the number of alt bases is (> min) or (< max)
        """  
        return cutoff_test(self.alt_count, direction, cutoff)
        
    def alt_percent_cutoff(direction, cutoff):
        """
        Tests if the percentage of alt bases is (> min) or (< max)
        """  
        return cutoff_test(self.alt_percent, direction, cutoff)
    

##############################################################################################################################################################################
##############################################################################################################################################################################


class VCFrow(BasicRowParser):
    
    """
    Each instance corresponds to one variant site (= one row in a .vcf file)
    """
    
    def __init__(self, row, startpos=None):
        

        BasicRowParser.__init__(self, row, startpos)
        
        self.alt = self.fields[4]
        self.qual = float(self.fields[5])
        self.fq = None
        self.dp = None
        
    def parse_info(self):  
        
        self.info = self.fields[7]
        if "DP=" in self.info:
            dp = re.search(("(DP=)([\d]+)"), self.info)
            self.dp = int(dp.groups()[1])
            
        if "FQ=" in self.info:
            fq = re.search(("(FQ=)([-]?[\d]+)"), self.info)
            self.fq = int(fq.groups()[1])
    
    ##########################################################
    # type questions
    

    def is_indel(self):
        """ test if "INDEL" is in the info field """
        self.info = self.fields[7]
        
        if "INDEL" in self.info:
            return True
        else:
            return False
    
    def get_indel_type(self):
    
        if self.is_indel():
            
            if len(self.get_alt()) > len(self.get_ref()):
                return "insertion"
            
            else:
                return "deletion"
            
        else:
            return None
        
        
        
    def is_snp(self):
        """ test if the alt field is anything other than "." """
   
        if self.alt == ".":
            return False
        else:
            return True
    
    def is_hom(self, max_fq = None):
        """ test if the FQ value is lower than the given number (-> a hom SNP), default is -40 \n False for Indels. """
        
        max_fq = -40.0 if max_fq is None else float(max_fq)
        
        self.parse_info()
        
        if self.fq <= max_fq and not self.is_indel():
            return True
        else:
            return False
            
    def is_het(self, min_fq = None):
        """ test if the FQ value is higher than the given number (-> a het SNP), default is 40 \n False for Indels. """
        
        min_fq = 40.0 if min_fq is None else float(min_fq)
        
        min_fq = float(min_fq)
        self.parse_info()
        
        if self.fq >= min_fq and not self.is_indel():
            return True
        else:
            return False
        
           
    ##########################################################
    # cutoffs
    
    def min_qual(self, min_qual = None):
        
        min_qual = 20 if min_qual == None else float(min_qual)
        
        if self.qual >= min_qual:
            return True
        else:
            return False
    
    def min_coverage(self, min_cov):
        
        min_cov = int(min_cov)
        self.parse_info()
        
        if self.dp >= min_cov:
            return True
        else:
            return False
    
    def max_coverage(self, max_cov):
        
        max_cov = int(max_cov)
        self.parse_info()
        
        if self.dp <= max_cov:
            return True
        else:
            return False
    
    def max_fq(self, max_fq):
        self.parse_info()

        if self.fq <= max_fq:
            return True
        else:
            return False
    
    ##########################################################
    # get values directly
    
    def get_dp(self):
        
        self.parse_info()
        return self.dp

    def get_fq(self):
        
        self.parse_info()
        return self.fq
                
                

##############################################################################################################################################################################
##############################################################################################################################################################################


class SNPeffectrow(BasicRowParser):
    
    """
    Each instance corresponds to one judged variant site (= output from judge_snps.3.1.py)
    
    Contig41        986094  T       A       Chrom4  intergenic
    Contig41        986091  T       C       Chrom4  intergenic
    Contig26        509649  C       T       Chrom1  intron  Contig26-snap.57
    Contig132       38932   C       T       Chrom5  intron  Contig132-snap.10
    Contig104       110245  G       T       Chrom4  exon    Contig104-snap.14       cel:-_PFAM:-    -       nonsynonymous   TGC->TGA        C->*    stop
    
    
    """
    
    def __init__(self, row, startpos=None):
        

        BasicRowParser.__init__(self, row, startpos)
        
        self.ref = self.fields[2]
        self.alt = self.fields[3]
        self.chrom = self.fields[4]
        self.location = self.fields[5]
        self.homolog = None
        self.snap = None
        
        if self.location == "exon":
            self.snap = self.fields[6]
            self.homolog = self.fields[7]
            
    
    def is_stop(self):
        if "stop" in row:
            return True
        else:
            return False
            

##############################################################################################################################################################################
##############################################################################################################################################################################

class CelegansSNPParser:
    
    """
    Parse the list of denovo SNP from the Bear 2009 PNAS paper
    
    e.g B574 I 4,441,240 A T IG
    """
    
    def __init__(self, row):

        row = row.rstrip("\n")
        
        self.fields = row.split("\t")
        
        self.chromosome = self.fields[1]
        self.chromosomename = "CHROMOSOME_" + self.fields[1]
        self.con = self.chromosomename
        
        self.pos = self.fields[2].replace(",", "")     
        self.pos = int(self.pos)
        
        self.ref = self.fields[3]
        self.alt = self.fields[4]
            