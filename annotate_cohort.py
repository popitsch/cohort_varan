#
#
# Copyright (C) 2020 XXX.  All rights reserved.
#
# This file is part of cohort_varan
# 
# See the file LICENSE for redistribution information.
#
# @author: niko.popitsch
"""
    Annotates a cohort of VCF files and writes TSV tables per pedigree
"""
from collections import OrderedDict
import copy
import csv, datetime, time, logging, sys, os, json
import gzip
import itertools
import math
from multiprocessing import Pool, Process, Manager
import warnings

import psutil
from tqdm import tqdm
import vcfpy

import pandas as pd
import pyranges as pr
from utils import check_file, files_exist, remove_file, bgzip, count_lines


if __name__ == '__main__':    
    pass

#===================== utility ===========================================

# @see https://stackoverflow.com/questions/2187269/print-only-the-message-on-warnings
def custom_formatwarning(msg, *args, **kwargs):
    """ a custom warning output format """
    # ignore everything except the message
    return str(msg) + '\n'
warnings.formatwarning = custom_formatwarning

# ==========================================================================
#    Variant filter tags
#    AF = maxPopAF too high
#    NC = no call in any affected sample
#    NI = no supported inheritance model (only if include_na_models is False)
#    LS = low score
#    ID = Incomplete data (CNVs only)
#    UT = unspecific variant type (CNVs only)
#    NP = no_pass (these variants are not included in the debug log!)
# ==========================================================================
filter_tags = ["AF", "NC", "NI", "LS", "ID", "UT", "NP" ]

# ==========================================================================
#    Scoring
# ==========================================================================
def getMax(values):
    """ get maximum value from passed variable. Supports int, float, string and lists of floats"""
    if values is None:
        return None
    if type(values) is int:
        return values
    if type(values) is float:
        return values
    if type(values) is str:
        return None
    if values is "NA":
        return None
    if not values:
        return None
    vmax = None
    for v in values:
        if v is not None and v is not "NA":
            if vmax is None:
                try:
                    vmax = float(v)
                except ValueError:
                    vmax = None
            else:
                try:
                    vmax = float(max(float(v), vmax)) 
                except ValueError:
                    vmax = None
    return vmax

def toFloat(s):
    """ convert to float. if a string containing comma-separated values is passes, the maximum value is chosen. """
    if s is None:
        return None
    if type(s) is float:
        return s
    if type(s) is int:
        return float(s)
    return getMax(s.split(","))

def fmt(v):
    """ format as string or return NA if None """
    if v is None:
        return "NA"
    return str(v)


def check_config( c ):
    """ check config for validity """
    # check existing sections
    for t in ["dataset_name", "input_data", "filters", "output_fields", "d_score_calc", "def_a_score", "so_term"]:
        if t not in c:
            logging.error("Error: config section ", t, "was not found! Exiting...")
            sys.exit(1)
    # check output_fields sections
    for t in ["included_info_fields", "included_info_factors", "included_info_flags"]:
        if t not in c["output_fields"]:
            logging.error("Error: config section ", t, "was not found! Exiting...")
            sys.exit(1)   
    # input data
    for t in ["pedigrees", "snv_vcf", "ped_pattern", "gado_pattern","exomiser_pattern", "reg_db", "gene_gff"]:
        if t not in c["input_data"]:
            logging.error("Error: input data section ", t, "was not found! Exiting...")
            sys.exit(1)    
    if len(c["input_data"]["pedigrees"])<=0:
            logging.error("Error: no pedigree id configured! Exiting...")
            sys.exit(1)  
    # check for existing input files
    for t in c["input_data"]:
        if t in ["snv_vcf","alias_table","gene_anno_table","reg_db","gene_anno"]:
            check_file(c["input_data"][t])


def calc_score(rec, conf, var_type, def_value):
    """ calculate d_score """
    scores=[]
    for i,f in enumerate(conf["d_score_calc"][var_type]["fields"]):
        if f.startswith("score_"):
            score = calc_score(rec, conf, f[6:], def_value)
        else:
            score = rec.INFO.get(f, None) 
        if score is not None:
            scores+=[toFloat(score) / conf["d_score_calc"][var_type]["norm"][i]]
    if len(scores)==0:
        return def_value
    #summarize
    if conf["d_score_calc"][var_type]["summarisation"] == "max":
        return getMax(scores+[def_value])
    elif conf["d_score_calc"][var_type]["summarisation"] == "mean":
        ssum=0
        for s in scores:
            ssum+=s 
        return float(ssum)/float(len(scores))
    sys.exit("Unknown summarisation method " + conf["d_score_calc"][var_type]["summarisation"])

# ==========================================================================
#    Util
# ==========================================================================
def serialize(rec):
    """ for use in map() operations """
    return rec.serialize()

def show_prog(q, max_value):
    """ show progress bar """
    prog = tqdm(total=max_value, unit=' vars', desc='Analyzing variants')
    while 1:
        try:
            to_add = q.get(timeout=1)
            prog.n += to_add
            prog.update(0)
            prog.set_description("Current mem: %i%%" % (psutil.virtual_memory().percent)) 
            q.task_done()
            if prog.n >= max_value:
                break
        except:
            continue 
# ==========================================================================
#    SO term
# ==========================================================================
class SoTerm():
    """ Represents a sequence ontology term """
    def __init__(self, term, d_score, a_score, ttype, so_id, description):
        self.term = term
        self.d_score = d_score
        self.a_score = a_score
        self.type = ttype
        self.so_id = so_id
        self.description = description
        
    def __repr__(self):
        return self.__str__()  
    def __str__(self):
        ret = ("%s [d:%f, a:%f]" % (self.term, self.d_score, self.a_score) )
        return (ret)   


# ==========================================================================
#    Pedigree data model
# ==========================================================================
class Person():
    """ Represents a person """
    def __init__(self, family, pid, sex, affected, mum, dad):
        self.family = family
        self.id = pid
        self.sex = sex
        self.affected = affected
        self.mum = mum
        self.dad = dad
        self.children = []
    def is_parent(self):
        return len(self.children) > 0
    def has_parent(self):
        return self.dad or self.mum
    def __repr__(self):
        return self.__str__()  
    def __str__(self):
        ret = ("[%s" % self.id)
        if self.mum:
            ret += " ,mum=%s" % self.mum.id
        if self.dad:
            ret += " ,dad=%s" % self.dad.id
        ret += "]"
        return (ret)   

    
class Family():
    """ Represents a family """
    def __init__(self, pid):
        self.id = pid
        self.members = OrderedDict()
    def affected_members(self):
        return [x for x in self.members.values() if x.affected == True]
    def unaffected_members(self):
        return [x for x in self.members.values() if x.affected == False]
    def __repr__(self):
        return self.__str__()  
    def __str__(self):
        return ('[%s: %s]' % (self.id, ", ".join(str(x) for x in self.members.values())))
    
# exomiser inheritance models
EXOMISER_IM = ["AR", "AD", "XR", "XD"]

    
class Pedigree():
    """ Represents a pedigree """
    def __init__(self, pedF):
        super().__init__()
        self.name = os.path.splitext(os.path.basename(pedF))[0] # filename w/o ext
        self.families = OrderedDict()
        self.readFromFile(pedF)
    def readFromFile(self, pedF):
        try:
            with open(pedF, 'r') as csvfile:
                reader = csv.reader(csvfile, delimiter='\t')
                for row in reader:
                    fid = row[0]
                    pid = row[1]
                    sex = row[4]
                    aff = row[5] == '2'
                    f = self.families.get(fid, Family(fid))       
                    p = Person(f, pid, sex, aff, None, None)
                    f.members[pid] = p
                    self.families[fid] = f
        except IndexError:
            sys.exit("Error parsing pedigree " + pedF)
        # update parent/child        
        with open(pedF, 'r') as csvfile:
            reader = csv.reader(csvfile, delimiter='\t')
            for row in reader:
                fid = row[0]
                pid = row[1]
                did = row[2]
                mid = row[3]
                f = self.families.get(fid)       
                p = f.members.get(pid)
                if mid in f.members.keys():
                    mum = f.members.get(mid)
                    p.mum = mum
                    mum.children.append(p)
                if did in f.members.keys():
                    dad = f.members.get(did)
                    p.dad = dad
                    dad.children.append(p)
    def all_ids(self):
        ret = []
        for f in self.families.values():
            for pid in f.members.keys():
                ret = ret + [pid]
        return (ret) 
    def affected_ids(self):
        ret = []
        for f in self.families.values():
            for m in f.affected_members():
                ret = ret + [m.id]
        return (ret) 
    def get_person(self, pid):
        for f in self.families.values():
            if pid in f.members.keys():
                return f.members[pid]
        return None
    def calc_inh_support(self, rec, min_gq, min_gq_dnm, min_aq_dnm=0):
        """ calculate support (=#individuals/GT that support the model-#individuals/GT that contradict the model) for all supported inheritance models.  
            The following rules are applied:
            * recessive: 
               -1    for all unaffected sample containing a HOM call
               +1    for all affected samples with HOM calls that are inherited from mum & dad 
               0     for all other samples
            * dominant:
               -1    for all unaffected samples with a HOM/HET call
               -1    for all affected samples containing a HOM call 
               +1    for all affected samples with a HET call that was inherited from an affected sample if available
               0     for all other samples
            * de_novo:
               -1    for all unaffected samples with a HET/HOM call
               -1    for all affected samples that inherit a call 
               +1    for all affected samples containing a HET call with high GQ that was not inherited from mum&dad 
                0    for all other samples
            NOTE that for missing data (including calls with low GQ) we assume that GT that supports the respective inheritance model
            independent of each other which may lead to the situation that different genotypes for the alleles are assumed per inh model. """
        sup_dom = 0
        sup_rec = 0
        sup_dnm = 0
        #sup_com = 0
        genotypes=[]
        genotype_quals=[]
        num_hq_calls = 0
        hom_aff = 0
        hom_unaff = 0
        het_aff = 0
        het_unaff = 0
        sup_comhet_mum=[]
        sup_comhet_dad=[]
        for pid in self.all_ids():
            p = self.get_person(pid)
            if p is None:
                print("Unknown pid %s [%s] " % (pid, self ) )
                continue
            if pid not in rec.call_for_sample: # no data for this sample
                gt = None
                gq = None
            else:
                call = rec.call_for_sample[pid]
                # get genotypes
                gt = call.gt_type # hom_ref = 0, het = 1, hom_alt = 2 (which alt is untracked), uncalled = None
                gq = getMax(call.data["GQ"])
            # calculate VAF if possible (FIXME: for 1st allele only)
            # vaf = ad[1] / dp if ( ad is not None and len(ad)>0 and dp is not None and dp>0 ) else None
            genotypes.append(str(gt) if gt is not None else "NA")
            genotype_quals.append(str(gq) if gq is not None else "NA")
            # get m/f genotypes if any
            good_gq = gq is not None and gq >= min_gq
            gt_m = rec.call_for_sample[p.mum.id].gt_type if ((p.mum is not None) and (p.mum.id in rec.call_for_sample)) else None
            gt_d = rec.call_for_sample[p.dad.id].gt_type if ((p.dad is not None) and (p.dad.id in rec.call_for_sample)) else None
            # check whether inherited/affected
            is_inherited = ( gt_m is not None and gt_m > 0 ) or ( gt_d is not None and gt_d > 0 )
            is_affected = pid in self.affected_ids()
            # skip calls with low GQ
            if not good_gq or gt is None:
                continue 
            # count calls
            num_hq_calls +=1
            if gt == 1:
                if is_affected:
                    het_aff += 1
                else:
                    het_unaff += 1
            if gt == 2:
                if is_affected:
                    hom_aff += 1
                else:
                    hom_unaff += 1
            # calc support for comphet (only if complete genotypes!)
            if gt == 1 and is_affected:
                gq_m = getMax(rec.call_for_sample[p.mum.id].data["GQ"] if p.mum and p.mum.id in rec.call_for_sample else None)
                gq_d = getMax(rec.call_for_sample[p.dad.id].data["GQ"] if p.dad and p.dad.id in rec.call_for_sample else None)
                if gq_m is not None and gq_m >= min_gq and gq_d is not None and gq_d >= min_gq:               
                    if gt_m == 1 and gt_d == 0:
                        sup_comhet_mum.append(pid)               
                    if gt_m == 0 and gt_d == 1:
                        sup_comhet_dad.append(pid)                                  
            # calc supp
            if is_affected:
                if gt==0: # homref. reduces support for dom/recessive inheritance
                    sup_dnm += 0
                    sup_dom -= 1
                    sup_rec -= 1
                elif gt==1: # het
                    sup_rec -= 1
                    if gt_m is not None and gt_d is not None:
                        # check gq of parent calls
                        gq_m = getMax(rec.call_for_sample[p.mum.id].data["GQ"] if p.mum and p.mum.id in rec.call_for_sample else None)
                        gq_d = getMax(rec.call_for_sample[p.dad.id].data["GQ"] if p.dad and p.dad.id in rec.call_for_sample else None)
                        if gq_m is not None and gq_m >= min_gq_dnm and gq_d is not None and gq_d >= min_gq_dnm and gq >= min_gq_dnm:
                            if gt_m == 0 and gt_d == 0:
                                # get additional metadata for DNM filtering
                                #dp = call.data["DP"] if "DP" in call.data else None
                                ad = call.data["AD"] if "AD" in call.data else None
                                ab = getMax(call.data["AB"]) if "AB" in call.data else None
                                aq = getMax(rec.INFO["AQ"]) if "AQ" in rec.INFO else None
                                filtered=False
                                # check for minimum AD of alt-allele if AD field set
                                if ad is not None and len(ad)>0 and ad[1]<4:
                                    filtered=True
                                # check for allelic balance
                                if ab is not None and ( ab < 0.2 or ab > 0.8):
                                    filtered=True
                                # check for allele quality
                                if aq is not None and aq < min_aq_dnm:
                                    filtered=True
                                if not filtered:
                                    sup_dnm += 1 # supports DNM
                            else:
                                sup_dnm -= 1 # inherited
                                
                            if ( gt_m > 0 and p.mum not in self.affected_ids() ) or ( gt_d > 0 and p.dad not in self.affected_ids() ):
                                sup_dom -= 1
                            elif ( gt_m > 0 and p.mum in self.affected_ids() ) or ( gt_d > 0 and p.dad in self.affected_ids() ):
                                sup_dom +=1
                elif gt==2: # hom
                    sup_dom -= 1
                    if is_inherited:
                        sup_dnm -= 1 
                    if ( gt_m is None or gt_m ==1 ) and ( gt_d is None or gt_d ==1 ): # hom inherited from mum&dad
                        sup_rec += 1
                    else:
                        sup_rec -= 1
            else:
                if gt==0: # homref
                    pass
                elif gt==1: # het
                    sup_dnm -= 1
                    sup_dom -= 1
                elif gt==2: # hom
                    sup_dnm -= 1
                    sup_rec -= 1
                    sup_dom -= 1
        gt_str = "\t".join(genotypes)
        gq_str = "\t".join(genotype_quals)
        return [sup_rec, sup_dom, sup_dnm, num_hq_calls,hom_aff,hom_unaff,het_aff,het_unaff, gt_str, gq_str, sup_comhet_mum, sup_comhet_dad]
    def __repr__(self):
        return self.__str__()  
    def __str__(self):
        return ('ped %s [%s]' % (self.name, ", ".join(str(x) for x in self.families.values())))

    
def add_link(ped2a2b2ids, ped_name, a, b, pid):
    """ add and entry to the passed dict """
    if ped_name not in ped2a2b2ids:
        ped2a2b2ids[ped_name]=OrderedDict()
    if a not in ped2a2b2ids[ped_name]:
        ped2a2b2ids[ped_name][a]=OrderedDict()
    if b not in ped2a2b2ids[ped_name][a]:
        ped2a2b2ids[ped_name][a][b]=[]
    ped2a2b2ids[ped_name][a][b].append(pid)
    return ped2a2b2ids

def add_link2(ped2gid2a2b2ids, ped_name, gid, a, b, pid):
    """ add and entry to the passed dict """
    if ped_name not in ped2gid2a2b2ids:
        ped2gid2a2b2ids[ped_name]=OrderedDict()
    if gid not in ped2gid2a2b2ids[ped_name]:
        ped2gid2a2b2ids[ped_name][gid]=OrderedDict()
    if a not in ped2gid2a2b2ids[ped_name][gid]:
        ped2gid2a2b2ids[ped_name][gid][a]=OrderedDict()
    if b not in ped2gid2a2b2ids[ped_name][gid][a]:
        ped2gid2a2b2ids[ped_name][gid][a][b]=[]
    ped2gid2a2b2ids[ped_name][gid][a][b].append(pid)
    return ped2gid2a2b2ids

def calc_region_coverage(a, b):
    """ calculates what fraction of a is covered by b. Works only for single pyranges intervals """
    if len(a.index)!=1 or len(b.index)!=1:
        raise Exception("multiple pr objects passed")
    if a.Chromosome[0] != b.Chromosome[0]:
        return 0.0
    len_a=a.End-a.Start+1
    d=a.subtract(b)
    if len(d)==0:
        return 1.0
    len_d=d.End-d.Start+1
    return (1-len_d/len_a)

def get_output_rec(inrec, include_info_fields, rfilter, out_id, info_keep=["SVTYPE", "SVLEN", "POS", "END"], info_add=OrderedDict()):
    """ return deepcopy of output record """
    rec=copy.deepcopy(inrec)
    if rfilter:
        rec.add_filter(rfilter)
    if not include_info_fields:
        for f in info_keep:
            if f in rec.INFO:
                info_add[f]=rec.INFO[f]
        rec.INFO.clear()
    if out_id:
        rec.INFO["var2reg_id"]=out_id
    if info_add:
        rec.INFO.update(info_add)
    return rec

def process_chrom(args):
    """ process individual chrom"""
    
    chrom,  snvF, cnvF, peds, gados, exos, gene_anno, gene_anno_columns, max_var_per_cat, include_na_models, chunksize, debug, outdir, queue, config = args 
    global alias2gid, terms, exon_int, utr_int, gene_int, id2reg_regions, reg_regions
    logging.info("Parsing %s SNVs from %s and CNVs from %s" % (chrom, snvF, cnvF))
    tqdm.write("Parsing %s SNVs from %s and CNVs from %s" % (chrom, snvF, cnvF))
     
    # for storing vars per inheritance model
    ped2gid2mod2vid=OrderedDict()
    # for storing comphet candidates
    ped2gid2vid2mumordad2pid=OrderedDict()
    # SO terms that were ignored
    ignored_SO_terms = set()
    # for keeping track of current variant
    wrote2results=OrderedDict() 
    # variant and record ids
    var_id_counter=0
    rec_id_counter=0
    known_rec_id_counter=0

    # output records
    records=OrderedDict()
    known_records=OrderedDict()
    for ped in peds:
        records[ped.name] = list()
        known_records[ped.name] = list()
        
    # output files
    recF=OrderedDict()
    recOut=OrderedDict()
    knownRecF=OrderedDict()
    knownRecOut=OrderedDict()
    for ped in peds:
        recF[ped.name] = outdir+"."+chrom+"." + ped.name + ".rec_tmp"
        recOut[ped.name] = open(recF[ped.name], 'w')
        knownRecF[ped.name] = outdir+"."+chrom+"." + ped.name + ".known_rec_tmp"
        knownRecOut[ped.name] = open(knownRecF[ped.name], 'w')
    logF=outdir+"."+chrom+".log_tmp"
    
    # filter counts
    filter_counts=OrderedDict()
    total_var_count=OrderedDict()
    for vcf_type in ["SNV", "CNV"]:
        filter_counts[vcf_type]=OrderedDict()
        total_var_count[vcf_type]=0
        for tag in filter_tags:
            filter_counts[vcf_type][tag]=0
    
    # iterate SNVs and CNVs
    vcf_types = ["SNV", "CNV"] if cnvF is not None else ["SNV"]
    # read samples we have data for from VCF header
    reader = vcfpy.Reader.from_path(snvF)
    reader2 = vcfpy.Reader.from_path(cnvF) if cnvF else None
    with open(logF, 'w') as log:
        for vcf_type in vcf_types:

            # for warning about slow pyranges queries
            did_warn_pyranges=False
        
            # load data
            infile = snvF if vcf_type == "SNV" else cnvF
            if not infile:
                tqdm.write("No VCF file configured for type %s." % infile )
                continue # no CNV file configured
            
            # read VCF header via vcfpy (# full file: 13Mio SNVs)
            with warnings.catch_warnings():
                warnings.simplefilter("ignore") # supress warnings from incomplete headers
                header = reader.header if vcf_type=="SNV" else reader2.header
                parser = vcfpy.parser.RecordParser(header, header.samples)
                if vcf_type=="CNV":
                        # fix issue with missing CN header in CNV vcf
                        header.add_format_line(vcfpy.OrderedDict([('ID', 'CN'), ('Description', 'Copynumber'), ('Number', '1'), ('Type', 'Float')]))
    
            # iterate variants
            names = ['CHROM', 'POS', 'ID', 'REF', 'ALT',
                       'QUAL', 'FILTER', 'INFO', 'FORMAT'] 
            dtype = {'CHROM': str, 'POS': int, 'ID': str, 'REF': str, 'ALT': str,
               'QUAL': str, 'FILTER': str, 'INFO': str, 'FORMAT':str}
            for s in header.samples.names:
                names = names + [s]
                dtype[s] = str
        
            var_calls = 0
            for d in pd.read_csv(infile,
                                 delimiter='\t',
                                 encoding='utf-8',
                                 header=None,
                                 names=names,
                                 dtype=dtype,
                                 comment="#",
                                 float_precision='high',
                                 chunksize=chunksize,
                                 quoting=csv.QUOTE_NONE,
                                 error_bad_lines=False): 
                # max var per cat
                if max_var_per_cat and var_calls >= max_var_per_cat:
                    break
            
                # iterate only pass variants on given chrom
                all_var = len(d.index)
                total_var_count[vcf_type]+=all_var
                d = d.query("(FILTER=='PASS' | FILTER=='.') & CHROM=='"+chrom+"'")
                filter_counts[vcf_type]["NP"]+=all_var-len(d.index)
                 
                for _, row in d.iterrows():
                    var_id_counter += 1
                    var_calls += 1
                    # progress bar
                    queue.put(1)
                    
                    # max var per cat
                    if max_var_per_cat and var_calls >= max_var_per_cat:
                        break
                    
                    is_filtered = False
                                
                    # parse VCF record
                    try:
                        row['INFO']=row['INFO'].replace(';;',';') # fix vcfpy problem with empty info entries (';;')
                        rec = parser.parse_line("\t".join(map(str, row)))
                    except:
                        e = sys.exc_info()[0]
                        logging.error(e)
                        logging.error("Cannot parse VCF record. Skipping...")
                        logging.error(row)
                        filter_counts[vcf_type]["ID"]+=1
                        is_filtered = True
                        continue
                        
                    # if this is a known variant, write in any case to known_vars table
                    is_known = len( rec.ID ) > 0
                        
                    # get pos/len
                    altstr=','.join(map(serialize, rec.ALT))
                    alt1=altstr.split(',')[0]
                    pos, end=None, None
                    if vcf_type == "SNV":
                        pos = rec.POS
                        end = pos + len(alt1) - 1 
                    elif vcf_type=="CNV":
                        if ("POS" in rec.INFO) and ("END" in rec.INFO):
                            pos = int(rec.INFO.get("POS"))
                            end = int(rec.INFO.get("END"))
                    if not pos or not end:
                        # INCOMPLETE_DATA. ignore.
                        if debug:
                            logging.error("ID\t%s:%s" % (rec.CHROM, rec.POS), file=log)
                        filter_counts[vcf_type]["ID"]+=1
                        is_filtered = True
                        # drop this variant in any case as it will cause the downstream code to break due to incomplete info
                        continue
                        
                    # get vartype
                    vartype = None
                    if vcf_type == "SNV":
                        if len(rec.REF)==len(alt1):
                            vartype = "SNV" if len(rec.REF)==1 else "MNV"
                        else:
                            vartype="INDEL"
                    elif vcf_type=="CNV":
                        # calculate variant type
                        is_unspecific=False
                        for a in rec.ALT:
                            if a.type == "BND":  # skip general breakpoint annotations as too unspecific 
                                is_unspecific = True
                            else:
                                vartype = a.value if a.value != "None" else a.type
                        if vartype is None or is_unspecific:
                            # UNSPECIFIC_TYPE. ignore
                            if debug:
                                logging.error("UT\t%s:%s" % (rec.CHROM, rec.POS), file=log)
                            filter_counts[vcf_type]["UT"]+=1
                            is_filtered = True
                            # drop this variant in any case as it will cause the downstream code to break due to incomplete info
                            continue
                        # this will map, e.g., DEL=>DEL, DEL:ME=>DEL:ME, DEL:ME:something=>DEL:ME
                        vartype = ":".join(vartype.split(":")[0:2])
                    
                    

                    # filter by max pop AF
                    #maxAF = calc_max_AF(rec, ["global_A1000G_AF", "global_ExAC_AF", "global_UK10K_AF", "global_gnomAD_AF"] if vcf_type =="SNV" else ["GD_POPMAX_AF", "1000g_max_AF"])
                    maxAF = calc_score(rec, config, "max_pop_af_snv", 0) if vcf_type =="SNV" else calc_score(rec, config,"max_pop_af_cnv",0)
                    if maxAF is not None and maxAF > config["filters"]["max_pop_af"]:
                        # TOO_COMMON. ignore
                        if debug:
                            logging.error("AF\t%s:%s" % (rec.CHROM, rec.POS), file=log)
                        filter_counts[vcf_type]["AF"]+=1
                        is_filtered = True
                        if not is_known:
                            continue
        
                    # cohort freq: FIXME (better impl for CNVS based on partial overlap?)
                    cohort_calls,cohort_het,cohort_hom = 0,0,0
                    for call in rec.calls:
                        gt = call.gt_type # hom_ref = 0, het = 1, hom_alt = 2 (which alt is untracked), uncalled = None
                        gq = call.data["GQ"]
                        min_gq = config["filters"]["min_gq"] if vcf_type == "SNV" else -1 
                        maxgq = getMax(gq)
                        if maxgq is None:
                            maxgq = 0    
                        if gt is not None and maxgq >= min_gq:
                            if gt == 1:
                                cohort_het +=1
                            elif gt == 2:
                                cohort_hom +=1
                            cohort_calls += 1
                    cohortAF = 0 if cohort_calls == 0 else float(cohort_het + cohort_hom * 2) / float(cohort_calls * 2) # allele frequency in cohort
                    
                    # iterate over all pedigrees
                    is_nocall = 0
                    is_ni = 0
                    is_low_score = 0
                    for idx, ped in enumerate(peds):
                        wrote2results[ped.name]=False
                    
                        # calc support for inh models
                        min_gq = config["filters"]["min_gq"] if vcf_type == "SNV" else -1 # GQ is not a reliable measure for CNVs
                        min_gq_dnm = (config["filters"]["min_gq_dnm"] if "min_gq_dnm" in config["filters"] else 20) if vcf_type == "SNV" else -1 # GQ is not a reliable measure for CNVs
                        sup_rec, sup_dom, sup_dnm, num_hq_calls, hom_aff, hom_unaff, het_aff, het_unaff, gt_str, gq_str, sup_comhet_mum, sup_comhet_dad = ped.calc_inh_support(rec, min_gq, min_gq_dnm )
                        max_sup = getMax([sup_rec, sup_dom, sup_dnm])
                        sup_comhet = len(sup_comhet_mum) + len(sup_comhet_dad)
                        num_affected_calls = hom_aff + het_aff
                        
                        # no call with GQ>thresh in affected samples. ignore.
                        if num_affected_calls == 0: 
                            is_nocall+=1
                            continue
                        # no support for any inh model or comphet. ignore.
                        if not include_na_models and sup_comhet == 0 and ( max_sup is None or max_sup <= 0):
                            is_ni+=1
                            continue
                        # collect INFO values   
                        infos=[]
                        for f in config["output_fields"]["included_info_fields"]:
                            infos.append(fmt(getMax(str(rec.INFO.get(f, "NA")).split(","))))
                        for f in config["output_fields"]["included_info_factors"]:
                            infos.append(str(rec.INFO.get(f, "NA")))
                        for f in config["output_fields"]["included_info_flags"]:
                            infos.append("1" if f in rec.INFO else "0")

                        #--------------------------------
                        # Find overlap with genic regions
                        #--------------------------------
                        # get term with max weight per gene
                        gene2term=OrderedDict()
                        gene2HGVSp=OrderedDict() 
                        if vcf_type == "SNV":
                            if "ANN" in rec.INFO:
                                for a in rec.INFO["ANN"]:
                                    # 'Allele | Annotation | Annotation_Impact | Gene_Name | Gene_ID | Feature_Type | Feature_ID | Transcript_BioType | Rank | HGVS.c | HGVS.p | cDNA.pos / cDNA.length | CDS.pos / CDS.length | AA.pos / AA.length | Distance | ERRORS / WARNINGS / INFO
                                    # ex: ANN=T|splice_donor_variant&intron_variant|HIGH|HDLBP|ENSG00000115677|transcript|ENST00000310931.8|protein_coding|14/27|c.1731+1G>A||||||,T|....
                                    for aa in rec.ALT:
                                        dat = a.split("|")
                                        gid=dat[3] # gene name
                                        gid = alias2gid[gid] if gid in alias2gid else gid
                                        #print(dat[0]+"/"+aa.serialize() + "-" + dat[1])
                                        if dat[0] == aa.serialize(): # check right alt-allele
                                            for imp in dat[1].split("&"): # iterate annotations
                                                if imp in terms:
                                                    if gid in gene2term:
                                                        # replace if more deleterious only (i.e., respect input order)
                                                        if terms[gene2term[gid]].d_score<terms[imp].d_score:
                                                            #print("%s < %s" %(gene2term[gid], imp))
                                                            gene2term[gid] = imp
                                                        #else:
                                                        #    print("%s >= %s" %(gene2term[gid], imp))
                                                    else:
                                                        gene2term[gid] = imp
                                                else:
                                                    ignored_SO_terms.add(imp)
                                            if dat[10]: # iterate HGVS.p
                                                hgvsp = gene2HGVSp[gid] if gid in gene2HGVSp else []
                                                hgvsp += [ (dat[6]+":" + dat[10]).replace(',', '_') ]
                                                gene2HGVSp[gid] = hgvsp
                            else:
                                warnings.warn("No ANN field found in SNV entry at %s:%i" % (rec.CHROM, pos) )
                        elif vcf_type == "CNV":
                                query = pr.PyRanges(chromosomes=[rec.CHROM], starts=[min(pos, end)], ends=[max(pos, end)])
                                join = query.join(exon_int)
                                # FIXME: select most severe per gene!
                                if 'gene_name' in join.columns:
                                    overlapping_exons = set(join.gene_name.values)   
                                    for gid in overlapping_exons:
                                        gene2term[gid]="exonic_sv" # must be in config file!
                                join = query.join(utr_int)
                                if 'gene_name' in join.columns:
                                    overlapping_utrs = set(join.gene_name.values)   
                                    for gid in overlapping_utrs:
                                        if gid not in gene2term:
                                            gene2term[gid]="utr_sv"
                                join = query.join(gene_int)
                                if 'gene_name' in join.columns:
                                    overlapping_introns = set(join.gene_name.values)   
                                    for gid in overlapping_introns:
                                        if gid not in gene2term:
                                            gene2term[gid]="intronic_sv" # we assume 'intron' here if not overlapping with exon/utr annotations of this gene.                     

                        # write known variants to own table
                        if is_known:
                            for gid in gene2term:
                                imp = gene2term[gid] # impact
                                term = terms[gene2term[gid]] # SO:term
                                aa_change = gene2HGVSp[gid] if gid in gene2HGVSp else ["NA"]
                                known_rec_id = chrom + ":" + str(known_rec_id_counter)
                                var_id = chrom + ":" + str(var_id_counter)
                                d_score = term.d_score # default
                                a_score = term.a_score
                                known_records[ped.name].append(  [known_rec_id, var_id, gid, rec.CHROM, pos, end, rec.REF, 
                                        altstr, vartype, ",".join(aa_change), ",".join(rec.ID) if len(rec.ID)>0 else "NA", fmt(maxAF), fmt(cohortAF), fmt(cohort_calls), fmt(cohort_het), fmt(cohort_hom), 
                                        "NA", "anno", term.type, "NA", term.term, fmt(d_score), fmt(a_score), 
                                        sup_rec, sup_dom, sup_dnm, 
                                        ",".join(sup_comhet_mum) if sup_comhet_mum else "NA", ",".join(sup_comhet_dad) if sup_comhet_dad else "NA",  
                                        num_hq_calls, hom_aff, hom_unaff, het_aff, het_unaff,
                                        gt_str, gq_str, "\t".join(infos) ] )
                                known_rec_id_counter += 1
                            if is_filtered:
                                continue
                        
        
                        handled_gids=[]
                        for gid in gene2term:
                            imp = gene2term[gid] # impact
                            term = terms[gene2term[gid]] # SO:term
                            aa_change = gene2HGVSp[gid] if gid in gene2HGVSp else ["NA"]
                            d_score = term.d_score # default
                            a_score = term.a_score
                            # FIXME: how to calc these scores for SVs?
                            if vcf_type == "SNV":
                                if term.type == "intron":
                                    # intronic variant, calculate d_score from NC scores or spliceAI et al.
                                    d_score=calc_score(rec, config, "intron", d_score)
                                elif term.type == "splicing":
                                    # splicing region variant according to VEP (3bp in exon, 8bp in intron). calc weight from spliceAI et al.
                                    d_score=calc_score(rec, config, "splicing", d_score)
                                elif term.type == "utr":
                                    # UTR variant. Calculate from NC scores
                                    d_score=calc_score(rec, config, "utr", d_score)
                            
                            if vcf_type == "CNV" or d_score >= config["filters"]["min_dscore"]:
                                handled_gids += gid
                                
                                rec_id = chrom + ":" + str(rec_id_counter)
                                var_id = chrom + ":" + str(var_id_counter)
                                records[ped.name].append(  [rec_id, var_id, gid, rec.CHROM, pos, end, rec.REF, 
                                    altstr, vartype, ",".join(aa_change), ",".join(rec.ID) if len(rec.ID)>0 else "NA", fmt(maxAF), fmt(cohortAF), fmt(cohort_calls), fmt(cohort_het), fmt(cohort_hom),
                                    "NA", "anno", term.type, "NA", term.term, fmt(d_score), fmt(a_score), 
                                    sup_rec, sup_dom, sup_dnm, 
                                    ",".join(sup_comhet_mum) if sup_comhet_mum else "NA", ",".join(sup_comhet_dad) if sup_comhet_dad else "NA",  
                                    num_hq_calls, hom_aff, hom_unaff, het_aff, het_unaff,
                                    gt_str, gq_str, "\t".join(infos) ] )
                                rec_id_counter += 1
                                
                                # add gene<>inh_model<>var links
                                if sup_rec > 0:
                                    ped2gid2mod2vid = add_link(ped2gid2mod2vid, ped.name, gid, "recessive", rec_id)
                                if sup_dom > 0:
                                    ped2gid2mod2vid = add_link(ped2gid2mod2vid, ped.name, gid, "dominant", rec_id)
                                if sup_dnm > 0:
                                    ped2gid2mod2vid = add_link(ped2gid2mod2vid, ped.name, gid, "denovo", rec_id)
                                if include_na_models:
                                    ped2gid2mod2vid = add_link(ped2gid2mod2vid, ped.name, gid, "NA", rec_id)
            
                                # comphet. sup_comhet_mum, sup_comhet_dad contain the pids of affected samples that inherit a het from mum/dad
                                for pid in sup_comhet_mum:
                                    ped2gid2vid2mumordad2pid = add_link2(ped2gid2vid2mumordad2pid, ped.name, gid, rec_id, "mum", pid)
                                for pid in sup_comhet_dad:
                                    ped2gid2vid2mumordad2pid = add_link2(ped2gid2vid2mumordad2pid, ped.name, gid, rec_id, "dad", pid)
            
                                wrote2results[ped.name]=True
                            
                        #--------------------------------
                        # Find overlap with regulatory regions. Filter for dscore only for SNVs
                        #--------------------------------    
                        if vcf_type == "SNV":
                            d_score=calc_score(rec, config, "noncoding", 0)
                        else:
                            d_score=None
                        if vcf_type == "CNV" or d_score >= config["filters"]["min_dscore"]:
                            gid2reg=OrderedDict()
                            gid2assoctype=OrderedDict()
                            if vcf_type == "SNV" and "snv_reg_id_info_field" in config["filters"]:
                                # try to get reg_id values from INFO fields
                                if rec.INFO.get(config["filters"]["snv_reg_id_info_field"], None):
                                    for rid in rec.INFO.get(config["filters"]["snv_reg_id_info_field"]).split(","):
                                        if rid in id2reg_regions:
                                            region=id2reg_regions[rid]
                                            has_assoc = region['controlled_genes'] and len(str(region['controlled_genes']))>0
                                            gids = str(region['controlled_genes']).split(",") if has_assoc else str(region['closestGene_symbol']).split(",")
                                            assoc_type = "reg_db" if has_assoc else "closest_gene"
                                            for gid in gids:
                                                gid = alias2gid[gid] if gid in alias2gid else gid
                                                if gid not in gid2reg:
                                                    gid2reg[gid]=[]
                                                gid2reg[gid].append(region)
                                                gid2assoctype[gid] = assoc_type
                            elif vcf_type == "CNV" and "cnv_reg_id_info_field" in config["filters"]:
                                # try to get reg_id values from INFO fields
                                if rec.INFO.get(config["filters"]["cnv_reg_id_info_field"], None):
                                    for rid in rec.INFO.get(config["filters"]["cnv_reg_id_info_field"]).split(","):
                                        if rid in id2reg_regions:
                                            region=id2reg_regions[rid]
                                            has_assoc = region['controlled_genes'] and len(str(region['controlled_genes']))>0
                                            assoc_type = "reg_db" if has_assoc else "closest_gene"
                                            gids = str(region['controlled_genes']).split(",") if has_assoc else str(region['closestGene_symbol']).split(",")
                                            for gid in gids:
                                                gid = alias2gid[gid] if gid in alias2gid else gid
                                                if gid not in gid2reg:
                                                    gid2reg[gid]=[]
                                                gid2reg[gid].append(region)
                                                gid2assoctype[gid] = assoc_type
                            else:
                                # pyranges query
                                if not did_warn_pyranges:
                                    warnings.warn("Could not extract pre-annotated reg ids from VCF. Searching for regulatory regions via pyranges which is much slower!")
                                    did_warn_pyranges=True
                                query=pr.PyRanges(chromosomes=rec.CHROM, starts=[pos], ends=[end])
                                overlapping_reg_regions = reg_regions.join(query)
                                df=overlapping_reg_regions.df
                                if df is not None:
                                    for _, region in df.iterrows():
                                        has_assoc = region['controlled_genes'] and len(str(region['controlled_genes']))>0
                                        assoc_type = "reg_db" if has_assoc else "closest_gene"
                                        gids = str(region['controlled_genes']).split(",") if has_assoc else str(region['closestGene_symbol']).split(",")
                                        for gid in gids:
                                            gid = alias2gid[gid] if gid in alias2gid else gid
                                            if gid not in gid2reg:
                                                gid2reg[gid]=[]
                                            gid2reg[gid].append(region)
                                            gid2assoctype[gid] = assoc_type
                            #--------------------------------
                            # write one var per gene/region
                            #--------------------------------    
                            for gid in gid2reg.keys():
                                for region in gid2reg[gid]: 
                                    a_score = 0
                                    if region['db_source'] in config["def_a_score"]:
                                        # read default a_score from config
                                        a_score = config["def_a_score"][region['db_source']]
                                    rec_id = chrom + ":" + str(rec_id_counter)
                                    var_id = chrom + ":" + str(var_id_counter)
                                    records[ped.name].append(  
                                        [rec_id, var_id, gid, rec.CHROM, pos, end, rec.REF, 
                                         altstr, vartype, "NA",",".join(rec.ID) if len(rec.ID)>0 else "NA", fmt(maxAF), fmt(cohortAF), fmt(cohort_calls), fmt(cohort_het), fmt(cohort_hom),
                                         region['id'], region['db_source'], gid2assoctype[gid], region['closestGene_dist'], region['std_type']+"_variant", 
                                         fmt(d_score), fmt(a_score), 
                                         sup_rec, sup_dom, sup_dnm,
                                         ",".join(sup_comhet_mum) if sup_comhet_mum else "NA", ",".join(sup_comhet_dad) if sup_comhet_dad else "NA",
                                         num_hq_calls, hom_aff, hom_unaff, het_aff, het_unaff,
                                         gt_str, gq_str, "\t".join(infos)] )
                                    rec_id_counter += 1
                                        
                                    # add gene<>inh_model<>var links
                                    if sup_rec > 0:
                                        ped2gid2mod2vid = add_link(ped2gid2mod2vid, ped.name, gid, "recessive", rec_id)
                                    if sup_dom > 0:
                                        ped2gid2mod2vid = add_link(ped2gid2mod2vid, ped.name, gid, "dominant", rec_id)
                                    if sup_dnm > 0:
                                        ped2gid2mod2vid = add_link(ped2gid2mod2vid, ped.name, gid, "denovo", rec_id)
                                    if include_na_models:
                                        ped2gid2mod2vid = add_link(ped2gid2mod2vid, ped.name, gid, "NA", rec_id)
            
                                    # comphet. sup_comhet_mum, sup_comhet_dad contain the pids of affected samples that inherit a het from mum/dad
                                    for pid in sup_comhet_mum:
                                        ped2gid2vid2mumordad2pid = add_link2(ped2gid2vid2mumordad2pid, ped.name, gid, rec_id, "mum", pid)
                                    for pid in sup_comhet_dad:
                                        ped2gid2vid2mumordad2pid = add_link2(ped2gid2vid2mumordad2pid, ped.name, gid, rec_id, "dad", pid)
            
                                    wrote2results[ped.name]=True
                        # LOW_SCORE variants. ignore
                        if not wrote2results[ped.name]:
                            is_low_score+=1
                    
                    if is_filtered and is_known: # avoid double-counting as this variant was only kept until here because it was known       
                        continue
                    
                    if is_nocall == len(peds):
                        # NO_CALL. ignore
                        if debug:
                            logging.error("NC\t%s:%s" % (rec.CHROM, rec.POS), file=log)
                        filter_counts[vcf_type]["NC"]+=1
                    if is_ni == len(peds):
                        # NO_INH_MODEL. ignore
                        if debug:
                            logging.error("NI\t%s:%s" % (rec.CHROM, rec.POS), file=log)
                        filter_counts[vcf_type]["NI"]+=1
                    if is_low_score == len(peds):
                        if debug:
                            logging.error("LS\t%s:%s" % (rec.CHROM, rec.POS), file=log)
                        filter_counts[vcf_type]["LS"]+=1
                    
                           
                # write records to disc
                if len(records) > 100:
                    for ped_id in records.keys():
                        for rdata in records[ped_id]:
                            print( "\t".join( str(x) for x in rdata), file=recOut[ped_id])
                        records[ped_id]=list()
                if len(known_records) > 100:
                    for ped_id in known_records.keys():
                        for rdata in known_records[ped_id]:
                            print( "\t".join( str(x) for x in rdata), file=knownRecOut[ped_id])
                        known_records[ped_id]=list()
                          
        # write remaining records       
        for ped_id in records.keys():
            for rdata in records[ped_id]:
                print( "\t".join( str(x) for x in rdata), file=recOut[ped_id])
            records[ped_id]=list()    

        # write remaining known records       
        for ped_id in known_records.keys():
            for rdata in known_records[ped_id]:
                print( "\t".join( str(x) for x in rdata), file=knownRecOut[ped_id])
            known_records[ped_id]=list()    
                
    # close streams
    for ped in peds:
        recOut[ped.name].close()
        knownRecOut[ped.name].close()
    #print("Done iterating %s" % chrom) 
            
    # write gene tables
    gene_ids = OrderedDict() 
    records_genes = OrderedDict()
    gene_count=0
    for idx, ped in enumerate(peds):
        records_genes[ped.name] = list()
        gene_ids[ped.name]=0
        if ped.name in ped2gid2mod2vid:
            for gid in ped2gid2mod2vid[ped.name]:
                for ms in ped2gid2mod2vid[ped.name][gid]:
                    gado_z = str(gados[idx][gid]) if gid in gados[idx] else "NA"
                    exo_ps = str(exos[idx][gid]) if gid in exos[idx] else "NA"
                    annostr=""
                    if gene_anno_columns:
                        anno=gene_anno[gid] if gid in gene_anno else ["NA"] * len(gene_anno_columns)
                        annostr="\t".join(str(x) for x in anno)
                    records_genes[ped.name].append([gene_ids[ped.name], gid, gado_z, exo_ps, ms, len(ped2gid2mod2vid[ped.name][gid][ms]), ped2gid2mod2vid[ped.name][gid][ms], annostr])
                    gene_count+=1
                gene_ids[ped.name] += 1

    # write comphet table
    comphet_ids = OrderedDict() 
    records_comphet = OrderedDict()
    for idx, ped in enumerate(peds):
        records_comphet[ped.name] = list()
        comphet_ids[ped.name]=0
        if ped.name in ped2gid2vid2mumordad2pid:
            for gid in ped2gid2vid2mumordad2pid[ped.name]:
                vids = ped2gid2vid2mumordad2pid[ped.name][gid]
                for _, v1 in enumerate(vids):
                    for _, v2 in enumerate(vids):
                        if v1 < v2:
                            # samples that inherit v1 from mum 
                            mum1 = ped2gid2vid2mumordad2pid[ped.name][gid][v1]["mum"] if "mum" in ped2gid2vid2mumordad2pid[ped.name][gid][v1] else []
                            dad1 = ped2gid2vid2mumordad2pid[ped.name][gid][v1]["dad"] if "dad" in ped2gid2vid2mumordad2pid[ped.name][gid][v1] else []  
                            mum2 = ped2gid2vid2mumordad2pid[ped.name][gid][v2]["mum"] if "mum" in ped2gid2vid2mumordad2pid[ped.name][gid][v2] else []
                            dad2 = ped2gid2vid2mumordad2pid[ped.name][gid][v2]["dad"] if "dad" in ped2gid2vid2mumordad2pid[ped.name][gid][v2] else []
                            m1d2 = list(set(mum1) & set(dad2)) # samples for which there is a hit inherited from mum and dad
                            m2d1 = list(set(mum2) & set(dad1))  
                            candidate_pid = set(m1d2 + m2d1)
                            if candidate_pid:
                                records_comphet[ped.name].append([comphet_ids[ped.name], 
                                                          gid, 
                                                          v1,
                                                          v2,
                                                          len(candidate_pid), 
                                                          ",".join(candidate_pid) if candidate_pid else "NA",
                                                          ",".join(mum1) if mum1 else "NA",
                                                          ",".join(dad1) if dad1 else "NA",
                                                          ",".join(mum2) if mum2 else "NA",
                                                          ",".join(dad2) if dad2 else "NA"]
                                                           )
                                comphet_ids[ped.name] += 1
    # write to disc
    comphetF=outdir+"."+chrom+".comphet_tmp"
    comphet_count=0
    with open(comphetF, 'w') as out:
        for ped_id in records_comphet.keys():
            for rdata in records_comphet[ped_id]:
                print("%s\t%s" % (ped_id, "\t".join( str(x) for x in rdata)), file=out)
                comphet_count=comphet_count+1
    del records_comphet
    tqdm.write("Finished %s with %i/%i SNV/CNV calls, %i genes and %i comphets" % (chrom, total_var_count["SNV"], total_var_count["CNV"], gene_count, comphet_count) )
    logging.info("Finished %s with %i/%i SNV/CNV calls, %i genes and %i comphets" % (chrom, total_var_count["SNV"], total_var_count["CNV"], gene_count, comphet_count))
    return( chrom, recF, knownRecF, records_genes, comphetF, ignored_SO_terms, total_var_count, filter_counts )


def annotate_cohort(config_file, outdir, threads ):
    """ Main method that orchestrates the annotation pipeline """
    startTime = time.time()
    global alias2gid, terms, exon_int, utr_int, gene_int, id2reg_regions, reg_regions
    
    # ignore warnings (mainly from vcfpy)
    if not sys.warnoptions:
        warnings.filterwarnings(action="ignore", module="vcfpy")
    #============================================================================
    
    # load + check config
    config = json.load(open(config_file), object_pairs_hook=OrderedDict)
    check_config(config)
    
    # ensure dirs
    if not os.path.exists(outdir):
            logging.info("Creating dir " + outdir)
            os.makedirs(outdir)
    if not outdir.endswith("/"):
        outdir += "/"
    log_file = os.path.join(outdir, "annotate_cohort.log")

    # start 
    print("Logging output to %s" % log_file)
    for handler in logging.root.handlers[:]:
        logging.root.removeHandler(handler)
    logging.basicConfig(filename=log_file, level=logging.DEBUG)                
    logging.info("==========================================================")
    logging.info("annotate_cohort")
    logging.info("==========================================================")
    logging.info("Started script at %s.", str(datetime.date.today()))

    # get chunk size
    chunksize = config["input_data"]["chunksize"] if "chunksize" in config["input_data"] else 1000000
    logging.info("Setting chunksize to %i" % (chunksize))
    
    # debug mode?
    debug = False
    if "debug" in config and config["debug"]:
        debug = True
    
    # write effective config
    outF=outdir + config["dataset_name"] + ".v2r.effective_conf.json"
    with open(outF, 'w') as outs:
        print(json.dumps(config, indent=4, sort_keys=False), file=outs)
    
    chroms = ['chr'+str(x) for x in range(1,23)]+["chrX", "chrY"]
    if "chrom" in config:
        chroms = config["filters"]["chrom"]
    logging.info("Handling chromosomes %s" % (", ".join(chroms)))
    
    # read gene name alias table
    alias2gid=OrderedDict()
    if "alias_table" in config["input_data"]:
        d = pd.read_csv(config["input_data"]["alias_table"],delimiter='\t',encoding='utf-8')
        alias2gid=pd.Series(d.symbol.values,index=d.alias).to_dict()
    logging.info("Loaded %i gene symbol aliases" % (len(alias2gid)))
    
    # check model output
    max_var_per_cat = config["input_data"]["max_var_per_cat"] if "max_var_per_cat" in config["input_data"] else None
    include_na_models = config["input_data"]["include_na_models"] if "include_na_models" in config["input_data"] else True
    if not include_na_models:
        logging.info("===========================================================================")
        logging.info("NOTE that NA model entries will be filtered with the current configuration!")
        logging.info("===========================================================================")
    
    # read SO terms and weights
    logging.info("Reading SO terms + weights")
    terms=OrderedDict()
    for t in config["so_term"]:
        d_score=config["so_term"][t]["d_score"]
        a_score=config["so_term"][t]["a_score"]
        a_type=config["so_term"][t]["a_type"]
        so_id=config["so_term"][t]["so_id"]
        descr=config["so_term"][t]["descr"]
        terms[t] = SoTerm(t,d_score,a_score,a_type,so_id,descr)
    logging.info("\t%s" % (terms))
    
    # read gene annotations
    gene_anno=OrderedDict()
    gene_anno_columns=None
    if "gene_anno_table" in config["input_data"]:
        nline = sum(1 for _ in gzip.open(config["input_data"]["gene_anno_table"], 'rb'))
        logging.info("Accessing file with %i additional gene annotations" % (nline))
        with tqdm(total=nline, unit=' annotations', desc='Loading gene annotations') as bar:  
            for d in pd.read_csv(config["input_data"]["gene_anno_table"],
                                     delimiter='\t',
                                     encoding='utf-8',
                                     comment="#",
                                     float_precision='high',
                                     chunksize=chunksize,
                                     quoting=csv.QUOTE_NONE,
                                     index_col=False,
                                     error_bad_lines=False): 
                gene_anno_columns=d.columns.tolist()
                del gene_anno_columns[0]
                for _, row in d.iterrows():
                    l = row.tolist()
                    del l[0]
                    l = ["NA" if math.isnan(x) else x for x in l]
                    gene_anno[row[0]]=l
                    bar.update(1)
        logging.info("\tread %i gene annotations" % (len(gene_anno)))        
    else:
        logging.info("No gene annotations configured!")
    
    
    # =================================================================================
    # read all regulatory regions
    # =================================================================================
    nline = sum(1 for _ in gzip.open(config["input_data"]["reg_db"], 'rb'))
    logging.info("Accessing file with %i regulatory regions" % (nline))
    with tqdm(total=nline, unit=' reg_regions', desc='Loading regulatory regions annotations') as bar:  
        names = ['Chromosome', 'Start', 'End', 'id', 
                 'std_type', 'db_source', 'N_sources', 'N_methods', 'constrain_pct', 'PhyloP100_median', 
                 'closestGene_symbol', 'closestGene_dist', 'closestProt_symbol', 'closestProt_dist', 'controlled_genes',
                 'closestGene_TSS_symbol', 'closestGene_TSS_dist', 'closestProt_TSS_symbol', 'closestProt_TSS_dist'] 
    
        dtype = {'Chromosome': str, 'Start': int, 'End': int, 'id': str, 'std_type': str, 'db_source': str, 
                 'constrain_pct': float, 'PhyloP100_median': float,
                 'closestGene_symbol': str, 'closestGene_dist': str, 
                 'closestProt_symbol': str, 'closestProt_dist': int,
                 'controlled_genes': str, 'N_methods':int
                 }
        reg_regions=[]
        id2reg_regions=OrderedDict()
        for d in pd.read_csv(config["input_data"]["reg_db"],
                                 delimiter='\t',
                                 encoding='utf-8',
                                 header=None,
                                 names=names,
                                 #dtype=dtype,
                                 comment="#",
                                 float_precision='high',
                                 chunksize=chunksize,
                                 quoting=csv.QUOTE_NONE,
                                 index_col=False,
                                 low_memory=False,
                                 error_bad_lines=False,
                                 na_values=['NA'],
                                 keep_default_na=False): # this is required to interpret empty strings as NA instead of NaN 
            if chroms:
                d = d[d['Chromosome'].isin(chroms)]
            d['Start']+=1 # 0-based coords!
            reg_regions.append(pr.PyRanges(d))
            bar.update(len(d.index))
    logging.info("\tbuilding index") 
    reg_regions=pr.concat(reg_regions)
    for _, region in reg_regions.df.iterrows():
        id2reg_regions[region['id']]=region 
    logging.info("\tread %i regulatory regions" % (len(reg_regions)))        
    with pd.option_context('display.max_rows',5):
        logging.info(reg_regions.df)
        
    # =================================================================================
    # read gene annotations from GFF
    # =================================================================================
    logging.info("Reading gene annotation intervals")
    with tqdm(total=6, desc='Loading gene annotation intervals') as bar:  
        gff = pr.read_gff3(config["input_data"]["gene_gff"])
        bar.update(1)
        exon_int = gff[gff.Feature == 'exon']
        bar.update(1)
        utr3_int = gff[gff.Feature == 'three_prime_UTR']
        bar.update(1)
        utr5_int = gff[gff.Feature == 'five_prime_UTR']
        bar.update(1)
        utr_int = pr.concat([utr3_int, utr5_int])
        bar.update(1)
        gene_int = gff[gff.Feature == 'gene']
        bar.update(1)
        logging.info("Loaded %i exon, %i utr and %i gene annotations from %s" % (len(exon_int), len(utr_int), len(gene_int), config["input_data"]["gene_gff"]))
        
    # read samples we have data for from VCF header
    reader = vcfpy.Reader.from_path(config["input_data"]["snv_vcf"])
    #reader2 = vcfpy.Reader.from_path(config["input_data"]["cnv_vcf"]) if 'cnv_vcf' in config["input_data"] else None # TODO: for now we assume same sample names in SNV and CNV VCF!
    
    
    # =================================================================================
    # read input pedigrees
    # =================================================================================
    logging.info("Reading pedigrees and gene rank data")
    peds=[]
    gados=[]
    exos=[]
    for pid in tqdm(config["input_data"]["pedigrees"], desc='Loading pedigrees'):
        pedf = config["input_data"]["ped_pattern"].format(id=pid)
        gadof = config["input_data"]["gado_pattern"].format(id=pid)
        anyexo = False
        for im in EXOMISER_IM:
            exof = config["input_data"]["exomiser_pattern"].format(IM=im, id=pid)
            anyexo = anyexo or files_exist(exof)
        if files_exist([pedf, gadof]) and anyexo:
            ped = Pedigree(pedf)
            # do we have data for at least one sample?
            if len(list(set(ped.all_ids()) & set(reader.header.samples.names))) > 0:
                peds += [ped]
                gid2rnk=OrderedDict()
                d = pd.read_csv(gadof, delimiter='\t', encoding='utf-8')
                for _, row in d.iterrows():
                    gid = row['Hgnc']
                    gid = alias2gid[gid] if gid in alias2gid else gid
                    gid2rnk[gid]=row['Zscore']
                gados+=[gid2rnk]
                exo2rnk=OrderedDict()
                # read exomiser tables
                for im in EXOMISER_IM:
                    exof = config["input_data"]["exomiser_pattern"].format(IM=im, id=pid)
                    if files_exist([exof]):
                        d = pd.read_csv(exof, delimiter='\t', encoding='utf-8')
                        for _, row in d.iterrows():
                            exo2rnk[row['#GENE_SYMBOL']]=row['EXOMISER_GENE_PHENO_SCORE']
                exos+=[exo2rnk]       
                logging.info("\t%s: %s, %i gado_rnk, %i exo_ps" % (ped.name, ','.join(ped.all_ids()), len(gid2rnk), len(exo2rnk)))
            else:
                logging.info("\tDid not find any data for the following samples in the snv VCF: %s. Ignoring ped %s." % (",".join(ped.all_ids()), pid))
        else:
            logging.info("\tCannot load all required data for id %s. Some files (ped, gado or exomiser results) do not exist. Ignoring (%s, %s, %s)" % (pid, pedf, gadof, config["input_data"]["exomiser_pattern"]) )

    logging.info("\tDone. Read %i pedigrees and %i gado rank lists" % (len(peds), len(gados)))
    if len(peds) == 0:
        sys.exit("Check input configuration / data...")


    # =================================================================================
    # split VCFs
    # =================================================================================
    logging.info("Splitting input VCFs per chrom")
    
    def split_vcf_by_chrom(vcf_file, chroms, out_prefix, chunksize=10000, threads=1):
        # check whether split files already exist (if previous run failed/stopped) and if yes do not recreate
        fout=OrderedDict()
        split_exists=True
        total_var_count=0 # for counting # of vars
        for c in chroms:
            fout[c] = out_prefix+"."+str(c)+".vcf.gz"
            if not files_exist(fout[c]):
                split_exists=False
                break
            total_var_count+=count_lines(fout[c])      
        if split_exists:
            print("found existing per-chromosome VCF files! Will *not* recreate...")
            return(fout, total_var_count)     
            
        # read VCF header via vcfpy (# full file: 13Mio SNVs)
        reader = vcfpy.Reader.from_path(vcf_file)
        fout=OrderedDict()
        total_var_count=0 # for counting # of vars
        for c in chroms:
            fout[c] = out_prefix+"."+str(c)+".vcf"
            writer = vcfpy.Writer.from_path(fout[c], reader.header)
            writer.close()
        try:
            for d in pd.read_csv(vcf_file,
                                 delimiter='\t',
                                 encoding='utf-8',
                                 header=None,
                                 dtype=str,
                                 comment="#",
                                 float_precision='high',
                                 chunksize=chunksize,
                                 quoting=csv.QUOTE_NONE,
                                 error_bad_lines=False): 
                found_chr = set(d[0])
                for c in found_chr:
                    if c in fout:
                        total_var_count+=len(d.index)
                        fil = d[d.iloc[:,0]==c]
                        fil.to_csv(fout[c], sep='\t', mode='a', header=None, index=False)
            for c in chroms:
                bgzip(fout[c], index=True, delinFile=True, maxthreads=threads)
                fout[c] = fout[c]+".gz"
        except pd.io.common.EmptyDataError:
            logging.info("No data in VCF file %s" % vcf_file)
        return(fout, total_var_count)
    
    with tqdm(total=2, desc='Splitting input data by chromosome') as bar:
        snvF_chr, total_var_count_snv = split_vcf_by_chrom(config["input_data"]["snv_vcf"], chroms, outdir+"snv", chunksize, threads)
        bar.update(1)
        if "cnv_vcf" in config["input_data"]:
            cnvF_chr, total_var_count_cnv = split_vcf_by_chrom(config["input_data"]["cnv_vcf"], chroms, outdir+"cnv", chunksize, threads)  
        else: 
            cnvF_chr, total_var_count_cnv = None, 0
        bar.update(1)
    
    print("Starting main loop")
    # .........................................................................
    #                Parallel calling
    # .........................................................................
    try:
        queue = Manager().Queue()
        progress = Process(target=show_prog, args=(queue, total_var_count_snv + total_var_count_cnv))
        progress.start()
        # get input files
        snvF_chr_f=[]
        cnvF_chr_f=[]
        for c in chroms:
            snvF_chr_f+=[snvF_chr[c]]
            cnvF_chr_f+=[cnvF_chr[c]] if "cnv_vcf" in config["input_data"] else [None]
        param = zip(
            chroms,  
            snvF_chr_f, 
            cnvF_chr_f, 
            itertools.repeat(peds), 
            itertools.repeat(gados), 
            itertools.repeat(exos), 
            itertools.repeat(gene_anno), 
            itertools.repeat(gene_anno_columns), 
            itertools.repeat(max_var_per_cat), 
            itertools.repeat(include_na_models), 
            itertools.repeat(chunksize),
            itertools.repeat(debug),
            itertools.repeat(outdir),
            itertools.repeat(queue),
            itertools.repeat(config)
            )
        with Pool(processes=threads) as pool:
            results = pool.map(process_chrom, param)
            logging.info("receiving results...")
            progress.terminate()
            data=OrderedDict()
            for c in results:
                data[c[0]]= c
    except:
        #print("Terminating progress thread.")
        progress.terminate()
        raise       
        
         
    # =================================================================================
    # merge results and write result file
    # =================================================================================
    logging.info("Finished, writing results.")
    
    # =================================================================================
    # Load all variant record
    # =================================================================================
    logging.info("Loading output records.")
    records=OrderedDict()
    for c in chroms:
        records[c]=OrderedDict()
        chrom, recF, knownRecF, records_gene, comphetF, chr_ignored_SO_terms, var_count, filter_counts = data[c]
        for ped_id in recF.keys():
            if os.stat(recF[ped_id]).st_size > 0:
                d = pd.read_csv(recF[ped_id],delimiter='\t',encoding='utf-8', header=None, dtype=str, na_values=[''], keep_default_na=False)
                for _, row in d.iterrows():
                    if not ped_id in records[c]:
                        records[c][ped_id]=list()
                    dat=[str(x) for x in row]
                    records[c][ped_id].append(dat)
            remove_file(recF[ped_id])
    
    logging.info("Writing variant records.")    


    # if true, ids in output tables are pedigree-dependent! 
    id_scope_pedigree=config["id_scope_pedigree"] if "id_scope_pedigree" in config else False
    if id_scope_pedigree:
        print("NOTE: id scope is set to PEDIGREE, i.e., var and rec ids are only record-wide unique")

    # SO terms that were ignored
    ignored_SO_terms = set()
    # for mapping rec_ids
    rec_id_map = OrderedDict()
    var_id_map = OrderedDict()
    rec_ids = OrderedDict() 
    var_ids = OrderedDict() 
    for _, ped in enumerate(peds):
        outF=outdir + config["dataset_name"]+ "." + ped.name + ".v2r.vars.tsv"
        with open(outF, 'w') as out:
            rec_ids[ped.name] = 0
            rec_id_map[ped.name] = OrderedDict()
            var_ids[ped.name] = 0
            var_id_map[ped.name] = OrderedDict()
            
            # write header
            header  = "# Samples (Affected are indicated by [*]):\n"
            for pid in ped.all_ids():
                if pid in ped.affected_ids():
                    pid+=" [*]"
                header+="#    "+pid + "\n"        
            header += "# Genotype codes: hom_ref = 0, het = 1, hom_alt = 2 (which alt is untracked), uncalled = NA\n"
            header += "# Supported inheritance models are not filtered for NA genotypes or genotype quality\n"
            header += "rec_id\tvar_id\tgene\tchr\tstart\tend\tref\talt\tvar_type\taa_change\tknown_ids\tmax_pop_af\tcohort_af\tcohort_calls\tcohort_het\tcohort_hom\treg_id\tdb_source\treg_type\tclosest_gene_dist\tconsequence\td_score\ta_score\tsup_rec\tsup_dom\tsup_dnm\tsup_comphet_mum\tsup_comphet_dad\tnum_hq_calls\thom_aff\thom_unaff\thet_aff\thet_unaff"
            for pid in ped.all_ids():
                header+="\tGT_"+pid        
            for pid in ped.all_ids():
                header+="\tGQ_"+pid        
            header += ("\t" + "\t".join(config["output_fields"]["included_info_fields"])  ) if len(config["output_fields"]["included_info_fields"])>0 else ""
            header += ("\t" + "\t".join(config["output_fields"]["included_info_factors"]) ) if len(config["output_fields"]["included_info_factors"])>0 else ""
            header += ("\t" + "\t".join(config["output_fields"]["included_info_flags"])   ) if len(config["output_fields"]["included_info_flags"])>0 else ""
            print(header, file=out)
        
            for c in chroms:
                chrom, recF, knownRecF, records_gene, comphetF, chr_ignored_SO_terms, var_count, filter_counts = data[c]
                ignored_SO_terms.update(chr_ignored_SO_terms)
                if ped.name in records[c]:
                    for rec in records[c][ped.name]:
                        # record ids
                        if id_scope_pedigree:
                            rec_id_map[ped.name][rec[0]] = "r" + str(rec_ids[ped.name])
                            rec[0] = rec_id_map[ped.name][rec[0]]
                        rec_ids[ped.name] += 1
                        # set variant ids
                        if rec[1] in var_id_map[ped.name]:
                            if id_scope_pedigree:
                                rec[1] = var_id_map[ped.name][rec[1]]
                        else:
                            var_id_map[ped.name][rec[1]] =  "v" + str(var_ids[ped.name])
                            if id_scope_pedigree:
                                rec[1] = var_id_map[ped.name][rec[1]]
                            var_ids[ped.name] += 1
                       
                        print("\t".join( str(x) for x in rec), file=out)
    del records
    
    # =================================================================================
    # Load all known variant record
    # =================================================================================
    logging.info("Loading known records.")
    records=OrderedDict()
    for c in chroms:
        records[c]=OrderedDict()
        chrom, recF, knownRecF, records_gene, comphetF, chr_ignored_SO_terms, var_count, filter_counts = data[c]
        for ped_id in knownRecF.keys():
            if os.stat(knownRecF[ped_id]).st_size > 0:
                d = pd.read_csv(knownRecF[ped_id],delimiter='\t',encoding='utf-8', header=None, dtype=str, na_values=[''], keep_default_na=False)
                for _, row in d.iterrows():
                    if not ped_id in records[c]:
                        records[c][ped_id]=list()
                    dat=[str(x) for x in row]
                    records[c][ped_id].append(dat)
            remove_file(knownRecF[ped_id])
    
    logging.info("Writing known variant records.")     
    known_rec_ids = OrderedDict() 
    for _, ped in enumerate(peds):
        outF=outdir + config["dataset_name"]+ "." + ped.name + ".v2r.known_vars.tsv"
        with open(outF, 'w') as out:
            known_rec_ids[ped.name] = 0
    
            # write header
            header  = "# Samples (Affected are indicated by [*]):\n"
            for pid in ped.all_ids():
                if pid in ped.affected_ids():
                    pid+=" [*]"
                header+="#    "+pid + "\n"        
            header += "# Genotype codes: hom_ref = 0, het = 1, hom_alt = 2 (which alt is untracked), uncalled = NA\n"
            header += "# Supported inheritance models are not filtered for NA genotypes or genotype quality\n"
            header += "rec_id\tvar_id\tgene\tchr\tstart\tend\tref\talt\tvar_type\taa_change\tknown_ids\tmax_pop_af\tcohort_af\tcohort_calls\tcohort_het\tcohort_hom\treg_id\tdb_source\treg_type\tclosest_gene_dist\tconsequence\td_score\ta_score\tsup_rec\tsup_dom\tsup_dnm\tsup_comphet_mum\tsup_comphet_dad\tnum_hq_calls\thom_aff\thom_unaff\thet_aff\thet_unaff"
            for pid in ped.all_ids():
                header+="\tGT_"+pid        
            for pid in ped.all_ids():
                header+="\tGQ_"+pid        
            header += ("\t" + "\t".join(config["output_fields"]["included_info_fields"])  ) if len(config["output_fields"]["included_info_fields"])>0 else ""
            header += ("\t" + "\t".join(config["output_fields"]["included_info_factors"]) ) if len(config["output_fields"]["included_info_factors"])>0 else ""
            header += ("\t" + "\t".join(config["output_fields"]["included_info_flags"])   ) if len(config["output_fields"]["included_info_flags"])>0 else ""
            print(header, file=out)
        
            for c in chroms:
                if ped.name in records[c]:
                    for rec in records[c][ped.name]:
                        # record ids
                        if id_scope_pedigree:
                            rec[0] = "r" + str(known_rec_ids[ped.name])
                        known_rec_ids[ped.name] += 1
                        # variant ids
                        if rec[1] in var_id_map[ped.name]:
                            if id_scope_pedigree:
                                rec[1] = var_id_map[ped.name][rec[1]]
                        print("\t".join( str(x) for x in rec), file=out)
    del records
    
    
    # =================================================================================
    # write genes table
    # =================================================================================
    logging.info("Writing gene records.")
    gene_ids = OrderedDict() 
    unique_genes = OrderedDict()
    for _, ped in enumerate(peds):
        gene_ids[ped.name]=0
        unique_genes[ped.name]=set()
        outF=outdir + config["dataset_name"] + "." + ped.name + ".v2r.genes.tsv"
        with open(outF, 'w') as out:
            header = "# "+ped.name+"\n"
            header += "rec_id\tgene\tgado_zscore\texomiser_gene_pheno_score\tinh_model\tvariants_n\tvariants"
            if gene_anno_columns:
                header+="\t"+"\t".join(gene_anno_columns)
            print(header, file=out)
            for c in chroms:
                chrom, recF, knownRecF, records_gene, comphetF, chr_ignored_SO_terms, var_count, filter_counts = data[c]
                if ped.name in records_gene:
                    for rec in records_gene[ped.name]:
                        if id_scope_pedigree:
                            linked_vars = ",".join(rec_id_map[ped.name][str(x)] for x in rec[6])
                            rec[6] = linked_vars
                            rec[0] = "g" + str(gene_ids[ped.name])
                        print("\t".join( str(x) for x in rec), file=out)
                        gene_ids[ped.name] += 1  
                        unique_genes[ped.name].update([rec[1]])        
    # =================================================================================
    # write comphet table
    # =================================================================================
    logging.info("Writing comphet records.")
    records_comphet=OrderedDict()
    for c in chroms:
        records_comphet[c]=OrderedDict()
        chrom, recF, knownRecF, records_gene, comphetF, chr_ignored_SO_terms, var_count, filter_counts = data[c]
        if os.stat(comphetF).st_size > 0:
            d = pd.read_csv(comphetF,delimiter='\t',encoding='utf-8', header=None, dtype=str, na_values=[''], keep_default_na=False)
            for _, row in d.iterrows():
                ped_id=row[0]
                if not ped_id in records_comphet[c]:
                    records_comphet[c][ped_id]=list()
                dat=[str(x) for x in row]
                records_comphet[c][ped_id].append(dat[1:])
        remove_file(comphetF)
    # write merged 
    comphet_ids = OrderedDict() 
    for _, ped in enumerate(peds):
        comphet_ids[ped.name]=0
        outF=outdir + config["dataset_name"] + "." + ped.name + ".v2r.comphet.tsv"
        with open(outF, 'w') as out:
            header = "# "+ped.name+"\n"
            header += "rec_id\tgene\tv1\tv2\tnum_aff\tcandidate\tmum1\tdad1\tmum2\tdad2"
            print(header, file=out)
            for c in chroms:
                chrom, recF, knownRecF, records_gene, comphetF, chr_ignored_SO_terms, var_count, filter_counts = data[c]
                if ped.name in records_comphet[c]:
                    for rec in records_comphet[c][ped.name]:
                        if id_scope_pedigree:
                            rec[0] = "c" + str(comphet_ids[ped.name])
                            rec[2] = rec_id_map[ped.name][rec[2]]
                            rec[3] = rec_id_map[ped.name][rec[3]]
                        print("\t".join( str(x) for x in rec), file=out)
                        comphet_ids[ped.name] += 1 
    del records_comphet
    
    # write index
    logging.info("Writing index.")
    outF=outdir + config["dataset_name"] + ".v2r.idx.tsv"
    with open(outF, 'w') as out:
        print("# V2 analysis date: "+time.strftime('%Y-%m-%d %H:%M:%S'), file=out)
        print("rec_id\tpedigree\tgene_records\tunique_genes\tvariant_records\tvariants\tcomphet_records\tn_all_samples\tn_affected_samples\tall_samples\taffected_samples\tgene_file\tvariant_file\tknown_variant_file\tcomphet_file\tpedigree_file\tgado_file\texomiser_files", file=out)
        rec_id = 0
        for _, ped in enumerate(peds):
            genf = os.path.abspath(outdir + config["dataset_name"] + "." + ped.name + ".v2r.genes.tsv")
            bgzip(genf, index=False, delinFile=True, maxthreads=threads)
            genf+=".gz"
            
            varf = os.path.abspath(outdir + config["dataset_name"] + "." + ped.name + ".v2r.vars.tsv")
            bgzip(varf, index=False, delinFile=True, maxthreads=threads)
            varf+=".gz"
            # TABIX; skip first 4+#affected_samples lines. FIXME: table not sorted as it will report SNV then CNV per chrom.
            #tabix(varf, additionalParameters=["-s4", "-b5", "-e6", "-S"+str(4+len(ped.affected_ids()))])
    
            known_varf = os.path.abspath(outdir + config["dataset_name"] + "." + ped.name + ".v2r.known_vars.tsv")
            bgzip(known_varf, index=False, delinFile=True, maxthreads=threads)
            known_varf+=".gz"
            # TABIX; skip first 4+#affected_samples lines
            #tabix(known_varf, additionalParameters=["-s4", "-b5", "-e6", "-S"+str(4+len(ped.affected_ids()))])
            
            comf = os.path.abspath(outdir + config["dataset_name"] + "." + ped.name + ".v2r.comphet.tsv")
            bgzip(comf, index=False, delinFile=True, maxthreads=threads)
            comf+=".gz"
            
            pedf = os.path.abspath(config["input_data"]["ped_pattern"].format(id=ped.name))
            gadof = os.path.abspath(config["input_data"]["gado_pattern"].format(id=ped.name))
            exofs=[]
            for im in EXOMISER_IM:
                exof = config["input_data"]["exomiser_pattern"].format(IM=im, id=ped.name)
                if files_exist(exof):
                    exofs+=[exof]
            rec=[rec_id,
                ped.name, 
                gene_ids[ped.name], 
                len(unique_genes[ped.name]),
                rec_ids[ped.name],        
                var_ids[ped.name],         
                comphet_ids[ped.name],             
                len(ped.all_ids()), len(ped.affected_ids()), 
                ",".join(ped.all_ids()),",".join(ped.affected_ids()),
                genf, varf, known_varf, comf, pedf, gadof, ",".join(exofs)]
            print("\t".join( str(x) for x in rec), file=out)
            rec_id += 1
    bgzip(outF, index=True, delinFile=True, maxthreads=threads)
        
    
    
    # write filtering stats             
    outF=outdir + config["dataset_name"]+".v2r.stats.tsv"
    with open(outF, 'w') as log:
        # write stats
        print("#=============================================", file=log)
        print("# Filtered variants statistics", file=log)
        print("# AF = maxPopAF too high", file=log)
        print("# NC = no call in any affected sample", file=log)
        print("# LS = low score", file=log)
        print("# NI = no supported inheritance model (only if include_na_models is False)", file=log)
        print("# ID = Incomplete data (CNVs only)", file=log)
        print("# UT = unspecific variant type (CNVs only)", file=log)
        print("# NP = no_pass (these variants are not included in the debug log!)", file=log)
        print("=============================================", file=log)
        all_var_count=OrderedDict()
        all_filter_sum=OrderedDict()
        all_filter_counts=OrderedDict()
        for vcf_type in ["SNV", "CNV"]:
            all_filter_counts[vcf_type]=OrderedDict()
            all_var_count[vcf_type]=0
            all_filter_sum[vcf_type]=0
            for tag in filter_tags:
                all_filter_counts[vcf_type][tag]=0
        print("var_type\tChr\tVariants\tfrac_filtered\t%s\t%s" % ("\t".join(filter_tags), "\t".join("frac_"+x for x in filter_tags) ), file=log )
        for vcf_type in ["SNV", "CNV"]:
            for c in chroms:
                chrom, recF, knownRecF, records_gene, comphetF, chr_ignored_SO_terms, var_count, filter_counts = data[c]
                fsum = sum(filter_counts[vcf_type].values())
                all_var_count[vcf_type] += var_count[vcf_type]
                all_filter_sum[vcf_type] += fsum
                for tag in filter_tags:
                    all_filter_counts[vcf_type][tag]+=filter_counts[vcf_type][tag]
                print("%s\t%s\t%i\t%f\t%s\t%s" % (vcf_type, 
                                                  chrom, 
                                                  var_count[vcf_type], 
                                                  float(fsum)/var_count[vcf_type] if var_count[vcf_type]>0 else 0, 
                                                  "\t".join(str(filter_counts[vcf_type][x]) for x in filter_tags), 
                                                  "\t".join(str(filter_counts[vcf_type][x]/float(fsum) if fsum > 0 else 0) for x in filter_tags) ), file=log)
            print("%s\t%s\t%i\t%f\t%s\t%s" % (vcf_type, "ALL",
                                      all_var_count[vcf_type], 
                                      float(all_filter_sum[vcf_type])/all_var_count[vcf_type] if all_var_count[vcf_type]>0 else 0, 
                                      "\t".join(str(all_filter_counts[vcf_type][x]) for x in filter_tags), 
                                      "\t".join(str(all_filter_counts[vcf_type][x]/float(all_filter_sum[vcf_type]) if all_filter_sum[vcf_type] > 0 else 0) for x in filter_tags) ), file=log)
        total_var_count=0
        total_filter_sum=0
        total_filter_counts=OrderedDict()
        for tag in filter_tags:
            total_filter_counts[tag]=0
        vcf_types = ["SNV", "CNV"] if cnvF_chr is not None else ["SNV"]
        for vcf_type in vcf_types:
            total_var_count+=all_var_count[vcf_type]
            total_filter_sum+=all_filter_sum[vcf_type]
            for tag in filter_tags:
                total_filter_counts[tag]+=all_filter_counts[vcf_type][tag]
        print("%s\t%s\t%i\t%f\t%s\t%s" % ("ALL", "ALL",
                                  total_var_count, 
                                  float(total_filter_sum)/total_var_count if total_var_count>0 else 0, 
                                  "\t".join(str(total_filter_counts[x]) for x in filter_tags), 
                                  "\t".join(str(total_filter_counts[x]/float(total_filter_sum) if total_filter_sum > 0 else 0) for x in filter_tags) ), file=log)
    bgzip(outF, index=True, delinFile=True, maxthreads=threads)
    
    # write merged debug log
    dfiles = [outdir+"."+c+".log_tmp" for c in chroms]
    if debug:
        with open(outdir + config["dataset_name"]+".v2r.debug.log", 'w') as outfile:
            for fname in dfiles:
                with open(fname) as infile:
                    for line in infile:
                        outfile.write(line)
    remove_file(dfiles)
    
    # =================================================================================
    # clean up
    # =================================================================================
    remove_file(snvF_chr_f)
    remove_file([x+".tbi" for x in snvF_chr_f if x is not None])
    remove_file(cnvF_chr_f)
    remove_file([x+".tbi" for x in cnvF_chr_f if x is not None])
    
    logging.warn("\nNOTE: the following SO terms were ignored: %s" % (ignored_SO_terms))
    logging.info("Finished in " + str(datetime.timedelta(seconds=time.time()-startTime)))
    print("Done.")

alias2gid, terms, exon_int, utr_int, gene_int, id2reg_regions, reg_regions = None,None,None,None,None,None,None