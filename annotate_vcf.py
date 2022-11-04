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
    Annotates a set of VCFs file using snpEff and vcfanno.
    Input:
        - sample sheet with 4 columns:
            - id: used as output filename prefix
            - genome: string that must match the respective entry in the VCF config file (see below)
            - source: 'platypus'|'deepvariant'|'gatk'. Used for VCF preprocessing.
            - vcf_file: path to VCF file
        - config file:
            - JSON file with the following structure:
                - prefix: in this section, variables can be defined that will be replacesd throughoput the remaining config file.
                    Example: "prefix" : { "@REF": "/path_prefix" } 
                - ref: reference genome definitions. The contained section names must match the entries in the sample sheet.
                    - subsections: FASTA (reference genomme fasta file), snpEff (snpEff annotation db version), chr (list of considered chromosome names, filter (snpEff filter string)
                - roi: regions of interest annotation BED files. Must have subsections for each configured reference genome.
                - af: population allele frequency VCF files.  Must have subsections for each configured reference genome and
                    beyond this 2 subsections "global" (global AF annotations) and "pop" (population specific AF annotations)
                    Each entry has 2 parts: VCF file and comma-separated list of info fields to be loaded from these vcf files
                - anno: general annotation files. Must have subsections for each configured reference genome.
                    Each entry has multiple parts, depending on file type:
                    BED files: [file_path]
                    VCF files: [file_path, 
                                comma-separated_list_of_info_field_names,
                                comma-separated_list_of_summarization_methods # see vcfanno docs; if omitted, 'self' will be used]                                
                    TSV files: [file_path, 
                                comma-separated_list_of_column_indices, 
                                comma-separated_list_of_variable_names, # if omitted, the TSV column names will be used 
                                comma-separated_list_of_summarization_methods # see vcfanno docs; if omitted, 'self' will be used]
                - known: VCF files with 'known' variants.  Must have subsections for each configured reference genome.
                    The respective IDs will be written to the annotated VCF iD section. It is possible to
                    configure a prefix string for these IDs
                - output: List of output fields for the created TSV file.
                    There must be a 'fields' section containing list entries with a column name and an optional operator: 
                        'max': split by comma, replace '.' values with nan (will be ignored) and select maximum value
                        'min': split by comma, replace '.' values with nan (will be ignored) and select minimum value
                        none: use value as is.
                - tools: Optional section for providing custom paths for 3rd party tools:
                    - vcf-sort
                    - snpSift
                    - snpEff
                    - vcfanno
                    If omitted, the pipeline will try to call the tool directly by name.
                - linux_temp_dir: optional, for configuring an alternative TEMP dir. 
                    
                    
"""

import csv, datetime, time, logging, sys, os, json
import re
from pandas.core.frame import DataFrame
import vcfpy
import numpy as np
import pandas as pd
from utils import pipelineStep, files_exist, remove_file, bgzip, move_id_to_info_field, add_contig_headers, get_tool, load_config

if __name__ == '__main__':    
    pass

#===================== utility ===========================================
class myset(set):
    """ overloads set __str__ method for nicer serialization """
    def __str__(self):
        return ",".join("{0}".format(n) for n in self)
    def __repr__(self):
        return self.__str__()  

def check_success(stage, success, start_time):
    """ check whether a particular processing stage was successful """
    elapsed_time = time.time() - start_time
    if not success:
        logging.error("Pipeline failed at stage: " + stage)
        sys.exit("Pipeline failed at stage: " + stage)
    else:
        logging.info("-----------------------------------------------------------------------------")
        logging.info("Finished stage " + stage + " in " + time.strftime("%H:%M:%S", time.gmtime(elapsed_time)))
        logging.info("-----------------------------------------------------------------------------")


def preprocessPlatypus(config, inVcf, ref, outVcf, tmpdir=None):
    """ Preprocess platypus VCF files. Sorts file and adds contig info headers """
    start_time = time.time()
    success = True
    logging.info("Adding CONTIG info to platypus VCF")
    vcfsorted = outVcf + ".sorted.vcf"    
    cmd = [get_tool(config, "vcf-sort")]
    if tmpdir:
        cmd += ["--temporary-directory", tmpdir]
    cmd+=[inVcf]
    pipelineStep(inVcf, vcfsorted, cmd, shell=True, stdout=vcfsorted)
    # fix CONTIG headers
    vcfsortedc = vcfsorted + ".CONTIG.vcf"
    add_contig_headers(vcfsorted, ref, outVcf)
    check_success("preprocessPlatypus", True, start_time)  
    return outVcf

def preprocessDeepvariant(inVcf):
    """ Preprocess DeepVariant files. Will remove ID entries. """
    start_time = time.time()
    success = True
    logging.info("Removing deepvariant entries from ID field")
    vcf_no_id = inVcf + ".no_id.vcf"
    move_id_to_info_field(inVcf, "old_deepvariant_ids", vcf_no_id)
    os.rename(vcf_no_id, inVcf)
    check_success("preprocessDeepvariant", True, start_time)  
    return inVcf

def preprocessGatk(inVcf):
    """ Preprocess GATK VCF files. Will remove ID entries.  """
    start_time = time.time()
    logging.info("Removing gatk entries from ID field")
    vcf_no_id = inVcf + ".no_id.vcf"
    move_id_to_info_field(inVcf, "old_gatk_ids", vcf_no_id)
    os.rename(vcf_no_id, inVcf)
    check_success("preprocessGatk", True, start_time)  
    return inVcf

def run_vcfanno(inVCF, outVCF, configF, luaF=None, threads=1, overwrite=True, additionalParameters=[], settings=[], exe='vcfanno'):
    """ run vcfanno """
    start_time = time.time()
    success = True
    if files_exist(outVCF) and not overwrite:
        print("VCF file " + inVCF + " already exists! SKIPPING re-creation...")
        logging.warn("VCF file " + inVCF + " already exists! SKIPPING re-creation...")
        return success
    iscompressed = outVCF.endswith('.gz')
    vcfraw = outVCF + ".raw.vcf"
    cmd = settings + [ exe ]
    if luaF is not None:
        cmd += ["-lua", luaF ]
    cmd += [ "-p", str(threads) ]
    cmd += [ configF ]
    cmd += [ inVCF ] + additionalParameters
    success = success and pipelineStep(inVCF, vcfraw, cmd, shell=True, stdout=vcfraw)
    check_success("vcfanno", success, start_time)  
    if iscompressed:
        bgzip(vcfraw, outFile=outVCF, index=True, delinFile=True)
    else:
        os.rename(vcfraw, outVCF)
    return success

def selectMax(df, colname, min2maxValues):
    """ selects the maximum category value """
    r=list(range(0,len(min2maxValues))) # [0, 1, ..., n]
    # split by comma, replace '.' value, map to number, select max, map back to name.
    return df[colname].str.split(",", expand=True).replace(".", np.nan).replace(min2maxValues, r).max(axis=1).replace(r, min2maxValues)


def postfilter(config, infile, outfile, genome, config_file, samples, overwrite=False):
    """ Postfilter annotated VCF and create TSV output file """
    start_time = time.time()
    success = True
    if files_exist(outfile) and not overwrite:
        logging.error("Outputfile exists!")
        return False
    # defines transitions and transversions
    TsTvMatrix = {"AG":"Ts", "GA": "Ts", 
              "CT": "Ts", "TC": "Ts",
              "AC": "Tv", "CA": "Tv",
              "AT": "Tv", "TA": "Tv",
              "CG": "Tv", "GC": "Tv",
              "GT": "Tv", "TG": "Tv"}  
    # to avoid "_csv.Error: field larger than field limit (131072)" errors
    # @see https://stackoverflow.com/questions/15063936/csv-error-field-larger-than-field-limit-131072
    maxInt = sys.maxsize
    while True:
        # decrease the maxInt value by factor 10 
        # as long as the OverflowError occurs.
        try:
            csv.field_size_limit(maxInt)
            break
        except OverflowError:
            maxInt = int(maxInt/10)           
    # expected input dataypes (autodetect others)
    dtypes={"CHROM":"object", "POS": "int64", "REF": "object", "ALT": "object", "ID": "object", "FILTER": "object", 
            "ANN[*].GENE": "object", "ANN[*].EFFECT": "object", "ANN[*].IMPACT": "object"}
    if "roi" in config and genome in config["roi"]:
        for roi in config["roi"][genome].keys():
            dtypes[roi]="bool"
    # used for column renaming   
    sampledict_gt={}
    sampledict_gq={}
    for i,v in enumerate(samples):
        sampledict_gt[i]="GT_" + v
        sampledict_gq[i]="GQ_" + v
        dtypes[sampledict_gt[i]]="object"
        dtypes[sampledict_gq[i]]="int64"
    writeHeader=True
    mode='w'
    # get list of AF fields
    AFS={}
    if "af" in config and genome in config["af"]:
        for k in config["af"][genome].keys():
            fields=[]
            for a, af in config["af"][genome][k].items():  
                if af[1].count(",") == 1:  # configured AC,AN and calculated AF
                    fields.append(k+"_"+a+"_AF")
                else:
                    fields.append(k+"_"+a+"_"+af[1])
            AFS[k]=fields     

    chunksize=1000000
    for df in pd.read_csv(infile,delimiter='\t',encoding='utf-8', comment="#", float_precision='high', chunksize=chunksize, na_values=["na","NA","."], dtype=dtypes): 
        # create output data frame
        o=DataFrame()
        o["CHROM"]=df["CHROM"]
        o["POS"]=df["POS"]
        o["REF"]=df["REF"]
        o["ALT"]=df["ALT"]
        o["FILTER"]=df["FILTER"]
        o["TYPE"] = (df["REF"].str.len()+df["ALT"].str.len()==2).map({ True: "SNV", False: "INDEL"})
        o["IsKnown"] = df["ID"].isna().replace({ True: "0", False: "1"})
        o["ID"]=df["ID"]
        try:
            o["Genes"]=np.array(list(map(myset, df["ANN[*].GENE"].astype('str').str.split(",", expand=False).values))) # get unique values
        except Exception as e: 
            logging.error(getattr(e, 'message', repr(e)))
            pos0=df.iloc[0]["CHROM"] + ":" + str(df.iloc[0]["POS"])
            pos1=df.iloc[len(df.index)-1]["CHROM"] + ":" + str(df.iloc[len(df.index)-1]["POS"])
            logging.error("Error converting Genes between %s and %s" % (pos0, pos1))
            o["Genes"]=None           
        try:
            o["Effects"]=list(map(myset, df["ANN[*].EFFECT"].astype('str').str.split(",", expand=False).values))
        except Exception as e:
            logging.error(getattr(e, 'message', repr(e)))
            pos0=df.iloc[0]["CHROM"] + ":" + str(df.iloc[0]["POS"])
            pos1=df.iloc[len(df.index)-1]["CHROM"] + ":" + str(df.iloc[len(df.index)-1]["POS"])
            logging.error("Error converting Effects between %s and %s" % (pos0, pos1))
            o["Effects"]=None        
        try:
            o["maxImpact"]=selectMax(df, "ANN[*].IMPACT", ["UNKNOWN","MODIFIER","LOW","MODERATE","HIGH"])
        except Exception as e:
            logging.error(getattr(e, 'message', repr(e)))
            pos0=df.iloc[0]["CHROM"] + ":" + str(df.iloc[0]["POS"])
            pos1=df.iloc[len(df.index)-1]["CHROM"] + ":" + str(df.iloc[len(df.index)-1]["POS"])
            logging.error("Error converting maxImpact between %s and %s" % (pos0, pos1))
            o["maxImpact"]=None     
        # create clean AFS columns
        if "af" in config and genome in config["af"]:
            for k in config["af"][genome].keys():
                clean=[]
                for a in AFS[k]:
                    if a in df.columns and not df[a].isnull().all(): # ignore Nan columns
                        if df[a].dtype==np.float64:
                            df[a+"_clean"]=df[a]
                        else:  # split comma-separated values and get max value
                            df[a+"_clean"]=df[a].astype(str).replace('nan',np.nan).str.split(",", expand=True).replace(["."], np.nan).fillna(value=np.nan).apply(pd.to_numeric).max(axis=1)
                        clean.append(a+"_clean")
                        if not "maxPopAF_"+k in o.columns:
                            o["maxPopAF_"+k]=(a+"="+df[a].astype(str))
                        else:
                            o["maxPopAF_"+k]=o["maxPopAF_"+k]+","+(a+"="+df[a].astype(str))
                # map to categories (or NaN if no value!)
                o["maxPopAF_cat_"+k]=pd.cut( df[clean].max(axis=1), bins=[ 0, 0.0001, 0.001, 0.02, 1], right=True, labels = [ "VERY_RARE", "RARE", "COMMON", "VERY_COMMON"])
        o["TsTv"]=(df["REF"].str.upper()+df["ALT"].str.upper()).map(TsTvMatrix)        
        # add configured output fields
        for f in config["output"]["fields"]:
            if f[0] in df:
                if f[1]=="max":
                    # split by comma, replace '.' values with nan (will be ignored) and select maximum value
                    # TODO: support other operations (min, mean, etc.)
                    o[f[0]]=df[f[0]].astype(str).str.split(",", expand=True).replace([".", "na", "NA"], [np.nan,np.nan,np.nan]).astype(float).max(axis=1)
                elif f[1]=="min":
                    # split by comma, replace '.' values with nan (will be ignored) and select minimum value
                    # TODO: support other operations (min, mean, etc.)
                    o[f[0]]=df[f[0]].astype(str).str.split(",", expand=True).replace([".", "na", "NA"], [np.nan,np.nan,np.nan]).astype(float).min(axis=1)
                else:
                    o[f[0]]=df[f[0]]
        # add ROIs
        if ("roi" in config and genome in config["roi"]):
            for roi in config["roi"][genome].keys():
                if roi in df:
                    o["ROI_" + roi]=df[roi].map({True: 1, False: 0})        
        # add GT and GQ fields
        gt=pd.DataFrame( df['GEN[*].GT'].astype(str).str.split(",", expand=True) ).rename(columns=sampledict_gt)
        gq=pd.DataFrame( df['GEN[*].GQ'].astype(str).str.split(",", expand=True) ).rename(columns=sampledict_gq)
        o=pd.concat([o,gt,gq], ignore_index=False, axis=1)
        # write output
        o.to_csv(outfile, sep='\t', mode=mode, header=writeHeader, index=False, na_rep=".")
        # from now on: append to file
        mode='a'
        writeHeader=False
    check_success("postfilter", success, start_time)
    return True
                
#===========================================================================

def annotate_vcf(config_file, infile, outdir, maxthreads=1, snpEffConfig=None, includeNonPass=False, noQC=False, overwrite=False):
    """ Annotates a VCF file with snpEff and vcfanno """    
    start_time=time.time()
    success=True
    # load config and merge with default tool config to enable system-specific program locations
    config = load_config(config_file)
    # confgured linux_temp_dir?
    linux_temp_dir = config["linux_temp_dir"] if "linux_temp_dir" in config else None
    # ensure dirs
    if not os.path.exists(outdir):
            logging.info("Creating dir " + outdir)
            os.makedirs(outdir)
    if not outdir.endswith("/"):
        outdir += "/"
    log_file = os.path.join(outdir, "annotate_vcf.log")   
    # start 
    print("Logging output to %s" % log_file)
    for handler in logging.root.handlers[:]:
        logging.root.removeHandler(handler)
    logging.basicConfig(filename=log_file, level=logging.DEBUG)                
    logging.info("==========================================================")
    logging.info("annotate_vcf")
    logging.info("==========================================================")
    logging.info("Effective configuration:")
    logging.info(json.dumps(config, indent=4, sort_keys=True))
    logging.info("Started script at %s.", str(datetime.date.today()))   
    with open(infile, 'rt') as tsvin:
        tsvin = csv.reader(tsvin, delimiter='\t')
        for row in tsvin:
            if (row[0].startswith("#")):
                continue
            PID = row[0]
            genome = row[1]
            source = row[2]
            inVcf = row[3]      
            sample_outdir = outdir + PID 
            if not files_exist(inVcf):
                logging.warn("Skipping configured input vcf as file not found: " + inVcf)
                continue
            if not files_exist(sample_outdir):
                os.makedirs(sample_outdir)
            tmpdir = sample_outdir + "/tmp"
            if not files_exist(tmpdir):
                os.makedirs(tmpdir)
            qcdir = sample_outdir + "/qc"
            if not files_exist(qcdir):
                os.makedirs(qcdir)
            ds = PID + "." + genome
            vcfa = tmpdir + "/" + ds + ".anno.vcf.gz"  # temporary unfiltered VCF
            vcfaf = sample_outdir + "/" + ds + ".anno+fil.vcf.gz"
            tsv = tmpdir + "/" + ds + ".anno+fil.tsv"  # temporary snpSift output
            final = sample_outdir + "/" + ds + ".final.tsv"
            logging.info("-------------------------")
            logging.info(ds)
            logging.info("-------------------------")
            if overwrite or not files_exist([vcfa]):
                tmp = tmpdir + "/" + ds + ".tmp1.vcf"
                tmp2 = tmpdir + "/" + ds + ".tmp2.vcf"
                # ###############################################################################
                #             NORMALIZE VCF (e.g., remove non-canonical chroms, non pass variants, etc.)
                # ###############################################################################
                norm_start_time = time.time()
                f = []
                for c in config["ref"][genome]["chr"]:    
                    f += ["( CHROM = '" + c + "' )"]
                filter_string = "|".join(f)
                if not includeNonPass:
                    filter_string = "(" + filter_string + ") & (( na FILTER ) | (FILTER = 'PASS'))"
                cmd = [get_tool(config, "snpSift"), "filter", "\"" + filter_string + "\"", "-f", inVcf ]
                success = success and pipelineStep(inVcf, tmp, cmd, shell=True, stdout=tmp)
                check_success("Normalization", success, norm_start_time)
                # ###############################################################################
                #             PREPROCESS VCF 
                # ###############################################################################
                if source.lower() == "platypus":
                    # ###############################################################################
                    # fix CONTIG headers               
                    # ###############################################################################
                    preprocessPlatypus(config, inVcf=tmp, ref=config["ref"][genome]["FASTA"], outVcf=tmp, tmpdir=linux_temp_dir)
                elif source.lower() == "deepvariant":
                    # ###############################################################################
                    # remove IDs as they are not informative        
                    # ###############################################################################
                    preprocessDeepvariant(inVcf=tmp)
                elif source.lower() == "gatk":
                    # ###############################################################################
                    # remove IDs as they are not informative        
                    # ###############################################################################
                    preprocessGatk(inVcf=tmp)                
                # ###############################################################################
                #             Annotate with snpEFF
                # ###############################################################################
                if (genome in config["ref"]):
                    snpEff_start_time=time.time()
                    cmd = [get_tool(config, "snpEff"), "ann", config["ref"][genome]["snpEff"]]
                    cmd += ["-i", "vcf"]
                    #cmd += ["-t"] # for some reason this does not work with Singularity?
                    cmd += ["-csvStats", tmpdir + "/" + ds + ".snpEff.ann.csvStats.txt"] 
                    if (snpEffConfig is not None):
                        cmd += ["-config", snpEffConfig]
                    cmd += ["-noLog"]
                    cmd += [tmp]
                    success = success and pipelineStep(tmp, tmp2, cmd, shell=True, stdout=tmp2)
                    check_success("snpEff", success, snpEff_start_time) 
                    os.rename(tmp2, tmp)
                else:
                    logging.warn("Skipped snpEff annotations as not defined for genome " + genome)
                # ###############################################################################
                #             Annotate with vcfanno
                # ###############################################################################  
                vcfaconf = tmpdir + "/" + ds + ".vcfanno.config"
                vcfalua = tmpdir + "/" + ds + ".vcfanno.lua"
                # create lua file
                with open(vcfalua, "w") as f: 
                    f.write("function setid(pre,...)\n")
                    f.write(" local t = {...}\n")
                    f.write(" local res = {}\n")
                    f.write(" local seen = {}\n")
                    f.write(" for idx, ids in pairs(t) do\n")
                    f.write("  local sep=\",\"\n")  # split comma-separated IDs
                    f.write("  for v in string.gmatch(ids, \"([^\"..sep..\"]+)\") do\n")
                    f.write("   for i, v in pairs(t) do\n")
                    f.write("    if v ~= \".\" and v ~= nil and v ~= \"\" then\n")
                    f.write("     if seen[v] == nil then\n")
                    f.write("      res[#res+1] = string.gsub(pre[i] .. v, \",\", \";\")\n")
                    f.write("      seen[v] = true\n")
                    f.write("     end\n")
                    f.write("    end\n")
                    f.write("   end\n")
                    f.write("  end\n")
                    f.write(" end\n")
                    f.write(" return table.concat(res, \";\")\n")
                    f.write("end\n")
                logging.info("Created LUA file at " + vcfalua)                    
                # create vcfanno config file     
                with open(vcfaconf, "w") as f: 
                    if ("roi" in config and genome in config["roi"]):
                        for r, rf in config["roi"][genome].items():
                            f.write("[[annotation]]\n")
                            f.write("file=\"" + rf + "\"\n")
                            f.write("columns = [3]\n")  # use column 3 as some BED do not contain 4 columns
                            f.write("ops = [\"flag\"]\n")
                            f.write("names = [\"" + r + "\"]\n")
                            f.write("\n")
                    else:
                        logging.warn("Skipped ROI annotation as not defined for genome " + genome)                       
                    if ("af" in config and genome in config["af"]):
                        file2block = {}
                        postanno = {}
                        for k in config["af"][genome].keys():
                            for a, af in config["af"][genome][k].items():
                                block = file2block[af[0]] if af[0] in file2block else [[], [], []]  # fields, ops, names
                                if af[1].count(",") == 1:  # configured AC,AN. We need a postannoblock
                                    postanno[k + "_" + a]=[]
                                    for x in af[1].split(","):
                                        postanno[k + "_" + a].append( "\"" + k + "_" + a + "_" + x +"\"" )
                                for x in af[1].split(","):
                                    block[0].append("\"{0}\"".format(x))
                                if (len(af) > 2):
                                    for x in af[1].split(","):
                                        block[1].append("\"{0}\"".format(x))
                                else:
                                    for x in af[1].split(","):
                                        block[1].append("\"max\"")
                                for x in af[1].split(","):
                                    block[2].append("\"" + k + "_" + a + "_{0}\"".format(x))
                                file2block[af[0]] = block
                        # write blocks
                        for file in file2block:    
                            block = file2block[file]
                            f.write("[[annotation]]\n")
                            f.write("file=\"" + file + "\"\n")
                            # select AF INFO fields
                            f.write("fields = [" + ", ".join(block[0]) + "]\n")
                            f.write("ops = [" + ", ".join(block[1]) + "]\n")
                            f.write("names = [" + ", ".join(block[2]) + "]\n")
                            f.write("\n")
                        # write postanno blocks
                        for k in postanno:
                            # divide ac/an
                            f.write("[[postannotation]]\n")
                            f.write("fields=[" + ", ".join(postanno[k]) + "]\n")
                            f.write("name=\"" + k + "_AF\"\n")
                            f.write("op=\"div2\"\n")
                            f.write("type=\"Float\"\n")
                            f.write("\n")
                    else:
                        logging.warn("Skipped AF annotation as not defined for genome " + genome)
                    if ("anno" in config and genome in config["anno"]):
                        for a, af in config["anno"][genome].items():
                            if af[0].endswith("bed") or af[0].endswith("bed.gz"):
                                # BED
                                f.write("[[annotation]]\n")
                                f.write("file=\"" + af[0] + "\"\n")
                                f.write("columns = [4]\n")  # copy value from 4th column (=BED name)
                                if (len(af) > 1):
                                    f.write("ops = [" + ", ".join("\"{0}\"".format(e) for e in af[1].split(",")) + "]\n")
                                else:
                                    f.write("ops = [\"self\"]\n")
                                f.write("names = [\"" + a + "\"]\n")
                                f.write("\n")
                            elif af[0].endswith("tsv") or af[0].endswith("tsv.gz") or af[0].endswith("txt") or af[0].endswith("txt.gz"):
                                # TSV: 0=file, 1=columns, 2=names, 3=ops
                                f.write("[[annotation]]\n")
                                f.write("file=\"" + af[0] + "\"\n")
                                # columns indices
                                f.write("columns = [" + ", ".join("{0}".format(e) for e in af[1].split(",")) + "]\n")
                                # column names
                                if (len(af) > 2):
                                    f.write("names = [" + ", ".join("\"" + a + "_{0}\"".format(e) for e in af[2].split(",")) + "]\n")
                                else: 
                                    # use column indices as names...
                                    f.write("names = [" + ", ".join("\"" + a + "_{0}\"".format(e) for e in af[1].split(",")) + "]\n")
                                # ops
                                if (len(af) > 3):
                                    f.write("ops = [" + ", ".join("\"{0}\"".format(e) for e in af[3].split(",")) + "]\n")
                                else:
                                    f.write("ops = [" + ", ".join("\"self\"".format(e) for e in af[1].split(",")) + "]\n")
                                f.write("\n")
                            else:
                                # VCF
                                f.write("[[annotation]]\n")
                                f.write("file=\"" + af[0] + "\"\n")
                                # select INFO fields
                                f.write("fields = [" + ", ".join("\"{0}\"".format(e) for e in af[1].split(",")) + "]\n")
                                if (len(af) > 2):
                                    f.write("ops = [" + ", ".join("\"{0}\"".format(e) for e in af[2].split(",")) + "]\n")
                                else:
                                    f.write("ops = [" + ", ".join("\"self\"".format(e) for e in af[1].split(",")) + "]\n")
                                f.write("names = [" + ", ".join("\"" + a + "_{0}\"".format(e) for e in af[1].split(",")) + "]\n")
                                f.write("\n")
                    else:
                        logging.warn("Skipped ANNO annotation as not defined for genome " + genome)                
                    if ("known" in config and genome in config["known"]):
                        knownFields=[]
                        knownFieldsQuote=[]
                        prefixes=[]
                        for a, af in config["known"][genome].items():
                            # annotate =INFO field
                            f.write("[[annotation]]\n")
                            f.write("file=\"" + af[0] + "\"\n")
                            f.write("fields=[\"ID\"]\n")
                            f.write("names=[\"" + a + "_ID\"]\n")
                            f.write("ops=[\"self\"]\n")
                            knownFieldsQuote+=["\"" + a + "_ID\""]
                            knownFields+=[a + "_ID"]
                            prefix = ""
                            if (len(af) > 1):
                                prefix = af[1]
                            prefixes+=["'" + prefix +"'"]
                        # postannotate: move from INFO fields to ID
                        f.write("[[postannotation]]\n")
                        f.write("name=\"ID\"\n")
                        f.write("fields=[" + ",".join(knownFieldsQuote)+"]\n")
                        #f.write("op=\"lua:setid(" + prefixes + "', " + a + "_ID)\"\n")
                        f.write("op=\"lua:setid({"+",".join(prefixes)+"}, "+",".join(knownFields)+")\"\n")
                        f.write("type=\"String\"\n")
                        # postannotate: delete INFO fields 
                        f.write("[[postannotation]]\n")
                        f.write("fields=[" + ",".join(knownFieldsQuote)+"]\n")
                        f.write("op=\"delete\"\n")
                        f.write("\n")
                logging.info("Created vcfanno configuration file at " + vcfaconf)   
                # run vcfanno
                run_vcfanno(inVCF=tmp, outVCF=tmp2, configF=vcfaconf, luaF=vcfalua, threads=maxthreads, 
                            settings=["GOGC=1000", "IRELATE_MAX_CHUNK=8000", "IRELATE_MAX_GAP=1000"],
                            exe=get_tool(config, "vcfanno"))
                os.rename(tmp2, tmp)
                # ###############################################################################
                #             compress results
                # ###############################################################################
                if vcfa.endswith(".gz"):
                    bgzip(tmp, vcfa, index=True, override=True, delinFile=True, maxthreads=maxthreads)
                else:
                    os.rename(tmp, vcfa)
            else:
                logging.warn("Annotated VCF exists, will not re-create: %s " % vcfa)    
            if overwrite or not files_exist([vcfaf]):
                # ###############################################################################
                #             Filter VCF (e.g., remove non-canonical chroms, non pass variants, etc.
                # ###############################################################################
                snpSift_start_time=time.time() 
                snpEffFilterF = tmpdir + "/" + ds + ".snpEffFilter.txt"
                with open(snpEffFilterF, "w") as f: 
                    f.write(config["ref"][genome]["filter"])
                tmp = vcfaf + ".tmp.vcf"
                cmd = [get_tool(config, "snpSift"), "filter"]
                cmd += ["-f", vcfa ]
                cmd += ["-e", snpEffFilterF]
                cmd += ["-i", "V2filter"]
                cmd += ["-p"]
                success = success and pipelineStep([vcfa, snpEffFilterF], tmp, cmd, shell=True, stdout=tmp)
                if vcfaf.endswith(".gz"):
                    bgzip(tmp, vcfaf, index=True, override=True, delinFile=True, maxthreads=maxthreads)
                else:
                    os.rename(tmp, vcfaf)
                check_success("Filtering", success, snpSift_start_time)        
            reader = vcfpy.Reader.from_path(vcfaf)
            samples = reader.header.samples.names
            logging.info("Samples: %s" % samples)
            if overwrite or not files_exist([tsv + ".gz"]):
                # ###############################################################################
                #             Extract fields to TSV
                # ###############################################################################
                extract_start_time=time.time()
                # "ANN[*].RANK", "ANN[*].CDNA_POS", "ANN[*].CDNA_LEN", "ANN[*].CDS_POS", "ANN[*].CDS_LEN","ANN[*].AA_POS", "ANN[*].AA_LEN", "ANN[*].DISTANCE",
                headers = []
                for h in ["CHROM", "POS", "ID", "REF", "ALT", "FILTER", "AF", "AC", "DP", "MQ",
                          "ANN[*].ALLELE", "ANN[*].EFFECT", "ANN[*].IMPACT", "ANN[*].GENE", "ANN[*].GENEID",
                          "ANN[*].FEATURE", "ANN[*].FEATUREID", "ANN[*].BIOTYPE", "ANN[*].HGVS_C",
                          "ANN[*].HGVS_P", "ANN[*].ERRORS", "GEN[*].GT", "GEN[*].GQ"]: 
                    headers += ["\"" + h + "\""]
                # add INFO fields from VCF
                pattern = re.compile("^[A-Za-z_][0-9A-Za-z_.]*\Z")
                for h in reader.header.info_ids():
                    # check compatibility
                    if pattern.match(h):
                        headers += ["\"" + h + "\""]
                    else:
                        logging.warn("Cannot extract INFO field " + h + " as not VCF compliant")
                cmd = [get_tool(config, "snpSift"), "extractFields"]
                cmd += ["-s", "\",\""]
                cmd += ["-e", "\".\""]
                cmd += [vcfaf]
                cmd += headers
                success = success and pipelineStep(vcfaf, tsv, cmd, shell=True, stdout=tsv)
                bgzip(tsv, delinFile=True, override=True, maxthreads=maxthreads)
                check_success("Extracting", success, extract_start_time)
            #
            # postprocess
            #
            if overwrite or not files_exist([final + ".gz"]):   
                postfilter(config, tsv + ".gz", final, genome, config_file, samples, overwrite=True)
                bgzip(final, delinFile=True, override=True, maxthreads=maxthreads)
    # check success
    check_success("pipeline", success, start_time)
    logging.info("Finished script successfully at %s.", str(datetime.date.today()))
    print("All done.")    
