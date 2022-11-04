#
#
# Copyright (C) 2020 XXX.  All rights reserved.
#
# This file is part of cohort_varan
# 
# See the file LICENSE for redistribution information.
#
# @author: niko.popitsch
if __name__ == '__main__':
    pass


from argparse import ArgumentTypeError
import sys, os, logging, gzip, json
from subprocess import check_output, STDOUT, CalledProcessError
import pysam


def get_tool(config, key):
    if (key not in config["tools"]):
        return key
    return config["tools"][key]

def load_config(confF):
    """ Load JSON config file and replace variables from 'prefix' block """
    config=json.load(open(confF))
    if "prefix" in config:
        s=json.dumps(config, separators=(',', ':'))
        for (k, v) in config["prefix"].items():
            s=s.replace(k, v)  
        config=json.loads(s)
    return config

def runTask(cmd, shell=False):
    """ run an external task """
    logging.info(cmd)
    if shell:
        out = check_output(" ".join(cmd), shell=True, stderr=STDOUT)
    else: 
        out = check_output(cmd, stderr=STDOUT)
    return out
        
def files_exist(files):
    """ check whether (a list of) files exist """
    if (type(files) is list) :
        for f in files:
            if f is None:
                return False
            if not os.path.exists(f):
                return False
    else:
        if files is None:
            return False
        if not os.path.exists(files):
            return False
    return True

def existing_file(files):
    """ use in argument parser """
    if files_exist(files):
        return files
    if (type(files) is list) :
        raise ArgumentTypeError("Not all files exist ["+",".join(files)+"]")
    else:
        raise ArgumentTypeError("Not all files exist ["+(files)+"]")
    
def remove_file(files):
    """ remove a (list of) file(s) (if it/they exists)"""
    if (type(files) is list) :
        for f in files:
            if f is not None and os.path.exists(f):
                os.remove(f)
    else:
        if files is not None and os.path.exists(files):
            os.remove(files)

def check_file(files):
    if (type(files) is list) :
        for f in files:
            check_file(f)
    else :
        if not files is None and not os.path.exists(files):
            print("Error: file", files, "was not found! Exiting...")
            sys.exit(1)       
                
def pipelineStep(inputfile, outFile, cmd, shell=False, stdout=None, append=False, logfile=None):
    try:
        if inputfile is not None:
            if (type(inputfile) is list):
                for f in inputfile:
                    check_file(f)
            else:
                check_file(inputfile) 
        
        if stdout is not None:
            if shell is False:
                raise ArgumentTypeError("When using the parameter stdout, the shell parameter must be set to True!")
            else: 
                if append:
                    cmd.append(">>")
                else:
                    cmd.append(">")
                cmd.append(stdout)                                                                                                                                                                                                                                                                 

                
        # Overwrite output file?!
        if outFile is not None:
            if (type(outFile) is list) :
                for f in outFile:
                    if files_exist(f):
                        logging.warn("Overwriting file %s", f)
                else:
                    if files_exist(outFile):
                        logging.warn("Overwriting file %s", outFile)

        out = runTask(cmd, shell)         
        if logfile is None:
            logging.info(out)
        else:
            with open(logfile, "a") as log:
                log.write(out)
        
        if outFile is not None:
            check_file(outFile)
        return True
    except CalledProcessError as e:
        logging.error(e.output)
        logging.error("ERROR %s - removing outputfile %s", e, outFile)
        if outFile is not None:
            if (type(outFile) is list) :
                for f in outFile:
                    remove_file([f]) 
            else:
                remove_file([outFile]) 
        return False     
    
    
# =================================================================
# specific tools
# =================================================================
def bgzip(inFile, outFile=None, index=False, override=False, delinFile=False, maxthreads=1, exe='bgzip'):
    """ bgzip a file """
    if outFile == None:
        outFile = inFile+".gz"
    success = True    
    if(not files_exist(outFile) or override):                 
        success = success and pipelineStep(inFile, outFile, [exe, "-@", str(maxthreads), "-c", inFile], shell=True, stdout=outFile)
        if(index):
            if outFile.endswith(".vcf"):
                success = success and index_vcf(outFile)
            else:
                success = success and tabix(outFile)
        if(success and delinFile):
            remove_file(inFile)
    else: 
        logging.info("Skipping bgzip. " + outFile + " already exists.");
    return success

def index_vcf(inFileVcf, override=False, exe='tabix'):
    """ index a VCF file """
    success = True
    idxFile = inFileVcf + ".tbi"
    if(not files_exist(idxFile) or override):
        success = success and pipelineStep(inFileVcf, idxFile, [exe, "-p", "vcf", inFileVcf], shell=True)
    return success

def tabix(inFileGz, override=False, additionalParameters=[], exe='tabix'):
    """ run tabix """
    success = True
    idxFile = inFileGz + ".tbi"
    if(not files_exist(idxFile) or override):
        success = success and pipelineStep(inFileGz, idxFile, [exe] + additionalParameters +[inFileGz], shell=True)
    return success

def count_lines(gzfile, comment_char='#'):
    """ count lines in a gzipped file, excluding header lines starting with comment_char """
    with gzip.open(gzfile, 'rb') as f:
        # header
        for line in f:
            if not line.decode("utf-8").startswith(comment_char):
                break
        for count, line in enumerate(f):
            pass
    return count+2

def move_id_to_info_field(vcf_in, info_field_name, vcf_out, desc=None):
    """ move all ID entries to a new info field """
    if desc is None:
        desc=info_field_name
    vcf = pysam.VariantFile(vcf_in,'r')
    header=vcf.header
    header.add_line("##cmd=chort_varan.moveId2InfoField()")
    header.add_line("##INFO=<ID=%s,Number=1,Type=String,Description=\"%s\">"%(info_field_name, desc))
    out=pysam.VariantFile(vcf_out, mode="w", header=header)
    for record in vcf.fetch():
        record.info[info_field_name]=','.join(record.id.split(';'))
        record.id='.'
        written=out.write(record)
    out.close()
 
 
def add_contig_headers(vcf_in, ref_fasta, vcf_out):
    """ add missing contig headers """
    vcf = pysam.VariantFile(vcf_in,'r')
    header=vcf.header
    fa=pysam.Fastafile(ref_fasta)
    for c in [c for c in fa.references if c not in header.contigs]:
        header.add_line("##contig=<ID=%s,length=%i" % (c, fa.get_reference_length(c)))
    out=pysam.VariantFile(vcf_out, mode="w", header=header)
    for record in vcf.fetch():
        written=out.write(record)
    out.close()
 

