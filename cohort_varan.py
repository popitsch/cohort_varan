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
    Main file that does all commandline handling
"""
from argparse import ArgumentParser, RawDescriptionHelpFormatter
import os, sys
from utils import existing_file
from annotate_vcf import annotate_vcf
from annotate_cohort import annotate_cohort

# Necessary for including python modules from a parent directory
sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
 
__version__ = '1.0'
logo = """
.--------------------------.
| c o h o r t _ v a r a n  |\\
|   __{    }__ __{    }__  | |
|  || ( .. ) ||  ( .. ) || | |
|  __  \__/  __   \__/  __ | |
| /..\  ''  /..\   ''  /..\| |
.--------------------------. |
 \__________v_1.0___________\|
"""

usage = '''                           

  Copyright (C) 2020 XXX.  All rights reserved.

  Distributed on an "AS IS" basis without warranties
  or conditions of any kind, either express or implied.

USAGE
'''
if __name__ == '__main__':
    MODULES = ["annotate_vcf", "annotate_cohort"]
    #============================================================================    
    if len(sys.argv) <= 1 or sys.argv[1] in ['-h', '--help']:
        print("usage: cohort_varan.py [-h] " + ",".join(MODULES))
        sys.exit(1)
    mod = sys.argv[1]  
    if mod not in MODULES:
        print("Invalid module %s selected. Please use one of " + ",".join(MODULES))
        sys.exit(1)
        
    parser = {}
    parser["annotate_vcf"] = ArgumentParser(description=usage, formatter_class=RawDescriptionHelpFormatter)
    parser["annotate_vcf"].add_argument("-i", "--in", type=existing_file, required=True, dest="infile", metavar="infile", help="Input TSV sample sheet. Columns: ID,genome,source(platypus|deepvariant|gatk),vcf_file")
    parser["annotate_vcf"].add_argument("-o", "--outdir", type=str, required=True, dest="outdir", metavar="outdir", help="output directory")
    parser["annotate_vcf"].add_argument("-t", "--threads", type=int, dest="threads", default=1, help="maximum threads [default: %(default)s]")
    parser["annotate_vcf"].add_argument("-c", "--config", type=existing_file, required=True, dest="config_file", help="Configuration file [default: %(default)s]")
    parser["annotate_vcf"].add_argument("-s", "--snpEffConfig", type=existing_file, required=False, dest="snpEffConfig", help="Optional snpEff configuration file [default: %(default)s]")
    parser["annotate_vcf"].add_argument("--includeNonPass", required=False, action="store_true", default=False, dest="includeNonPass", help="If set, non-PASS variants will not be filtered")
    parser["annotate_vcf"].add_argument("--overwrite", required=False, action="store_true", default=False, dest="overwrite", help="If set, existing files will be overwritten")

    parser["annotate_cohort"] = ArgumentParser(description=usage, formatter_class=RawDescriptionHelpFormatter)
    parser["annotate_cohort"].add_argument("-c", "--config", type=existing_file, required=True, dest="config_file", metavar="config_file", help="JSON config file")
    parser["annotate_cohort"].add_argument("-o", "--outdir", type=str, required=True, dest="outdir", metavar="outdir", help="output directory")
    parser["annotate_cohort"].add_argument("-t", "--threads", type=int, dest="threads", default=1, help="maximum threads [default: %(default)s]")

    print(logo)
    print("module: " + mod)
    args = parser[mod].parse_args(sys.argv[2:])
    #============================================================================
    if mod == "annotate_vcf":
        annotate_vcf(args.config_file, args.infile, args.outdir, int(args.threads), args.snpEffConfig, args.includeNonPass, args.noQC, args.overwrite)
    if mod == "annotate_cohort":
        annotate_cohort(args.config_file, args.outdir, int(args.threads))

