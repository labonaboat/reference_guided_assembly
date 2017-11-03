#!/usr/bin/env python

import os
import sys
from shutil import copyfile
import glob
import argparse
import textwrap

working_directory_at_start = os.getcwd()

parser = argparse.ArgumentParser(prog='PROG', formatter_class=argparse.RawDescriptionHelpFormatter, description=textwrap.dedent('''\

---------------------------------------------------------
rgAssembly --> *Provide Usage Instructions*

'''), epilog='''---------------------------------------------------------''')

parser.add_argument('-m', '--email', action='store', dest='email', help='[**APHIS only**, specify own SMTP address for functionality] email options: all, s, tod, jess, suelee, chris, email_address')
parser.add_argument('-s', '--species', action='store', dest='in_species', help='OPTIONAL: Seed using species designation')
parser.add_argument('-f', '--file', action='store', dest='in_file', help='OPTIONAL: Seed using single file designation')
parser.add_argument('-d', '--dir', action='store', dest='in_dir', help='OPTIONAL: Seed using all FASTAs within directory')

args = parser.parse_args()
print ("\nSET ARGUMENTS: ")
print (args)
print("")

if args.email == "all":
    email_list = "tod.p.stuber@aphis.usda.gov, Jessica.A.Hicks@aphis.usda.gov, Christine.R.Quance@aphis.usda.gov, Suelee.Robbe-Austerman@aphis.usda.gov, patrick.m.camp@aphis.usda.gov, David.T.Farrell@aphis.usda.gov, Robin.L.Swanson@aphis.usda.gov, hannah.m.tharp@aphis.usda.gov, Doris.M.Bravo@aphis.usda.gov, eto3@cdc.gov"
elif args.email == "tod":
    email_list = "tod.p.stuber@aphis.usda.gov"
elif args.email == "jess":
    email_list = "Jessica.A.Hicks@aphis.usda.gov"
elif args.email == "suelee":
    email_list = "tod.p.stuber@aphis.usda.gov, Jessica.A.Hicks@aphis.usda.gov, Suelee.Robbe-Austerman@aphis.usda.gov"
elif args.email == "chris":
    email_list = "tod.p.stuber@aphis.usda.gov, Jessica.A.Hicks@aphis.usda.gov, Christine.R.Quance@aphis.usda.gov, Suelee.Robbe-Austerman@aphis.usda.gov"
else:
    email_list = "tod.p.stuber@aphis.usda.gov"

############
def get_seed_list(dir_to_seed):
    if not dir_to_seed.endswith('/'):
        dir_to_seed = dir_to_seed + "/"
    all_fna = glob.glob(dir_to_seed + "*.fna")
    all_fa = glob.glob(dir_to_seed + "*.fa")
    all_fas = glob.glob(dir_to_seed + "*.fas")
    all_fasta = glob.glob(dir_to_seed + "*.fasta")
    seed_list = all_fna + all_fa + all_fas + all_fasta
    return seed_list

if args.in_species:
    if args.in_species == "bvd":
        dir_to_seed = "/bioinfo11/MKillian/Analysis/script_dependents/bvd/"
        seed_list = get_seed_list(dir_to_seed)
        print("Seed using species designation %s" % args.in_species)
        print("There are this many FASTA: %s" % len(seed_list))
    else:
        print("Incorrect species option given, unable to find argument listed as a species")
        print("Argument given %s:  " % args.in_species)
        sys.exit(0)
elif args.in_file:
    seed_list = []
    seed_list.append(args.in_file)
    print("Seed using single file designation %s" % args.in_file)
    print("There are this many FASTA: %s" % len(seed_list))
elif args.in_dir:
    print(args.in_dir)
    seed_list = get_seed_list(args.in_dir)
    print("Seed using all FASTAs within directory %s" % args.in_dir)
    print("There are this many FASTA: %s" % len(seed_list))
else:
    run_assembly = True
    print("run_assembly == %s" % run_assembly)
    print("for testing though... exiting")
    sys.exit(0)

for i in seed_list:
    print(i)

print("\nDONE\n")

# class to align reads, determine best match



# Created November 2017 by Tod Stuber
