import os
import sys
import glob
import __init__
import argparse
import subprocess




trim_galore = '/nethome/sravishankar/tools/trim_galore'
bwa_meth = '/nethome/sravishankar/tools/bwa-meth-0.10/bwameth.py'
picard = '/nethome/sravishankar/tools/picard-tools-1.119/'
arg = picard+ 'AddOrReplaceReadGroups.jar'
mbam = picard + 'MergeSamFiles.jar'
md = picard + 'MarkDuplicates.jar'
bissnp = '/nethome/sravishankar/tools/BisSNP-0.82.2.jar'
java = 'java'
