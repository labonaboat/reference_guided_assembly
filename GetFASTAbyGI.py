#!/usr/bin/env python

from Bio import SeqIO
from Bio import Entrez
from sys import argv

#provide a GI or accession number
script, acc, writeFile = argv

#print ("This is the GI number: %s " % gi)
#print ("This is the file to write too: %s " % writeFile)
#print ("script: %s " % script)

entrezDbName = 'nucleotide'
Entrez.email = 'tod.p.stuber@usda.gov'

entryData = Entrez.efetch(db=entrezDbName, id=acc, retmode="text", rettype='fasta')

local_file=open(writeFile,"w")
local_file.write(entryData.read())
entryData.close()
local_file.close()

