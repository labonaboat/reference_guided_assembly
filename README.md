## reference_guided_assembly
USDA APHIS Veterinary Services (VS) virus, mainly influenza, Kraken / Krona and reference guided assembly.

## OVERVIEW
**idvirus.sh** takes input as single or paired FASTQ reads.  Kraken is ran using either standard, host or flu specific database.  Kraken results are converted to a Krona graph.  Reads are then aligned against reference genomes selected based on an input argument.  Consensus sequence is BLAST against NT database with the best hit used in another alignment.  This is repeated three times.  The resulting consensus is formated, and a LaTex report is output with analysis summary.   
**idemail.sh** run takes in multiple samples at once and provides an additional summary report of all genomes.
