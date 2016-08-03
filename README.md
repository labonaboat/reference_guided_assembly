## reference_guided_assembly
USDA APHIS Veterinary Services (VS) virus, mainly influenza, Kraken / Krona and reference guided assembly.

## OVERVIEW
**idvirus.sh** takes input as single or paired FASTQ reads.  Kraken ran using either standard, host or flu specific database.  Kraken results are converted to Krona graph.  Following Kraken pipeline aligns reads against select genomes based on given argument.  Consensus sequence is BLAST against NT database with the best hit used to alignment reads again.  This is repeated three times.  The resulting consensus is formated, and a LaTex pdf report is output with analysis summary.   

**idemail.sh** takes multiple samples at once and provides an additional summary report of all genomes.
