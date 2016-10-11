#!/bin/sh

alias pause='read -p "$LINENO Enter"'
root=`pwd`

genotypingcodes="/scratch/report/flu_genotyping_codes.txt"
krakenDatabase="/home/shared/databases/kraken/flu_jhu/fludb_20150820_with_hosts"
sampleName=`ls *.fastq* | head -1 | sed 's/[_.].*//'`
printf "\n Working on $sampleName\n\n"

# NCBI downloaded reference location
mydb="/data/mydb"
#Delete files in local database that may be empty
dzdo chmod 755 * ${mydb}/*
for i in ${mydb}/*; do
    if [ -s $i ]; then
        echo "" > /dev/null
    else
        echo "file $i is empty and has been deleted"
        rm -f $i
    fi
done

#Establish Read Files
if [ -f *R2* ]; then
    export sampleType="paired"
    r1=${root}/`ls *_R1*`
    echo "R1 reads: $r1"
    r2=${root}/`ls *_R2*`
    echo "R2 reads: $r2"
    mkdir original_reads
    cp $r1 $r2 original_reads

    echo "*** Trimming"
    #######################
    #      Trimming       #
    #######################

    #Trim the reads with bbmap tool kit (bbduk plugin)
    #about twice as fast as trimmomatic

    strain=$(echo $revReads | sed 's/_.*//' | sed 's/\..*//')
    echo -e "Quality trimming sample "$strain""

    bbduk.sh -Xmx80g \
    in1=$r1 \
    in2=$r2 \
    ref="/usr/local/bin/bbmap/resources/nextera_flu.fa.gz" \
    ktrim=r k=23 mink=11 hdist=1 \
    qtrim=lr trimq=5 \
    minlen=36 \
    out1=${root}/${sampleName}_Trimmed_R1.fastq.gz \
    out2=${root}/${sampleName}_Trimmed_R2.fastq.gz \
    stats=trim_stats.txt \
    qchist=qc_by_base.txt \
    threads=auto \
    showspeed=f

    r1t="${root}/${sampleName}_Trimmed_R1.fastq.gz"
    echo "R1 reads to be used after trimmed: $r1t"
    r2t="${root}/${sampleName}_Trimmed_R2.fastq.gz"
    echo "R2 reads to be used after trimmed:: $r2t"

else
    echo "Just a single read present"
    export sampleType="single"
    single_read=${root}/*fastq*
    echo "Forward Reads:  $single_read"
fi

# make unzipped files available and keep zipped
pigz -dk *gz

#######################
#       Kraken        #
#######################

echo "Kraken database selected is: $krakenDatabase"
date
echo "*** Kraken is finding reads"
#Run Kraken
if [[ $sampleType == "paired" ]]; then
kraken --db ${krakenDatabase} --threads ${NR_CPUS} --paired ${r1t%.gz} ${r2t%.gz} > $sampleName-output.txt && kraken-report --db ${krakenDatabase} $sampleName-output.txt > $sampleName-report.txt
else
kraken --db ${krakenDatabase} --threads ${NR_CPUS}  ${single_read%.gz} > $sampleName-output.txt && kraken-report --db ${krakenDatabase} $sampleName-output.txt > $sampleName-report.txt
fi

date
echo "*** Transforming Kraken output to Krona graph"

echo "------> jhu Building Krona Graph..."
kraken2krona.sh -i $sampleName-output.txt -k ${krakenDatabase} -o $sampleName-jhu-output.txt -r $sampleName-jhu-Krona_id_graphic.html

mkdir ${root}/kraken
mv $sampleName-output.txt $sampleName-report.txt $sampleName-jhu-Krona_id_graphic.html ${root}/kraken

# Set variables and paths
koutput=${root}/kraken/$sampleName-output.txt
kreport=${root}/kraken/$sampleName-report.txt
khtml=${root}/kraken/$sampleName-jhu-Krona_id_graphic.html
echo "kraken output: $koutput"
echo "kraken report: $kreport"
echo "krona html: $khtml"

########################

cat >${root}/pick_fasta.pl <<'EOL'
#!/usr/bin/env perl

use strict;
use warnings;

use Bio::Seq;
use Bio::SeqIO;

my $in_file = $ARGV[0];
if (not defined $in_file) {
    die "ERROR!!! --> An input file was not provided.  Provide FASTA!"
}

print "my in_file: $in_file\n";
(my $file_name = $in_file) =~ s/_scaffolds//;
print "File name: $file_name\n";

my $seqin_obj = Bio::SeqIO->new(-file => $in_file, -format => 'fasta');
my $seqout_obj = Bio::SeqIO->new(-file => ">$file_name", -fomat => 'fasta');

# print out first (and largest) seq obj
# isolate just one contig if multiple
my $inseq = $seqin_obj->next_seq;
$seqout_obj->write_seq($inseq);

`blastn -query $file_name -out strand.$file_name -db /data/BLAST/db/nt -max_target_seqs 1 -num_threads 16 -outfmt "6 sstrand"`;

# Put sequence into correct orientation
open(my $sh, '<', "strand.$file_name") or die "Can't open strand.$file_name";
my $first_line = <$sh>;
chomp $first_line;
my $seqbackin_obj = Bio::SeqIO->new(-file => "$file_name", -format => 'fasta');
print "This is the first line $first_line\n";

while (my $seq_obj = $seqbackin_obj->next_seq ){
    my $seqout_correct_orientation_obj = Bio::SeqIO->new(-file => ">$file_name", -fomat => 'fasta');
    if ($first_line eq "minus"){
        print "The reverse complement will be used\n";
        $seqout_correct_orientation_obj->write_seq($seq_obj->revcom);
    }elsif($first_line eq "plus"){
        print "Strand orientation is correct\n";
        $seqout_correct_orientation_obj->write_seq($seq_obj);
    }else{
        die "ERROR!!! --> Correct strand orientation not found";
    }
}
EOL

chmod 755 ${root}/pick_fasta.pl
#######################

function parse_kraken_report() {

echo "perl here-doc to parse report"
cat >${root}/parse_kraken_report.pl <<'EOL'
#!/usr/bin/env perl

my $kraken_report = $ARGV[0];
open(my $kraken_report_handle, '<', $kraken_report) or die "Can't open $kraken_report";

my $kraken_output = $ARGV[1];
open(my $kraken_output_handle, '<', $kraken_output) or die "Can't open $kraken_output";

# create hash
#my %taxon_segment;
# iterate each line of file
my @gene_taxon;
my $segment = "no_segment_info";
while ( <$kraken_report_handle> ) {
    chomp;
    # remove spaces at beginning of line
    s/^ +//;
    # ids columned on two spaces
    # find/replace two spaces with tabs
    s/ {2}/\t/g;
    # split file on tabs
    my @line = split(/\t/);
    # turn warnings off for the undef columns
    no warnings 'uninitialized';
    # column 11 contains segment id
    if ($line[12] eq "") {
        # if column empty use last segment id

        #print "$segment $line[4] \n";
        push @gene_taxon, "$segment $line[4]";
    } else {
        # if column not empty update segment id
        $segment = $line[12];
        push @gene_taxon, "$segment $line[4]";
        #print "$segment $line[4] \n";
    }
}

@gene_taxon = grep /M1,M2|NS1,NS2|NA|HA|NP|PA,PA-X|PB2|PB1/, @gene_taxon;

#foreach (@gene_taxon) {
#    print "$_ \n";
#}

# array to hash
my %gene_taxon_hash;
foreach (@gene_taxon) {
    my @gene_taxon_split = split(/ /);
    # add each line to hash
    # key:gene $gene_taxon_split[0], value:reference to taxon array $gene_taxon_split[1]
    push @ { $gene_taxon_hash {$gene_taxon_split[0]} }, $gene_taxon_split[1];
}
#print Dumper \%gene_taxon_hash;

# place columns 1 and 2 into array, delimited with tab
# column1: header
# column2: taxon identification
my @header_array;
while ( <$kraken_output_handle> ) {
    chomp;
    my @line = split(/\t/);
    push @header_array, "$line[1] \t $line[2]";
}

print "\nParsing reads by gene\n\n";
# for each $key/gene make a file holding headers
# these header can than be used to fgrep FASTQ files
# %gene_taxon_hash -> key:gene => value:reference to array of taxon
while (my ($key, $value) = each %gene_taxon_hash) {
    # for each value of each key
    print "\nWorking on $key\n";
    $key =~ s/,/-/g;
    my $outfile = "$path" . "$key" . "_headers.txt";
    foreach my $taxon_by_gene (@$value) {
        # look for a match, taxon to taxon
        # $header_array -> list of read header \t taxon id
        foreach (@header_array) {
            my @header_split = split /\t/;
            # $header_split[0]=read header, $header_split[1]=taxon id
            # put header finding
            # if taxon with gene label matches to taxon with header label
            #print "1: $taxon_by_gene, 2: $header_split[1]\n";
            if ($taxon_by_gene == $header_split[1]) {
                open (my $outfileh, '>>', $outfile) or die "Could not open file '$outfile' $!";
                my $header_line = $header_split[0];
                $header_line =~ s/ $//;
                $header_line =~ s/^ //;
                print $outfileh "$header_line\n";
                # check line --> print "$header_split[0] \t $header_split[1] \t $key \t $taxon_by_gene\n";
            }
        }
    }
}

EOL

chmod 755 ${root}/parse_kraken_report.pl

${root}/parse_kraken_report.pl $kreport $koutput

}

parse_kraken_report

mkdir header_files
mv *headers.txt header_files
cd header_files

function parse_reads() {
header_name=${each_header%_headers.txt}
 
if [[ $sampleType == "paired" ]]; then
    echo "Get the R1 reads for $header_name"
    fgrep --no-group-separator -A3 -h -f ${each_header} ${r1t%.gz} > ${header_name}_R1.fastq
    echo "Get the R2 reads for $header_name"
    fgrep --no-group-separator -A3 -h -f ${each_header} ${r2t%.gz} > ${header_name}_R2.fastq 
    mkdir -p ${root}/isolated_reads
    mkdir ${root}/isolated_reads/${header_name}
    mv ${header_name}_R1.fastq ${header_name}_R2.fastq ${root}/isolated_reads/${header_name}
    pigz ${root}/isolated_reads/${header_name}/*fastq
    r1i="${root}/isolated_reads/${header_name}/${header_name}_R1.fastq.gz"
    r2i="${root}/isolated_reads/${header_name}/${header_name}_R2.fastq.gz"
    
    cd ${root}/isolated_reads/${header_name}
    spades.py -t 32 -k 127 --careful -1 $r1i -2 $r2i -o ${root}/isolated_reads/${header_name}/

    cp scaffolds.fasta ${sampleName}_${header_name}_scaffolds.fasta
    
    ${root}/pick_fasta.pl ${sampleName}_${header_name}_scaffolds.fasta
    # update the scaffold file with single contig with correct orientation
    mv ${sampleName}_${header_name}.fasta ${sampleName}_${header_name}_scaffolds.fasta

    ref="${root}/isolated_reads/${header_name}/${sampleName}_${header_name}_scaffolds.fasta"

    r=`basename $ref | sed 's/\.fasta//'`
    n="${sampleName}_${header_name}"

    echo "n: $n"
    echo "r: $r"
    echo "r1i: $r1i"
    echo "r2i: $r2i"
    echo "ref: $ref"

    samtools faidx $ref
    java -Xmx4g -jar ${picard} CreateSequenceDictionary REFERENCE=${ref} OUTPUT=${r}.dict
    bwa index $ref
    bwa mem -M -t 16 -R @RG"\t"ID:"$n""\t"PL:ILLUMINA"\t"PU:"$n"_RG1_UNIT1"\t"LB:"$n"_LIB1"\t"SM:"$n" $ref $r1i $r2i > $n.sam
    samtools view -bh -T $ref $n.sam > $n.all.bam
    samtools view -h -F4 $n.all.bam > $n.mappedReads.sam
    samtools view -bh -F4 -T $ref  $n.mappedReads.sam > $n.raw.bam
    samtools sort $n.raw.bam -o $n.sorted.bam
    samtools index $n.sorted.bam
    java -Xmx4g -jar  ${picard} MarkDuplicates INPUT=$n.sorted.bam OUTPUT=$n.dup.bam METRICS_FILE=$n.FilteredReads.xls ASSUME_SORTED=true REMOVE_DUPLICATES=true
    samtools index $n.dup.bam

    java -Xmx4g -jar ${gatk} -R $ref -T HaplotypeCaller -I $n.dup.bam -o $n.hapreadyAll.vcf -dontUseSoftClippedBases -allowNonUniqueKmersInRef
    grep '^#' $n.hapreadyAll.vcf > header
    grep -v '^#' $n.hapreadyAll.vcf | awk 'BEGIN{FS="\t"}{if ($6 > 750) print $0}' > body
    mv $n.hapreadyAll.vcf $n.original_hapreadyAll.vcf
    cat header body > $n.hapreadyAll.vcf

    java -jar ${gatk} -T FastaAlternateReferenceMaker -R $ref -o ${n}.readreference.fasta -V $n.hapreadyAll.vcf

else
    echo "get the single read setup at line $LINENO"
fi

}

for each_header in *headers.txt; do
    cd ${root}/header_files
    parse_reads #&
done
wait
pause
cd ${root}/isolated_reads

pwd
date
printf "\nDONE\n\n";
# 2016-10-10 stuber

