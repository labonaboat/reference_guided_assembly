#!/bin/sh

starttime=$(date +%s)
alias pause='read -p "$LINENO Enter"'
root=`pwd`
pythonGetFasta=`which GetFASTAbyGI.py`

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

function getfluname () {

echo "using xlrd to get flu genotyping codes from MiSeq_Samples.xlsx"
date

cat >./excelcolumnextract.py <<EOL
#!/usr/bin/env python

import xlrd
from sys import argv

script, input = argv

wb = xlrd.open_workbook(input)

sh = wb.sheet_by_index(0)
for rownum in range(sh.nrows):
    row = sh.row_values(rownum)
    # use the join/map/str to output comma delimited list only
    print (', '.join(map(str, row)))

EOL

chmod 755 ./excelcolumnextract.py

./excelcolumnextract.py /bioinfo11/MKillian/MiSeq\ samples/MiSeq_Samples.xlsx | sed 's/, /,/g' | awk 'BEGIN{FS=","; OFS="\t"} {print $2, $4, $5}' | sed -e 's/[.*:()/\?]/_/g' -e 's/ /_/g' -e 's/_-/_/' -e 's/-_/_/' -e 's/__/_/g' -e 's/[_-]$//' > /scratch/report/flu_genotyping_codes.txt 

rm ./excelcolumnextract.py

}
getfluname

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
# BLAST1
    echo "short BLAST"
    blastn -query ${n}.readreference.fasta -db /data/BLAST/db/nt -num_threads 20 -out ${n}-blast1.txt -max_target_seqs 1 -outfmt "6 saccver"
    #rm ${refname}.readreference.fasta

    head -1 ${n}-blast1.txt > ${n}-blast1.txt.temp
    mv ${n}-blast1.txt.temp ${n}-blast1.txt

    if [ -s ${n}-blast1.txt ]; then
        echo "Something was found"
    else
        echo "No matches" > ERROR-REPORT.txt
        #exit 1
    fi

    acc=`head -1 ${n}-blast1.txt | sed 's/\..*//'`
    writeFile="${acc}.fasta"
    echo "acc $acc"
    echo "writeFile $writeFile"
    python -u "${pythonGetFasta}" $acc $writeFile
    ls ${mydb} >  ${mydb}/list.txt
    p=`grep "${acc}" ${mydb}/list.txt`

    if [[ -z "$p" ]]; then
        echo "Downloading from NCBI"
        echo "This is the file to write to:  $writeFile"
        echo "Running python script to grab fasta file"
        grep -v "Resource temporarily unavailable" $writeFile > $writeFile.temp; mv $writeFile.temp $writeFile
        if [ -s $writeFile ]; then
            echo "Downloaded from NCBI, Good to continue."
        else
            echo "Try downloading again"
            sleep 20
            echo "Running python script to grab fasta file"
            python -u "${pythonGetFasta}" $acc $writeFile
            grep -v "Resource temporarily unavailable" $writeFile > $writeFile.temp; mv $writeFile.temp $writeFile
            sleep 5
            if [ ! -s $writeFile ]; then
                echo "Downloaded from NCBI, Good to continue."
            else
                echo "Try downloading again"
                sleep 120
                echo "Running python script to grab fasta file"
                python -u "${pythonGetFasta}" $acc $writeFile
                grep -v "Resource temporarily unavailable" $writeFile > $writeFile.temp; mv $writeFile.temp $writeFile
                sleep 5
                if [ ! -s $writeFile ]; then
                    echo "Downloaded from NCBI, Good to continue."
                else
                    echo "Try downloading again"
                    sleep 320
                    echo "Running python script to grab fasta file"
                    python -u "${pythonGetFasta}" $acc $writeFile
                    grep -v "Resource temporarily unavailable" $writeFile > $writeFile.temp; mv $writeFile.temp $writeFile
                    sleep 5
                    if [ ! -s $writeFile ]; then
                        read -p "Fasta file ${acc} failed to download from NCBI, Manually download if desired to proceed.  Line: $LINENO"
                    fi
                fi
            fi
        fi
        cp $writeFile ${mydb}
    else
        echo "File is local"
        cp ${mydb}/${p} ./
        writeFile=${p}
    fi

    ref="$writeFile"
    r=`basename $ref | sed 's/\.fasta//'`
    n="${sampleName}_${header_name}"
echo "writefile $writeFile"

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
# BLAST2
    echo "short BLAST"
    blastn -query ${n}.readreference.fasta -db /data/BLAST/db/nt -num_threads 20 -out ${n}-blast1.txt -max_target_seqs 1 -outfmt "6 saccver"
    #rm ${refname}.readreference.fasta

    head -1 ${n}-blast1.txt > ${n}-blast1.txt.temp
    mv ${n}-blast1.txt.temp ${n}-blast1.txt

    if [ -s ${n}-blast1.txt ]; then
        echo "Something was found"
    else
        echo "No matches" > ERROR-REPORT.txt
        #exit 1
    fi

    acc=`head -1 ${n}-blast1.txt | sed 's/\..*//'`
    writeFile="${acc}.fasta"
    echo "acc $acc"
    echo "writeFile $writeFile"
    python -u "${pythonGetFasta}" $acc $writeFile
    ls ${mydb} >  ${mydb}/list.txt
    p=`grep "${acc}" ${mydb}/list.txt`

    if [[ -z "$p" ]]; then
        echo "Downloading from NCBI"
        echo "This is the file to write to:  $writeFile"
        echo "Running python script to grab fasta file"
        grep -v "Resource temporarily unavailable" $writeFile > $writeFile.temp; mv $writeFile.temp $writeFile
        if [ -s $writeFile ]; then
            echo "Downloaded from NCBI, Good to continue."
        else
            echo "Try downloading again"
            sleep 20
            echo "Running python script to grab fasta file"
            python -u "${pythonGetFasta}" $acc $writeFile
            grep -v "Resource temporarily unavailable" $writeFile > $writeFile.temp; mv $writeFile.temp $writeFile
            sleep 5
            if [ ! -s $writeFile ]; then
                echo "Downloaded from NCBI, Good to continue."
            else
                echo "Try downloading again"
                sleep 120
                echo "Running python script to grab fasta file"
                python -u "${pythonGetFasta}" $acc $writeFile
                grep -v "Resource temporarily unavailable" $writeFile > $writeFile.temp; mv $writeFile.temp $writeFile
                sleep 5
                if [ ! -s $writeFile ]; then
                    echo "Downloaded from NCBI, Good to continue."
                else
                    echo "Try downloading again"
                    sleep 320
                    echo "Running python script to grab fasta file"
                    python -u "${pythonGetFasta}" $acc $writeFile
                    grep -v "Resource temporarily unavailable" $writeFile > $writeFile.temp; mv $writeFile.temp $writeFile
                    sleep 5
                    if [ ! -s $writeFile ]; then
                        read -p "Fasta file ${acc} failed to download from NCBI, Manually download if desired to proceed.  Line: $LINENO"
                    fi
                fi
            fi
        fi
        cp $writeFile ${mydb}
    else
        echo "File is local"
        cp ${mydb}/${p} ./
        writeFile=${p}
        fi

    ref="$writeFile"
    r=`basename $ref | sed 's/\.fasta//'`
    n="${sampleName}_${header_name}"

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

    java -jar ${gatk} -T FastaAlternateReferenceMaker -R $ref -o ${n}.final-readreference.fasta -V $n.hapreadyAll.vcf

##################
else
    echo "get the single read setup at line $LINENO"
fi

}

for each_header in *headers.txt; do
    cd ${root}/header_files
    parse_reads &
done
wait
cd ${root}
mkdir final_assemblies
find . -name "*final-readreference.fasta" -exec cp {} final_assemblies \;
cd ${root}/final_assemblies

#####################
#Format fasta headers for NCBI
genotypingcodes="/scratch/report/flu_genotyping_codes.txt"
$subtype="###UPDATE###"
grep `echo $sampleName | sed 's/_denovo//'` $genotypingcodes | head -1 > ${sampleName}.information
if [[ -s ${sampleName}.information ]]; then
    echo "file exists and is greater than zero, continue"
    #column 1: sample
    sample=`awk 'BEGIN{FS="\t"}{print $1}' ${sampleName}.information`
    echo "sample $sample"

    #column 2: species
    species=`awk 'BEGIN{FS="\t"}{print $2}' ${sampleName}.information`
    speciesspace=`awk 'BEGIN{FS="\t"}{print $2}' ${sampleName}.information | sed 's/_/ /g'`
    echo "species $species"
    echo "speciesspace $speciesspace"

    #column 3: state
    state=`awk 'BEGIN{FS="\t"}{print $3}' ${sampleName}.information`
    echo "state $state"
    statespace=`awk 'BEGIN{FS="\t"}{print $3}' ${sampleName}.information | sed 's/_/ /g'`
    syear=`echo "$sample" | sed 's/-.*//'`
    echo "syear $syear"
    sampleyear=`echo "20${syear}"`
    echo "sampleyear $sampleyear"

    sed  's/>.*/>Seq1/' *PB2.final-readreference.fasta > ${sampleName}.temp
    sed  's/>.*/>Seq2/' *PB1.final-readreference.fasta >> ${sampleName}.temp
    sed  's/>.*/>Seq3/' *PA-PA-X.final-readreference.fasta >> ${sampleName}.temp
    sed  's/>.*/>Seq4/' *HA.final-readreference.fasta >> ${sampleName}.temp
    sed  's/>.*/>Seq5/' *NP.final-readreference.fasta >> ${sampleName}.temp
    sed  's/>.*/>Seq6/' *NA.final-readreference.fasta >> ${sampleName}.temp
    sed  's/>.*/>Seq7/' *M1-M2.final-readreference.fasta >> ${sampleName}.temp
    sed  's/>.*/>Seq8/' *NS1-NS2.final-readreference.fasta >> ${sampleName}.temp

# Create "here-document"
cat >./param.txt <<EOL

>Seq1-test [organism=Influenza A virus](A/${species}/${state}/${sample}/${sampleyear}(${subtype})) segment 1, polymerase PB2 (PB2) gene, complete cds.
>Seq2-test [organism=Influenza A virus](A/${species}/${state}/${sample}/${sampleyear}(${subtype})) segment 2, polymerase PB1 (PB1) gene, complete cds.
>Seq3-test [organism=Influenza A virus](A/${species}/${state}/${sample}/${sampleyear}(${subtype})) segment 3, polymerase PA (PA) gene, complete cds.
>Seq4-test [organism=Influenza A virus](A/${species}/${state}/${sample}/${sampleyear}(${subtype})) segment 4, hemagglutinin (HA) gene, complete cds.
>Seq5-test [organism=Influenza A virus](A/${species}/${state}/${sample}/${sampleyear}(${subtype})) segment 5, nucleoprotein (NP) gene, complete cds.
>Seq6-test [organism=Influenza A virus](A/${species}/${state}/${sample}/${sampleyear}(${subtype})) segment 6, neuraminidase (NA) gene, complete cds.
>Seq7-test [organism=Influenza A virus](A/${species}/${state}/${sample}/${sampleyear}(${subtype})) segment 7, matrix protein 2 (M2) and matrix protein 1 (M1) genes, complete cds.
>Seq8-test [organism=Influenza A virus](A/${species}/${state}/${sample}/${sampleyear}(${subtype})) segment 8, non-structural protein NS1 and non-structural protein NS2 (NS) gene, complete cds.

EOL
    awk 'NR==FNR{a[$1]=$0;next} ($1 in a){ print a[$1]; next}1' param.txt ${sampleName}.temp > ${sampleName}-test-submissionfile.fasta
    cp ${sampleName}-test-submissionfile.fasta $root
else
    echo "metadata not available"
    sed  's/>.*/>Seq1-test/' *PB2.final-readreference.fasta > ${sampleName}-test_nometa.fasta
    sed  's/>.*/>Seq2-test/' *PB1.final-readreference.fasta >> ${sampleName}-test_nometa.fasta
    sed  's/>.*/>Seq3-test/' *PA-PA-X.final-readreference.fasta >> ${sampleName}-test_nometa.fasta
    sed  's/>.*/>Seq4-test/' *HA.final-readreference.fasta >> ${sampleName}-test_nometa.fasta
    sed  's/>.*/>Seq5-test/' *NP.final-readreference.fasta >> ${sampleName}-test_nometa.fasta
    sed  's/>.*/>Seq6-test/' *NA.final-readreference.fasta >> ${sampleName}-test_nometa.fasta
    sed  's/>.*/>Seq7-test/' *M1-M2.final-readreference.fasta >> ${sampleName}-test_nometa.fasta
    sed  's/>.*/>Seq8-test/' *NS1-NS2.final-readreference.fasta >> ${sampleName}-test_nometa.fasta

    cp ${sampleName}-test_nometa.fasta $root
fi
pwd

#rm *temp
rm *information
rm param.txt
##################
cd $root
rm *fastq
rm -r final_assemblies
rm -r header_files
rm -r isolated_reads
rm *pl
rm *txt
rm *Trimmed*

endtime=`date +%s`
runtime=$((endtime-starttime))
printf 'Runtime: %dh:%dm:%ds\n' $(($runtime/3600)) $(($runtime%3600/60)) $(($runtime%60))
printf "\nDONE\n\n";
# 2016-10-10 stuber

