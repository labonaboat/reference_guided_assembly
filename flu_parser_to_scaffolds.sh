#/bin/sh

cdir=`pwd`
mkdir isolated_reads
mv NS1-NS2*.fastq isolated_reads/
mv M1-M2*.fastq isolated_reads/
mv PB1*.fastq isolated_reads/
mv PA-PA-X*.fastq isolated_reads/
mv HA*.fastq isolated_reads/
mv NP*.fastq isolated_reads/
mv NA*.fastq isolated_reads/
mv PB2*.fastq isolated_reads/
cd isolated_reads/
pigz *fastq

packagefastqs.sh

printf "\n\nSPAdes is running in the backgroud\n"
date

currentdir=`pwd`
for f in *; do
    cd $currentdir
    echo $f; cd ./$f
    if [ -f *R2* ]; then
        #trimreads.sh
        #cd trimmedreads
        # possibly k77 to complete a contig that is truncated, but appears k127 is best overall
        # the smaller k value sometimes creates false ends 
        (spades.py -t 32 -k 127 --careful -1 *_R1.fastq.gz -2 *_R2.fastq.gz -o ./
        #####
        mkdir alignment; cp *gz alignment; cp scaffolds.fasta alignment
        cd alignment

        ref=`ls | grep .fasta`
        forReads=`ls | grep _R1`
        revReads=`ls | grep _R2`

        r=`echo $ref | sed 's/\..*//'`
        n=`echo $revReads | sed 's/_.*//' | sed 's/\..*//'`

        samtools faidx $ref
        java -Xmx4g -jar ${picard} CreateSequenceDictionary REFERENCE=${ref} OUTPUT=${r}.dict
        bwa index $ref
        bwa mem -M -t 16 -R @RG"\t"ID:"$n""\t"PL:ILLUMINA"\t"PU:"$n"_RG1_UNIT1"\t"LB:"$n"_LIB1"\t"SM:"$n" $ref $forReads $revReads > $n.sam
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

        mv ${n}.readreference.fasta ../
        cd ..
        mv scaffolds.fasta original_scaffolds.fasta
        mv ${n}.readreference.fasta scaffolds.fasta) &
        #####
    else
        #trimreads-iontorrent.sh
        #cd trimmedreads 
        spades.py -k 45,47,49,51,53,57,59,61,63,65,77,99,127 --mismatch-correction -s *.fastq.gz -o ./ &
    fi
done; wait

cd $cdir/isolated_reads

mkdir scaffolds
for i in *; do
    cp $i/scaffolds.fasta scaffolds/${i}_scaffolds.fasta
done

mv scaffolds $cdir/scaffolds_from_isolated_reads
cd $cdir/scaffolds_from_isolated_reads

for i in *fasta; do
    fasta_size_selection.pl $i &
done
wait

rm out*
rm strand*

mkdir scaffolds notables
mv *scaffolds.fasta scaffolds/
mv *notable.fasta notables

subtype=Need_to_update

cd notables
name=`pwd | sed 's:/scaffolds_from_isolated_reads/notables::'`; sampleName=`basename $name`

#Format fasta headers for NCBI
genotypingcodes="/scratch/report/flu_genotyping_codes.txt"
echo "genotyping codes not given"
echo "sampleName: $sampleName"

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

    sed  's/>.*/>Seq1/' PB2_notable.fasta > ${sampleName}.temp
    sed  's/>.*/>Seq2/' PB1_notable.fasta >> ${sampleName}.temp
    sed  's/>.*/>Seq3/' PA-PA-X_notable.fasta >> ${sampleName}.temp
    sed  's/>.*/>Seq4/' HA_notable.fasta >> ${sampleName}.temp
    sed  's/>.*/>Seq5/' NP_notable.fasta >> ${sampleName}.temp
    sed  's/>.*/>Seq6/' NA_notable.fasta >> ${sampleName}.temp
    sed  's/>.*/>Seq7/' M1-M2_notable.fasta >> ${sampleName}.temp
    sed  's/>.*/>Seq8/' NS1-NS2_notable.fasta >> ${sampleName}.temp

# Create "here-document"
cat >./param.txt <<EOL

>Seq1 [organism=Influenza A virus](A/${species}/${state}/${sample}/${sampleyear}(${subtype})) segment 1, polymerase PB2 (PB2) gene, complete cds.
>Seq2 [organism=Influenza A virus](A/${species}/${state}/${sample}/${sampleyear}(${subtype})) segment 2, polymerase PB1 (PB1) gene, complete cds.
>Seq3 [organism=Influenza A virus](A/${species}/${state}/${sample}/${sampleyear}(${subtype})) segment 3, polymerase PA (PA) gene, complete cds.
>Seq4 [organism=Influenza A virus](A/${species}/${state}/${sample}/${sampleyear}(${subtype})) segment 4, hemagglutinin (HA) gene, complete cds.
>Seq5 [organism=Influenza A virus](A/${species}/${state}/${sample}/${sampleyear}(${subtype})) segment 5, nucleoprotein (NP) gene, complete cds.
>Seq6 [organism=Influenza A virus](A/${species}/${state}/${sample}/${sampleyear}(${subtype})) segment 6, neuraminidase (NA) gene, complete cds.
>Seq7 [organism=Influenza A virus](A/${species}/${state}/${sample}/${sampleyear}(${subtype})) segment 7, matrix protein 2 (M2) and matrix protein 1 (M1) genes, complete cds.
>Seq8 [organism=Influenza A virus](A/${species}/${state}/${sample}/${sampleyear}(${subtype})) segment 8, non-structural protein NS1 and non-structural protein NS2 (NS) gene, complete cds.

EOL
    awk 'NR==FNR{a[$1]=$0;next} ($1 in a){ print a[$1]; next}1' param.txt ${sampleName}.temp > ${sampleName}-submissionfile.fasta
else
    echo "metadata not available"
    sed  's/>.*/>Seq1/' PB2_notable.fasta > ${sampleName}_nometa.fasta
    sed  's/>.*/>Seq2/' PB1_notable.fasta >> ${sampleName}_nometa.fasta
    sed  's/>.*/>Seq3/' PA-PA-X_notable.fasta >> ${sampleName}_nometa.fasta
    sed  's/>.*/>Seq4/' HA_notable.fasta >> ${sampleName}_nometa.fasta
    sed  's/>.*/>Seq5/' NP_notable.fasta >> ${sampleName}_nometa.fasta
    sed  's/>.*/>Seq6/' NA_notable.fasta >> ${sampleName}_nometa.fasta
    sed  's/>.*/>Seq7/' M1-M2_notable.fasta >> ${sampleName}_nometa.fasta
    sed  's/>.*/>Seq8/' NS1-NS2_notable.fasta >> ${sampleName}_nometa.fasta
fi
pwd

#rm *temp
rm *information
rm param.txt

# 2016-09-17 stuber
