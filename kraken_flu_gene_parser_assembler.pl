#!/usr/bin/env perl

use strict;
use warnings;

use File::Basename;
use Data::Dumper qw(Dumper);
use Cwd;

# Usage:
    # virus_segment_finder.pl < full path to fastq files >
    # Files can be zipped or unzipped, single or paired reads
    # Paired files must be labed R1 and R2

print "----- START -----\n\n";



# one argument expected , which must be path to fastq files 
#my $path = $ARGV[0];
my $path = getcwd;
if (not defined $path) {
    die "ERROR!!! --> Argument to directory containing FASTQ/s was not provided"
}

# place a "/" at the end of path if one not provided
$path =~ s|/?$|/|;
#print "Path to files: $path\n";

# Determine if single or paired read
my @files = < ${path}*fastq* >;
my $read_type;
my $filename;
# check if array is empty
if (@files) {
    foreach (@files) {
        # get basename
        my $basename = basename($_);
        $filename=$basename;
        chomp $filename;
        # try to find R[1-2] from basename
        if ($basename =~ /R[1-2]{1}/g) {
            $read_type = "paired";
        } else {
            $read_type = "single";
        }
    }
} else {
    die "ERROR!!! FASTQ not found at $path Look";
}

print "Read type found: $read_type\n\n";

# make zipped and unzipped reads available
use IO::Uncompress::Gunzip qw(gunzip $GunzipError);
use IO::Compress::Gzip qw(gzip $GzipError);
my @files = < ${path}*fastq* >;
foreach my $file (@files) {
    # -B  File is a "binary" file
    # -T  File is an ASCII or UTF-8 text file
    # -e File exists
    # Test if FASTQs are zipped or not.  Make zip and unzip available
    # Binary and text file not present
    # strip .gz
    my $no_gz = substr($file, 0, -3);
        if ( -B $file && ! ( -e $no_gz )) {
        (my $file_noext = $file) =~ s/\.[^.]+$//;
        gunzip $file => $file_noext or die "ERROR!!! gunzip failed for $file";
        }
    my $with_gz = $file . ".gz";
        # Text file and binary not present
    if ( -T $file && ! ( -e $with_gz )) {
        gzip $file => "${file}.gz" or die "gzip failed: $GzipError";
        }
}

# allow to be globally available
# paired reads
my $input_R1_zip;
my $input_R2_zip;
my $input_R1_unzip;
my $input_R2_unzip;
# single reads
my $input_zip;
my $input_unzip;

if ($read_type eq "paired") {
    # place zip and unzipped file into variables
    my @input_R1_zip = < ${path}*R1*fastq.gz >;
    my @input_R2_zip = < ${path}*R2*fastq.gz >;
    my @input_R1_unzip = < ${path}*R1*fastq >;
    my @input_R2_unzip = < ${path}*R2*fastq >;

    foreach (@input_R1_zip) {
        $input_R1_zip = $_;
        }
    foreach (@input_R2_zip) {
        $input_R2_zip = $_;
        }
    foreach (@input_R1_unzip) {
        $input_R1_unzip = $_;
        }
    foreach (@input_R2_unzip) {
        $input_R2_unzip = $_;
        }

    print "Zipped Files:\n";
    print "$input_R1_zip\n";
    print "$input_R2_zip\n\n";

    print "Unzipped Files:\n";
    print "$input_R1_unzip\n";
    print "$input_R2_unzip\n\n";
} else {
    
    my @input_zip = < ${path}*fastq.gz >;
    my @input_unzip = < ${path}*fastq >;
    
    foreach (@input_zip) {
        $input_zip = $_;
    }
    foreach (@input_unzip) {
        $input_unzip = $_;
    }
    
    print "Zipped File:\n";
    print "$input_zip\n";
    
    print "Unzipped File:\n";
    print "$input_unzip\n\n";
    
}

my $krakenDatabase="/home/shared/databases/kraken/flu_jhu/fludb_20150820_with_hosts";

$filename =~ s/[_.].*//;
my $kraken_output = $path . $filename . "-output.txt";
my $kraken_report = $path . $filename . "-report.txt";
print "Working on file: $filename\n";
if ($read_type eq "paired") {
    print "Kraken has started with paired reads\n\n";
    `kraken --db $krakenDatabase --threads 50 --paired $input_R1_unzip $input_R2_unzip > $kraken_output`;
    `kraken-report --db $krakenDatabase $kraken_output > $kraken_report`;
} else {
    print "Kraken has started with single read\n\n";
    `kraken --db $krakenDatabase --threads 50 $input_unzip > $kraken_output`;
    `kraken-report --db $krakenDatabase $kraken_output > $kraken_report`;
}

my $jhu_out = $path . $filename . "-jhu-output.txt";
my $html_out = $path . $filename . ".html";

print "\nBuilding Krona Graph... using JHU kraken2krona.sh\n";
`kraken2krona.sh -i $kraken_output -k $krakenDatabase -o $jhu_out -r $html_out`;

###

#my $kraken_report = $ARGV[0];
open(my $kraken_report_handle, '<', $kraken_report) or die "Can't open $kraken_report";

#my $kraken_output = $ARGV[1];
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
    if ($read_type eq "paired") {
        my $R1out = "$path" . "$key" . "_R1.fastq";
        my $R2out = "$path" . "$key" . "_R2.fastq";
        print "fgrep R1\n"; 
        `fgrep --no-group-separator -A3 -h -f $outfile $input_R1_unzip > "$R1out"`;
        print "fgrep R2\n";
        `fgrep --no-group-separator -A3 -h -f $outfile $input_R2_unzip > "$R2out"`;
    } else {
        my $readout = "$path" . "$key" . ".fastq";
        print "fgrep reads\n";
        `fgrep --no-group-separator -A3 -h -f $outfile $input_unzip > "$readout"`;
    }
}

`flu_parser_to_scaffolds.sh`;



print "\nDONE\n\n";
### created 2016-08-24 tstuber
