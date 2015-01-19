#!/usr/bin/perl -w
use strict;

my $geneWiseInFile = "";
my $geneWiseProteinFile = "";
my $geneWiseDNAFile = "";
my $geneWise3LevelGFF3File = "";
my $geneWisePseudoGeneDNAFile = "";
my $summaryOutputFile = "";
my %possiblyTruncated = ();
foreach my $arg (@ARGV)
{
        ($arg =~ m/^-(h|help)$/) and usage();
        ($arg =~ m/^-g=(.+)$/) and $geneWiseInFile = $1;
}

if (-e $geneWiseInFile)
{
        $geneWiseProteinFile = $geneWiseInFile . ".faa";
        $geneWiseDNAFile = $geneWiseInFile . ".fna";
        $geneWise3LevelGFF3File = $geneWiseInFile . ".gff3";
        $geneWisePseudoGeneDNAFile = $geneWiseInFile . ".pseudo.fna";
}
($summaryOutputFile) or $summaryOutputFile = $geneWiseInFile . ".sum";

my $oldSep = $/;
$/= "genewise ";
open(IN, "<$geneWiseInFile") or die "Can not open genewise file $geneWiseInFile $!\n";
open(FAA, ">$geneWiseProteinFile") or die "Can not open genewise faa outfile $geneWiseProteinFile $!\n";
open(FNA, ">$geneWiseDNAFile") or die "Can not open genewise fna outfile $geneWiseDNAFile $!\n";
open(FNA, ">$geneWiseDNAFile") or die "Can not open genewise fna outfile $geneWiseDNAFile $!\n";
open(PSEUDO, ">$geneWisePseudoGeneDNAFile") or die "Can not open genewise fna outfile $geneWisePseudoGeneDNAFile $!\n";
open(GFF3, ">$geneWise3LevelGFF3File") or die "Can not open genewise fna outfile $geneWise3LevelGFF3File $!\n";
open(SUM, ">$summaryOutputFile") or die "Can not open genewise fna outfile $summaryOutputFile $!\n";
my $recordCounter = 0;
while ( my $record = <IN>)
{
        next if ($record eq "genewise ");
        $recordCounter++;
        if ($recordCounter % 100 == 0)
        {
                select STDERR; print STDERR "$recordCounter\n";
        }
        my @lines = split/\n/, $record;
        my $getTranslation = 0;
        my $getDNASequence = 0;
        my $getGFF3 = 0;
        my $summaryLine = "";
        my $faaSeq = "";
        my $fnaSeq = "";
        my $isPseudoGene = 0;
        my $lastmRNAParentId = "";
        my $lastGeneId = "";

        if ($record =~ m/is a pseudo gene/s)
        {
                $isPseudoGene = 1;
        }
        my @gff3 =();
        my $lineCount = 0;
        my $summaryTxt = "";
        foreach my $i (0 .. $#lines)
        {
                my $line = $lines[$i];
                chomp $line;
                my @lineData = split/\t/, $line;
                if ($line =~ m/^Bits\tQuery/)
                {
                        $summaryLine = $lines[$i + 1];
                        next;
                }
                elsif ($line =~ m/^>.+\.sp\.tr$/)
                {
                        #>scaffold00001.[695800:696851].sp.tr
                        $getTranslation = 1;
                        my $nameLine = $line;
                        #$nameLine =~ s/\.tr$//g; # replace terminal .tr
                        $nameLine =~ s/\./_/g; # replace .
                        $nameLine =~ s/\:/_/g; # replace :
                        $nameLine =~ s/\[//g; # replace [
                        $nameLine =~ s/\]//g; # replace ]
                        my @nameData = split /\s+/, $nameLine;
                        $nameData[0] = $nameData[0] . "_" . $recordCounter;
                        $nameLine = join(" ", @nameData);
                        $faaSeq = $nameLine . " " . $summaryLine . "\n";
                        next;
                }
                elsif ($line =~ m/^>.+\.sp$/)
                {
                        $getDNASequence = 1;
                        my $nameLine = $line;
                        $nameLine =~ s/\./_/g; # replace .
                        $nameLine =~ s/\:/_/g; # replace :
                        $nameLine =~ s/\[//g; # replace [
                        $nameLine =~ s/\]//g; # replace ]
                        my @nameData = split /\s+/, $nameLine;
                        $nameData[0] = $nameData[0] . "_" . $recordCounter;
                        $lastGeneId = $nameData[0];
                        $lastGeneId =~ s/^>//;

                        $nameLine = join(" ", @nameData);
                        $fnaSeq = $nameLine . " " . $summaryLine . "\n";
                        next;
                }
                elsif (($line =~ m!^//$!) and ($i < $#lines))
                {
                        if (($isPseudoGene == 0) and ($lines[$i + 1] =~ m/^>.+\.sp\.tr$/))
                        {
                                $getTranslation = 1;
                                $getDNASequence = 0;
                                $getGFF3 = 0;
                                next;
                        }
                        elsif ($lines[$i + 1] =~ m/^>.+\.sp$/)
                        {
                                $getTranslation = 0;
                                $getDNASequence = 1;
                                $getGFF3 = 0;
                                next;
                        }
                        elsif ($getDNASequence)
                        {
                                $getTranslation = 0;
                                $getDNASequence = 0;
                                $getGFF3 = 1;
                                next;
                        }
                        else
                        {
                                $getDNASequence = 0;
                                $getTranslation = 0;
                                $getGFF3 = 0;
                                next;
                        }
                        next;
                }
                elsif ($getTranslation)
                {
                        $faaSeq .= $line;
                        next;
                }
                elsif ($getDNASequence)
                {
                        $fnaSeq .= $line;
                        next;
                }
                elsif ($getGFF3 and ($#lineData == 8))
                {
                        # GFF3 line
                        if ($lineData[2] eq 'match')
                        {
                                # scaffold00001   GeneWise        match   695800  696851  140.81  +       .       scaffold00001-genewise-prediction-1

                                #my $geneIdAttr = "ID=" . $lineData[0] . "_" . $lineData[3] . "_" . $lineData[4] . "_sp_" . $recordCounter;
                                my $geneIdAttr = "ID=$lastGeneId";
                                #my $mrnaIdAttr = "ID=" . $lineData[0] . "_" . $lineData[3] . "_" . $lineData[4] . ".mrna";
                                my $mrnaIdAttr = "ID=$lastGeneId" . "_mrna";
                                $lastmRNAParentId = $mrnaIdAttr;
                                my $nameAttr = "Name=$lastGeneId";
                                my $parentAttr = "Parent=$geneIdAttr";
                                my $pseudoGeneAttr = "";
                                if ($isPseudoGene == 1)
                                {
                                        $pseudoGeneAttr = ";IsAPseudogene=1";
                                }
                                push(@gff3, join("\t", $lineData[0], $lineData[1], 'gene', $lineData[3], $lineData[4], $lineData[5], $lineData[6], $lineData[7], "$geneIdAttr;$nameAttr$pseudoGeneAttr") );
                                push(@gff3, join("\t", $lineData[0], $lineData[1], 'mRNA', $lineData[3], $lineData[4], $lineData[5], $lineData[6], $lineData[7], "$mrnaIdAttr;$nameAttr;$parentAttr") );
                                $summaryTxt = join("\t", $lastGeneId, $lineData[0], $lineData[3], $lineData[4], $lineData[6]);
                        }
                        elsif ($lineData[2] eq 'cds')
                        {
                                my $parentAttr = "Parent=$lastmRNAParentId";
                                push(@gff3, join("\t", $lineData[0], $lineData[1], 'CDS', $lineData[3], $lineData[4], $lineData[5], $lineData[6], $lineData[7], "$parentAttr") );

                        }
                }


        }
        if ($isPseudoGene)
        {
                select PSEUDO; print PSEUDO "$fnaSeq\n";
        }
        else
        {
                if ($faaSeq)
                {
                        select FAA; print FAA "$faaSeq\n";
                        my @fastaLine = split /\n/, $faaSeq;
                        if ($fastaLine[1] !~ m/^M/)
                        {
                                $possiblyTruncated{$fastaLine[0]}++;
                        }
                }
                else
                {
                        die "No protein sequence for \n#$record#\n";
                }
                if ($fnaSeq)
                {
                        select FNA; print FNA "$fnaSeq\n";

                }
                else
                {
                        die "No DNA sequence for \n$record\n";
                }
        }
        if (scalar(@gff3))
        {
                my $gff3 = join("\n", @gff3);
                select GFF3; print GFF3 "$gff3\n";
                select SUM; print SUM "$summaryTxt\n";
        }
        else
        {
                die "No DNA sequence for \n$record\n";
        }
}
close(SUM);
close(GFF3);
close(PSEUDO);
close(FNA);
close(FAA);
close(IN);

open(OUT, ">possiblyTruncated.list") or die "Can not open possiblyTruncated.list $!\n";
foreach my $line (keys %possiblyTruncated)
{
        print OUT "$line\n";
}

close(OUT);
exit(0);

sub usage
{
        print "$0 -g=geneWiseInFile\n";
        exit(0);
}


__END__

=head1 DESCRIPTION

Takes genewise records (all concatenated into single file
and parses out the proteins and DNA sequence and reformats
the GFF3 to GBrowse 3 level system

=head1 SYNOPSIS

$0 -g=geneWiseInFile

Where geneWiseInFile is output file from genewise and/or
several genewise results file concatenated into single 
results file.

=head1 DISCLAIMER

THIS SOFTWARE IS PROVIDED "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, 
BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A 
PARTICULAR PURPOSE ARE DISCLAIMED. 

IN NO EVENT SHALL THE AUTHORS, COPYRIGHT HOLDER, OR CONTRIBUTORS BE LIABLE 
FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL 
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR 
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER 
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, 
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE 
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

=head1 GeneWise Record Structure

Genewise records start with:
genewise $Name:  $ (unreleased release)
This program is freely distributed under a GPL. See source directory
Copyright (c) GRL limited: portions of the code are from separate copyright

Query protein:       gene00217-v1.0-hybrid
Comp Matrix:         BLOSUM62.bla
Gap open:            12
Gap extension:       2
Start/End            default
Target Sequence      scaffold00001
Strand:              reverse
Start/End (protein)  default
Gene Parameter file: gene.stat
Splice site model:   GT/AG only
GT/AG bits penalty   -9.95
Codon Table:         codon.table
Subs error:          1e-06
Indel error:         1e-06
Null model           syn
Algorithm            623
Bits   Query         start end Target      start end   idels introns
330.06 gene00217-v1.0-hybrid    1  164 scaffold00001 195114 193933    0    1
//
Gene 1
Gene 195114 193933
  Exon 195114 194878 phase 0
  Exon 194187 193933 phase 0
//
>scaffold00001.[195114:193931].sp.tr
MPAGHGLRSRTRDLFSRAFRKKGYIPLSTYLKTYRIGDHVDVKVNGAVHKGMPHKFYHGR
TGRVWNVTKRAIGVEINKQVGNRIIRKRIHVRVEHVQPSRCTEEFRLRKIKNDQLKAEAK
AKGEVISTKRQPKGPKPGFRVEGAVMETVTPIPYDVVNDLKGGY
//
>scaffold00001.[195114:193931].sp
ATGCCGGCTGGACACGGTCTCCGTTCCCGGACTCGCGATCTGTTCTCGAGGGCCTTCAGG
AAGAAGGGGTACATACCACTCTCCACCTACCTCAAGACCTACAGGATCGGCGACCATGTC
GACGTCAAGGTTAACGGCGCCGTCCACAAGGGTATGCCCCACAAGTTCTACCACGGCCGT
ACCGGTCGCGTCTGGAACGTTACCAAGCGCGCCATCGGTGTCGAGATCAACAAGCAGGTT
GGAAACAGAATCATCAGGAAGAGGATCCACGTGCGTGTGGAGCATGTGCAGCCATCAAGG
TGCACCGAAGAATTCCGTTTGAGGAAGATTAAGAATGACCAGCTGAAGGCCGAGGCTAAA
GCAAAGGGTGAGGTCATTAGCACCAAGAGGCAGCCAAAGGGACCAAAGCCTGGTTTCAGG
GTGGAAGGTGCCGTCATGGAAACTGTTACTCCCATTCCTTATGATGTCGTCAACGATCTT
AAAGGAGGGTAC
//
scaffold00001   GeneWise        match   195114  193933  330.06  -       .       scaffold00001-genewise-prediction-1
scaffold00001   GeneWise        cds     195114  194878  0.00    -       0       scaffold00001-genewise-prediction-1
scaffold00001   GeneWise        intron  194877  194188  0.00    -       .       scaffold00001-genewise-prediction-1
scaffold00001   GeneWise        cds     194187  193933  0.00    -       0       scaffold00001-genewise-prediction-1
//
genewise $Name:  $ (unreleased release)
This program is freely distributed under a GPL. See source directory
Copyright (c) GRL limited: portions of the code are from separate copyright

Query protein:       gene04071-v1.0-hybrid
Comp Matrix:         BLOSUM62.bla
Gap open:            12
Gap extension:       2
Start/End            default
Target Sequence      scaffold00001
Strand:              forward
Start/End (protein)  default
Gene Parameter file: gene.stat
Splice site model:   GT/AG only
GT/AG bits penalty   -9.95
Codon Table:         codon.table
Subs error:          1e-06
Indel error:         1e-06
Null model           syn
Algorithm            623
Bits   Query         start end Target      start end   idels introns
135.90 gene04071-v1.0-hybrid    1  209 scaffold00001 695800 696851    0    2
//
Gene 1
Gene 695800 696851
  Exon 695800 696045 phase 0
  Exon 696166 696245 phase 0
  Exon 696575 696851 phase 2
//
>scaffold00001.[695800:696851].sp.tr
VASAEEAEDTAAPLSSPNSTVSSFQMDFGIRSGGRSSKRDLEVDADRASDDENGSTRKKL
RLSKDQSAFLEESFKEHSTLNPKQKLALAKQLNLRPRQVEVWFQNRRARTKLKQTEVDCE
YLKRCCETLTEENRRLQKELQELRALKTSQPFYMQLPATTLTMCPSCERVVTTASANTST
TTTTTTTNNHIKSVQKKQPRL
//
>scaffold00001.[695800:696851].sp
GTGGCCTCGGCTGAGGAGGCAGAAGATACGGCGGCGCCGCTATCATCTCCGAACAGCACA
GTTTCATCGTTTCAGATGGATTTTGGAATAAGAAGTGGAGGAAGATCGAGTAAAAGAGAT
TTGGAGGTTGATGCCGACAGAGCGAGTGATGACGAGAACGGGTCGACTCGAAAGAAACTC
AGGCTCTCTAAAGATCAATCGGCTTTTCTTGAGGAGAGCTTCAAAGAGCACAGCACTCTC
AATCCTAAGCAAAAACTTGCTCTGGCTAAACAGTTGAATCTTCGTCCTCGCCAAGTGGAA
GTGTGGTTTCAGAATCGAAGAGCAAGGACCAAGCTGAAGCAGACAGAAGTAGATTGCGAG
TACTTAAAGAGATGCTGTGAGACACTGACAGAAGAGAATAGGAGGTTACAAAAAGAACTG
CAAGAATTAAGAGCTTTGAAGACCTCTCAGCCTTTCTACATGCAGCTGCCTGCCACCACA
CTCACCATGTGCCCCTCATGTGAGCGCGTGGTCACAACTGCCTCAGCCAACACCTCCACC
ACCACCACCACCACCACCACCAACAACCACATAAAGTCTGTTCAGAAGAAGCAGCCAAGG
TTA
//
scaffold00001   GeneWise        match   695800  696851  140.81  +       .       scaffold00001-genewise-prediction-1
scaffold00001   GeneWise        cds     695800  696045  0.00    +       0       scaffold00001-genewise-prediction-1
scaffold00001   GeneWise        intron  696046  696165  0.00    +       .       scaffold00001-genewise-prediction-1
scaffold00001   GeneWise        cds     696166  696245  0.00    +       0       scaffold00001-genewise-prediction-1
scaffold00001   GeneWise        intron  696246  696574  0.00    +       .       scaffold00001-genewise-prediction-1
scaffold00001   GeneWise        cds     696575  696851  0.00    +       1       scaffold00001-genewise-prediction-1
//

