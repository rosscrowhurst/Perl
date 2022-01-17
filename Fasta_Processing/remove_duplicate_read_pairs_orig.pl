#!/usr/bin/env perl
# original script date: 27/03/2013
# revised 18-09-2013
use strict;
use warnings;
my $read1Infile = "";
my $read2Infile = "";
my $read1Outfile = "";
my $read2Outfile = "";
my $step = 100000;
foreach my $arg(@ARGV)
{
    ($arg =~ m/^-h$/) and usage();
    ($arg =~ m/^-f1=(.+)$/) and $read1Infile = $1;
    ($arg =~ m/^-f2=(.+)$/) and $read2Infile = $1;
    ($arg =~ m/^-o1=(.+)$/) and $read1Outfile = $1;
    ($arg =~ m/^-o2=(.+)$/) and $read2Outfile = $1;
    ($arg =~ m/^-step=(\d+)$/) and $step = $1;
}
($read1Infile eq "") and usage();
($read1Outfile eq "") and usage();
my $pair = 0;
if ($read2Infile ne "")
{
        ($read2Outfile) or usage(error => "Must supply outfile name for read2");
        $pair = 1;
}

my %uniqueReads = ();
my $uniques = 0;
my $dups = 0;
my $total = 0;
open(R1_IN, "<$read1Infile") or die "Can not open read 1 input file: $read1Infile $!\n";
open(R1_OUT, ">$read1Outfile") or die "Can not open read 1 out file: $read1Outfile $!\n";
if ($pair)
{
        open(R2_IN, "<$read2Infile") or die "Can not open read 2 input file: $read2Infile $!\n";
        open(R2_OUT, ">$read2Outfile") or die "Can not open read 2 out file: $read2Outfile $!\n";
}
while (my $read1Name = <R1_IN>)
{
    $total++;
    my $read1Sequence = <R1_IN>;
    my $read1Line3 = <R1_IN>;
    my $read1Qual = <R1_IN>;

    my $read2Name = "";
    my $read2Sequence = "";
    my $read2Line3 = "";
    my $read2Qual = "";

    if ($pair)
    {
        $read2Name = <R2_IN>;
        $read2Sequence = <R2_IN>;
        $read2Line3 = <R2_IN>;
        $read2Qual = <R2_IN>;
    }
    if ($pair)
    {
        my $pairSeq = $read1Sequence . $read2Sequence;

        #if ((!exists $uniqueReads{$read1Sequence}) and (!exists $uniqueReads{$read2Sequence}))
        if (!exists $uniqueReads{$pairSeq})
        {
                $uniqueReads{$pairSeq}++;
                #$uniqueReads{$read1Sequence}++;
                #$uniqueReads{$read2Sequence}++;
                $uniques++;
                select R1_OUT;
                print R1_OUT "$read1Name$read1Sequence$read1Line3$read1Qual";
                select R2_OUT;
                print R2_OUT "$read2Name$read2Sequence$read2Line3$read2Qual";
        }
        else
        {
                 $dups++;
        }

    }
    else
    {
        if (!exists $uniqueReads{$read1Sequence})
        {
                $uniqueReads{$read1Sequence}++;
                $uniques++;
                select R1_OUT;
                print R1_OUT "$read1Name$read1Sequence$read1Line3$read1Qual";
        }
        else
        {
                $dups++;
        }
    }
    if ($total % $step == 0)
    {
        my $percentdup = sprintf("%.2f",  ($dups / $total)*100);
        my $percentUniq = sprintf("%.2f",  ($uniques / $total)*100);
        select STDERR;
        print STDERR "TOTAL: $total     UNIQUE: $uniques ($percentUniq%)        DUPS: $dups ($percentdup%)\n";
    }

}

close(R1_OUT);
close(R1_IN);
if ($pair)
{
        close(R2_OUT);
        close(R2_IN);
}

my $percentdup = sprintf("%.2f", ($dups / $total)*100);
my $percentUniq = sprintf("%.2f", ($uniques / $total)*100);
select STDERR;
print STDERR "TOTAL: $total     UNIQUE: $uniques ($percentUniq%)        DUPS: $dups ($percentdup%)\n";



exit(0);

sub usage
{
    my %params = (error => "", @_);
    select STDERR;
    if (exists $params{'error'})
    {
        print STDERR "Error: $params{'error'}\n";
    }
    print STDERR "USAGE: \n";
    print STDERR "PAIR: $0 -f1=read1_infile -f2=read2_infile -o1=read1_outfile -o2=read2_outfile\n";
    print STDERR "SINGLE: $0 -f1=read1_infile -o1=read1_outfile \n";
    exit(0);
}
