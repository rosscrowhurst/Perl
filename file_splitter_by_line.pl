#!/usr/bin/perl -w
use strict;
use FileHandle;

my $usage = "USAGE: $0 file2split number_of_parts\n";
(-r $ARGV[0]) or die $usage;
($ARGV[1] =~ m/^\d+$/) or die $usage;

my $filehandles = {};
foreach my $key (0 .. ($ARGV[1] - 1))
{
        open( my $fh, '>', "myfile_$key") or die "Cannot open myfile_$key for writing $!\n";
        $filehandles->{$key} = $fh;
}
open(IN, "<$ARGV[0]") or die "can not open $ARGV[0] $!\n";
my $inputLineCount = 0;
while (my $line = <IN>)
{
        my $idx = $inputLineCount % $ARGV[1];
        print { $filehandles->{$idx} } $line;
        $inputLineCount++;
}
close(IN);

foreach my $key (0 .. ($ARGV[1] - 1))
{
        close  { $filehandles->{$key} };
}
