#!/usr/bin/perl
use strict;
use warnings;
use 5.010;
use Getopt::Long qw(GetOptions);

my $fname='';
my $fraglen = 200;
my $mtemp = 60;
my $primerlen = 20;
my $help;

GetOptions(
    'from|f=s' => \$fname,
    'fraglen|l=s' => \$fraglen,
    'mtemp|t=s' => \$mtemp,
    'primerlen|p=s' => \$primerlen,
    'help' => \$help
    ) or die "Usage: $0 --from NAME\n Please use --help to check the instructions.";
if ($help) {
	print "\-f, \-\-from \n  Name of input file REQUIRED\n\-l, \-\-fraglen \n  Length of the fragment Default : 200\n\-t, \-\-mtemp \n  Melting temperature Default : 60\n\-p, \-\-primerlen \n  Primer length 15 to 25 Default : 20\n\-\-help \n  Help information\n";
	print "Example : perl rsun22project.pl \-f gene_1.fna \-l 200 \-t 56 \-p 18\n"
}
if  ($fname eq '') {
	die "Please use \-\-help to check the instructions!";
}
my $longseq;
open FILE, "$fname";
my $line1 = <FILE>;
while (my $line2 = <FILE>) {
        chomp $line2;
        $longseq = $longseq . $line2;
}
if ($primerlen<15 || $primerlen>25) {
	die "Invalid primer length! Please run the program again\n";
}
my $lennumber;
$lennumber = length($longseq);
my $maxlennumber;
$maxlennumber = $lennumber - 2*$primerlen;
if ($fraglen > $maxlennumber) {
	die "Fragment length is longer than limit $maxlennumber. Please run the program again!\n";
}
if ($fraglen <100 ) {
	die "Fragment length is shorter than limit 100. Please run the program again!\n";
}
my $position1 = 0;
my $length = $lennumber - $fraglen +1;
while ($position1<$length) {
	my $position2 = $position1 + $fraglen -$primerlen;
	my $primer1 = substr ($longseq, $position1, $primerlen);
	my $fragment = substr ($longseq, $position1, $fraglen);
	my $primer2 = substr ($longseq, $position2, $primerlen);
	my @sub1 = split /$primer1/, $longseq;
	my @sub2 = split /$primer2/, $longseq;
	my $sub11 = join (' ',@sub1);
	my $sub22 = join (' ',@sub2);
	my $rprimer2 = $primer2;
	$rprimer2 =~ tr/ACGT/TGCA/;
	if ($primer1 ne $sub11 && $primer2 ne $sub22 && tempdiff($primer1, $mtemp)<5 && tempdiff($rprimer2, $mtemp)<5) {
		my $temp1=melttemp($primer1);
		my $temp2=melttemp($rprimer2);
		print "Forward primer is $primer1 at melting point $temp1\nFragment sequence is $fragment \nReverse primer is $rprimer2 at melting point $temp2\n\n";
	}
	$position1 += 1;
}

sub tempdiff {
	my $primer = $_[0];
	my $g = ($primer =~ tr/G//);
	my $c = ($primer =~ tr/C//);
	my $temp = 64.9 + 41*($g+$c-16.4)/(length($primer));
	my $diff = abs($temp-$_[1]);
	return $diff;
}
sub melttemp {
	my $primer = $_[0];
	my $g = ($primer =~ tr/G//);
	my $c = ($primer =~ tr/C//);
	my $temp = 64.9 + 41*($g+$c-16.4)/(length($primer));
}
close FILE;

