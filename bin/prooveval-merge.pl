#!/usr/bin/env perl
use warnings;
use strict;

use Pod::Usage;
use Getopt::Long;

=head1 NAME

prooveval-merge.pl

=head1 DESCRIPTION

Merge prooveval statistic tables (obtained from threaded prooveval execution).

=head1 CHANGELOG

see git log

=cut
	
=head1 SYNOPSIS

  prooveval-merge.pl STATS1 STATS2 [STATS3 ...] > STATS

=cut

my %opt;

GetOptions( 	# NOTE: defaults are set in new();
	\%opt, qw(
		help|h
	)
) or exit(255);

pod2usage(1) if $opt{help};


my @i;
my @o;
my $file = 0;
$i[$file]=[];
while(<>){
	chomp();
	push @{$i[$file]}, [split("\t", $_)];
	if(eof){
		$file++;
		$i[$file]=[];
	};
}

pop @i; # last one is a stub

$o[0] = $i[0][0]; # header

foreach(@i){
	for(my $i=1; $i<@$_; $i++){
		my $r = $_->[$i];
		$o[$i][0] = $r->[0] || '';
		next unless @$r > 1; # more than 1 column
		for(my $j=2; $j<@$r; $j++){
			my $c = $r->[$j];
			$o[$i][$j]+=$c;
		}
		$o[$i][1] = acc(@{$r}[3..$#{$r}]);
	}
}

foreach(@o){
	print join("\t", @$_),"\n";
}

sub acc{
	my $bases = 0;
	$bases+=$_ for @_;
	return -1 unless $bases;
	return $_[0] * 100 / $bases;
}

