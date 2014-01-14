#!/usr/bin/env perl
use warnings;
use strict;

use Pod::Usage;
use Getopt::Long;
use Log::Log4perl qw(:no_extra_logdie_message);
use Log::Log4perl::Level;
my $L = Log::Log4perl::get_logger();
use Data::Dumper;

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

# init a root logger in exec mode
Log::Log4perl->init(\q(
	log4perl.rootLogger					= INFO, Screen
	log4perl.logger.Phrap				= WARN, Screen
	log4perl.appender.Screen			= Log::Log4perl::Appender::Screen
	log4perl.appender.Screen.stderr		= 1
	log4perl.appender.Screen.layout		= PatternLayout
	log4perl.appender.Screen.layout.ConversionPattern = [%d{yy-MM-dd HH:mm:ss}] [prooveval-merge.pl] %m%n
));

GetOptions( 	# NOTE: defaults are set in new();
	\%opt, qw(
		help|h
	)
) or exit(255);

pod2usage(1) if $opt{help};

$L->info('Merging '.@ARGV." files (\n\t".join("\n\t",@ARGV).")");


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
		$o[$i][1] = sprintf("%0.3f", acc(@{$r}[3..$#{$r}]));
	}
}

$L->debug(Dumper(\@i, \@o));

foreach(@o){
	print join("\t", @$_),"\n";
}

sub acc{
	my $bases = 0;
	$bases+=$_ for @_;
	return -1 unless $bases;
	return $_[0] * 100 / $bases;
}

