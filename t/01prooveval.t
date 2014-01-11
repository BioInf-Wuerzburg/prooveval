#!/usr/bin/env perl

use strict;
use warnings;

use Test::More;
use Data::Dumper;

use FindBin qw($RealBin);
use lib "$RealBin/../lib/";

use Log::Log4perl qw(:easy :levels);
Log::Log4perl->init(\q(
	log4perl.rootLogger					= DEBUG, Screen
	log4perl.appender.Screen			= Log::Log4perl::Appender::Screen
	log4perl.appender.Screen.stderr		= 1
	log4perl.appender.Screen.layout		= PatternLayout
	log4perl.appender.Screen.layout.ConversionPattern = [%d{yy-MM-dd HH:mm:ss}] [%C] %m%n
));

#--------------------------------------------------------------------------#
=head2 load module

=cut

BEGIN { use_ok('prooveval'); }

my $class = 'prooveval';

#--------------------------------------------------------------------------#
=head2 sample data

=cut


# create data file names from name of this <file>.t
my $Ref_file = $FindBin::RealBin.'/../sample/F.antasticus_genome.fa';
my $Qry_org_file = $FindBin::RealBin.'/../sample/F.antasticus_long_org.fa';
my $Qry_err_file = $FindBin::RealBin.'/../sample/F.antasticus_long_err.fa';
my $Qry_crr_un_file = $FindBin::RealBin.'/../sample/F.antasticus_long_crr.untrimmed.fa';
my $Qry_crr_file = $FindBin::RealBin.'/../sample/F.antasticus_long_crr.trimmed.fa';
my $Sry_crr_file = $FindBin::RealBin.'/../sample/F.antasticus_long_crr.trimmed.fa.sry';

(my $Psl_file = $FindBin::RealScript) =~ s/t$/psl/;
(my $Dmp_file = $FindBin::RealScript) =~ s/t$/dmp/;
(my $Mer_file = $FindBin::RealScript) =~ s/t$/mer/;


#--------------------------------------------------------------------------#
=head2 Modulino _Main

=cut

subtest '_Main' => sub{
	can_ok($class, '_Main');
};

my $self;
my $x;
subtest "$class->new" => sub{
	$self = new_ok($class, [
		ref => [$Ref_file],
		qry => [$Qry_crr_file],
	]);
};


$x = "prep_ref";
subtest '$o->'.$x => sub{
	can_ok($self, $x);
	$self->$x;
};

$x = "prep_qry";
subtest '$o->'.$x => sub{
	can_ok($self, $x);
	$self->$x;
};


# not yet implemented
#$x = "run_gmap";
#subtest '$o->'.$x => sub{
#	can_ok($self, $x);
#	$self->$x;
#};


$x = "parse_gmap";
subtest '$o->'.$x => sub{
	$self->{gmap_sry} = $Sry_crr_file;
	can_ok($self, $x);
	$self->$x;
};

print Dumper($self->{stats});

done_testing();	
__END__

$x = "cluster_ovl";
subtest '$o->'.$x => sub{
	can_ok($self, $x);
	$self->$x;
};

$x = "stats";
subtest '$o->'.$x => sub{
	can_ok($self, $x);
	print $self->$x;
	is($self->$x, "         1 singleton
         2 skipped
         0 failed
         2 success
");
};






















