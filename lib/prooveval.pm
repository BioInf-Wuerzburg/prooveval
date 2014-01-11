#!/usr/bin/perl
package prooveval;

use warnings;
no warnings 'qw';
use strict;


use Carp;
use Getopt::Long;
use Pod::Usage;
use Log::Log4perl qw(:no_extra_logdie_message);
use Log::Log4perl::Level;
my $L = Log::Log4perl::get_logger();

use Data::Dumper;
$Data::Dumper::Sortkeys = 1;

use File::Path;

use Bio::DB::Fasta;
use SeqStore;
use Gmap::Parser;
use Gmap::Record;

use Fasta::Parser;
use Fasta::Seq;

=head1 Changelog

see git log

=cut

##------------------------------------------------------------------------##

our $VERSION = '0.01';

my %opt = (
	refine_partials => 1
);

# run _Main if modulino is used as executable
__PACKAGE__->_Main(@ARGV) unless caller;


##------------------------------------------------------------------------##

sub _Main{
	my $class = shift;
	
	# init a root logger in exec mode
	Log::Log4perl->init(\q(
		log4perl.rootLogger					= INFO, Screen
		log4perl.logger.Phrap				= WARN, Screen
		log4perl.appender.Screen			= Log::Log4perl::Appender::Screen
		log4perl.appender.Screen.stderr		= 1
		log4perl.appender.Screen.layout		= PatternLayout
		log4perl.appender.Screen.layout.ConversionPattern = [%d{yy-MM-dd HH:mm:ss}] [%C] %m%n
	));
	
	
	GetOptions( 	# NOTE: defaults are set in new();
		\%opt, qw(
			ref=s@{1,}
			qry=s@{1,}
			gmap_sry|gmap-sry=s
			tmp_dir=s
			gmap_bin|gmap-bin=s
			exonerate_bin|exonerate-bin=s
			allow_edge_mappings|allow-edge-mappings!
			refine_partials|refine-partials!
			refine_all|refine-all!
			keep_tmp|keep-tmp!
			debug!
		)
	) or exit(255);
	
	$L->level($DEBUG) if $opt{debug};
	$L->debug('Verbose level set to DEBUG');
	$L->debug(Dumper(\%opt));	
	my $self = $class->new(%opt);

	# check input
	#	qry, ref, gmap-sry
	$self->prep_ref();
	
	$self->prep_qry();
	
	
	## run
		
	# run gmap (unless gmap-sry)
		# ref: db
		# mapping
	unless ($self->{gmap_sry}){
		$self->run_gmap();
	}
	
	# parse gmap-sry
	$self->parse_gmap;
	
	# exonerate
		# ref: db, cut using gmap coords
		# qry: SeqStore fetch
		# align
		# parse
	
	
	print Dumper($self->{stats});
	
	my $exo = $self->{stats}{exo};
	printf "%0.5f\n", $exo->{ma} * 100 / ($exo->{ma}+$exo->{mm}+$exo->{de}+$exo->{in}+$exo->{dr});
}








##------------------------------------------------------------------------##

=head2 new

Instantiate prooveval modulino handle

=cut

sub new{
	my $proto = shift;
	my $self;
	my $class;
	
	# object method -> clone + overwrite
	if($class = ref $proto){ 
		return bless ({%$proto, @_}, $class);
	}

	# class method -> construct + overwrite
	# init empty obj
	$self = {
		gmap_bin => 'gmap',				# assume exported
		exonerate_bin => 'exonerate', 	# assume exported
		keep_tmp => 0,
		tmp_dir => '/tmp', 					# set to RAM-disk to increase speed
		ref => [],
		qry => [],
		gmap_sry => undef,
		exo_sry => 'prooveval.exonerate.aln',
		# overwrite with params
		@_,
		# protected
		tmp => undef,
		tmp_qry => undef,
		tmp_ref => undef,
		ref_store => undef,
		qry_store => undef,
		stats => {
			record_count => 0, 
			gmap => {
				unmapped_count => 0, 
				chimera_count => 0,
				edge_mapped_count => 0,
				p0 => [], # 0 == chimeric
			},
			exo => {
				
			},
		},
		_exo_sry_fh => undef,
	};
	
	$self->{tmp} = $self->{tmp_dir}.'/prooveval/';
	$self->{tmp_qry} = $self->{tmp}.'/qry.fa';
	$self->{tmp_ref} = $self->{tmp}.'/ref.fa';
	
	
	bless $self, $proto;
	$L->debug(ref $self, ' instantiated');
	
	# create /tmp/dir
	$L->debug('Setting up tmp file dir "', $self->{tmp}, '"');
	unless(
		File::Path::make_path($self->{tmp})
		||
		-d $self->{tmp}
	){
		$L->logdie('Cannot create tmp data in '.$self->{tmp})
	}
	
	#
	open(EXO, '>', $self->{exo_sry}) 
		or $L->logdie($!, ": ", $self->{exo_sry});
	close EXO;
	open($self->{_exo_sry_fh}, '>>', $self->{exo_sry}) 
		or $L->logdie($!, ": ", $self->{exo_sry});
	
	# return
	return $self;
}


sub DESTROY{
	my $self = shift;
	unless($self->{keep_tmp}){
		File::Path::remove_tree($self->{tmp}, {result => \(my $rm)} );
		$L->debug('Removed tmp data: ', @$rm); 
	};
}

=head2 prep_ref

Prepare reference data (Index, db)

=cut

sub prep_ref{
	my ($self) = @_;
	$L->info("Preparing reference db");
	$self->{ref_store} = SeqStore->new(path => $self->{ref}[0]);
}



=head2 prep_qry

Prepare query data (Index, db)

=cut

sub prep_qry{
	my ($self) = @_;
	$L->info("Preparing query db");
	$self->{qry_store} = SeqStore->new(path => $self->{qry}[0]);
}


=head2 run_gmap

=cut

sub run_gmap{
	$L->logdie('Running gmap within pipeline not yet implemented!');
}


=head2 parse_gmap

=cut

sub parse_gmap{
	my ($self) = @_;

	$L->info('Parsing gmap summary and running exonerate refinements');

	$L->logdie("Gmap summary file required") unless $self->{gmap_sry} && (-e $self->{gmap_sry} || $self->{gmap_sry} eq '-');
	
	my $gmap_sry = $self->{gmap_sry} eq '-' ? undef : $self->{gmap_sry};
	
	my $gp = Gmap::Parser->new(file => $gmap_sry)->check_format() 
	|| $L->exit("$self->{gmap_sry} does not look like gmap summary"); # defaults to &STDIN
	
	my %S = %{$self->{stats}{gmap}};
	$self->{stats}{gmap} = \%S;
	
	while(my $r = $gp->next_record){
		$self->{stats}{record_count}++;
		my $pc = @{$r->{paths}};
		unless($pc){
			$S{unmapped_count}++;
		}else{
			if ($r->{chimera}){
				$pc = 0;
			};
			
			# analyse primary path: path "0""
			my $p0 = $r->{paths}[0];
			
			## exclude
			# edge mappings
			unless($self->{allow_edge_mappings}){
				# check if mapping is to close to an end
				if(
					$p0->{qry_len} > $p0->{ref_hit_to}
					|| 
					$p0->{qry_len} > $p0->{ref_len} - $p0->{ref_hit_from}
				){
					$S{edge_mapped_count}++;
					next; 
				}
			}
			# N's (ref/qry)

			my ($dr, $ma, $mm, $in, $de) = (0,0,0,0,0);
			
			$ma +=$p0->{aln_ma};
			$mm +=$p0->{aln_mm};
			$in +=$p0->{aln_ins};
			$de +=$p0->{aln_del};

=pod	
	#		if($opt{ref_cov}){
	#			# init ref cov
	#			unless(exists $REF_COV{$p0->{ref_id}}){
	#				$REF_COV{$p0->{ref_id}} = [(0)x $p0->{ref_len}];
	#			}
	#			map{$_++}@{$REF_COV{$p0->{ref_id}}}[$p0->{ref_hit_from}..$p0->{ref_hit_to}];
	#		}
=cut

			# include secondary path in analysis if chimera
			if($r->{chimera}){
				my $p1 = $r->{paths}[1];
				$S{chimera_count}++;
				$ma +=$p1->{aln_ma};
				$mm +=$p1->{aln_mm};
				$in +=$p1->{aln_ins};
				$de +=$p1->{aln_del};
	
				$dr = ($p0->{qry_len}-$p0->{qry_hit_len}-$p1->{qry_hit_len});

=pod				
	#			if($opt{ref_cov}){
	#				# init ref cov
	#				unless(exists $REF_COV{$p1->{ref_id}}){
	#					$REF_COV{$p1->{ref_id}} = [(0)x $p1->{ref_len}];
	#				}
	#				map{$_++}@{$REF_COV{$p1->{ref_id}}}[$p1->{ref_hit_from}..$p1->{ref_hit_to}];
	#			}
=cut

			}else{
				$dr = ($p0->{qry_len}-$p0->{qry_hit_len});
			}
			
			$S{p0}[$pc]{dr} += $dr;
			$S{p0}[$pc]{ma} += $ma;
			$S{p0}[$pc]{mm} += $mm;
			$S{p0}[$pc]{in} += $in;
			$S{p0}[$pc]{de} += $de;
	
			
			if($opt{refine_all} || ($opt{refine_partials} && $dr)){
				# exonerate
				$self->gmap2seqs($r);
				$self->run_exonerate();
			}else{
				$L->debug('Skipping refinement');
				# use unrefined stats like exo stats
				$self->{stats}{exo}{dr} += $dr;
				$self->{stats}{exo}{ma} += $ma;
				$self->{stats}{exo}{mm} += $mm;
				$self->{stats}{exo}{in} += $in;
				$self->{stats}{exo}{de} += $de;
			}
			
			
			#last if $r->{id} eq "long_error_7_0";
			
		}
	}
	print "\n";
	
}


=head2 gmap2seqs

=cut

sub gmap2seqs{
	my ($self, $r) = @_;
	my $p0 = $r->{paths}[0];
	my ($from, $to) = sort{$a <=> $b}($p0->{ref_hit_from}, $p0->{ref_hit_to});
	$from-= $p0->{qry_len}*1.5;
	$to+= $p0->{qry_len}*1.5;
	my $qry = $self->{qry_store}->get_seq($r->id);
	my $ref = $self->{ref_store}->get_seq($p0->{ref_id}, $from, $to);
	
	open (QRY, '>', $self->{tmp_qry}) or $L->logdie($!,": ", $self->{tmp_qry});
	open (REF, '>', $self->{tmp_ref}) or $L->logdie($!,": ", $self->{tmp_ref});

	print QRY $qry;
	print REF $ref;
	
	close QRY;
	close REF;
	


}

=head2 run_exonerate

Run exonrate for a mapped read to refine alignment results

=cut

sub run_exonerate{
	my ($self) = @_;
	my $exo_sry_fh = $self->{_exo_sry_fh};
	print "\r".$self->{stats}{record_count} unless $self->{stats}{record_count}%100;
	
	#last if $self->{stats}{record_count} > 25;

=pod
#		--frameshift -5
#		--model est2genome


realignment of GENOMIC reads using exonerate
	--querytype dna
	--targettype dna

	--model affine:bestfit  semi-global, query cov 100% ensured
	--exhaustive			required for affine:bestfit model, which isnt heuristic

	--frameshift -5			reduced from -28 : gaps are not frame depended 
	--gapopen -5 			reduced	from -12 : gaps are quite likely 
	--gapextend -3			reduced from -4
	--codongapopen -5		reduced from -18 : gaps are not codon depended 
	--codongapextend -3		reduced from -8

	--bestn 1				we expect a single result
	--subopt n				no interest in subopt

	--showalignment 1 
	--alignmentwidth 150 
	--showvulgar 0
	--showcigar 1	
	--ryo ">%ti\t%tl\t%tab\t%tae\t%tal\t%qi\t%ql\t%qab\t%qae\t%qal\t%qS\t%ps\t%et\t%ei\n"


realignment of TRANSCRIPTOMIC reads using exonerate
	--querytype dna
	--targettype dna

	--model est2genome		intron modeling, query terminal base drops possible => compensated somehow with low gap costs

	--frameshift -5			reduced from -28 : gaps are not frame depended 
	--gapopen -5 			reduced	from -12 : gaps are quite likely
	--gapextend -3			reduced from -4
	--codongapopen -5		reduced from -18 : gaps are not codon depended 
	--codongapextend -3		reduced from -8
	
	--bestn 1				we expect a single result
	--subopt n				no interest in subopt

=cut

	my $exo_cmd=join(" ", qw(
		exonerate 
		--querytype dna
		--targettype dna

		--exhaustive 1
		--model affine:bestfit
		--bestn 1
		--subopt n
		--frameshift -5

		--showalignment 1 
		--alignmentwidth 150 
		--showvulgar 0
		--showcigar 1
	),
	("--ryo" => '">'.join('\t', 
		"%qi",		# query id
		"%ql",		# query length
#		"%qab",		# query alignment start
#		"%qae",		# query alignment end
		"%qal",		# query alignment length
#		"%qS",		# query strand
#		"%ps", 		# query similarity
#		"%et",		# equivalenced total (m+mm)
		"%ei",		# equivalenced matches
		"%em",		# equivalenced mismatches
	).'\n"'),
	, $self->{tmp_qry}, $self->{tmp_ref});
	
	$L->debug($exo_cmd);
	
	my @exo_sry = qx($exo_cmd); 
	print $exo_sry_fh @exo_sry;
	chomp @exo_sry;
	
	my @exo_ryo = grep {/^>/}@exo_sry;
	my @exo_cigar = grep {/^cigar: /}@exo_sry;

	#print Dumper(\@exo_aln);
	
	foreach (@exo_ryo){
		s/^>//;
		my ($qid,$qlen,$qalnlen,$et,$em)= split("\t");
		#printf ("$_\t%0.2f\n", $qalnlen/$qlen*100);
		#$qsim/=100;
		$self->{stats}{exo}{ma}+=$et;
		$self->{stats}{exo}{mm}+=$em;
		$self->{stats}{exo}{dr}+=($qlen - $qalnlen);
	}
	
	
	foreach (@exo_cigar){
		my ($cigar, $query_id, $query_start, $query_end, $query_strand,
		$target_id, $target_start, $target_end, $target_strand, $score,
		@cigar) = split(/\s+/, $_);

		$L->debug(@cigar);

		my %cigar;
		for(my $i=0; $i<@cigar; $i+=2 ){
			$cigar{$cigar[$i]}+= $cigar[$i+1];
		}
		
		my @unknown_cigar = grep{/[^DIM]/}keys %cigar;
		$L->logdie("Unknown cigars: @unknown_cigar") if @unknown_cigar; 
		
		$self->{stats}{exo}{de}+= $cigar{'D'} || 0;
		$self->{stats}{exo}{in}+= $cigar{'I'} || 0;
	}

	
	
#	exit 1;
}






1;

=head1 AUTHOR

Thomas Hackl S<thomas.hackl@uni-wuerzburg.de>

=cut


