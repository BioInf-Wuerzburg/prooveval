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
use File::Basename;
use File::Spec;

use Bio::DB::Fasta;
use SeqStore;
use Gmap::Parser;
use Gmap::Record;

use Fasta::Parser;
use Fasta::Seq;

=head1 NAME

prooveval

=head1 DESCRIPTION

Analyse proovread correction results using gmap mappings (and exonerate 
refinements).

=head1 CHANGELOG

see git log

=cut
	
=head1 SYNOPSIS

  prooveval [<OPTIONS>] --ref <genome.fa> --qry <reads.fa> --gmap-sry <reads.sry>
  cat <reads.sry> | prooveval [<OPTIONS>] --ref <genome.fa> --qry <reads.fa> --gmap-sry -

=cut

=head1 OPTIONS

=over

=item -r|--ref=<FASTA>

Reference file, usually genome.

=item -q|--qry=<FASTA>

Query read file.

=item -g|--gmap-sry=<GMAPSRY|->

Gmap mapping results for --qry against --ref in gmap summary format. '-' 
 for STDIN. 

=item [-o|--out=<PATH|->] [`basename --qry`]

Output prefix. '-' for STDOUT.

=item [-t|transcriptomic] [OFF]

Expect spliced mappings in gmap and exonerate.

=item [--[no]-refine-partials] [ON]

Refine partial gmap mappings (coverage < 100%) with exonerate.

=item [--[no]-refine-all] [OFF]

Refine all gmap mappings with exonerate.


=item [--include-edge-mappings] [OFF]

Do not exclude mappings overlapping reference boundaries in statistics.

=item [-T|--tmp-root=<PATH>] [/tmp]

Root location for temporary files.

=item [-k|--keep-tmp] [OFF]

Do not remove temporary files after run.

=item [--gmap-bin=<PATH>] [gmap]

Path to gmap binary.

=item [--exonerate-bin=<PATH>] [exonerate]

Path to exonerate binary.

=back

=cut

##------------------------------------------------------------------------##

=head1 GLOBALS

=cut

our $VERSION = '0.01';

my %opt = (
	refine_partials => 1,
);

# run _Main if modulino is used as executable
__PACKAGE__->_Main(@ARGV) unless caller;


##------------------------------------------------------------------------##

=head1 Main

=cut

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
			ref|r=s@{1,}
			qry|q=s@{1,}
			out|o=s
			gmap_sry|gmap-sry|g=s
			transcriptomic|s!
			refine_partials|refine-partials!
			refine_all|refine-all!
			include_edge_mappings|include-edge-mappings!
			gmap_bin|gmap-bin=s
			exonerate_bin|exonerate-bin=s
			tmp_root|tmp-root|T=s
			keep_tmp|keep-tmp|k!
			debug!
			help|h
		)
	) or exit(255);
	
	pod2usage(1) if $opt{help};
	
	# required	
	for(qw(ref qry gmap_sry)){
		pod2usage("required: --$_") unless defined ($opt{$_}) 
	};
	
	
	$opt{out} = basename($opt{qry}[0]).'.stats' unless $opt{out};
	
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
	# exonerate
		# ref: db, cut using gmap coords
		# qry: SeqStore fetch
		# align
		# parse
	$self->parse_gmap;
	
	$self->stats

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

	my %empty_stats = (
		nr => 0,
		ma => 0,
		mm => 0,
		de => 0,
		in => 0,
		dr => 0,
	);

	# class method -> construct + overwrite
	# init empty obj
	$self = {
		gmap_bin => 'gmap',				# assume exported
		exonerate_bin => 'exonerate', 	# assume exported
		keep_tmp => 0,
		tmp_root => '/tmp',				# set to RAM-disk to increase speed
		ref => [],
		qry => [],
		gmap_sry => undef,
		# overwrite with params
		@_,
		# protected
		tmp => undef,
		tmp_qry => undef,
		tmp_ref => undef,
		tmp_exo => undef,
		ref_store => undef,
		qry_store => undef,
		stats => {
			record_count => 0, 
			gmap => {
				excluded => {
					chimera_count => 0,
					unmapped_count => 0, 
					edge_mapped_count => 0,
					multi_exon_count => 0,
				},
				p0 => {
					1 => {%empty_stats},
					2 => {%empty_stats},
					3 => {%empty_stats},
					4 => {%empty_stats},
					5 => {%empty_stats},
					'1-5' => {%empty_stats},
					'chimera' => {%empty_stats},
				}, # 0 == chimeric
			},
			exo => {
				refined => {%empty_stats},
				bypass => {%empty_stats},
				preref => {%empty_stats},
			},
		},
		_tmp_exo_fh => undef,
	};

	$self->{tmp} = File::Temp->newdir(
			"prooveval-XXXXXXXXXX", 
			DIR => $opt{tmp_root}, 
			CLEANUP => 0
	);
	
	$self->{tmp_qry} = File::Spec->catfile($self->{tmp},'qry.fa');		
	$self->{tmp_ref} = File::Spec->catfile($self->{tmp},'ref.fa');		
	$self->{tmp_exo} = File::Spec->catfile($self->{tmp},'exonerate.aln');		
	
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
	open(EXO, '>', $self->{tmp_exo}) 
		or $L->logdie($!, ": ", $self->{tmp_exo});
	close EXO;
	open($self->{_tmp_exo_fh}, '>>', $self->{tmp_exo}) 
		or $L->logdie($!, ": ", $self->{tmp_exo});
	
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
		#last if $self->{stats}{record_count} >= 25;
		
		$self->{stats}{record_count}++;
		my $pc = @{$r->{paths}};
		unless($pc){
			$S{excluded}{unmapped_count}++;
		}else{
			if ($r->{chimera}){
				$pc = 'chimera';
			};
			
			# analyse primary path: path "0""
			my $p0 = $r->{paths}[0];
			
			## exclude
			# edge mappings
			unless($self->{include_edge_mappings}){
				# check if mapping is to close to an end
				if(
					$p0->{qry_len} > $p0->{ref_hit_to}
					|| 
					$p0->{qry_len} > $p0->{ref_len} - $p0->{ref_hit_from}
				){
					$S{excluded}{edge_mapped_count}++;
					next; 
				}
			}
			# N's (ref/qry)
			# TODO
			
			# multi exon genomic
			unless($opt{transcriptomic}){
				if($p0->{number_of_exons} > 1){
					$S{excluded}{multi_exon_count}++;
					next;
				}
			}
			

			my ($dr, $ma, $mm, $in, $de) = (0,0,0,0,0);
			
			$ma +=$p0->{aln_ma};
			$mm +=$p0->{aln_mm};
			$in +=$p0->{aln_ins};
			$de +=$p0->{aln_del};
			
			
			################################################
			# What shall we do with crappy alignments !!!! #
			################################################
			

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

				if($p1){ # for some reasons there are single exon chimeras
					# multi exon genomic
					unless($opt{transcriptomic}){
						if($p1->{number_of_exons} > 1){
							$S{excluded}{multi_exon_count}++;
							next;
						}
					}
	
					$S{excluded}{chimera_count}++;
					$ma +=$p1->{aln_ma};
					$mm +=$p1->{aln_mm};
					$in +=$p1->{aln_ins};
					$de +=$p1->{aln_del};
		
					$dr = ($p0->{qry_len}-$p0->{qry_hit_len}-$p1->{qry_hit_len});
				}else{
					# single path chimera
				}


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
			
			
			# store gmap stats
			# chimera a stored pc=0 -> individually
			
			$S{p0}{$pc}{nr}++;
			$S{p0}{$pc}{dr} += $dr;
			$S{p0}{$pc}{ma} += $ma;
			$S{p0}{$pc}{mm} += $mm;
			$S{p0}{$pc}{in} += $in;
			$S{p0}{$pc}{de} += $de;
			# run (and store) exonerate stats
			next if $r->{chimera};
			
			if($opt{refine_all} || ($opt{refine_partials} && $dr)){
				# exonerate
				$self->{stats}{exo}{preref}{nr}++;
				$self->{stats}{exo}{preref}{dr} += $dr;
				$self->{stats}{exo}{preref}{ma} += $ma;
				$self->{stats}{exo}{preref}{mm} += $mm;
				$self->{stats}{exo}{preref}{in} += $in;
				$self->{stats}{exo}{preref}{de} += $de;
				
				$self->gmap2seqs($r);
				$self->run_exonerate();
			}else{
				$L->debug('Skipping refinement');
				# use unrefined stats as exo stats
				$self->{stats}{exo}{bypass}{nr}++;
				$self->{stats}{exo}{bypass}{dr} += $dr;
				$self->{stats}{exo}{bypass}{ma} += $ma;
				$self->{stats}{exo}{bypass}{mm} += $mm;
				$self->{stats}{exo}{bypass}{in} += $in;
				$self->{stats}{exo}{bypass}{de} += $de;
			}
			
			
			#last if $r->{id} eq "long_error_7_0";
			
		}
	}
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
	my $tmp_exo_fh = $self->{_tmp_exo_fh};
	#print "\r".$self->{stats}{record_count} unless $self->{stats}{record_count}%100;
	
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
	),
	$opt{transcriptomic} ? qw(--model e2g) : qw(--model a:b),
	qw(
		--exhaustive 1
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
	
	my @tmp_exo = qx($exo_cmd); 
	print $tmp_exo_fh @tmp_exo;
	chomp @tmp_exo;
	
	my @exo_ryo = grep {/^>/}@tmp_exo;
	my @exo_cigar = grep {/^cigar: /}@tmp_exo;

	#print Dumper(\@tmp_exo);
	
	foreach (@exo_ryo){
		s/^>//;
		my ($qid,$qlen,$qalnlen,$et,$em)= split("\t");
		#printf ("$_\t%0.2f\n", $qalnlen/$qlen*100);
		#$qsim/=100;
		$self->{stats}{exo}{refined}{nr}++;
		$self->{stats}{exo}{refined}{ma}+=$et;
		$self->{stats}{exo}{refined}{mm}+=$em;
		$self->{stats}{exo}{refined}{dr}+=($qlen - $qalnlen);
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
		
		$self->{stats}{exo}{refined}{de}+= $cigar{'D'} || 0;
		$self->{stats}{exo}{refined}{in}+= $cigar{'I'} || 0;
	}

	
	
#	exit 1;
}



sub stats{
	my $self = shift;
	
		
	# combined stats
	my $exo_stats = $self->{stats}{exo};
	foreach my $type(qw(refined bypass)){
		foreach (keys %{$exo_stats->{$type}}){
			$exo_stats->{'re+by'}{$_} += $exo_stats->{$type}{$_}
		} 
	}
	
	my $path_stats = $self->{stats}{gmap}{p0};
	foreach my $type(keys %$path_stats){
		next if $type eq 'chimera';
		foreach (keys %{$path_stats->{$type}}){
			$path_stats->{'1-5'}{$_} += $path_stats->{$type}{$_}
		} 
	}
	
	
	$L->debug(Dumper($self->{stats}));
	
	
	my $ofh = \*STDOUT;

	if($opt{out} && $opt{out} ne '-'){
		open($ofh, ">", $opt{out});
	}

	my $pat = "%s:%s\t%0.3f".("\t%d"x6)."\n";


	my $stats = $self->{stats}{gmap}{p0};
	my $type;
	my @cs = qw(nr ma mm de in dr);

	# header
	printf $ofh "%s".("\t%s"x7)."\n", '# source', qw( idy reads match mismatch deletion insertion dropped);
	print $ofh ('-'x100)."\n";

	foreach $type (sort keys %{$stats}){
		next if $type eq 'chimera' or $type eq '1-5';
		printf $ofh $pat, 'gmap:paths', $type, _acc($stats->{$type}) || -1, @{$stats->{$type}}{@cs};
	}
	print $ofh ('-'x100)."\n";
	$type = '1-5';
	printf $ofh $pat, 'gmap:paths', $type,  _acc($stats->{$type}) || -1, @{$stats->{$type}}{@cs};
	
	print $ofh "\n";
	$type = 'chimera';
	printf $ofh $pat, 'gmap', $type,  _acc($stats->{$type}) || -1, @{$stats->{$type}}{@cs};
	print $ofh "\n";
	
	$stats = $self->{stats}{exo};
	foreach $type(qw(bypass preref refined)){
		printf $ofh $pat, 'exo', $type,  _acc($stats->{$type}) || -1, @{$stats->{$type}}{@cs};
	}

	print $ofh ('-'x100)."\n";
	$type = 're+by';
	printf $ofh $pat, 'exo', $type,  _acc($stats->{$type}) || -1, @{$stats->{$type}}{@cs};

	if($opt{out} && $opt{out} ne '-'){
		close $ofh;
	}
}


sub _acc{
	my $s = shift;
	my $bases = $s->{ma}+$s->{mm}+$s->{de}+$s->{in}+$s->{dr};
	return undef unless $bases;
	return $s->{ac} = $s->{ma} * 100 / $bases;
}

1;

=head1 AUTHOR

Thomas Hackl S<thomas.hackl@uni-wuerzburg.de>

=cut



