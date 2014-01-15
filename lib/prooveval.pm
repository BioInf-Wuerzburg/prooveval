#!/usr/bin/perl
package prooveval;

use warnings;
no warnings 'qw';
use strict;


use Carp;
use Getopt::Long qw(:config no_ignore_case bundling);
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

use threads;
use Thread::Queue;

=head1 NAME

prooveval

=head1 DESCRIPTION

Analyse proovread correction results using gmap mappings (and exonerate 
refinements).

=cut

=head1 CHANGELOG

see git log

=cut

=head1 TODO

=over

=item N50

=item TP

=item Unmapped

=back

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

=item [-t|threads]

Number of threads for running exonerate.

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

=item [-M|--max-records=<INT>] [0]

Stop evaluation after this many records. "0" reads all records.

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
			ref|r=s
			qry|q=s
			out|o=s
			gmap_sry|gmap-sry|g=s
			transcriptomic|s!
			thread_num|threads|t=i
			refine_partials|refine-partials!
			refine_all|refine-all!
			include_edge_mappings|include-edge-mappings!
			gmap_bin|gmap-bin=s
			exonerate_bin|exonerate-bin=s
			tmp_root|tmp-root|T=s
			keep_tmp|keep-tmp|k!
			max_records|max-records|M=i
			debug!
			help|h
		)
	) or exit(255);
	
	pod2usage(1) if $opt{help};
	
	# required	
	for(qw(ref qry gmap_sry)){
		pod2usage("required: --$_") unless defined ($opt{$_}) 
	};
	
	
	$opt{out} = basename($opt{qry}).'.stats' unless $opt{out};
	
	$L->level($DEBUG) if $opt{debug};
	$L->debug('Verbose level set to DEBUG');
	$L->debug(Dumper(\%opt));	
	my $self = $class->new(%opt);

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
	
	$self->stats;
	
	unless($self->{keep_tmp}){
		$L->debug("Removing tmp dir");
		File::Path::remove_tree($self->{tmp}, {result => \(my $rm)} );
		$L->debug("Removed tmp data: @$rm"); 
	};
	

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
		tmp_root => '/tmp',				# set to RAM-disk to increase speed
		ref => [],
		qry => [],
		gmap_sry => undef,
		thread_num => 1,
		# overwrite with params
		@_,
		# protected
		tmp => undef,
		ref_store => undef,
		qry_store => undef,
		record_count => 0, 
		stats => {
			gmap_unmapped => empty_stats(),
			gmap_chimera => empty_stats(),
			gmap_edge_mapped => empty_stats(),
			gmap_multi_exon => empty_stats(),
			'gmap_p0:1' => empty_stats(),
			'gmap_p0:2' => empty_stats(),
			'gmap_p0:3' => empty_stats(),
			'gmap_p0:4' => empty_stats(),
			'gmap_p0:5' => empty_stats(),
			exo_refined => empty_stats(),
			exo_bypass => empty_stats(),
			exo_preref => empty_stats(),
		},
		_thread_queue_length => 50,
	};


	bless $self, $proto;

	# tmp dir
	$L->debug('Setting up tmp...');

	$self->{tmp} = File::Temp->newdir(
			"prooveval-XXXXXXXXXX", 
			DIR => $self->{tmp_root}, 
			CLEANUP => 0
	);

	unless(
		File::Path::make_path($self->{tmp})
		||
		-d $self->{tmp}
	){
		$L->logdie('Cannot create tmp data in '.$self->{tmp})
	}

	$L->debug("tmp dir: ", $self->{tmp});

	# prepare thread
	$self->{tmp_qry_patt} = File::Spec->catfile($self->{tmp},"qry.%d.fa");
	$self->{tmp_ref_patt} = File::Spec->catfile($self->{tmp},"ref.%d.fa");
	$self->{tmp_exo_patt} = File::Spec->catfile($self->{tmp},"exonerate#%d.fa");
	
	
	# exonerate parameter
	$self->{exo_opt} = join(" ", qw(
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
	);
	
	$L->debug("Exonerate options: ", $self->{exo_opt});
	
	# thread launching
	my $Q = Thread::Queue->new();
	$L->debug("Initialized thread queue");

	$L->debug("Launching $self->{thread_num} threads ...");
	for (1..$self->{thread_num}){
		threads->create({context => 'list'}, sub {
		
			my $tmp_exo = File::Spec->catfile($self->{tmp},'exonerate#'.threads->tid().'.aln');		
			
			$L->debug("Launched thread#",threads->tid());
			
			# loop
			while (my $c = $Q->dequeue(1)){
				$self->run_exonerate($c);
			}
			
			# return
			return %{$self->{stats}{exo_refined}}; # cant return $self ref, gets destroyed
		})
	};

	sleep 1;

	$self->{_thread_queue} = $Q;

	$L->debug("Launched $self->{thread_num} threads");


	# check input, create/load dbs
	$L->debug('Setting up ref and qry db...');

	$self->prep_ref();
	$self->prep_qry();

	$L->debug('Loaded/created ref and qry db...');

	
	# return
	$L->debug(ref $self, ' object fully invoked');
	return $self;
}

#''
#sub DESTROY{
#	my $self = shift;
#	
#	# thread safe
#	unless(threads->tid()){
#		# needs to be thread safe - only main thread is allowed to remove
#		unless($self->{keep_tmp}){
#			$L->debug("Removing tmp dir");
#			File::Path::remove_tree($self->{tmp}, {result => \(my $rm)} );
#			$L->debug("Removed tmp data: @$rm"); 
#		};
#	}
#}



=head2 prep_ref

Prepare reference data (Index, db)

=cut

sub prep_ref{
	my ($self) = @_;
	$L->info("Preparing reference db");
	$self->{ref_store} = SeqStore->new(path => $self->{ref});
}



=head2 prep_qry

Prepare query data (Index, db)

=cut

sub prep_qry{
	my ($self) = @_;
	$L->info("Preparing query db");
	$self->{qry_store} = SeqStore->new(path => $self->{qry});
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

	$L->debug('Setting up gmap summary parser');
	
	$L->logdie("Gmap summary file required") unless $self->{gmap_sry} && (-e $self->{gmap_sry} || $self->{gmap_sry} eq '-');
	
	my $gmap_sry = $self->{gmap_sry} eq '-' ? undef : $self->{gmap_sry};
	
	my $gp = Gmap::Parser->new(file => $gmap_sry)->check_format() 
	|| $L->exit("$self->{gmap_sry} does not look like gmap summary"); # defaults to &STDIN
	
	$L->debug("Instanciated gmap parser");
	
	my %S = %{$self->{stats}};
	$self->{stats} = \%S;
	
	while(1){
		next unless $self->{_thread_queue}->pending < $self->{_thread_queue_length};

		$L->debug("Processing record");
		
		my $r = $gp->next_record;
		last unless defined $r;
		
		last if $self->{max_records} && $self->{record_count} >= $self->{max_records};
		
		$self->{record_count}++;
		my $pc = @{$r->{paths}};
		unless($pc){
			$L->debug("Excluded: unmapped");
			$S{gmap_unmapped}{nr}++;
			# TODO: from db!!
			push @{$S{gmap_unmapped}{ql}}, length($self->{qry_store}->get_seq($r->id)->seq);
			
		}else{
			if ($r->{chimera}){
				# we do stats, but store them separately, no refinement
				$pc = 'gmap_chimera'; 
			}else{
				$pc = 'gmap_p0:'.$pc;
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
					$L->debug("Excluded: edge_mapped");
					$S{gmap_edge_mapped}{nr}++;
					push @{$S{gmap_edge_mapped}{ql}}, $p0->{qry_len};
					next; 
				}
			}
			# N's (ref/qry)
			# TODO
			
			# multi exon genomic
			unless($opt{transcriptomic}){
				if($p0->{number_of_exons} > 1){
					$L->debug("Excluded: multi_exon");
					$S{gmap_multi_exon}{nr}++;
					push @{$S{gmap_multi_exon}{ql}}, $p0->{qry_len};
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
							$L->debug("Excluded: multi_exon (chimera secondary path)");
							$S{gmap_multi_exon}{nr}++;
							push @{$S{gmap_multi_exon}{ql}}, $p0->{qry_len};
							next;
						}
					}
	
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
			
			$S{$pc}{nr}++;
			$S{$pc}{bp}{dr} += $dr;
			$S{$pc}{bp}{ma} += $ma;
			$S{$pc}{bp}{mm} += $mm;
			$S{$pc}{bp}{in} += $in;
			$S{$pc}{bp}{de} += $de;
			push @{$S{$pc}{ql}}, $p0->{qry_len};
			
			# run (and store) exonerate stats
			if($r->{chimera}){
				$L->debug("Excluded: chimera");
				next;
			};
			
			
			if($opt{refine_all} || ($opt{refine_partials} && $dr)){
				# exonerate
				$L->debug("Enqueuing record for refinement");
				$self->{stats}{exo_preref}{nr}++;
				$self->{stats}{exo_preref}{bp}{dr} += $dr;
				$self->{stats}{exo_preref}{bp}{ma} += $ma;
				$self->{stats}{exo_preref}{bp}{mm} += $mm;
				$self->{stats}{exo_preref}{bp}{in} += $in;
				$self->{stats}{exo_preref}{bp}{de} += $de;
				push @{$S{exo_preref}{ql}}, $p0->{qry_len};
				
				$self->gmap2seqs($r, $self->{record_count});
				$self->{_thread_queue}->enqueue($self->{record_count});
				
			}else{
				# use unrefined stats as exo stats
				$L->debug('Bypassing refinement');
				$self->{stats}{exo_bypass}{nr}++;
				$self->{stats}{exo_bypass}{bp}{dr} += $dr;
				$self->{stats}{exo_bypass}{bp}{ma} += $ma;
				$self->{stats}{exo_bypass}{bp}{mm} += $mm;
				$self->{stats}{exo_bypass}{bp}{in} += $in;
				$self->{stats}{exo_bypass}{bp}{de} += $de;
				push @{$S{exo_bypass}{ql}}, $p0->{qry_len}
			}
			#last if $r->{id} eq "long_error_7_0";
		}
	}
	

	$L->debug("Joining threads");
	# send undef to worker threads as "done" signal
	$self->{_thread_queue}->enqueue((undef)x $self->{thread_num});


	# join threads and integrate refined exonerate stats
	foreach (threads->list()){
		$L->debug("Joining thread#".$_->tid());
		append_stats($self->{stats}{exo_refined}, {$_->join()});

	}
	
	$L->debug('Done parsing and refining');
}


=head2 gmap2seqs

=cut

sub gmap2seqs{
	my ($self, $r, $c) = @_;
	my $p0 = $r->{paths}[0];
	my ($from, $to) = sort{$a <=> $b}($p0->{ref_hit_from}, $p0->{ref_hit_to});
	$from-= $p0->{qry_len}*1.5;
	$to+= $p0->{qry_len}*1.5;
	my $qry = $self->{qry_store}->get_seq($r->id);
	my $ref = $self->{ref_store}->get_seq($p0->{ref_id}, $from, $to);
	
	open (QRY, '>', $self->tmp_qry($c)) or $L->logdie($!,": ", $self->tmp_qry($c));
	open (REF, '>', $self->tmp_ref($c)) or $L->logdie($!,": ", $self->tmp_ref($c));

	print QRY $qry;
	print REF $ref;
	
	close QRY;
	close REF;
	


}

=head2 run_exonerate

Run exonrate for a mapped read to refine alignment results

=cut

sub run_exonerate{
	my ($self, $c) = @_;
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

	
	my $tmp_ref = $self->tmp_ref($c);
	my $tmp_qry = $self->tmp_qry($c);
	my $tmp_exo = $self->tmp_exo(threads->tid());
	
	my @tmp_exo = qx($self->{exonerate_bin} $self->{exo_opt} $tmp_qry $tmp_ref | tee -a $tmp_exo); # TODO: use tee tmp_file  
	chomp @tmp_exo;
	
	my @exo_ryo = grep {/^>/}@tmp_exo;
	my @exo_cigar = grep {/^cigar: /}@tmp_exo;

	#print Dumper(\@tmp_exo);
	
	foreach (@exo_ryo){
		s/^>//;
		my ($qid,$qlen,$qalnlen,$et,$em)= split("\t");
		#printf ("$_\t%0.2f\n", $qalnlen/$qlen*100);
		#$qsim/=100;
		$self->{stats}{exo_refined}{nr}++;
		$self->{stats}{exo_refined}{bp}{ma}+=$et;
		$self->{stats}{exo_refined}{bp}{mm}+=$em;
		$self->{stats}{exo_refined}{bp}{dr}+=($qlen - $qalnlen);
		
	}
	
	
	foreach (@exo_cigar){
		my ($cigar, $query_id, $query_start, $query_end, $query_strand,
		$target_id, $target_start, $target_end, $target_strand, $score,
		@cigar) = split(/\s+/, $_);

		my %cigar;
		for(my $i=0; $i<@cigar; $i+=2 ){
			$cigar{$cigar[$i]}+= $cigar[$i+1];
		}
		
		my @unknown_cigar = grep{/[^DIM]/}keys %cigar;
		$L->logdie("Unknown cigars: @unknown_cigar") if @unknown_cigar; 
		
		$self->{stats}{exo_refined}{bp}{de}+= $cigar{'D'} || 0;
		$self->{stats}{exo_refined}{bp}{in}+= $cigar{'I'} || 0;
	}

	unlink $tmp_ref, $tmp_qry unless $self->{keep_tmp};
	
#	exit 1;
}



sub stats{
	my $self = shift;
	
	$L->info("compiling statistics");
	
	my $S = $self->{stats};
	# combined stats
	
	$S->{exo_refined}{ql} = $S->{exo_preref}{ql}; 
	
	# exo_re+by
	$S->{'exo_re+by'} = empty_stats();
	append_stats($S->{'exo_re+by'}, @{$S}{qw(exo_refined exo_bypass)});

	# gmap_p0:1-5
	$S->{'gmap_p0:1-5'} = empty_stats();
	append_stats($S->{'gmap_p0:1-5'}, @{$S}{qw(gmap_p0:1 gmap_p0:2 gmap_p0:3 gmap_p0:4 gmap_p0:5)});
	
	foreach (values %$S){
		if($_->{ql}){
			my %nx = Nx(lengths => $_->{ql}, N => [50]);
			$_->{N50} = $nx{50};
		}
	}
	
	$L->debug(Dumper($S));
	
	
	my $ofh = \*STDOUT;

	if($opt{out} && $opt{out} ne '-'){
		open($ofh, ">", $opt{out});
	}

	print $ofh join("\t", qw(category R_total R_category R_N50 B_idy B_match B_mismatch B_deletion B_insertion B_dropped)),"\n";
	print $ofh ('-'x100)."\n";

	foreach my $k(sort keys %$S){
		print join("\t", 
			substr($k, 0, 14), 
			$self->{record_count}, 
			$S->{$k}{nr}, 
			$S->{$k}{N50} || '-NA-', 
			@{$S->{$k}{bp}}{qw(ma mm de in dr)}
		), "\n";
	}

	if($opt{out} && $opt{out} ne '-'){
		close $ofh;
	}
}


sub append_stats{
	my $s = shift;
	
	foreach my $a (@_){
		foreach (qw(ma mm de in dr)){
			$s->{bp}{$_} += $a->{bp}{$_};
		}
		$s->{nr} += $a->{nr};
		$s->{ql} = [] unless defined $s->{ql};
		push @{$s->{ql}}, @{$a->{ql}};
	}
}

sub _acc{
	my $s = shift;
	my $bases = $s->{bp}{ma}+$s->{bp}{mm}+$s->{bp}{de}+$s->{bp}{in}+$s->{bp}{dr};
	return undef unless $bases;
	return $s->{ac} = $s->{bp}{ma} * 100 / $bases;
}

sub tmp_qry{
	my ($self, $c) = @_;
	return sprintf $self->{tmp_qry_patt}, $c;
}

sub tmp_ref{
	my ($self, $c) = @_;
	return sprintf $self->{tmp_ref_patt}, $c;
}

sub tmp_exo{
	my ($self, $c) = @_;
	return sprintf $self->{tmp_exo_patt}, $c;
}

sub Nx{
	my $p = {
		lengths => [],
		N => [50],
		@_
	};
	
	my $total_length = 0;
	$total_length+=$_ for @{$p->{lengths}};
	my %Ns;
	my @L = sort{$b <=> $a}@{$p->{lengths}};
	my @ls = map{$total_length * ($_/100)}@{$p->{N}};
	
	my $lc = 0;
	foreach my $l (@L){
		$lc+= $l;
		if($lc >= $ls[0]){
			shift @ls;
			$Ns{shift @{$p->{N}}} = $l;
			last unless @ls;
		}
	}
	return %Ns;
}	


sub empty_stats{
	return {
		ql => [],
		nr => 0,
		bp => {
			ma => 0,
			mm => 0,
			de => 0,
			in => 0,
			dr => 0,
	 	}
	};
}

1;

=head1 AUTHOR

Thomas Hackl S<thomas.hackl@uni-wuerzburg.de>

=cut



