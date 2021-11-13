#! /usr/bin/env perl

use strict;
use warnings FATAL => 'all';
#use Carp::Always;

use FindBin;
my $PROGDIR = $FindBin::Bin;

# use lib "$FindBin::Bin";
# use Xyzzy;

use constant { TRUE => 1, FALSE => 0 };

# ------------------------------------------------------------------------

use File::Basename;
use Getopt::Std;

our $opt_h;
our $opt_v;

sub usage {
  my $progname = basename($0);
  print STDERR "Usage: $progname [options] ...\n";
  print STDERR "-h - print help\n";
  print STDERR "-v - verbose\n";
  exit(@_);
}

my $stat = getopts('hv');
if (!$stat) {
  usage(1);
}
if ($opt_h) {
  usage();
}
if (scalar(@ARGV) != 6) {
  usage(1);
}

my ($workdir,$proteome_gff,$firsts_blast,$best_breaks_txt,$orig_suffix,$new_suffix) = @ARGV;


# ------------------------------------------------------------------------

use constant { NO_VALUE => ";no-value;" };

sub parse_gff_attributes {
  my ($raw_attributes) = @_;
  my $attributes = {};
  foreach my $key_val (split(/; */,$raw_attributes)) {
    my ($key,$val);
    if ( $key_val =~ /^([^=]+)=(.*)/ ) {
      ($key,$val) = ($1,$2);
    } else {
      ($key,$val) = ($key_val, NO_VALUE);
    }
    $attributes->{$key} = $val;
  }
  return $attributes;
}

my $orf_details = {};

open (my $fh, "<", $proteome_gff) || die;
while (<$fh>) {
  chomp;
  if (/^#/) {
    next;
  }
  my ($accession,$source,$feature,$start,$end,
      $score,$strand,$frame,$raw_attributes) = split(/\t/,$_);
  my $attributes = parse_gff_attributes($raw_attributes);
  my $orf_name = $attributes->{name};
  my $crosses_origin = $attributes->{cross} || FALSE;
  (!defined($orf_details->{$orf_name})) || die;
  $orf_details->{$orf_name} = {
			       accession => $accession,
			       start => $start,
			       end => $end,
			       strand => $strand,
			       crosses => $crosses_origin
			      };
}
close $fh;

# ------------------------------------------------------------------------

my $bit_scores = {};

open ($fh, "<", $firsts_blast) || die;
while (<$fh>) {
  chomp;
  my ($query_id, $subj_id, $pers_id, $align_len, $mismatches,
      $gap_opens, $query_start, $query_end, $subj_start, $subj_end,
      $evalue, $bit_score) = split(/\t/);

  ($subj_start < $subj_end) || die;

  if (!defined($bit_scores->{$subj_id})) {
    $bit_scores->{$subj_id} = $bit_score;
  } elsif ($bit_scores->{$subj_id} < $bit_score) {
    $bit_scores->{$subj_id} = $bit_score;
  }

}
close $fh;

my $firsts = {};
foreach my $orf_name ( keys %{$bit_scores} ) {
  my $details = $orf_details->{$orf_name};
  (defined($details)) || die;
  my $accession = $details->{accession};
  if (!defined($firsts->{$accession})) {
    $firsts->{$accession} = [];
  }
  push @{$firsts->{$accession}}, $orf_name;
}

# ------------------------------------------------------------------------

my $best_breaks = {};

open ($fh, "<", $best_breaks_txt) || die;
while (<$fh>) {
  chomp;
  (/(.+): [0-9]+ \((.+)\)/) || die;
  $best_breaks->{$1} = [split(/,/,$2)];
}
close $fh;

# ------------------------------------------------------------------------

my @new_fastas;

use File::Basename;
foreach my $orig_fasta ( sort (<$workdir/*$orig_suffix>) ) {
  my $accession = basename($orig_fasta,$orig_suffix);
  my $new_fasta = "$workdir/$accession$new_suffix";
  if ($opt_v) {
    print STDERR "# $accession\n";
  }

  my ($pos,$strand);

  if (defined($firsts->{$accession})) {
    my @first_orf_names = @{$firsts->{$accession}};
    @first_orf_names = sort { $bit_scores->{$b} <=> $bit_scores->{$a} } @first_orf_names;
    if ($opt_v) {
      print STDERR "## first: ",join(", ",map { $_.":".$bit_scores->{$_} } @first_orf_names),"\n";
    }
    my $orf_name = $first_orf_names[0];
    my $details = $orf_details->{$orf_name};
    (defined($details)) || die;
    $strand = $details->{strand};
    if ( $strand eq "+" ) {
      $pos = $details->{start};
    } else {
      $pos = $details->{end};
    }
  } else {
    (defined($best_breaks->{$accession})) || die;
    my @best_breaks = @{$best_breaks->{$accession}};
    if ($opt_v) {
      print STDERR "## best_breaks: ",join(", ",@best_breaks),"\n";
    }
    $pos = $best_breaks[-1];
    $strand = "+";
  }

  if ($opt_v) {
    print STDERR "## new index 1: $pos($strand)\n";
  }

  print "$PROGDIR/reindex-fasta $strand $pos $orig_fasta > $new_fasta\n";
  push @new_fastas, $new_fasta;
}

print "cat ",join(" ",@new_fastas),"\n";
