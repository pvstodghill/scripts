#! /usr/bin/env perl

use strict;
use warnings FATAL => 'all';
use Carp::Always;

# use FindBin;
# use lib "$FindBin::Bin";
# use Xyzzy;

use constant { TRUE => 1, FALSE => 0 };

# ------------------------------------------------------------------------
# Process the command line
# ------------------------------------------------------------------------

use File::Basename;
use Getopt::Std;

our $opt_h;

sub usage {
  my $progname = basename($0);
  print STDERR "Usage: $progname [options] a.afa b.afa ...\n";
  print STDERR "-h - print help\n";
  exit(@_);
}

my $stat = getopts('h');
if (!$stat) {
  usage(1);
}
if ($opt_h) {
  usage();
}
if (scalar(@ARGV) == 0) {
  usage(1);
}

# ------------------------------------------------------------------------

my %strains;
my @sequences;

sub emit {
  my ($seqs,$seq,$defline) = @_;
  if (!defined($defline)) {
    return;
  }
  ($defline =~ /^>([^ ]*)/) || die "defline=\"$defline\",";
  my $strain = $1;
  $strains{$strain} = TRUE;
  (!defined($seqs->{$strain})) || die "more than one sequence: $strain,";
  $seqs->{$strain} = $seq;
}

foreach my $afa_name ( @ARGV ) {
  my $defline = undef;
  my $seq = "";
  my $seqs = {};
  open(my $afa_fh,"<", $afa_name) || die "cannot open <<$afa_name>>,";
  while (<$afa_fh>) {
    chomp;
    if ( /^>/ ) {
      emit($seqs,$seq,$defline);
      $defline = $_;
      $seq = "";
    } else {
      s/ //g;
      $seq .= $_;
    }
  }
  emit($seqs,$seq,$defline);
  close $afa_fh;
  push @sequences, $seqs;
}

# ------------------------------------------------------------------------

use constant { WIDTH => 70 };

strain:
foreach my $strain (sort (keys %strains)) {
  my $seq = "";
  foreach my $seqs ( @sequences ) {
    my $subseq = $seqs->{$strain};
    if (!defined($subseq)) {
      print STDERR "# Dropping $strain\n";
      next strain;
    }
    $seq .= $subseq;
  }
  print ">$strain\n";
  while (length($seq) > WIDTH) {
    my $s = substr($seq,0,WIDTH);
    $seq = substr($seq,WIDTH);
    print "$s\n";
  }
  if (length($seq) > 0) {
    print "$seq\n";
  }
}
