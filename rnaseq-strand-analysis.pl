#! /usr/bin/env perl

use strict;
use warnings FATAL => 'all';
#use Carp::Always;

# use FindBin;
# use lib "$FindBin::Bin";
# use Xyzzy;

use constant { TRUE => 1, FALSE => 0 };


# ------------------------------------------------------------------------

use File::Basename;
use Getopt::Std;

our $opt_h;
our $opt_o;

my $usage_str = "";
my $progname = basename($0);
$usage_str .= "Usage: $progname [options] forward.txt reverse.txt [locus_tag1 ...]\n";
$usage_str .= "-h - print help\n";
$usage_str .= "-o FILE - output to FILE\n";

sub usage {
  print STDERR $usage_str;
  exit(@_);
}

my $stat = getopts('ho:');
if (!$stat) {
  usage(1);
}
if ($opt_h) {
  usage();
}
if (scalar(@ARGV) < 2) {
  usage(1);
}

my ($forward_fn,$reverse_fn,@locus_tags) = @ARGV;

# ------------------------------------------------------------------------

my %all_locus_tags;

sub read_file {
  my ($fn) = @_;
  my $counts = {};

  open (my $fh,"<",$fn) || die "Cannot open <<$fn>>,";
  while (<$fh>) {
    if (/^#/) {
      next;
    }
    chomp;
    my ($Geneid,$Chr,$Start,$End,$Strand,$Length,$count) = split(/\t/);
    if ($Geneid eq "Geneid") {
      next;
    }
    $all_locus_tags{$Geneid} = TRUE;
    $counts->{$Geneid} = $count;
  }
  close $fh;

  return $counts;
}

my $forward_counts = read_file($forward_fn);
my $reverse_counts = read_file($reverse_fn);

# ------------------------------------------------------------------------

my $forward_vote = 0;
my $reverse_vote = 0;

my $pers_margin = .1;

sub vote {
  my ($tag,$forward_count,$reverse_count) = @_;
  my $result;
  if ( $forward_count > (1+$pers_margin)*$reverse_count ) {
    $forward_vote++;
    $result = "forward";
  } elsif ( $reverse_count > (1+$pers_margin)*$forward_count ) {
    $reverse_vote++;
    $result = "reverse";
  }
  print join("\t",$tag,$forward_count,$reverse_count,$result),"\n";
}

# ------------------------------------------------------------------------

my $all_forward = 0;
my $all_reverse = 0;

foreach my $locus_tag ( keys %all_locus_tags ) {
  $all_forward += $forward_counts->{$locus_tag};
  $all_reverse += $reverse_counts->{$locus_tag};
}

print join("\t","","forward","reverse","result"),"\n";
vote("all",$all_forward,$all_reverse);

# ------------------------------------------------------------------------

foreach my $locus_tag ( @locus_tags ) {
  (defined($forward_counts->{$locus_tag}))
    || die "No such locus tag <<$locus_tag>>\n";
  vote($locus_tag,$forward_counts->{$locus_tag},
       $reverse_counts->{$locus_tag});
}

# ------------------------------------------------------------------------

my $orientation;
if ($forward_vote > $reverse_vote) {
  $orientation = "forward";
} elsif ($forward_vote < $reverse_vote) {
  $orientation = "reverse";
} else {
  die "Tie.\n";
}

print "\n";
print "ORIENTATION=$orientation\n";
if ( $opt_o ) {
  open(my $out_fh,">",$opt_o) || die "Cannot open for writing <<$opt_o>>,";
  print $out_fh "ORIENTATION=$orientation\n";
  close $out_fh;
}

