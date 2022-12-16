#! /usr/bin/env perl

use strict;
use feature 'postderef';
use warnings FATAL => 'all';
no warnings "experimental::postderef";
#use Carp::Always;

# use FindBin;
# use lib "$FindBin::Bin";
# use Xyzzy;

use constant { TRUE => 1, FALSE => 0 };

# ------------------------------------------------------------------------

use File::Basename;
use Getopt::Std;

our $opt_d;
our $opt_C;
our $opt_S;
our $opt_D;
our $opt_h;
our $opt_v;

my $usage_str = "";
my $progname = basename($0);
$usage_str .= "Usage: $progname [options]\n";

# C:100.0%[S:99.8%,D:0.2%],F:0.0%,M:0.0%,n:440
# 440	Complete BUSCOs (C)
# 439	Complete and single-copy BUSCOs (S)
# 1	Complete and duplicated BUSCOs (D)
# 0	Fragmented BUSCOs (F)
# 0	Missing BUSCOs (M)
# 440	Total BUSCO groups searched


$usage_str .= "-C PERS - Only print C>=PERS\n";
$usage_str .= "-D PERS - Only print D<=PERS\n";
$usage_str .= "-S PERS - Only print S>=PERS\n";
$usage_str .= "-d DB - Limit to lineage database, DB\n";
$usage_str .= "-h - print help\n";
$usage_str .= "-v - print stats at the end\n";

sub usage {
  print STDERR $usage_str;
  exit(@_);
}

my $stat = getopts('C:D:S:d:hv');
if (!$stat) {
  usage(1);
}
if ($opt_h) {
  usage();
}
if (scalar(@ARGV) != 0) {
  usage(1);
}

# ------------------------------------------------------------------------

my $num_in = 0;
my $num_C_fail = 0;
my $num_D_fail = 0;
my $num_S_fail = 0;
my $num_d_fail = 0;
my $num_pass = 0;

while (<STDIN>) {
  chomp;
  my ($Name,$db,$C,$S,$D,$F,$M,$n) = split(/\t/);
  if ($Name eq "Name") {
    print "$_\n";
    next;
  }

  $num_in++;

  $C =~ s/\%$//;
  $S =~ s/\%$//;
  $D =~ s/\%$//;
  $F =~ s/\%$//;
  $M =~ s/\%$//;

  my $failed = FALSE;

  if ( defined($opt_C) && !($C>=$opt_C) ) {
    $failed = TRUE;
    $num_C_fail++;
  }
  if ( defined($opt_D) && !($D<=$opt_D) ) {
    $failed = TRUE;
    $num_D_fail++;
  }
  if ( defined($opt_S) && !($S>=$opt_S) ) {
    $failed = TRUE;
    $num_S_fail++;
  }
  if ( defined($opt_d) && !($db eq $opt_d) ) {
    $failed = TRUE;
    $num_d_fail++;
  }

  if ($failed) {
    next;
  }

  print "$_\n";
  $num_pass++;

}

if ( $opt_v ) {
  print STDERR "## Input genomes: $num_in\n";
  if ( $num_in != $num_pass ) {
    print STDERR "## Failed:\n";
  }
  if ( $num_C_fail > 0 ) {
    print STDERR "## - C>$opt_C: $num_C_fail\n";
  }
  if ( $num_S_fail > 0 ) {
    print STDERR "## - S>$opt_S: $num_S_fail\n";
  }
  if ( $num_D_fail > 0 ) {
    print STDERR "## - D<$opt_D: $num_D_fail\n";
  }
  if ( $num_d_fail > 0 ) {
    print STDERR "## - db!=$opt_d: $num_d_fail\n";
  }
  print STDERR "## Passed genomes: $num_pass\n";
}
