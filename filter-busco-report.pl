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

our $opt_C;
our $opt_S;
our $opt_h;

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
$usage_str .= "-S PERS - Only print S>=PERS\n";
$usage_str .= "-h - print help\n";

sub usage {
  print STDERR $usage_str;
  exit(@_);
}

my $stat = getopts('C:S:h');
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

while (<STDIN>) {
  chomp;
  my ($Name,$db,$C,$S,$D,$F,$M,$n) = split(/\t/);
  if ($Name eq "Name") {
    print "$_\n";
    next;
  }
  $C =~ s/\%$//;
  $S =~ s/\%$//;
  $D =~ s/\%$//;
  $F =~ s/\%$//;
  $M =~ s/\%$//;
  if ( defined($opt_C) && !($C>=$opt_C) ) { next; }
  if ( defined($opt_S) && !($S>=$opt_S) ) { next; }
  print "$_\n";
}
