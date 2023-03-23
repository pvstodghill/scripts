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

# 'wide character' warning.
binmode(STDIN, ":utf8");
binmode(STDOUT, ":utf8");

# ------------------------------------------------------------------------

use File::Basename;
use Getopt::Std;

our $opt_F;
our $opt_M;
our $opt_c = 95.0;
our $opt_d;
our $opt_h;

my $usage_str = "";
my $progname = basename($0);
$usage_str .= "Usage: $progname [options]\n";

$usage_str .= "-F FASTANI.txt - FastANI formatted input\n";
$usage_str .= "-M MATRIX.txt - Matrix formatted input\n";
$usage_str .= "-c CUTOFF - cutoff for clade equivalence [$opt_c]\n";
$usage_str .= "-d OUTPUT.dot - graph of clades\n";
$usage_str .= "-h - print help\n";
$usage_str .= "\n";
$usage_str .= "exactly one required: -F -M\n";

sub usage {
  print STDERR $usage_str;
  exit(@_);
}

my $stat = getopts('F:M:c:d:h');
if (!$stat) {
  usage(1);
}
if ($opt_h) {
  usage();
}
if (scalar(@ARGV) != 0) {
  usage(1);
}
if (!$opt_F && !$opt_M) {
  usage(1);
}

# ------------------------------------------------------------------------
# Unison/Find functions transliterated from
# - https://en.wikipedia.org/wiki/Disjoint-set_data_structure
# ------------------------------------------------------------------------

sub uf_create {
  return {};
}

sub uf_ensure {
  my ($uf,$label) = @_;
  my $x = $uf->{$label};
  if (!defined($x)) {
    $x = {};
    $x->{label} = $label;
    $x->{parent} = $x;
    $x->{size} = 1;
    $uf->{$label} = $x;
  }
  return $x;
}

sub uf_find2 {
  my ($x) = @_;
  if ($x->{parent} != $x) {
    $x->{parent} = uf_find2($x->{parent});
    return $x->{parent};
  } else {
    return $x;
  }
}

sub uf_find {
  my ($uf,$label) = @_;
  my $x = uf_ensure($uf,$label);
  my $y = uf_find2($x);
  return $y->{label};
}

sub uf_union2 {
  my ($x,$y) = @_;
  $x = uf_find2($x);
  $y = uf_find2($y);
  if ($x == $y) {
    return;
  }
  if ($x->{size} < $y->{size}) {
    ($x,$y) = ($y,$x);
  }
  $y->{parent} = $x;
  $x->{size} += $y->{size};
}

sub uf_union {
  my ($uf,$label1,$label2) = @_;
  my $x = uf_ensure($uf,$label1);
  my $y = uf_ensure($uf,$label2);
  uf_union2($x,$y);
}

# ------------------------------------------------------------------------

sub is_T {
  my ($s) = @_;
  if ( $s =~ /^T__/ ) {
    return TRUE;
  } elsif ( $s =~ /\(T\)$/ ) {
    return TRUE;
  } else {
    return FALSE;
  }
}

sub is_TR {
  my ($s) = @_;
  if ( $s =~ /^TR__/ ) {
    return TRUE;
  } elsif ( $s =~ /\(T,R\)$/ ) {
    return TRUE;
  } else {
    return FALSE;
  }
}

sub is_R {
  my ($s) = @_;
  if ( $s =~ /^R__/ ) {
    return TRUE;
  } elsif ( $s =~ /\(R\)$/ ) {
    return TRUE;
  } else {
    return FALSE;
  }
}

sub fix_name {
  my ($s) = @_;
  $s =~ s/^TR__(.*)/$1(T,R)/;
  $s =~ s/^T__(.*)/$1(T)/;
  $s =~ s/^R__(.*)/$1(R)/;
  return $s;
}

# ------------------------------------------------------------------------

sub read_fastani {
  my ($filename) = @_;

  my $nodes = {};
  my $edges = {};
  my $uf = uf_create();

  open(my $in_fh, "<", $filename) || die "Cannot open for reading: <<$filename>>,";
  while (<$in_fh>) {
    chomp;
    my ($a,$b,$score,@rest) = split("\t");

    # $a = fix_name($a);
    # $b = fix_name($b);

    $nodes->{$a} = TRUE;
    $nodes->{$b} = TRUE;

    if ( $a eq $b ) {
      next;
    }
    if ( $score < $opt_c ) {
      next;
    }

    my $tag;
    if ( $a lt $b ) {
      $tag = "$a&&&$b";
    } else {
      $tag = "$b&&&$a";
    }

    my $l = $edges->{$tag};
    if (!defined($l)) {
      $l = $edges->{$tag} = [];
    }
    push $l->@*, $score;

    uf_union($uf,$a,$b);
  }

  return ($nodes,$edges,$uf);
}

# ------------------------------------------------------------------------

sub read_pyani {
  my ($filename) = @_;

  my $nodes = {};
  my $edges = {};
  my $uf = uf_create();

  my @col_names;
  open(my $in_fh, "<", $filename) || die "Cannot open for reading: <<$filename>>,";
  while (<$in_fh>) {
    chomp;
    my @fields = split(/\t/);

    if (scalar(@col_names) == 0) {
      shift @fields;
      @col_names = @fields;
      for (my $j=0; $j<=$#col_names; $j++) {
	my $name = $col_names[$j];
	# $name = fix_name($name);
	$nodes->{$name} = TRUE;
	$col_names[$j] = $name;
      }
      next;
    }

    my ($a,@col_values) = @fields;
    # $a = fix_name($a);

    for (my $j=0; $j<=$#col_values; $j++) {
      my $b = $col_names[$j];
      my $score = $col_values[$j];

      if ( $a eq $b ) {
	next;
      }
      if ( $score < $opt_c ) {
	next;
      }

      my $tag;
      if ( $a lt $b ) {
	$tag = "$a&&&$b";
      } else {
	$tag = "$b&&&$a";
      }

      my $l = $edges->{$tag};
      if (!defined($l)) {
	$l = $edges->{$tag} = [];
      }
      push $l->@*, $score;

      uf_union($uf,$a,$b);
    }
  }

  return ($nodes,$edges,$uf);
}

# ------------------------------------------------------------------------

my ($nodes,$edges,$uf);

if ( $opt_F ) {
  ($nodes,$edges,$uf) = read_fastani($opt_F);
} elsif ( $opt_M ) {
  ($nodes,$edges,$uf) = read_pyani($opt_M);
} else {
  die;
}

# ------------------------------------------------------------------------

my $clade_names = {};
my $clade_members = {};

my $suppress_member = {};


foreach my $node (keys($nodes->%*)) {
  my $color = uf_find($uf,$node);
  my $members = $clade_members->{$color};
  if (!defined($members)) {
    $members = $clade_members->{$color} = [];
  }
  push $members->@*, $node;
}

my $num_unnamed=0;

foreach my $color ( keys($clade_members->%*) ) {
  my @types;
  foreach my $member ( $clade_members->{$color}->@* ) {
    if ( is_T($member) || is_TR($member) ) {
      push @types, $member;
      $suppress_member->{$member} = TRUE;
    }
  }
  if (scalar(@types) == 0) {
    $num_unnamed++;
    $clade_names->{$color} = sprintf("~Unnamed #%d",$num_unnamed);
  } else {
    $clade_names->{$color} = join(", ", sort {$a cmp $b} @types);
  }
}

foreach my $color (sort { $clade_names->{$a} cmp $clade_names->{$b} } (keys($clade_names->%*))) {
  my $name = $clade_names->{$color};
  $name =~ s/^~//;
  print "### $name\n\n";
  my @members = $clade_members->{$color}->@*;
  my $skip = FALSE;
  foreach my $member ( sort {$a cmp $b} @members ) {
    if ($suppress_member->{$member}) { next; }
    print " - $member\n";
    $skip = TRUE;
  }
  if ($skip) {
    print "\n";
  }
}

# ------------------------------------------------------------------------

sub min {
  my ($x,$y) = @_;
  if (!defined($x)) {
    return $y;
  } elsif (!defined($y)) {
    return $x;
  } elsif ($x < $y) {
    return $x;
  } else {
    return $y;
  }
}

# ------------------------------------------------------------------------

if (!$opt_d) { exit; }

# dot language ref: https://graphviz.org/doc/info/lang.html

open(my $out_fh, ">", $opt_d) || die "Cannot open for writing: <<$opt_d>>,\n";

print $out_fh "strict graph {\n";

foreach my $edge ( sort {$a cmp $b} (keys($edges->%*)) ) {

  my ($a,$b) = split(/&&&/,$edge);
  my @scores = $edges->{$edge}->@*;
  my $weight = sprintf("%.2f",min(@scores));
  print $out_fh "$a -- $b [color=\"black\",weight=$weight]\n";

}

print $out_fh "\n";

foreach my $node ( sort {$a cmp $b} (keys($nodes->%*)) ) {

  my @attrs = sprintf("label=\"%s\"",fix_name($node));
  if (is_T($node)) {
    push @attrs, "color=\"green\"";
  } elsif (is_R($node)) {
    push @attrs, "color=\"red\"";
  } elsif (is_TR($node)) {
    push @attrs, "color=\"yellow\"";
  }

  my $attrs = join(", ", @attrs);

  print $out_fh "$node [$attrs]\n";
}

print $out_fh "}\n";
close($out_fh) || die;
