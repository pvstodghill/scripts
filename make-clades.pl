#! /usr/bin/env perl

use strict;
use feature 'postderef';
use warnings FATAL => 'all';
no warnings "experimental::postderef";
#use Carp::Always;

use File::Basename;

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

our $opt_A;
our $opt_F;
our $opt_L;
our $opt_M;
our $opt_N;
our $opt_U = "";
our $opt_c = 95.0;
our $opt_d;
our $opt_h;
our $opt_s;
our $opt_t;

my $usage_str = "";
my $progname = basename($0);
$usage_str .= "Usage: $progname [options]\n";

$usage_str .= "-A MASH.out - raw Mash output\n";
$usage_str .= "-F FASTANI.txt - FastANI formatted input\n";
$usage_str .= "-L LENS - replicon lengths\n";
$usage_str .= "-M MATRIX.txt - Matrix formatted input\n";
$usage_str .= "-N CLADE_NAMES.txt - tsv: REPLICON_NAME CLADE_NAME\n";
$usage_str .= "-U TAG - prefix for unnamed clades\n";
$usage_str .= "-c CUTOFF - cutoff for clade equivalence [$opt_c]\n";
$usage_str .= "-d OUTPUT.dot - graph of clades\n";
$usage_str .= "-h - print help\n";
$usage_str .= "-s - don't print singletons\n";
$usage_str .= "-t - output tsv\n";

$usage_str .= "\n";
$usage_str .= "exactly one required: -A -F -M\n";

sub usage {
  print STDERR $usage_str;
  exit(@_);
}

my $stat = getopts('A:F:L:M:N:U:c:d:hts');
if (!$stat) {
  usage(1);
}
if ($opt_h) {
  usage();
}
if (scalar(@ARGV) != 0) {
  usage(1);
}
if ((defined($opt_A) ? 1 : 0)
    + (defined($opt_F) ? 1 : 0)
    + (defined($opt_M) ? 1 : 0)
    != 1) {
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

sub read_mash {
  my ($filename) = @_;

  my $nodes = {};
  my $edges = {};
  my $uf = uf_create();

  open(my $in_fh, "<", $filename) || die "Cannot open for reading: <<$filename>>,";
  while (<$in_fh>) {
    chomp;
    my ($a,$b,$score,@rest) = split("\t");

    $a = basename($a,".fna",".fasta");
    $b = basename($b,".fna",".fasta");
    # $a = fix_name($a);
    # $b = fix_name($b);

    $nodes->{$a} = TRUE;
    $nodes->{$b} = TRUE;

    if ( $a eq $b ) {
      next;
    }
    if ( $score > $opt_c ) {
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

if ( $opt_A ) {
  ($nodes,$edges,$uf) = read_mash($opt_A);
} elsif ( $opt_F ) {
  ($nodes,$edges,$uf) = read_fastani($opt_F);
} elsif ( $opt_M ) {
  ($nodes,$edges,$uf) = read_pyani($opt_M);
} else {
  die;
}

# ------------------------------------------------------------------------

my $replicon_lengths = {};

if ( $opt_L ) {
  open(my $len_fh,"<",$opt_L) || die "Cannot open: <<$opt_L>>,";
  while (<$len_fh>) {
    chomp;
    if ($_ eq "") { next; }
    my ($replicon_name,$replicon_length) = split(/\t/);
    (!defined($replicon_lengths->{$replicon_name})) || die;
    $replicon_lengths->{$replicon_name} = $replicon_length;
  }
  close $len_fh;
}

# ------------------------------------------------------------------------

my $clade_names = {};

if ( $opt_N ) {
  open(my $names_fh,"<",$opt_N) || die "Cannot open: <<$opt_N>>,";
  while (<$names_fh>) {
    chomp;
    if ($_ eq "") { next; }
    my ($replicon_name,$clade_name) = split(/\t/);
    (!defined($clade_names->{$replicon_name})) || die;
    $clade_names->{$replicon_name} = $clade_name;
  }
  close($names_fh);
}

# ------------------------------------------------------------------------

my $color_members = {};

foreach my $node (keys($nodes->%*)) {
  my $color = uf_find($uf,$node);
  my $members = $color_members->{$color};
  if (!defined($members)) {
    $members = $color_members->{$color} = [];
  }
  push $members->@*, $node;
}

# ------------------------------------------------------------------------

my @all_colors = keys($color_members->%*);

my $clade_size = {};
my $clade_key = {};

foreach my $color ( @all_colors ) {
  my $len = 0;
  my $key;
  foreach my $member ( $color_members->{$color}->@* ) {
    my $n = $replicon_lengths->{$member};
    if (defined($n)) {
      $len += $n;
    } else {
      $len += 1;
    }
    if (!defined($key) || $member lt $key) {
      $key = $member;
    }
  }
  $clade_size->{$color} = $len;
  $clade_key->{$color} = $key;
}

my @sorted_colors = sort {
  $clade_size->{$b} <=> $clade_size->{$a}
    || $clade_key->{$a} cmp $clade_key->{$b}
  } @all_colors;

# ------------------------------------------------------------------------

my $color_names = {};
my $num_unnamed=0;

foreach my $color ( @sorted_colors ) {
  my $color_name;
  foreach my $member ( $color_members->{$color}->@* ) {
    my $name = $clade_names->{$member};
    if (!defined($name)) {
      ;
    } elsif ( !defined($color_name) ) {
      $color_name = $name;
    } elsif ( $color_name eq $name ) {
      ;
    } else {
      die "Multiple names for clade: <<$color_name>>, <<$name>>,";
    }
  }

  if (!defined($color_name)) {
    $num_unnamed++;
    $color_name = sprintf("Unnamed %s%d",$opt_U,$num_unnamed);
  }
  $color_names->{$color} = $color_name;
}

# ------------------------------------------------------------------------

my $num_colors = 0;

foreach my $color ( @sorted_colors ) {
  my $name = $color_names->{$color};
  my @members = $color_members->{$color}->@*;

  if ($opt_s && scalar(@members) == 1) { next; }

  if ( $opt_t ) {

    my $name2 = $name;
    $name2 =~ s/^Unnamed //;
    foreach my $member ( sort {$a cmp $b} @members ) {
      my $len = $replicon_lengths->{$member};
      if (defined($len)) {
	print join("\t",$name2,$member,$len),"\n";
      } else {
	print join("\t",$name2,$member),"\n";
      }
    }

  } else {

    print "### $name\n\n";
    my $skip = FALSE;
    foreach my $member ( sort {$a cmp $b} @members ) {
      my $len = $replicon_lengths->{$member};
      if (defined($len)) {
	print " - $member [$len]\n";
      } else {
	print " - $member\n";
      }
      $skip = TRUE;
    }
    if ($skip) {
      print "\n";
    }
  }

  $num_colors++;

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
