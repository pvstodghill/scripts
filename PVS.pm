package PVS;

use strict;
use warnings;

# Cribbed from *Man perlmod*
BEGIN {
  require Exporter;

  # Inherit from Exporter to export functions and variables
  our @ISA         = qw(Exporter);

  # Functions and variables which are exported by default
  our @EXPORT      = qw(min max);

  # Functions and variables which can be optionally exported
  # our @EXPORT_OK   = qw($Var1 %Hashit func3);
}

########################################################################
# The correct way to include this file is as follows,
#
# use FindBin;
# use lib "$FindBin::Bin";
# use PVS;
########################################################################

########################################################################

sub read_fasta {
    my ($filename) = @_;
    my %h;
    my $defline;
    my $seq;
    my ($fh,$stdin_p);
    if ($filename eq "-") {
      $stdin_p = 1;
      $fh = *STDIN;
    } elsif ($filename =~ /(.*\.tgz):(.*)/) {
      my ($tarfile,$filename2) = ($1,$2);
      $stdin_p = 0;
      (open $fh,"tar zxOf $tarfile $filename2 |") 
	|| die "Can't open $filename\n";
    } else {
      $stdin_p = 0;
      (open $fh,"<$filename")
	|| die "Can't open $filename\n";
    }
    while (<$fh>) {
        s/[\n\r]+$//;
	if ( $_ =~ /^>(.*)/ ) {
	    if ( defined($defline) ) {
		$h{$defline} = $seq;
	    }
	    $defline = $1;
	    $seq = "";
	} else {
	    $seq .= $_;
	}
    }
    if ( defined($defline) ) {
	$h{$defline} = $seq;
    }
    if ( !$stdin_p ) {
      close $fh;
    }
    if (wantarray) {
      return %h;
    } else {
      my @seqs = values(%h);
      return $seqs[0];
    }
}

sub read_fasta2 {
    my ($filename) = @_;
    my @results;
    my $defline;
    my $seq;
    (open my $fh,"<$filename") || die;
    while (<$fh>) {
	chomp $_;
	if ( $_ =~ /^>(.*)/ ) {
	    if ( defined($defline) ) {
		push @results, $defline;
		push @results, $seq;
	    }
	    $defline = $1;
	    $seq = "";
	} else {
	    $seq .= $_;
	}
    }
    if ( defined($defline) ) {
		push @results, $defline;
		push @results, $seq;
    }
    close $fh;
    return @results;
}

########################################################################

my $fasta_width = 70;

sub write_fasta {
  my ($fh, @args) = @_;

  while ($#args >= 0) {
    my $defline = shift @args;
    my $seq = shift @args;
    (defined($defline) && defined($seq)) || die;
    print $fh ">$defline\n";
    my $len = length($seq);
    my $i;
    for ($i=0; $i<$len; $i += $fasta_width) {
      print $fh substr($seq,$i,$fasta_width),"\n";
    }
  }
}

########################################################################

# From, http://www.animalgenome.org/blast/docs/fasta.html
# 
# Sequence:   ACGTMRWSYKVHDBN.
# Complement: TGCAKYWSRMBDHVN.

sub complement {
    my ($s) = @_;
    ( $s !~ /([^ACGTMRWSYKVHDBN.])/i ) || die "Can't complement $1\n";
    $s =~ tr/acgtmrwsykvhdbnACGTMRWSYKVHDBN./tgcakywsrmbdhvnTGCAKYWSRMBDHVN./;
    return $s;
}

sub reverse {
    my ($s) = @_;
    $s = CORE::reverse($s);
    return $s;
}

sub reverse_complement {
    my ($s) = @_;
    $s = complement($s);
    $s = &reverse($s);
    return $s;
}

########################################################################

sub average {
  my $count = 0;
  my $sum = 0.0;
  foreach my $x ( @_ ) {
    if ( !defined($x) ) { next; }
    $count++;
    $sum += $x;
  }
  if ($count > 0) {
    return $sum / $count;
  } else {
    return undef;
  }
}

########################################################################

# Radix trees.

sub radix_tree_init {
  my ($alphabet) = @_;
  my $n = length($alphabet);
  my @decode;
  for (my $i=0; $i<$n; $i++) {
    my $c = substr($alphabet,$i,1);
    $decode[ord(lc($c))] = $i+1;
    $decode[ord(uc($c))] = $i+1;
  }
  return [undef,$n,[@decode]];
}

sub radix_tree_insert {
  my ($tree,$string) = @_;
  my ($bins,$n,$decode) = @$tree;
  ${$tree}[0] = xradix_tree_insert($string, 0, $bins, $n, $decode);
}

sub xradix_tree_insert {
  my ($string, $i, $bins, $n, $decode) = @_;
  if (!defined($bins)) {
    $bins = make_bins($n);
  }
  ${$bins}[0]++;
  if ($i == length($string)) { 
    ;
  } else {
    my $c = substr($string,$i,1);
    my $j = ${$decode}[ord($c)];
    ${$bins}[$j] = xradix_tree_insert($string,$i+1,${$bins}[$j], $n, $decode);
  }
  return $bins;
}

sub make_bins {
  my ($n) = @_;
  my @bins = (0);
#   while ($n > 0) {
#     push @bins, undef;
#     $n--;
#   }
  return [@bins];
}

sub radix_tree_search {
  my ($tree,$string) = @_;
  my ($bins,$n,$decode) = @$tree;
  return xradix_tree_search($string, 0, $bins, $n, $decode);
}

sub xradix_tree_search {
  my ($string, $i, $bins, $n, $decode) = @_;
  (defined($bins)) || die;
  if ($i == length($string)) {
    return ($i, ${$bins}[0]);
  } else {    
    my $c = substr($string,$i,1);
    my $j = ${$decode}[ord($c)];
    if (!defined(${$bins}[$j])) {
      return ($i, ${$bins}[0]);
    } else {
      my $bins2 = ${$bins}[$j];
      return xradix_tree_search($string,$i+1,$bins2,$n,$decode);
    }
  }
}




########################################################################

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


sub max {
  my ($x,$y) = @_;
  if (!defined($x)) {
    return $y;
  } elsif (!defined($y)) {
    return $x;
  } elsif ($x > $y) {
    return $x;
  } else {
    return $y;
  }
}

########################################################################

sub logn {
  my ($m,$n) = @_;
  return log($m)/log($n);
}

sub log2 {
  my ($m) = @_;
  return log($m)/log(2.0);
}
sub log10 {
  my ($m) = @_;
  return log($m)/log(10.0);
}

########################################################################

# Make arrays of 0's of the lengths provided.

sub make_replicon_arrays {
  my ($init,%lengths) = @_;
  my %replicons;
  foreach my $key (keys %lengths) {
    my @a = ($init);
    my $len = $lengths{$key};
    for (my $i=1; $i<=$len+1; $i++) {
      push @a, $init;
    }
    $replicons{$key} = [ @a ];
    # There is an extra 0 at the beginning and end of each replicon array
    # (i.e., index 0 and index length+1). This simplifies certain operations.
  }
  return %replicons;
}

our %DC3000_lengths = (
		      'NC_004578', 6397126,
		      'NC_004632', 67473,
		      'NC_004633', 73661,
		     );

########################################################################

sub compute_pers_cg {
  my $cg = 0;
  my $at = 0;
  foreach my $s ( @_ ) {

    for (my $i = 0; $i < length($s); $i++) {
      my $c = lc(substr($s,$i,1));
      if ( $c eq "c" || $c eq "g" ) {
	$cg++;
      } elsif ( $c eq "a" || $c eq "t" ) {
	$at++;
      } else {
	die "c = $c\n";
      }
    }


  }

  my $total = $cg + $at;
  my $result = 100.0*((1.0*$cg) / (1.0*$total));

  return $result;
}

########################################################################
# Transform gff coordinates by formulae
#
# old_transform_gff_coords("start-5","end+15",10,20,"+")
# -> (5,35)
# old_transform_gff_coords("start-5","end+15",100,200,"-")
# -> (85,205)
########################################################################

sub parse_gff_xform_formula {
  my ($s) = @_;
  $s =~ s/[ \t]+//g;
  my ($from_start,$delta);
  if ( $s =~ /^start(\W*.*)/ ) {
    $from_start = 1;
    $delta = $1;
  } elsif ( $s =~ /^end(\W*.*)/ ) {
    $from_start = 0;
    $delta = $1;
  } else {
    print STDERR "Ill-formed coordinate: $s\n";
    exit(1);
  }
  $delta =~ s/^(-?)[+](-?)/$1$2/;
  if ( $delta eq "" ) { $delta = 0; }
  return ($from_start,$delta);
}

sub old_transform_gff_coord {
  my ($from_start,$delta,$start,$end) = @_;
  my $x;
  if ( $from_start ) {
    $x = $start;
  } else {
    $x = $end;
  }
  $x += $delta;
  return $x;
}

sub old_transform_gff_coords {
  my ($start_formula,$end_formula,$lb1,$ub1,$strand) = @_;
  my $factor = 1;
  my ($feature_start,$feature_end) = ($lb1,$ub1);
  if ($strand eq "-") {
    $factor = -1;
    ($feature_start,$feature_end) = ($ub1,$lb1);
  }
  my ($start_from_start,$start_delta) = parse_gff_xform_formula($start_formula);
  my ($end_from_start,$end_delta) = parse_gff_xform_formula($end_formula);
  my $xxx_start = old_transform_gff_coord($start_from_start,$factor*$start_delta,
			    $feature_start,$feature_end);
  my $xxx_end = old_transform_gff_coord($end_from_start,$factor*$end_delta,
			  $feature_start,$feature_end);

  my ($lb2,$ub2) = ($xxx_start,$xxx_end);
  if ($strand eq "-") {
    ($lb2,$ub2) = ($xxx_end,$xxx_start);
  }
  # fixme: check for empty seq or error?
  ($lb2 <= $ub2) || die;
  # fixme: circular wrap-around
  return ($lb2,$ub2);
}

########################################################################

sub eval_gff_coords {
  my ($formula,$start,$end,$length) = @_;
  my $expr = $formula;
  $expr =~ s/\bstart\b/ ($start) /g;
  $expr =~ s/\bend\b/ ($end) /g;
  $expr =~ s/\blength\b/ ($length) /g;
  my $result = eval $expr;
  $result = int($result); # ick
  return $result;
}

sub transform_gff_coords {
  my ($start_formula,$end_formula,$flip_strand,$lb1,$ub1,$strand1) = @_;

  my ($start1,$end1);
  if ( $strand1 eq "+" || $strand1 eq "." ) {
    ($start1,$end1) = ($lb1,$ub1);
  } elsif ( $strand1 eq "-") {
    ($start1,$end1) = (-$ub1,-$lb1);
  } else {
    die "transform_gff_coords: strand = <$strand1>,"
  }

  my $length1 = $ub1 - $lb1 + 1;

  my $start2 = eval_gff_coords($start_formula,$start1,$end1,$length1);
  my $end2 = eval_gff_coords($end_formula,$start1,$end1,$length1);

  my ($lb2,$ub2);
  if ( $strand1 eq "+" || $strand1 eq "." ) {
    ($lb2,$ub2) = ($start2,$end2);
  } elsif ( $strand1 eq "-") {
    ($lb2,$ub2) = (-$end2,-$start2);
  } else {
    die "transform_gff_coords: strand = <$strand1>,"
  }
  # fixme: check for empty seq or error?
  ($lb2 <= $ub2) || die;
  # fixme: circular wrap-around

  my $strand2;
  if (!$flip_strand) {
    $strand2 = $strand1;
  } elsif ( $strand1 eq "+" ) {
    $strand2 = "-";
  } elsif ( $strand1 eq "-" ) {
    $strand2 = "+";
  } elsif ( $strand1 eq "." ) {
    die;
  } else {
    die;
  }
  return ($lb2,$ub2,$strand2);
}

########################################################################

1;
