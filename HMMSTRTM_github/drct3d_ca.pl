#!/usr/bin/perl -w
# drct3d.pl -- add pdb coordinates to drct file (binay sequence profile).
#  by Yu Shao
## BUG!!!! PDB resSeq numbers are not being used, but they should be.
##          C.B Wed Aug 11 15:02:23 EDT 2004

use strict;

# global variables:
my $drctfile;
my $pdbfile;
my $chainID;
my $pdbcode;

# drct i/o
my $template;
my $RECORDSIZE;
my $record;
my @fields;

### start of drct record structure ##
my $drctLen = 0;
my @drctSeq = ();
my @drctSeq_list = ();
my @drctResSeq = ();
my @drctStructure_name = ();
my @drctChainID = ();
my @drctSeq_resid = ();
my @drctSS = ();
my @drctRotamer = ();
my @drctXYZca = ();
my @drctXYZcb = ();
my @drctAng = ();
my @drctCsum = ();
my @drctVall_pro = ();
my @drctGap_freq = ();
my @drctIns_freq = ();
my @drctDomainID = ();
my @drctCtx = ();
my @drctClusterID = ();
### end of drct record structure ##

### pdb seq records
my $pdbLen = 0;
my $pdbLenAtom = 0;
my @pdbSeq = ();
my @pdbResSeq = ();
my @pdbSeqAtom = ();
my @pdbCoords = ();
my @pdbCoordsN = ();
my @pdbCoordsC = ();
my @pdbCoordsCB = ();
my @pdbCoordsCG = ();
my @pdbCoordsAtom = ();

my @dih = ();



die ("Usage: drct3d.pl drct pdb chainID pdbCode\n") unless scalar(@ARGV) >= 4;

$drctfile = shift;
$pdbfile = shift;

$chainID = shift;
if ($chainID eq "_" || $chainID eq "0") {
  $chainID = " ";
}
$pdbcode = shift;

print "drct3d.pl v.0.1 Yu Shao (shaoy\@rpi.edu) Oct 2003\n";
&readDRCT;
&readPDB;
&getdihedral;
&dynprog;
&writeDRCT;
print "done.\n";


# read the binary drct file
sub readDRCT {

  my $i;
  my $k;

  open (DRCT, $drctfile) || die ("$drctfile: file open failed.\n");
  binmode (DRCT);

  # 2 2 4 1 1 1 1 12 12 4 80 4 4 1 1 6
  # record size: 136 bytes
  # (0)seq_list (1)resSeq (2)structure_name (3)chainID (4)seq_resid (5)ss (6)rotamer (7-9)xyzca (10-12)ang (13)csum (14-33)vall_pro (34)gap_freq (35)ins_freq (36)domainID (37)ctx (38-40)clusterID
  # $template = "s s A4 A A A A f f f f f f f f f f f f f f f f f f f f f f f f f f f f f A A s3";
  $template = "s s A4 A A A A f3 f3 f f20 f f A A s3";
  $RECORDSIZE =  136;

  $drctLen = 0;
  until (eof(DRCT)) {
    read (DRCT, $record, $RECORDSIZE) == $RECORDSIZE || die ("short read\n");
    @fields = unpack ($template, $record);

#    print "$#fields  ";
#   for ($i = 0; $i <= $#fields; $i++) {
#     print "\"$fields[$i]\"  ";
#   }
#   print "\n";

    $drctSeq_list[$drctLen] = $fields[0];
    $drctResSeq[$drctLen] = $fields[1];
    $drctStructure_name[$drctLen] = $fields[2];
    $drctChainID[$drctLen] = $fields[3];
    $drctSeq[$drctLen] = $fields[4];
    $drctSS[$drctLen] = $fields[5];
    $drctRotamer[$drctLen] = $fields[6];
    for ($i = 7; $i < 10; $i++) {
      $drctXYZca[$drctLen][$i-7] = $fields[$i];
    }
    for ($i = 10; $i < 13; $i++) {
      $drctAng[$drctLen][$i-10] = $fields[$i];
    }
    $drctCsum[$drctLen] = $fields[13];
    for ($i = 14; $i < 34; $i++) {
      $drctVall_pro[$drctLen][$i-14] = $fields[$i];
    }
    $drctGap_freq[$drctLen] = $fields[34];
    $drctIns_freq[$drctLen] = $fields[35];
    $drctDomainID[$drctLen] = $fields[36];
    $drctCtx[$drctLen] = $fields[37];
    for ($i = 38; $i < 41; $i++) {
      $drctClusterID[$drctLen][$i-38] = $fields[$i];
    }

    $drctLen ++;

  # last;
  }
  close (DRCT);

# print drctSeq to output
  print "**************************************************\n";
  print "$drctfile (len = $drctLen):\n";
  print "--------------------------------------------------\n";
  $k = 0;
  for ($i = 0; $i < $drctLen; $i++) {
    print "$drctSeq[$i]";
    $k++;
    if ($k == 50) {
      print "\n";
      $k = 0;
    }
  }
  print "\n";
  print "**************************************************\n";
}

# read the pdb seq and coords for the given chain
sub readPDB {
  my @file = ();
  my $i;
  my $k;
  my ($j, $curSeq, $resSeq) = 0;
  my $atom;

  open (PDB, $pdbfile) || die ("$pdbfile: file open failed.\n");
  @file = <PDB>;
  close (PDB);

  $pdbLen = 0;
  $pdbLenAtom = 0;
  @pdbSeq = ();
#  cb
  @pdbResSeq = ();
  @pdbSeqAtom = ();
  @pdbCoords = ();
  @pdbCoordsAtom = ();

  $curSeq = -9999999999;
  for $i (0 .. 9999) {
    for $j (0 .. 2) {
      $pdbCoords[$i][$j] = 99999999;
      $pdbCoordsN[$i][$j] = 99999999;
      $pdbCoordsC[$i][$j] = 99999999;
      $pdbCoordsCB[$i][$j] = 99999999;
      $pdbCoordsCG[$i][$j] = 99999999;
    }
  }
  $j = -1;
  for ($i = 0; $i < $#file; $i++) {
    chomp ($file[$i]);
    if ($file[$i] =~ /^TER/ && $j > 0) {
      last;
    }
    if ($file[$i] !~ /^ATOM/ && $file[$i] !~ /^HETATM/) {
      next;
    }
    if (substr ($file[$i], 21, 1) ne $chainID) {
      next;
    }
#   if (substr ($file[$i], 12, 4) ne " CA ") {
#     next;
#   }

    # include the insertion code!!!!
    $resSeq = substr ($file[$i], 22, 5);
    if ($resSeq ne $curSeq) {
      $curSeq = $resSeq;
      $j++;
      # print "$file[$i] $resSeq $j\n";
    }
    $atom = substr ($file[$i], 12, 4);
    # note: some residue may have less or even only one atom
    #       in the pdb file. so read the residue name from
    #       each of them to avoid the uninitialized problem.
    if ($atom eq " CA ") {
      $pdbSeq[$j] = &three2one (substr($file[$i], 17, 3));
# cb 11-aug-04
      $pdbResSeq[$j] = (substr($file[$i], 22, 4));
       # print "$pdbSeq[$j] $pdbResSeq[$j] $j\n";
      $pdbCoords[$j][0] = substr ($file[$i], 30, 8);
      $pdbCoords[$j][1] = substr ($file[$i], 38, 8);
      $pdbCoords[$j][2] = substr ($file[$i], 46, 8);
#  cb  Use CA coordinates if CB is not found (i.e. for Gly)
      $pdbCoordsCB[$j][0] = substr ($file[$i], 30, 8);
      $pdbCoordsCB[$j][1] = substr ($file[$i], 38, 8);
      $pdbCoordsCB[$j][2] = substr ($file[$i], 46, 8);
    } elsif ($atom eq " N  ") {
      $pdbSeq[$j] = &three2one (substr($file[$i], 17, 3));
      $pdbCoordsN[$j][0] = substr ($file[$i], 30, 8);
      $pdbCoordsN[$j][1] = substr ($file[$i], 38, 8);
      $pdbCoordsN[$j][2] = substr ($file[$i], 46, 8);
    } elsif ($atom eq " C  ") {
      $pdbSeq[$j] = &three2one (substr($file[$i], 17, 3));
      $pdbCoordsC[$j][0] = substr ($file[$i], 30, 8);
      $pdbCoordsC[$j][1] = substr ($file[$i], 38, 8);
      $pdbCoordsC[$j][2] = substr ($file[$i], 46, 8);
    } elsif ($atom eq " CB ") {
      $pdbSeq[$j] = &three2one (substr($file[$i], 17, 3));
      $pdbCoordsCB[$j][0] = substr ($file[$i], 30, 8);
      $pdbCoordsCB[$j][1] = substr ($file[$i], 38, 8);
      $pdbCoordsCB[$j][2] = substr ($file[$i], 46, 8);
    } elsif ($atom eq " CG ") {
      $pdbSeq[$j] = &three2one (substr($file[$i], 17, 3));
      $pdbCoordsCG[$j][0] = substr ($file[$i], 30, 8);
      $pdbCoordsCG[$j][1] = substr ($file[$i], 38, 8);
      $pdbCoordsCG[$j][2] = substr ($file[$i], 46, 8);
    }
  }
  $pdbLen = $j + 1;
  
# for ($i = 0; $i < $pdbLen; $i++) {
#   print ("$pdbSeq[$i]    $pdbCoords[$i][0]   $pdbCoords[$i][1]    $pdbCoords[$i][2]\n"); 
# }

# print pdbSeq to output
  print "$pdbfile (len = $pdbLen):\n";
  print "--------------------------------------------------\n";
  $k = 0;
  for ($i = 0; $i < $pdbLen; $i++) {
    print "$pdbSeq[$i]";
    $k++;
    if ($k == 50) {
      print "\n";
      $k = 0;
    }
  }
  print "\n";
  print "**************************************************\n";

}

sub three2one {
  my $res = shift;
  if ($res eq "ALA") {
    return "A";
  } elsif ($res eq "CYS") {
    return "C";
  } elsif ($res eq "ASP") {
    return "D";
  } elsif ($res eq "GLU") {
    return "E";
  } elsif ($res eq "PHE") {
    return "F";
  } elsif ($res eq "GLY") {
    return "G";
  } elsif ($res eq "HIS") {
    return "H";
  } elsif ($res eq "ILE") {
    return "I";
  } elsif ($res eq "LYS") {
    return "K";
  } elsif ($res eq "LEU") {
    return "L";
  } elsif ($res eq "MET") {
    return "M";
  } elsif ($res eq "ASN") {
    return "N";
  } elsif ($res eq "PRO") {
    return "P";
  } elsif ($res eq "GLN") {
    return "Q";
  } elsif ($res eq "ARG") {
    return "R";
  } elsif ($res eq "SER") {
    return "S";
  } elsif ($res eq "THR") {
    return "T";
  } elsif ($res eq "VAL") {
    return "V";
  } elsif ($res eq "TRP") {
    return "W";
  } elsif ($res eq "TYR") {
    return "Y";
  } else {
    return "X";
  }
}

sub getdihedral {
  my ($i, $j) = 0;

  for $i (0 .. $pdbLen - 1) {
    for $j (0 .. 2) {
      $dih[$i][$j] = 0;
    }
  }

  for $i (0 .. $pdbLen - 1) {
    # phi (2 .. N)
    if ($i > 0) {
      if (!(
            ($pdbCoordsC[$i-1][0] == 99999999  && 
             $pdbCoordsC[$i-1][1] == 99999999  && 
             $pdbCoordsC[$i-1][2] == 99999999)
         || ($pdbCoordsN[$i][0] == 99999999    &&
             $pdbCoordsN[$i][1] == 99999999    &&
             $pdbCoordsN[$i][2] == 99999999)
         || ($pdbCoords[$i][0] == 99999999     &&
             $pdbCoords[$i][1] == 99999999     &&
             $pdbCoords[$i][2] == 99999999)
         || ($pdbCoordsC[$i][0] == 99999999    &&
             $pdbCoordsC[$i][1] == 99999999    &&
             $pdbCoordsC[$i][2] == 99999999)
         )) {
        $dih[$i][0] = &dihedral ($pdbCoordsC[$i-1], $pdbCoordsN[$i], 
                                 $pdbCoords[$i], $pdbCoordsC[$i]);
      }
    }
    # psi (1 .. N-1)
    if ($i < ($pdbLen - 1)) {
      if (!(
            ($pdbCoordsN[$i][0] == 99999999    && 
             $pdbCoordsN[$i][1] == 99999999    && 
             $pdbCoordsN[$i][2] == 99999999)
         || ($pdbCoords[$i][0] == 99999999     &&
             $pdbCoords[$i][1] == 99999999     &&
             $pdbCoords[$i][2] == 99999999)
         || ($pdbCoordsC[$i][0] == 99999999    &&
             $pdbCoordsC[$i][1] == 99999999    &&
             $pdbCoordsC[$i][2] == 99999999)
         || ($pdbCoordsN[$i+1][0] == 99999999  &&
             $pdbCoordsN[$i+1][1] == 99999999  &&
             $pdbCoordsN[$i+1][2] == 99999999)
         )) {
        $dih[$i][1] = &dihedral ($pdbCoordsN[$i], $pdbCoords[$i], 
                                 $pdbCoordsC[$i], $pdbCoordsN[$i+1]);
      }
    }
    # omega (1 .. N-1)
    if ($i < ($pdbLen - 1)) {
      if (!(
            ($pdbCoords[$i][0] == 99999999    && 
             $pdbCoords[$i][1] == 99999999    && 
             $pdbCoords[$i][2] == 99999999)
         || ($pdbCoordsC[$i][0] == 99999999   &&
             $pdbCoordsC[$i][1] == 99999999   &&
             $pdbCoordsC[$i][2] == 99999999)
         || ($pdbCoordsN[$i+1][0] == 99999999 &&
             $pdbCoordsN[$i+1][1] == 99999999 &&
             $pdbCoordsN[$i+1][2] == 99999999)
         || ($pdbCoords[$i+1][0] == 99999999  &&
             $pdbCoords[$i+1][1] == 99999999  &&
             $pdbCoords[$i+1][2] == 99999999)
         )) {
        $dih[$i][2] = &dihedral ($pdbCoords[$i], $pdbCoordsC[$i], 
                                 $pdbCoordsN[$i+1], $pdbCoords[$i+1]);
      }
    }
  }
}

# simple local dynamic programming using linear gap penalty and
# identity scoring matrix.
sub dynprog {
  my @score = ();
  my @ali = ();
  my @trace = ();
  my $match = 1;
  my $mismatch = -1;
  my $gap = -1;
  my $i;
  my $j;
  my $width;
  my $maxscore;
  my $maxposi;
  my $maxposj;
  my $alisize;
  my $nummis;
  my $numgap;

# initialize the scoring matrix.
  for ($i = 0; $i <= $drctLen; $i++) {
    $score[$i][0] = 0;
  }
  for ($i = 0; $i <= $pdbLen; $i++) {
    $score[0][$i] = 0;
  }
  for ($i = 1; $i <= $drctLen; $i++) {
    for ($j = 1; $j <= $pdbLen; $j++) {
      if ($drctSeq[$i-1] eq $pdbSeq[$j-1]) {
        $score[$i][$j] = $match;
      } else {
        $score[$i][$j] = $mismatch;
      }
    }
  }

# initialize the traceback matrix.
# [i-1][j]  -> 0
# [i][j]    -> 1
# [i][j-1]  -> 2
# terminate -> 3
  for ($i = 0; $i <= $drctLen; $i++) {
    for ($j = 0; $j <= $pdbLen; $j++) {
      $trace[$i][$j] = 0;
    }
  }

# local dynamic programming
  for ($i = 1; $i <= $drctLen; $i++) {
    for ($j = 1; $j <=$pdbLen; $j++) {
      if (($score[$i-1][$j-1] + $score[$i][$j]) >= ($score[$i-1][$j] + $gap) && ($score[$i-1][$j-1] + $score[$i][$j]) >= ($score[$i][$j-1] + $gap)) {
        $score[$i][$j] = $score[$i-1][$j-1] + $score[$i][$j];
        $trace[$i][$j] = 1;
      } elsif (($score[$i-1][$j] + $gap) >= ($score[$i-1][$j-1] + $score[$i][$j]) && ($score[$i-1][$j] + $gap) >= ($score[$i][$j-1] + $gap)) {
        $score[$i][$j] = $score[$i-1][$j] + $gap;
        $trace[$i][$j] = 0;
      } else {
        $score[$i][$j] = $score[$i][$j-1] + $gap;
        $trace[$i][$j] = 2;
      }
      if ($score[$i][$j] < 0) {
        $score[$i][$j] = 0;
        $trace[$i][$j] = 3;
      }
    }
  }

# traceback
  $maxscore = -1;
  $maxposi = 0;
  $maxposj = 0;
  for ($i = 0; $i <= $drctLen; $i++) {
    for ($j = 0; $j <= $pdbLen; $j++) {
      if ($score[$i][$j] >= $maxscore) {
        $maxscore = $score[$i][$j];
        $maxposi = $i;
        $maxposj = $j;
      }
    }
  }

  $alisize = 0;
  $i = $maxposi;
  $j = $maxposj;
  while ($i && $j) {
    if ($score[$i][$j] == 0) {
      last;
    }
    if ($trace[$i][$j] == 0) {
      $ali[$alisize][0] = $i - 1;
      $ali[$alisize][1] = -1;
      $i--;
    } elsif ($trace[$i][$j] == 1) {
      $ali[$alisize][0] = $i - 1;
      $ali[$alisize][1] = $j - 1;
      $i--;
      $j--;
    } elsif ($trace[$i][$j] == 2) {
      $ali[$alisize][0] = -1;
      $ali[$alisize][1] = $j - 1;
      $j--;
    }
    $alisize++;
  }

# for ($i = ($alisize - 1); $i >= 0; $i--) {
#   if ($ali[$i][0] == -1) {
#     print "_  ";
#   } else {
#     print "$drctSeq[$ali[$i][0]]  ";
#   }
#   if ($ali[$i][1] == -1) {
#     print "_  ";
#   } else {
#     print "$pdbSeq[$ali[$i][1]]  ";
#   }
#   print "\n";
# }
# print "$alisize\n";

# replace the dcrt coords and dihedral angles with those of pdb 
# based on the alignment
  for ($i = 0; $i < $alisize; $i++) {
    if ($ali[$i][0] != -1 && $ali[$i][1] != -1) {
# cb 11-aug-04
      $drctResSeq[$ali[$i][0]] = $pdbResSeq[$ali[$i][1]];
      for ($j = 0; $j < 3; $j++) {
        $drctXYZca[$ali[$i][0]][$j] = $pdbCoords[$ali[$i][1]][$j];
        $drctXYZcb[$ali[$i][0]][$j] = $pdbCoordsCB[$ali[$i][1]][$j];
        $drctAng[$ali[$i][0]][$j] = $dih[$ali[$i][1]][$j];
      }
    }
  }

# number of mismatches in the alignment
  $nummis = 0;
  $numgap = 0;
  my @mistag = ();
  for ($i = 0; $i < $alisize; $i++) {
    $mistag[$i] = 0;
  }
  for ($i = 0; $i < $alisize; $i++) {
    if ($drctSeq[$ali[$i][0]] eq "X" || $pdbSeq[$ali[$i][1]] eq "X") {
      next;
    }
    if ($ali[$i][0] == -1 || $ali[$i][1] == -1) {
      $numgap ++;
      next;
    }
    if ($drctSeq[$ali[$i][0]] ne $pdbSeq[$ali[$i][1]]) {
      $nummis++;
      $mistag[$i] = 1;
    }
  }

  print ("Alignment size: $alisize, score: $maxscore, mismatches: $nummis, gaps: $numgap\n");
  print "--------------------------------------------------\n";

# output the alignment to screen
  $i = $alisize - 1;
  while ($i >= 0) {
    if ($i >= 49) {
      $width = 50;
    } else {
      $width = $i + 1;
    }

    printf ("%6s%6d  ", "DRCT:", $ali[$i][0] + 1);
    for ($j = $i; $j > $i - $width; $j--) {
      if ($ali[$j][0] == -1 ) {
        printf ("%s", "-");
      } else {
        printf ("%s", $drctSeq[$ali[$j][0]]);
      }
    }
    for ($j = 0; $j < 50 - $width; $j++) {
      printf ("%s", " ");
    }
    printf ("  %6d\n", $ali[$i - $width + 1][0] + 1);

    print ("              ");
    for ($j = $i; $j > $i - $width ; $j--) {
      if ($mistag[$j] == 1) {
        print "+";
      } else {
        print " ";
      }
    }
    print "\n";

    printf ("%6s%6d  ", "PDB:", $ali[$i][1] + 1);
    for ($j = $i; $j > $i - $width ; $j--) {
      if ($ali[$j][1] == -1 ) {
        printf ("%s", "-");
      } else {
        printf ("%s", $pdbSeq[$ali[$j][1]]);
      }
    }
    for ($j = 0; $j < 50 - $width; $j++) {
      printf ("%s", " ");
    }
    printf ("  %6d\n\n", $ali[$i - $width + 1][1] + 1);

    $i = $i - $width;
  }
  print "**************************************************\n";
}


sub writeDRCT {
  
  my $i;
  my $j;
  my @entry = ();

  open (DRCT, ">$drctfile") || die ("$drctfile: output file open failed.\n");
  binmode (DRCT);
  
  for ($i = 0; $i < $drctLen; $i++) {
    ## @entry = ($drctSeq_list[$i], $drctResSeq[$i], $drctStructure_name[$i], $drctChainID[$i], $drctSeq[$i], $drctSS[$i], $drctRotamer[$i]);
    ## NEW VERSION uses pdbcode and chainID from command line  12-AUG-04
    @entry = ($drctSeq_list[$i], $drctResSeq[$i], $pdbcode, $chainID, $drctSeq[$i], $drctSS[$i], $drctRotamer[$i]);
    for ($j = 0; $j < 3; $j++) {
      @entry = (@entry, $drctXYZca[$i][$j]);
      ##  @entry = (@entry, $drctXYZcb[$i][$j]);
    }
    for ($j = 0; $j < 3; $j++) {
      @entry = (@entry, $drctAng[$i][$j]);
    }
    @entry = (@entry, $drctCsum[$i]);
    for ($j = 0; $j < 20; $j++) {
      @entry = (@entry, $drctVall_pro[$i][$j]);
    }
    @entry = (@entry, $drctGap_freq[$i], $drctIns_freq[$i], $drctDomainID[$i], $drctCtx[$i]);
    for ($j = 0; $j < 3; $j++) {
      @entry = (@entry, $drctClusterID[$i][$j]);
    }
    $record = pack ($template, @entry);
    print (DRCT $record);

#   print ("$#entry  ");
#   for ($j = 0; $j <= $#entry; $j++) {
#     print ("\"$entry[$j]\"  ");
#   }
#   print "\n";
  }

  close (DRCT);

}

##########################################
# dihedral angle calculating subroutines #
##########################################

# calculate the dihedral angle determined by 4 atoms in the space
sub dihedral {
  my ($vec1, $vec2, $vec3, $vec4) = @_;
  my (@vec21, @vec23, @vec34, @vecx, @vecy, @uvecx, @uvecy) = ();
  
  &subvec ($vec3, $vec2, \@vec23);
  &subvec ($vec1, $vec2, \@vec21);

  &crossprod (\@vec23, \@vec21, \@vecy);
  &unitvec (\@vecy, \@uvecy);

  &crossprod (\@uvecy, \@vec23, \@vecx);
  &unitvec (\@vecx, \@uvecx);

  &subvec ($vec4, $vec3, \@vec34);
# $cx = &dotprod (\@vec34, \@uvecx);
# $cy = &dotprod (\@vec34, \@uvecy);

  return atan2 (&dotprod (\@vec34, \@uvecy), &dotprod (\@vec34, \@uvecx)) / 3.1415926 * 180;
}

# vec3 = vec1 - vec2
sub subvec {
  my ($vec1, $vec2, $vec3) = @_;
  my $i;
  for ($i = 0; $i < 3; $i++) {
    $vec3->[$i] = $vec1->[$i] - $vec2->[$i];
  }
}

# vec3 is the cross-product of vec1 and vec2
sub crossprod {
  my ($vec1, $vec2, $vec3) = @_;
  $vec3->[0] = $vec1->[1] * $vec2->[2] - $vec1->[2] * $vec2->[1];
  $vec3->[1] = $vec1->[2] * $vec2->[0] - $vec1->[0] * $vec2->[2];
  $vec3->[2] = $vec1->[0] * $vec2->[1] - $vec1->[1] * $vec2->[0];
}

# return the dot-product of vec1 and vec2
sub dotprod {
  my ($vec1, $vec2) = @_;
  return $vec1->[0]*$vec2->[0] + $vec1->[1]*$vec2->[1] + $vec1->[2]*$vec2->[2]; 
}


# calculate the unit vector of vec1 and store it in vet2
sub unitvec {
  my ($vec1, $vec2) = @_;
  my $mod = sqrt (&dotprod ($vec1, $vec1));
  my $i;
  for $i (0 .. 2) {
    $vec2->[$i] = $vec1->[$i] / $mod;
  }
}
