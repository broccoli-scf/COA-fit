#!/usr/bin/env perl
# 
# Parameterization program using pdl's lmfit (Levenberg-Marquardt non-linear least squares fitting).
# Code by Johannes Margraf and Duminda Ranasinghe (DFT version) (2016).
#
# IP Syntax: <Filename> IP <OrbitalNumber> <Value> <StdDev>
# AT Syntax: <Filename> AT <Number of Fragment Types> [<FragName> <Occurance>](<- repeat as necessary) <Value> <StdDev>
#
  use PDL;
  use PDL::Math;
  use PDL::Fit::LM;
  use strict;
  use warnings;

  our @namelist;
  our @iMO;
  our @refType;
  our @nFrag;
  our @FragName;
  our @FragQuan;
  our $lastname = "bacalala";
  my $scratch;
  my @temp;
  my $i;
  my $j;

  print "Starting parameterization... \n\n";

  open(INPUT,"<","ref.dat");

  # Read number of datapoints from first line
  our $nDat=<INPUT>;
  my  $xdata = pdl sequence $nDat;
  my  $ydata = zeroes($nDat);
  my  $wt    = ones($nDat);

  
  # Read data labels and reference data
  $i=0;
  while(<INPUT>){ 
    $scratch = $_;
    chomp $scratch;
    @temp =split / +/, $scratch;
    $namelist[$i] = $temp[0];
    $refType[$i]  = $temp[1];
    if($refType[$i] eq "IP"){
      $iMO[$i]           = $temp[2];    
      $ydata->slice($i) .= $temp[3];
      $wt->slice($i)    .= $temp[4];
      $scratch = $wt->slice($i);
    }elsif($refType[$i] eq "AT"){
      $nFrag[$i]          = $temp[2];
      for($j=0;$j<$nFrag[$i];$j++){
        $FragName[$i][$j] = $temp[3+(2*$j)];
        $FragQuan[$i][$j] = $temp[3+(2*$j)+1];
      }
      $ydata->slice($i)  .= $temp[3+2*$nFrag[$i]];  
      $wt->slice($i)     .= $temp[4+2*$nFrag[$i]];
      $scratch = $wt->slice($i);
    }else{
      print "<!> Unknown data type $refType[$i] <!>";
      exit;
    }
    $i++;
  }

  close INPUT;
  
  # Print Ref Data and StdDev
  print "Reference Data: \n";
  print "$ydata\n\n";
  print "Std. Devs: \n";
  print "$wt\n\n";

  # set initial prameters in a pdl (order in accord with fit function below)
  my $initp = pdl [0.0041757957,0.017762461,-0.0145164318];

  # Use lmfit. Fourth input argument is reference to user-defined
  # subroutine (\&costfn) detailed below.
  my ($yf,$pf,$cf,$if) = lmfit $xdata, $ydata, $wt, \&costfn, $initp, {Maxiter => 25, Eps => 1e-5};

  # Standard output
  print "\nXDATA\n$xdata\nY DATA\n$ydata\n\nY DATA FIT\n$yf\n\n";
  print "Fitted Parameters\n$pf\n\nCOVARIANCE MATRIX\n$cf\n\n";
  print "NUMBER ITERATIONS\n$if\n\n";



#=======================================================================================================
  sub costfn {
#=======================================================================================================
# Costfunction for LM optimizer. Use &runcalc to calculate current results and partial derivatives.
  # leave this line as is
    my ($x,$par,$ym,$dyda) = @_;
  # stepsize for numerical gradients
    my $step=0.0005;

  # name fit parameters for use in the calculation         
  # replace (0..1) with (0..x) where x is equal to your number of fit parameters minus 1
    my ($b1s,$b2s,$b2p) = map { $par->slice("($_)") } (0..2);
   
  # Write function with dependent variable $ym,
  # independent variable $x, and fit parameters as specified above.
  # Use the .= (dot equals) assignment operator to express the equality
  # (not just a plain equals)
    print "running stat calculations\n";
    $ym .= &runcalc($x,$b1s,$b2s,$b2p,1) ;
    print " $ym \n";

  # Edit only the (0..1) part to (0..x) as above
    my (@dy) = map {$dyda -> slice(",($_)") } (0..2);

#   TODO: write as loop!
    print "running grad calculations par1\n";
    my $dplus  = &runcalc($x,$b1s+$step,$b2s,$b2p,0);
    my $dminus = &runcalc($x,$b1s-$step,$b2s,$b2p,0);
    $dy[0] .= ($dplus-$dminus)/(2.0*$step);
#    print " $dy[0] \n";

    print "running grad calculations par2\n";
    $dplus  = &runcalc($x,$b1s,$b2s+$step,$b2p,0);
    $dminus = &runcalc($x,$b1s,$b2s-$step,$b2p,0);
    $dy[1] .= ($dplus-$dminus)/(2.0*$step);
#    print " $dy[1] \n";

    print "running grad calculations par3\n";
    $dplus  = &runcalc($x,$b1s,$b2s,$b2p+$step,0);
    $dminus = &runcalc($x,$b1s,$b2s,$b2p-$step,0);
    $dy[2] .= ($dplus-$dminus)/(2.0*$step);
#    print " $dy[2] \n";

#    $dplus  = &runcalc($x,$b1s,$b2s,$b2p,$expo+0.1,0);
#    $dminus = &runcalc($x,$b1s,$b2s,$b2p,$expo-0.1,0);
#    $dy[3] .= ($dplus-$dminus)/(2.0*0.1);
#    print " $dy[3] \n";
  }

#=======================================================================================================
  sub runcalc {
#=======================================================================================================
# Run calculations with current parameter set. Call parseout to get results.
    my ($x,$b1s,$b2s,$b2p,$print) = @_;
    my $results = zeroes($nDat);
    my $i;
    my $j;
    my $scratch;
    my @temp;
    my $name;
    my $res;
    my $docalc;
#    my $prog;

    if($print){print "$b1s $b2s $b2p\n"};
    for($i=0;$i<$nDat;$i++){
#      $prog=$i/$nDat*100.0;
      print "$i ";

      $name = $namelist[$i] ;
      if($name ne $lastname){
        $docalc = "true";
        $lastname = $name;
      }else{
        $docalc = "";
      }
#      open(PAR,">","PAR.TXT");
#      print PAR "$b1s $b2s $b2p 4.00 \n";
#      close PAR;
      # run molecular calculation
#      `/home/jmargraf/scripts/acesingen.sh $name`;
#      `/home/jmargraf/scripts/runCOA.sh $name`;
      if($docalc){
        `perl /Users/hans/prog/samsa/samsa.pl $name 19 $b1s $b2s $b2p`;
      }  
      if($refType[$i] eq "AT"){ 
      #  loop over number of fragments
        for($j=0;$j<$nFrag[$i];$j++){
      #    run fragment calculations
          `perl /Users/hans/prog/samsa/samsa.pl $FragName[$i][$j] 19 $b1s $b2s $b2p`;

#          `/home/jmargraf/scripts/acesingen.sh $FragName[$i][$j]`;
#          `/home/jmargraf/scripts/runCOA.sh $FragName[$i][$j]`;
        }
      }
      $res = &parseout($name,$i);
      $results->slice($i) .= $res;
    }
    print "\n";
    return $results;  
  }

#=======================================================================================================
  sub parseout {
#=======================================================================================================
# Parse output and get result of interest.
    my $name = $_[0];
    my $iDat = $_[1];
    my $scratch;
    my @temp;
    my $nHOMO;
    my $i;
    my $result;

    if($refType[$iDat] eq "IP"){
    # get orbital eigenvalue
      open(OUTPUT,"<","$name.out");
      while(<OUTPUT>){
        $scratch = $_;
        chomp $scratch;
        @temp =split / +/, $scratch ;
        if($temp[1] and ($temp[1] eq "No") and ($temp[3] eq "occ")){ # No of occ orbitals
          $nHOMO=$temp[6];
#          print "nHOMO = $nHOMO \n";
        }
#  Orbital Energies:
#    n      E[Ha]        E[eV]

#      0    -8.800059    -239.449611
        if($temp[1] and ($temp[1] eq "Orbital") and ($temp[2] eq "Energies:")){
          for($i=0;$i<($nHOMO+2);$i++){
            $scratch = <OUTPUT>;
            if($i==($iMO[$iDat]+2)){
              chomp $scratch;
              @temp =split / +/, $scratch ;
              $result = $temp[3];
#              print "E_HOMO($name) = $result \n";
            }
          }
          last;
        }
      }

      close(OUTPUT);
      return $result;

    }elsif($refType[$iDat] eq "AT"){
    # get moleculear energy
    # E_final = -144.016720857077
      my $Emol = 0.0;
      my $Eat  = 0.0;
      open(OUTPUT,"<","$name.out");
      while(<OUTPUT>){
        $scratch = $_;
        chomp $scratch;
        @temp =split / +/, $scratch ;
        if($temp[1] and ($temp[1] eq "E_final") ){
          $Emol = $temp[3];
          last;
        }
      }
      close(OUTPUT);
    # get atomic energies      
      for($i=0;$i<$nFrag[$iDat];$i++){
        open(OUTPUT,"<","$FragName[$iDat][$i].out");
        while(<OUTPUT>){
          $scratch = $_;
          chomp $scratch;
          @temp =split / +/, $scratch ;
          if($temp[1] and ($temp[1] eq "E_final") ){
            $Eat = $Eat + $FragQuan[$iDat][$i]*$temp[3];
#            print "$FragName[$iDat][$i] $FragQuan[$iDat][$i] $Eat  \n";
            last;
          }
        }
        close(OUTPUT);
      }
#      print "Emol = $Emol / Eat = $Eat \n";
      $result = -1.0*($Emol - $Eat)*627.5 ; 
      return $result;
    }

  }

