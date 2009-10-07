#!/usr/bin/perl -w

use strict;

#
# EDIT firstRunNumber
# EDIT enter VALUEs
# EDIT macro template to place values
# RUN script
#

my $firstRunNumber = 1000;

my $numberOfEventsPerRun = 10000;

my $minEnergy = 1;
my $maxEnergy = 1; # GeV
my $energyStep = 1;
my $minAngle = 0.1;
my $maxAngle = 1.0;
my $angleStep = 0.1;

my $condor_dir = "/home/home4/institut_1b/beischer/src/geant4/resolution/condor";
my $executable = "/home/home4/institut_1b/beischer/geant4/8.1.p02/slc4_ia32_gcc34/bin/Linux-g++/resolution";

################################################################################################################

my $currentRun = $firstRunNumber;

for (my $angle = $minAngle; $angle <= $maxAngle; $angle += $angleStep) {
    for (my $energy = $minEnergy; $energy <= $maxEnergy; $energy += $energyStep) {

    print "Run $currentRun:\n";
    my $condorfile = &make_condor_file($currentRun);
    print "$condorfile\n";
    my $macrofile = &make_macro_file($currentRun, $energy, $numberOfEventsPerRun);
    print "$macrofile\n";

    unlink("/home/home4/institut_1b/beischer/src/geant4/resolution/results/perdaix_${energy}_GeV_${angle}_deg.root");
    system "condor_submit", "$condorfile";
    ++$currentRun;
}

sub make_macro_file{
  # make a macro file, first argument is run number, second is current value (to be placed in macro template below)

  my ($energy, $angle, $nEvents) = @_;

  my $rotation = $angle/2.;

  open MACROFILE, ">perdaix_${energy}_GeV_${angle}_deg.mac" or die "Error: Cannot make macro file: $!";

  print MACROFILE <<EOF;
/control/verbose 1
/run/verbose 0
/tracking/verbose 0
/RES/Fit/Verbose 0

/RES/Gun/RandomOrigin
/RES/Gun/RandomDirection
/gun/energy ${energy} GeV

/control/alias moduleRot -${rotation}
/control/alias moduleInternalRot ${angle}

/RES/Det/AddModule 0. 0. 18.5 cm
/RES/Det/SetModuleRotation 0 {moduleRot}
/RES/Det/SetModuleInternalRotation 0 {moduleInternalRot}
/RES/Det/SetModuleWidth 0 20.736

/RES/Det/AddModule 0. 0. 4.5 cm
/RES/Det/SetModuleRotation 1 {moduleRot}
/RES/Det/SetModuleInternalRotation 1 {moduleInternalRot}
/RES/Det/SetModuleWidth 1 13.824

/RES/Det/AddModule 0. 0. -4.5 cm
/RES/Det/SetModuleRotation 2 {moduleRot}
/RES/Det/SetModuleInternalRotation 2 {moduleInternalRot}
/RES/Det/SetModuleWidth 2 13.824

/RES/Det/AddModule 0. 0. -18.5 cm
/RES/Det/SetModuleRotation 3 {moduleRot}
/RES/Det/SetModuleInternalRotation 3 {moduleInternalRot}
/RES/Det/SetModuleWidth 3 20.736

#/RES/Field/SetInhomFieldFrom tables/perdaix_07_jul_2009.table
/RES/Field/SetDummyField 0.27 0.0 0.0 tesla
#/RES/Field/SetUniformField 0.3 0.0 0.0 tesla

/RES/Fit/Method blobel

/RES/Data/OverWriteFile true
/RES/Data/SetFileName /home/home4/institut_1b/beischer/src/geant4/resolution/results/perdaix_${energy}_GeV_${angle}_deg.root
/RES/Run/StoreResults

/run/initialize

/RES/Run/Generate ${nEvents}
/RES/Run/Reconstruct
EOF

  close MACROFILE;

  # return value: macro file name
  "perdaix_${energy}_GeV_${angle}_deg.mac";

}

sub make_condor_file{
  # make a condor file, first ( and only ) argument is run number

  my $ith = shift;
  my $currentRun = $ith;

  open CONDORFILE, ">${condor_dir}/condor/res_${currentRun}.condor" or die "Error: Cannot make condor file: $!";

  print CONDORFILE <<EOF;
universe        = vanilla
executable      = $executable
arguments	= ${condor_dir}/mac/res_$currentRun.mac
initialdir      = ${condor_dir}
EOF

  print CONDORFILE <<'EOF';
getenv          = true
requirements    = ( OpSys == "LINUX" )
rank            = KFlops
EOF

  print CONDORFILE <<EOF;
error           = ${condor_dir}/ERR/res_${currentRun}.ERR
output          = ${condor_dir}/STD/res_${currentRun}.STD
log             = ${condor_dir}/LOG/res_${currentRun}.LOG
EOF

  print CONDORFILE <<'EOF';
notification    = error
notify_user	= beischer@physik.rwth-aachen.de
coresize        = 1000
Queue
EOF

  close CONDORFILE;

  # return value: condor file name
  "${condor_dir}/condor/res_${currentRun}.condor";

}
