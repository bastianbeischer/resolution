#!/usr/bin/perl -w

use strict;

#
# EDIT firstRunNumber
# EDIT enter VALUEs
# EDIT macro template to place values
# RUN script
#

my $firstRunNumber = 1600;

my $numberOfEventsPerRun = 10000;

my $minMomentum = 0.5;
my $maxMomentum = 10.0;
my $momentumStep = 0.5;
my $minAngle = 5.0;
my $maxAngle = 5.0;
my $angleStep = 0.1;

my $project_dir = "/home/home4/institut_1b/beischer/src/geant4/resolution";
my $condor_dir  = "${project_dir}/condor";
my $table_dir   = "${project_dir}/tables";
my $result_dir  = "${project_dir}/results";
my $executable  = "$ENV{G4BIN}/$ENV{G4SYSTEM}/resolution";

################################################################################################################

my $currentRun = $firstRunNumber;

for (my $angle = $minAngle; $angle <= $maxAngle; $angle += $angleStep) {
    for (my $momentum = $minMomentum; $momentum <= $maxMomentum; $momentum += $momentumStep) {

    my $momentumString = sprintf("%.1f", $momentum);
    my $angleString = sprintf("%.2f", $angle);
    my $filename = "${result_dir}/perdaix_${momentumString}_GeV_${angleString}_deg_inhom_msc.root";

    print "Run $currentRun:\n";
    my $condorfile = &make_condor_file($currentRun);
    print "$condorfile\n";
    my $macrofile = &make_macro_file($currentRun, $momentum, $angle, $numberOfEventsPerRun, $filename);
    print "$macrofile\n";

    unlink($filename);
    system "condor_submit", "$condorfile";
    ++$currentRun;
  }
}

sub make_macro_file{
  # make a macro file, first argument is run number, second is current value (to be placed in macro template below)

  my ($currentRun, $momentum, $angle, $nEvents, $filename) = @_;

  my $rotation = $angle/2.;

  my $momentumString = sprintf("%.1f", $momentum);
  my $angleString = sprintf("%.2f", $angle);

  open MACROFILE, ">${condor_dir}/mac/res_${currentRun}.mac" or die "Error: Cannot make macro file: $!";

  print MACROFILE <<EOF;
/control/verbose 1
/run/verbose 0
/tracking/verbose 0
/RES/Fit/Verbose 0

/RES/Gun/RandomOrigin
/RES/Gun/RandomDirection
/gun/particle e-
/gun/momentumAmp ${momentum} GeV

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

/RES/Field/SetInhomFieldFrom ${table_dir}/perdaix_07_jul_2009.table
#/RES/Field/SetDummyField 0.27 0.0 0.0 tesla
#/RES/Field/SetUniformField 0.3 0.0 0.0 tesla

/RES/Fit/Method blobel

/RES/Data/OverWriteFile true
/RES/Data/SetFileName ${filename}
/RES/Run/StoreResults

/run/initialize

/process/activate msc
/RES/Run/Generate ${nEvents}
/process/inactivate msc
/RES/Run/Reconstruct
EOF

  close MACROFILE;

  # return value: macro file name
  "${condor_dir}/mac/res_${currentRun}.mac";

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

