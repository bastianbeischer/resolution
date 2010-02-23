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

my $minMomentum = 50.0;
my $maxMomentum = 50.0;
my $momentumStep = 0.1;
my $minAngle = 0.1;
my $maxAngle = 10.0;
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
    my $filename = "${result_dir}/pebs01_${momentumString}_GeV_${angleString}_deg_inhom_msc.root";

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

/RES/Gun/StartZ 120 cm
/RES/Gun/RandomOrigin
/RES/Gun/RandomDirection
/gun/energy ${momentum} GeV
/gun/particle e-

/control/alias moduleRot -${rotation}
/control/alias moduleInternalRot ${angle}
/control/alias fiberEfficiency 1.0
/control/alias fiberResolution 35

# panel 1_1
/RES/Det/AddModule 0. 0. 116.0 cm
/RES/Det/SetModuleRotation 0 {moduleRot}
/RES/Det/SetModuleInternalRotation 0 {moduleInternalRot}
/RES/Det/SetModuleLength 0 152.0
/RES/Det/SetModuleWidth 0 152.0
/RES/Det/SetModuleUpperSigmaV 0 {fiberResolution}
/RES/Det/SetModuleLowerSigmaV 0 {fiberResolution}
/RES/Det/SetModuleUpperEfficiency 0 {fiberEfficiency}
/RES/Det/SetModuleLowerEfficiency 0 {fiberEfficiency}

# panel 1_2
/RES/Det/AddModule 0. 0. 114.0 cm
/RES/Det/SetModuleRotation 1 {moduleRot}
/RES/Det/SetModuleInternalRotation 1 {moduleInternalRot}
/RES/Det/SetModuleLength 1 152.0
/RES/Det/SetModuleWidth 1 152.0
/RES/Det/SetModuleUpperSigmaV 1 {fiberResolution}
/RES/Det/SetModuleLowerSigmaV 1 {fiberResolution}
/RES/Det/SetModuleUpperEfficiency 1 {fiberEfficiency}
/RES/Det/SetModuleLowerEfficiency 1 {fiberEfficiency}

# panel 2_1
/RES/Det/AddModule 0. 0. 28.0 cm
/RES/Det/SetModuleRotation 2 {moduleRot}
/RES/Det/SetModuleInternalRotation 2 {moduleInternalRot}
/RES/Det/SetModuleLength 2 86.0
/RES/Det/SetModuleWidth 2 86.0
/RES/Det/SetModuleUpperSigmaV 2 {fiberResolution}
/RES/Det/SetModuleLowerSigmaV 2 {fiberResolution}
/RES/Det/SetModuleUpperEfficiency 2 {fiberEfficiency}
/RES/Det/SetModuleLowerEfficiency 2 {fiberEfficiency}

# panel 2_2
/RES/Det/AddModule 0. 0. 26.0 cm
/RES/Det/SetModuleRotation 3 {moduleRot}
/RES/Det/SetModuleInternalRotation 3 {moduleInternalRot}
/RES/Det/SetModuleLength 3 86.0
/RES/Det/SetModuleWidth 3 86.0
/RES/Det/SetModuleUpperSigmaV 3 {fiberResolution}
/RES/Det/SetModuleLowerSigmaV 3 {fiberResolution}
/RES/Det/SetModuleUpperEfficiency 3 {fiberEfficiency}
/RES/Det/SetModuleLowerEfficiency 3 {fiberEfficiency}

# panel 3_1
/RES/Det/AddModule 0. 0. -25.0 cm
/RES/Det/SetModuleRotation 4 {moduleRot}
/RES/Det/SetModuleInternalRotation 4 {moduleInternalRot}
/RES/Det/SetModuleLength 4 86.0
/RES/Det/SetModuleWidth 4 86.0
/RES/Det/SetModuleUpperSigmaV 4 {fiberResolution}
/RES/Det/SetModuleLowerSigmaV 4 {fiberResolution}
/RES/Det/SetModuleUpperEfficiency 4 {fiberEfficiency}
/RES/Det/SetModuleLowerEfficiency 4 {fiberEfficiency}

# panel 3_2
/RES/Det/AddModule 0. 0. -27.0 cm
/RES/Det/SetModuleRotation 5 {moduleRot}
/RES/Det/SetModuleInternalRotation 5 {moduleInternalRot}
/RES/Det/SetModuleLength 5 86.0
/RES/Det/SetModuleWidth 5 86.0
/RES/Det/SetModuleUpperSigmaV 5 {fiberResolution}
/RES/Det/SetModuleLowerSigmaV 5 {fiberResolution}
/RES/Det/SetModuleUpperEfficiency 5 {fiberEfficiency}
/RES/Det/SetModuleLowerEfficiency 5 {fiberEfficiency}

# panel 4_1
/RES/Det/AddModule 0. 0. -70.0 cm
/RES/Det/SetModuleRotation 6 {moduleRot}
/RES/Det/SetModuleInternalRotation 6 {moduleInternalRot}
/RES/Det/SetModuleLength 6 86.0
/RES/Det/SetModuleWidth 6 86.0
/RES/Det/SetModuleUpperSigmaV 6 {fiberResolution}
/RES/Det/SetModuleLowerSigmaV 6 {fiberResolution}
/RES/Det/SetModuleUpperEfficiency 6 {fiberEfficiency}
/RES/Det/SetModuleLowerEfficiency 6 {fiberEfficiency}

# panel 4_2
/RES/Det/AddModule 0. 0. -72.0 cm
/RES/Det/SetModuleRotation 7 {moduleRot}
/RES/Det/SetModuleInternalRotation 7 {moduleInternalRot}
/RES/Det/SetModuleLength 7 86.0
/RES/Det/SetModuleWidth 7 86.0
/RES/Det/SetModuleUpperSigmaV 7 {fiberResolution}
/RES/Det/SetModuleLowerSigmaV 7 {fiberResolution}
/RES/Det/SetModuleUpperEfficiency 7 {fiberEfficiency}
/RES/Det/SetModuleLowerEfficiency 7 {fiberEfficiency}

/RES/Field/SetInhomFieldFrom ${table_dir}/pebs01_25_jan_2010.table
#/RES/Field/SetDummyField 0.27 0.0 0.0 tesla
#/RES/Field/SetUniformField 0.3 0.0 0.0 tesla

/RES/Data/OverWriteFile true
/RES/Data/SetFileName ${filename}
/RES/Run/StoreResults

/run/initialize

/process/inactivate eIoni
/process/inactivate eBrem

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

