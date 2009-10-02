#!/usr/bin/perl -w

use strict;

#
# EDIT firstRunNumber
# EDIT enter VALUEs
# EDIT macro template to place values
# RUN script
#

my $numberOfEventsPerRun = 10000;

my $minEnergy = 1;
my $maxEnergy = 10; # GeV
my $minAngle = 5.0;
my $maxAngle = 5.0;

my $executable = "resolution";

################################################################################################################



for (my $angle = $minAngle; $angle <= $maxAngle; ++$angle) {
    for (my $energy = $minEnergy; $energy <= $maxEnergy; ++$energy) {
        my $macrofile = &make_macro_file($energy, $angle, $numberOfEventsPerRun);
        print "$macrofile\n";

        unlink("/home/home4/institut_1b/beischer/src/geant4/resolution/results/perdaix_${energy}_GeV_${angle}_deg.root");
        #system "resolution", "";
    }
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

#/RES/Gun/RandomOrigin
/RES/Gun/RandomDirection
/gun/energy ${energy} GeV
#/gun/position 0 0 50 cm
#/gun/direction 0 -2 -400

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
/RES/Data/SetFileName results/perdaix_${energy}_GeV_${angle}_deg.root
/RES/Run/StoreResults

/run/initialize

/RES/Run/Generate ${nEvents}
/RES/Run/Reconstruct

EOF

  close MACROFILE;

  # return value: macro file name
  "perdaix_${energy}_GeV_${angle}_deg.mac";

}
