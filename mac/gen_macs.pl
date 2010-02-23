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
my $maxEnergy = 1; # GeV
my $energyStep = 1;
my $minAngle = 0.1;
my $maxAngle = 1.0;
my $angleStep = 0.1;

my $executable = "resolution";

################################################################################################################



for (my $angle = $minAngle; $angle <= $maxAngle; $angle += $angleStep) {
    for (my $energy = $minEnergy; $energy <= $maxEnergy; $energy += $energyStep) {
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

/RES/Gun/StartZ 120 cm
/RES/Gun/RandomOrigin
/RES/Gun/RandomDirection
/gun/energy ${energy} GeV
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

/RES/Field/SetInhomFieldFrom tables/pebs01_25_jan_2010.table
#/RES/Field/SetDummyField 0.27 0.0 0.0 tesla
#/RES/Field/SetUniformField 0.3 0.0 0.0 tesla

/RES/Data/OverWriteFile true
/RES/Data/SetFileName results/pebs01_${energy}_GeV_${angle}_deg_inhom_msc.root
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
  "perdaix_${energy}_GeV_${angle}_deg.mac";

}
