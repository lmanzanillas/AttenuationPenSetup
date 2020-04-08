# MaterialOpt
GEANT4 Simulation for optimizing thickness of shielding.

Simulates various sources of radiation incident on a target. Target thickness and material can be controlled via macro, as can the source type and energy. Includes geometry for detection of scintillation photons in a PMT coupled to the base of the tile, or two SiPMs coupled to the side of the tile.

The following histograms are produced:

Light Output - number of photons detected at the chosen detector type.
Energy deposited in target, in MeV.
Light Yield - number of photons produced in an event.
Ratio of Light Output to Light Yield.
Number of photons leaving the target. Geometry independent version of light output.

To use: Navigate to build folder and enter:

cmake -DGeant4_DIR=/opt/geant4/lib64/Geant4-10.3.0 ..

This creates the makefile as needed. Then, run make to create the simulation. ./PEN runs the program.

Additional macro commands and cases:

/PEN/det/setTargetMat [material]

Set the target material to supplied material. Material must be loaded in Geant4 prior to calling. Default to air.

Materials provided:

G4_TEFLON
G4_Al
G4_Si
G4_Cu
Water
G4_Pyrex_Glass
PEN
Galactic - vacuum

/PEN/det/setWorldMat [material]

Set the material for the world volume. Default to air.

Tested cases:

air
Galactic - vacuum

/PEN/det/setSize [number with unit]

Set target thickness to supplied value.

NB: The simulation generates the name of the file automatically using the BestUnit function. Using a thickness that is converted to a decimal value (15 mm becomes 1.5 cm) causes issues with file naming.

/PEN/det/setDetectorType [int 0-1]

Set detector type. Default 0.

  0 - PMT coupled to base of target.
  1 - 2 SiPMs, coupled to the side of the target.

/PEN/gun/sourceType [int 0-6]

Choose the type of source used in the simulation. Allowed cases:

  0 - perpendicular mono-energetic gamma with fixed position - use with /PEN/gun/sourceEnergy
  1 - 60Co source
  2 - 137Cs source
  3 - 90Sr source
	4 - 241Am source
	5 - 106Ru source

/PEN/gun/sourceEnergy [number with unit]

Set energy for the mono-energetic gamma. Automatically sets source type to 0.
