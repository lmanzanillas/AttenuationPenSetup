#ifndef SILICONPLATECONSTUCTION_H
#define SILICONPLATECONSTUCTION_H

#include "globals.hh"
#include "G4VUserDetectorConstruction.hh"
#include "DetectorMessenger.hh"
#include "G4SystemOfUnits.hh"


class SiliconPlateConstruction
{
  public:
    G4VSolid* ConstructPlate();

  private:
    G4double fSiliconPlate_h;
    G4double fHolderWidth;
};
#endif
