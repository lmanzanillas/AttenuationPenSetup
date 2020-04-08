#ifndef LIGHTGUIDECONSTUCTION_H
#define LIGHTGUIDECONSTUCTION_H

#include "globals.hh"
#include "G4VUserDetectorConstruction.hh"
#include "DetectorMessenger.hh"
#include "G4SystemOfUnits.hh"
#include "G4LogicalVolume.hh"

class LightGuideConstruction
{
  public:
    G4VSolid* ConstructPlate();
    G4LogicalVolume* ConstructGuideLog();

  private:
    G4double fSiliconPlate_h;
    G4double fHolderWidth;
};
#endif
