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

    G4double GetLighGuideSizeX()  {return LighGuideSizeX;};
    G4double GetLighGuideSizeY()  {return LighGuideSizeY;};
    G4double GetLighGuideSizeZ()  {return LighGuideSizeZ;};
    G4double GetDepthGreaseHole()  {return depthGreaseHole;};
    G4double GetDiameterGreaseHole()  {return diameterGreaseHole;};
    G4double GetSlitToPlaceThinSintillatorOffset()  {return slitToPlaceThinSintillatorOffset;};
    G4double GetSlitToPlaceThinSintillatorDepth()  {return slitToPlaceThinSintillatorDepth;};
    G4double GetSlitToPlaceThinSintillatorHeight()  {return slitToPlaceThinSintillatorHeight;};

  private:
    //Hole dimensions for grease hole
    G4double TubsStartAngle ;
    G4double TubsSpanningAngle ;
    G4double diameterGreaseHole;
    G4double depthGreaseHole ;

    //Initial block dimensions
    G4double LighGuideSizeX ;
    G4double LighGuideSizeY ;
    G4double LighGuideSizeZ ;
    G4double chamderSize ;
    G4double slitToPlaceThinSintillatorHeight ;
    G4double slitToPlaceThinSintillatorOffset ;
    G4double slitToPlaceThinSintillatorDepth ;


    G4double fSiliconPlate_h;
    G4double fHolderWidth;
};
#endif
