#include "G4UnionSolid.hh"
#include "G4SubtractionSolid.hh"
#include "G4IntersectionSolid.hh"
#include "G4UnionSolid.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"

#include "G4Material.hh"
#include "G4Element.hh"
#include "G4LogicalBorderSurface.hh"
#include "G4LogicalSkinSurface.hh"
#include "G4OpticalSurface.hh"

#include "G4SystemOfUnits.hh"
#include "G4UnitsTable.hh"
#include "G4NistManager.hh"

#include "LightGuideConstruction.hh"

// Class to define logical / virtual volume for the light guide

G4VSolid* LightGuideConstruction::ConstructPlate(){

  //Hole dimensions for grease hole
  TubsStartAngle =  0;
  TubsSpanningAngle = 360 * deg;
  diameterGreaseHole = 25.5*mm;
  depthGreaseHole = 1.*mm;

  //Initial block dimensions
  LighGuideSizeX = 20.*mm;
  LighGuideSizeY = 14.*mm;
  LighGuideSizeZ = 14.*mm;
  chamderSize = 2.25*mm;
  slitToPlaceThinSintillatorHeight = 0.125*mm;
  slitToPlaceThinSintillatorOffset = 1.75*mm;
  slitToPlaceThinSintillatorDepth = 0.5*mm;


  G4RotationMatrix* rm = new G4RotationMatrix();
  rm->rotateX(45.*deg);
  G4RotationMatrix* rm1 = new G4RotationMatrix();
  rm1->rotateZ(-35.*deg);
  G4RotationMatrix* rm2 = new G4RotationMatrix();
  rm2->rotateZ(-35.*deg);
  rm2->rotateX(45.*deg);
  G4RotationMatrix* rm3 = new G4RotationMatrix();
  rm3->rotateY(90.*deg);

  //char solidname[100];

  G4VSolid* initialBlock = new G4Box("lg_block", LighGuideSizeX, LighGuideSizeY, LighGuideSizeZ);
  G4VSolid* angleBlock = new G4Box("angleBlock", 35.*mm, 20.*mm, 15.*mm);

  G4VSolid* chamferEdge = new G4Box("chamfer", LighGuideSizeX+1.*mm, chamderSize, chamderSize);
  G4VSolid* chamferTwo = new G4Box("chamferTwo", 35.*mm, chamderSize+1.*mm, chamderSize+1.*mm);

  G4VSolid* boxSlitToPlaceThinSintillator = new G4Box("boxSlitToPlaceThinSintillator", slitToPlaceThinSintillatorDepth, slitToPlaceThinSintillatorHeight, LighGuideSizeY);

  G4VSolid* holeForOpticalGrease = new G4Tubs("dip", 0, diameterGreaseHole/2., depthGreaseHole/2., TubsStartAngle, TubsSpanningAngle);

  G4SubtractionSolid* oneEdge = new G4SubtractionSolid("oneEdge", initialBlock, chamferEdge, rm, G4ThreeVector(0*mm,LighGuideSizeY, LighGuideSizeY));
  oneEdge = new G4SubtractionSolid("oneEdge", oneEdge, chamferEdge, rm, G4ThreeVector(0*mm,-LighGuideSizeY, LighGuideSizeY));
  oneEdge = new G4SubtractionSolid("oneEdge", oneEdge, chamferEdge, rm, G4ThreeVector(0*mm,-LighGuideSizeY, -LighGuideSizeY));
  oneEdge = new G4SubtractionSolid("oneEdge", oneEdge, chamferEdge, rm, G4ThreeVector(0*mm, LighGuideSizeY, -LighGuideSizeY));

  G4SubtractionSolid* twoEdge = new G4SubtractionSolid("twoEdge", oneEdge,angleBlock,rm1, G4ThreeVector(LighGuideSizeX,-LighGuideSizeY, 0));

  G4SubtractionSolid* threeEdge = new G4SubtractionSolid("threeEdge", twoEdge, chamferTwo, rm2, G4ThreeVector(+chamderSize, -chamderSize, chamderSize+LighGuideSizeY));
  threeEdge = new G4SubtractionSolid("threeEdge", threeEdge, chamferTwo, rm2, G4ThreeVector(+chamderSize, -chamderSize, -chamderSize-LighGuideSizeY));
  threeEdge = new G4SubtractionSolid("threeEdge", threeEdge, holeForOpticalGrease, rm3, G4ThreeVector(-LighGuideSizeX+0.5*mm, 0, 0));

  G4SubtractionSolid* guide = new G4SubtractionSolid("final_guide", threeEdge, boxSlitToPlaceThinSintillator, 0, G4ThreeVector(LighGuideSizeX, LighGuideSizeY-slitToPlaceThinSintillatorOffset, 0));

  G4VSolid* finalPlate = guide;
  return finalPlate;

}

// G4LogicalVolume* LightGuideConstruction::ConstructGuideLog(){
//   G4VSolid* light_guide = ConstructPlate();
//
//   G4NistManager* man = G4NistManager::Instance();
//
//   G4Material* PMMA = man->FindOrBuildMaterial("G4_PLEXIGLASS");
//
//   G4double refractive_index[] = {1.49, 1.49, 1.49, 1.49, 1.49};
//   G4double abs[] = {0.5*m, 0.5*m, 0.5*m, 0.5*m, 0.5*m};
//   G4double refl[] = {0.9, 0.9, 0.9, 0.9, 0.9};
//   G4double energy[] = {2.48*eV, 2.58*eV, 2.68*eV, 2.78*eV, 3.1*eV};
//
//   G4MaterialPropertiesTable* pmmaMPT = new G4MaterialPropertiesTable();
//   pmmaMPT->AddProperty("RINDEX", energy, refractive_index, 5);
//   pmmaMPT->AddProperty("ABSLENGTH", energy, abs, 5);
//   pmmaMPT->AddProperty("REFLECTIVITY", energy, refl, 5);
//   PMMA->SetMaterialPropertiesTable(pmmaMPT);
//
//   G4LogicalVolume* guide_log = new G4LogicalVolume(light_guide, PMMA, "Ligh_guide_log");
//   return guide_log;
// }
