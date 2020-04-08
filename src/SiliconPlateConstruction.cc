#include "G4UnionSolid.hh"
#include "G4SubtractionSolid.hh"
#include "G4IntersectionSolid.hh"
#include "G4UnionSolid.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "SiliconPlateConstruction.hh"

G4VSolid* SiliconPlateConstruction::ConstructPlate(){
  fSiliconPlate_h = 1.5*mm;
  fHolderWidth=90.00*mm;

  G4double TubsStartAngle =  0;
  G4double TubsSpanningAngle = 360 * deg;

  G4double trg_b = (fHolderWidth/sqrt(2.))*mm;
  G4double trg_h = 45.*mm;
  G4double rect_x = 2.*(trg_h + 9.5*mm);
  G4double rect_y = fHolderWidth;
  G4double rect_bc_x = 28.*mm;
  G4double cut_x = 20.*mm;
  G4double cut_y = 31.*mm;

  G4RotationMatrix* rm = new G4RotationMatrix();
  rm->rotateZ(45.*deg);
  G4RotationMatrix* rm1 = new G4RotationMatrix();
  rm1->rotateZ(-120.*deg);
  G4RotationMatrix* rm2 = new G4RotationMatrix();
  rm2->rotateZ(120.*deg);

  char solidname[100];

  G4VSolid* solid_LowerPlate_Rect  = new G4Box(solidname,
                                      rect_x/2,
                                      rect_y/2,
                                      0.5*fSiliconPlate_h);
  // build rectangular cuts
  G4VSolid* solid_LowerPlate_CutRect0  = new G4Box(solidname,
                                         (0.5*rect_x- rect_bc_x+1*mm)*mm/2,
                                         (rect_y+1*mm)/2,
                                         0.5*(fSiliconPlate_h+.1*mm));
  G4VSolid* solid_LowerPlate_CutRect1  = new G4Box(solidname,
                                                   trg_b/2,
                                                   trg_b/2,
                                                   0.5*(fSiliconPlate_h+.1*mm));

  G4VSolid* solid_LowerPlate_CutRect2  = new G4Box(solidname,
                                                   cut_x/2,
                                                   cut_y/2,
                                                   0.5*(fSiliconPlate_h+.1*mm));

  //build rectangular cut plates (to imitate rounded corners, important to fit into MS)
  G4VSolid* solid_LowerPlate_CutRect3  = new G4Box(solidname,
                                                    9*mm/2,
                                                    20*mm/2,
                                                    0.5*(fSiliconPlate_h+.1*mm));

  G4VSolid* solid_LowerPlate_CutRect4  = new G4Box(solidname,
                                                   6*mm/2,
                                                   20*mm/2,
                                                   0.5*(fSiliconPlate_h+.1*mm));

  // cut final silicon plate

  G4VSolid* solid_LowerPlate_Base0 = new G4SubtractionSolid(solidname,
                                                            solid_LowerPlate_Rect,
                                                            solid_LowerPlate_CutRect0,
                                                            0,
                                                            G4ThreeVector((-.5*rect_bc_x-.5*.5*rect_x-.5*1*mm),0.,0.));


  G4VSolid* solid_LowerPlate_Base1 = new G4SubtractionSolid(solidname,
                                                            solid_LowerPlate_Base0,
                                                            solid_LowerPlate_CutRect1,
                                                            rm,
                                                            G4ThreeVector(rect_x/2,trg_h,0.));


  G4VSolid* solid_LowerPlate_Base2 = new G4SubtractionSolid(solidname,
                                                            solid_LowerPlate_Base1,
                                                            solid_LowerPlate_CutRect1,
                                                            rm,
                                                            G4ThreeVector(rect_x/2,-trg_h,0.));


  G4VSolid* solid_LowerPlate_Base3 = new G4SubtractionSolid(solidname,
                                                            solid_LowerPlate_Base2,
                                                            solid_LowerPlate_CutRect2,
                                                            0,
                                                            G4ThreeVector(0.,(0.5*cut_y-3.5*mm),0.));


  G4VSolid* solid_LowerPlate_Base4 = new G4SubtractionSolid(solidname,
                                                            solid_LowerPlate_Base3,
                                                            solid_LowerPlate_CutRect3,
                                                            0,
                                                            G4ThreeVector(rect_x/2,0,0.));


  G4VSolid* solid_LowerPlate_Base5 = new G4SubtractionSolid(solidname,
                                                            solid_LowerPlate_Base4,
                                                            solid_LowerPlate_CutRect4,
                                                            rm1,
                                                            G4ThreeVector(-rect_bc_x,trg_h,0.));


  G4VSolid* solid_LowerPlate_Base = new G4SubtractionSolid(solidname,
                                                            solid_LowerPlate_Base5,
                                                            solid_LowerPlate_CutRect4,
                                                            rm2,
                                                            G4ThreeVector(-rect_bc_x,-trg_h,0.));


  G4Tubs* solid_LowerPlate_hole  = new G4Tubs(solidname,
                0*mm,
                (1.5*mm+0.1*mm),
                0.5*fSiliconPlate_h,
                TubsStartAngle, TubsSpanningAngle);

  //if(fCrystalRadius>41.*mm){fMaximalAllowedCrystalRadius=46.*mm;}
  //radius condition identifies GTFs
  G4double tmpR = 46.*mm + 1*mm + 1.5*mm;
  G4VSolid* current_solid_LowerPlate = solid_LowerPlate_Base;

  for (int i = 0; i < 3; i++) {
    const G4double VertBarAngle = ((G4double) i) * 120.*deg;
    const G4ThreeVector VertBarTranslation (/* x */ tmpR * std::cos(VertBarAngle),  /* y */ tmpR * std::sin(VertBarAngle), /* z */ 0.);

    G4VSolid* solid_LowerPlate_nHoles = new G4SubtractionSolid(solidname,
                     current_solid_LowerPlate,
                     solid_LowerPlate_hole,
                     0,
                     VertBarTranslation);
    current_solid_LowerPlate = solid_LowerPlate_nHoles;
  }
  G4VSolid* final_plate = current_solid_LowerPlate;
  return final_plate;

}
