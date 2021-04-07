//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
/// \file OpNovice/include/OpNoviceDetectorConstruction.hh
/// \brief Definition of the OpNoviceDetectorConstruction class
//
//
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef DetectorConstruction_h
#define DetectorConstruction_h 1

#include "globals.hh"
#include "G4VUserDetectorConstruction.hh"
#include "DetectorMessenger.hh"
#include "G4SystemOfUnits.hh"
#include "G4MaterialPropertiesTable.hh"
#include "PenMaterials.hh"
#include "G4GDMLParser.hh"
class DetectorMessenger;
class G4LogicalVolume;
class G4Material;
class PenMaterials;
class G4Box;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class DetectorConstruction : public G4VUserDetectorConstruction
{
  public:
    DetectorConstruction();
    virtual ~DetectorConstruction();
    static G4VPhysicalVolume *
    GetPhysicalVolumeByName(const G4String &name);

  public:
    virtual G4VPhysicalVolume* Construct();
    //static G4VPhysicalVolume* GetPhysicalVolumeByName(const G4String &name);
    void SetSize  (G4double);
    void SetTargetMaterial(G4String);
    void SetWorldMaterial(G4String);
    void SetPMTPlacement(G4ThreeVector);
    virtual G4VPhysicalVolume* UpdateGeometry();
    void SetPropertyTable(G4Material* mat, G4MaterialPropertiesTable* tab);
    void SetDetectorType(G4int);
    void SetTargetSampleLength(G4double);
    void SetTargetSampleThickness(G4double);
    void SetTargetSampleWidth(G4double);
    void SetDetectorCollimatorX(G4double);
    void SetDetectorCollimatorThickness(G4double);
    void SetNumberOfTargetSamples(G4int);
    void SetLY(G4double);
    void SetRes(G4double);
    void SetABS(G4double);
    void SetSigAlpha(G4double);
    void SetSigAlphaSides(G4double);
    void SetSigAlphaBottom(G4double);
    void SetPMTReflectivitySides(G4double);
    void SetPMTReflectivityBottom(G4double);
    void MaterialPropertiesScintillator();
    void SetVolName(G4ThreeVector);
    void SetRI(G4double);
    void SetReflectorOn(G4bool );     

    void SetDetectorName(G4String);
    void SetABSFile(G4String);



  public:
    //const G4VPhysicalVolume* GetWorld() {return physicPenSampleBox;};
    const G4VPhysicalVolume* GetWorld() {return physicWorldBox;};
    G4int GetDetectorType(){return fDetectorType;};
    G4double GetDetectorCollimatorX(){return fDetectorCollimatorX;};
    G4double GetDetectorCollimatorY(){return fDetectorCollimatorY;};
    G4double GetSourceContainerY(){return fSourceContainerY;};
    G4double GetLY(){return fLY;};
    G4double GetRes(){return fRES;};
    G4double GetABS(){return AbsorptionLength;};
    G4double GetSigAlpha(){return fSigAlphaPENAirTopBottom;};
    G4double GetSigAlphaSides(){return fSigAlphaSides;};
    G4double GetSigAlphaBottom(){return fSigAlphaBottom;};
    G4double GetPMTReflectivitySides(){return pmtReflectivitySides;};
    G4double GetPMTReflectivityBottom(){return pmtReflectivityBottom;};
    G4String GetDetectorName(){return fDetectorName;};
    G4String GetVolName(){return fVolName;};
    G4double GetRI(){return fRI;};

    G4double GetCollimatorThickness()  {return halfCollimatorThickness;};
    G4double GetTargetSampleLength()  {return halfPenSampleLength;};
    G4double GetTargetSampleThickness()  {return halfPenSampleThickness;};
    G4double GetTargetSampleWidth()  {return halfPenSampleWidth;};
    G4int GetNumberOfTargetSamples()  {return nSamples;};
    G4String GetTargetMaterialName(){return fTargetName;};
    G4Material*        GetWorldMaterial()   {return fWorldMaterial;};
    G4Material*        GetTargetMaterial()   {return fTargetMaterial;};
    G4ThreeVector* fSourceVector;
    G4bool GetReflectorOn(){return reflectorOn;}

    void               DefineMaterials();

  private:

    G4double halfSizeDarkBoxX;
    G4double halfSizeDarkBoxY;
    G4double halfSizeDarkBoxZ;
    G4double SpaceBetweenSamples;
    G4double halfGreaseThickness;

    G4double fLY;
    G4double fRES;
    G4double AbsorptionLength;
    G4double fSigAlphaPENAirTopBottom;
    G4double fSigAlphaSides;
    G4double fSigAlphaBottom;
    G4double pmtReflectivitySides;
    G4double pmtReflectivityBottom;
    G4double fRI;
    //add position
    G4double halfPenSampleLength;
    G4double halfPenSampleThickness;
    G4double halfPenSampleWidth;
    G4int nSamples;

    G4double userActivePhotoCathodeLength;
    G4double userActivePhotoCathodeWidth;
    G4double activePhotoCathodePMTLength;
    G4double activePhotoCathodePMTWidth;

    G4double fSiliconPlate_h;
    G4double fHolderWidth;

    G4Material* fWorldMaterial;
    G4Material* fTargetMaterial;
    G4Material* fGlassMaterialPMT;
    G4Material* fPhotoCathodeMaterial;
    G4String fTargetName;
    G4double halfSourceContainerThickness;
    G4double halfCollimatorThickness;
    G4Box* fWorldBox;
    G4Box* penSampleBox;

    G4LogicalVolume* penLogicBox;
    G4LogicalVolume* logicWorldBox;

    G4VSolid* boxGreaseBetweenSamples;
    G4VSolid* boxGreaseBetweenSamplesAndPMT;
    G4LogicalVolume* logicGreasePEN;
    G4LogicalVolume* logicGreasePENPMT;
    
    G4ThreeVector fPMTPlacement = G4ThreeVector(0,0,63.5*mm+halfPenSampleThickness);

    G4VPhysicalVolume* physicWorldBox;
    G4VPhysicalVolume* physicPenSampleBox;
    G4VPhysicalVolume* physicPenGrease;
    G4VPhysicalVolume* physicSourceContainer;
    G4VPhysicalVolume* physicCollimator;
    G4VPhysicalVolume* physicBoxPMTShell;
    G4VPhysicalVolume* physicBoxPhotoCathodeSupport;
    G4VPhysicalVolume* physicActivePhotoCathodeTriggerPMT;
    G4VPhysicalVolume* physicBoxEmptyInsidePMT;

    G4VPhysicalVolume* physicActivePhotoCathodePMT1;
    G4VPhysicalVolume* physicActivePhotoCathodePMT2;
    G4VPhysicalVolume* physicActivePhotoCathodePMT3;
    G4VPhysicalVolume* physicActivePhotoCathodePMT4;
    G4VPhysicalVolume* physicActivePhotoCathodePMT5;
    G4VPhysicalVolume* physicPhotoCathodeSupportPMT1;
    G4VPhysicalVolume* physicPhotoCathodeSupportPMT2;
    G4VPhysicalVolume* physicPhotoCathodeSupportPMT3;
    G4VPhysicalVolume* physicPhotoCathodeSupportPMT4;
    G4VPhysicalVolume* physicPhotoCathodeSupportPMT5;

    G4VPhysicalVolume* physicReflectorFoilBoxOverPEN;
    G4VPhysicalVolume* vacPlacement;

    G4VPhysicalVolume* physicPenStackedSamples;

    G4VPhysicalVolume* physicTriggerFoilEJ212;
    G4VPhysicalVolume* physicReflectorFoilAroundEJ212Foil;
    G4VPhysicalVolume* physicLightGuide;
    G4VPhysicalVolume* physicOpticalGrease;

    G4bool reflectorOn;


    G4OpticalSurface* surfaceAirTargetTopBottom;
    G4OpticalSurface* surfaceAirTargetLongSide;
    G4OpticalSurface* surfaceAirTargetShortSide;
    G4OpticalSurface* surfaceCathodeSupport;
    G4OpticalSurface* surfaceCathodeSupportBottom;
    G4OpticalSurface* surfaceGreaseTargetSides;
    G4OpticalSurface* surfaceGreaseTargetBottom;
    G4OpticalSurface* surfaceGreaseTargetMiddle;

    G4MaterialPropertiesTable* MPT_PEN;
    G4MaterialPropertiesTable* MPT_GlassPMT;
    G4MaterialPropertiesTable* MPT_Target;
    G4MaterialPropertiesTable* SMPT_surfaceAirTargetTopBottom;
    G4MaterialPropertiesTable* MPT_World;
    G4MaterialPropertiesTable* MPT_SurfaceBottom;
    G4MaterialPropertiesTable* MPT_SurfaceBetween;
    G4MaterialPropertiesTable* MPT_SurfaceSides;

    G4LogicalBorderSurface* logicSurfaceAirTarget;
    G4LogicalBorderSurface* logicSurfaceTargetAir;
    G4LogicalBorderSurface* logicSurfaceAirTargetStacked;
    G4LogicalBorderSurface* logicSurfacePENStackedAir;
    G4LogicalBorderSurface* logicSurfacePENStackedPENStacked;
    G4LogicalBorderSurface* logicSurfacePENStackedPENStackedCathodeSupport;
    G4LogicalBorderSurface* logicSurfaceEJ212PMMA;
    G4LogicalBorderSurface* logicSurfacePMMAEJ212;
    G4LogicalBorderSurface* logicSurfaceAirPMMA;
    G4LogicalBorderSurface* logicSurfacePMMAAir;
    G4LogicalBorderSurface* logicSurfacePENStackedReflectorFoilBoxOverPEN;
    G4LogicalBorderSurface* logicSurfaceReflectorFoilBoxOverPENPENStacked;
    G4LogicalBorderSurface* logicSurfaceEJ212Reflector;
    G4LogicalBorderSurface* logicSurfaceReflectorEJ212;
    G4LogicalBorderSurface* logicSurfacePenReflectorFoilBoxOverPEN;

    G4Material* PenMaterial;
    G4Material* materialBialkali;
    G4Material* materialSi;
    G4Material* fGe;
    G4Material* materialAir;
    G4Material* fVacuum;
    G4Material* materialTriggerFoilEJ212;
    G4Material* Pstyrene;
    G4Material* fGlass;
    G4Material* fPOM;
    G4Material* fABS;
    G4Material* materialPMMA;
    G4Material* materialGreaseEJ550;
    G4Material* materialTeflon;
    G4Material* materialVikuiti;
    G4Material* materialPolyethylene;
    G4Material* materialTitanium;

    G4int fDetectorType;
    G4double fDetectorCollimatorX;
    G4double fDetectorCollimatorY;
    G4double fSourceContainerY;
    G4String fDetectorName;
    G4String fVolName;
    G4String fABSFile;

    PenMaterials* materialConstruction;
    DetectorMessenger* fDetectorMessenger;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif /*DetectorConstruction_h*/
