#include "DetectorConstruction.hh"
#include "DetectorMessenger.hh"
#include "SiliconPlateConstruction.hh"
#include "LightGuideConstruction.hh"
#include "PrimaryGeneratorAction.hh"

#include "G4Material.hh"
#include "G4Element.hh"
#include "G4LogicalBorderSurface.hh"
#include "G4LogicalSkinSurface.hh"
#include "G4OpticalSurface.hh"

#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4Sphere.hh"
#include "G4Cons.hh"
#include "G4Torus.hh"
#include "G4Hype.hh"

#include "G4Transform3D.hh"
#include "G4LogicalVolume.hh"
#include "G4ThreeVector.hh"
#include "G4PVPlacement.hh"
#include "G4SystemOfUnits.hh"
#include "G4UnitsTable.hh"
#include "G4NistManager.hh"

#include "G4MultiFunctionalDetector.hh"
#include "G4SDManager.hh"
#include "G4PSEnergyDeposit.hh"
#include <G4VPrimitiveScorer.hh>

#include "G4UnionSolid.hh"
#include "G4SubtractionSolid.hh"
#include "G4IntersectionSolid.hh"
#include "G4UnionSolid.hh"
#include "G4VoxelLimits.hh"

#include "G4MTRunManager.hh"
#include "G4PhysicalConstants.hh"

#include "G4GeometryManager.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4SolidStore.hh"

#include "G4Navigator.hh"
#include "G4TransportationManager.hh"

#include "G4GDMLParser.hh"

#include <G4VisAttributes.hh>
#include <iostream>
#include <fstream>

using namespace std;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

/*
Constructs DetectorConstruction, defines default values.
*/

DetectorConstruction::DetectorConstruction()
: G4VUserDetectorConstruction(),  penSampleBox(nullptr), penLogicBox(nullptr), 
physicPenSampleBox(nullptr),physicCollimator(nullptr),
AirPEN(nullptr),	
MPT_PEN(nullptr)
{
  fDetectorMessenger = new DetectorMessenger(this);
  //MPT_PEN = new G4MaterialPropertiesTable();
  halfSizeDarkBoxX = halfSizeDarkBoxY = halfSizeDarkBoxZ = 1.*m;
  fTargetName = "holder";
  reflectorOn = false;
  halfPenSampleLength = 37.*mm;
  halfPenSampleThickness = 0.85*mm;
  halfPenSampleWidth = 4.5*mm;
  halfCollimatorThickness = 5.0*mm;
  nSamples = 4;
  fDetectorType = 1;
  fDetectorCollimatorX = 0*mm;
  AbsorptionLength = 20.;//value at 400 nm
  fRES = 1.0;
  fLY = 5500./MeV;
  fDetectorName = "PenAttenuationSetup";
  fABSFile = "Exp4_long";
  fVolName = "World";
  fSigAlpha = 0.50;
  DefineMaterials();
  SetTargetMaterial("PEN");
  SetWorldMaterial("Air");
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::~DetectorConstruction(){
  delete physicWorldBox;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

//Set the position of the collimator in X
void DetectorConstruction::SetDetectorCollimatorX(G4double value){
  fDetectorCollimatorX  = value;
  G4RunManager::GetRunManager()->ReinitializeGeometry();
}
void DetectorConstruction::SetLY(G4double value){
  fLY=value;
  MPT_PEN->AddConstProperty("SCINTILLATIONYIELD",fLY/MeV);
}

void DetectorConstruction::SetReflectorOn(G4bool b) {   
  reflectorOn = b;  
  G4RunManager::GetRunManager()->ReinitializeGeometry(); 
}

void DetectorConstruction::SetRes(G4double value){
  fRES=value;
  DefineMaterials();
  UpdateGeometry();
  G4MTRunManager::GetRunManager()->PhysicsHasBeenModified();
}

/*
Sets which detector geometry is used.
*/
void DetectorConstruction::SetDetectorType(G4int value){
  fDetectorType=value;
  
  if(value == 0){
	nSamples =1;
	halfPenSampleLength = 15. * mm;
	halfPenSampleWidth = 15. * mm;
  }
  if(value == 2){
	halfPenSampleLength = 15. * mm;
	halfPenSampleWidth = 15. * mm;
  }
  G4RunManager::GetRunManager()->ReinitializeGeometry();
  //UpdateGeometry();
  //G4MTRunManager::GetRunManager()->PhysicsHasBeenModified();
}

void DetectorConstruction::SetNumberOfTargetSamples(G4int value){
  nSamples = value;
  G4RunManager::GetRunManager()->ReinitializeGeometry();
}

//Sets dimmensions of target, thickness corresponds to the Y coordinate, Length to x.
void DetectorConstruction::SetTargetSampleLength(G4double value){
  halfPenSampleLength = (value/2.)*mm;
  G4RunManager::GetRunManager()->ReinitializeGeometry();
}

void DetectorConstruction::SetTargetSampleThickness(G4double value){
  halfPenSampleThickness = (value/2.)*mm;
  G4RunManager::GetRunManager()->ReinitializeGeometry();
}

void DetectorConstruction::SetTargetSampleWidth(G4double value){
  halfPenSampleWidth = (value/2.)*mm;
  G4RunManager::GetRunManager()->ReinitializeGeometry();
}
void DetectorConstruction::SetABS(G4double value){
  AbsorptionLength=value;
  //read file and add the value given by the user
  G4double wavelength;
  char filler;
  G4double varAbsorLength;
  G4double emission;
  G4double rindex;

  G4double wlPhotonEnergy[102]  = {0};
  G4double ABSORPTION_PEN[102] = {0};

  G4int absEntries = 0;
  ifstream ReadAbs;

  G4String absFile = "../input_files/"+fABSFile+".csv";
  ReadAbs.open(absFile);
  if(ReadAbs.is_open())
  {
    while(!ReadAbs.eof())
    {
      ReadAbs>>wavelength>>filler>>varAbsorLength>>filler>>emission>>filler>>rindex;
      if(ReadAbs.eof()){
        break;
      }
      wlPhotonEnergy[absEntries] = (1240./wavelength)*eV;
      ABSORPTION_PEN[absEntries] = (varAbsorLength + AbsorptionLength)*mm;
      absEntries++;
    }
  }

  else G4cout<<"Error opening file: " <<absFile<<G4endl;
  ReadAbs.close();
  absEntries--;

  const G4int nEntries1 = sizeof(wlPhotonEnergy)/sizeof(G4double);
  assert(sizeof(ABSORPTION_PEN) == sizeof(wlPhotonEnergy));
  MPT_PEN->AddProperty("ABSLENGTH",    wlPhotonEnergy, ABSORPTION_PEN, nEntries1)->SetSpline(true); // *
}

void DetectorConstruction::SetSigAlpha(G4double value){
  fSigAlpha=value;
  AirPEN -> SetPolish(fSigAlpha);
  AirPEN -> SetMaterialPropertiesTable(MPT_PEN);
}

void DetectorConstruction::SetDetectorName(G4String name){
  fDetectorName=name;
}

void DetectorConstruction::SetABSFile(G4String fileName){
  fABSFile = fileName;
  DefineMaterials();
  UpdateGeometry();
  G4MTRunManager::GetRunManager()->PhysicsHasBeenModified();
  DefineMaterials();
}


/*
Sets material of target.
*/
void DetectorConstruction::SetTargetMaterial(G4String materialChoice)
{
  // search the material by its name
  G4Material* pttoMaterial = G4NistManager::Instance()->FindOrBuildMaterial(materialChoice);

  if (pttoMaterial) {
    fTargetMaterial = pttoMaterial;
    fTargetName = fTargetMaterial->GetName();
    if ( penLogicBox ) { penLogicBox->SetMaterial(fTargetMaterial); }
  } else {
    G4cout << "\n--> warning from DetectorConstruction::SetMaterial : "
           << materialChoice << " not found" << G4endl;
  }
  G4MTRunManager::GetRunManager()->PhysicsHasBeenModified();
}

void DetectorConstruction::SetRI(G4double value){
  fRI = value;
  UpdateGeometry();
  G4MTRunManager::GetRunManager()->PhysicsHasBeenModified();
  DefineMaterials();
}

/*
Sets material of world volume.
*/
void DetectorConstruction::SetWorldMaterial(G4String materialChoice)
{
  // search the material by its name
  G4Material* pttoMaterial =
     G4NistManager::Instance()->FindOrBuildMaterial(materialChoice);
  if (pttoMaterial) {
    fWorldMaterial = pttoMaterial;
    if ( logicWorldBox ) { logicWorldBox->SetMaterial(fWorldMaterial); }
  } else {
    G4cout << "\n--> warning from DetectorConstruction::SetMaterial : "
           << materialChoice << " not found" << G4endl;
  }
  G4MTRunManager::GetRunManager()->PhysicsHasBeenModified();
}

/*
Defines materials used in simulation. Sets material properties for PEN and other optical components.
*/
void DetectorConstruction::DefineMaterials(){
  // ============================================================= Materials =============================================================
  G4double a, z, density;
  G4int nelements;

  // materialAir
  G4Element* N = new G4Element("Nitrogen", "N", z=7 , a=14.01*g/mole);
  G4Element* O = new G4Element("Oxygen"  , "O", z=8 , a=16.00*g/mole);

  materialAir = new G4Material("Air", density=1.29*mg/cm3, nelements=2);
  materialAir->AddElement(N, 70.*perCent);
  materialAir->AddElement(O, 30.*perCent);

  G4NistManager* man = G4NistManager::Instance();
  // Water
  G4Element* H = new G4Element("Hydrogen", "H", z=1 , a=1.01*g/mole);

  G4Material* water = new G4Material("Water", density= 1.0*g/cm3, nelements=2);
  water->AddElement(H, 2);
  water->AddElement(O, 1);

  G4Element* C = new G4Element("Carbon", "C", z=12, a=12*g/mole);
  //G4Element* Pb = new G4Element("Lead", "Pb", z=87, a=207*g/mole);
  fGlass = man->FindOrBuildMaterial("G4_Pyrex_Glass");
  fPOM = new G4Material("POM",density=1.41*g/cm3,nelements=3);
  fPOM->AddElement(O,1);
  fPOM->AddElement(C,1);
  fPOM->AddElement(H,2);

  fABS = new G4Material("ABS",density=1.07*g/cm3,nelements=3);
  fABS->AddElement(C,15);
  fABS->AddElement(H,17);
  fABS->AddElement(N,1);

  // Scintillators
  G4int number_of_atoms;
  PenMaterial = new G4Material("PEN", density= 1.3*g/cm3, nelements=3);
  PenMaterial->AddElement(O, number_of_atoms=4);
  PenMaterial->AddElement(H, number_of_atoms=10);
  PenMaterial->AddElement(C, number_of_atoms=14);

  G4double wavelength;
  char filler;
  G4double varAbsorLength;
  G4double emission;
  G4double rindex;

  G4double wlPhotonEnergy[102]  = {0};
  G4double ABSORPTION_PEN[102] = {0};
  G4double EMISSION_PEN[102] = {0};
  G4double RINDEX_PEN[102] = {0};
  G4double RINDEX_Air[102] = {0};
  //G4double emsAbs[102] = {0};

  G4int absEntries = 0;
  ifstream ReadAbs;

  G4String absFile = "../input_files/"+fABSFile+".csv";
  ReadAbs.open(absFile);
  if(ReadAbs.is_open())
  {
    while(!ReadAbs.eof())
    {
      ReadAbs>>wavelength>>filler>>varAbsorLength>>filler>>emission>>filler>>rindex;
      if(ReadAbs.eof()){
        break;
      }
      wlPhotonEnergy[absEntries] = (1240/wavelength)*eV;
      //ABSORPTION_PEN[absEntries] = (varAbsorLength+AbsorptionLength)*mm;
      ABSORPTION_PEN[absEntries] = (varAbsorLength)*mm;
      EMISSION_PEN[absEntries] = emission;
      RINDEX_PEN[absEntries] = rindex;
      RINDEX_Air[absEntries] = 1.0;
      //RINDEX_PEN[absEntries] = 1.65;
      //emsAbs[absEntries] = 0.02;
      //emission_fibre[absEntries] = 1.0;
      absEntries++;
    }
  }

  else G4cout<<"Error opening file: " <<absFile<<G4endl;
  ReadAbs.close();
  absEntries--;

  const G4int nEntries1 = sizeof(wlPhotonEnergy)/sizeof(G4double);
  assert(sizeof(RINDEX_PEN) == sizeof(wlPhotonEnergy));
  assert(sizeof(ABSORPTION_PEN) == sizeof(wlPhotonEnergy));
  assert(sizeof(EMISSION_PEN) == sizeof(wlPhotonEnergy));
  assert(sizeof(RINDEX_Air == sizeof(wlPhotonEnergy)));
  
  MPT_PEN = new G4MaterialPropertiesTable();

  MPT_PEN->AddProperty("RINDEX",       wlPhotonEnergy, RINDEX_PEN, nEntries1)->SetSpline(true);
  MPT_PEN->AddProperty("ABSLENGTH",    wlPhotonEnergy, ABSORPTION_PEN, nEntries1)->SetSpline(true); // *
  MPT_PEN->AddProperty("FASTCOMPONENT",wlPhotonEnergy, EMISSION_PEN, nEntries1)->SetSpline(true);
  MPT_PEN->AddProperty("SLOWCOMPONENT",wlPhotonEnergy, EMISSION_PEN, nEntries1)->SetSpline(true);

  MPT_PEN->AddConstProperty("SCINTILLATIONYIELD",fLY/MeV); // * 2.5 * PEN = PS, 10*PEN=PS
  MPT_PEN->AddConstProperty("RESOLUTIONSCALE",fRES); // * 1, 4, 8
  MPT_PEN->AddConstProperty("FASTTIMECONSTANT", 5.198*ns);
  MPT_PEN->AddConstProperty("SLOWTIMECONSTANT",24.336*ns);
  MPT_PEN->AddConstProperty("YIELDRATIO",1.0);

  PenMaterial->SetMaterialPropertiesTable(MPT_PEN);

  AirPEN = new G4OpticalSurface("AirPEN",glisur, ground, dielectric_dielectric);
  AirPEN -> SetPolish(fSigAlpha);
  AirPEN -> SetMaterialPropertiesTable(MPT_PEN);

  density = universe_mean_density;    //from PhysicalConstants.h
  fVacuum = new G4Material("Galactic", z=1., a=1.008*g/mole, density,
                           kStateGas,2.73*kelvin,3.e-18*pascal);
  //
  // fAir
  G4MaterialPropertiesTable* MPT_Air = new G4MaterialPropertiesTable();
  MPT_Air->AddProperty("RINDEX", wlPhotonEnergy, RINDEX_Air, nEntries1)->SetSpline(true);

  materialAir->SetMaterialPropertiesTable(MPT_Air);
  fVacuum->SetMaterialPropertiesTable(MPT_Air);
  materialSi = man->FindOrBuildMaterial("G4_Si");

  materialTriggerFoilEJ212 = new G4Material("EJ212", density= 1.03*g/cm3, 2);
  materialTriggerFoilEJ212->AddElement(C, 0.475);
  materialTriggerFoilEJ212->AddElement(H, 0.525);

  Pstyrene = new G4Material("Polystyrene", density= 1.03*g/cm3, 2);
  Pstyrene->AddElement(C, 8);
  Pstyrene->AddElement(H, 8);

  G4double rindexEnergy[500] = {0};
  G4double scintIndex[500] = {0};

  G4int rindexEntries = 0;
  ifstream ReadRindex;

  G4String rindex_file="../input_files/rindexScint.txt";
  ReadRindex.open(rindex_file);

  if(ReadRindex.is_open())
  {
      while(!ReadRindex.eof())
      {
          ReadRindex>>wavelength>>filler>>scintIndex[rindexEntries];
          rindexEnergy[rindexEntries] = (1240./wavelength)*eV;
          rindexEntries++;
      }
  }else G4cout<<"Error opening file: "<<rindex_file<<G4endl;
  ReadRindex.close();
  rindexEntries--;

  G4double scintEnergy[501] = {0};
  G4double scintEmit[501] = {0};
  G4double scintEmitSlow[501] = {0};

  G4int scintEntries = 0;
  ifstream ReadScint;

  G4String Scint_file="../input_files/pTP_emission.txt";
  ReadScint.open(Scint_file);

  if(ReadScint.is_open())
  {
      while(!ReadScint.eof())
      {
          ReadScint>>wavelength>>filler>>scintEmit[scintEntries];
          //convert wavelength to eV:
          scintEnergy[scintEntries] = (1240./wavelength)*eV;
          scintEmitSlow[scintEntries] = scintEmit[scintEntries];
          scintEntries++;
      }
  }else G4cout<<"Error opening file: "<<Scint_file<<G4endl;
  ReadScint.close();
  scintEntries--;

  G4int absorbEntries = 0;
  G4double varAbsorbLength;
  G4double absorbEnergy[501] = {0};
  G4double Absorb[501] = {0};

  ifstream ReadAbsorb;
  G4String ReadAbsorbLength="../input_files/PlasticBulkAbsorb2.cfg";

  ReadAbsorb.open(ReadAbsorbLength);
  if (ReadAbsorb.is_open())
  {
      while(!ReadAbsorb.eof())
      {
          ReadAbsorb>>wavelength>>filler>>varAbsorbLength;
          absorbEnergy[absorbEntries]=(1240/wavelength)*eV;
          Absorb[absorbEntries]=(varAbsorbLength+10.)*m;
          absorbEntries++;
      }
  }else G4cout<<"Error opening file: "<<ReadAbsorb<<G4endl;
  ReadAbsorb.close();
  absorbEntries--;

  G4double wlsEnergy[501] = {0};
  G4double wlsEmit[501] = {0};

  G4int wlsScintEntries = 0;
  ifstream ReadWLSScint;

  G4String wls_Scint_file="../input_files/full_popop_emission.cfg";
  ReadWLSScint.open(wls_Scint_file);

  if(ReadWLSScint.is_open())
  {
      while(!ReadWLSScint.eof())
      {
          ReadWLSScint>>wavelength>>filler>>wlsEmit[500-wlsScintEntries];
          //convert wavelength to eV:
          wlsEnergy[500 - wlsScintEntries] = (1240/wavelength)*eV;
          wlsScintEntries++;
      }
  }else G4cout<<"Error opening file: "<<wls_Scint_file<<G4endl;
  ReadWLSScint.close();
  wlsScintEntries--;

  G4int wlsAbsorbEntries = 0;
  G4double wlsAbsorbEnergy[501] = {0};
  G4double wlsAbsorb[501] = {0};

  ifstream ReadWLSAbsorb;
  G4String ReadWLSAbsorbLength="../input_files/scintAbsLen.txt";

  ReadWLSAbsorb.open(ReadWLSAbsorbLength);
  if (ReadWLSAbsorb.is_open())
  {
      while(!ReadWLSAbsorb.eof())
      {
          ReadWLSAbsorb>>wavelength>>filler>>varAbsorbLength;
          wlsAbsorbEnergy[wlsAbsorbEntries]=(1240./wavelength)*eV;
          wlsAbsorb[wlsAbsorbEntries]=varAbsorbLength*m;
          wlsAbsorbEntries++;
      }
  }else G4cout<<"Error opening file: "<<ReadWLSAbsorb<<G4endl;
  ReadWLSAbsorb.close();
  wlsAbsorbEntries--;

  G4MaterialPropertiesTable* MPT_FoilEJ212 = new G4MaterialPropertiesTable();

  MPT_FoilEJ212->AddProperty("WLSABSLENGTH",wlsAbsorbEnergy,wlsAbsorb,wlsAbsorbEntries);
  MPT_FoilEJ212->AddProperty("WLSCOMPONENT",wlsEnergy,wlsEmit,wlsScintEntries);
  MPT_FoilEJ212->AddConstProperty("WLSTIMECONSTANT", 12*ns);

  MPT_FoilEJ212->AddProperty("RINDEX",        rindexEnergy,  scintIndex, rindexEntries);
  MPT_FoilEJ212->AddProperty("ABSLENGTH",     absorbEnergy, Absorb,     absorbEntries);
  MPT_FoilEJ212->AddProperty("FASTCOMPONENT", scintEnergy,  scintEmit,  scintEntries);
  MPT_FoilEJ212->AddProperty("SLOWCOMPONENT",scintEnergy, scintEmitSlow,     scintEntries);

  MPT_FoilEJ212->AddConstProperty("SCINTILLATIONYIELD",11520./MeV);
  //MPT_FoilEJ212->AddConstProperty("SCINTILLATIONYIELD",1000./MeV);
  MPT_FoilEJ212->AddConstProperty("RESOLUTIONSCALE",4.0);
  MPT_FoilEJ212->AddConstProperty("FASTTIMECONSTANT", 2.1*ns);
  MPT_FoilEJ212->AddConstProperty("SLOWTIMECONSTANT",14.2*ns);
  MPT_FoilEJ212->AddConstProperty("YIELDRATIO",1.0);

  materialTriggerFoilEJ212->SetMaterialPropertiesTable(MPT_FoilEJ212);

  // G4double density, a;
  // G4int nelements, z;
  // G4Element* O = new G4Element("Oxygen"  , "O", z=8 , a=16.00*g/mole);
  // G4Element* C = new G4Element("Carbon", "C", z=12, a=12*g/mole);
  // G4Element* H = new G4Element("Hydrogen", "H", z=1, a=1.01*g/mole);

  materialPMMA = new G4Material("PMMA", density = 1.18*g/cm3, nelements = 3);
  materialPMMA->AddElement(C, 5);
  materialPMMA->AddElement(O, 2);
  materialPMMA->AddElement(H, 8);

  G4double refractive_index[] = {1.49, 1.49, 1.49, 1.49, 1.49, 1.49};
  G4double absPMMA[] = {5*m, 5*m, 5*m, 5*m, 5*m, 5*m};
  G4double reflPMMA[] = {0.1, 0.1, 0.1, 0.1, 0.1, 0.1};
  G4double energyPMMA[] = {2.18*eV, 2.48*eV, 2.58*eV, 2.68*eV, 2.78*eV, 4.1*eV};
  const G4int nEntries3 = sizeof(energyPMMA)/sizeof(G4double);

  G4MaterialPropertiesTable* MPT_PMMA = new G4MaterialPropertiesTable();
  MPT_PMMA->AddProperty("RINDEX", energyPMMA, refractive_index, nEntries3);
  MPT_PMMA->AddProperty("ABSLENGTH", energyPMMA, absPMMA, nEntries3)->SetSpline(true);
  MPT_PMMA->AddProperty("REFLECTIVITY", energyPMMA, reflPMMA, nEntries3)->SetSpline(true);
  materialPMMA->SetMaterialPropertiesTable(MPT_PMMA);

  ej_550 = man->FindOrBuildMaterial("G4_WATER");
  G4MaterialPropertiesTable* MPT_EJ550Grease = new G4MaterialPropertiesTable();
  G4double ej_refractive_index[] = {1.46, 1.46, 1.46, 1.46, 1.46, 1.46};
  MPT_EJ550Grease->AddProperty("ABSLENGTH", energyPMMA, absPMMA, nEntries3);
  MPT_EJ550Grease->AddProperty("RINDEX", energyPMMA, ej_refractive_index, nEntries3)->SetSpline(true);
  ej_550->SetMaterialPropertiesTable(MPT_EJ550Grease);
}

void DetectorConstruction::SetVolName(G4ThreeVector thePoint){
  G4Navigator* theNavigator = G4TransportationManager::GetTransportationManager()->GetNavigatorForTracking();
  G4VPhysicalVolume* myVolume= theNavigator->LocateGlobalPointAndSetup(thePoint);
  fVolName =  myVolume->GetName();
}

void DetectorConstruction::SetPropertyTable(G4Material* material, G4MaterialPropertiesTable* table){
  material->SetMaterialPropertiesTable(table);
}

void DetectorConstruction::UpdateGeometry(){
  // clean-up previous geometry
  G4SolidStore::GetInstance()->Clean();
  G4LogicalVolumeStore::GetInstance()->Clean();
  G4PhysicalVolumeStore::GetInstance()->Clean(); 
  
  //define new one
  G4RunManager::GetRunManager()->DefineWorldVolume(Construct());
  G4RunManager::GetRunManager()->GeometryHasBeenModified();

  //G4MTRunManager::GetRunManager()->DefineWorldVolume(Construct());
}

/*
Clears stored geometry, then constructs all volumes that can be used in the simulation.

Builds and places volumes in world.

Defines detector sensitivities and properties.
*/
G4VPhysicalVolume* DetectorConstruction::Construct()
{
  G4NistManager* man = G4NistManager::Instance();

  G4Material* teflon = man->FindOrBuildMaterial("G4_TEFLON");
  G4Material* polyethylene = man->FindOrBuildMaterial("G4_POLYETHYLENE");
  //G4Material* air = man->FindOrBuildMaterial("G4_AIR");

// ============================================================= Define Volumes =============================================================

  //The experimental Dark Box
  fWorldBox = new G4Box("World",halfSizeDarkBoxX,halfSizeDarkBoxY,halfSizeDarkBoxZ);
  logicWorldBox = new G4LogicalVolume(fWorldBox,materialAir,"World",0,0,0);
  physicWorldBox = new G4PVPlacement(0,G4ThreeVector(),logicWorldBox,"World",0,false,0);

  penSampleBox = new G4Box("target", halfPenSampleLength, halfPenSampleThickness, halfPenSampleWidth);
  penLogicBox = new G4LogicalVolume(penSampleBox,PenMaterial, "target",0,0,0);


  //EJ212 foil scintillator used for trigger
  double halfThicknessTriggerFoilEJ212 = 50.*um;//100 um thickness EJ foil for trigger
  double halfWidthTriggerFoilEJ212 = 7.5*mm;
  double halfLengthTriggerFoilEJ212 = 15.*mm;
  G4Box* boxTriggerFoilEJ212 = new G4Box("triggerFoilEJ212", halfLengthTriggerFoilEJ212, halfThicknessTriggerFoilEJ212, halfWidthTriggerFoilEJ212);
  G4LogicalVolume* logicBoxTriggerFoilEJ212 = new G4LogicalVolume(boxTriggerFoilEJ212, materialTriggerFoilEJ212, "triggerFoilEJ212", 0, 0, 0);

  //Reflector foil around EJ212 scintillator
  //Define a box with larger dimmensions than EJ212 box and then subract EJ212 scintillator
  double halfThicknessReflectorFoil = 100.*um;
  double halfWidththReflectorFoil = 7.8*mm;
  double halfLengthReflectorFoil = 14.0*mm;
  G4VSolid* boxReflectorFoil = new G4Box("foil", halfLengthReflectorFoil, halfThicknessReflectorFoil, halfWidththReflectorFoil);
  G4SubtractionSolid* reflectorFoilAroundEJ212Foil = new G4SubtractionSolid("reflectorFoilAroundEJ212Foil", boxReflectorFoil, boxTriggerFoilEJ212, 0, G4ThreeVector(0, 0, 0));
  G4LogicalVolume* logicReflectorFoilAroundEJ212Foil = new G4LogicalVolume(reflectorFoilAroundEJ212Foil, teflon, "foil", 0, 0, 0);

  //Reflector foil to be placed over pen samples, not used in attenuation
  double totalReflectorFoilThickness = 100.*um;
  double halfReflectorBoxOverPEN = 3.0*mm;

  G4VSolid* boxReflectorPEN = new G4Box("foilPEN",halfPenSampleLength, halfReflectorBoxOverPEN, halfPenSampleWidth);
  G4Box* reflectorFoilBox = new G4Box("foil", halfPenSampleLength-totalReflectorFoilThickness, halfReflectorBoxOverPEN, halfPenSampleWidth-totalReflectorFoilThickness);
  G4SubtractionSolid* reflectorFoilOverPEN = new G4SubtractionSolid("reflectorFoilOverPEN",boxReflectorPEN,reflectorFoilBox,0,G4ThreeVector(0, totalReflectorFoilThickness, 0));
  //G4LogicalVolume* logicReflectorFoilBox = new G4LogicalVolume(reflectorFoilBox, teflon, "foil", 0, 0, 0);
  G4LogicalVolume* logicReflectorFoilBox = new G4LogicalVolume(reflectorFoilOverPEN, teflon, "foil", 0, 0, 0);

  //Passive collimator
  double innerRadiusCollimator = 1.*mm; 
  double externalRadiusCollimator = 12.*mm; 
  halfCollimatorThickness= 5.*mm; 
  
  G4Tubs* collimatorTube = new G4Tubs("collimator",innerRadiusCollimator,externalRadiusCollimator,halfCollimatorThickness,0.,360.*deg);
  G4LogicalVolume* logicCollimator = new G4LogicalVolume(collimatorTube,polyethylene,"collimator",0,0,0); 


  G4ThreeVector point = G4ThreeVector(0,0,5*cm);
  G4Navigator* pointNavigator = new G4Navigator();
  pointNavigator->SetWorldVolume(physicWorldBox);
  pointNavigator->LocateGlobalPointAndSetup(point);

  // ============================================================= Detectors =============================================================

  char filler;
  G4double wavelength;
  G4double cathodeEfficiency;
  G4double photocathEnergy[37];
  G4double photoCathodeQuantumEfficiency[37];
  G4double perfectEfficiency[37];
  G4double photoCathodeRelectivity[37];
  //G4String pmtFile = "../input_files/pmtQE.csv";
  G4String pmtFile = "../input_files/QE_H11934_300.csv";

  ifstream ReadEff;
  G4int effCounter = 0;
  ReadEff.open(pmtFile);

  if(ReadEff.is_open())
  {
    while(!ReadEff.eof())
    {
      ReadEff>>wavelength>>filler>>cathodeEfficiency;
      if(ReadEff.eof()){
        break;
      }
      photocathEnergy[effCounter] = (1240./wavelength)*eV;
      photoCathodeQuantumEfficiency[effCounter] = cathodeEfficiency/100.;
      perfectEfficiency[effCounter] = 1.;
      photoCathodeRelectivity[effCounter] = 0.;
      effCounter++;
    }
  }

  else G4cout<<"Error opening file: " <<pmtFile<<G4endl;
  ReadEff.close();
  effCounter--;

  const G4int nPMT_EFF = sizeof(photocathEnergy)/sizeof(G4double);

  G4OpticalSurface* perfectOptSurf = new G4OpticalSurface("perfect",glisur,polished, dielectric_metal);
  G4MaterialPropertiesTable* MPT_Detector = new G4MaterialPropertiesTable();
  MPT_Detector->AddProperty("EFFICIENCY", photocathEnergy, perfectEfficiency,nPMT_EFF);
  MPT_Detector->AddProperty("REFLECTIVITY", photocathEnergy, photoCathodeRelectivity,nPMT_EFF);
  perfectOptSurf->SetMaterialPropertiesTable(MPT_Detector);

  G4OpticalSurface* pmtOpticalSurface = new G4OpticalSurface("pmt", glisur, polished, dielectric_metal);
  G4MaterialPropertiesTable* MPT_PMT = new G4MaterialPropertiesTable();
  MPT_PMT->AddProperty("EFFICIENCY", photocathEnergy, photoCathodeQuantumEfficiency,nPMT_EFF);
  MPT_PMT->AddProperty("REFLECTIVITY", photocathEnergy, photoCathodeRelectivity,nPMT_EFF);
  pmtOpticalSurface->SetMaterialPropertiesTable(MPT_PMT);

  G4LogicalVolume* tileDetectorLog = new G4LogicalVolume(penSampleBox,materialSi, "tile_sensor");
  new G4LogicalSkinSurface("pen_det_surf",tileDetectorLog,perfectOptSurf);

  LightGuideConstruction guide;
  G4VSolid* LightGuidePMMA = guide.ConstructPlate();
  //G4cout <<  G4BestUnit(logicLightGuidePMMA->GetMass(true),"Mass") << G4endl;


  //Define PMTs
  G4double casingPMTLength = 30.*mm / 2.;
  G4double casingPMTWidth = 30.*mm / 2.;
  G4double casingPMTHeight = 32.5*mm / 2.;
  G4double casingPMTHoleLength = 29.*mm / 2.;
  G4double casingPMTHoleWidth = 29.*mm / 2.;
  G4double casingPMTHoleHeight = 31.5*mm / 2.;

  G4double inactivePhotoCathodePMTLength = 26.2*mm / 2.;
  G4double inactivePhotoCathodePMTWidth = 26.2*mm / 2.;
  G4double activePhotoCathodePMTLength = 23.*mm / 2.;
  G4double activePhotoCathodePMTWidth = 23.*mm / 2.;
  G4double activePhotoCathodePMTThickness = 0.8*mm;

  G4Box* boxCasingPMT = new G4Box("case",casingPMTLength,casingPMTWidth,casingPMTHeight);
  G4Box* boxEmptyInsidePMT = new G4Box("hole",casingPMTHoleLength,casingPMTHoleWidth,casingPMTHoleHeight);
  G4SubtractionSolid* boxPMTShell = new G4SubtractionSolid("boxPMTShell",boxCasingPMT,boxEmptyInsidePMT);

  G4LogicalVolume* logicBoxEmptyInsidePMT = new G4LogicalVolume(boxEmptyInsidePMT, fVacuum,"vacuum");
  G4LogicalVolume* logicBoxPMTShell = new G4LogicalVolume(boxPMTShell,fPOM,"logicBoxPMTShell");

  G4Box* boxTotalPhotoCathodeSupport = new G4Box("pmt_cathode",
		inactivePhotoCathodePMTLength,
		activePhotoCathodePMTThickness,
		inactivePhotoCathodePMTWidth);
  G4Box* boxActivePhotoCathodePMT = new G4Box("pmt_active",
		activePhotoCathodePMTLength,
		activePhotoCathodePMTThickness,
		activePhotoCathodePMTWidth);
  G4SubtractionSolid* inactiveBoxPhotoCathodeSupport = new G4SubtractionSolid("inactiveBoxPhotoCathodeSupportode",
		boxTotalPhotoCathodeSupport,
		boxActivePhotoCathodePMT);
  G4LogicalVolume* logicBoxPhotoCathodeSupport = new G4LogicalVolume(inactiveBoxPhotoCathodeSupport,
		fGlass,
		"logicBoxPhotoCathodeSupport");
  G4LogicalVolume* logicBoxActivePhotoCathodePMT = new G4LogicalVolume(boxActivePhotoCathodePMT,
		fGlass,
		"logicBoxActivePhotoCathodePMT");
  new G4LogicalSkinSurface("pmt_surf", logicBoxActivePhotoCathodePMT,pmtOpticalSurface);
  //new G4LogicalSkinSurface("spec_surf",logicBoxActivePhotoCathodePMT,perfectOptSurf);

  //Define orientation of PMTs
  G4RotationMatrix* rotationMatrix = new G4RotationMatrix(0,0,0);
  rotationMatrix->rotateZ(90*deg);
  G4RotationMatrix* rotationMatrix1 = new G4RotationMatrix(0,0,0);
  rotationMatrix1->rotateZ(90*deg);
  rotationMatrix1->rotateX(180*deg);
  G4RotationMatrix* rotationMatrix2 = new G4RotationMatrix(0,0,0);
  rotationMatrix2->rotateX(90*deg);
  G4RotationMatrix* rotationMatrix3 = new G4RotationMatrix(0,0,0);
  rotationMatrix3->rotateX(90*deg);
  rotationMatrix3->rotateZ(180*deg);
  G4RotationMatrix* rotationMatrix4 = new G4RotationMatrix(0,0,0);
  rotationMatrix4->rotateX(180*deg);
  G4RotationMatrix* rotationMatrix5 = new G4RotationMatrix(0,0,0);
  rotationMatrix5->rotateZ(90*deg);
  rotationMatrix5->rotateX(90*deg);
  G4RotationMatrix* rotationMatrixCollimator = new G4RotationMatrix(0,0,0);
  rotationMatrixCollimator->rotateX(90*deg);

  
  // Set Draw G4VisAttributes
  G4VisAttributes* visAttr = new G4VisAttributes();
  visAttr->SetVisibility(false);
  logicWorldBox->SetVisAttributes(visAttr);
  logicBoxEmptyInsidePMT->SetVisAttributes(visAttr);

  //PEN and EJ212 
  G4VisAttributes* visualAttributesScintillators = new G4VisAttributes(G4Colour::Blue());
  visualAttributesScintillators->SetVisibility(true);
  penLogicBox->SetVisAttributes(visualAttributesScintillators);
  logicBoxTriggerFoilEJ212->SetVisAttributes(visualAttributesScintillators);

  //Collimator
  G4VisAttributes* visAttributesCollimator = new G4VisAttributes(G4Colour::White());
  visAttributesCollimator->SetVisibility(true);
  visAttributesCollimator->SetForceSolid(true);
  logicCollimator->SetVisAttributes(visAttributesCollimator);

  //reflector foils
  G4VisAttributes* visulaAttributesReflectors = new G4VisAttributes(G4Colour::Red());
  visulaAttributesReflectors->SetVisibility(true);
  logicReflectorFoilAroundEJ212Foil->SetVisAttributes(visulaAttributesReflectors);
  logicReflectorFoilBox->SetVisAttributes(visulaAttributesReflectors);

  // Active Detectors
  G4VisAttributes* visulaAttributesDetectors = new G4VisAttributes(G4Colour::Yellow());
  visulaAttributesDetectors->SetVisibility(true);
  visulaAttributesDetectors->SetForceSolid(true);
  logicBoxActivePhotoCathodePMT->SetVisAttributes(visulaAttributesDetectors);

  // Inactive volumes
  G4VisAttributes* visualAttributesInactiveMaterials = new G4VisAttributes(G4Colour::Gray());
  logicBoxPhotoCathodeSupport->SetVisAttributes(visualAttributesInactiveMaterials);
  logicBoxPMTShell->SetVisAttributes(visualAttributesInactiveMaterials);

  man = G4NistManager::Instance();

  G4LogicalVolume* logicLightGuidePMMA = new G4LogicalVolume(LightGuidePMMA, materialPMMA, "Ligh_guideLog");

  G4double greaseHoleDiameter = 23.5*mm;
  G4double greaseHoleDepth = 1*mm;
  G4Tubs* tubeSiGrease = new G4Tubs("SiGrease", 0, greaseHoleDiameter/2., greaseHoleDepth/2., 0, 360*deg);

  G4LogicalVolume* logicTubeSiGrease = new G4LogicalVolume(tubeSiGrease, ej_550, "logicTubeSiGrease");
  logicTubeSiGrease->SetVisAttributes(visulaAttributesReflectors);


  //  ============================================================= Place volumes =============================================================
  
  // Place main tile at centre of world volume
  switch(fDetectorType){

  case 0:
    for(int iSample = 0; iSample < nSamples; iSample++){
  	physicPenStackedSamples = new G4PVPlacement(0, 
				G4ThreeVector(0,iSample*2*halfPenSampleThickness + iSample*5.*um,0),
				penLogicBox,
				"target_"+std::to_string(iSample+1),
				logicWorldBox,false,iSample,false);
     }

     //physicPenSampleBox = new G4PVPlacement(0, 
     //			G4ThreeVector(0,0,0),
     //			penLogicBox,
     //			"target",
     //			logicWorldBox,false,0,false);
     // Trigger and light guide placement, with triggerFoilEJ212 PMT
     if(reflectorOn){
     		physicReflectorFoilAroundEJ212Foil = new G4PVPlacement(0, 
			G4ThreeVector(0,18.*mm,0), 
			logicReflectorFoilAroundEJ212Foil, 
			"foil", 
			logicWorldBox, false, 0, false);
     }
     physicTriggerFoilEJ212 = new G4PVPlacement(0, 
			//G4ThreeVector(0,0,0), 
			G4ThreeVector(0,18.*mm,0), 
			logicBoxTriggerFoilEJ212, 
			"triggerFoilEJ212", 
			logicWorldBox, false, 0, false);
			//logicReflectorFoilAroundEJ212Foil, false, 0, false);
     physicLightGuide = new G4PVPlacement(rotationMatrix4, 
			G4ThreeVector(-(halfLengthTriggerFoilEJ212 + 20*mm - 0.4*mm), (18*mm + 14*mm - 1.75*mm), 0),
			logicLightGuidePMMA,
			"LightGuide",
			logicWorldBox,false,0,false);
     physicOpticalGrease =  new G4PVPlacement(rotationMatrix5, 
			G4ThreeVector(-19*mm - (halfLengthTriggerFoilEJ212 + 20*mm), (18*mm + 14*mm - 1.75*mm),0), 
			logicTubeSiGrease, 
			"grease", 
			logicWorldBox, false, 0, false);
     physicActivePhotoCathodeTriggerPMT = new G4PVPlacement(rotationMatrix,
			G4ThreeVector(-(20*mm + activePhotoCathodePMTThickness),0,0),
			logicBoxActivePhotoCathodePMT,
			"trigger_pmt",
			logicLightGuidePMMA,false,0,false);
     physicBoxPhotoCathodeSupport = new G4PVPlacement(0,
			G4ThreeVector(0,0,0),
			logicBoxPhotoCathodeSupport,
			"inactive_detector1",
			logicBoxActivePhotoCathodePMT,false,0,false);
     physicBoxPMTShell = new G4PVPlacement(0,
			G4ThreeVector(0,-(casingPMTLength + activePhotoCathodePMTThickness),0),
			logicBoxPMTShell,
			"pmt1", 
			logicBoxActivePhotoCathodePMT,false,0,false);
     physicBoxEmptyInsidePMT = new G4PVPlacement(0,
			G4ThreeVector(), 
			logicBoxEmptyInsidePMT, 
			"void1", 
			logicBoxPMTShell, false, 0, false);

     // Main PMT placements
     physicActivePhotoCathodePMT1 = new G4PVPlacement(rotationMatrix, 
			G4ThreeVector(-(halfPenSampleLength + activePhotoCathodePMTThickness),0,0), 
			logicBoxActivePhotoCathodePMT, 
			"main_pmt_1", 
			logicWorldBox, false, 0, false);
     physicActivePhotoCathodePMT2 = new G4PVPlacement(rotationMatrix1, 
			G4ThreeVector((halfPenSampleLength + activePhotoCathodePMTThickness),0,0), 
			logicBoxActivePhotoCathodePMT, 
			"main_pmt_2", 
			logicWorldBox, false, 0, false);
     physicActivePhotoCathodePMT3 = new G4PVPlacement(0, 
     			G4ThreeVector(0,-(halfPenSampleThickness+activePhotoCathodePMTThickness),0), 
     			logicBoxActivePhotoCathodePMT, 
     			"main_pmt_3", 
     			logicWorldBox, false, 0, false);
     physicActivePhotoCathodePMT4 = new G4PVPlacement(rotationMatrix2, 
     			G4ThreeVector(0,0,(halfPenSampleWidth+activePhotoCathodePMTThickness)),
     			logicBoxActivePhotoCathodePMT,
     			"main_pmt_4", 
     			logicWorldBox, false, 0, false);
     physicActivePhotoCathodePMT5 = new G4PVPlacement(rotationMatrix3, 
     			G4ThreeVector(0,0,-(halfPenSampleWidth+activePhotoCathodePMTThickness)),
     			logicBoxActivePhotoCathodePMT, 
     			"main_pmt_5", 
     			logicWorldBox, false, 0, false);
     break;

  case 1:
     //physicPenSampleBox = new G4PVPlacement(0, G4ThreeVector(0,0,0),penLogicBox,"target",logicWorldBox,false,0,false);
     for(int iSample = 0; iSample < nSamples; iSample++){
  	physicPenStackedSamples = new G4PVPlacement(0, 
				G4ThreeVector(0,iSample*2*halfPenSampleThickness+iSample*5.*um,0),
				penLogicBox,
				"target_"+std::to_string(iSample+1),
				logicWorldBox,false,iSample,false);
     }
     //Reflector foil over PEN samples, uncomment if needed
     //physicReflectorFoilBoxOverPEN = new G4PVPlacement(0, 
     //				G4ThreeVector(0,totalReflectorFoilThickness+halfPenSampleThickness*10,0),
     //				logicReflectorFoilBox,
     //				"target",
     //				logicWorldBox,false,0,false);
     //Collimator, Trigger and light guide placement, with triggerFoilEJ212 PMT
     physicTriggerFoilEJ212 = new G4PVPlacement(0,
			G4ThreeVector(fDetectorCollimatorX,18*mm,0), 
			logicBoxTriggerFoilEJ212, 
			"triggerFoilEJ212",
			logicWorldBox, false, 0, false);
     if(reflectorOn){
     		physicReflectorFoilAroundEJ212Foil = new G4PVPlacement(0, 
			G4ThreeVector(fDetectorCollimatorX,18*mm,0), 
			logicReflectorFoilAroundEJ212Foil, 
			"foil", 
			logicWorldBox, false, 0, false);
     }
     physicCollimator = new G4PVPlacement(rotationMatrixCollimator, 
			G4ThreeVector(fDetectorCollimatorX,19*mm+halfCollimatorThickness,0), 
			logicCollimator, 
			"collimator", 
			logicWorldBox, false, 0, false);
     physicLightGuide = new G4PVPlacement(rotationMatrix4, 
			G4ThreeVector(-(halfLengthTriggerFoilEJ212+20*mm-0.4*mm)+fDetectorCollimatorX, (18*mm+14*mm-1.75*mm), 0),
			logicLightGuidePMMA,
			"LightGuide",
			logicWorldBox,false,0,false);
     physicOpticalGrease =  new G4PVPlacement(rotationMatrix5, 
			G4ThreeVector(-19*mm-(halfLengthTriggerFoilEJ212+20*mm)+fDetectorCollimatorX,(18*mm+14*mm-1.75*mm),0), 
			logicTubeSiGrease, 
			"grease", 
			logicWorldBox, false, 0, false);
     physicActivePhotoCathodeTriggerPMT = new G4PVPlacement(rotationMatrix,
			G4ThreeVector(-(20*mm+activePhotoCathodePMTThickness),0,0),
			logicBoxActivePhotoCathodePMT,
			"trigger_pmt",
			logicLightGuidePMMA,false,0,false);
     physicBoxPhotoCathodeSupport = new G4PVPlacement(0,
			G4ThreeVector(0,0,0),
			logicBoxPhotoCathodeSupport,
			"inactive_detector1",
			logicBoxActivePhotoCathodePMT,false,0,false);
     physicBoxPMTShell = new G4PVPlacement(0,
			G4ThreeVector(0,-(casingPMTLength + activePhotoCathodePMTThickness),0),
			logicBoxPMTShell,
			"pmt1", 
			logicBoxActivePhotoCathodePMT,false,0,false);
     physicBoxEmptyInsidePMT = new G4PVPlacement(0,
			G4ThreeVector(), 
			logicBoxEmptyInsidePMT, 
			"void1", 
			logicBoxPMTShell, false, 0, false);

     // Main PMT placements
     physicActivePhotoCathodePMT1 = new G4PVPlacement(rotationMatrix, 
			G4ThreeVector(-(halfPenSampleLength+activePhotoCathodePMTThickness), 0, 0), 
			logicBoxActivePhotoCathodePMT, 
			"main_pmt_1", 
			logicWorldBox, false, 0, false);
     physicActivePhotoCathodePMT2 = new G4PVPlacement(rotationMatrix1, 
			G4ThreeVector((halfPenSampleLength+activePhotoCathodePMTThickness), 0, 0), 
			logicBoxActivePhotoCathodePMT, 
			"main_pmt_2", 
			logicWorldBox, false, 0, false);
     G4MTRunManager::GetRunManager()->GeometryHasBeenModified();
     //physicActivePhotoCathodePMT3 = new G4PVPlacement(0, 
     //			G4ThreeVector(0,-(halfPenSampleThickness+activePhotoCathodePMTThickness),0), 
     //			logicBoxActivePhotoCathodePMT, 
     //			"main_pmt_3", 
     //			logicWorldBox, false, 0, false);
     //physicActivePhotoCathodePMT4 = new G4PVPlacement(rotationMatrix2, 
     //			G4ThreeVector(0,0,(halfPenSampleWidth+activePhotoCathodePMTThickness)),
     //			logicBoxActivePhotoCathodePMT, 
     //			"main_pmt_4", 
     //			logicWorldBox, false, 0, false);
     //physicActivePhotoCathodePMT5 = new G4PVPlacement(rotationMatrix3, 
     //			G4ThreeVector(0,0,-(halfPenSampleWidth+activePhotoCathodePMTThickness)),
     //			logicBoxActivePhotoCathodePMT, 
     //			"main_pmt_5", 
     //			logicWorldBox, false, 0, false);
     break;

     case 2:
     //physicPenSampleBox = new G4PVPlacement(0, G4ThreeVector(0,0,0),penLogicBox,"target",logicWorldBox,false,0,false);
     for(int iSample = 0; iSample < nSamples; iSample++){
  	physicPenStackedSamples = new G4PVPlacement(0, 
				G4ThreeVector(0,iSample*2*halfPenSampleThickness + iSample*5.*um,0),
				penLogicBox,
				"target_"+std::to_string(iSample+1),
				logicWorldBox,false,iSample,false);
     }
     //Reflector foil over PEN samples, uncomment if needed
     physicReflectorFoilBoxOverPEN = new G4PVPlacement(0, 
     				G4ThreeVector(0,halfReflectorBoxOverPEN + 2*halfPenSampleThickness*(nSamples)-halfPenSampleThickness + (nSamples)*5.*um,0),
     				logicReflectorFoilBox,
     				"reflector",
     				logicWorldBox,false,0,false);
     //Collimator, Trigger and light guide placement, with triggerFoilEJ212 PMT
     physicTriggerFoilEJ212 = new G4PVPlacement(0,
			G4ThreeVector(fDetectorCollimatorX,18*mm,0), 
			logicBoxTriggerFoilEJ212, 
			"triggerFoilEJ212",
			logicWorldBox, false, 0, false);
     if(reflectorOn){
     		physicReflectorFoilAroundEJ212Foil = new G4PVPlacement(0, 
			G4ThreeVector(fDetectorCollimatorX,18*mm,0), 
			logicReflectorFoilAroundEJ212Foil, 
			"foil", 
			logicWorldBox, false, 0, false);
     }
     physicCollimator = new G4PVPlacement(rotationMatrixCollimator, 
			G4ThreeVector(fDetectorCollimatorX,19*mm+halfCollimatorThickness,0), 
			logicCollimator, 
			"collimator", 
			logicWorldBox, false, 0, false);
     physicLightGuide = new G4PVPlacement(rotationMatrix4, 
			G4ThreeVector(-(halfLengthTriggerFoilEJ212+20*mm-0.4*mm)+fDetectorCollimatorX, (18*mm+14*mm-1.75*mm), 0),
			logicLightGuidePMMA,
			"LightGuide",
			logicWorldBox,false,0,false);
     physicOpticalGrease =  new G4PVPlacement(rotationMatrix5, 
			G4ThreeVector(-19*mm-(halfLengthTriggerFoilEJ212+20*mm)+fDetectorCollimatorX,(18*mm+14*mm-1.75*mm),0), 
			logicTubeSiGrease, 
			"grease", 
			logicWorldBox, false, 0, false);
     physicActivePhotoCathodeTriggerPMT = new G4PVPlacement(rotationMatrix,
			G4ThreeVector(-(20*mm+activePhotoCathodePMTThickness),0,0),
			logicBoxActivePhotoCathodePMT,
			"trigger_pmt",
			logicLightGuidePMMA,false,0,false);
     physicBoxPhotoCathodeSupport = new G4PVPlacement(0,
			G4ThreeVector(0,0,0),
			logicBoxPhotoCathodeSupport,
			"inactive_detector1",
			logicBoxActivePhotoCathodePMT,false,0,false);
     physicBoxPMTShell = new G4PVPlacement(0,
			G4ThreeVector(0,-(casingPMTLength + activePhotoCathodePMTThickness),0),
			logicBoxPMTShell,
			"pmt1", 
			logicBoxActivePhotoCathodePMT,false,0,false);
     physicBoxEmptyInsidePMT = new G4PVPlacement(0,
			G4ThreeVector(), 
			logicBoxEmptyInsidePMT, 
			"void1", 
			logicBoxPMTShell, false, 0, false);

     // Main PMT placements
     physicActivePhotoCathodePMT1 = new G4PVPlacement(rotationMatrix, 
			G4ThreeVector(-(halfPenSampleLength+activePhotoCathodePMTThickness), 0, 0), 
			logicBoxActivePhotoCathodePMT, 
			"main_pmt_1", 
			logicWorldBox, false, 0, false);
     physicActivePhotoCathodePMT2 = new G4PVPlacement(rotationMatrix1, 
			G4ThreeVector((halfPenSampleLength+activePhotoCathodePMTThickness), 0, 0), 
			logicBoxActivePhotoCathodePMT, 
			"main_pmt_2", 
			logicWorldBox, false, 0, false);
     physicActivePhotoCathodePMT3 = new G4PVPlacement(0, 
     			G4ThreeVector(0,-(halfPenSampleThickness+activePhotoCathodePMTThickness),0), 
     			logicBoxActivePhotoCathodePMT, 
     			"main_pmt_3", 
     			logicWorldBox, false, 0, false);
     physicActivePhotoCathodePMT4 = new G4PVPlacement(rotationMatrix2, 
     			G4ThreeVector(0,0,(halfPenSampleWidth+activePhotoCathodePMTThickness)),
     			logicBoxActivePhotoCathodePMT, 
     			"main_pmt_4", 
     			logicWorldBox, false, 0, false);
     physicActivePhotoCathodePMT5 = new G4PVPlacement(rotationMatrix3, 
     			G4ThreeVector(0,0,-(halfPenSampleWidth+activePhotoCathodePMTThickness)),
     			logicBoxActivePhotoCathodePMT, 
     			"main_pmt_5", 
     			logicWorldBox, false, 0, false);
     G4MTRunManager::GetRunManager()->GeometryHasBeenModified();
     break;



  }

  //============================================================= Surfaces =============================================================

  G4OpticalSurface* AirEJ212 = new G4OpticalSurface("AirEJ212", 
				glisur, 
				polished, 
				dielectric_dielectric);
  AirEJ212->SetPolish(0.99);

  G4OpticalSurface* AirPMMA = new G4OpticalSurface("AirPMMA", 
				glisur, 
				polished, 
				dielectric_dielectric);
  AirPMMA->SetPolish(1.0);

  G4OpticalSurface* PMMATrigger = new G4OpticalSurface("PMMATrigger", 
				glisur, 
				ground, 
				dielectric_dielectric);
  PMMATrigger->SetPolish(0.1);

  logicSurfaceEJ212Reflector = new G4LogicalBorderSurface("EJ212Reflector", 
				physicTriggerFoilEJ212, 
				physicReflectorFoilAroundEJ212Foil, 
				AirEJ212);

  logicSurfaceReflectorEJ212 = new G4LogicalBorderSurface("ReflectorEJ212", 
				physicWorldBox, 
				physicReflectorFoilAroundEJ212Foil, 
				AirEJ212);

  logicSurfacePenReflectorFoilBoxOverPEN = new G4LogicalBorderSurface("PenReflectorFoilBoxOverPEN", 
				physicPenSampleBox, 
				physicReflectorFoilBoxOverPEN, 
				AirEJ212);

  logicSurfacePENStackedReflectorFoilBoxOverPEN = new G4LogicalBorderSurface("PENStackedReflectorFoilBoxOverPEN", 
				physicPenStackedSamples, 
				physicReflectorFoilBoxOverPEN, 
				AirEJ212);

  logicSurfaceReflectorFoilBoxOverPENPENStacked = new G4LogicalBorderSurface("ReflectorFoilBoxOverPENPENStacked", 
				physicReflectorFoilBoxOverPEN, 
				physicPenStackedSamples, 
				AirEJ212);
  
  logicSurfacePENStackedAir = new G4LogicalBorderSurface("PENStackedAir", 
				physicPenStackedSamples, 
				physicWorldBox, 
				AirPEN);

  logicSurfaceAirPENStacked = new G4LogicalBorderSurface("AirPENStacked", 
				physicWorldBox, 
				physicPenStackedSamples, 
				AirPEN);

  logicSurfaceAirPMMA = new G4LogicalBorderSurface("AirPMMA", 
				physicLightGuide, 
				physicWorldBox, 
				AirPMMA);

  logicSurfacePMMAAir = new G4LogicalBorderSurface("PMMAAir", 
  				physicWorldBox, 
  				physicLightGuide, 
  				AirPMMA);

  logicSurfaceEJ212PMMA = new G4LogicalBorderSurface("TriggerPMMA", 
				physicLightGuide, 
				physicTriggerFoilEJ212, 
				PMMATrigger);

   
  logicSurfacePMMAEJ212 = new G4LogicalBorderSurface("PMMATrigger", 
				physicTriggerFoilEJ212, 
				physicLightGuide, 
				PMMATrigger);

  logicSurfaceAirPEN = new G4LogicalBorderSurface("AirPEN",
  				physicWorldBox,
  				physicPenSampleBox,
  				AirPEN);

  logicSurfacePENAir = new G4LogicalBorderSurface("PENAir",
                                physicPenSampleBox,
                                physicWorldBox,
                                AirPEN);

  return physicWorldBox;
}
