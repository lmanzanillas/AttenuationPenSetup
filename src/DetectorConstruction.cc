#include "DetectorConstruction.hh"
#include "DetectorMessenger.hh"
#include "SiliconPlateConstruction.hh"
#include "LightGuideConstruction.hh"
#include "PrimaryGeneratorAction.hh"

#include "G4MaterialTable.hh"
#include "PenMaterials.hh"
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
#include <iterator>

using namespace std;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

/*
Constructs DetectorConstruction, defines default values.
*/

DetectorConstruction::DetectorConstruction()
: G4VUserDetectorConstruction(),  penSampleBox(nullptr), penLogicBox(nullptr), 
physicPenSampleBox(nullptr),physicCollimator(nullptr),
AirTarget(nullptr),	
surfaceCathodeSupport(nullptr),	
MPT_PEN(nullptr)
{
  fDetectorMessenger = new DetectorMessenger(this);
  MPT_PEN = new G4MaterialPropertiesTable();
  halfSizeDarkBoxX = halfSizeDarkBoxY = halfSizeDarkBoxZ = 1.*m;
  fTargetName = "PEN";
  reflectorOn = false;
  //halfPenSampleLength = 15.*mm;
  halfPenSampleLength = 37.*mm;
  halfPenSampleThickness = 0.85*mm;
  halfPenSampleWidth = 4.5*mm;
  //halfPenSampleThickness = 5.0*mm;
  //halfPenSampleWidth = 15.*mm;
  halfCollimatorThickness = 5.0*mm;
  nSamples = 1;
  fDetectorType = 1;
  fDetectorCollimatorX = 0*mm;
  AbsorptionLength = 1.5;//value at 400 nm
  fRES = 1.0;
  fLY = 5500./MeV;
  userActivePhotoCathodeLength = 23.0/2.*mm;
  userActivePhotoCathodeWidth = 23.0/2.*mm;
  fDetectorName = "PenAttenuationSetup";
  fABSFile = "PEN_ABS";
  fVolName = "World";
  fSigAlpha = 1.00;
  fSigAlphaSides = 0.90;
  pmtReflectivity = 0.50;
  materialConstruction = new PenMaterials;
  DefineMaterials();
  fTargetMaterial = G4Material::GetMaterial("PVT_structure");
  fGlassMaterialPMT = G4Material::GetMaterial("BorosilicateGlass");
  MPT_GlassPMT = fGlassMaterialPMT->GetMaterialPropertiesTable();
  SetABS(AbsorptionLength);
  SetWorldMaterial("Air");
  //Construct();
  SetTargetMaterial("PVT_structure");
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
  G4RunManager::GetRunManager()->ReinitializeGeometry(); 
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
        userActivePhotoCathodeLength = 23.0/2. * mm;
        userActivePhotoCathodeWidth = 23.0/2. *mm;
  }
  else if(value == 1){
        userActivePhotoCathodeLength = 23.0/2. * mm;
        userActivePhotoCathodeWidth = 23.0/2. *mm;
  }
  else if(value == 2){
	halfPenSampleLength = 15. * mm;
	halfPenSampleWidth = 15. * mm;
        userActivePhotoCathodeLength = 23.0/2. * mm;
        userActivePhotoCathodeWidth = 23.0/2. *mm;
  }
  else if(value == 3){
  	nSamples = 1;
  	userActivePhotoCathodeLength = 6. /2. * mm;
  	userActivePhotoCathodeWidth = 4. /2. *mm;
  }
  else if(value == 4){
	halfPenSampleLength = 15. * mm;
	halfPenSampleWidth = 15. * mm;
        userActivePhotoCathodeLength = 23.0/2. * mm;
        userActivePhotoCathodeWidth = 23.0/2. *mm;
  }
  UpdateGeometry();
  G4RunManager::GetRunManager()->ReinitializeGeometry();
  //G4MTRunManager::GetRunManager()->PhysicsHasBeenModified();
}

void DetectorConstruction::SetNumberOfTargetSamples(G4int value){
  nSamples = value;
  //UpdateGeometry();
  G4RunManager::GetRunManager()->ReinitializeGeometry();
  //G4MTRunManager::GetRunManager()->PhysicsHasBeenModified();
}

//Sets dimmensions of target, thickness corresponds to the Y coordinate, Length to x.
void DetectorConstruction::SetTargetSampleLength(G4double value){
  halfPenSampleLength = (value/2.)*mm;
  //UpdateGeometry();
  G4RunManager::GetRunManager()->ReinitializeGeometry();
}

void DetectorConstruction::SetTargetSampleThickness(G4double value){
  halfPenSampleThickness = (value/2.)*mm;
  //UpdateGeometry();
  G4RunManager::GetRunManager()->ReinitializeGeometry();
}

void DetectorConstruction::SetTargetSampleWidth(G4double value){
  halfPenSampleWidth = (value/2.)*mm;
  //UpdateGeometry();
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
      ABSORPTION_PEN[absEntries] = (varAbsorLength*AbsorptionLength)*mm; //use measured value of attenuation to constrain curve and then change values multiplying the curve for a given factor
      absEntries++;
    }
  }

  else G4cout<<"Error opening file: " <<absFile<<G4endl;
  ReadAbs.close();
  absEntries--;

  const G4int nEntries1 = sizeof(wlPhotonEnergy)/sizeof(G4double);
  assert(sizeof(ABSORPTION_PEN) == sizeof(wlPhotonEnergy));
  MPT_PEN->AddProperty("ABSLENGTH",    wlPhotonEnergy, ABSORPTION_PEN, nEntries1)->SetSpline(true); // *
  G4MTRunManager::GetRunManager()->PhysicsHasBeenModified();
}

void DetectorConstruction::SetSigAlpha(G4double value){
  fSigAlpha=value;
  AirTarget -> SetSigmaAlpha(fSigAlpha);
  AirTarget -> SetMaterialPropertiesTable(MPT_Target);
   
  G4RunManager::GetRunManager()->ReinitializeGeometry();
  G4MTRunManager::GetRunManager()->PhysicsHasBeenModified();
}

void DetectorConstruction::SetSigAlphaSides(G4double value){
  fSigAlphaSides=value;
  surfaceCathodeSupport->SetSigmaAlpha(fSigAlphaSides);
  surfaceCathodeSupport->SetMaterialPropertiesTable(MPT_GlassPMT);
   
  G4RunManager::GetRunManager()->ReinitializeGeometry();
  G4MTRunManager::GetRunManager()->PhysicsHasBeenModified();
}
void DetectorConstruction::SetPMTReflectivity(G4double value){
  pmtReflectivity=value;
  G4RunManager::GetRunManager()->ReinitializeGeometry();
}
void DetectorConstruction::SetDetectorName(G4String name){
  fDetectorName=name;
}

void DetectorConstruction::SetABSFile(G4String fileName){
  fABSFile = fileName;
  G4MTRunManager::GetRunManager()->PhysicsHasBeenModified();
  DefineMaterials();
}


/*
Sets material of target.
*/
void DetectorConstruction::SetTargetMaterial(G4String materialChoice)
{
  // search the material by its name
  G4Material* pttoMaterial = G4Material::GetMaterial(materialChoice);

  if (pttoMaterial) {
    fTargetMaterial = pttoMaterial;
    fTargetName = fTargetMaterial->GetName();
    if(penLogicBox)penLogicBox->SetMaterial(fTargetMaterial);
    G4cout<<" light yield: "<<fTargetMaterial->GetMaterialPropertiesTable()->GetConstProperty("SCINTILLATIONYIELD")<<" photons/MeV"<<G4endl;  
    //G4cout<<" abs: "<<fTargetMaterial->GetMaterialPropertiesTable()->GetConstProperty("ABSLENGTH")<<" mm"<<G4endl;  
    //G4cout<<" surface: "<<fTargetMaterial->GetMaterialPropertiesTable()->GetConstProperty("ABSLENGTH")<<" mm"<<G4endl;  
  } else {
    G4cout << "\n--> warning from DetectorConstruction::SetMaterial : "
           << materialChoice << " not found" << G4endl;
  }
  //if(MPT_Target)MPT_Target->AddConstProperty("SCINTILLATIONYIELD",10000./MeV);
  G4RunManager::GetRunManager()->ReinitializeGeometry();
  //G4MTRunManager::GetRunManager()->PhysicsHasBeenModified();
}

void DetectorConstruction::SetRI(G4double value){
  fRI = value;
  UpdateGeometry();
  G4RunManager::GetRunManager()->ReinitializeGeometry();
  //G4MTRunManager::GetRunManager()->PhysicsHasBeenModified();
  //DefineMaterials();
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
  G4RunManager::GetRunManager()->ReinitializeGeometry();
  //G4MTRunManager::GetRunManager()->PhysicsHasBeenModified();
}

/*
Defines materials used in simulation. Sets material properties for PEN and other optical components.
*/
void DetectorConstruction::DefineMaterials(){
  // ============================================================= Materials =============================================================
  //materialConstruction = new PenMaterials;
  materialConstruction-> Construct();
  materialAir = G4Material::GetMaterial("Air");
  materialBialkali = G4Material::GetMaterial("Bialkali");
  fGlass = G4Material::GetMaterial("BorosilicateGlass");
  fPOM = G4Material::GetMaterial("POM");
  fABS = G4Material::GetMaterial("ABS");
  PenMaterial = G4Material::GetMaterial("PEN");
  materialSi = G4Material::GetMaterial("G4_Si");
  materialTriggerFoilEJ212 = G4Material::GetMaterial("EJ212");
  Pstyrene = G4Material::GetMaterial("Polystyrene");
  materialPMMA = G4Material::GetMaterial("PMMA");
  fVacuum = G4Material::GetMaterial("Vacuum");
  materialGreaseEJ550 = G4Material::GetMaterial("Grease");
  materialTeflon = G4Material::GetMaterial("G4_TEFLON");
  materialVikuiti = G4Material::GetMaterial("Vikuiti");
  materialTitanium = G4Material::GetMaterial("titanium");
  materialPolyethylene = G4Material::GetMaterial("G4_POLYETHYLENE");

  G4cout<<" materials ok "<<G4endl;

  G4double wavelength;
  char filler;
  G4double varAbsorLength;
  G4double emission;
  G4double rindex;

  G4double wlPhotonEnergy[102]  = {0};
  G4double ABSORPTION_PEN[102] = {0};
  G4double RINDEX_PEN[102] = {0};

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
      ABSORPTION_PEN[absEntries] = (varAbsorLength)*mm;
      RINDEX_PEN[absEntries] = rindex;
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
  
  MPT_PEN = new G4MaterialPropertiesTable();

  MPT_PEN->AddProperty("RINDEX",       wlPhotonEnergy, RINDEX_PEN, nEntries1)->SetSpline(true);
  MPT_PEN->AddProperty("ABSLENGTH",    wlPhotonEnergy, ABSORPTION_PEN, nEntries1)->SetSpline(true); // *

  // Read primary emission spectrum from PEN
  // Measurements from MPP Munich
  G4double pWavelength;
  G4String  Scint_file ="../properties/PEN_EM_SPECTRUM.dat";
  std::ifstream ReadScint2(Scint_file),ReadScintPEN;
  //count number of entries
  ReadScint2.unsetf(std::ios_base::skipws);
  //unsigned line_count = std::count(
  int line_count = std::count(
        std::istream_iterator<char>(ReadScint2),
        std::istream_iterator<char>(), 
        '\n');
  std::cout << "Lines: " << line_count << "\n";
  ReadScint2.close();
  G4double PEN_EMISSION[500]; 
  G4double PEN_WL_ENERGY[500]; 
  G4int nEntriesPEN = 0;
  ReadScintPEN.open(Scint_file);
  if(ReadScintPEN.is_open()){
        while(!ReadScintPEN.eof()){
                 ReadScintPEN>>pWavelength>>PEN_EMISSION[nEntriesPEN];
                 PEN_WL_ENERGY[nEntriesPEN] = (1240./pWavelength)*eV;//convert wavelength to eV
 	         //G4cout<<nEntriesPEN<<" wl "<<PEN_WL_ENERGY[nEntriesPEN]<<" "<<PEN_EMISSION[nEntriesPEN]<<G4endl;
                 nEntriesPEN++;
	         if(nEntriesPEN > (line_count-1)){ G4cout << " entries completed " << G4endl; break;}
        }
  }
  else
       G4cout << "Error opening file: " << Scint_file << G4endl;
  ReadScintPEN.close();
  G4cout<<" nEntriesPEN "<<nEntriesPEN<<G4endl;

  MPT_PEN->AddProperty("FASTCOMPONENT",PEN_WL_ENERGY, PEN_EMISSION, line_count)->SetSpline(true);
  MPT_PEN->AddProperty("SLOWCOMPONENT",PEN_WL_ENERGY, PEN_EMISSION, line_count)->SetSpline(true);

  MPT_PEN->AddConstProperty("SCINTILLATIONYIELD",fLY/MeV); // * 2.5 * PEN = PS, 10*PEN=PS
  MPT_PEN->AddConstProperty("RESOLUTIONSCALE",fRES); // * 1, 4, 8
  MPT_PEN->AddConstProperty("FASTTIMECONSTANT", 5.198*ns);
  MPT_PEN->AddConstProperty("SLOWTIMECONSTANT",24.336*ns);
  MPT_PEN->AddConstProperty("YIELDRATIO",1.);

  PenMaterial->SetMaterialPropertiesTable(MPT_PEN);
  //pvt_structure->SetMaterialPropertiesTable(MPT_PEN);

  
  G4cout<<" pen ok "<<G4endl;


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

  Scint_file="../input_files/pTP_emission.txt";
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
          Absorb[absorbEntries]=(varAbsorbLength)*m;
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
          wlsAbsorb[wlsAbsorbEntries]=varAbsorbLength*mm;
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

  //MPT_FoilEJ212->AddConstProperty("SCINTILLATIONYIELD",11520./MeV);
  MPT_FoilEJ212->AddConstProperty("SCINTILLATIONYIELD",10./MeV);//set low LY to make it faster, intead use Edep for coincidences
  MPT_FoilEJ212->AddConstProperty("RESOLUTIONSCALE",4.0);
  MPT_FoilEJ212->AddConstProperty("FASTTIMECONSTANT", 2.1*ns);
  MPT_FoilEJ212->AddConstProperty("SLOWTIMECONSTANT",14.2*ns);
  MPT_FoilEJ212->AddConstProperty("YIELDRATIO",1.0);

  materialTriggerFoilEJ212->SetMaterialPropertiesTable(MPT_FoilEJ212);

  G4cout<<" EJ212 ok "<<G4endl;

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
// ============================================================= Define Volumes =============================================================
  //The experimental Dark Box
  fWorldBox = new G4Box("World",halfSizeDarkBoxX,halfSizeDarkBoxY,halfSizeDarkBoxZ);
  logicWorldBox = new G4LogicalVolume(fWorldBox,materialAir,"World",0,0,0);
  physicWorldBox = new G4PVPlacement(0,G4ThreeVector(),logicWorldBox,"World",0,false,0);

  penSampleBox = new G4Box("target", halfPenSampleLength, halfPenSampleThickness, halfPenSampleWidth);
  penLogicBox = new G4LogicalVolume(penSampleBox,fTargetMaterial, "target",0,0,0);


  //EJ212 foil scintillator used for trigger
  double halfThicknessTriggerFoilEJ212 = 45.*um;//100 um thickness EJ foil for trigger
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
  G4LogicalVolume* logicReflectorFoilAroundEJ212Foil = new G4LogicalVolume(reflectorFoilAroundEJ212Foil, materialVikuiti, "foil", 0, 0, 0);

  //source holder for Bi source
  G4double innerRadiusSourceContainer = 0.*mm;
  G4double externalRadiusSourceContainer = 5.*mm;
  halfSourceContainerThickness = 53*um;
  G4Tubs* SourceContainerDisk = new G4Tubs("sourceContainer",innerRadiusSourceContainer,externalRadiusSourceContainer,halfSourceContainerThickness,0.,360.*deg);
  G4LogicalVolume* logicSourceContainer = new G4LogicalVolume(SourceContainerDisk,materialTitanium,"sourceContainer",0,0,0); 
 
  //Passive collimator
  double innerRadiusCollimator = 1.*mm; 
  double externalRadiusCollimator = 12.*mm; 
  halfCollimatorThickness= 10.*mm; 
  
  G4Tubs* collimatorTube = new G4Tubs("collimator",innerRadiusCollimator,externalRadiusCollimator,halfCollimatorThickness,0.,360.*deg);
  G4LogicalVolume* logicCollimator = new G4LogicalVolume(collimatorTube,materialPolyethylene,"collimator",0,0,0); 


  G4ThreeVector point = G4ThreeVector(0,0,5*cm);
  G4Navigator* pointNavigator = new G4Navigator();
  pointNavigator->SetWorldVolume(physicWorldBox);
  pointNavigator->LocateGlobalPointAndSetup(point);

  // ============================================================= Detectors =============================================================

  LightGuideConstruction guide;
  G4VSolid* LightGuidePMMA = guide.ConstructPlate();
  G4double halfLightGuideSizeX = guide.GetLighGuideSizeX();
  G4double halfLightGuideSizeY = guide.GetLighGuideSizeY();
  //G4double halfLightGuideSizeZ = guide.GetLighGuideSizeZ();
  G4double depthGreaseHole = guide.GetDepthGreaseHole();
  G4double diameterGreaseHole = guide.GetDiameterGreaseHole();
  G4double slitToPlaceThinSintillatorOffset = guide.GetSlitToPlaceThinSintillatorOffset();
  G4double slitToPlaceThinSintillatorDepth = guide.GetSlitToPlaceThinSintillatorDepth();
  //G4double slitToPlaceThinSintillatorHeight = guide.GetSlitToPlaceThinSintillatorHeight();
  G4cout<<" guide size x "<<guide.GetLighGuideSizeX()<<" slid depth "<<slitToPlaceThinSintillatorDepth<<G4endl;
  //G4cout <<  G4BestUnit(logicLightGuidePMMA->GetMass(true),"Mass") << G4endl;


  //Define PMTs geometry
  G4double casingPMTLength = 30.0*mm / 2.;
  G4double casingPMTWidth = 30.0*mm / 2.;
  G4double casingPMTHeight = 30.0*mm / 2.;
  G4double casingPMTHoleLength = 28.0*mm / 2.;
  G4double casingPMTHoleWidth = 28.0*mm / 2.;
  G4double casingPMTHoleHeight = 28.0*mm / 2.;

  G4double inactivePhotoCathodePMTLength = 26.2*mm / 2.;
  G4double inactivePhotoCathodePMTWidth = 26.2*mm / 2.;
  activePhotoCathodePMTLength = userActivePhotoCathodeLength;
  activePhotoCathodePMTWidth = userActivePhotoCathodeWidth;

  G4double inactivePhotoCathodePMTThickness = 0.8*mm / 2.;
  G4double activePhotoCathodePMTThickness = 0.1*mm / 2.;

  //PMT casing
  G4Box* boxCasingPMT = new G4Box("case",casingPMTLength,casingPMTWidth,casingPMTHeight);
  G4Box* boxEmptyInsidePMT = new G4Box("hole",casingPMTHoleLength,casingPMTHoleWidth,casingPMTHoleHeight);
  G4SubtractionSolid* boxPMTShell = new G4SubtractionSolid("boxPMTShell",boxCasingPMT,boxEmptyInsidePMT,0,G4ThreeVector(0,0,0));

  G4LogicalVolume* logicBoxEmptyInsidePMT = new G4LogicalVolume(boxEmptyInsidePMT, fVacuum,"vacuum");
  G4LogicalVolume* logicBoxPMTShell = new G4LogicalVolume(boxPMTShell,fPOM,"logicBoxPMTShell");


  //Photocathode 
  G4Box* boxTotalPhotoCathodeSupport = new G4Box("pmt_cathode",
		inactivePhotoCathodePMTLength,
		inactivePhotoCathodePMTThickness,
		inactivePhotoCathodePMTWidth);
  G4Box* boxActivePhotoCathodePMT = new G4Box("pmt_active",
		activePhotoCathodePMTLength,
		activePhotoCathodePMTThickness,
		activePhotoCathodePMTWidth);

  G4LogicalVolume* logicBoxPhotoCathodeSupport = new G4LogicalVolume(boxTotalPhotoCathodeSupport,
		fGlass,
		"logicBoxPhotoCathodeSupport");
  G4LogicalVolume* logicBoxActivePhotoCathodePMT = new G4LogicalVolume(boxActivePhotoCathodePMT,
		materialBialkali,
		"logicBoxActivePhotoCathodePMT");

  

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
  
   //Reflector foil to be placed over pen samples, not used in attenuation
  //Placement of trigger foil and others
  G4double offSetCollimator = 5.*mm;
  G4double offSetTriggerSetup = 5.*mm;
  double totalReflectorFoilThickness = 150.*um;
  double halfReflectorBoxOverPEN = (casingPMTHeight - 2*halfPenSampleThickness*nSamples + halfPenSampleThickness)/2.;
  //if( halfReflectorBoxOverPEN < 0.1){halfReflectorBoxOverPEN = casingPMTHeight - halfPenSampleThickness;}
  G4VSolid* boxReflectorPEN = new G4Box("foilPEN",halfPenSampleLength, halfReflectorBoxOverPEN, halfPenSampleWidth);
  G4VSolid* reflectorFoilBox = new G4Box("foil", halfPenSampleLength-totalReflectorFoilThickness, halfReflectorBoxOverPEN, halfPenSampleWidth-totalReflectorFoilThickness);
  G4SubtractionSolid* reflectorFoilOverPEN = new G4SubtractionSolid("reflectorFoilOverPEN",boxReflectorPEN,reflectorFoilBox,0,G4ThreeVector(0, totalReflectorFoilThickness, 0));
  G4LogicalVolume* logicReflectorFoilBox = new G4LogicalVolume(reflectorFoilOverPEN, materialVikuiti, "foil", 0, 0, 0);


  //-(halfLengthTriggerFoilEJ212+20*mm-0.4*mm)+fDetectorCollimatorX
  G4double xPositionLightGuide = fDetectorCollimatorX - halfLengthTriggerFoilEJ212 - halfLightGuideSizeX + slitToPlaceThinSintillatorDepth;
  G4double yPositionLightGuide = casingPMTHeight  + offSetTriggerSetup + halfLightGuideSizeY;
  G4double xPositionGrease = xPositionLightGuide - halfLightGuideSizeX + depthGreaseHole/2.;
  G4double positionYTriggerFoil = yPositionLightGuide - halfLightGuideSizeY + slitToPlaceThinSintillatorOffset;
  G4cout<<" positionYTriggerFoil "<<positionYTriggerFoil<<G4endl;
  G4double xPositionPMTPhotoCathode = xPositionLightGuide-halfLightGuideSizeX -2*inactivePhotoCathodePMTThickness - activePhotoCathodePMTThickness;
  //position of PMT with respect to Light Guide
  G4double xPositionInactivePMTPhotoCathode = -halfLightGuideSizeX -inactivePhotoCathodePMTThickness;
  fDetectorCollimatorY = positionYTriggerFoil + halfCollimatorThickness + offSetCollimator;
  

  G4double SpaceBetweenSamples = 5.0*um;

  if(fDetectorType == 0){fSourceContainerY = positionYTriggerFoil + offSetCollimator + halfSourceContainerThickness;}
  else if (fDetectorType == 1){fSourceContainerY = positionYTriggerFoil + 2*halfCollimatorThickness + offSetCollimator + halfSourceContainerThickness;} 
  else if (fDetectorType == 2){fSourceContainerY = positionYTriggerFoil + offSetCollimator + halfSourceContainerThickness;} 
  else if (fDetectorType == 3){fSourceContainerY = positionYTriggerFoil + offSetCollimator + halfSourceContainerThickness;} 
  else if (fDetectorType == 4){fSourceContainerY = positionYTriggerFoil + 2*halfCollimatorThickness + offSetCollimator + halfSourceContainerThickness;} 

  
  // Set Draw G4VisAttributes
  G4VisAttributes* visAttr = new G4VisAttributes();
  visAttr->SetVisibility(false);
  logicWorldBox->SetVisAttributes(visAttr);

  //PEN and EJ212 
  G4VisAttributes* visualAttributesScintillators = new G4VisAttributes(G4Colour::Blue());
  visualAttributesScintillators->SetVisibility(true);
  penLogicBox->SetVisAttributes(visualAttributesScintillators);
  logicBoxTriggerFoilEJ212->SetVisAttributes(visualAttributesScintillators);

  //Collimator
  G4VisAttributes* visAttributesCollimator = new G4VisAttributes(G4Colour::White());
  visAttributesCollimator->SetVisibility(true);
  visAttributesCollimator->SetForceSolid(true);
  visAttributesCollimator->SetForceAuxEdgeVisible(true);
  visAttributesCollimator->SetLineWidth(3.0);
  visAttributesCollimator->SetForceLineSegmentsPerCircle(50);
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
  //Source container
  G4VisAttributes* visualAttributesSource = new G4VisAttributes(G4Colour::Gray());
  visualAttributesSource->SetVisibility(true);
  visualAttributesSource->SetForceSolid(true);
  visualAttributesSource->SetForceAuxEdgeVisible(true);
  logicSourceContainer->SetVisAttributes(visualAttributesSource);
  
   
  // Inactive volumes
  G4VisAttributes* visualAttributesInactiveMaterials = new G4VisAttributes(G4Colour::Gray());
  visualAttributesInactiveMaterials->SetVisibility(true);
  visualAttributesInactiveMaterials->SetForceSolid(false);
  visualAttributesInactiveMaterials->SetForceAuxEdgeVisible(true);
  //logicBoxEmptyInsidePMT->SetVisAttributes(visualAttributesInactiveMaterials);
  logicBoxPhotoCathodeSupport->SetVisAttributes(visualAttributesInactiveMaterials);
  G4VisAttributes* visualAttributesInactive = new G4VisAttributes(G4Colour::Black());
  visualAttributesInactive->SetVisibility(true);
  visualAttributesInactive->SetForceSolid(true);
  visualAttributesInactive->SetForceAuxEdgeVisible(true);
  logicBoxPMTShell->SetVisAttributes(visualAttributesInactive);


  G4LogicalVolume* logicLightGuidePMMA = new G4LogicalVolume(LightGuidePMMA, materialPMMA, "Ligh_guideLog");
  G4VisAttributes* visualAttributesLightGuide = new G4VisAttributes(G4Colour::Gray());
  visualAttributesLightGuide->SetVisibility(true);
  //visualAttributesLightGuide->SetForceSolid(true);
  logicLightGuidePMMA->SetVisAttributes(visualAttributesLightGuide);

  G4double greaseHoleDiameter = diameterGreaseHole*mm;
  G4double greaseHoleDepth = depthGreaseHole*mm;
  G4Tubs* tubeSiGrease = new G4Tubs("SiGrease", 0, greaseHoleDiameter/2., greaseHoleDepth/2., 0, 360*deg);

  G4LogicalVolume* logicTubeSiGrease = new G4LogicalVolume(tubeSiGrease, materialGreaseEJ550, "logicTubeSiGrease");
  logicTubeSiGrease->SetVisAttributes(visulaAttributesReflectors);


  //  ============================================================= Place volumes =============================================================
  
  // Place main tile at centre of world volume
  switch(fDetectorType){

  case 0:
    for(int iSample = 0; iSample < nSamples; iSample++){
  	physicPenStackedSamples = new G4PVPlacement(0, 
				G4ThreeVector(0,iSample*2*halfPenSampleThickness + iSample*SpaceBetweenSamples,0),
				penLogicBox,
				"target_"+std::to_string(iSample+1),
				logicWorldBox,false,iSample,false);
     }

     // Trigger and light guide placement, with triggerFoilEJ212 PMT
     if(reflectorOn){
     		physicReflectorFoilAroundEJ212Foil = new G4PVPlacement(0, 
			G4ThreeVector(0.*mm + fDetectorCollimatorX,positionYTriggerFoil,0), 
			logicReflectorFoilAroundEJ212Foil, 
			"foil", 
			logicWorldBox, false, 0, false);
     }
     //Trigger foil
     physicSourceContainer = new G4PVPlacement(rotationMatrixCollimator,
                        G4ThreeVector(0*mm + fDetectorCollimatorX,fSourceContainerY,0),
                        logicSourceContainer,
                        "collimator",
                        logicWorldBox, false, 0, false);
     physicTriggerFoilEJ212 = new G4PVPlacement(0, 
			G4ThreeVector(0.*mm + fDetectorCollimatorX,positionYTriggerFoil,0), 
			logicBoxTriggerFoilEJ212, 
			"triggerFoilEJ212", 
			logicWorldBox, false, 0, false);
     physicLightGuide = new G4PVPlacement(rotationMatrix4, 
			G4ThreeVector(xPositionLightGuide,yPositionLightGuide, 0),
			logicLightGuidePMMA,
			"LightGuide",
			logicWorldBox,false,0,false);
     physicOpticalGrease =  new G4PVPlacement(rotationMatrix5, 
			G4ThreeVector(xPositionGrease,yPositionLightGuide,0), 
			logicTubeSiGrease, 
			"grease", 
			logicWorldBox, false, 0, false);
     //PMT for trigger
     //Use the log active photocathode as mother of the PMT and then just place it according to the setup configuration
      physicBoxPhotoCathodeSupport = new G4PVPlacement(rotationMatrix,
                        G4ThreeVector(xPositionInactivePMTPhotoCathode,0,0),
                        logicBoxPhotoCathodeSupport,
                        "inactive_detector1",
                        logicLightGuidePMMA,false,0,false);

      physicActivePhotoCathodeTriggerPMT = new G4PVPlacement(rotationMatrix,
                        G4ThreeVector(xPositionPMTPhotoCathode,yPositionLightGuide,0),
                        logicBoxActivePhotoCathodePMT,
                        "trigger_pmt",
                        logicWorldBox,false,0,false);

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
     // PMT1
     physicActivePhotoCathodePMT1 = new G4PVPlacement(rotationMatrix, 
			G4ThreeVector(-(halfPenSampleLength + activePhotoCathodePMTThickness + 2*inactivePhotoCathodePMTThickness), 0, 0), 
			logicBoxActivePhotoCathodePMT, 
			"main_pmt_1", 
			logicWorldBox, false, 0, false);
     physicPhotoCathodeSupportPMT1 = new G4PVPlacement(rotationMatrix,
			G4ThreeVector(-(halfPenSampleLength + inactivePhotoCathodePMTThickness), 0, 0), 
			logicBoxPhotoCathodeSupport, 
			"support1", 
			logicWorldBox, false, 0, false);
     //PMT2
     physicActivePhotoCathodePMT2 = new G4PVPlacement(rotationMatrix1, 
			G4ThreeVector((halfPenSampleLength + activePhotoCathodePMTThickness + 2*inactivePhotoCathodePMTThickness), 0, 0), 
			logicBoxActivePhotoCathodePMT, 
			"main_pmt_2", 
			logicWorldBox, false, 0, false);
     physicPhotoCathodeSupportPMT2 = new G4PVPlacement(rotationMatrix1,
			G4ThreeVector((halfPenSampleLength + inactivePhotoCathodePMTThickness), 0, 0), 
			logicBoxPhotoCathodeSupport, 
			"support2", 
			logicWorldBox, false, 0, false);
     //PMT3
     physicActivePhotoCathodePMT3 = new G4PVPlacement(0, 
     			G4ThreeVector(0,-(halfPenSampleThickness + activePhotoCathodePMTThickness + 2*inactivePhotoCathodePMTThickness),0), 
     			logicBoxActivePhotoCathodePMT, 
     			"main_pmt_3", 
     			logicWorldBox, false, 0, false);
     physicPhotoCathodeSupportPMT3 = new G4PVPlacement(0,
			G4ThreeVector(0,-(halfPenSampleThickness + inactivePhotoCathodePMTThickness), 0), 
			logicBoxPhotoCathodeSupport, 
			"support3", 
			logicWorldBox, false, 0, false);
     //PMT4
     physicActivePhotoCathodePMT4 = new G4PVPlacement(rotationMatrix2, 
     			G4ThreeVector(0,0,(halfPenSampleWidth + activePhotoCathodePMTThickness + 2*inactivePhotoCathodePMTThickness)),
     			logicBoxActivePhotoCathodePMT, 
     			"main_pmt_4", 
     			logicWorldBox, false, 0, false);
     physicPhotoCathodeSupportPMT4 = new G4PVPlacement(rotationMatrix2,
			G4ThreeVector(0,0,(halfPenSampleWidth + inactivePhotoCathodePMTThickness)), 
			logicBoxPhotoCathodeSupport, 
			"support4", 
			logicWorldBox, false, 0, false);
     //PMT5
     physicActivePhotoCathodePMT5 = new G4PVPlacement(rotationMatrix3, 
     			G4ThreeVector(0,0,-(halfPenSampleWidth + activePhotoCathodePMTThickness + 2*inactivePhotoCathodePMTThickness)),
     			logicBoxActivePhotoCathodePMT, 
     			"main_pmt_5", 
     			logicWorldBox, false, 0, false);
     physicPhotoCathodeSupportPMT5 = new G4PVPlacement(rotationMatrix3,
			G4ThreeVector(0,0,-(halfPenSampleWidth + inactivePhotoCathodePMTThickness)), 
			logicBoxPhotoCathodeSupport, 
			"support5", 
			logicWorldBox, false, 0, false);
     break;

    

  case 1:
     //Attenuation setup consisting on pen samples coupled to two PMTs and the trigger setup that can be moved
     for(int iSample = 0; iSample < nSamples; iSample++){
  	physicPenStackedSamples = new G4PVPlacement(0, 
				G4ThreeVector(0,iSample*2*halfPenSampleThickness+iSample*SpaceBetweenSamples,0),
				penLogicBox,
				"target_"+std::to_string(iSample+1),
				logicWorldBox,false,iSample,false);
     }
     //Collimator, Trigger and light guide placement, with triggerFoilEJ212 PMT
     physicTriggerFoilEJ212 = new G4PVPlacement(0,
			G4ThreeVector(0.*mm + fDetectorCollimatorX,positionYTriggerFoil,0), 
			logicBoxTriggerFoilEJ212, 
			"triggerFoilEJ212",
			logicWorldBox, false, 0, false);
     if(reflectorOn){
     		physicReflectorFoilAroundEJ212Foil = new G4PVPlacement(0, 
			G4ThreeVector(0.*mm + fDetectorCollimatorX,positionYTriggerFoil,0), 
			logicReflectorFoilAroundEJ212Foil, 
			"foil", 
			logicWorldBox, false, 0, false);
     }
     physicCollimator = new G4PVPlacement(rotationMatrixCollimator, 
			G4ThreeVector(0*mm + fDetectorCollimatorX,fDetectorCollimatorY,0), 
			logicCollimator, 
			"collimator", 
			logicWorldBox, false, 0, false);
     physicSourceContainer = new G4PVPlacement(rotationMatrixCollimator,
                        G4ThreeVector(0*mm + fDetectorCollimatorX,fSourceContainerY,0),
                        logicSourceContainer,
                        "collimator",
                        logicWorldBox, false, 0, false);
     physicLightGuide = new G4PVPlacement(rotationMatrix4, 
			G4ThreeVector(xPositionLightGuide,yPositionLightGuide, 0),
			logicLightGuidePMMA,
			"LightGuide",
			logicWorldBox,false,0,false);
     physicOpticalGrease =  new G4PVPlacement(rotationMatrix5, 
			G4ThreeVector(xPositionGrease,yPositionLightGuide,0), 
			logicTubeSiGrease, 
			"grease", 
			logicWorldBox, false, 0, false);
     //PMT for trigger
     physicBoxPhotoCathodeSupport = new G4PVPlacement(rotationMatrix,
                        G4ThreeVector(xPositionInactivePMTPhotoCathode,0,0),
                        logicBoxPhotoCathodeSupport,
                        "inactive_detector1",
                        logicLightGuidePMMA,false,0,false);

      physicActivePhotoCathodeTriggerPMT = new G4PVPlacement(rotationMatrix,
                        G4ThreeVector(xPositionPMTPhotoCathode,yPositionLightGuide,0),
                        logicBoxActivePhotoCathodePMT,
                        "trigger_pmt",
                        logicWorldBox,false,0,false);

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
     // PMT1
     physicActivePhotoCathodePMT1 = new G4PVPlacement(rotationMatrix,
                        G4ThreeVector(-(halfPenSampleLength + activePhotoCathodePMTThickness + 2*inactivePhotoCathodePMTThickness), 0, 0),
                        logicBoxActivePhotoCathodePMT,
                        "main_pmt_1",
                        logicWorldBox, false, 0, false);
     physicPhotoCathodeSupportPMT1 = new G4PVPlacement(rotationMatrix,
                        G4ThreeVector(-(halfPenSampleLength + inactivePhotoCathodePMTThickness), 0, 0),
                        logicBoxPhotoCathodeSupport,
                        "support1",
                        logicWorldBox, false, 0, false);
     //PMT2
     physicActivePhotoCathodePMT2 = new G4PVPlacement(rotationMatrix1,
                        G4ThreeVector((halfPenSampleLength + activePhotoCathodePMTThickness + 2*inactivePhotoCathodePMTThickness), 0, 0),
                        logicBoxActivePhotoCathodePMT,
                        "main_pmt_2",
                        logicWorldBox, false, 0, false);
     physicPhotoCathodeSupportPMT2 = new G4PVPlacement(rotationMatrix1,
                        G4ThreeVector((halfPenSampleLength + inactivePhotoCathodePMTThickness), 0, 0),
                        logicBoxPhotoCathodeSupport,
                        "support2",
                        logicWorldBox, false, 0, false);

     break;

     case 2:
     for(int iSample = 0; iSample < nSamples; iSample++){
  	physicPenStackedSamples = new G4PVPlacement(0, 
				G4ThreeVector(0,iSample*2*halfPenSampleThickness + iSample*SpaceBetweenSamples,0),
				penLogicBox,
				"target_"+std::to_string(iSample+1),
				logicWorldBox,false,iSample,false);
     }
     //Reflector foil over PEN samples, uncomment if needed
     
     physicReflectorFoilBoxOverPEN = new G4PVPlacement(0, 
     				G4ThreeVector(0,halfReflectorBoxOverPEN + 2*halfPenSampleThickness*(nSamples)-halfPenSampleThickness + (nSamples)*SpaceBetweenSamples,0),
     				logicReflectorFoilBox,
     				"reflector",
     				logicWorldBox,false,0,false);
     
     //Collimator, Trigger and light guide placement, with triggerFoilEJ212 PMT
     physicTriggerFoilEJ212 = new G4PVPlacement(0,
			G4ThreeVector(0.*mm + fDetectorCollimatorX,positionYTriggerFoil,0), 
			logicBoxTriggerFoilEJ212, 
			"triggerFoilEJ212",
			logicWorldBox, false, 0, false);
     if(reflectorOn){
     		physicReflectorFoilAroundEJ212Foil = new G4PVPlacement(0, 
			G4ThreeVector(0.*mm + fDetectorCollimatorX,positionYTriggerFoil,0), 
			logicReflectorFoilAroundEJ212Foil, 
			"foil", 
			logicWorldBox, false, 0, false);
     }
     //physicCollimator = new G4PVPlacement(rotationMatrixCollimator, 
     //			G4ThreeVector(fDetectorCollimatorX,19*mm+halfCollimatorThickness,0), 
     //			logicCollimator, 
     //			"collimator", 
     //			logicWorldBox, false, 0, false);
     physicSourceContainer = new G4PVPlacement(rotationMatrixCollimator,
                        G4ThreeVector(0*mm + fDetectorCollimatorX,fSourceContainerY,0),
                        logicSourceContainer,
                        "collimator",
                        logicWorldBox, false, 0, false);
     physicLightGuide = new G4PVPlacement(rotationMatrix4, 
			G4ThreeVector(xPositionLightGuide,yPositionLightGuide, 0),
			logicLightGuidePMMA,
			"LightGuide",
			logicWorldBox,false,0,false);
     physicOpticalGrease =  new G4PVPlacement(rotationMatrix5, 
			G4ThreeVector(xPositionGrease,yPositionLightGuide,0), 
			logicTubeSiGrease, 
			"grease", 
			logicWorldBox, false, 0, false);
     //trigger pmt
     physicBoxPhotoCathodeSupport = new G4PVPlacement(rotationMatrix,
                        G4ThreeVector(xPositionInactivePMTPhotoCathode,0,0),
                        logicBoxPhotoCathodeSupport,
                        "inactive_detector1",
                        logicLightGuidePMMA,false,0,false);

     physicActivePhotoCathodeTriggerPMT = new G4PVPlacement(rotationMatrix,
                        G4ThreeVector(xPositionPMTPhotoCathode,yPositionLightGuide,0),
                        logicBoxActivePhotoCathodePMT,
                        "trigger_pmt",
                        logicWorldBox,false,0,false);

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
     // PMT1
     physicActivePhotoCathodePMT1 = new G4PVPlacement(rotationMatrix, 
			G4ThreeVector(-(halfPenSampleLength + activePhotoCathodePMTThickness + 2*inactivePhotoCathodePMTThickness), 0, 0), 
			logicBoxActivePhotoCathodePMT, 
			"main_pmt_1", 
			logicWorldBox, false, 0, false);
     physicPhotoCathodeSupportPMT1 = new G4PVPlacement(rotationMatrix,
			G4ThreeVector(-(halfPenSampleLength + inactivePhotoCathodePMTThickness), 0, 0), 
			logicBoxPhotoCathodeSupport, 
			"support1", 
			logicWorldBox, false, 0, false);
     //PMT2
     physicActivePhotoCathodePMT2 = new G4PVPlacement(rotationMatrix1, 
			G4ThreeVector((halfPenSampleLength + activePhotoCathodePMTThickness + 2*inactivePhotoCathodePMTThickness), 0, 0), 
			logicBoxActivePhotoCathodePMT, 
			"main_pmt_2", 
			logicWorldBox, false, 0, false);
     physicPhotoCathodeSupportPMT2 = new G4PVPlacement(rotationMatrix1,
			G4ThreeVector((halfPenSampleLength + inactivePhotoCathodePMTThickness), 0, 0), 
			logicBoxPhotoCathodeSupport, 
			"support2", 
			logicWorldBox, false, 0, false);
     //PMT3
     physicActivePhotoCathodePMT3 = new G4PVPlacement(0, 
     			G4ThreeVector(0,-(halfPenSampleThickness + activePhotoCathodePMTThickness + 2*inactivePhotoCathodePMTThickness),0), 
     			logicBoxActivePhotoCathodePMT, 
     			"main_pmt_3", 
     			logicWorldBox, false, 0, false);
     physicPhotoCathodeSupportPMT3 = new G4PVPlacement(0,
			G4ThreeVector(0,-(halfPenSampleThickness + inactivePhotoCathodePMTThickness), 0), 
			logicBoxPhotoCathodeSupport, 
			"support3", 
			logicWorldBox, false, 0, false);
     //PMT4
     physicActivePhotoCathodePMT4 = new G4PVPlacement(rotationMatrix2, 
     			G4ThreeVector(0,0,(halfPenSampleWidth + activePhotoCathodePMTThickness + 2*inactivePhotoCathodePMTThickness)),
     			logicBoxActivePhotoCathodePMT, 
     			"main_pmt_4", 
     			logicWorldBox, false, 0, false);
     physicPhotoCathodeSupportPMT4 = new G4PVPlacement(rotationMatrix2,
			G4ThreeVector(0,0,(halfPenSampleWidth + inactivePhotoCathodePMTThickness)), 
			logicBoxPhotoCathodeSupport, 
			"support4", 
			logicWorldBox, false, 0, false);
     //PMT5
     physicActivePhotoCathodePMT5 = new G4PVPlacement(rotationMatrix3, 
     			G4ThreeVector(0,0,-(halfPenSampleWidth + activePhotoCathodePMTThickness + 2*inactivePhotoCathodePMTThickness)),
     			logicBoxActivePhotoCathodePMT, 
     			"main_pmt_5", 
     			logicWorldBox, false, 0, false);
     physicPhotoCathodeSupportPMT5 = new G4PVPlacement(rotationMatrix3,
			G4ThreeVector(0,0,-(halfPenSampleWidth + inactivePhotoCathodePMTThickness)), 
			logicBoxPhotoCathodeSupport, 
			"support5", 
			logicWorldBox, false, 0, false);
     break;

    case 3:
    for(int iSample = 0; iSample < nSamples; iSample++){
  	physicPenStackedSamples = new G4PVPlacement(0, 
				G4ThreeVector(0,iSample*2*halfPenSampleThickness + iSample*5.*um,0),
				penLogicBox,
				"target_"+std::to_string(iSample+1),
				logicWorldBox,false,iSample,false);
     }

     //Setup to reproduce the spectrometer measurements, it consit of only 1 PMT

     physicBoxPhotoCathodeSupport = new G4PVPlacement(0,
                        G4ThreeVector(0,(activePhotoCathodePMTThickness + inactivePhotoCathodePMTThickness),0),
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
			G4ThreeVector(-(halfPenSampleLength + activePhotoCathodePMTThickness + 2*inactivePhotoCathodePMTThickness),0,0), 
			logicBoxActivePhotoCathodePMT, 
			"main_pmt_1", 
			logicWorldBox, false, 0, false);
     break;

    case 4:
    for(int iSample = 0; iSample < nSamples; iSample++){
  	physicPenStackedSamples = new G4PVPlacement(0, 
				G4ThreeVector(0,iSample*2*halfPenSampleThickness + iSample*SpaceBetweenSamples,0),
				penLogicBox,
				"target_"+std::to_string(iSample+1),
				logicWorldBox,false,iSample,false);
     }
     //Reflector foil over PEN samples, uncomment if needed
     
     physicReflectorFoilBoxOverPEN = new G4PVPlacement(0, 
     				G4ThreeVector(0,halfReflectorBoxOverPEN + 2*halfPenSampleThickness*(nSamples)-halfPenSampleThickness + (nSamples)*SpaceBetweenSamples,0),
     				logicReflectorFoilBox,
     				"reflector",
     				logicWorldBox,false,0,false);
     
     //Collimator, Trigger and light guide placement, with triggerFoilEJ212 PMT
     physicTriggerFoilEJ212 = new G4PVPlacement(0,
			G4ThreeVector(0.*mm + fDetectorCollimatorX,positionYTriggerFoil,0), 
			logicBoxTriggerFoilEJ212, 
			"triggerFoilEJ212",
			logicWorldBox, false, 0, false);
     if(reflectorOn){
     		physicReflectorFoilAroundEJ212Foil = new G4PVPlacement(0, 
			G4ThreeVector(0.*mm + fDetectorCollimatorX,positionYTriggerFoil,0), 
			logicReflectorFoilAroundEJ212Foil, 
			"foil", 
			logicWorldBox, false, 0, false);
     }
     physicCollimator = new G4PVPlacement(rotationMatrixCollimator, 
     			G4ThreeVector(0*mm + fDetectorCollimatorX,fDetectorCollimatorY,0), 
     			logicCollimator, 
     			"collimator", 
     			logicWorldBox, false, 0, false);
    
     physicSourceContainer = new G4PVPlacement(rotationMatrixCollimator, 
     			G4ThreeVector(0*mm + fDetectorCollimatorX,fSourceContainerY,0), 
     			logicSourceContainer, 
     			"collimator", 
     			logicWorldBox, false, 0, false);


     physicLightGuide = new G4PVPlacement(rotationMatrix4, 
			G4ThreeVector(xPositionLightGuide,yPositionLightGuide, 0),
			logicLightGuidePMMA,
			"LightGuide",
			logicWorldBox,false,0,false);
     physicOpticalGrease =  new G4PVPlacement(rotationMatrix5, 
			G4ThreeVector(xPositionGrease,yPositionLightGuide,0), 
			logicTubeSiGrease, 
			"grease", 
			logicWorldBox, false, 0, false);
     //Trigger pmt
     physicBoxPhotoCathodeSupport = new G4PVPlacement(rotationMatrix,
                        G4ThreeVector(xPositionInactivePMTPhotoCathode,0,0),
                        logicBoxPhotoCathodeSupport,
                        "inactive_detector1",
                        logicLightGuidePMMA,false,0,false);

     physicActivePhotoCathodeTriggerPMT = new G4PVPlacement(rotationMatrix,
                        G4ThreeVector(xPositionPMTPhotoCathode,yPositionLightGuide,0),
                        logicBoxActivePhotoCathodePMT,
                        "trigger_pmt",
                        logicWorldBox,false,0,false);

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
     // PMT1
     physicActivePhotoCathodePMT1 = new G4PVPlacement(rotationMatrix, 
			G4ThreeVector(-(halfPenSampleLength + activePhotoCathodePMTThickness + 2*inactivePhotoCathodePMTThickness), 0, 0), 
			logicBoxActivePhotoCathodePMT, 
			"main_pmt_1", 
			logicWorldBox, false, 0, false);
     physicPhotoCathodeSupportPMT1 = new G4PVPlacement(rotationMatrix,
			G4ThreeVector(-(halfPenSampleLength + inactivePhotoCathodePMTThickness), 0, 0), 
			logicBoxPhotoCathodeSupport, 
			"support1", 
			logicWorldBox, false, 0, false);
     //PMT2
     physicActivePhotoCathodePMT2 = new G4PVPlacement(rotationMatrix1, 
			G4ThreeVector((halfPenSampleLength + activePhotoCathodePMTThickness + 2*inactivePhotoCathodePMTThickness), 0, 0), 
			logicBoxActivePhotoCathodePMT, 
			"main_pmt_2", 
			logicWorldBox, false, 0, false);
     physicPhotoCathodeSupportPMT2 = new G4PVPlacement(rotationMatrix1,
			G4ThreeVector((halfPenSampleLength + inactivePhotoCathodePMTThickness), 0, 0), 
			logicBoxPhotoCathodeSupport, 
			"support2", 
			logicWorldBox, false, 0, false);
     //PMT3
     physicActivePhotoCathodePMT3 = new G4PVPlacement(0, 
     			G4ThreeVector(0,-(halfPenSampleThickness + activePhotoCathodePMTThickness + 2*inactivePhotoCathodePMTThickness),0), 
     			logicBoxActivePhotoCathodePMT, 
     			"main_pmt_3", 
     			logicWorldBox, false, 0, false);
     physicPhotoCathodeSupportPMT3 = new G4PVPlacement(0,
			G4ThreeVector(0,-(halfPenSampleThickness + inactivePhotoCathodePMTThickness), 0), 
			logicBoxPhotoCathodeSupport, 
			"support3", 
			logicWorldBox, false, 0, false);
     //PMT4
     physicActivePhotoCathodePMT4 = new G4PVPlacement(rotationMatrix2, 
     			G4ThreeVector(0,0,(halfPenSampleWidth + activePhotoCathodePMTThickness + 2*inactivePhotoCathodePMTThickness)),
     			logicBoxActivePhotoCathodePMT, 
     			"main_pmt_4", 
     			logicWorldBox, false, 0, false);
     physicPhotoCathodeSupportPMT4 = new G4PVPlacement(rotationMatrix2,
			G4ThreeVector(0,0,(halfPenSampleWidth + inactivePhotoCathodePMTThickness)), 
			logicBoxPhotoCathodeSupport, 
			"support4", 
			logicWorldBox, false, 0, false);
     //PMT5
     physicActivePhotoCathodePMT5 = new G4PVPlacement(rotationMatrix3, 
     			G4ThreeVector(0,0,-(halfPenSampleWidth + activePhotoCathodePMTThickness + 2*inactivePhotoCathodePMTThickness)),
     			logicBoxActivePhotoCathodePMT, 
     			"main_pmt_5", 
     			logicWorldBox, false, 0, false);
     physicPhotoCathodeSupportPMT5 = new G4PVPlacement(rotationMatrix3,
			G4ThreeVector(0,0,-(halfPenSampleWidth + inactivePhotoCathodePMTThickness)), 
			logicBoxPhotoCathodeSupport, 
			"support5", 
			logicWorldBox, false, 0, false);
     break;
  }

  //============================================================= Surfaces =============================================================

  //Main target, set the sigma alpha here, get the material table of the target and include the sigma alpha
  // ==== AirTarget surfaces ====================

  AirTarget = new G4OpticalSurface("AirTarget",unified, ground, dielectric_dielectric);
  AirTarget -> SetSigmaAlpha(fSigAlpha);
  MPT_Target = new G4MaterialPropertiesTable();
  MPT_Target = fTargetMaterial->GetMaterialPropertiesTable();
  //G4MaterialPropertyVector *check_values = MPT_Target->GetProperty("RAYLEIGH");
  //check_values->DumpValues ();
  AirTarget -> SetMaterialPropertiesTable(MPT_Target);

  logicSurfacePENStackedAir = new G4LogicalBorderSurface("AirTarget", 
				physicPenStackedSamples, 
				physicWorldBox, 
				AirTarget);

  logicSurfaceAirTargetStacked = new G4LogicalBorderSurface("AirTarget", 
				physicWorldBox, 
				physicPenStackedSamples, 
				AirTarget);

  logicSurfacePENStackedPENStacked = new G4LogicalBorderSurface("AirTarget",
                                physicPenStackedSamples, 
                                physicPenStackedSamples,
                                AirTarget);

  logicSurfaceAirTarget = new G4LogicalBorderSurface("AirTarget",
  				physicWorldBox,
  				physicPenSampleBox,
  				AirTarget);

  logicSurfaceTargetAir = new G4LogicalBorderSurface("AirTarget",
                                physicPenSampleBox,
                                physicWorldBox,
                                AirTarget);

  //Reflector surface, air supposed to be between reflector and samples 
  //reflectivity is set to 98%
  G4OpticalSurface* scintWrap = new G4OpticalSurface("ScintWrap");

  scintWrap->SetModel(unified);
  scintWrap->SetType(dielectric_metal);
  scintWrap->SetFinish(polished);

  G4double energyReflector[] = {2.0*eV, 3.0*eV, 4.0*eV, 5.0*eV, 6.0*eV, 7.0*eV};
  const G4int num = sizeof(energyReflector)/sizeof(G4double);
  G4double reflectivity[] = {0.98, 0.98, 0.98, 0.98, 0.98, 0.98};
  assert(sizeof(reflectivity) == sizeof(energyReflector));

  G4MaterialPropertiesTable* scintWrapProperty = materialVikuiti->GetMaterialPropertiesTable();

  scintWrapProperty->AddProperty("REFLECTIVITY",energyReflector,reflectivity,num);
  scintWrap->SetMaterialPropertiesTable(scintWrapProperty); 

  //This should have no effect, since aire between scinitllator and samples
  logicSurfacePenReflectorFoilBoxOverPEN = new G4LogicalBorderSurface("ScintWrap", 
				physicPenSampleBox, 
				physicReflectorFoilBoxOverPEN, 
				scintWrap);
  
  logicSurfacePENStackedReflectorFoilBoxOverPEN = new G4LogicalBorderSurface("ScintWrap", 
				physicPenStackedSamples , 
				physicReflectorFoilBoxOverPEN, 
				scintWrap);


  //this surface is the one working, the other should work in case no air in between
  new G4LogicalBorderSurface("ScintWrap",
                                physicWorldBox,
                                physicReflectorFoilBoxOverPEN,
                                scintWrap);
 
  //Define surface in case a reflector is used with trigger foil
  logicSurfaceEJ212Reflector = new G4LogicalBorderSurface("ScintWrap", 
				physicTriggerFoilEJ212, 
				physicReflectorFoilAroundEJ212Foil, 
				scintWrap);

  logicSurfaceReflectorEJ212 = new G4LogicalBorderSurface("ScintWrap", 
				physicWorldBox, 
				physicReflectorFoilAroundEJ212Foil, 
				scintWrap);

  //Surface between air and thin scintillator
  G4OpticalSurface* AirEJ212 = new G4OpticalSurface("AirEJ212", 
				unified, 
				polished, 
				dielectric_dielectric);
  AirEJ212->SetPolish(1.0);

  new G4LogicalBorderSurface("AirEJ212",
                             physicTriggerFoilEJ212,
                             physicWorldBox,
                             AirEJ212);

  new G4LogicalBorderSurface("AirEJ212",
                             physicWorldBox,
                             physicTriggerFoilEJ212,
                             AirEJ212);

  //Surface between Air and PMMA light guide
  G4OpticalSurface* AirPMMA = new G4OpticalSurface("AirPMMA", 
				unified, 
				polished, 
				dielectric_dielectric);
  AirPMMA->SetPolish(1.0);

  logicSurfaceAirPMMA = new G4LogicalBorderSurface("AirPMMA", 
				physicLightGuide, 
				physicWorldBox, 
				AirPMMA);

  logicSurfacePMMAAir = new G4LogicalBorderSurface("PMMAAir", 
  				physicWorldBox, 
  				physicLightGuide, 
  				AirPMMA);

  //Surface between Light guide and trigger foil
  G4OpticalSurface* PMMATrigger = new G4OpticalSurface("PMMATrigger", 
				unified, 
				ground, 
				dielectric_dielectric);
  PMMATrigger->SetPolish(0.1);

  logicSurfaceEJ212PMMA = new G4LogicalBorderSurface("PMMATrigger", 
				physicLightGuide, 
				physicTriggerFoilEJ212, 
				PMMATrigger);
   
  logicSurfacePMMAEJ212 = new G4LogicalBorderSurface("PMMATrigger", 
				physicTriggerFoilEJ212, 
				physicLightGuide, 
				PMMATrigger);


  //surface between light guide and physicOpticalGrease
  G4OpticalSurface* PMMAGrease = new G4OpticalSurface("PMMAGrease", 
				unified, 
				ground, 
				dielectric_dielectric);
  PMMAGrease->SetPolish(0.1);
  
  new G4LogicalBorderSurface("PMMAGrease",
                            physicLightGuide,
                            physicOpticalGrease,
                            PMMAGrease);

  new G4LogicalBorderSurface("PMMAGrease",
                            physicOpticalGrease,
                            physicLightGuide,
                            PMMAGrease);

  G4OpticalSurface* PMMAPhotoCathodeSupport = new G4OpticalSurface("PMMACathodeSupport", 
				unified, 
				ground, 
				dielectric_dielectric);
  PMMAGrease->SetPolish(0.1);
  
  new G4LogicalBorderSurface("PMMACathodeSupport",
                           physicOpticalGrease,
                           physicBoxPhotoCathodeSupport,
                           PMMAPhotoCathodeSupport);
                           
  new G4LogicalBorderSurface("PMMACathodeSupport",
                           physicBoxPhotoCathodeSupport,
                           physicOpticalGrease,
                           PMMAPhotoCathodeSupport);
  
  new G4LogicalBorderSurface("PMMACathodeSupport",
                           physicBoxPhotoCathodeSupport,
                           physicLightGuide,
                           PMMAPhotoCathodeSupport);

  new G4LogicalBorderSurface("PMMACathodeSupport",
                           physicLightGuide,
                           physicBoxPhotoCathodeSupport,
                           PMMAPhotoCathodeSupport);
  

  //Photocathode Support in contact with the samples  
  surfaceCathodeSupport = new G4OpticalSurface("surfaceCathodeSupport", unified, ground, dielectric_dielectric);
  const G4int NUM = 2;
  G4double pp[NUM] = {2.038*eV, 4.144*eV};
  G4double reflectivityGlass[NUM] = {pmtReflectivity, pmtReflectivity};
  G4double efficiency[NUM] = {0.0, 0.0};
  MPT_GlassPMT->AddProperty("TRANSMITTANCE",pp,reflectivityGlass,NUM);
  MPT_GlassPMT->AddProperty("EFFICIENCY",pp,efficiency,NUM);
  //G4MaterialPropertiesTable* cathodeSupportProperty = new G4MaterialPropertiesTable();
  //G4MaterialPropertiesTable* cathodeSupportProperty = fGlass -> GetMaterialPropertiesTable();
  surfaceCathodeSupport -> SetSigmaAlpha(fSigAlphaSides);
  surfaceCathodeSupport -> SetMaterialPropertiesTable(MPT_GlassPMT); 

  new G4LogicalBorderSurface("surfaceCathodeSupport", 
				physicPenStackedSamples,
                                physicBoxPhotoCathodeSupport,
                                surfaceCathodeSupport);

  new G4LogicalBorderSurface("surfaceCathodeSupport", 
                                physicBoxPhotoCathodeSupport,
				physicPenStackedSamples,
                                surfaceCathodeSupport);

  new G4LogicalBorderSurface("surfaceCathodeSupport", 
                                physicBoxPhotoCathodeSupport,
				physicWorldBox,
                                surfaceCathodeSupport);

  new G4LogicalBorderSurface("surfaceCathodeSupport", 
				physicWorldBox,
                                physicBoxPhotoCathodeSupport,
                                surfaceCathodeSupport);

  new G4LogicalBorderSurface("surfaceCathodeSupport", 
				physicPenStackedSamples,
                                physicPhotoCathodeSupportPMT1,
                                surfaceCathodeSupport);
  new G4LogicalBorderSurface("surfaceCathodeSupport", 
                                physicPhotoCathodeSupportPMT1,
				physicPenStackedSamples,
                                surfaceCathodeSupport);

  new G4LogicalBorderSurface("surfaceCathodeSupport", 
				physicPenStackedSamples,
                                physicPhotoCathodeSupportPMT2,
                                surfaceCathodeSupport);
  new G4LogicalBorderSurface("surfaceCathodeSupport", 
                                physicPhotoCathodeSupportPMT2,
				physicPenStackedSamples,
                                surfaceCathodeSupport);

  new G4LogicalBorderSurface("surfaceCathodeSupport", 
				physicPenStackedSamples,
                                physicPhotoCathodeSupportPMT3,
                                surfaceCathodeSupport);
  new G4LogicalBorderSurface("surfaceCathodeSupport", 
                                physicPhotoCathodeSupportPMT3,
				physicPenStackedSamples,
                                surfaceCathodeSupport);

  new G4LogicalBorderSurface("surfaceCathodeSupport", 
				physicPenStackedSamples,
                                physicPhotoCathodeSupportPMT4,
                                surfaceCathodeSupport);
  new G4LogicalBorderSurface("surfaceCathodeSupport", 
                                physicPhotoCathodeSupportPMT4,
				physicPenStackedSamples,
                                surfaceCathodeSupport);

  new G4LogicalBorderSurface("surfaceCathodeSupport", 
				physicPenStackedSamples,
                                physicPhotoCathodeSupportPMT5,
                                surfaceCathodeSupport);
  new G4LogicalBorderSurface("surfaceCathodeSupport", 
                                physicPhotoCathodeSupportPMT5,
				physicPenStackedSamples,
                                surfaceCathodeSupport);

  // PMT photocathode
  // Quantum efficiency
  char filler;
  G4double wavelength;
  G4double cathodeEfficiency;
  G4double photocathEnergy[37];
  G4double photoCathodeQuantumEfficiency[37];
  //G4double perfectEfficiency[37];
  G4String pmtFile = "../input_files/QE_H11934_300.csv";

  ifstream ReadEff;
  G4int effCounter = 0;
  ReadEff.open(pmtFile);
  if(ReadEff.is_open())
  {
    while(!ReadEff.eof())
    {
      ReadEff>>wavelength>>filler>>cathodeEfficiency;
      photocathEnergy[effCounter] = (1240./wavelength)*eV;
      photoCathodeQuantumEfficiency[effCounter] = cathodeEfficiency/100.;
      //perfectEfficiency[effCounter] = 1.0;
      effCounter++;
      if(effCounter > 36){break;}
    }
  }
  else G4cout<<"Error opening file: " <<pmtFile<<G4endl;
  ReadEff.close();

  //Read real part of photocatode RINDEX
  G4double RealRefl,ImRefl, rindexEnergy;
  G4int nEntriesReR = 348;
  G4double photocath_ReR[348];
  G4double photocath_ReR_Energy[348];
  pmtFile = "../properties/Rindex_bialkali_real_energy.txt";
  
  effCounter = 0;
  ReadEff.open(pmtFile);
  if(ReadEff.is_open())
  {
    while(!ReadEff.eof())
    {
      ReadEff>>rindexEnergy>>RealRefl;
      photocath_ReR[effCounter] = RealRefl;
      photocath_ReR_Energy[effCounter] = rindexEnergy;
      G4cout<<" "<<effCounter<<" PMT reflectivity: energy ReR"<<photocath_ReR_Energy[effCounter]<<" "<<photocath_ReR[effCounter]<<G4endl;
      effCounter++;
      if(effCounter > (nEntriesReR-1)){break;}
    }
  }
  else G4cout<<"Error opening file: " <<pmtFile<<G4endl;
  ReadEff.close();


  G4int nEntriesImR = 496;
  G4double photocath_ImR[496];
  G4double photocath_ImR_Energy[496];
  pmtFile = "../properties/Rindex_bialkali_im_energy.txt";

  effCounter = 0;
  ReadEff.open(pmtFile);
  if(ReadEff.is_open())
  {
    while(!ReadEff.eof())
    {
      ReadEff>>rindexEnergy>>ImRefl;
      photocath_ImR[effCounter] = ImRefl;
      photocath_ImR_Energy[effCounter] = rindexEnergy;
      G4cout<<" "<<effCounter<<" PMT reflectivity im : Energy ImR "<<photocath_ImR_Energy[effCounter]<<" im "<<photocath_ImR[effCounter]<<G4endl;
      effCounter++;
      if(effCounter > (nEntriesImR-1)){break;}
    }
  }
  else G4cout<<"Error opening file: " <<pmtFile<<G4endl;
  ReadEff.close();

  //G4OpticalSurface*  phcath_opsurf = new G4OpticalSurface("phcath_opsurf",
    // unified, polished, dielectric_dielectric);
  //new G4LogicalBorderSurface("phcath_surf",physicActivePhotoCathodePMT1,physicPhotoCathodeSupportPMT1,  phcath_opsurf);
  //new G4LogicalBorderSurface("phcath_surf",physicActivePhotoCathodePMT2,physicPhotoCathodeSupportPMT2,  phcath_opsurf);
  //new G4LogicalBorderSurface("phcath_surf",physicActivePhotoCathodePMT3,physicPhotoCathodeSupportPMT3,  phcath_opsurf);
  //new G4LogicalBorderSurface("phcath_surf",physicActivePhotoCathodePMT4,physicPhotoCathodeSupportPMT4,  phcath_opsurf);
  //new G4LogicalBorderSurface("phcath_surf",physicActivePhotoCathodePMT5,physicPhotoCathodeSupportPMT5,  phcath_opsurf);
  //new G4LogicalBorderSurface("phcath_surf",physicPhotoCathodeSupportPMT2, physicActivePhotoCathodePMT2, phcath_opsurf);
  //new G4LogicalBorderSurface("phcath_surf",physicPhotoCathodeSupportPMT3, physicActivePhotoCathodePMT3, phcath_opsurf);
  //new G4LogicalBorderSurface("phcath_surf",physicPhotoCathodeSupportPMT4, physicActivePhotoCathodePMT4, phcath_opsurf);
  //new G4LogicalBorderSurface("phcath_surf",physicPhotoCathodeSupportPMT5, physicActivePhotoCathodePMT5, phcath_opsurf);


  const G4int nPMT_EFF = sizeof(photocathEnergy)/sizeof(G4double);
  G4OpticalSurface* pmtOpticalSurface = new G4OpticalSurface("pmt", glisur, polished, dielectric_metal);
  G4MaterialPropertiesTable* MPT_PMT = new G4MaterialPropertiesTable();
  MPT_PMT->AddProperty("EFFICIENCY", photocathEnergy, photoCathodeQuantumEfficiency,nPMT_EFF);
  MPT_PMT->AddProperty("REALRINDEX", photocath_ReR_Energy, photocath_ReR,nEntriesReR);
  MPT_PMT->AddProperty("IMAGINARYRINDEX", photocath_ImR_Energy, photocath_ImR, nEntriesImR);
  pmtOpticalSurface->SetMaterialPropertiesTable(MPT_PMT);
  new G4LogicalSkinSurface("pmt_surf", logicBoxActivePhotoCathodePMT,pmtOpticalSurface);

      
   return physicWorldBox;
}
