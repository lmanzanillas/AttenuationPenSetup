// Make this appear first!
#include "G4Timer.hh"

#include "RunAction.hh"
#include "RunActionMessenger.hh"

#include "G4Run.hh"
#include "G4Types.hh"
#include "G4Material.hh"

#include "Analysis.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"

#include "EventAction.hh"

#include <string>
#include <ctime>
#include <sys/stat.h>

#include "DetectorConstruction.hh"
#include "PrimaryGeneratorAction.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RunAction::RunAction(DetectorConstruction* det, PrimaryGeneratorAction* kin)
 : G4UserRunAction(),fDetector(det),fPrimary(kin),
   fTimer(0)
{
  fTimer = new G4Timer;
  fMan = G4AnalysisManager::Instance();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RunAction::~RunAction()
{
  delete fTimer;
}

std::string datetime()
{
    time_t rawtime;
    struct tm * timeinfo;
    char buffer[80];

    time (&rawtime);
    timeinfo = localtime(&rawtime);

    strftime(buffer,80,"%d_%m_%Y-%H-%M-%S",timeinfo);
    return std::string(buffer);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::BeginOfRunAction(const G4Run* aRun)
{
  fMan = G4AnalysisManager::Instance();

  G4String targetName = fDetector->GetTargetMaterialName();
  G4String targetThickness = G4BestUnit(fDetector->GetTargetSampleThickness(),"Length");
  G4int sourceName = fPrimary->GetSourceType();
  G4String sourceString;

  switch (sourceName){
    case 0:
      sourceString = "Cs137";
      break;
    case 1:
      sourceString = "Bi207";
      break;
    case 2:
      sourceString = "Sr90";
      break;
    case 3:
      sourceString = "Mono-beta";
      break;
  }

  G4cout << sourceString << G4endl;

  //G4String positionZ = G4BestUnit((fDetector->GetTargetSampleThickness()-fPrimary->GetSourcePositionZ()),"Length");
  G4String positionX ="_x_"+std::to_string((int)fPrimary->GetSourcePositionX())+"_mm";


  G4cout << "### Run " << aRun->GetRunID() << " start." << G4endl;
  fTimer->Start();
  std::stringstream stream_abs;
  stream_abs << std::fixed << std::setprecision(2) <<fDetector->GetABS();
  G4String abs_string = stream_abs.str();
  G4String s_LY = "_LY_"+std::to_string((int)fDetector->GetTargetMaterial()->GetMaterialPropertiesTable()->GetConstProperty("SCINTILLATIONYIELD"))+"ph_MeV_ABS_";
  G4String s_detector_type = "_Det_"+std::to_string(fDetector->GetDetectorType());
  //SigmaAlpha sides
  std::stringstream stream_SigmaAlphaSides;
  stream_SigmaAlphaSides << std::fixed << std::setprecision(2) <<fDetector->GetSigAlphaSides();
  G4String s_SigAlphaSides ="_SigAlphaSides_"+ stream_SigmaAlphaSides.str();
  //SigmaAlpha top bottom
  std::stringstream stream_SigmaAlphaBottom;
  stream_SigmaAlphaBottom << std::fixed << std::setprecision(2) <<fDetector->GetSigAlphaBottom();
  G4String s_SigAlphaBottom ="_SigAlphaBottom_"+ stream_SigmaAlphaBottom.str();
  //transmittance side pmts
  std::stringstream stream_PMT_reflectivitySides;
  stream_PMT_reflectivitySides << std::fixed << std::setprecision(2) <<fDetector->GetPMTReflectivitySides();
  G4String s_PMT_reflectivitySides ="_PmtReflecSide_"+ stream_PMT_reflectivitySides.str();
  //transmittance bottom pmt
  std::stringstream stream_PMT_reflectivityBottom;
  stream_PMT_reflectivityBottom << std::fixed << std::setprecision(2) <<fDetector->GetPMTReflectivityBottom();
  G4String s_PMT_reflectivityBottom ="_PmtReflecBottom_"+ stream_PMT_reflectivityBottom.str();
  //Target material
  G4String s_Target_Material ="_"+ fDetector->GetTargetMaterialName();
  //n Samples
  G4String s_n_targets = "_"+std::to_string(fDetector->GetNumberOfTargetSamples())+"_samples";
  //Thickness
  std::stringstream stream_thickness;
  stream_thickness << std::fixed << std::setprecision(2) <<2*fDetector->GetTargetSampleThickness();
  G4String s_target_thickness ="_"+ stream_thickness.str()+"mm";


  G4String directorName = "../output/"+sourceString+"_"
	  +datetime()
	  +positionX 
	  +s_LY
	  +abs_string
	  +s_detector_type
          +s_n_targets
          +s_target_thickness
	  +s_Target_Material+"/";

  mkdir(directorName, 0777);
  
  fFileName = directorName+fDetector->GetDetectorName()
	  +s_SigAlphaSides
	  +s_SigAlphaBottom
	  +s_PMT_reflectivitySides
	  +s_PMT_reflectivityBottom
          +".csv";

  fMan->SetVerboseLevel(0);
  fMan->OpenFile(fFileName);
  if(fDetector->GetDetectorType()==3){
 	 fMan->CreateNtuple("PEN","Detailed MC Info");
	 fMan->CreateNtupleDColumn("PhotonWL");
	 fMan->CreateNtupleDColumn("PhotonDistance");
	 fMan->CreateNtupleDColumn("IsPhotonDetected");
	fMan->FinishNtuple();
  }else{
  	fMan->CreateNtuple("PEN","Detailed MC Info");
	  fMan->CreateNtupleDColumn("Event_Number");
	  fMan->CreateNtupleDColumn("EDep_Trigger_FoilEJ212");
	  fMan->CreateNtupleDColumn("EDep_PEN_Stacked1");
	  fMan->CreateNtupleDColumn("EDep_PEN_Stacked2");
	  fMan->CreateNtupleDColumn("EDep_PEN_Stacked3");
	  fMan->CreateNtupleDColumn("EDep_PEN_Stacked4");
	  fMan->CreateNtupleDColumn("EDep_InactiveMaterial");
	  fMan->CreateNtupleDColumn("N_Trigger");
	  fMan->CreateNtupleDColumn("N_Left");
	  fMan->CreateNtupleDColumn("N_Right");
	  fMan->CreateNtupleDColumn("N_Bottom");
	  fMan->CreateNtupleDColumn("N_Front");
	  fMan->CreateNtupleDColumn("N_Back");
	  fMan->FinishNtuple();
  }

}

void RunAction::SetFileName(G4String fileName)
{
  fFileName = fileName;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::EndOfRunAction(const G4Run* aRun)
{
  fTimer->Stop();
  G4cout << "number of event = " << aRun->GetNumberOfEvent()
         << " " << *fTimer << G4endl;
    G4cout << "End Run" << G4endl;

    fMan->Write();
    fMan->CloseFile();
    delete fMan;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
