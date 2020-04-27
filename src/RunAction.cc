// Make this appear first!
#include "G4Timer.hh"

#include "RunAction.hh"
#include "RunActionMessenger.hh"

#include "G4Run.hh"
#include "G4Types.hh"

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
  G4String positionX = std::to_string((int)fPrimary->GetSourcePositionX());


  G4cout << "### Run " << aRun->GetRunID() << " start." << G4endl;
  fTimer->Start();
  G4String abs_string = std::to_string((int)fDetector->GetABS());
  //replace(abs_string.begin(),abs_string.end(),'.',',');

  //G4String alpha_string = std::to_string(fDetector->GetSigAlpha());
  //replace(alpha_string.begin(),alpha_string.end(),'.',',');

  G4String directorName = "../output/"+sourceString+"_" + datetime()+"_x_"+positionX+"_mm_LY_" + std::to_string((int)fDetector->GetLY())+"ph_MeV_ABS_"+abs_string+"_SigAlpha_"+std::to_string(fDetector->GetSigAlpha())+"/";
  mkdir(directorName, 0777);
  fFileName = directorName+fDetector->GetDetectorName()+".csv";

  fMan->SetVerboseLevel(0);
  fMan->OpenFile(fFileName);

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
