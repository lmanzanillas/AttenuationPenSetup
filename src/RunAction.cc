// Make this appear first!
#include "G4Timer.hh"

#include "RunAction.hh"
#include "RunActionMessenger.hh"

#include "G4Run.hh"

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

    strftime(buffer,80,"%d-%m-%Y %H-%M-%S",timeinfo);
    return std::string(buffer);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::BeginOfRunAction(const G4Run* aRun)
{
  fMan = G4AnalysisManager::Instance();

  G4String targetName = fDetector->GetTargetMaterialName();
  G4String targetThickness = G4BestUnit(fDetector->GetTargetSize(),"Length");
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

  G4String directorName = "../output/"+sourceString+"-" + datetime()+"/";
  mkdir(directorName, 0777);

  G4String position = G4BestUnit((fDetector->GetTargetSize()-fPrimary->GetSourcePosition()),"Length");

  G4cout << "### Run " << aRun->GetRunID() << " start." << G4endl;
  fTimer->Start();
  G4String abs_string = std::to_string(fDetector->GetABS());
  replace(abs_string.begin(),abs_string.end(),'.',',');

  G4String alpha_string = std::to_string(fDetector->GetSigAlpha());
  replace(alpha_string.begin(),alpha_string.end(),'.',',');

  fFileName = directorName+"_"+fDetector->GetDetectorName()+".csv";

  fMan->SetVerboseLevel(0);
  fMan->OpenFile(fFileName);

  fMan->CreateNtuple("DetectedPhotons","Escapes");
  fMan->CreateNtupleDColumn("EDep_Trigger");
  fMan->CreateNtupleDColumn("N_Trigger");
  fMan->CreateNtupleDColumn("EDep_PEN");
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
