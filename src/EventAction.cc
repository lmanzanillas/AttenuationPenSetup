#include "EventAction.hh"
#include "RunAction.hh"
#include "Analysis.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"
#include "DetectorConstruction.hh"

#include "G4Event.hh"
#include "G4RunManager.hh"

#include "G4LogicalVolumeStore.hh"
#include "G4LogicalVolume.hh"
#include "G4Box.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

EventAction::EventAction(RunAction* runAction)
: G4UserEventAction(),
fRunAction(runAction)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

EventAction::~EventAction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void EventAction::BeginOfEventAction(const G4Event* myEvent)
{
	fDetectedPhotons = 0;
	fDepositedEnergy = 0;
	fDepositedEnergyTrigger = 0;

	fFrontPhoton =0;
	fBackPhoton = 0;
	fBottomPhoton=0;
	fLeftPhoton=0;
	fRightPhoton=0;
	if (myEvent->GetEventID() % 1000 == 0)
		G4cout << "event no.: " << myEvent->GetEventID() << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void EventAction::AddWavelength(G4double wavelength){
	auto analysisManager = G4AnalysisManager::Instance();

	analysisManager->FillH1(6, wavelength);
}

void EventAction::EndOfEventAction(const G4Event* myEvent)
{

	auto analysisManager = G4AnalysisManager::Instance();

	 if (fDetectedPhotons >= 1){
		//G4cout <<myEvent->GetEventID()<<": " << fDetectedPhotons << G4endl;
		analysisManager->FillNtupleDColumn(0,fDetectedPhotons);
		analysisManager->FillNtupleDColumn(1,fDepositedEnergyTrigger);

		if(fLeftPhoton > 5 || fRightPhoton > 5 || fBottomPhoton > 5 || fFrontPhoton > 5 || fBackPhoton > 5){
			analysisManager->FillNtupleDColumn(2, fDepositedEnergy);
			analysisManager->FillNtupleDColumn(3, fLeftPhoton);
			analysisManager->FillNtupleDColumn(4, fRightPhoton);
			analysisManager->FillNtupleDColumn(5, fBottomPhoton);
			analysisManager->FillNtupleDColumn(6, fFrontPhoton);
			analysisManager->FillNtupleDColumn(7, fBackPhoton);
			analysisManager->AddNtupleRow(0);
		}

	 }

}
