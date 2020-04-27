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
	fDetectedPhotons = 0.;
	depositedEnergyPENStackedSample1 = 0.;
	depositedEnergyPENStackedSample2 = 0.;
	depositedEnergyPENStackedSample3 = 0.;
	depositedEnergyPENStackedSample4 = 0.;
	depositedEnergyTriggerFoilEJ212 = 0.;
	depositedEnergyInactiveMaterial = 0.;

	fFrontPhoton = 0.;
	fBackPhoton = 0.;
	fBottomPhoton = 0.;
	fLeftPhoton = 0.;
	fRightPhoton = 0.;
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

	//if (fDetectedPhotons >= 0){
	if (depositedEnergyTriggerFoilEJ212 > 0.){
		//G4cout <<myEvent->GetEventID()<<": " << fDetectedPhotons << G4endl;
		analysisManager->FillNtupleDColumn(0, myEvent->GetEventID());
		analysisManager->FillNtupleDColumn(1, depositedEnergyTriggerFoilEJ212);
		analysisManager->FillNtupleDColumn(2, depositedEnergyPENStackedSample1);
		analysisManager->FillNtupleDColumn(3, depositedEnergyPENStackedSample2);
		analysisManager->FillNtupleDColumn(4, depositedEnergyPENStackedSample3);
		analysisManager->FillNtupleDColumn(5, depositedEnergyPENStackedSample4);
		analysisManager->FillNtupleDColumn(6, depositedEnergyInactiveMaterial);
		analysisManager->FillNtupleDColumn(7, fDetectedPhotons);
		analysisManager->FillNtupleDColumn(8, fLeftPhoton);
		analysisManager->FillNtupleDColumn(9, fRightPhoton);
		analysisManager->FillNtupleDColumn(10, fBottomPhoton);
		analysisManager->FillNtupleDColumn(11, fFrontPhoton);
		analysisManager->FillNtupleDColumn(12, fBackPhoton);
		analysisManager->AddNtupleRow(0);

	 }

}
