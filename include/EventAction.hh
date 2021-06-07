#ifndef EventAction_h
#define EventAction_h 1

#include "G4UserEventAction.hh"
#include "globals.hh"
#include "G4ThreeVector.hh"

class DetectorConstruction;
class RunAction;

/// Event action class
///

class EventAction : public G4UserEventAction
{
  public:
    EventAction(DetectorConstruction* det, RunAction* runAction);
    virtual ~EventAction();

    virtual void BeginOfEventAction(const G4Event* event);
    virtual void EndOfEventAction(const G4Event* event);
    void    AddDetectedPhoton(void){fDetectedPhotons++;};
    G4double   GetNumberDetectedPhotons(void){return fDetectedPhotons;};
    void    AddDepositedEnergyPENStackedSample1(G4double newEnergy){depositedEnergyPENStackedSample1 += newEnergy;};
    void    AddDepositedEnergyPENStackedSample2(G4double newEnergy){depositedEnergyPENStackedSample2 += newEnergy;};
    void    AddDepositedEnergyPENStackedSample3(G4double newEnergy){depositedEnergyPENStackedSample3 += newEnergy;};
    void    AddDepositedEnergyPENStackedSample4(G4double newEnergy){depositedEnergyPENStackedSample4 += newEnergy;};
    void    AddDepositedEnergyEJ212TriggerFoil(G4double newEnergy){depositedEnergyTriggerFoilEJ212 += newEnergy;};
    void    AddDepositedEnergyInactiveMaterial(G4double newEnergy){depositedEnergyInactiveMaterial += newEnergy;};

    void    AddFirstPenX(G4double new_x){xFirstPen = new_x;};
    void    AddFirstPenY(G4double new_y){yFirstPen = new_y;};
    void    AddFirstPenZ(G4double new_z){zFirstPen = new_z;};

    void    AddWaveLength(G4double PhotonWL);
    void    AddPhotonTravelledDistance(G4double PhotonTravelledDistance);
    void    AddIsPhotonDetected(G4int IsPhotonDetected);

    void AddFrontPhoton(void){fFrontPhoton++;};
    void AddBackPhoton(void){fBackPhoton++;};
    void AddBottomPhoton(void){fBottomPhoton++;};
    void AddLeftPhoton(void){fLeftPhoton++;};
    void AddRightPhoton(void){fRightPhoton++;};

  private:

    DetectorConstruction*   fDetector;
    RunAction* 	fRunAction;

    G4double depositedEnergyPENStackedSample1;
    G4double depositedEnergyPENStackedSample2;
    G4double depositedEnergyPENStackedSample3;
    G4double depositedEnergyPENStackedSample4;
    G4double depositedEnergyTriggerFoilEJ212;
    G4double depositedEnergyInactiveMaterial;

    G4double fDetectedPhotons;
    G4double fFrontPhoton;
    G4double fBackPhoton;
    G4double fBottomPhoton;
    G4double fLeftPhoton;
    G4double fRightPhoton;
     
    //position of first interaction
    G4double xFirstPen;
    G4double yFirstPen;
    G4double zFirstPen;
     
    //Spectrometer information
    G4double PhotonWL;
    G4double PhotonTravelledDistance;
    G4int IsPhotonDetected;

};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
