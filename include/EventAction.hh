#ifndef EventAction_h
#define EventAction_h 1

#include "G4UserEventAction.hh"
#include "globals.hh"
#include "G4ThreeVector.hh"

class RunAction;

/// Event action class
///

class EventAction : public G4UserEventAction
{
  public:
    EventAction(RunAction* runAction);
    virtual ~EventAction();

    virtual void BeginOfEventAction(const G4Event* event);
    virtual void EndOfEventAction(const G4Event* event);
    void    AddDetectedPhoton(void){fDetectedPhotons++;};
    G4double   GetNumberDetectedPhotons(void){return fDetectedPhotons;};
    void    AddDepositedEnergy(G4double newEnergy){fDepositedEnergy += newEnergy;};
    void    AddDepositedEnergyTrigger(G4double newEnergy){fDepositedEnergyTrigger += newEnergy;};
    void    AddWavelength(G4double newWavelength);
    void    AddIWavelength(G4double startWavelength);

    void AddFrontPhoton(void){fFrontPhoton++;};
    void AddBackPhoton(void){fBackPhoton++;};
    void AddBottomPhoton(void){fBottomPhoton++;};
    void AddLeftPhoton(void){fLeftPhoton++;};
    void AddRightPhoton(void){fRightPhoton++;};

  private:
    RunAction* 	fRunAction;

    G4double fDepositedEnergy;
    G4double fDepositedEnergyTrigger;

    G4double fDetectedPhotons;
    G4double fFrontPhoton;
    G4double fBackPhoton;
    G4double fBottomPhoton;
    G4double fLeftPhoton;
    G4double fRightPhoton;

};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
