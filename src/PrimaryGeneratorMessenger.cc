#include "PrimaryGeneratorMessenger.hh"

#include "PrimaryGeneratorAction.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4SystemOfUnits.hh"
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PrimaryGeneratorMessenger::
  PrimaryGeneratorMessenger(PrimaryGeneratorAction* Gun)
  : G4UImessenger(),
    fAction(Gun),fGunDir(0)
{

  fGunDir = new G4UIdirectory("/PEN/gun/");
  fGunDir->SetGuidance("PrimaryGenerator control");

  fSourceType = new G4UIcmdWithAnInteger("/PEN/gun/sourceType",this);
  fSourceType->SetGuidance("Choose the type of source");
  fSourceType->SetParameterName("sourceType",true);
  fSourceType->SetDefaultValue(0);
  fSourceType->AvailableForStates(G4State_Idle);

  fSourceEnergy = new G4UIcmdWithADoubleAndUnit("/PEN/gun/sourceEnergy", this);
  fSourceEnergy->SetGuidance("Choose source energy");
  fSourceEnergy->SetParameterName("sourceEnergy",true);
  fSourceEnergy->SetDefaultValue(60*keV);
  fSourceType->AvailableForStates(G4State_Idle);

  fSourcePosition = new G4UIcmdWithADoubleAndUnit("/PEN/gun/sourcePosition",this);
  fSourcePosition->SetGuidance("Set Source z position");
  fSourcePosition->SetParameterName("fPosition",true);
  fSourcePosition->SetDefaultValue(2*cm);
  fSourcePosition->AvailableForStates(G4State_Idle);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PrimaryGeneratorMessenger::~PrimaryGeneratorMessenger()
{
  delete fSourceType;
  delete fGunDir;
  delete fSourcePosition;
  delete fSourceEnergy;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PrimaryGeneratorMessenger::SetNewValue(
                                        G4UIcommand* command, G4String newValue)
{
  if( command == fSourceType ) {
      fAction->SetSourceType(fSourceType->GetNewIntValue(newValue));
  }

  if (command == fSourceEnergy){
    fAction->SetSourceType(0);
    fAction->SetSourceEnergy(fSourceEnergy->GetNewDoubleValue(newValue));
  }
  if (command == fSourcePosition){
    fAction->SetSourcePosition(fSourcePosition->GetNewDoubleValue(newValue));
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
