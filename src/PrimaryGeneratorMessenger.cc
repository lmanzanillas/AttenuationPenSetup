#include "PrimaryGeneratorMessenger.hh"

#include "PrimaryGeneratorAction.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4SystemOfUnits.hh"
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PrimaryGeneratorMessenger::PrimaryGeneratorMessenger(PrimaryGeneratorAction* Gun)
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
  fSourceEnergy->SetDefaultValue(2000.*keV);
  fSourceType->AvailableForStates(G4State_Idle);

  fSourcePositionX = new G4UIcmdWithADoubleAndUnit("/PEN/gun/sourcePositionX",this);
  fSourcePositionX->SetGuidance("Set Source x position");
  fSourcePositionX->SetParameterName("fPositionX",true);
  fSourcePositionX->SetDefaultValue(0.*mm);
  fSourcePositionX->AvailableForStates(G4State_Idle);

  fSourcePositionY = new G4UIcmdWithADoubleAndUnit("/PEN/gun/sourcePositionY",this);
  fSourcePositionY->SetGuidance("Set Source y position");
  fSourcePositionY->SetParameterName("fPositionY",true);
  fSourcePositionY->SetDefaultValue(0.*mm);
  fSourcePositionY->AvailableForStates(G4State_Idle);

  fSourcePositionZ = new G4UIcmdWithADoubleAndUnit("/PEN/gun/sourcePositionZ",this);
  fSourcePositionZ->SetGuidance("Set Source z position");
  fSourcePositionZ->SetParameterName("fPositionZ",true);
  fSourcePositionZ->SetDefaultValue(30.*mm);
  fSourcePositionZ->AvailableForStates(G4State_Idle);

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PrimaryGeneratorMessenger::~PrimaryGeneratorMessenger()
{
  delete fSourceType;
  delete fSourceEnergy;
  delete fGunDir;
  delete fSourcePositionX;
  delete fSourcePositionY;
  delete fSourcePositionZ;
  //delete fAction;//debug
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PrimaryGeneratorMessenger::SetNewValue(G4UIcommand* command, G4String newValue)
{
  if(command == fSourceType) {
  	fAction->SetSourceType(fSourceType->GetNewIntValue(newValue));
  }

  if(command == fSourceEnergy){
  	fAction->SetSourceType(3);
  	fAction->SetSourceEnergy(fSourceEnergy->GetNewDoubleValue(newValue));
  }
  if(command == fSourcePositionX){
  	fAction->SetSourcePositionX(fSourcePositionX->GetNewDoubleValue(newValue));
  }
  if(command == fSourcePositionY){
  	fAction->SetSourcePositionY(fSourcePositionY->GetNewDoubleValue(newValue));
  }
  if(command == fSourcePositionZ){
  	fAction->SetSourcePositionZ(fSourcePositionZ->GetNewDoubleValue(newValue));
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
