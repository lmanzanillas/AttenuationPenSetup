#include "DetectorMessenger.hh"

#include "DetectorConstruction.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithADouble.hh"
#include "G4UIcmdWithoutParameter.hh"
#include "G4UIcmdWithAnInteger.hh"

DetectorMessenger::DetectorMessenger(DetectorConstruction * Det)
:G4UImessenger(),fDetector(Det),
 fPENDir(0),
 fDetDir(0),
 fTMaterCMD(0),
 fWMaterCMD(0),
 fSizeCMD(0),
 fTypeCMD(0),
 fLYCMD(0),fResCMD(0),fDetNameCMD(0),fABSCMD(0), fAbsFileCMD(0)
{
  fDetDir = new G4UIdirectory("/PEN/det/");
  fDetDir->SetGuidance("detector construction commands");

  fTMaterCMD = new G4UIcmdWithAString("/PEN/det/setTargetMat",this);
  fTMaterCMD->SetGuidance("Select material of the target.");
  fTMaterCMD->SetParameterName("choice",false);
  fTMaterCMD->AvailableForStates(G4State_PreInit,G4State_Idle);
  fTMaterCMD->SetToBeBroadcasted(false);

  fWMaterCMD = new G4UIcmdWithAString("/PEN/det/setWorldMat",this);
  fWMaterCMD->SetGuidance("Select material of the world.");
  fWMaterCMD->SetParameterName("choice",false);
  fWMaterCMD->AvailableForStates(G4State_PreInit,G4State_Idle);
  fWMaterCMD->SetToBeBroadcasted(false);

  fSizeCMD = new G4UIcmdWithADoubleAndUnit("/PEN/det/setSize",this);
  fSizeCMD->SetGuidance("Set size of the box");
  fSizeCMD->SetParameterName("Size",false);
  fSizeCMD->SetRange("Size>0.");
  fSizeCMD->SetUnitCategory("Length");
  fSizeCMD->AvailableForStates(G4State_PreInit,G4State_Idle);
  fSizeCMD->SetToBeBroadcasted(false);

  fTypeCMD = new G4UIcmdWithAnInteger("/PEN/det/setDetectorType",this);
  fTypeCMD->SetGuidance("Set detector type");
  fTypeCMD->AvailableForStates(G4State_PreInit,G4State_Idle);
  fTypeCMD->SetToBeBroadcasted(false);

  fLYCMD = new G4UIcmdWithADouble("/PEN/det/setLY",this);
  fLYCMD->SetGuidance("Set scint LY");
  fLYCMD->SetParameterName("Light Yield",false);
  fLYCMD->AvailableForStates(G4State_PreInit,G4State_Idle);
  fLYCMD->SetToBeBroadcasted(false);

  fResCMD = new G4UIcmdWithADouble("/PEN/det/setRes",this);
  fResCMD->SetGuidance("Set Res");
  fResCMD->AvailableForStates(G4State_PreInit,G4State_Idle);
  fResCMD->SetToBeBroadcasted(false);

  fABSCMD = new G4UIcmdWithADouble("/PEN/det/setABS",this);
  fABSCMD->SetGuidance("Set ABS");
  fABSCMD->AvailableForStates(G4State_Init,G4State_Idle);
  fABSCMD->SetToBeBroadcasted(false);

  fAlphaCMD = new G4UIcmdWithADouble("/PEN/det/setSigAlpha",this);
  fAlphaCMD->SetGuidance("Set fSigAlpha");
  fAlphaCMD->AvailableForStates(G4State_Init,G4State_Idle);
  fAlphaCMD->SetToBeBroadcasted(false);

  fDetNameCMD = new G4UIcmdWithAString("/PEN/det/setDetName",this);
  fDetNameCMD->SetGuidance("Set detetector file name.");
  fDetNameCMD->SetParameterName("choice",false);
  fDetNameCMD->AvailableForStates(G4State_PreInit,G4State_Idle);
  fDetNameCMD->SetToBeBroadcasted(false);

  fAbsFileCMD = new G4UIcmdWithAString("/PEN/det/setAbsFile",this);
  fAbsFileCMD->SetGuidance("Set abs file name.");
  fAbsFileCMD->SetParameterName("choice",false);
  fAbsFileCMD->AvailableForStates(G4State_PreInit,G4State_Idle);
  fAbsFileCMD->SetToBeBroadcasted(false);

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorMessenger::~DetectorMessenger()
{
  delete fTMaterCMD;
  delete fWMaterCMD;
  delete fSizeCMD;
  delete fTypeCMD;
  delete fDetDir;
  delete fPENDir;
  delete fLYCMD;
  delete fResCMD;
  delete fABSCMD;
  delete fDetNameCMD;
  delete fAbsFileCMD;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorMessenger::SetNewValue(G4UIcommand* command,G4String newValue)
{
   if( command == fTMaterCMD )
    { fDetector->SetTargetMaterial(newValue);}

   if( command == fWMaterCMD )
    { fDetector->SetWorldMaterial(newValue);}

    if( command == fDetNameCMD )
     { fDetector->SetDetectorName(newValue);}

   if( command == fSizeCMD )
    { fDetector->SetSize(fSizeCMD->GetNewDoubleValue(newValue));}

  if( command == fTypeCMD )
    { fDetector->SetDetectorType(fTypeCMD->GetNewIntValue(newValue));}

  if(command == fLYCMD)
  {fDetector->SetLY(fLYCMD->GetNewDoubleValue(newValue));}

  if(command==fResCMD)
  {fDetector->SetRes(fResCMD->GetNewDoubleValue(newValue));}

  if(command==fABSCMD)
  {fDetector->SetABS(fABSCMD->GetNewDoubleValue(newValue));}

  if(command==fAlphaCMD)
  {fDetector->SetSigAlpha(fAlphaCMD->GetNewDoubleValue(newValue));}

  if(command==fAbsFileCMD)
  {fDetector->SetABSFile(newValue);}
}
