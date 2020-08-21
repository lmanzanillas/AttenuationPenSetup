#include "DetectorMessenger.hh"

#include "DetectorConstruction.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithADouble.hh"
#include "G4UIcmdWithoutParameter.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcmdWithABool.hh"

#include "G4RunManager.hh"

DetectorMessenger::DetectorMessenger(DetectorConstruction * Det)
:G4UImessenger(),
 fDetector(Det),
 fPENDir(0),
 fDetDir(0),
 commandSetWorldMaterial(0),
 commandSetDetectorName(0),
 commandSetDetectorType(0),
 commandSetNumberOfTargetSamples(0),
 commandSetTargetSampleLength(0),
 commandSetTargetSampleThickness(0),
 commandSetTargetSampleWidth(0),
 commandSetTargetMaterial(0),
 commandSetAlphaSigma(0),
 commandSetAlphaSigmaSides(0),
 commandSetAlphaSigmaBottom(0),
 commandSetPMTReflectivity(0),
 commandSetCollimatorPositionX(0),
 commandSetLY(0),
 commandSetResolutionLY(0),
 commandSetPenAbsorption(0), 
 commandSetAbsorptionFile(0)
 //commandSetReflectorOn(false),
 {
  fDetDir = new G4UIdirectory("/PEN/det/");
  fDetDir->SetGuidance("detector construction commands");

  commandSetTargetMaterial = new G4UIcmdWithAString("/PEN/det/setTargetMat",this);
  commandSetTargetMaterial->SetGuidance("Select material of the target.");
  commandSetTargetMaterial->SetParameterName("choice",false);
  commandSetTargetMaterial->AvailableForStates(G4State_PreInit,G4State_Idle);
  commandSetTargetMaterial->SetToBeBroadcasted(false);

  commandSetWorldMaterial = new G4UIcmdWithAString("/PEN/det/setWorldMat",this);
  commandSetWorldMaterial->SetGuidance("Select material of the world.");
  commandSetWorldMaterial->SetParameterName("choice",false);
  commandSetWorldMaterial->AvailableForStates(G4State_PreInit,G4State_Idle);
  commandSetWorldMaterial->SetToBeBroadcasted(false);

  commandSetTargetSampleLength = new G4UIcmdWithADoubleAndUnit("/PEN/det/setTargetLength",this);
  commandSetTargetSampleLength->SetGuidance("Set length of target samples");
  commandSetTargetSampleLength->SetParameterName("SampleLength",false);
  commandSetTargetSampleLength->SetRange("SampleLength>0.");
  commandSetTargetSampleLength->SetUnitCategory("Length");
  commandSetTargetSampleLength->AvailableForStates(G4State_PreInit,G4State_Idle);
  commandSetTargetSampleLength->SetToBeBroadcasted(false);

  commandSetTargetSampleThickness = new G4UIcmdWithADoubleAndUnit("/PEN/det/setTargetThickness",this);
  commandSetTargetSampleThickness->SetGuidance("Set thickness of target samples");
  commandSetTargetSampleThickness->SetParameterName("SampleThickness",false);
  commandSetTargetSampleThickness->SetRange("SampleThickness>0.");
  commandSetTargetSampleThickness->SetUnitCategory("Length");
  commandSetTargetSampleThickness->AvailableForStates(G4State_PreInit,G4State_Idle);
  commandSetTargetSampleThickness->SetToBeBroadcasted(false);

  commandSetTargetSampleWidth = new G4UIcmdWithADoubleAndUnit("/PEN/det/setTargetWidth",this);
  commandSetTargetSampleWidth->SetGuidance("Set width of target samples");
  commandSetTargetSampleWidth->SetParameterName("SampleWidth",false);
  commandSetTargetSampleWidth->SetRange("SampleWidth>0.");
  commandSetTargetSampleWidth->SetUnitCategory("Length");
  commandSetTargetSampleWidth->AvailableForStates(G4State_PreInit,G4State_Idle);
  commandSetTargetSampleWidth->SetToBeBroadcasted(false);

  commandSetDetectorType = new G4UIcmdWithAnInteger("/PEN/det/setDetectorType",this);
  commandSetDetectorType->SetGuidance("Set detector type");
  commandSetDetectorType->AvailableForStates(G4State_PreInit,G4State_Idle);
  commandSetDetectorType->SetToBeBroadcasted(false);

  commandSetNumberOfTargetSamples = new G4UIcmdWithAnInteger("/PEN/det/setNTargetSamples",this);
  commandSetNumberOfTargetSamples->SetGuidance("Set number of target Samples");
  commandSetNumberOfTargetSamples->AvailableForStates(G4State_PreInit,G4State_Idle);
  commandSetNumberOfTargetSamples->SetToBeBroadcasted(false);

  
  commandSetCollimatorPositionX = new G4UIcmdWithADoubleAndUnit("/PEN/det/setCollimatorPositionX",this);
  commandSetCollimatorPositionX ->SetGuidance("Set position of collimator in x");
  commandSetCollimatorPositionX ->SetParameterName("CollimatorPositionX",false);
  commandSetCollimatorPositionX ->SetRange("CollimatorPositionX>-37.4&&CollimatorPositionX<37.4");//limits according to the size of the samples
  commandSetCollimatorPositionX ->SetUnitCategory("Length");
  commandSetCollimatorPositionX ->AvailableForStates(G4State_PreInit,G4State_Idle);
  commandSetCollimatorPositionX ->SetToBeBroadcasted(false);

  commandSetLY = new G4UIcmdWithADouble("/PEN/det/setLY",this);
  commandSetLY->SetGuidance("Set scint LY");
  commandSetLY->SetParameterName("Light Yield",false);
  commandSetLY->AvailableForStates(G4State_PreInit,G4State_Idle);
  commandSetLY->SetToBeBroadcasted(false);

  commandSetResolutionLY = new G4UIcmdWithADouble("/PEN/det/setRes",this);
  commandSetResolutionLY->SetGuidance("Set Res");
  commandSetResolutionLY->AvailableForStates(G4State_PreInit,G4State_Idle);
  commandSetResolutionLY->SetToBeBroadcasted(false);

  commandSetPenAbsorption = new G4UIcmdWithADouble("/PEN/det/setABS",this);
  commandSetPenAbsorption->SetGuidance("Set ABS");
  commandSetPenAbsorption->AvailableForStates(G4State_Init,G4State_Idle);
  commandSetPenAbsorption->SetToBeBroadcasted(false);

  commandSetAlphaSigma = new G4UIcmdWithADouble("/PEN/det/setSigAlpha",this);
  commandSetAlphaSigma->SetGuidance("Set surf SigAlpha");
  commandSetAlphaSigma->AvailableForStates(G4State_Init,G4State_Idle);
  commandSetAlphaSigma->SetToBeBroadcasted(false);

  commandSetAlphaSigmaSides = new G4UIcmdWithADouble("/PEN/det/setSigAlphaSides",this);
  commandSetAlphaSigmaSides->SetGuidance("Set surf SigAlpha in contact with PMT");
  commandSetAlphaSigmaSides->AvailableForStates(G4State_Init,G4State_Idle);
  commandSetAlphaSigmaSides->SetToBeBroadcasted(false);

  commandSetAlphaSigmaBottom = new G4UIcmdWithADouble("/PEN/det/setSigAlphaBottom",this);
  commandSetAlphaSigmaBottom->SetGuidance("Set surf SigAlpha in contact with bottom PMT");
  commandSetAlphaSigmaBottom->AvailableForStates(G4State_Init,G4State_Idle);
  commandSetAlphaSigmaBottom->SetToBeBroadcasted(false);

  commandSetPMTReflectivity = new G4UIcmdWithADouble("/PEN/det/setPMTReflectivity",this);
  commandSetPMTReflectivity->SetGuidance("Set PMT reflectivity");
  commandSetPMTReflectivity->AvailableForStates(G4State_Init,G4State_Idle);
  commandSetPMTReflectivity->SetToBeBroadcasted(false);

  commandSetDetectorName = new G4UIcmdWithAString("/PEN/det/setDetName",this);
  commandSetDetectorName->SetGuidance("Set detetector file name.");
  commandSetDetectorName->SetParameterName("choice",false);
  commandSetDetectorName->AvailableForStates(G4State_PreInit,G4State_Idle);
  commandSetDetectorName->SetToBeBroadcasted(false);

  commandSetAbsorptionFile = new G4UIcmdWithAString("/PEN/det/setAbsFile",this);
  commandSetAbsorptionFile->SetGuidance("Set abs file name.");
  commandSetAbsorptionFile->SetParameterName("choice",false);
  commandSetAbsorptionFile->AvailableForStates(G4State_PreInit,G4State_Idle);
  commandSetAbsorptionFile->SetToBeBroadcasted(false);

  commandSetReflectorOn = new G4UIcmdWithABool("/PEN/det/setReflectorOn",this);
  commandSetReflectorOn->SetGuidance("Enable/Disable reflector.");
  commandSetReflectorOn->AvailableForStates(G4State_PreInit,G4State_Idle);
  commandSetReflectorOn->SetToBeBroadcasted(false);

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorMessenger::~DetectorMessenger()
{
  delete fDetDir;
  delete fPENDir;
  delete commandSetWorldMaterial;
  delete commandSetDetectorName;
  delete commandSetDetectorType;
  delete commandSetNumberOfTargetSamples;
  delete commandSetTargetSampleLength;
  delete commandSetTargetSampleThickness;
  delete commandSetTargetSampleWidth;
  delete commandSetTargetMaterial;
  delete commandSetAlphaSigma;
  delete commandSetAlphaSigmaSides;
  delete commandSetAlphaSigmaBottom;
  delete commandSetPMTReflectivity;
  delete commandSetCollimatorPositionX;
  delete commandSetLY;
  delete commandSetResolutionLY;
  delete commandSetPenAbsorption;
  delete commandSetAbsorptionFile;
  delete commandSetReflectorOn;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorMessenger::SetNewValue(G4UIcommand* command,G4String newValue)
{
   if( command == commandSetTargetMaterial ){
	fDetector->SetTargetMaterial(newValue);
   }

   if( command == commandSetWorldMaterial ){
	fDetector->SetWorldMaterial(newValue);
   }

   if( command == commandSetDetectorName ){
	fDetector->SetDetectorName(newValue);
   }

   if( command == commandSetTargetSampleLength ){
	fDetector->SetTargetSampleLength(commandSetTargetSampleLength->GetNewDoubleValue(newValue));
   }

   if( command == commandSetTargetSampleThickness ){
	fDetector->SetTargetSampleThickness(commandSetTargetSampleThickness->GetNewDoubleValue(newValue));
   }

   if( command == commandSetTargetSampleWidth ){
	fDetector->SetTargetSampleWidth(commandSetTargetSampleWidth->GetNewDoubleValue(newValue));
   }

   if( command == commandSetDetectorType ){
	fDetector->SetDetectorType(commandSetDetectorType->GetNewIntValue(newValue));
   }
  
   if( command == commandSetNumberOfTargetSamples ){
	fDetector->SetNumberOfTargetSamples(commandSetNumberOfTargetSamples->GetNewIntValue(newValue));
   }
  
   if( command == commandSetCollimatorPositionX ){
	fDetector->SetDetectorCollimatorX(commandSetCollimatorPositionX->GetNewDoubleValue(newValue));
	//G4RunManager::GetRunManager()->ReinitializeGeometry();
   }  

   if(command == commandSetLY){
	fDetector->SetLY(commandSetLY->GetNewDoubleValue(newValue));
   }

   if(command == commandSetResolutionLY){
	fDetector->SetRes(commandSetResolutionLY->GetNewDoubleValue(newValue));
   }

   if(command == commandSetPenAbsorption){
	fDetector->SetABS(commandSetPenAbsorption->GetNewDoubleValue(newValue));
   }

   if(command == commandSetAlphaSigma){
	fDetector->SetSigAlpha(commandSetAlphaSigma->GetNewDoubleValue(newValue));
   }

   if(command == commandSetAlphaSigmaSides){
	fDetector->SetSigAlphaSides(commandSetAlphaSigmaSides->GetNewDoubleValue(newValue));
   }

   if(command == commandSetAlphaSigmaBottom){
	fDetector->SetSigAlphaBottom(commandSetAlphaSigmaBottom->GetNewDoubleValue(newValue));
   }

   if(command == commandSetPMTReflectivity){
	fDetector->SetPMTReflectivity(commandSetPMTReflectivity->GetNewDoubleValue(newValue));
   }

   if(command == commandSetAbsorptionFile){
	fDetector->SetABSFile(newValue);
   }
   
   if(command == commandSetReflectorOn){
	fDetector->SetReflectorOn(newValue);
   }
}
