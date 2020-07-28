#ifndef DetectorMessenger_h
#define DetectorMessenger_h 1

#include "G4UImessenger.hh"
#include "globals.hh"

class DetectorConstruction;
class G4UIdirectory;
class G4UIcmdWithAString;
class G4UIcmdWithADoubleAndUnit;
class G4UIcmdWithADouble;
class G4UIcmdWithoutParameter;
class G4UIcmdWithAnInteger;
class G4UIcmdWithABool;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class DetectorMessenger: public G4UImessenger
{
  public:

    DetectorMessenger(DetectorConstruction* );
   ~DetectorMessenger();

    virtual void SetNewValue(G4UIcommand*, G4String);

  private:

    DetectorConstruction*      fDetector;


    G4UIdirectory*             fPENDir;
    G4UIdirectory*             fDetDir;

    G4UIcmdWithAString*         commandSetWorldMaterial;
    G4UIcmdWithAString* 	commandSetDetectorName;
    G4UIcmdWithAnInteger*       commandSetDetectorType;
    G4UIcmdWithAnInteger*       commandSetNumberOfTargetSamples;
    G4UIcmdWithADoubleAndUnit*  commandSetTargetSampleLength;
    G4UIcmdWithADoubleAndUnit*  commandSetTargetSampleThickness;
    G4UIcmdWithADoubleAndUnit*  commandSetTargetSampleWidth;
    G4UIcmdWithAString*         commandSetTargetMaterial;
    G4UIcmdWithADouble* 	commandSetAlphaSigma;
    G4UIcmdWithADouble* 	commandSetAlphaSigmaSides;
    G4UIcmdWithADouble* 	commandSetPMTReflectivity;
    G4UIcmdWithADoubleAndUnit*  commandSetCollimatorPositionX;
    G4UIcmdWithADouble*	        commandSetLY;
    G4UIcmdWithADouble* 	commandSetResolutionLY;
    G4UIcmdWithADouble* 	commandSetPenAbsorption;
    G4UIcmdWithAString* 	commandSetAbsorptionFile;
    G4UIcmdWithABool*           commandSetReflectorOn;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
