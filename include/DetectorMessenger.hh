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
    G4UIcmdWithAString*        fTMaterCMD;
    G4UIcmdWithAString*        fWMaterCMD;
    G4UIcmdWithADoubleAndUnit* fSizeCMD;
    G4UIcmdWithAnInteger* fTypeCMD;
    G4UIcmdWithADouble* fLYCMD;
    G4UIcmdWithADouble* fResCMD;
    G4UIcmdWithADouble* fABSCMD;
    G4UIcmdWithADouble* fAlphaCMD;
    G4UIcmdWithAString* fDetNameCMD;
    G4UIcmdWithAString* fAbsFileCMD;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
