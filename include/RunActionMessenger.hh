// #ifndef RunActionMessenger_h
// #define RunActionMessenger_h 1
//
// #include "G4UImessenger.hh"
// #include "globals.hh"
//
// class RunAction;
// class G4UIdirectory;
// class G4UIcmdWithAString;
// class G4UIcmdWithADoubleAndUnit;
// class G4UIcmdWithADouble;
// class G4UIcmdWithoutParameter;
// class G4UIcmdWithAnInteger;
//
// class RunActionMessenger: public G4UImessenger
// {
//   public:
//
//     RunActionMessenger(G4UserRunAction );
//    ~RunActionMessenger();
//
//     virtual void SetNewValue(G4UIcommand*, G4String);
//
//   private:
//
//     G4Run*      fRunAction;
//
//     G4UIdirectory*             fPENDir;
//     G4UIdirectory*             fRunDir;
//     G4UIcmdWithAString*        fRunNameCMD;
// };
//
// //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//
// #endif
