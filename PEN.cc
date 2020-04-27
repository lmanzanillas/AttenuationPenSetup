#include "G4RunManager.hh"

 #ifdef G4MULTITHREADED
 #include "G4MTRunManager.hh"
 #endif

#include "G4UImanager.hh"

#include "PhysicsList.hh"
#include "DetectorConstruction.hh"

#include "ActionInitialization.hh"

#ifdef G4VIS_USE
#include "G4VisExecutive.hh"
#endif

#ifdef G4UI_USE
#include "G4UIExecutive.hh"
#endif

#include "unistd.h"
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
namespace {
  void PrintUsage() {
    G4cerr << " Usage: " << G4endl;
    G4cerr << " OpNovice [-m macro ] [-u UIsession] [-t nThreads] [-r seed] "
           << G4endl;
    G4cerr << "   note: -t option is available only for multi-threaded mode."
           << G4endl;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

int main(int argc,char** argv)
{
    G4cout << " " << G4endl
    << "================================="<< G4endl
    << "                                 " << G4endl
    << "    #####   #####  #      #      " << G4endl
    << "   #     #  #      # #    #      " << G4endl
    << "   #     #  #      #  #   #      " << G4endl
    << "   # ####   #####  #   #  #      " << G4endl
    << "   #        #      #    # #      " << G4endl
    << "   #        #      #     ##      " << G4endl
    << "   #        #####  #      #      " << G4endl
    << "                                 " << G4endl
    << "================================="<< G4endl
    << "Luis Manzanillas, Connor Hayward "<< G4endl
    << "================================="<< G4endl
    << G4endl << G4endl;
  // Evaluate arguments
  //
  if ( argc > 9 ) {
    PrintUsage();
    return 1;
  }

  G4String macro;
  G4String session;
#ifdef G4MULTITHREADED
  G4int nThreads = 1;
#endif

  G4long myseed = 345354;
  for ( G4int i=1; i<argc; i=i+2 ) {
     if      ( G4String(argv[i]) == "-m" ) macro   = argv[i+1];
     else if ( G4String(argv[i]) == "-u" ) session = argv[i+1];
     else if ( G4String(argv[i]) == "-r" ) myseed  = atoi(argv[i+1]);
#ifdef G4MULTITHREADED
     else if ( G4String(argv[i]) == "-t" ) {
                    nThreads = G4UIcommand::ConvertToInt(argv[i+1]);
    }
#endif
    else {
      PrintUsage();
      return 1;
    }
  }

  // Choose the Random engine
  //
  G4Random::setTheEngine(new CLHEP::RanecuEngine);

  // Construct the default run manager
  //
#ifdef G4MULTITHREADED
   G4MTRunManager * runManager = new G4MTRunManager;
   if ( nThreads > 0 ) runManager->SetNumberOfThreads(nThreads);
#else
  G4RunManager * runManager = new G4RunManager;
#endif

  // Seed the random number generator manually
  G4Random::setTheSeed(myseed);

  // Initialize G4 kernel
  //

  // Set mandatory initialization classes
  //
  // Detector construction
  DetectorConstruction* det = new DetectorConstruction;
  runManager->SetUserInitialization(det);

  //
  // Physics list
  runManager-> SetUserInitialization(new PhysicsList());

  // User action initialization
  runManager->SetUserInitialization(new ActionInitialization(det));
  //runManager->SetUserInitialization(new ActionInitialization());
  //G4cout<<"All init"<<G4endl;
  runManager->Initialize();


#ifdef G4VIS_USE
  // Initialize visualization
  //
//  G4VisManager* visManager = new G4VisExecutive;
  // G4VisExecutive can take a verbosity argument - see /vis/verbose guidance.
  G4VisManager* visManager = new G4VisExecutive("Quiet");
  visManager->Initialize();
#endif

  // Get the pointer to the User Interface manager
  //
  G4UImanager* UImanager = G4UImanager::GetUIpointer();

  if ( macro.size() ) {
     // Batch mode
     G4String command = "/control/execute ";
     UImanager->ApplyCommand(command+macro);
  }
  else // Define UI session for interactive mode
  {
#ifdef G4UI_USE
    G4cout << "g4UIExecutice!" << G4endl;
     G4UIExecutive * ui = new G4UIExecutive(argc,argv,session);
#ifdef G4VIS_USE
     G4cout << "G4VisExecutice!" << G4endl;
     UImanager->ApplyCommand("/control/execute vis.mac");
#else
     UImanager->ApplyCommand("/control/execute OpNovice.in");
#endif
     if (ui->IsGUI())
        UImanager->ApplyCommand("/control/execute gui.mac");
        //UImanager->ApplyCommand("/tracking/verbose 1");
     ui->SessionStart();
     delete ui;
#endif
  }

  // Job termination
  // Free the store: user actions, physics_list and detector_description are
  //                 owned and deleted by the run manager, so they should not
  //                 be deleted in the main() program !

#ifdef G4VIS_USE
  delete visManager;
#endif
  delete runManager;

  return 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
