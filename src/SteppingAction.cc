#include "SteppingAction.hh"
#include "EventAction.hh"
#include "Analysis.hh"
#include "G4SteppingManager.hh"
#include "G4SDManager.hh"
#include "G4EventManager.hh"
#include "G4ProcessManager.hh"
#include "G4Track.hh"
#include "G4Step.hh"
#include "G4Event.hh"
#include "G4StepPoint.hh"
#include "G4TrackStatus.hh"
#include "G4VPhysicalVolume.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTypes.hh"
#include "G4SystemOfUnits.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

SteppingAction::SteppingAction(EventAction* eventAction)
: fOneStepPrimaries(false), fEventAction(eventAction)
{

	fExpectedNextStatus = Undefined;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

SteppingAction::~SteppingAction() {}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void SteppingAction::UserSteppingAction(const G4Step * theStep)
{
	//auto analysisManager = G4AnalysisManager::Instance();
	G4Track* theTrack = theStep->GetTrack();
	//G4cout << "Begin Stepping" << G4endl;
	fExpectedNextStatus = Undefined;

	G4StepPoint* thePrePoint = theStep->GetPreStepPoint();
	G4VPhysicalVolume* thePrePV = thePrePoint->GetPhysicalVolume();
	//G4TouchableHistory* theTouchable = (G4TouchableHistory*)(thePrePoint->GetTouchable());
	//G4int copyNo = theTouchable->GetVolume()->GetCopyNo();


	if ( thePrePV->GetName()=="target_1" ){
		fEventAction->AddDepositedEnergyPENStackedSample1(theStep->GetTotalEnergyDeposit()/keV);
	}
	else if ( thePrePV->GetName()=="target_2" ){
		fEventAction->AddDepositedEnergyPENStackedSample2(theStep->GetTotalEnergyDeposit()/keV);
        }
	else if ( thePrePV->GetName()=="target_3" ){
		fEventAction->AddDepositedEnergyPENStackedSample3(theStep->GetTotalEnergyDeposit()/keV);
        }
	else if ( thePrePV->GetName()=="target_4" ){
		fEventAction->AddDepositedEnergyPENStackedSample4(theStep->GetTotalEnergyDeposit()/keV);
        }
	else if ( thePrePV->GetName()=="triggerFoilEJ212"){
		fEventAction->AddDepositedEnergyEJ212TriggerFoil(theStep->GetTotalEnergyDeposit()/keV);
	}else{
		fEventAction->AddDepositedEnergyInactiveMaterial(theStep->GetTotalEnergyDeposit()/keV);
	}

	G4StepPoint* thePostPoint = theStep->GetPostStepPoint();
	G4VPhysicalVolume* thePostPV = thePostPoint->GetPhysicalVolume();

	G4OpBoundaryProcessStatus boundaryStatus=Undefined;
	static G4ThreadLocal G4OpBoundaryProcess* boundary=NULL;



  	//find the boundary process only once
	if(!boundary){
		G4ProcessManager* pm	= theStep->GetTrack()->GetDefinition()->GetProcessManager();
		G4int nprocesses = pm->GetProcessListLength();
		G4ProcessVector* pv = pm->GetProcessList();
		G4int i;
		for( i=0;i<nprocesses;i++){
			if((*pv)[i]->GetProcessName()=="OpBoundary"){
				boundary = (G4OpBoundaryProcess*)(*pv)[i];
				break;
			}
		}
	}

  	if(!thePostPV)
  	{//out of world
  		fExpectedNextStatus=Undefined;
  		return;
  	}

  	G4ParticleDefinition* particleType = theTrack->GetDefinition();
  	if(particleType==G4OpticalPhoton::OpticalPhotonDefinition())
  	{
	  	boundaryStatus=boundary->GetStatus();

	    	//Check to see if the partcile was actually at a boundary
	    	//Otherwise the boundary status may not be valid
	    	//Prior to Geant4.6.0-p1 this would not have been enough to check
	  	if(thePostPoint->GetStepStatus()==fGeomBoundary){
	  		if(fExpectedNextStatus==StepTooSmall){
	  			if(boundaryStatus!=StepTooSmall){
	  				G4ExceptionDescription ed;
	  				ed << "SteppingAction::UserSteppingAction(): "
		  			<< "No reallocation step after reflection!"
		  			<< G4endl;
		  			G4Exception("SteppingAction::UserSteppingAction()", "PENExpl01",
	  				FatalException,ed,
	  				"Something is wrong with the surface normal or geometry");
	  			}
	  		}
	  		fExpectedNextStatus=Undefined;
	  		switch(boundaryStatus){
		  		case Absorption:{
					// fEventAction->AddAbsorbedPhoton();
				}
		  		break;
		      	case Detection:
		      	{
				G4String pvName = thePostPV->GetName();
				if(pvName == "trigger_pmt"){
 					fEventAction->AddDetectedPhoton();
				}
				else if(pvName == "main_pmt_1"){
					fEventAction->AddLeftPhoton();
				}
				else if(pvName == "main_pmt_2"){
					fEventAction->AddRightPhoton();
				}
				else if(pvName == "main_pmt_3"){
					fEventAction->AddBottomPhoton();
				}
				else if(pvName == "main_pmt_4"){
					fEventAction->AddFrontPhoton();
				}
				else if(pvName == "main_pmt_5"){
					fEventAction->AddBackPhoton();
				}
			      	break;
		      	}
		      	case FresnelReflection:
		      	case TotalInternalReflection:
		      	case LambertianReflection:
		      	case LobeReflection:
		      	case SpikeReflection:
		      	case BackScattering:
		      	//trackInformation->IncReflections();
		      	fExpectedNextStatus=StepTooSmall;
		      		break;
		      	default:
		      		break;
	  		}
		}
	}

}
