#include "PrimaryGeneratorAction.hh"
#include "PrimaryGeneratorMessenger.hh"

#include "G4MTRunManager.hh"

#include "Randomize.hh"
#include "DetectorConstruction.hh"
#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4SystemOfUnits.hh"
#include "G4IonTable.hh"
#include "G4DecayTable.hh"
#include "G4OpticalPhoton.hh"
#include "G4VPhysicalVolume.hh"

#include "G4Navigator.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PrimaryGeneratorAction::PrimaryGeneratorAction(DetectorConstruction* det)
: G4VUserPrimaryGeneratorAction()
{
        fDetector = det;
	fPrimaryMessenger = new PrimaryGeneratorMessenger(this);
	G4int n_particle = 1;
	fParticleGun = new G4ParticleGun(n_particle);
	//fParticleGun = new G4ParticleGun();
	fPositionZ = 30.*mm;
	//m_detectorConstruction = det;
	//m_detectorConstruction->GetMinDetectorLimits();
	//fPositionX = fDetector -> GetDetectorCollimatorX();
	fPositionX = 0. *mm;
	fSourceType = 1;
 	//default kinematic
	//G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
	//G4ParticleDefinition* particle = particleTable->FindParticle("e-");
	fSourceEnergy = 1000*keV;
	fPhotonWavelength = 0;
	fParticleName = "void";
	// fPoint = G4ThreeVector();
	//DefineParticle();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PrimaryGeneratorAction::~PrimaryGeneratorAction()
{
	delete fParticleGun;
	//delete fDetector;
}

/*
void PrimaryGeneratorAction::DefineParticle(){


	G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();


	fParticleGun->SetParticleDefinition(particleTable->FindParticle("e-"));
	fParticleGun->SetParticleEnergy(1000.*keV);

	G4ThreeVector position = G4ThreeVector(fPositionX,30.*mm, 0.*mm);
	fParticleGun->SetParticlePosition(position);
	fParticleGun->SetParticleMomentumDirection(G4ThreeVector(0,-1.,0));
}
*/
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

// Randomises placement and momentum vectors for required sources.

void PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
	G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
	//G4cout << "Gen Primaries" <<G4endl;
	G4double ionCharge   = 0.*eplus;
	G4double excitEnergy = 0.*eV;
	G4int Z=0, A=0;
	G4ParticleDefinition* ion;
	G4ThreeVector position = G4ThreeVector(fPositionX, 30*mm, 0*mm);
	//G4ThreeVector position = G4ThreeVector(fPositionX, 30*mm, 0*mm);
        //G4cout<<" position x "<<GetSourcePositionX()<<G4endl;
	//G4String name;

	switch (fSourceType) {
		case 0:
			Z = 55;
			A = 137;
			ion = G4IonTable::GetIonTable()->GetIon(Z,A,excitEnergy);
			fParticleGun->SetParticleEnergy(0.*eV);
			fParticleGun->SetParticlePosition(position);
			fParticleGun->SetParticleMomentumDirection(G4ThreeVector(0.,0.,0.));
			fParticleGun->SetParticleDefinition(ion);
			fParticleGun->SetParticleCharge(ionCharge);
			break;
		case 1:
			Z = 83;
			A = 207;
			ion = G4IonTable::GetIonTable()->GetIon(Z,A,excitEnergy);
			fParticleGun->SetParticleEnergy(0.*eV);
			fParticleGun->SetParticlePosition(position);
			fParticleGun->SetParticleMomentumDirection(G4ThreeVector(0.,0.,0.));
			fParticleGun->SetParticleDefinition(ion);
			fParticleGun->SetParticleCharge(ionCharge);
			break;
		case 2:
			Z = 38;
			A = 90;
			ion = G4IonTable::GetIonTable()->GetIon(Z,A,excitEnergy);
			fParticleGun->SetParticleEnergy(0.*eV);
			fParticleGun->SetParticlePosition(position);
			fParticleGun->SetParticleMomentumDirection(G4ThreeVector(0.,0.,0.));
			fParticleGun->SetParticleDefinition(ion);
			fParticleGun->SetParticleCharge(ionCharge);
			break;
		case 3:
			fParticleGun->SetParticleDefinition(particleTable->FindParticle("e-"));
			fParticleGun->SetParticleEnergy(2000.*keV);
			fParticleGun->SetParticlePosition(position);
			//rx = static_cast <float> (rand()) / static_cast <float> (RAND_MAX);
			//ry = static_cast <float> (rand()) / static_cast <float> (RAND_MAX);
			//rz = static_cast <float> (rand()) / static_cast <float> (RAND_MAX);
			//fParticleGun->SetParticleMomentumDirection(G4ThreeVector(2*(rx-1),2*(ry-1),2*(rz-1)));
			fParticleGun->SetParticleMomentumDirection(G4ThreeVector(0,-1,0));
 			break;

	}
	fParticleGun->GeneratePrimaryVertex(anEvent);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PrimaryGeneratorAction::SetParticleName(G4int Z, G4int A, G4double excitEnergy)
{
	fParticleName = G4IonTable::GetIonTable()->GetIonName(Z,A,excitEnergy);
}

void PrimaryGeneratorAction::SetSourcePositionZ(G4double newPosition){
	fPositionZ = newPosition;
}
void PrimaryGeneratorAction::SetSourcePositionX(G4double newPosition){
	fPositionX = newPosition;
}

void PrimaryGeneratorAction::SetSourceType(G4int newType)
{
	if (newType <= 5 && newType >= 0)
	{
		fSourceType = newType;
	}
	else
	{
		G4cerr << "The option is out of the possible values (0-5)!" << G4endl;
		G4cerr << "The default option (0) is set!" << G4endl;
		fSourceType = 0;
	}
	//DefineParticle();
}

void PrimaryGeneratorAction::SetSourceEnergy(G4double newEnergy)
{
	if (newEnergy>0)
	{
		fSourceEnergy = newEnergy;
	}
	else{
		G4cerr << "New energy is < 0." << G4endl;
		G4cerr << "The default option 60 keV is set!" << G4endl;
		fSourceEnergy = 60.*keV;
	}
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PrimaryGeneratorAction::SetPhotonWavelength(G4double newValue)
{
	if (newValue > 200 && newValue < 700)
	{
		fPhotonWavelength = newValue;
	}
	else
	{
		G4cerr << "The new desired wavelength is out of range (200-700 nm)!" << G4endl;
		G4cerr << "The photon wavelength is set to default value (420 nm)!" << G4endl;
		fPhotonWavelength = 420;
	}
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
