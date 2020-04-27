#ifndef PrimaryGeneratorAction_h
#define PrimaryGeneratorAction_h 1

#include "G4VUserPrimaryGeneratorAction.hh"
#include "PrimaryGeneratorMessenger.hh"
#include "globals.hh"

class G4ParticleGun;
class G4Event;
class PrimaryGeneratorMessenger;
class DetectorConstruction;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class PrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{

public:
	PrimaryGeneratorAction(DetectorConstruction*);
	~PrimaryGeneratorAction();

// public:
// 	PrimaryGeneratorAction();
// 	virtual ~PrimaryGeneratorAction();

public:
	virtual void GeneratePrimaries(G4Event*);
	//void DefineParticle();
	void  SetSourceType(G4int newType);
	void SetSourceEnergy(G4double newEnergy);
	G4int GetSourceType(void){return fSourceType;};
	G4double GetSourceEnergy(void){return fSourceEnergy;};
	void SetPhotonWavelength(G4double newValue);
	G4String GetParticleName(void){return fParticleName;};
	void SetParticleName(G4int Z, G4int A, G4double excitEnergy);
	G4double GetSourcePositionX(void){return fPositionX;};
	void SetSourcePositionX(G4double newValue);
	G4double GetSourcePositionY(void){return fPositionY;};
	void SetSourcePositionY(G4double newValue);
	G4double GetSourcePositionZ(void){return fPositionZ;};
	void SetSourcePositionZ(G4double newValue);
private:
	G4double fPositionX;
	G4double fPositionY;
	G4double fPositionZ;
	G4String fParticleName;
	PrimaryGeneratorMessenger* fPrimaryMessenger;
	G4ParticleGun*  fParticleGun;
	G4int           fSourceType;
	G4double	fSourceEnergy;
	DetectorConstruction*      fDetector;
	G4double		fPhotonWavelength;
	PrimaryGeneratorMessenger* fGunMessenger;
	G4ThreeVector fPoint;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif /*PrimaryGeneratorAction_h*/
