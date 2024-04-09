#ifndef DETECTOR_HH

#define DETECTOR_HH

#include "G4VSensitiveDetector.hh"
#include "G4RunManager.hh"
#include "G4AnalysisManager.hh"

class MuonSensitiveDetector : public G4VSensitiveDetector{
public:
	MuonSensitiveDetector(G4String);
	~MuonSensitiveDetector();
private:
	virtual G4bool ProcessHits(G4Step*, G4TouchableHistory*);
};

#endif
