#include "detector.hh"

MuonSensitiveDetector::MuonSensitiveDetector(G4String name) : G4VSensitiveDetector(name){
}

MuonSensitiveDetector::~MuonSensitiveDetector(){}

G4bool MuonSensitiveDetector::ProcessHits(G4Step *aStep, G4TouchableHistory *ROhist){
  G4Track *track = aStep->GetTrack();
  
  
  G4StepPoint *preStepPoint = aStep->GetPreStepPoint();
  G4StepPoint *postStepPoint = aStep->GetPostStepPoint();
  G4ThreeVector prePosition = preStepPoint->GetPosition();
  G4ThreeVector postPosition = postStepPoint->GetPosition();

  G4double dx = postPosition[0] - prePosition[0];
  G4double dy = postPosition[1] - prePosition[1];
  G4double dz = postPosition[2] - prePosition[2];

  G4String particleName = track->GetDefinition()->GetParticleName(); 
  
  if(particleName != "mu-" && particleName != "mu+") return true;

  G4cout << "Particle position: " << prePosition << " " << particleName << G4endl;

  G4cout << postPosition << G4endl;

  track->SetTrackStatus(fStopAndKill); // stop the track from further propagation to reduce duplicate detections 
  
  G4int evt = G4RunManager::GetRunManager()->GetCurrentEvent()->GetEventID();

  G4AnalysisManager *man = G4AnalysisManager::Instance();
  man->FillNtupleIColumn(0, evt);
  man->FillNtupleDColumn(1, prePosition[0]);
  man->FillNtupleDColumn(2, prePosition[1]);
  man->FillNtupleDColumn(3, prePosition[2]);
  
  man->FillNtupleDColumn(4, dx);
  man->FillNtupleDColumn(5, dy);
  man->FillNtupleDColumn(6, dz);

  man->AddNtupleRow(0);

}
