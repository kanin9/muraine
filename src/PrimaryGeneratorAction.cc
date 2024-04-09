//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
//
/// \file B1/src/PrimaryGeneratorAction.cc
/// \brief Implementation of the B1::PrimaryGeneratorAction class
#include <cmath>
#include "PrimaryGeneratorAction.hh"

#include "G4LogicalVolumeStore.hh"
#include "G4LogicalVolume.hh"
#include "G4Box.hh"
#include "G4RunManager.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"


#define PI 3.14159265

  

namespace B1
{

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PrimaryGeneratorAction::PrimaryGeneratorAction()
{
  //G4int n_particle = 1;
  fEcoMug = new EcoMug();
  fParticleGun = new G4ParticleGun();
  fParticleGun->SetParticleEnergy(0.25 * GeV);
  // default particle kinematic
  //G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
  //G4String particleName;
  //G4ParticleDefinition* particle
  //  = particleTable->FindParticle(particleName="mu-");
  //fParticleGun->SetParticleDefinition(particle);
  //fParticleGun->SetParticlePosition(G4ThreeVector(0., 1.0, 0.0));
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PrimaryGeneratorAction::~PrimaryGeneratorAction()
{
  delete fParticleGun;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
  //this function is called at the begining of ecah event
  //

  // In order to avoid dependence of PrimaryGeneratorAction
  // on DetectorConstruction class we get Envelope volume
  // from G4LogicalVolumeStore.

  G4double envSizeXY = 10;
  G4double envSizeZ = 10;
  G4double size = 1.0;

  if (!fEnvelopeBox)
  {
    G4LogicalVolume* envLV
      = G4LogicalVolumeStore::GetInstance()->GetVolume("Envelope");
    if ( envLV ) fEnvelopeBox = dynamic_cast<G4Box*>(envLV->GetSolid());
  }

  if ( fEnvelopeBox ) {
    envSizeXY = fEnvelopeBox->GetXHalfLength()*2.;
    envSizeZ = fEnvelopeBox->GetZHalfLength()*2.;
  }
  else  {
    G4ExceptionDescription msg;
    msg << "Envelope volume of box shape not found.\n";
    msg << "Perhaps you have changed geometry.\n";
    msg << "The gun will be place at the center.";
    G4Exception("PrimaryGeneratorAction::GeneratePrimaries()",
     "MyCode0002",JustWarning,msg);
  }
  
  G4int intensity = 100;
  
  fEcoMug->SetUseSky();

  auto skySize = size * envSizeXY;
  fEcoMug->SetSkyCenterPosition({0, 0, 8});
  fEcoMug->SetSkySize({skySize, skySize});
  
  float fSpherePositionThetaMin = 0 * deg;
  float fSpherePositionThetaMax = 90 * deg;
  float fSpherePositionPhiMin = 0 * deg;
  float fSpherePositionPhiMax = 360 * deg;

  float fMomentumMin = 0 * GeV;
  float fMomentumMax = 1 * TeV;
  float fThetaMin = 0 * deg;
  float fThetaMax = 1 * deg;
  float fPhiMin = 0 * deg;
  float fPhiMax = 360 * deg;

  fEcoMug->SetMinimumMomentum(fMomentumMin / GeV);
  fEcoMug->SetMaximumMomentum(fMomentumMax / GeV);
  fEcoMug->SetMinimumTheta(fThetaMin / rad);
  fEcoMug->SetMaximumTheta(fThetaMax / rad);
  fEcoMug->SetMinimumPhi(fPhiMin / rad);
  fEcoMug->SetMaximumPhi(fPhiMax / rad);
  fEcoMug->SetHSphereMinPositionTheta(fSpherePositionThetaMin / rad);
  fEcoMug->SetHSphereMaxPositionTheta(fSpherePositionThetaMax / rad);
  fEcoMug->SetHSphereMinPositionPhi(fSpherePositionPhiMin / rad);
  fEcoMug->SetHSphereMaxPositionPhi(fSpherePositionPhiMax / rad);
  
  fEcoMug->Generate();
  
  fParticleGun->SetNumberOfParticles(1);
  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
  G4ParticleDefinition* muonMinus = particleTable->FindParticle("mu-");
  G4ParticleDefinition* muonPlus = particleTable->FindParticle("mu+");

  fEcoMug->GetCharge() < 0 ? fParticleGun->SetParticleDefinition(muonMinus)
                           : fParticleGun->SetParticleDefinition(muonPlus);
  
  const auto& pos = fEcoMug->GetGenerationPosition();
  //G4cout << "Generated particle at " << pos[0] << " " << pos[1] << " " << pos[2] << G4endl;
  //fParticleGun->SetParticlePosition({pos[0] * m, pos[1] * m, pos[2] * m});
  
  G4double fX = (G4UniformRand() - 0.5) * 0.1;
  G4double fZ = (G4UniformRand() - 0.5) * 0.1;
  fParticleGun->SetParticlePosition({fX*m, 0.21*m, fZ*m});

  G4ThreeVector d_cart(1, 1, 1);
  d_cart.setTheta(fEcoMug->GetGenerationTheta()); // in rad
  d_cart.setPhi(fEcoMug->GetGenerationPhi());     // in rad
  d_cart.setMag(1 * m);

  //G4cout << "Momentum is " << d_cart[0] << " " << d_cart[1] << " " << d_cart[2] << G4endl;

  fParticleGun->SetParticleMomentumDirection({0, -1,  0}); // y and z axis are swapped, somehow

  const auto& p_tot = fEcoMug->GetGenerationMomentum() * GeV;
  const auto& mu_mass = muonMinus->GetPDGMass();

  //G4double particleEnergy = 0.25 * GeV;//std::sqrt(p_tot*p_tot + mu_mass * mu_mass) - mu_mass;
  //G4cout << particleEnergy << G4endl;
  //fParticleGun->SetParticleEnergy(particleEnergy);
  fParticleGun->GeneratePrimaryVertex(anEvent);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

}


