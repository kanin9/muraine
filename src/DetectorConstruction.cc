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
/// \file B1/src/DetectorConstruction.cc
/// \brief Implementation of the B1::DetectorConstruction class

#include "DetectorConstruction.hh"

#include "G4RunManager.hh"
#include "G4NistManager.hh"
#include "G4Box.hh"
#include "G4Colour.hh"
#include "G4Cons.hh"
#include "G4Orb.hh"
#include "G4Sphere.hh"
#include "G4Trd.hh"
#include "G4Ellipsoid.hh"
#include "G4ExtrudedSolid.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4SystemOfUnits.hh"
#include "G4VisAttributes.hh"
#include <math.h>

namespace B1
{

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* DetectorConstruction::Construct()
{
  // Get nist material manager
  G4NistManager* nist = G4NistManager::Instance();

  // Envelope parameters
  //
  G4double env_sizeXY = 0.2*m, env_sizeZ = 0.2*m;
  //G4Material* env_mat = nist->FindOrBuildMaterial("G4_WATER");

  // Option to switch on/off checking of volumes overlaps
  //
  G4bool checkOverlaps = true;

  //
  // World
  //
  G4double world_sizeXY = 1.2*env_sizeXY;
  G4double world_sizeZ  = 1.2*env_sizeZ;
  G4Material* env_mat = nist->FindOrBuildMaterial("G4_AIR");
  //G4Material* world_mat = new G4Material("G4_ICE", 0.9167*g/cm3, 2);
  G4Material* world_mat = nist->FindOrBuildMaterial("G4_AIR");
  //world_mat->AddElement(nist->FindOrBuildElement("H"), 2);
  //world_mat->AddElement(nist->FindOrBuildElement("O"), 1);
  
  
  auto solidWorld = new G4Box("World",                           // its name
    0.5 * world_sizeXY, world_sizeXY, 0.5 * world_sizeZ);  // its size

  auto logicWorld = new G4LogicalVolume(solidWorld,  // its solid
    world_mat,                                       // its material
    "World");                                        // its name

  auto physWorld = new G4PVPlacement(nullptr,  // no rotation
    G4ThreeVector(),                           // at (0,0,0)
    logicWorld,                                // its logical volume
    "World",                                   // its name
    nullptr,                                   // its mother  volume
    false,                                     // no boolean operation
    0,                                         // copy number
    checkOverlaps);                            // overlaps checking

  //
  // Envelope
  //
  auto solidEnv = new G4Box("Envelope",                    // its name
    0.5 * env_sizeXY, 0.5 * env_sizeXY, 0.5 * env_sizeZ);  // its size

  auto logicEnv = new G4LogicalVolume(solidEnv,  // its solid
    env_mat,                                     // its material
    "Envelope");                                 // its name

  new G4PVPlacement(nullptr,  // no rotation
    G4ThreeVector(),          // at (0,0,0)
    logicEnv,                 // its logical volume
    "Envelope",               // its name
    logicWorld,               // its mother  volume
    false,                    // no boolean operation
    0,                        // copy number
    checkOverlaps);           // overlaps checking

  //
  // Shape 1
  //
  /*
  G4Material* shape1_mat = nist->FindOrBuildMaterial("G4_A-150_TISSUE");
  G4ThreeVector pos1 = G4ThreeVector(0, 2*cm, -7*cm);

  // Conical section shape
  G4double shape1_rmina =  0.*cm, shape1_rmaxa = 2.*cm;
  G4double shape1_rminb =  0.*cm, shape1_rmaxb = 4.*cm;
  G4double shape1_hz = 3.*cm;
  G4double shape1_phimin = 0.*deg, shape1_phimax = 360.*deg;
  auto solidShape1 = new G4Cons("Shape1", shape1_rmina, shape1_rmaxa, shape1_rminb, shape1_rmaxb,
    shape1_hz, shape1_phimin, shape1_phimax);

  auto logicShape1 = new G4LogicalVolume(solidShape1,  // its solid
    shape1_mat,                                        // its material
    "Shape1");                                         // its name

  new G4PVPlacement(nullptr,  // no rotation
    pos1,                     // at position
    logicShape1,              // its logical volume
    "Shape1",                 // its name
    logicEnv,                 // its mother  volume
    false,                    // no boolean operation
    0,                        // copy number
    checkOverlaps); */          // overlaps checking

  //
  // Shape 2
  //
  G4Material* shape2_mat = nist->FindOrBuildMaterial("G4_SILICON_DIOXIDE");
  // Trapezoid shape
  //G4double shape2_dxa = 12*cm, shape2_dxb = 12*cm;
  //G4double shape2_dya = 0*cm, shape2_dyb = 16*cm;
  //G4double shape2_dz  = 6*cm;

  //std::vector<G4TwoVector> triangle(3);
  //triangle[0] = G4TwoVector(0*m, 0*m);
  //triangle[1] = G4TwoVector(8*m, 0*m);
  //triangle[2] = G4TwoVector(0*m, 8*m);
  //G4ExtrudedSolid* solidShape2 = new G4ExtrudedSolid("Bedrock", triangle, 4*m, G4TwoVector(0,0), 1., G4TwoVector(0,0), 1.);
   
  //auto solidShape2 = new G4Box("Bedrock", 0.49 * env_sizeXY, 0.25 * env_sizeZ, 0.49 * env_sizeXY);
  
  G4double semiX = 0.03 * env_sizeXY / m;
  G4double semiY = 0.03  * env_sizeZ  / m;
  G4double semiZ = 0.02 * env_sizeXY / m;

  //auto solidShape2 = new G4Box("Bedrock", semiX * m, semiY * m, semiZ * m);
  
  G4double waterSemiX = 0.5 * env_sizeXY / m;
  G4double waterSemiY = 0.5 * env_sizeZ / m;
  G4double waterSemiZ = 0.5 * env_sizeXY / m;
  auto waterShape = new G4Box("Water", waterSemiX * m, waterSemiY * m, waterSemiZ * m);
  
  G4Material* waterMaterial = nist->FindOrBuildMaterial("G4_WATER");
  
  auto logicWater = new G4LogicalVolume(waterShape, waterMaterial, "Water");

  logicWater->SetVisAttributes(G4VisAttributes(G4Colour(0, 0, 1, 0.3)));

  new G4PVPlacement(nullptr, G4ThreeVector(0*m, -0.1 * m + waterSemiY * m, 0*m), logicWater, "Water", logicEnv, false, 555555, checkOverlaps);

  auto solidShape2 = new G4Ellipsoid("Bedrock", semiX * m, semiY * m, semiZ * m, 0, 0);
  // auto solidShape2 = new G4Trd("Shape2",  // its name
  //  0.5 * shape2_dxa, 0.0 * shape2_dxb, 0.5 * shape2_dya, 0.5 * shape2_dyb,
  //  0.5 * shape2_dz);  // its size
	


  auto logicShape2 = new G4LogicalVolume(solidShape2,  // its solid
    shape2_mat,                                        // its material
    "Bedrock");                                         // its name
  
  logicShape2->SetVisAttributes(G4VisAttributes(
			  	G4Colour(0.32, 0.33, 0.3)
			  ));
  
  /*G4ThreeVector pos2 = G4ThreeVector(0*m, 0*m, 0*m);
  new G4PVPlacement(nullptr,
		  pos2,
		  logicShape2,
		  "Bedrock",
		  logicEnv,
		  false,
		  12,
		  checkOverlaps);
	*/
  
  G4int stoneCount = 0;
  G4double stoneVolume = semiX * semiY * semiZ * 4/3 * CLHEP::pi;
  G4double containerVolume = env_sizeXY / m * env_sizeXY / m * env_sizeZ / m;
  
  G4double longAxis = 2 * semiY;
  G4int verticalBoundary = std::floor(2 * waterSemiY / longAxis);
  //G4int verticalBoundary = 0;
  G4int boundary = std::floor(env_sizeXY / m / longAxis);
  
  /*new G4PVPlacement(nullptr,
		  	G4ThreeVector(0, -0.1*m + semiY*m , 0),
			logicShape2,
			"Bedrock",
			logicEnv,
			false,
			13,
			checkOverlaps);
  */

  G4double minX = 999, minY = 999, minZ = 999;
  G4double maxX = -1, maxY = -1, maxZ = -1;

  for(G4int i = 0; i < boundary; i++){
  	  for(G4int j = 1 + i%2; j < boundary - 1; j++){
	  	for(G4int k = 1 + i%2; k < boundary - 1; k++){
			stoneCount++;
			G4double offset = 0;
			if(i % 2){
				offset = longAxis / 2;
			}
			
			G4double ppX = 0.2 * 0.5 - semiX + offset - k * longAxis;
			G4double ppZ = 0.2 * 0.5 - semiZ + offset - j * longAxis;
			
			G4double ppY = -(waterSemiY - semiY - i * longAxis);
			
			if(i >= verticalBoundary){
				ppY = -(0.2 * 0.5 - semiY - i * longAxis);
			}

			G4cout << "ggs: " << ppX << " " << ppY << " " << ppZ << G4endl;
			
 			maxX = std::max(ppX, maxX);
			maxY = std::max(ppY, maxY);
			maxZ = std::max(ppZ, maxZ);
			minX = std::min(ppX, minX);
			minY = std::min(ppY, minY);
			minZ = std::min(ppZ, minZ);

			G4ThreeVector pos2 = G4ThreeVector(ppX * m, (ppY + 0.64 * semiY) * m, ppZ * m);
			
			G4RotationMatrix* rotation = new G4RotationMatrix();
  			rotation->rotateX(360 * (G4UniformRand() - 0.5) * deg);
  			rotation->rotateY(360 * (G4UniformRand() - 0.5) * deg);
  			rotation->rotateZ(360 * (G4UniformRand() - 0.5) * deg);
			
			auto motherVolume = logicWater;

			if(i >= verticalBoundary){
			  motherVolume = logicEnv;
			}

	  		new G4PVPlacement(rotation,  // no rotation
    				pos2,                     // at position
    				logicShape2,              // its logical volume
    				"Bedrock",                 // its name
    				motherVolume,                 // its mother  volume
    				false,                    // no boolean operation
    				i + j*100 + k*10000,                        // copy number
    				checkOverlaps);           // overlaps checking
	  	}
	}
  }
  
  containerVolume = (maxX - minX + longAxis / 2) * (maxY - minY + longAxis / 2) * (maxZ - minZ + longAxis / 2);

  G4cout << "Overall stone count is " << stoneCount << G4endl;
  G4cout << "Volume taken by stones is " << stoneCount * stoneVolume << G4endl;
  G4cout << "Volume overall is " << containerVolume << G4endl;
  G4cout << "Phantom density is " << (stoneCount * stoneVolume * 2.6 + (containerVolume - stoneVolume * stoneCount)) / containerVolume << " g/cm3" << G4endl;

  G4Box *solidMuonDetector = new G4Box("solidMuonDetector", 0.25 * env_sizeXY, 0.01 * env_sizeZ, 0.25 * env_sizeXY);
  G4LogicalVolume *logicMuonDetector = new G4LogicalVolume(solidMuonDetector, world_mat, "logicMuonDetector");
  fScoringVolume = logicMuonDetector;
  new G4PVPlacement(nullptr, G4ThreeVector(0.0, -0.52 * env_sizeZ, 0*m), logicMuonDetector, "physMuonDetector", logicWorld, false, 11, checkOverlaps);

  //
  //always return the physical World
  //
  return physWorld;
}


void DetectorConstruction::ConstructSDandField(){
  DetectorManager = G4SDManager::GetSDMpointer();

  G4String SDname = "muDetector";
  G4VSensitiveDetector *muDetector = DetectorManager->FindSensitiveDetector(SDname, false);
  
  if(!muDetector){
    muDetector = new MuonSensitiveDetector("muDetector");
    DetectorManager->AddNewDetector(muDetector);
    SetSensitiveDetector("logicMuonDetector", muDetector);
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

}
