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
/// \file DetectorConstruction.cc
/// \brief Implementation of the B1::DetectorConstruction class

#include "DetectorConstruction.hh"

#include "G4RunManager.hh"
#include "G4NistManager.hh"
#include "G4Box.hh"
#include "G4Cons.hh"
#include "G4Orb.hh"
#include "G4Sphere.hh"
#include "G4Trd.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4SystemOfUnits.hh"
#include "G4Tubs.hh"
#include "G4SubtractionSolid.hh"
#include "G4Types.hh"
#include "G4VisAttributes.hh"


namespace B1
{

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::DetectorConstruction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::~DetectorConstruction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* DetectorConstruction::Construct()
{
  // Get nist material manager
  G4NistManager* nist = G4NistManager::Instance();
  
  G4Material* world_mat = nist->FindOrBuildMaterial("G4_AIR"); // the world volume is Air
  G4Material* box_mat = nist->FindOrBuildMaterial("G4_GLASS_PLATE"); // box material is water
  G4Material* box_Argon = nist->FindOrBuildMaterial("G4_Ar"); // box material is water
  G4Material* box_Lead = nist->FindOrBuildMaterial("G4_Pb"); // box material is water
  G4Material* Copper =  nist->FindOrBuildMaterial("G4_Cu"); 
  G4Material* honeyCombe =  nist->FindOrBuildMaterial("G4_PVC"); 
  // Option to switch on/off checking of volumes overlaps
  //
  G4bool checkOverlaps = true;


  // construct a world volume box of dimensions 30 cm x 30 cm x 30 cm
  G4Box* solidWorld =  new G4Box("World", 50*cm, 50*cm, 50*cm);     //Solid volume, Note convention: Half length is used
  G4LogicalVolume* logicWorld =   new G4LogicalVolume(solidWorld, world_mat,"World"); // logic volume
  // physical volume
  G4VPhysicalVolume* physWorld =
    new G4PVPlacement(0,                     //no rotation
                      G4ThreeVector(),       //at (0,0,0)
                      logicWorld,            //its logical volume
                      "World",               //its name
                      0,                     //its mother  volume
                      false,                 //no boolean operation
                      0,                     //copy number
                      checkOverlaps);        //overlaps checking


  // Big box: construct a water box of dimensions 20 cm x 20 cm x 20 cm
  G4double l1 = 26.7;
  G4double l2 = 23.7;

  G4Box* GlassPlate =  new G4Box("Box", (l1/2)*cm, (l1/2)*cm, (3/2)*mm);     //Note convention: Half length is used
  G4Box* LeadCoat =  new G4Box("Box", (23.7/2)*cm, (23.7/2)*cm, (0.3/2)*mm);     //Note convention: Half length is used
  G4Box* cornerPlate1 =  new G4Box("Box", (3.7/2)*cm, (3.7/2)*cm, (8/2)*mm);     //Note convention: Half length is used
  
  

  G4RotationMatrix* rotation = new G4RotationMatrix();
  // G4RotationMatrix* rotation2 = new G4RotationMatrix();
  // rotation2->rotateZ(90*deg)
  
  
  G4RotationMatrix* rotation1 = new G4RotationMatrix();
  rotation1->rotateZ(45*deg);

    
  G4RotationMatrix* rotation3 = new G4RotationMatrix();
  rotation3->rotateZ(90*deg);
  
  // Glass corner  cut

  G4ThreeVector trans(l1/2*cm, l1/2* cm, 0);
  G4ThreeVector trans2(-l1/2*cm, l1/2* cm, 0);
  G4ThreeVector trans3(-l1/2*cm, -l1/2* cm, 0);
  G4ThreeVector trans4(l1/2*cm, -l1/2* cm, 0);

  G4SubtractionSolid *subSolid1 = new G4SubtractionSolid("Box-Cylinder", GlassPlate, cornerPlate1, rotation1, trans);
  G4SubtractionSolid *subSolid2 = new G4SubtractionSolid("Box-Cylinder", subSolid1, cornerPlate1, rotation1, trans2);
  G4SubtractionSolid *subSolid3 = new G4SubtractionSolid("Box-Cylinder", subSolid2, cornerPlate1, rotation1, trans3);
  G4SubtractionSolid *finalPlate = new G4SubtractionSolid("Box-Cylinder", subSolid3, cornerPlate1, rotation1, trans4);

// # Lead coating 
  G4ThreeVector trans_l(l2/2*cm, l2/2* cm, 0);
  G4ThreeVector trans2_l(-l2/2*cm, l2/2* cm, 0);
  G4ThreeVector trans3_l(-l2/2*cm, -l2/2* cm, 0);
  G4ThreeVector trans4_l(l2/2*cm, -l2/2* cm, 0);

  G4SubtractionSolid *subSolid1_l = new G4SubtractionSolid("Box-Cylinder", LeadCoat, cornerPlate1, rotation1, trans_l);
  G4SubtractionSolid *subSolid2_l = new G4SubtractionSolid("Box-Cylinder", subSolid1_l, cornerPlate1, rotation1, trans2_l);
  G4SubtractionSolid *subSolid3_l = new G4SubtractionSolid("Box-Cylinder", subSolid2_l, cornerPlate1, rotation1, trans3_l);
  G4SubtractionSolid *LeadPlate = new G4SubtractionSolid("Box-Cylinder", subSolid3_l, cornerPlate1, rotation1, trans4_l);



// Definig logical volume


  G4LogicalVolume* logicBox =   new G4LogicalVolume(finalPlate, box_mat,"Glass plate");            //its name
  G4LogicalVolume* logicArgon =   new G4LogicalVolume(finalPlate, box_Argon,"Argon"); 
  G4LogicalVolume* logicLead =   new G4LogicalVolume(LeadPlate, box_Lead,"Lead plate");  
  

  G4Color yellow(1.0, 1.0, 0.0);
  G4Color red(1.0, 0.0, 0.0);
  G4Color blue(0.0, 0.0, 1.0);
  G4Color green(0.0, 1.0, 0.0);
  G4Color copper(0.7, 0.45, 0.20);

  G4VisAttributes* visAttributes1 = new G4VisAttributes();
  visAttributes1->SetColor(red);
  visAttributes1->SetVisibility(true);
  visAttributes1->SetForceSolid(true);


  G4VisAttributes* visAttributes2 = new G4VisAttributes();
  visAttributes2->SetColor(blue);
  visAttributes2->SetVisibility(true);
  visAttributes2->SetForceSolid(true);


  G4VisAttributes* visAttributes3 = new G4VisAttributes();
  visAttributes3->SetColor(yellow);
  visAttributes3->SetVisibility(true);
  visAttributes3->SetForceSolid(true);


  G4VisAttributes* visAttributes4 = new G4VisAttributes();
  visAttributes4->SetColor(green);
  visAttributes4->SetVisibility(true);
  visAttributes4->SetForceSolid(true);


  G4VisAttributes* visAttributes5 = new G4VisAttributes();
  visAttributes5->SetColor(copper);
  visAttributes5->SetVisibility(true);
  visAttributes5->SetForceSolid(true);

  logicBox->SetVisAttributes(visAttributes3);
  logicArgon ->SetVisAttributes(visAttributes2) ;
  logicLead ->SetVisAttributes(visAttributes1) ;

  // G4LogicalVolume* logicBox2 = G4LogicalVolume(logicBox);

  // Upper Glass
  G4VPhysicalVolume* physBox =
                      new G4PVPlacement(rotation,     // rotation matrix
                      G4ThreeVector(0,0,3*mm),      //at // translation
                      logicBox,             //its logical volume
                      "Glass plate",                 //its name
                      logicWorld,            //its mother  volume
                      false,                 //no boolean operation
                      0,                     //copy number
                      checkOverlaps);        //overlaps checking
// Lower glass
  G4VPhysicalVolume* physBox2 =
                      new G4PVPlacement(rotation,     // rotation matrix
                      G4ThreeVector(0,0,-3*mm),      //at // translation
                      logicBox,             //its logical volume
                      "Glass plate2",                 //its name
                      logicWorld,            //its mother  volume
                      false,                 //no boolean operation
                      0,                     //copy number
                      checkOverlaps);        //overlaps checking

// Lead upper
    G4VPhysicalVolume* LeadCoat1 =
                      new G4PVPlacement(rotation,     // rotation matrix
                      G4ThreeVector(0,0,4.65*mm),      //at // translation
                      logicLead,             //its logical volume
                      "Lead Coat ",                 //its name
                      logicWorld,            //its mother  volume
                      false,                 //no boolean operation
                      0,                     //copy number
                      checkOverlaps);        //overlaps checking

// Lead Lower
    G4VPhysicalVolume* LeadCoat2 =
                      new G4PVPlacement(rotation,     // rotation matrix
                      G4ThreeVector(0,0,-4.65*mm),      //at // translation
                      logicLead,             //its logical volume
                      "Lead Coat ",                 //its name
                      logicWorld,            //its mother  volume
                      false,                 //no boolean operation
                      0,                     //copy number
                      checkOverlaps);        //overlaps checking

// Argon
    G4VPhysicalVolume* physArgon =
                      new G4PVPlacement(rotation,     // rotation matrix
                      G4ThreeVector(0,0,0),      //at // translation
                      logicArgon,             //its logical volume
                      "Lead Coat ",                 //its name
                      logicWorld,            //its mother  volume
                      false,                 //no boolean operation
                      0,                     //copy number
                      checkOverlaps);        //overlaps checking



//  Strip and honey combe
  G4Box *pickup_pannel = new G4Box("top_plate", 13.4 * cm, 13.4* cm, 0.4 * cm);
  G4Box *hole_Box = new G4Box("hole_box", 0.8 * cm, 0.8 * cm, 0.5* cm);
  G4RotationMatrix* y_rot = new G4RotationMatrix;
  y_rot->rotateY(0*deg);

  G4ThreeVector zTrans(0, 0, 0);
    G4SubtractionSolid *subtraction1 = new G4SubtractionSolid("Box-Cylinder", pickup_pannel, hole_Box, y_rot, zTrans);
     
      for(int i = -3; i <= 3; i++)
      {
          for(int j = -3; j <= 3; j++)
          {
              G4ThreeVector zTrans(4*i*cm, 4*j*cm, 0*mm);
              subtraction1 = new G4SubtractionSolid("Box-Cylinder2", subtraction1, hole_Box, y_rot, zTrans);
          }
      }
    


      G4LogicalVolume* logicHoney =   new G4LogicalVolume(subtraction1, box_mat ,"Box-Cylinder2");            //its name
        // logicHoney ->SetVisAttributes(visAttributes4) ;

    G4VPhysicalVolume *physHoney1=
        new G4PVPlacement(0,                      // rotation matrix
                          G4ThreeVector(0, 0, 10.2*mm), // at // translation
                          logicHoney,               // its logical volume
                          "Box-Honey",         // its name
                          logicWorld,             // its mother  volume
                          false,                  // no boolean operation
                          0,                      // copy number
                          checkOverlaps);         // overlaps checking

  G4VPhysicalVolume *physHoney2 =
        new G4PVPlacement(0,                      // rotation matrix
                          G4ThreeVector(0, 0, -10.2*mm), // at // translation
                          logicHoney,               // its logical volume
                          "Box-Honey",         // its name
                          logicWorld,             // its mother  volume
                          false,                  // no boolean operation
                          0,                      // copy number
                          checkOverlaps);         // overlaps checking                        



  G4Box *box = new G4Box("Box", 13.4 * cm, 13.4 * cm, 1.4 * mm);
  // G4Tubs *cyl = new G4Tubs("Cylinder", 0 * mm, 2 * cm, 2 * cm, 0, 2 * CLHEP::pi); // rmin=0;rmax=15 mm; len=40mm and phi from 0 to 2pi
  G4Box *box1 = new G4Box("Box1", 13.5 * cm, 1 * mm, 1.5 * mm);

  G4RotationMatrix* yRot = new G4RotationMatrix;  // Rotates X and Z axes only
  yRot->rotateY(0*deg);                     // Rotates 45 degrees
  
  G4ThreeVector zzTrans(0, 0, 0);
  G4SubtractionSolid* subtraction =   new G4SubtractionSolid("Box-Cylinder", box, box1, yRot, zzTrans);
  
  for(int i=-3; i <= 3; i++){
  G4ThreeVector zzTrans(0, 3.1*i*cm, 0);
  subtraction =   new G4SubtractionSolid("Box-Cylinder", subtraction, box1, yRot, zzTrans);
         }
         
  G4LogicalVolume* logicBox1 =   new G4LogicalVolume(subtraction, box_mat,"Box-Cylinder");
  logicBox1 ->SetVisAttributes(visAttributes5) ;
  G4VPhysicalVolume *physStrip1 =
      new G4PVPlacement(0,                      // rotation matrix
                        G4ThreeVector(0, 0, -6.2*mm), // at // translation
                        logicBox1,               // its logical volume
                        "Box-Cylinder",         // its name
                        logicWorld,             // its mother  volume
                        false,                  // no boolean operation
                        0,                      // copy number
                        checkOverlaps);         // overlaps checking
  
  G4VPhysicalVolume *physStrip2 =
      new G4PVPlacement(rotation3,                      // rotation matrix
                        G4ThreeVector(0, 0, 6.2*mm), // at // translation
                        logicBox1,               // its logical volume
                        "Box-Cylinder",         // its name
                        logicWorld,             // its mother  volume
                        false,                  // no boolean operation
                        0,                      // copy number
                        checkOverlaps);         // overlaps checking
  // Set Shape2 as scoring volume
  //
  fScoringVolume = logicBox;

  //
  //always return the physical World
  //
  return physWorld;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

}
