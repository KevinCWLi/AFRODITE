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
//      ----------------------------------------------------------------
//                          AFRODITE (iThemba Labs)
//      ----------------------------------------------------------------
//
//      Github repository: https://www.github.com/KevinCWLi/AFRODITE
//
//      Main Author:    K.C.W. Li
//
//      email: likevincw@gmail.com
//

#include "DetectorConstruction.hh"

#include "G4NistManager.hh"
#include "G4Box.hh"
#include "G4Sphere.hh"
#include "G4LogicalVolume.hh"
#include "G4RotationMatrix.hh"
#include "G4Transform3D.hh"

#include "G4SubtractionSolid.hh"
#include "G4UnionSolid.hh"
#include "G4IntersectionSolid.hh"
#include "G4Trd.hh"

#include "G4Material.hh"
#include "G4MaterialPropertiesTable.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4SDManager.hh"
#include "G4GeometryTolerance.hh"
#include "G4GeometryManager.hh"
#include "G4Isotope.hh"
#include "G4OpticalSurface.hh"
#include "G4LogicalBorderSurface.hh"
#include "G4LogicalSkinSurface.hh"

#include "G4UserLimits.hh"
#include "G4VisAttributes.hh"
#include "G4Colour.hh"
#include "G4SystemOfUnits.hh"
#include "globals.hh"

#include "G4GlobalMagFieldMessenger.hh"
#include "G4AutoDelete.hh"

#include "G4PhysicalVolumeStore.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4SolidStore.hh"
#include "G4PhysicalConstants.hh"

#include "G4UniformMagField.hh"
#include "G4QuadrupoleMagField.hh"
#include "G4Mag_UsualEqRhs.hh"
#include "G4ClassicalRK4.hh"

#include "G4MagneticField.hh"
#include "G4FieldManager.hh"
#include "G4TransportationManager.hh"
#include "G4Mag_UsualEqRhs.hh"
#include "G4AutoDelete.hh"

#include "CADMesh.hh"
#include "MagneticFieldMapping.hh"
//#include "G4BlineTracer.hh"

#include <fstream>
#include <string>

using namespace std;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4ThreadLocal G4GlobalMagFieldMessenger* DetectorConstruction::fMagFieldMessenger = 0;

G4ThreadLocal G4QuadrupoleMagField* DetectorConstruction::MagneticField_AFRODITE_Q = 0;
G4ThreadLocal G4UniformMagField* DetectorConstruction::MagneticField_AFRODITE_D1 = 0;

G4ThreadLocal G4FieldManager* DetectorConstruction::fieldManagerMagneticField_AFRODITE_Q = 0;
G4ThreadLocal G4FieldManager* DetectorConstruction::fieldManagerMagneticField_AFRODITE_D1 = 0;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::DetectorConstruction()
 : G4VUserDetectorConstruction(),
fAbsorberPV(0), fGapPV(0), fCheckOverlaps(false), PhysiCLOVER_HPGeCrystal(0), PhysiCLOVER_Shield_BGOCrystal(0), PhysiCLOVER_Shield_PMT(0), PhysiPlasticScint(0), PhysiHAGAR_NaICrystal(0), PhysiHAGAR_Annulus(0), PhysiHAGAR_FrontDisc(0), PhysiAFRODITE_Dipole1(0), PhysiAFRODITE_Quadrupole(0)
{
    WorldSize = 15.*m;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::~DetectorConstruction()
{ 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* DetectorConstruction::Construct()
{
    
    //////////////////////////////////////////////////////////
    ////                                                  ////
    ////            DETECTOR ARRAY POSITIONING            ////
    ////                                                  ////
    //////////////////////////////////////////////////////////
    
    ////////////////////////////
    ////    CLOVER SETUP
    
    CLOVER_AllPresent_Override = false;
    CLOVER_AllAbsent_Override = true;
    
    CLOVER_Shield_AllPresent_Override = false;
    CLOVER_Shield_AllAbsent_Override = true;
    
    
    //  CLOVER 1
    CLOVER_Presence[0] = true;
    CLOVER_Shield_Presence[0] = true;
    CLOVER_Distance[0] = 0*cm;
    CLOVER_phi[0] = 90*deg;
    CLOVER_theta[0] = 135*deg;
    CLOVER_rotm[0].rotateX(45.*deg);
    
    //  CLOVER 2
    CLOVER_Presence[1] = false;
    CLOVER_Shield_Presence[1] = false;
    CLOVER_Distance[1] = 10*cm;
    CLOVER_phi[1] = 0*deg;
    CLOVER_theta[1] = 135*deg;
    CLOVER_rotm[1].rotateY(-45.0*deg);
    
    //  CLOVER 3
    CLOVER_Presence[2] = false;
    CLOVER_Shield_Presence[2] = false;
    CLOVER_Distance[2] = 10*cm;
    CLOVER_phi[2] = 270*deg;
    CLOVER_theta[2] = 135*deg;
    CLOVER_rotm[2].rotateX(-45.0*deg);
    
    //  CLOVER 4
    CLOVER_Presence[3] = false;
    CLOVER_Shield_Presence[3] = false;
    CLOVER_Distance[3] = 10*cm;
    CLOVER_phi[3] = 180*deg;
    CLOVER_theta[3] = 135*deg;
    CLOVER_rotm[3].rotateY(45.0*deg);
    
    //  CLOVER 5
    CLOVER_Presence[4] = false;
    CLOVER_Shield_Presence[4] = false;
    CLOVER_Distance[4] = 10*cm;
    CLOVER_phi[4] = 45*deg;
    CLOVER_theta[4] = 90*deg;
    CLOVER_rotm[4].rotateY(90.0 *deg);
    CLOVER_rotm[4].rotateZ(-135.0*deg);
    
    //  CLOVER 6
    CLOVER_Presence[5] = false;
    CLOVER_Shield_Presence[5] = false;
    CLOVER_Distance[5] = 10*cm;
    CLOVER_phi[5] = 0*deg;
    CLOVER_theta[5] = 90*deg;
    CLOVER_rotm[5].rotateY(-90.0*deg);
    
    //  CLOVER 7
    CLOVER_Presence[6] = false;
    CLOVER_Shield_Presence[6] = false;
    CLOVER_Distance[6] = 10*cm;
    CLOVER_phi[6] = 180*deg;
    CLOVER_theta[6] = 90*deg;
    CLOVER_rotm[6].rotateY(90.0*deg);
    
    //  CLOVER 8
    CLOVER_Presence[7] = false;
    CLOVER_Shield_Presence[7] = false;
    CLOVER_Distance[7] = 10*cm;
    CLOVER_phi[7] = 135*deg;
    CLOVER_theta[7] = 90*deg;
    CLOVER_rotm[7].rotateY(90.0*deg);
    CLOVER_rotm[7].rotateZ(-45.0*deg);
    
    
    //  CLOVER 9
    CLOVER_Presence[8] = false;
    CLOVER_Shield_Presence[8] = false;
    CLOVER_Distance[8] = 10*cm;
    CLOVER_phi[8] = -45*deg;
    CLOVER_theta[8] = 90*deg;
    CLOVER_rotm[8].rotateY(-90.0 *deg);
    CLOVER_rotm[8].rotateZ(-45.0 *deg);
    
    
    for (G4int i=0; i<numberOf_CLOVER; i++)
    {
        if(CLOVER_AllPresent_Override) CLOVER_Presence[i] = true;
        if(CLOVER_AllAbsent_Override) CLOVER_Presence[i] = false;
        if(CLOVER_AllPresent_Override && CLOVER_AllAbsent_Override) CLOVER_Presence[i] = false;
        
        if(CLOVER_Shield_AllPresent_Override) CLOVER_Shield_Presence[i] = true;
        if(CLOVER_Shield_AllAbsent_Override) CLOVER_Shield_Presence[i] = false;
        if(CLOVER_Shield_AllPresent_Override && CLOVER_Shield_AllAbsent_Override) CLOVER_Shield_Presence[i] = false;
    }
    
    
    /////////////////////////////////////////
    ////        PlasticScint SETUP
    
    PlasticScint_AllPresent_Override = true;
    PlasticScint_AllAbsent_Override = false;
    
    
    //  PlasticScint 1
    PlasticScint_Presence[0] = true;
    PlasticScint_CentrePositionX[0] = 0.; // cm
    PlasticScint_CentrePositionY[0] = 28.75; // cm
    PlasticScint_CentrePositionZ[0] = 180.; // cm
    PlasticScint_CentrePosition[0] = G4ThreeVector(PlasticScint_CentrePositionX[0]*cm, PlasticScint_CentrePositionY[0]*cm, PlasticScint_CentrePositionZ[0]*cm);
    //PlasticScint_RotationY[0] = -14.03; // deg
    //PlasticScint_rotm[0].rotateY(PlasticScint_RotationY[0]*deg);
    
    //  PlasticScint 2
    PlasticScint_Presence[1] = true;
    PlasticScint_CentrePositionX[1] = 0.; // cm
    PlasticScint_CentrePositionY[1] = 17.25; // cm
    PlasticScint_CentrePositionZ[1] = 180.; // cm
    PlasticScint_CentrePosition[1] = G4ThreeVector(PlasticScint_CentrePositionX[1]*cm, PlasticScint_CentrePositionY[1]*cm, PlasticScint_CentrePositionZ[1]*cm);
    //PlasticScint_RotationY[1] = -14.13; // deg
    //PlasticScint_rotm[1].rotateY(PlasticScint_RotationY[1]*deg);
    
    //  PlasticScint 3
    PlasticScint_Presence[2] = true;
    PlasticScint_CentrePositionX[2] = 0.; // cm
    PlasticScint_CentrePositionY[2] = 5.75; // cm
    PlasticScint_CentrePositionZ[2] = 180.; // cm
    PlasticScint_CentrePosition[2] = G4ThreeVector(PlasticScint_CentrePositionX[2]*cm, PlasticScint_CentrePositionY[2]*cm, PlasticScint_CentrePositionZ[2]*cm);
    //PlasticScint_RotationY[2] = -14.23; // deg
    //PlasticScint_rotm[2].rotateY(PlasticScint_RotationY[2]*deg);
    
    //  PlasticScint 4
    PlasticScint_Presence[3] = true;
    PlasticScint_CentrePositionX[3] = 0.; // cm
    PlasticScint_CentrePositionY[3] = -5.75; // cm
    PlasticScint_CentrePositionZ[3] = 180.; // cm
    PlasticScint_CentrePosition[3] = G4ThreeVector(PlasticScint_CentrePositionX[3]*cm, PlasticScint_CentrePositionY[3]*cm, PlasticScint_CentrePositionZ[3]*cm);
    //PlasticScint_RotationY[3] = -14.33; // deg
    //PlasticScint_rotm[3].rotateY(PlasticScint_RotationY[3]*deg);
    
    //  PlasticScint 5
    PlasticScint_Presence[4] = true;
    PlasticScint_CentrePositionX[4] = 0.; // cm
    PlasticScint_CentrePositionY[4] = -17.25; // cm
    PlasticScint_CentrePositionZ[4] = 180.; // cm
    PlasticScint_CentrePosition[4] = G4ThreeVector(PlasticScint_CentrePositionX[4]*cm, PlasticScint_CentrePositionY[4]*cm, PlasticScint_CentrePositionZ[4]*cm);
    //PlasticScint_RotationY[4] = -14.43; // deg
    //PlasticScint_rotm[4].rotateY(PlasticScint_RotationY[4]*deg);
    
    //  PlasticScint 6
    PlasticScint_Presence[5] = true;
    PlasticScint_CentrePositionX[5] = 0.; // cm
    PlasticScint_CentrePositionY[5] = -28.75; // cm
    PlasticScint_CentrePositionZ[5] = 180.; // cm
    PlasticScint_CentrePosition[5] = G4ThreeVector(PlasticScint_CentrePositionX[5]*cm, PlasticScint_CentrePositionY[5]*cm, PlasticScint_CentrePositionZ[5]*cm);
    //PlasticScint_RotationY[5] = -14.53; // deg
    //PlasticScint_rotm[5].rotateY(PlasticScint_RotationY[5]*deg);
    
    //  PlasticScint 7
    PlasticScint_Presence[6] = true;
    PlasticScint_CentrePositionX[6] = 0.; // cm
    PlasticScint_CentrePositionY[6] = 28.75; // cm
    PlasticScint_CentrePositionZ[6] = 180. + 11.0; // cm
    PlasticScint_CentrePosition[6] = G4ThreeVector(PlasticScint_CentrePositionX[6]*cm, PlasticScint_CentrePositionY[6]*cm, PlasticScint_CentrePositionZ[6]*cm);
    //PlasticScint_RotationY[6] = -14.63; // deg
    //PlasticScint_rotm[6].rotateY(PlasticScint_RotationY[6]*deg);
    
    //  PlasticScint 8
    PlasticScint_Presence[7] = true;
    PlasticScint_CentrePositionX[7] = 0.; // cm
    PlasticScint_CentrePositionY[7] = 17.25; // cm
    PlasticScint_CentrePositionZ[7] = 180. + 11.0; // cm
    PlasticScint_CentrePosition[7] = G4ThreeVector(PlasticScint_CentrePositionX[7]*cm, PlasticScint_CentrePositionY[7]*cm, PlasticScint_CentrePositionZ[7]*cm);
    //PlasticScint_RotationY[7] = -14.73; // deg
    //PlasticScint_rotm[7].rotateY(PlasticScint_RotationY[7]*deg);
    
    //  PlasticScint 9
    PlasticScint_Presence[8] = true;
    PlasticScint_CentrePositionX[8] = 0.; // cm
    PlasticScint_CentrePositionY[8] = 5.75; // cm
    PlasticScint_CentrePositionZ[8] = 180. + 11.0; // cm
    PlasticScint_CentrePosition[8] = G4ThreeVector(PlasticScint_CentrePositionX[8]*cm, PlasticScint_CentrePositionY[8]*cm, PlasticScint_CentrePositionZ[8]*cm);
    //PlasticScint_RotationY[8] = -14.83; // deg
    //PlasticScint_rotm[8].rotateY(PlasticScint_RotationY[8]*deg);
    
    //  PlasticScint 10
    PlasticScint_Presence[9] = true;
    PlasticScint_CentrePositionX[9] = 0.; // cm
    PlasticScint_CentrePositionY[9] = -5.75; // cm
    PlasticScint_CentrePositionZ[9] = 180. + 11.0; // cm
    PlasticScint_CentrePosition[9] = G4ThreeVector(PlasticScint_CentrePositionX[9]*cm, PlasticScint_CentrePositionY[9]*cm, PlasticScint_CentrePositionZ[9]*cm);
    //PlasticScint_RotationY[9] = -14.93; // deg
    //PlasticScint_rotm[9].rotateY(PlasticScint_RotationY[9]*deg);
    
    //  PlasticScint 11
    PlasticScint_Presence[10] = true;
    PlasticScint_CentrePositionX[10] = 0.; // cm
    PlasticScint_CentrePositionY[10] = -17.25; // cm
    PlasticScint_CentrePositionZ[10] = 180. + 11.0; // cm
    PlasticScint_CentrePosition[10] = G4ThreeVector(PlasticScint_CentrePositionX[10]*cm, PlasticScint_CentrePositionY[10]*cm, PlasticScint_CentrePositionZ[10]*cm);
    //PlasticScint_RotationY[10] = -14.103; // deg
    //PlasticScint_rotm[10].rotateY(PlasticScint_RotationY[10]*deg);
    
    //  PlasticScint 12
    PlasticScint_Presence[11] = true;
    PlasticScint_CentrePositionX[11] = 0.; // cm
    PlasticScint_CentrePositionY[11] = -28.75; // cm
    PlasticScint_CentrePositionZ[11] = 180. + 11.0; // cm
    PlasticScint_CentrePosition[11] = G4ThreeVector(PlasticScint_CentrePositionX[11]*cm, PlasticScint_CentrePositionY[11]*cm, PlasticScint_CentrePositionZ[11]*cm);
    //PlasticScint_RotationY[11] = -14.113; // deg
    //PlasticScint_rotm[11].rotateY(PlasticScint_RotationY[11]*deg);
    
    for(G4int i=0; i<numberOf_PlasticScint; i++)
    {
        if(PlasticScint_AllPresent_Override) PlasticScint_Presence[i] = true;
        if(PlasticScint_AllAbsent_Override) PlasticScint_Presence[i] = false;
        if(PlasticScint_AllPresent_Override && PlasticScint_AllAbsent_Override) PlasticScint_Presence[i] = false;
    }

    
    ////////////////////////////////
    ////        LEPS SETUP
    
    LEPS_AllPresent_Override = true;
    LEPS_AllAbsent_Override = false;
    
    
    //  LEPS 1
    LEPS_Presence[0] = true;
    LEPS_Distance[0] = 4.5*cm;
    LEPS_phi[0] = 0.*deg;
    LEPS_theta[0] = 0.*deg;
    LEPS_rotm[0].rotateX(180.*deg);
    
    //  LEPS 2
    LEPS_Presence[1] = true;
    LEPS_Distance[1] = 4.5*cm;
    LEPS_phi[1] = 0.*deg;
    LEPS_theta[1] = 90.*deg;
    LEPS_rotm[1].rotateY(-90.*deg);
    
    //  LEPS 3
    LEPS_Presence[2] = true;
    LEPS_Distance[2] = 4.5*cm;
    LEPS_phi[2] = 90.*deg;
    LEPS_theta[2] = 90.*deg;
    LEPS_rotm[2].rotateX(90.*deg);
    
    //  LEPS 4
    LEPS_Presence[3] = true;
    LEPS_Distance[3] = 4.5*cm;
    LEPS_phi[3] = 180.*deg;
    LEPS_theta[3] = 90*deg;
    LEPS_rotm[3].rotateY(90.*deg);
    
    //  LEPS 5
    LEPS_Presence[4] = true;
    LEPS_Distance[4] = 4.5*cm;
    LEPS_phi[4] = 270.*deg;
    LEPS_theta[4] = 90*deg;
    LEPS_rotm[4].rotateX(-90.*deg);
    
    
    //  LEPS 6
    LEPS_Presence[5] = true;
    LEPS_Distance[5] = 4.5*cm;
    LEPS_phi[5] = 0*deg;
    LEPS_theta[5] = 180*deg;
    LEPS_rotm[5].rotateY(0.*deg);
    
    
    
    for (G4int i=0; i<numberOf_LEPS; i++)
    {
        if( LEPS_AllPresent_Override == true ) LEPS_Presence[i] = true;
        if( LEPS_AllAbsent_Override == true ) LEPS_Presence[i] = false;
        if( LEPS_AllPresent_Override == true && LEPS_AllAbsent_Override == true ) LEPS_Presence[i] = false;
    }

    
    
    ////////////////////////////
    ////    HAGAR SETUP
    
    //  HAGAR NaI Crystal
    HAGAR_NaICrystal_Presence = false;
    HAGAR_NaICrystal_CentrePosition = G4ThreeVector(-100.*cm, 0.*cm, 0.*cm);
    HAGAR_rotm.rotateY(-90*deg);
    
    //  HAGAR Annulus
    HAGAR_Annulus_Presence = false;
    HAGAR_Annulus_CentrePosition = G4ThreeVector(-105.*cm, 0.*cm, 0.*cm);
    
    //  HAGAR Front
    HAGAR_FrontDisc_Presence = false;
    HAGAR_FrontDisc_CentrePosition = HAGAR_Annulus_CentrePosition + G4ThreeVector((61./2 +8/2)*cm, 0.*cm, 0.*cm);
    
    

    ////////////////////////////////////////////////////
    ////                                            ////
    ////            AFRODITE VAULT SETUP            ////
    ////                                            ////
    ////////////////////////////////////////////////////
    
    //  AFRODITE Quadrupole
    AFRODITE_Quadrupole = false;
    Ideal_Quadrupole = true;
    Mapped_Quadrupole = false;
    //dBdr_gradient_AFRODITE_Q = 0.030*tesla/cm;  // gradient = dB/dr for Ideal Quadrupole
    dBdr_gradient_AFRODITE_Q = 0.001*tesla/cm;  // gradient = dB/dr for Ideal Quadrupole
    AFRODITE_Quadrupole_CentrePosition = G4ThreeVector(0.*cm, 0.*cm,100.*cm);
    //AFRODITE_Quadrupole_rotm.rotateZ(90.*deg);
    //AFRODITE_Quadrupole_rotm.rotateY(90*deg);
    //AFRODITE_Quadrupole_rotm.rotateX(90*deg);
    
    //  AFRODITE Dipole 1
    AFRODITE_Dipole1 = false;
    AFRODITE_Dipole1_BZ = -2.30*tesla;
    //AFRODITE_Dipole1_BZ = -3.30*tesla;
    AFRODITE_Dipole1_CentrePosition = G4ThreeVector(75*cm, 0.*cm,250*cm);
    AFRODITE_Dipole1_rotm.rotateX(-90.*deg);
    AFRODITE_Dipole1_rotm.rotateY(180.*deg);
    
    if(Ideal_Quadrupole && Mapped_Quadrupole || !Ideal_Quadrupole && !Mapped_Quadrupole)
    {
        Ideal_Quadrupole = false;
        Mapped_Quadrupole = false;
    }
    
    ////////////////////////////////////////////////////////
    ////                                                ////
    ////        AFRODITE VAULT - STRUCTURE SETUP        ////
    ////                                                ////
    ////////////////////////////////////////////////////////
    
    
    ////////////////////////////////////////////////
    ////    New AFRODITE Target Chamber by Mathis
    AFRODITE_MathisTC_Presence = false;
    
    /////////////////////////////////////
    ////    AFRODITE Target
    AFRODITE_Target_Presence = true;
    
    
  // Define materials 
  DefineMaterials();
  
  // Define volumes
  return DefineVolumes();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::DefineMaterials()
{
    G4NistManager* nistManager = G4NistManager::Instance();
    
    //  NIST Material Database - Materials
    nistManager->FindOrBuildMaterial("G4_Galactic");
    nistManager->FindOrBuildMaterial("G4_AIR");
    nistManager->FindOrBuildMaterial("G4_BGO");
    nistManager->FindOrBuildMaterial("G4_Ge");
    nistManager->FindOrBuildMaterial("G4_Al");
    nistManager->FindOrBuildMaterial("G4_Si");
    nistManager->FindOrBuildMaterial("G4_PLASTIC_SC_VINYLTOLUENE");
    nistManager->FindOrBuildMaterial("G4_MYLAR");
    nistManager->FindOrBuildMaterial("G4_W");
    nistManager->FindOrBuildMaterial("G4_Ar");
    nistManager->FindOrBuildMaterial("G4_Be");
    nistManager->FindOrBuildMaterial("G4_CARBON_DIOXIDE");
    nistManager->FindOrBuildMaterial("G4_SODIUM_IODIDE");
    nistManager->FindOrBuildMaterial("G4_LITHIUM_CARBONATE");

    //  NIST Elementary Material Database - ELEMENTS
    nistManager->FindOrBuildElement("H");
    nistManager->FindOrBuildElement("C");
    nistManager->FindOrBuildElement("N");
    nistManager->FindOrBuildElement("O");
    nistManager->FindOrBuildElement("Fe");
    nistManager->FindOrBuildElement("Co");
    nistManager->FindOrBuildElement("Ni");
    nistManager->FindOrBuildElement("Cu");
    nistManager->FindOrBuildElement("Pb");
    nistManager->FindOrBuildElement("W");
    nistManager->FindOrBuildElement("Li");


    // Vacuum
    //new G4Material("Galactic", z=1., a=1.01*g/mole,density= universe_mean_density,kStateGas, 2.73*kelvin, 3.e-18*pascal);
    
    
    // Print materials
  G4cout << *(G4Material::GetMaterialTable()) << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* DetectorConstruction::DefineVolumes()
{
    //////////////////////////////////////
    //          Get Elements            //
    //////////////////////////////////////
    
    G4Element* Hydrogen = G4Element::GetElement("H");
    G4Element* Carbon = G4Element::GetElement("C");
    G4Element* Nitrogen = G4Element::GetElement("N");
    G4Element* Oxygen = G4Element::GetElement("O");
    G4Element* Iron = G4Element::GetElement("Fe");
    G4Element* Cobalt = G4Element::GetElement("Co");
    G4Element* Nickel = G4Element::GetElement("Ni");
    G4Element* Copper = G4Element::GetElement("Cu");
    G4Element* Lead = G4Element::GetElement("Pb");
    G4Element* Tungsten = G4Element::GetElement("W");
    G4Element* Lithium = G4Element::GetElement("Li");

    
    //////////////////////////////////////
    //          Get Materials           //
    //////////////////////////////////////
    
    ////    NIST Defined Elemental Material
    G4Material* G4_Ge_Material = G4Material::GetMaterial("G4_Ge");
    G4Material* G4_Al_Material  = G4Material::GetMaterial("G4_Al");
    G4Material* G4_Si_Material = G4Material::GetMaterial("G4_Si");
    G4Material* G4_W_Material = G4Material::GetMaterial("G4_W");
    G4Material* G4_Ar_Material = G4Material::GetMaterial("G4_Ar");
    G4Material* G4_Be_Material = G4Material::GetMaterial("G4_Be");

    ////    NIST Defined Materials and Compounds
    G4Material* G4_Galactic_Material = G4Material::GetMaterial("G4_Galactic");
    G4Material* G4_AIR_Material = G4Material::GetMaterial("G4_AIR");
    G4Material* G4_BGO_Material = G4Material::GetMaterial("G4_BGO");
    G4Material* G4_PLASTIC_SC_VINYLTOLUENE_Material = G4Material::GetMaterial("G4_PLASTIC_SC_VINYLTOLUENE");
    G4Material* G4_Mylar_Material = G4Material::GetMaterial("G4_MYLAR");
    G4Material* G4_CARBON_DIOXIDE_Material = G4Material::GetMaterial("G4_CARBON_DIOXIDE");
    G4Material* G4_SODIUM_IODIDE_Material = G4Material::GetMaterial("G4_SODIUM_IODIDE");
    G4Material* G4_LITHIUM_CARBONATE_Material = G4Material::GetMaterial("G4_LITHIUM_CARBONATE");

    ////    CLOVER Detector Shield, HEAVIMET Material
    G4Material* Heavimet_Material = new G4Material("Heavimet_Material",19.25*g/cm3, 5);
    Heavimet_Material->AddElement( Tungsten, 94.20*perCent);
    Heavimet_Material->AddElement( Nickel, 4.35*perCent);
    Heavimet_Material->AddElement( Iron, 0.85*perCent);
    Heavimet_Material->AddElement( Cobalt, 0.50*perCent);
    Heavimet_Material->AddElement( Copper, 0.10*perCent);
    
    
    ////    PlasticScint Detector Materials
    ////    This Plastic is identical in composition to NE102, which is the plastic utilized in the Neutron Wall plastic Scintillators which have been employed within the AFRODITE vault for numerous (3He,n) experiments
    G4Material* BC408_Material = new G4Material("BC408_Material",1.032*g/cm3, 2);
    BC408_Material->AddElement( Hydrogen, 8.4748*perCent);
    BC408_Material->AddElement( Carbon, 91.5252*perCent);
    
    ////    VDC Detector Materials
    G4Material* VDC_SR_Gas_Material = new G4Material("VDC_SR_Gas_Material", 0.00184212*g/cm3, 2);
    VDC_SR_Gas_Material->AddMaterial( G4_Ar_Material, 90.*perCent);
    VDC_SR_Gas_Material->AddMaterial( G4_CARBON_DIOXIDE_Material, 10.*perCent);
    
    
    
    //////////////////////////////////////////////////////////
    //                      WORLD                           //
    //////////////////////////////////////////////////////////
    
    G4Box* SolidWorld = new G4Box("World", WorldSize/2, WorldSize/2, WorldSize/2);
    
    G4LogicalVolume*
    LogicWorld = new G4LogicalVolume(SolidWorld,                //its solid
                                     G4_Galactic_Material,      //its material
                                     "World");                  //its name
    G4VPhysicalVolume*
    PhysiWorld = new G4PVPlacement(0,                        //no rotation
                                   G4ThreeVector(),          //at (0,0,0)
                                   LogicWorld,               //its logical volume
                                   "World",                  //its name
                                   0,                        //its mother  volume
                                   false,                    //no boolean operation
                                   0);                       //copy number

    
    
    //////////////////////////////////////////////////////////
    //                      VACUUM CHAMBER
    //////////////////////////////////////////////////////////
    
    G4ThreeVector positionVacuumChamber = G4ThreeVector(0,0,0);
    
    G4Box* SolidVacuumChamber = new G4Box("VacuumChamber", (100./2)*cm, (100./2)*cm, (100./2)*cm);
    
    G4LogicalVolume* LogicVacuumChamber = new G4LogicalVolume(SolidVacuumChamber, G4_Galactic_Material,"VacuumChamber",0,0,0);
    
    
    new G4PVPlacement(0,               // no rotation
                      positionVacuumChamber, // at (x,y,z)
                      LogicVacuumChamber,       // its logical volume
                      "VacuumChamber",       // its name
                      LogicWorld,         // its mother  volume
                      false,           // no boolean operations
                      0,               // copy number
                      fCheckOverlaps); // checking overlaps
    
    
    
    
    
    //////////////////////////////////////////////////////////
    //              Scattering Chamber - CADMesh
    //////////////////////////////////////////////////////////
    
    if(AFRODITE_MathisTC_Presence)
    {
        G4ThreeVector offset_MathisTC = G4ThreeVector(0*cm, 0*cm, 0*cm);
        
        CADMesh * mesh_MathisTC = new CADMesh("../AFRODITE/Mesh-Models/STRUCTURES/MathisTC/MathisTC.ply", "PLY", mm, offset_MathisTC, false);
        
        G4VSolid * SolidMathisTC = mesh_MathisTC->TessellatedMesh();
        
        G4LogicalVolume* LogicMathisTC = new G4LogicalVolume(SolidMathisTC, G4_Al_Material, "BACTAR", 0, 0, 0);
        
        new G4PVPlacement(0,               // no rotation
                          G4ThreeVector(), // at (x,y,z)
                          LogicMathisTC,       // its logical volume
                          "BACTAR",       // its name
                          LogicVacuumChamber,         // its mother  volume
                          false,           // no boolean operations
                          0,               // copy number
                          fCheckOverlaps); // checking overlaps
        
        
        //  Visualisation
        G4VisAttributes* AFRODITE_MathisTC_VisAtt = new G4VisAttributes(G4Colour(0.8,0.8,0.8));
        AFRODITE_MathisTC_VisAtt->SetForceSolid(true);
        LogicMathisTC->SetVisAttributes(AFRODITE_MathisTC_VisAtt);
    }

    
    //////////////////////////////////////////////////////////
    //                      TARGET
    //////////////////////////////////////////////////////////
    
    G4ThreeVector positionTarget = G4ThreeVector(0,0,0);
    
    G4Box* SolidTarget = new G4Box("Target", (20./2)*mm, (20./2)*mm, (2.42/2)*um);
    
    G4LogicalVolume* LogicTarget = new G4LogicalVolume(SolidTarget, G4_LITHIUM_CARBONATE_Material,"Target",0,0,0);
    
    if(AFRODITE_Target_Presence)
    {
        new G4PVPlacement(0,               // no rotation
                          positionTarget, // at (x,y,z)
                          LogicTarget,       // its logical volume
                          "Target",       // its name
                          LogicVacuumChamber,         // its mother  volume
                          false,           // no boolean operations
                          0,               // copy number
                          fCheckOverlaps); // checking overlaps
    }
    

    
    
    /////////////////////////////////////////////
    //             CLOVER DEFINITION           //
    /////////////////////////////////////////////
    
    G4double CLOVERtoShield_displacement = 10.;  // cm
    
    G4ThreeVector offset_CLOVERInternalVacuum = G4ThreeVector(0*cm, 0*cm, -CLOVERtoShield_displacement*cm);
    G4ThreeVector offset_CLOVEREncasement = G4ThreeVector(0*cm, 0*cm, -CLOVERtoShield_displacement*cm);
    G4ThreeVector offset_CLOVERHPGeCrystal1 = G4ThreeVector(0*cm, 0*cm, -CLOVERtoShield_displacement*cm);
    G4ThreeVector offset_CLOVERHPGeCrystal2 = G4ThreeVector(0*cm, 0*cm, -CLOVERtoShield_displacement*cm);
    G4ThreeVector offset_CLOVERHPGeCrystal3 = G4ThreeVector(0*cm, 0*cm, -CLOVERtoShield_displacement*cm);
    G4ThreeVector offset_CLOVERHPGeCrystal4 = G4ThreeVector(0*cm, 0*cm, -CLOVERtoShield_displacement*cm);
    
    G4LogicalVolume * Logic_CLOVER_InternalVacuum[numberOf_CLOVER];
    G4LogicalVolume * Logic_CLOVER_Encasement;
    G4LogicalVolume * Logic_CLOVER_HPGeCrystal[4];
    
    if( CLOVER_Presence[0] || CLOVER_Presence[1] || CLOVER_Presence[2] || CLOVER_Presence[3] || CLOVER_Presence[4] || CLOVER_Presence[5] || CLOVER_Presence[6] || CLOVER_Presence[7] )
    {
        //////////////////////////////////////////////////////////
        //              CLOVER Internal Vacuum - CADMesh
        
        CADMesh * mesh_CLOVERInternalVacuum = new CADMesh("../AFRODITE/Mesh-Models/DETECTORS/CLOVER/CLOVER-InternalVacuum/CloverInternalVacuum.ply", "PLY", mm, offset_CLOVERInternalVacuum, false);
        
        G4VSolid * Solid_CLOVERInternalVacuum = mesh_CLOVERInternalVacuum->TessellatedMesh();
        
        for(G4int i=0; i<numberOf_CLOVER; i++)
        {
            Logic_CLOVER_InternalVacuum[i] = new G4LogicalVolume(Solid_CLOVERInternalVacuum, G4_Galactic_Material, "LogicCLOVERInternalVacuum", 0, 0, 0);
        }
        
        
        ///////////////////////////////////////////////////////
        //              CLOVER Encasement - CADMesh
        
        CADMesh * mesh_CLOVEREncasement = new CADMesh("../AFRODITE/Mesh-Models/DETECTORS/CLOVER/CloverEncasement/CloverEncasement.ply", "PLY", mm, offset_CLOVEREncasement, false);
        
        G4VSolid * Solid_CLOVEREncasement = mesh_CLOVEREncasement->TessellatedMesh();
        
        Logic_CLOVER_Encasement = new G4LogicalVolume(Solid_CLOVEREncasement, G4_Al_Material, "LogicCLOVERCloverEncasement", 0, 0, 0);
        
        
        //////////////////////////////////////////////////////////
        //              CLOVER HPGeCrystals - CADMesh
        
        CADMesh * mesh_CLOVERHPGeCrystal1 = new CADMesh("../AFRODITE/Mesh-Models/DETECTORS/CLOVER/HPGeCrystals/HPGeCrystal1.ply", "PLY", mm, offset_CLOVERHPGeCrystal1, false);
        CADMesh * mesh_CLOVERHPGeCrystal2 = new CADMesh("../AFRODITE/Mesh-Models/DETECTORS/CLOVER/HPGeCrystals/HPGeCrystal2.ply", "PLY", mm, offset_CLOVERHPGeCrystal2, false);
        CADMesh * mesh_CLOVERHPGeCrystal3 = new CADMesh("../AFRODITE/Mesh-Models/DETECTORS/CLOVER/HPGeCrystals/HPGeCrystal3.ply", "PLY", mm, offset_CLOVERHPGeCrystal3, false);
        CADMesh * mesh_CLOVERHPGeCrystal4 = new CADMesh("../AFRODITE/Mesh-Models/DETECTORS/CLOVER/HPGeCrystals/HPGeCrystal4.ply", "PLY", mm, offset_CLOVERHPGeCrystal4, false);
        
        G4VSolid * Solid_HPGeCrystal1 = mesh_CLOVERHPGeCrystal1->TessellatedMesh();
        G4VSolid * Solid_HPGeCrystal2 = mesh_CLOVERHPGeCrystal2->TessellatedMesh();
        G4VSolid * Solid_HPGeCrystal3 = mesh_CLOVERHPGeCrystal3->TessellatedMesh();
        G4VSolid * Solid_HPGeCrystal4 = mesh_CLOVERHPGeCrystal4->TessellatedMesh();
        
        Logic_CLOVER_HPGeCrystal[0] = new G4LogicalVolume(Solid_HPGeCrystal1, G4_Ge_Material,"LogicCLOVERHPGeCrystal",0,0,0);
        Logic_CLOVER_HPGeCrystal[1] = new G4LogicalVolume(Solid_HPGeCrystal2, G4_Ge_Material,"LogicCLOVERHPGeCrystal",0,0,0);
        Logic_CLOVER_HPGeCrystal[2] = new G4LogicalVolume(Solid_HPGeCrystal3, G4_Ge_Material,"LogicCLOVERHPGeCrystal",0,0,0);
        Logic_CLOVER_HPGeCrystal[3] = new G4LogicalVolume(Solid_HPGeCrystal4, G4_Ge_Material,"LogicCLOVERHPGeCrystal",0,0,0);
        
    }
    
    
    ////////////////////////////////////////////
    //              CLOVER SHIELD
    ////////////////////////////////////////////
    
    G4ThreeVector offset_CLOVER_Shield_Body = G4ThreeVector(0*cm, 0*cm, 0*cm);
    G4ThreeVector offset_CLOVER_Shield_Heavimet = G4ThreeVector(0*cm, 0*cm, 0*cm);
    G4ThreeVector offset_CLOVER_Shield_BGOCrystals = G4ThreeVector(0*cm, 0*cm, 0*cm);
    G4ThreeVector offset_CLOVER_Shield_PMT = G4ThreeVector(0*cm, 0*cm, 0*cm);
    G4ThreeVector offset_CLOVER_Shield_PMTConArray = G4ThreeVector(0*cm, 0*cm, 0*cm);
    
    G4LogicalVolume* Logic_CLOVER_Shield_Body;
    G4LogicalVolume* Logic_CLOVER_Shield_Heavimet;
    G4LogicalVolume* Logic_CLOVER_Shield_BGOCrystal[16];
    G4LogicalVolume* Logic_CLOVER_Shield_PMT[16];
    
    if( CLOVER_Shield_Presence[0] || CLOVER_Shield_Presence[1] || CLOVER_Shield_Presence[2] || CLOVER_Shield_Presence[3] || CLOVER_Shield_Presence[4] || CLOVER_Shield_Presence[5] || CLOVER_Shield_Presence[6] || CLOVER_Shield_Presence[7] || CLOVER_Shield_Presence[8])
    {
        
        ///////////////////////////////////////////////////////
        //              CLOVER Shield Body - CADMesh
        ///////////////////////////////////////////////////////
        
        CADMesh * mesh_CLOVER_Shield_Body = new CADMesh("../AFRODITE/Mesh-Models/DETECTORS/CLOVER/Shield/Body/Body.ply", "PLY", mm, offset_CLOVER_Shield_Body, false);
        
        G4VSolid * Solid_CLOVER_Shield_Body = mesh_CLOVER_Shield_Body->TessellatedMesh();
        
        Logic_CLOVER_Shield_Body = new G4LogicalVolume(Solid_CLOVER_Shield_Body, G4_Al_Material, "LogicCLOVERShieldBody", 0, 0, 0);


        ///////////////////////////////////////////////////////
        //              CLOVER Shield Heavimet - CADMesh
        ///////////////////////////////////////////////////////
        
        CADMesh * mesh_CLOVER_Shield_Heavimet = new CADMesh("../AFRODITE/Mesh-Models/DETECTORS/CLOVER/Shield/Heavimet/Heavimet.ply", "PLY", mm, offset_CLOVER_Shield_Heavimet, false);
        
        G4VSolid * Solid_CLOVER_Shield_Heavimet = mesh_CLOVER_Shield_Heavimet->TessellatedMesh();
        
        Logic_CLOVER_Shield_Heavimet = new G4LogicalVolume(Solid_CLOVER_Shield_Heavimet, Heavimet_Material, "LogicCLOVERShieldHeavimet", 0, 0, 0);

        
        
        ///////////////////////////////////////////////////////
        //      CLOVER Shield BGO Crystals - CADMesh
        ///////////////////////////////////////////////////////
        
        G4Box* Solid_BGOCrystal_template = new G4Box("BGOCrystal_template", (200./2)*mm, (200./2)*mm, (600./2)*mm);
        G4Box* Solid_BGOCrystal_razor = new G4Box("BGOCrystal_razor", (500./2)*mm, (500./2)*mm, (500./2)*mm);
        G4Tubs* Solid_BGOCrystal_PMTrazor = new G4Tubs("BGOCrystal_PMTrazor", 0.*mm, 500.*mm, 57.*mm, 0.*deg, 360.*deg);

        G4VSolid* Solid_BGOCrystal_shear[16][7];
        G4ThreeVector   position_BGOCrystal_razor[16][7];
        G4RotationMatrix* rm_BGOCrystal_razor[16][7];
        
        for(G4int i=0; i<16; i++)
        {
            for(G4int j=0; j<7; j++)
            {
                rm_BGOCrystal_razor[i][j] = new G4RotationMatrix();
            }
        }

        ////    The components of the positions used by the razor to shear each BGOCrystal shape
        ////    There are two "unique crystals" (up to rotations/reflections), therefore only two sets of positions are needed.
        ////    There are seven shears needed per BGO Crystal
        G4double    cP_BGOCrystal_razor[2][7][3];
        G4double    sign_BGOCrystal_razor[8][3];
        G4double    crm_BGOCrystal_razor[2][7];
        G4int       signrm_BGOCrystal_quad[4][3];

        ////    The components of positions for the razor
        ////    cP_BGOCrystal1_razor
        cP_BGOCrystal_razor[0][0][0] = 0., cP_BGOCrystal_razor[0][0][1] = 0., cP_BGOCrystal_razor[0][0][2] = (250.-42.);
        cP_BGOCrystal_razor[0][1][0] = -0.2000000, cP_BGOCrystal_razor[0][1][1] = -189.5684723, cP_BGOCrystal_razor[0][1][2] = -213.7959824;
        cP_BGOCrystal_razor[0][2][0] = -0.2000000, cP_BGOCrystal_razor[0][2][1] = 326.0592535, cP_BGOCrystal_razor[0][2][2] = -150.4848585;
        cP_BGOCrystal_razor[0][3][0] = -0.2000000, cP_BGOCrystal_razor[0][3][1] = 312.5280186, cP_BGOCrystal_razor[0][3][2] = -74.2400764;
        cP_BGOCrystal_razor[0][4][0] = (-0.20 + 250.), cP_BGOCrystal_razor[0][4][1] = 74.3269833, cP_BGOCrystal_razor[0][4][2] = -151.6673049;
        cP_BGOCrystal_razor[0][5][0] = -287.9700263, cP_BGOCrystal_razor[0][5][1] = 61.7954410, cP_BGOCrystal_razor[0][5][2] = -209.6135091;
        cP_BGOCrystal_razor[0][6][0] = 0., cP_BGOCrystal_razor[0][6][1] = 83.3486382, cP_BGOCrystal_razor[0][6][2] = -310.4623689;

        ////    cP_BGOCrystal2_razor
        cP_BGOCrystal_razor[1][0][0] = 0., cP_BGOCrystal_razor[1][0][1] = 0., cP_BGOCrystal_razor[1][0][2] = (250.-42.);
        cP_BGOCrystal_razor[1][1][0] = -0.2000000, cP_BGOCrystal_razor[1][1][1] = -189.5684723, cP_BGOCrystal_razor[1][1][2] = -213.7959824;
        cP_BGOCrystal_razor[1][2][0] = -0.2000000, cP_BGOCrystal_razor[1][2][1] = 326.0592535, cP_BGOCrystal_razor[1][2][2] = -150.4848585;
        cP_BGOCrystal_razor[1][3][0] = -0.2000000, cP_BGOCrystal_razor[1][3][1] = 312.5280186, cP_BGOCrystal_razor[1][3][2] = -74.2400764;
        cP_BGOCrystal_razor[1][4][0] = 208.0290511, cP_BGOCrystal_razor[1][4][1] = 69.0956904, cP_BGOCrystal_razor[1][4][2] = -269.0694506;
        cP_BGOCrystal_razor[1][5][0] = -252.0980447, cP_BGOCrystal_razor[1][5][1] = -101.2517474, cP_BGOCrystal_razor[1][5][2] = -172.2364983;
        cP_BGOCrystal_razor[1][6][0] = 0., cP_BGOCrystal_razor[1][6][1] = 83.3486382, cP_BGOCrystal_razor[1][6][2] = -310.4623689;

        ////    Reflections (for a total of 8 versions of the aforementioned two "unique" BGO Crystals
        sign_BGOCrystal_razor[0][0] = 1, sign_BGOCrystal_razor[0][1] = 1, sign_BGOCrystal_razor[0][2] = 1;
        sign_BGOCrystal_razor[1][0] = -1, sign_BGOCrystal_razor[1][1] = -1, sign_BGOCrystal_razor[1][2] = 1;
        sign_BGOCrystal_razor[2][0] = -1, sign_BGOCrystal_razor[2][1] = 1, sign_BGOCrystal_razor[2][2] = 1;
        sign_BGOCrystal_razor[3][0] = 1, sign_BGOCrystal_razor[3][1] = -1, sign_BGOCrystal_razor[3][2] = 1;
        sign_BGOCrystal_razor[4][0] = -1, sign_BGOCrystal_razor[4][1] = -1, sign_BGOCrystal_razor[4][2] = 1;
        sign_BGOCrystal_razor[5][0] = 1, sign_BGOCrystal_razor[5][1] = 1, sign_BGOCrystal_razor[5][2] = 1;
        sign_BGOCrystal_razor[6][0] = 1, sign_BGOCrystal_razor[6][1] = -1, sign_BGOCrystal_razor[6][2] = 1;
        sign_BGOCrystal_razor[7][0] = -1, sign_BGOCrystal_razor[7][1] = 1, sign_BGOCrystal_razor[7][2] = 1;

        
        ////    A sign per rotation (about 3 axes), per quadrant
        ////    0->X, 1->Y, 2->Z
        signrm_BGOCrystal_quad[0][0] = 1, signrm_BGOCrystal_quad[0][1] = 1, signrm_BGOCrystal_quad[0][2] = 1;
        signrm_BGOCrystal_quad[1][0] = -1, signrm_BGOCrystal_quad[1][1] = 1, signrm_BGOCrystal_quad[1][2] = 1;
        signrm_BGOCrystal_quad[2][0] = -1, signrm_BGOCrystal_quad[2][1] = -1, signrm_BGOCrystal_quad[2][2] = 1;
        signrm_BGOCrystal_quad[3][0] = 1, signrm_BGOCrystal_quad[3][1] = -1, signrm_BGOCrystal_quad[3][2] = 1;

        ////    The two "unique" sets of rotations;
        crm_BGOCrystal_razor[0][0] = 0.;
        crm_BGOCrystal_razor[0][1] = -7.;
        crm_BGOCrystal_razor[0][2] = -7.;
        crm_BGOCrystal_razor[0][3] = -14.5;
        crm_BGOCrystal_razor[0][4] = 0.;
        crm_BGOCrystal_razor[0][5] = -7.;
        crm_BGOCrystal_razor[0][6] = -3.;

        crm_BGOCrystal_razor[1][0] = 0.;
        crm_BGOCrystal_razor[1][1] = -7.;
        crm_BGOCrystal_razor[1][2] = -7.;
        crm_BGOCrystal_razor[1][3] = -14.5;
        crm_BGOCrystal_razor[1][4] = -7.;
        crm_BGOCrystal_razor[1][5] = -45.;
        crm_BGOCrystal_razor[1][6] = -3.;
        
        
        //  ph suffix-> place holder
        G4double    cP_BGOCrystal_ph[3];
        G4double    signP_BGOCrystal_ph[3];

        G4int       index1, index2;
        G4double    swap[2];
        
        for(G4int i=0; i<16; i++)
        {
            ////    Determining which copy number (0->7) of pair of "two unique" BGO Crystals we are working with
            if( ((i/2)%2)==0 )
            {
                ////    Depending upon the copy number, the order of the crystals shall change
                if((i%2)==0)
                    index1 = 0;
                else
                    index1 = 1;
            }
            else
            {
                ////    Depending upon the copy number, the order of the crystals shall change
                if((i%2)==0)
                    index1 = 1;
                else
                    index1 = 0;
            }

            ////    Positions of the shearing
            for(G4int j=0; j<7; j++)
            {
                cP_BGOCrystal_ph[0] = cP_BGOCrystal_razor[index1][j][0];
                cP_BGOCrystal_ph[1] = cP_BGOCrystal_razor[index1][j][1];
                cP_BGOCrystal_ph[2] = cP_BGOCrystal_razor[index1][j][2];
                
                ////    Reflection
                if((i>1 && i<6) || (i>9 && i<14))
                {
                    swap[0] = cP_BGOCrystal_ph[0];
                    swap[1] = cP_BGOCrystal_ph[1];

                    cP_BGOCrystal_ph[0] = swap[1];
                    cP_BGOCrystal_ph[1] = swap[0];
                }
                
                signP_BGOCrystal_ph[0] = sign_BGOCrystal_razor[i/2][0];
                signP_BGOCrystal_ph[1] = sign_BGOCrystal_razor[i/2][1];
                signP_BGOCrystal_ph[2] = sign_BGOCrystal_razor[i/2][2];

                position_BGOCrystal_razor[i][j] = G4ThreeVector(signP_BGOCrystal_ph[0]*cP_BGOCrystal_ph[0]*mm,
                                                                signP_BGOCrystal_ph[1]*cP_BGOCrystal_ph[1]*mm,
                                                                signP_BGOCrystal_ph[2]*cP_BGOCrystal_ph[2]*mm);

            }
            
            ////    Rotations of the shearing
            ////    Determining which quadrant (eg. (X-, Y+)) of BGO Crystals we are working with
            ////    This determines between the two "unique" planes of rotation for the shearing
            ////    The "index1" is also neccesary to determine which of the two crystals we are working with
            
            
            ////    Determining the quadrant
            index2 = i - ((i/4)*4);
            
            if((i/4)==0 || (i/4)==2) index2 = index2;
            if((i/4)==1 || (i/4)==3) index2 = 3 - index2;
 
            if(index2==0)
            {
                rm_BGOCrystal_razor[i][0]->rotateZ(signrm_BGOCrystal_quad[(i/4)][2]*crm_BGOCrystal_razor[0][0]*deg);
                rm_BGOCrystal_razor[i][1]->rotateX(signrm_BGOCrystal_quad[(i/4)][0]*crm_BGOCrystal_razor[0][1]*deg);
                rm_BGOCrystal_razor[i][2]->rotateX(signrm_BGOCrystal_quad[(i/4)][0]*crm_BGOCrystal_razor[0][2]*deg);
                rm_BGOCrystal_razor[i][3]->rotateX(signrm_BGOCrystal_quad[(i/4)][0]*crm_BGOCrystal_razor[0][3]*deg);
                rm_BGOCrystal_razor[i][4]->rotateX(signrm_BGOCrystal_quad[(i/4)][0]*crm_BGOCrystal_razor[0][4]*deg);
                rm_BGOCrystal_razor[i][5]->rotateY(signrm_BGOCrystal_quad[(i/4)][1]*crm_BGOCrystal_razor[0][5]*deg);
                rm_BGOCrystal_razor[i][6]->rotateX(signrm_BGOCrystal_quad[(i/4)][0]*crm_BGOCrystal_razor[0][6]*deg);
            }
            if(index2==1)
            {
                rm_BGOCrystal_razor[i][0]->rotateZ(signrm_BGOCrystal_quad[(i/4)][2]*crm_BGOCrystal_razor[1][0]*deg);
                rm_BGOCrystal_razor[i][1]->rotateX(signrm_BGOCrystal_quad[(i/4)][0]*crm_BGOCrystal_razor[1][1]*deg);
                rm_BGOCrystal_razor[i][2]->rotateX(signrm_BGOCrystal_quad[(i/4)][0]*crm_BGOCrystal_razor[1][2]*deg);
                rm_BGOCrystal_razor[i][3]->rotateX(signrm_BGOCrystal_quad[(i/4)][0]*crm_BGOCrystal_razor[1][3]*deg);
                rm_BGOCrystal_razor[i][4]->rotateY(signrm_BGOCrystal_quad[(i/4)][1]*crm_BGOCrystal_razor[1][4]*deg);
                rm_BGOCrystal_razor[i][5]->rotateZ(signrm_BGOCrystal_quad[(i/4)][2]*crm_BGOCrystal_razor[1][5]*deg);
                rm_BGOCrystal_razor[i][6]->rotateX(signrm_BGOCrystal_quad[(i/4)][0]*crm_BGOCrystal_razor[1][6]*deg);
            }
            if(index2==2)
            {
                rm_BGOCrystal_razor[i][0]->rotateZ(signrm_BGOCrystal_quad[(i/4)][2]*crm_BGOCrystal_razor[1][0]*deg);
                rm_BGOCrystal_razor[i][1]->rotateY(signrm_BGOCrystal_quad[(i/4)][1]*crm_BGOCrystal_razor[1][1]*deg);
                rm_BGOCrystal_razor[i][2]->rotateY(signrm_BGOCrystal_quad[(i/4)][1]*crm_BGOCrystal_razor[1][2]*deg);
                rm_BGOCrystal_razor[i][3]->rotateY(signrm_BGOCrystal_quad[(i/4)][1]*crm_BGOCrystal_razor[1][3]*deg);
                rm_BGOCrystal_razor[i][4]->rotateX(signrm_BGOCrystal_quad[(i/4)][0]*crm_BGOCrystal_razor[1][4]*deg);
                rm_BGOCrystal_razor[i][5]->rotateZ(signrm_BGOCrystal_quad[(i/4)][2]*crm_BGOCrystal_razor[1][5]*deg);
                rm_BGOCrystal_razor[i][6]->rotateY(signrm_BGOCrystal_quad[(i/4)][1]*crm_BGOCrystal_razor[1][6]*deg);
            }
            if(index2==3)
            {
                rm_BGOCrystal_razor[i][0]->rotateZ(signrm_BGOCrystal_quad[(i/4)][2]*crm_BGOCrystal_razor[0][0]*deg);
                rm_BGOCrystal_razor[i][1]->rotateY(signrm_BGOCrystal_quad[(i/4)][1]*crm_BGOCrystal_razor[0][1]*deg);
                rm_BGOCrystal_razor[i][2]->rotateY(signrm_BGOCrystal_quad[(i/4)][1]*crm_BGOCrystal_razor[0][2]*deg);
                rm_BGOCrystal_razor[i][3]->rotateY(signrm_BGOCrystal_quad[(i/4)][1]*crm_BGOCrystal_razor[0][3]*deg);
                rm_BGOCrystal_razor[i][4]->rotateY(signrm_BGOCrystal_quad[(i/4)][1]*crm_BGOCrystal_razor[0][4]*deg);
                rm_BGOCrystal_razor[i][5]->rotateX(signrm_BGOCrystal_quad[(i/4)][0]*crm_BGOCrystal_razor[0][5]*deg);
                rm_BGOCrystal_razor[i][6]->rotateY(signrm_BGOCrystal_quad[(i/4)][1]*crm_BGOCrystal_razor[0][6]*deg);
            }
            
            
        }

        
        char nameChar[512];
        G4String nameG4String;

        for(G4int i=0; i<16; i++)
        {
            for(G4int j=0; j<6; j++)
            {
                sprintf(nameChar,"Solid_BGOCrystal%d_shear%d", i, j);
                nameG4String = string(nameChar);
                
                if(j==0)
                {
                    Solid_BGOCrystal_shear[i][j] = new G4SubtractionSolid(nameG4String, Solid_BGOCrystal_template, Solid_BGOCrystal_razor, rm_BGOCrystal_razor[i][j], position_BGOCrystal_razor[i][j]);
                }
                else
                {
                    Solid_BGOCrystal_shear[i][j] = new G4SubtractionSolid(nameG4String, Solid_BGOCrystal_shear[i][j-1], Solid_BGOCrystal_razor, rm_BGOCrystal_razor[i][j], position_BGOCrystal_razor[i][j]);
                }
                
            }
            
            ////    The final shear for each crystal (corresponding to the surface that the PMT attaches) is executed in such a way that the resulting surface is well defined for simple PMT coupling
            ////    This scintillation still needs to be fully implemented
            
            sprintf(nameChar,"Solid_BGOCrystal%d", i);
            nameG4String = string(nameChar);
            
            Solid_BGOCrystal_shear[i][6] = new G4SubtractionSolid(nameChar, Solid_BGOCrystal_shear[i][5], Solid_BGOCrystal_PMTrazor, rm_BGOCrystal_razor[i][6], position_BGOCrystal_razor[i][6]);

        }
        
        for(G4int i=0; i<16; i++)
        {
            Logic_CLOVER_Shield_BGOCrystal[i] = new G4LogicalVolume(Solid_BGOCrystal_shear[i][6], G4_BGO_Material, "LogicCLOVERShieldBGOCrystal", 0, 0, 0);
        }
        
        
        
        ////////////////////////////////////
        //      CLOVER Shield PMT's
        ////////////////////////////////////
        
        G4VSolid* Solid_CLOVER_Shield_PMT[16];
        
        ////    Positioning of the centre of the PMT's
        ////    Once again, there are only two "unique" positions (the rest can be obtained through reflections about the cartesian axes)
        G4ThreeVector   position_Shield_PMT[16];
        G4double        cP_Shield_PMT[2][3];
        
        cP_Shield_PMT[0][0] = -22., cP_Shield_PMT[0][1] = 83.3486382, cP_Shield_PMT[0][2] = -310.4623689;
        cP_Shield_PMT[1][0] = -58.5, cP_Shield_PMT[1][1] = 83.3486382, cP_Shield_PMT[1][2] = -310.4623689;

        
        ////    A unique rotation of the PMT's within each quadrant (X+, X-, Y+, Y-)
        G4RotationMatrix* rm_Shield_PMT[4];
        
        for(G4int i=0; i<4; i++)
        {
            rm_Shield_PMT[i] = new G4RotationMatrix();
        }
        
        rm_Shield_PMT[0]->rotateX(-3.*deg);
        rm_Shield_PMT[1]->rotateY(-3.*deg);
        rm_Shield_PMT[2]->rotateX(3.*deg);
        rm_Shield_PMT[3]->rotateY(3.*deg);
        
        G4Tubs* Solid_Shield_PMTtemplate = new G4Tubs("Shield_PMT", 0.*mm, 14.3*mm, 57.*mm, 0.*deg, 360.*deg);
        G4Box* Solid_Shield_PMTRaw = new G4Box("Shield_PMTtemplate", (1000./2)*mm, (1000./2)*mm, (1000./2)*mm);

        
        //  placeholders
        G4double    cP_Shield_PMT_ph[3];
        G4double    signP_Shield_PMT_ph[3];

        for(G4int i=0; i<16; i++)
        {
            ////    Determining which copy number (0->7) of pair of "two unique" PMT positions we are working with
            if( ((i/2)%2)==0 )
            {
                ////    Depending upon the copy number, the order of the PMT's shall change
                if((i%2)==0)
                    index1 = 0;
                else
                    index1 = 1;
            }
            else
            {
                ////    Depending upon the copy number, the order of the PMT's shall change
                if((i%2)==0)
                    index1 = 1;
                else
                    index1 = 0;
            }
            
            ////    Positions of the Intersection
            cP_Shield_PMT_ph[0] = cP_Shield_PMT[index1][0];
            cP_Shield_PMT_ph[1] = cP_Shield_PMT[index1][1];
            cP_Shield_PMT_ph[2] = cP_Shield_PMT[index1][2];
            
            if((i>1 && i<6) || (i>9 && i<14))
            {
                swap[0] = cP_Shield_PMT_ph[0];
                swap[1] = cP_Shield_PMT_ph[1];
                
                cP_Shield_PMT_ph[0] = swap[1];
                cP_Shield_PMT_ph[1] = swap[0];
            }
            
            signP_Shield_PMT_ph[0] = sign_BGOCrystal_razor[i/2][0];
            signP_Shield_PMT_ph[1] = sign_BGOCrystal_razor[i/2][1];
            signP_Shield_PMT_ph[2] = sign_BGOCrystal_razor[i/2][2];
            
            position_Shield_PMT[i] = G4ThreeVector(signP_Shield_PMT_ph[0]*cP_Shield_PMT_ph[0]*mm,
                                                   signP_Shield_PMT_ph[1]*cP_Shield_PMT_ph[1]*mm,
                                                   signP_Shield_PMT_ph[2]*cP_Shield_PMT_ph[2]*mm);
            
            ////////////////////////////////////////////////////////////////////
            
            sprintf(nameChar,"Solid_CLOVER_Shield_PMT%d", i);
            nameG4String = string(nameChar);
            
            ////    Selecting which quadrant of PMT's (X+, X-, Y+, Y-) that we are working with
            if(i==0 || i==1 || i==14 || i==15)  index2 = 0;
            if(i>=2 && i<=5)  index2 = 1;
            if(i>=6 && i<=9)  index2 = 2;
            if(i>=10 && i<=13)  index2 = 3;

            Solid_CLOVER_Shield_PMT[i] = new G4IntersectionSolid(nameG4String, Solid_Shield_PMTRaw, Solid_Shield_PMTtemplate, rm_Shield_PMT[index2], position_Shield_PMT[i]);
            
            Logic_CLOVER_Shield_PMT[i] = new G4LogicalVolume(Solid_CLOVER_Shield_PMT[i], G4_AIR_Material,"LogicCLOVERShieldPMT",0,0,0);

        }
    }
    
    
    ////////////////////////////////////////////////////
    //               CLOVER INITIALIZATION            //
    ////////////////////////////////////////////////////
    
    
    for(G4int i=0; i<numberOf_CLOVER; i++)
    {
        CLOVER_position[i] = (CLOVER_Distance[i])*G4ThreeVector( sin(CLOVER_theta[i]) * cos(CLOVER_phi[i]), sin(CLOVER_theta[i]) * sin(CLOVER_phi[i]), cos(CLOVER_theta[i]));
        
        CLOVER_transform[i] = G4Transform3D(CLOVER_rotm[i],CLOVER_position[i]);
        
        CLOVER_Shield_position[i] = (CLOVER_Distance[i])*G4ThreeVector( sin(CLOVER_theta[i]) * cos(CLOVER_phi[i]), sin(CLOVER_theta[i]) * sin(CLOVER_phi[i]), cos(CLOVER_theta[i]));
        
        CLOVER_Shield_transform[i] = CLOVER_transform[i];
        
        /////////////////////////////
        //          CLOVER
        if(CLOVER_Presence[i])
        {
            new G4PVPlacement(CLOVER_transform[i],   // transformation matrix
                              Logic_CLOVER_Encasement,       // its logical volume
                              "CLOVER_Encasement",       // its name
                              LogicVacuumChamber,         // its mother  volume
                              false,           // no boolean operations
                              i,               // copy number
                              fCheckOverlaps); // checking overlaps
            
            
            
            
            
            for (int j=0; j<4; j++)
            {
                PhysiCLOVER_HPGeCrystal = new G4PVPlacement(0,               // no rotation
                                                            G4ThreeVector(0,0,0), // at (x,y,z)
                                                            Logic_CLOVER_HPGeCrystal[j],
                                                            "CLOVER_HPGeCrystal", // its name
                                                            Logic_CLOVER_InternalVacuum[i],
                                                            false,           // no boolean operations
                                                            i*4 + j,               // copy number
                                                            fCheckOverlaps); // checking overlaps
                
            }
            
            
            new G4PVPlacement(CLOVER_transform[i],
                              Logic_CLOVER_InternalVacuum[i],
                              "CLOVER_InternalVacuum",       // its name
                              LogicVacuumChamber,         // its mother  volume
                              false,           // no boolean operations
                              i,               // copy number
                              fCheckOverlaps); // checking overlaps
            
        }
        
        
        /////////////////////////////
        //      CLOVER SHIELD
        if(CLOVER_Shield_Presence[i])
        {
            
            new G4PVPlacement(CLOVER_Shield_transform[i],
                              Logic_CLOVER_Shield_Body,
                              "CLOVER_Shield_Body",       // its name
                              LogicVacuumChamber,         // its mother  volume
                              false,           // no boolean operations
                              i,               // copy number
                              fCheckOverlaps); // checking overlaps
            
            for (int j=0; j<16; j++)
            {
                PhysiCLOVER_Shield_BGOCrystal = new G4PVPlacement(CLOVER_Shield_transform[i],
                                                                  Logic_CLOVER_Shield_BGOCrystal[j],
                                                                  "CLOVER_Shield_BGOCrystal",
                                                                  LogicVacuumChamber,
                                                                  false,
                                                                  i*16 + j,
                                                                  fCheckOverlaps);
                
                PhysiCLOVER_Shield_PMT = new G4PVPlacement(CLOVER_Shield_transform[i],
                                                           Logic_CLOVER_Shield_PMT[j],
                                                           "CLOVER_Shield_PMT", // its name
                                                           LogicVacuumChamber,
                                                           false, // no boolean operations
                                                           i*16 + j,  // copy number
                                                           fCheckOverlaps); // checking overlaps
                
            }
            

            new G4PVPlacement(CLOVER_Shield_transform[i],
                              Logic_CLOVER_Shield_Heavimet,
                              "CLOVER_Shield_Heavimet",       // its name
                              LogicVacuumChamber,         // its mother  volume
                              false,           // no boolean operations
                              i,               // copy number
                              fCheckOverlaps); // checking overlaps
        
        }
    }
    

    //////////////////////////////////////
    //      PlasticScint DEFINITION     //
    //////////////////////////////////////

    G4LogicalVolume* Logic_PlasticScint[numberOf_PlasticScint];
    G4Box* Solid_PlasticScint = new G4Box("Scintillator", (600/2)*mm, (100/2)*mm, (100/2)*mm);
    
    for(G4int i=0; i<numberOf_PlasticScint; i++)
    {
        Logic_PlasticScint[i] = new G4LogicalVolume(Solid_PlasticScint, BC408_Material,"Scintillator",0,0,0);
    }

    
    //////////////////////////////////////////
    //      PlasticScint INITIALISATION     //
    //////////////////////////////////////////
    
    for(G4int i=0; i<numberOf_PlasticScint; i++)
    {
        if(PlasticScint_Presence[i])
        {
            PlasticScint_transform[i] = G4Transform3D(PlasticScint_rotm[i],PlasticScint_CentrePosition[i]);
            
            PhysiPlasticScint = new G4PVPlacement(PlasticScint_transform[i],
                                            Logic_PlasticScint[i],    // its logical volume
                                            "PlasticScint",           // its name
                                            LogicWorld,         // its mother  volume
                                            false,              // no boolean operations
                                            i,                  // copy number
                                            fCheckOverlaps);    // checking overlaps
            
        }
    }
    

    
    
    ////////////////////////////////////////
    ////        LEPS DEFINITION         ////
    ////////////////////////////////////////
    
    //////////////////////////////////////////////////////////
    //              LEPS Internal Vacuum - CADMesh
    //////////////////////////////////////////////////////////
    
    G4Tubs* Solid_LEPS_InternalVacuum = new G4Tubs("Solid_LEPSInternalVacuum", 0.*mm, 38.4*mm, 45.0*mm, 0.*deg, 360*deg);
    G4LogicalVolume* Logic_LEPS_InternalVacuum[numberOf_LEPS];
    
    for(G4int i=0; i<numberOf_LEPS; i++)
    {
        Logic_LEPS_InternalVacuum[i] = new G4LogicalVolume(Solid_LEPS_InternalVacuum, G4_Galactic_Material, "LogicLEPSInternalVacuum", 0, 0, 0);
    }
    
    ///////////////////////////////////////////////////////
    //              LEPS Encasement - CADMesh
    ///////////////////////////////////////////////////////
    
    G4Tubs* Solid_LEPS_Encasement = new G4Tubs("Solid_LEPSEncasement", 38.5*mm, 40.0*mm, (90./2)*mm, 0.*deg, 360*deg);
    
    G4LogicalVolume* Logic_LEPS_Encasement = new G4LogicalVolume(Solid_LEPS_Encasement, G4_Al_Material, "LogicLEPSLEPSEncasement", 0, 0, 0);
    
    
    ///////////////////////////////////////////////////////
    //              LEPS Beryllium Window - CADMesh
    ///////////////////////////////////////////////////////
    
    G4Tubs* Solid_LEPS_Window = new G4Tubs("Solid_LEPSWindow", 0.*mm, 38.5*mm, (0.3/2)*mm, 0.*deg, 360*deg);
    
    G4LogicalVolume* Logic_LEPS_Window = new G4LogicalVolume(Solid_LEPS_Window, G4_Be_Material, "Logic_LEPS_Window", 0, 0, 0);
    
    
    //////////////////////////////////////////////////////////
    //              LEPS HPGeCrystals - CADMesh
    //////////////////////////////////////////////////////////
    
    G4Tubs* Solid_HPGeCrystal = new G4Tubs("Solid_HPGeCrystal1", 0.*mm, 33.0*mm, 5.5*mm, 0.*deg, 90.*deg);
    
    G4LogicalVolume* Logic_LEPS_HPGeCrystal;
    Logic_LEPS_HPGeCrystal = new G4LogicalVolume(Solid_HPGeCrystal, G4_Ge_Material,"LogicLEPSHPGeCrystal",0,0,0);
    
    LEPS_HPGeCrystal_rotm[0].rotateZ(0.*deg);
    LEPS_HPGeCrystal_rotm[1].rotateZ(90.*deg);
    LEPS_HPGeCrystal_rotm[2].rotateZ(180.*deg);
    LEPS_HPGeCrystal_rotm[3].rotateZ(270.*deg);

    for(G4int i=0; i<4; i++)
    {
        LEPS_HPGeCrystal_transform[i] = G4Transform3D(LEPS_HPGeCrystal_rotm[i], G4ThreeVector(0,0,(29.0-0.5)*mm));
    }
    
    ////////////////////////////////////////////////////
    //               LEPS INITIALIZATION
    ////////////////////////////////////////////////////
    
    
    for(G4int i=0; i<numberOf_LEPS; i++)
    {
        LEPS_position[i] = (LEPS_Distance[i] + 4.5*cm)*G4ThreeVector( std::sin(LEPS_theta[i]) * std::cos(LEPS_phi[i]), std::sin(LEPS_theta[i]) * std::sin(LEPS_phi[i]), std::cos(LEPS_theta[i]));
        
        LEPS_transform[i] = G4Transform3D(LEPS_rotm[i],LEPS_position[i]);
        
        LEPS_InternalVacuum_position[i] = (LEPS_Distance[i]+ 4.5*cm + 0.5*mm)*G4ThreeVector( std::sin(LEPS_theta[i]) * std::cos(LEPS_phi[i]), std::sin(LEPS_theta[i]) * std::sin(LEPS_phi[i]), std::cos(LEPS_theta[i]));
        LEPS_InternalVacuum_transform[i] = G4Transform3D(LEPS_rotm[i],LEPS_InternalVacuum_position[i]);
        
        LEPS_Window_position[i] = (LEPS_Distance[i] + 4.5*cm +(-45.0+0.15)*mm)*G4ThreeVector( std::sin(LEPS_theta[i]) * std::cos(LEPS_phi[i]), std::sin(LEPS_theta[i]) * std::sin(LEPS_phi[i]), std::cos(LEPS_theta[i]));
        LEPS_Window_transform[i] = G4Transform3D(LEPS_rotm[i],LEPS_Window_position[i]);
        
        /////////////////////////////
        //          LEPS
        if(LEPS_Presence[i] == true)
        {
            
            new G4PVPlacement(LEPS_transform[i],   // transformation matrix
                              Logic_LEPS_Encasement,       // its logical volume
                              "LEPSEncasement",       // its name
                              LogicVacuumChamber,         // its mother  volume
                              false,           // no boolean operations
                              i,               // copy number
                              fCheckOverlaps); // checking overlaps
            
            new G4PVPlacement(LEPS_Window_transform[i],   // transformation matrix
                              Logic_LEPS_Window,       // its logical volume
                              "LEPSWindow",       // its name
                              LogicVacuumChamber,         // its mother  volume
                              false,           // no boolean operations
                              i,               // copy number
                              fCheckOverlaps); // checking overlaps
            
            
            
            for (int j=0; j<4; j++)
            {
                Physical_LEPS_HPGeCrystal = new G4PVPlacement(LEPS_HPGeCrystal_transform[j],
                                                              Logic_LEPS_HPGeCrystal,       // its logical volume
                                                              "LEPSHPGeCrystal",       // its name
                                                              Logic_LEPS_InternalVacuum[i],    // its mother  volume
                                                              false,           // no boolean operations
                                                              j + (i*4),               // copy number
                                                              fCheckOverlaps); // checking overlaps
                
            }
            
            new G4PVPlacement(LEPS_InternalVacuum_transform[i],
                              Logic_LEPS_InternalVacuum[i],
                              "LEPSInternalVacuum",       // its name
                              LogicVacuumChamber,         // its mother  volume
                              false,           // no boolean operations
                              i,               // copy number
                              fCheckOverlaps); // checking overlaps
            
            
        }
        
    }
    
    
    
    //////////////////////////////////
    //      HAGAR DEFINITION        //
    //////////////////////////////////
    
    ////////////////////////////////
    //      HAGAR - NaI Crystal
    G4Tubs* Solid_HAGAR_NaICrystal = new G4Tubs("HAGAR_NaICrystal", 0.*cm, (23.8/2)*cm, (35.6/2)*cm, 0.*deg, 360.*deg);
    
    G4LogicalVolume* Logic_HAGAR_NaICrystal = new G4LogicalVolume(Solid_HAGAR_NaICrystal, G4_SODIUM_IODIDE_Material, "Logic_HAGAR_NaICrystal", 0, 0, 0);
    
    ////////////////////////////////
    //      HAGAR - Annulus
    G4Tubs* Solid_HAGAR_Anulus = new G4Tubs("HAGAR_Anulus", (28.86/2)*cm, ((28.86/2) + 9.84)*cm, (61/2)*cm, 0.*deg, 360.*deg);
    
    G4LogicalVolume* Logic_HAGAR_Annulus = new G4LogicalVolume(Solid_HAGAR_Anulus, BC408_Material, "Logic_HAGAR_Annulus", 0, 0, 0);
    
    ////////////////////////////////
    //      HAGAR - Front Disc
    G4Tubs* Solid_HAGAR_FrontDisc = new G4Tubs("HAGAR_FrontDisc", 0.*cm, (48.58/2)*cm, (8/2)*cm, 0.*deg, 360.*deg);
    
    G4LogicalVolume* Logic_HAGAR_FrontDisc = new G4LogicalVolume(Solid_HAGAR_FrontDisc, BC408_Material, "Logic_HAGAR_FrontDisc", 0, 0, 0);
    
    
    
    ////////////////////////////////////////////////////
    //               HAGAR INITIALIZATION             //
    ////////////////////////////////////////////////////
    
    
    ///////////////////////////////
    //      HAGAR - NaI Crystal
    if(HAGAR_NaICrystal_Presence)
    {
        HAGAR_transform = G4Transform3D(HAGAR_rotm, HAGAR_NaICrystal_CentrePosition);
        
        PhysiHAGAR_NaICrystal = new G4PVPlacement(HAGAR_transform,   // transformation matrix
                                                  Logic_HAGAR_NaICrystal,       // its logical volume
                                                  "HAGAR_NaICrystal",       // its name
                                                  LogicWorld,         // its mother  volume
                                                  false,           // no boolean operations
                                                  0,               // copy number
                                                  fCheckOverlaps); // checking overlaps
    }
    
    /////////////////////////////
    //      HAGAR - Annulus
    if(HAGAR_Annulus_Presence)
    {
        HAGAR_transform = G4Transform3D(HAGAR_rotm, HAGAR_Annulus_CentrePosition);
        
        PhysiHAGAR_Annulus = new G4PVPlacement(HAGAR_transform,   // transformation matrix
                                               Logic_HAGAR_Annulus,       // its logical volume
                                               "HAGAR_Annulus",       // its name
                                               LogicWorld,         // its mother  volume
                                               false,           // no boolean operations
                                               0,               // copy number
                                               fCheckOverlaps); // checking overlaps
    }
    
    ///////////////////////////////
    //      HAGAR - Front Disc
    if(HAGAR_FrontDisc_Presence)
    {
        HAGAR_transform = G4Transform3D(HAGAR_rotm, HAGAR_FrontDisc_CentrePosition);
        
        PhysiHAGAR_FrontDisc = new G4PVPlacement(HAGAR_transform,   // transformation matrix
                                                 Logic_HAGAR_FrontDisc,       // its logical volume
                                                 "HAGAR_FrontDisc",       // its name
                                                 LogicWorld,         // its mother  volume
                                                 false,           // no boolean operations
                                                 0,               // copy number
                                                 fCheckOverlaps); // checking overlaps
    }

    
    
    
    

    
    //////////////////////////////////////////////////
    //      AFRODITE SPECTROMETER INITIIALIZATION       //
    //////////////////////////////////////////////////
    
    
    // A field object is held by a field manager
    // Find the global Field Manager
    G4TransportationManager* tmanagerMagneticField = G4TransportationManager::GetTransportationManager();
    tmanagerMagneticField->GetPropagatorInField()->SetLargestAcceptableStep(1*mm);
    G4double minStepMagneticField = 0.0025*mm ;
    
    
    //////////////////////////////////////////////////////
    //              AFRODITE - QUADRUPOLE
    //////////////////////////////////////////////////////
    
    if(AFRODITE_Quadrupole)
    {
        
        G4RotationMatrix* AFRODITE_Q_MagField_rotm = new G4RotationMatrix;
        //AFRODITE_Q_MagField_rotm->rotateX(90.*deg);
        
        AFRODITE_Quadrupole_transform = G4Transform3D(AFRODITE_Quadrupole_rotm, AFRODITE_Quadrupole_CentrePosition);
        
        G4Box* Solid_AFRODITE_Quadrupole = new G4Box("Solid_AFRODITE_Quadrupole", (50./2)*cm, (50./2)*cm, (30./2)*cm);
        
        G4LogicalVolume* Logic_AFRODITE_Quadrupole = new G4LogicalVolume(Solid_AFRODITE_Quadrupole, G4_Galactic_Material,"Logic_AFRODITE_Quadrupole",0,0,0);
        
        ////    IDEAL MAGNETIC FIELD for QUADRUPOLE
        if(Ideal_Quadrupole)
        {
            MagneticField_AFRODITE_Q = new G4QuadrupoleMagField(dBdr_gradient_AFRODITE_Q, AFRODITE_Quadrupole_CentrePosition, AFRODITE_Q_MagField_rotm);
            fEquationMagneticField_AFRODITE_Q = new G4Mag_UsualEqRhs(MagneticField_AFRODITE_Q);
            
            fieldManagerMagneticField_AFRODITE_Q = new G4FieldManager(MagneticField_AFRODITE_Q);
            
            stepperMagneticField_AFRODITE_Q = new G4ClassicalRK4( fEquationMagneticField_AFRODITE_Q );
            fieldManagerMagneticField_AFRODITE_Q -> SetDetectorField(MagneticField_AFRODITE_Q);
            
            fChordFinder_AFRODITE_Q = new G4ChordFinder( MagneticField_AFRODITE_Q, minStepMagneticField, stepperMagneticField_AFRODITE_Q);
            
            Logic_AFRODITE_Quadrupole -> SetFieldManager(fieldManagerMagneticField_AFRODITE_Q, true) ;
            
            
            //G4BlineTracer* theBlineTool = new G4BlineTracer();
            //theBlineTool->ComputeBlines();
        }
        
        ////    MAPPED MAGNETIC FIELD for QUADRUPOLE
        if(Mapped_Quadrupole)
        {
            
            G4double z_Q_Offset = 4.4*mm+ 100*cm;
            
            G4MagneticField* PurgMagField = new MagneticFieldMapping("../AFRODITE/MagneticFieldMaps/Quadrupole_MagneticFieldMap.TABLE", z_Q_Offset);
            fEquationMagneticField_AFRODITE_Q = new G4Mag_UsualEqRhs(PurgMagField);
            
            fieldManagerMagneticField_AFRODITE_Q = new G4FieldManager(PurgMagField);
            
            stepperMagneticField_AFRODITE_Q = new G4ClassicalRK4( fEquationMagneticField_AFRODITE_Q );
            fieldManagerMagneticField_AFRODITE_Q -> SetDetectorField(PurgMagField);
            
            fChordFinder_AFRODITE_Q = new G4ChordFinder( PurgMagField, minStepMagneticField, stepperMagneticField_AFRODITE_Q);
            
            Logic_AFRODITE_Quadrupole -> SetFieldManager(fieldManagerMagneticField_AFRODITE_Q, true) ;
            
            
            //G4BlineTracer* theBlineTool = new G4BlineTracer();
            
        }
        
        PhysiAFRODITE_Quadrupole = new G4PVPlacement(AFRODITE_Quadrupole_transform,
                                                 Logic_AFRODITE_Quadrupole,       // its logical volume
                                                 "AFRODITE_Quadrupole",       // its name
                                                 LogicWorld,         // its mother  volume
                                                 false,           // no boolean operations
                                                 0,               // copy number
                                                 fCheckOverlaps); // checking overlaps
        
    }
    
    //////////////////////////////////////////////////////
    //              AFRODITE - DIPOLE 1
    //////////////////////////////////////////////////////
    
    /*
    // magnetic field ----------------------------------------------------------
    MagneticField_AFRODITE_D1 = new G4UniformMagField(G4ThreeVector(0., AFRODITE_Dipole1_BZ, 0.));
    fieldManagerMagneticField_AFRODITE_D1 = new G4FieldManager();
    fieldManagerMagneticField_AFRODITE_D1->SetDetectorField(MagneticField_AFRODITE_D1);
    fieldManagerMagneticField_AFRODITE_D1->CreateChordFinder(MagneticField_AFRODITE_D1);
    G4bool forceToAllDaughters = true;
    */
    
    if(AFRODITE_Dipole1)
    {
        ////    MAGNETIC FIELD for DIPOLE 1
        MagneticField_AFRODITE_D1 = new G4UniformMagField(G4ThreeVector(0., AFRODITE_Dipole1_BZ, 0.));
        fEquationMagneticField_AFRODITE_D1 = new G4Mag_UsualEqRhs(MagneticField_AFRODITE_D1);
        
        fieldManagerMagneticField_AFRODITE_D1 = new G4FieldManager(MagneticField_AFRODITE_D1);
        
        stepperMagneticField_AFRODITE_D1 = new G4ClassicalRK4( fEquationMagneticField_AFRODITE_D1 );
        fieldManagerMagneticField_AFRODITE_D1 -> SetDetectorField(MagneticField_AFRODITE_D1);
        
        fChordFinder_AFRODITE_D1 = new G4ChordFinder( MagneticField_AFRODITE_D1, minStepMagneticField, stepperMagneticField_AFRODITE_D1);

        
        /////////////////////////////////////////////
        
        AFRODITE_Dipole1_transform = G4Transform3D(AFRODITE_Dipole1_rotm, AFRODITE_Dipole1_CentrePosition);
        
        //G4Box* Solid_AFRODITE_Dipole1 = new G4Box("Solid_AFRODITE_Dipole1", (50./2)*cm, (50./2)*cm, (30./2)*cm);
        //G4Tubs* Solid_AFRODITE_Dipole1 = new G4Tubs("Solid_AFRODITE_Dipole1", 50.*cm, 100.0*cm, 30.*cm, 0.*deg, 40.*deg);
        G4Tubs* Solid_AFRODITE_Dipole1 = new G4Tubs("Solid_AFRODITE_Dipole1", 30.*cm, 150.0*cm, 30.*cm, 0.*deg, 40.*deg);

        G4LogicalVolume* Logic_AFRODITE_Dipole1 = new G4LogicalVolume(Solid_AFRODITE_Dipole1, G4_Galactic_Material,"Logic_AFRODITE_Dipole1",0,0,0);
        Logic_AFRODITE_Dipole1 -> SetFieldManager(fieldManagerMagneticField_AFRODITE_D1, true) ;
        
        // Register the field and its manager for deleting
        //G4AutoDelete::Register(MagneticField_AFRODITE_D1);
        //G4AutoDelete::Register(fieldManagerMagneticField_AFRODITE_D1);
        
        PhysiAFRODITE_Dipole1 = new G4PVPlacement(AFRODITE_Dipole1_transform,
                                              Logic_AFRODITE_Dipole1,       // its logical volume
                                              "AFRODITE_Dipole1",       // its name
                                              LogicWorld,         // its mother  volume
                                              false,           // no boolean operations
                                              0,               // copy number
                                              fCheckOverlaps); // checking overlaps
        
    }
    
    
    

    
    
    
    //////////////////////////////////////////////////////////
    //                      VISUALISATION
    //////////////////////////////////////////////////////////
    
    G4VisAttributes* World_VisAtt= new G4VisAttributes(G4Colour(0., 0., 0.));
    World_VisAtt->SetVisibility(true);
    LogicWorld->SetVisAttributes(World_VisAtt);
    
    //G4VisAttributes* Target_VisAtt= new G4VisAttributes(G4Colour(1.0, 1.0, 1.0));
    //LogicTarget->SetVisAttributes(Target_VisAtt);
    
    
    ////////////////////////////
    //  Vacuum Chamber
    
    G4VisAttributes* VacuumChamber_VisAtt= new G4VisAttributes(G4Colour(1.0, 1.0, 1.0));
    //VacuumChamber_VisAtt->SetVisibility(false);
    LogicVacuumChamber->SetVisAttributes(VacuumChamber_VisAtt);
    
    
    //////////////////////////////////
    //      PlasticScint VISUALIZATION
    //////////////////////////////////
    
    G4VisAttributes* scintillatorVisAtt = new G4VisAttributes(G4Colour(1.0, 1.0, 0.0));
    for(G4int i=0; i<numberOf_PlasticScint; i++)
    {Logic_PlasticScint[i] -> SetVisAttributes(scintillatorVisAtt);}
    
    
    //////////////////////////////////
    //      HAGAR VISUALIZATION
    //////////////////////////////////
    
    //  HAGAR - NaI Crystal
    G4VisAttributes* HAGAR_NaICrystal_VisAtt = new G4VisAttributes(G4Colour(1.0, 1.0, 0.0));
    HAGAR_NaICrystal_VisAtt->SetForceSolid(true);
    Logic_HAGAR_NaICrystal->SetVisAttributes(HAGAR_NaICrystal_VisAtt);
    
    //  HAGAR - Annulus
    G4VisAttributes* HAGAR_Annulus_VisAtt = new G4VisAttributes(G4Colour(0.0, 0.7, 0.0));
    Logic_HAGAR_Annulus->SetVisAttributes(HAGAR_Annulus_VisAtt);
    
    //  HAGAR - Front Disc
    G4VisAttributes* HAGAR_FrontDisc_VisAtt = new G4VisAttributes(G4Colour(1.0, 1.0, 0.0));
    HAGAR_FrontDisc_VisAtt->SetForceSolid(true);
    Logic_HAGAR_FrontDisc->SetVisAttributes(HAGAR_FrontDisc_VisAtt);
    
    
    
    //
    //always return the physical World
    //
    return PhysiWorld;
    
    
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::ConstructField()
{ 
    // Create global magnetic field messenger.
    // Uniform magnetic field is then created automatically if
    // the field value is not zero.
    G4ThreeVector fieldValue = G4ThreeVector();
    //G4ThreeVector fieldValue = G4ThreeVector(0., 5*tesla, 0.);

    fMagFieldMessenger = new G4GlobalMagFieldMessenger(fieldValue);
    fMagFieldMessenger->SetVerboseLevel(1);
  
    // Register the field messenger for deleting
    G4AutoDelete::Register(fMagFieldMessenger);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
