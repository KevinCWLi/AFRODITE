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

#ifndef DetectorConstruction_h
#define DetectorConstruction_h 1

#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"
#include "G4RotationMatrix.hh"
#include "G4Transform3D.hh"

#include "G4UniformMagField.hh"
#include "G4QuadrupoleMagField.hh"
#include "G4Mag_UsualEqRhs.hh"
#include "G4MagIntegratorStepper.hh"
#include "G4ChordFinder.hh"

//#include "TabulatedField3D.hh"
#include "G4PropagatorInField.hh"
#include "G4PropagatorInField.hh"
#include "G4FieldManager.hh"


class G4VPhysicalVolume;
class G4GlobalMagFieldMessenger;

///////////////////////////////
class G4LogicalVolume;
class G4Material;
class G4UserLimits;
class G4Isotope;
class G4Element;
class G4LogicalBorderSurface;
class G4LogicalSkinSurface;
class G4Mag_UsualEqRhs;
class G4MagIntegratorStepper;
class G4ChordFinder;
class G4PropagatorInField;
class G4FieldManager;
class G4UniformMagField;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

//////////////////////////////////////////////////////////
//                  DETECTOR ARRAY SETUP                //
//////////////////////////////////////////////////////////

///////////////     CLOVER DETECTORS     ///////////////////
const G4int     numberOf_CLOVER = 9;
const G4int     numberOf_CLOVER_Shields = 9;

///////////////     TIGRESS DETECTORS     ///////////////////
const G4int     numberOf_TIGRESS = 1;
const G4int     numberOf_TIGRESS_BGO = 1;

///////////////     PLASTIC SCINTILLATOR DETECTORS     ///////////////////
const G4int     numberOf_PlasticScint = 12;

///////////////     LEPS DETECTORS     ///////////////////
const G4int     numberOf_LEPS = 8;


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class DetectorConstruction : public G4VUserDetectorConstruction
{
  public:
    DetectorConstruction();
    virtual ~DetectorConstruction();

  public:
    virtual G4VPhysicalVolume* Construct();
    virtual void ConstructField();

    // get methods
    //
    //const G4VPhysicalVolume* GetAbsorberPV() const;
    //const G4VPhysicalVolume* GetGapPV() const;
     
  private:
    // methods
    //
    void DefineMaterials();
    G4VPhysicalVolume* DefineVolumes();
  
    // data members
    //
    static G4ThreadLocal G4GlobalMagFieldMessenger*  fMagFieldMessenger;
                                      // magnetic field messenger
     
    G4VPhysicalVolume*   fAbsorberPV; // the absorber physical volume
    G4VPhysicalVolume*   fGapPV;      // the gap physical volume
    
    G4bool  fCheckOverlaps; // option to activate checking of volumes overlaps
    
    /////////////////////////////
    //          WORLD
    G4double WorldSize;
    
    
    /////////////////////////////////////
    //          CLOVER DETECTORS
    /////////////////////////////////////
    
    G4bool              CLOVER_AllPresent_Override;
    G4bool              CLOVER_AllAbsent_Override;
    G4bool              CLOVER_Presence[numberOf_CLOVER];
    G4double            CLOVER_Distance[numberOf_CLOVER];
    G4RotationMatrix    CLOVER_rotm[numberOf_CLOVER];
    G4Transform3D       CLOVER_transform[numberOf_CLOVER];
    G4ThreeVector       CLOVER_position[numberOf_CLOVER];
    G4double            CLOVER_phi[numberOf_CLOVER];
    G4double            CLOVER_theta[numberOf_CLOVER];
    
    //  CLOVER HPGe Crystals
    G4VPhysicalVolume*  PhysiCLOVER_HPGeCrystal;
    
    ///////////////////////////////////////////////////////////////
    //          CLOVER - BGO Shield   (Manufacturer: Cyberstar)
    ///////////////////////////////////////////////////////////////
    
    G4bool              CLOVER_Shield_AllPresent_Override;
    G4bool              CLOVER_Shield_AllAbsent_Override;
    G4bool              CLOVER_Shield_Presence[numberOf_CLOVER_Shields];
    G4ThreeVector       CLOVER_Shield_position[numberOf_CLOVER_Shields];
    G4Transform3D       CLOVER_Shield_transform[numberOf_CLOVER_Shields];
    
    //  Shield BGO Crystal Scintillators
    G4VPhysicalVolume*  PhysiCLOVER_Shield_BGOCrystal;
    
    //  Shield PMT Tubes
    G4VPhysicalVolume*  PhysiCLOVER_Shield_PMT;

    
    ///////////////////////////////
    //      PlasticScint DETECTORS
    ///////////////////////////////
    
    G4bool              PlasticScint_AllPresent_Override;
    G4bool              PlasticScint_AllAbsent_Override;
    G4bool              PlasticScint_Presence[numberOf_PlasticScint];
    G4RotationMatrix    PlasticScint_rotm[numberOf_PlasticScint];
    G4Transform3D       PlasticScint_transform[numberOf_PlasticScint];
    G4ThreeVector       PlasticScint_CentrePosition[numberOf_PlasticScint];
    G4double            PlasticScint_CentrePositionX[numberOf_PlasticScint];
    G4double            PlasticScint_CentrePositionY[numberOf_PlasticScint];
    G4double            PlasticScint_CentrePositionZ[numberOf_PlasticScint];
    G4double            PlasticScint_RotationY[numberOf_PlasticScint];
    
    G4VPhysicalVolume*  PhysiPlasticScint;
    
    /////////////////////////////////////
    //              HAGAR
    /////////////////////////////////////
    
    G4bool              HAGAR_NaICrystal_Presence;
    G4bool              HAGAR_Annulus_Presence;
    G4bool              HAGAR_FrontDisc_Presence;
    
    G4ThreeVector       HAGAR_NaICrystal_CentrePosition;
    G4ThreeVector       HAGAR_Annulus_CentrePosition;
    G4ThreeVector       HAGAR_FrontDisc_CentrePosition;
    G4RotationMatrix    HAGAR_rotm;
    G4Transform3D       HAGAR_transform;
    
    //  HAGAR NaI Crystal
    G4VPhysicalVolume*  PhysiHAGAR_NaICrystal;
    
    //  HAGAR Annulus
    G4VPhysicalVolume*  PhysiHAGAR_Annulus;
    
    //  HAGAR Front Disc
    G4VPhysicalVolume*  PhysiHAGAR_FrontDisc;
    
    
    
    
    //////////////////////////////////////
    //          K600 SPECTROMETER
    //////////////////////////////////////
    
    //////////////////////////////////////
    //          K600 - QUADRUPOLE
    G4bool              Ideal_Quadrupole;
    G4bool              Mapped_Quadrupole;
    G4bool              AFRODITE_Quadrupole;
    
    G4VPhysicalVolume*  PhysiAFRODITE_Quadrupole;
    
    G4ThreeVector       AFRODITE_Quadrupole_CentrePosition;
    G4RotationMatrix    AFRODITE_Quadrupole_rotm;
    G4Transform3D       AFRODITE_Quadrupole_transform;
    
    ////    MAGNETIC FIELD for QUADRUPOLE
    static G4ThreadLocal G4FieldManager* fieldManagerMagneticField_AFRODITE_Q;
    static G4ThreadLocal G4QuadrupoleMagField* MagneticField_AFRODITE_Q;
    
    G4double                dBdr_gradient_AFRODITE_Q;   // gradient = dB/dr
    G4Mag_UsualEqRhs*       fEquationMagneticField_AFRODITE_Q;
    G4MagIntegratorStepper* stepperMagneticField_AFRODITE_Q;
    G4ChordFinder*          fChordFinder_AFRODITE_Q;
    
    
    //////////////////////////////////////
    //          AFRODITE - DIPOLE 1
    G4bool              AFRODITE_Dipole1;
    G4VPhysicalVolume*  PhysiAFRODITE_Dipole1;
    
    G4ThreeVector       AFRODITE_Dipole1_CentrePosition;
    G4RotationMatrix    AFRODITE_Dipole1_rotm;
    G4Transform3D       AFRODITE_Dipole1_transform;
    
    ////    MAGNETIC FIELD for DIPOLE 1
    static G4ThreadLocal G4FieldManager* fieldManagerMagneticField_AFRODITE_D1;
    static G4ThreadLocal G4UniformMagField* MagneticField_AFRODITE_D1;
    
    G4double                AFRODITE_Dipole1_BZ;
    G4Mag_UsualEqRhs*       fEquationMagneticField_AFRODITE_D1;
    G4double                minStepMagneticField;
    G4MagIntegratorStepper* stepperMagneticField_AFRODITE_D1;
    G4ChordFinder*          fChordFinder_AFRODITE_D1;
    
    
    
    ////////////////////////////////
    ////        STRUCTURES      ////
    ////////////////////////////////
    
    /////////////////////////////////////
    //  K600 BACTAR Target Chamber
    G4bool      AFRODITE_MathisTC_Presence;
    
    /////////////////////////////////////
    //  K600 Target
    G4bool      AFRODITE_Target_Presence;
    
    
};



//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

