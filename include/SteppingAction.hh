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

#ifndef SteppingAction_h
#define SteppingAction_h 1

#include "G4UserSteppingAction.hh"
#include "globals.hh"
#include "G4ThreeVector.hh"

class DetectorConstruction;
class EventAction;

/// Stepping action class.
///
/// In UserSteppingAction() there are collected the energy deposit and track 
/// lengths of charged particles in Absober and Gap layers and
/// updated in EventAction.

class SteppingAction : public G4UserSteppingAction
{
public:
  SteppingAction(const DetectorConstruction* detectorConstruction,
                    EventAction* eventAction);
  virtual ~SteppingAction();

  virtual void UserSteppingAction(const G4Step* step);
    
private:
  const DetectorConstruction* fDetConstruction;
  EventAction*  fEventAction;
    
    G4double    fCharge;
    G4double    fMass;
    G4ThreeVector worldPosition;
    G4ThreeVector localPosition;
    
    
    ////    PADDLE DETECTOR - Plastic Scintillator
    G4double    edepPADDLE;
    G4int       PADDLE_ITS;  // Interaction Time Sample for the Neutron Wall Plastic Scintillators
    G4int       PADDLENo;
    
    ////    VDC DETECTOR - Plastic Scintillator
    G4double    edepVDC;
    G4int       VDC_ITS;  // Interaction Time Sample for the Neutron Wall Plastic Scintillators
    G4int       WireChamberNo;
    G4int       hit_StoredChannelNo;

    //  Local Position
    G4double    xPosL;
    G4double    yPosL;
    G4double    zPosL;
    
    //  World Position
    G4double    xPosW;
    G4double    yPosW;
    G4double    zPosW;
    
    G4double    xShift = 4*(cos(40) + tan(40)*cos(50));
    G4double    xOffset;
    
    ////    TIARA DETECTOR
    G4double    edepTIARA_AA;
    G4int       TIARA_AA_ITS;             // Interaction Time Sample for the TIARA Detectors
    G4int       TIARANo;
    G4int       TIARA_RowNo;
    G4int       TIARA_SectorNo;
    
    ////    CLOVER DETECTOR
    G4double    edepCLOVER_HPGeCrystal;
    G4int       CLOVER_HPGeCrystal_ITS;   // Interaction Time Sample for the CLOVER Detectors
    G4int       CLOVERNo;
    G4int       CLOVER_HPGeCrystalNo;
    
    ////    CLOVER BGO-Crystal, Compton Supression Shield
    G4double    edepCLOVER_BGOCrystal;
    G4int       CLOVER_BGO_ITS; // Interaction Time Sample for the BGO Crystals within the CLOVER Shield
    
    G4double    interactiontime;
    G4int       channelID;
    G4String    volumeName;
    
    
    //////////////////////////////////
    //      GEOMETRY ANALYSIS
    //////////////////////////////////
    
    G4double normVector, theta, phi;
    
    
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
