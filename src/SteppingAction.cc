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

#include "SteppingAction.hh"
#include "EventAction.hh"
#include "DetectorConstruction.hh"
#include "G4SystemOfUnits.hh"

#include "G4Step.hh"
#include "G4RunManager.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

SteppingAction::SteppingAction(
                      const DetectorConstruction* detectorConstruction,
                      EventAction* eventAction)
  : G4UserSteppingAction(),
    fDetConstruction(detectorConstruction),
    fEventAction(eventAction)
{
    
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

SteppingAction::~SteppingAction()
{ 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void SteppingAction::UserSteppingAction(const G4Step* aStep)
{
    G4StepPoint* preStepPoint = aStep->GetPreStepPoint();
    G4TouchableHandle theTouchable = preStepPoint->GetTouchableHandle();
    
    // get particle name/definition
    //G4String particleName = aStep->GetTrack()->GetDefinition()->GetParticleName();
    //G4ParticleDefinition* particle = aStep->GetTrack()->GetDefinition();
    
    // get particle lifetime
    //G4double lifetime = aStep->GetTrack()->GetDefinition()->GetIonLifeTime()/ns;
    
    // get interaction time of the current step
    interactiontime = preStepPoint->GetGlobalTime()/ns;
    
    // get volume of the current step
    G4VPhysicalVolume* volume = theTouchable->GetVolume();
    
    // get volume name of the current step
    volumeName = volume->GetName();


    //G4cout << "Here is the particleName    "<< particleName << G4endl;
    //G4cout << "Here is the lifetime    "<< lifetime << G4endl;
    //G4cout << "Here is the interactiontime    "<< interactiontime << G4endl;
    //G4cout << "Here is the TEST    "<< TEST << G4endl;
    //G4cout << "Here is the interactiontime    "<< interactiontime << G4endl;
    //G4cout << "                                "<< G4endl;
    

    //G4ParticleDefinition* particleOI = G4Gamma::Gamma();
    
    
    ////////////////////////////////////////////
    //              TIARA ARRAY
    ////////////////////////////////////////////
    
    if(interactiontime < TIARA_TotalSampledTime && volumeName == "TIARA_AA_RS")
    {
        edepTIARA_AA = aStep->GetTotalEnergyDeposit()/MeV;

        if(edepTIARA_AA != 0)
        {
            channelID = volume->GetCopyNo();
            
            TIARANo = channelID/128;
            TIARA_RowNo = (channelID - (TIARANo*128))/8;
            TIARA_SectorNo = (channelID - (TIARANo*128))%8;

            TIARA_AA_ITS = interactiontime/TIARA_SamplingTime;
            edepTIARA_AA = aStep->GetTotalEnergyDeposit()/MeV;
            
            if(fEventAction->GetVar_TIARA_AA(TIARANo, TIARA_RowNo, TIARA_SectorNo, 0, TIARA_AA_ITS)==0)
            {
                worldPosition = preStepPoint->GetPosition();
                
                xPosW = worldPosition.x()/m;
                yPosW = worldPosition.y()/m;
                zPosW = worldPosition.z()/m;
                
                normVector = pow(pow(xPosW,2) + pow(yPosW,2) + pow(zPosW,2) , 0.5);
                theta = acos(zPosW/normVector)/deg;
                
                if(xPosW==0)
                {
                    if(yPosW==0) phi = 0;
                    if(yPosW>0) phi = 90;
                    if(yPosW<0) phi = 270;
                }
                else
                {
                    phi = atan(yPosW/xPosW)/deg;
                    
                    if(xPosW>0 && yPosW>0) phi = phi; // deg
                    if(xPosW<0 && yPosW>0) phi = phi + 180.; // deg
                    if(xPosW<0 && yPosW<0) phi = phi + 180.; // deg
                    if(xPosW>0 && yPosW<0) phi = phi + 360.; // deg
                }
                
                fEventAction->SetVar_TIARA_AA(TIARANo, TIARA_RowNo, TIARA_SectorNo, 1, TIARA_AA_ITS, theta);
                fEventAction->SetVar_TIARA_AA(TIARANo, TIARA_RowNo, TIARA_SectorNo, 2, TIARA_AA_ITS, phi);

            }
            
            fEventAction->FillVar_TIARA_AA(TIARANo, TIARA_RowNo, TIARA_SectorNo, 0, TIARA_AA_ITS, edepTIARA_AA);
        }
    }
    
    
    ////////////////////////////////////////////
    //          PlasticScint DETECTORS
    ////////////////////////////////////////////
    
    if (interactiontime < PlasticScint_TotalSampledTime)
    {
        if (volumeName == "PlasticScint")
        {
            channelID = volume->GetCopyNo();
            
            PlasticScintNo = channelID;
            
            PlasticScint_ITS = interactiontime/PlasticScint_SamplingTime;
            edepPlasticScint = aStep->GetTotalEnergyDeposit()/MeV;
            
            worldPosition = preStepPoint->GetPosition();
            localPosition = theTouchable->GetHistory()->GetTopTransform().TransformPoint(worldPosition);
            
            fEventAction->AddEnergy_PlasticScint( PlasticScintNo, PlasticScint_ITS, edepPlasticScint);
            fEventAction->TagTOF_PlasticScint(PlasticScintNo, PlasticScint_ITS, interactiontime);
            fEventAction->AddEWpositionX_PlasticScint( PlasticScintNo, PlasticScint_ITS, edepPlasticScint*localPosition.x());
            fEventAction->AddEWpositionY_PlasticScint( PlasticScintNo, PlasticScint_ITS, edepPlasticScint*localPosition.y());
            
            //if(fEventAction->Get_PlasticScint_Trig(i) == false) fEventAction->Set_PlasticScint_Trig(i, true);
        }
    }
    
    
    
    ////////////////////////////////////////////////
    //                  CLOVERS
    ////////////////////////////////////////////////
    
    if(interactiontime < CLOVER_TotalSampledTime)
    {
        if (volumeName == "CLOVER_HPGeCrystal")
        {
            channelID = volume->GetCopyNo();
            
            CLOVERNo = channelID/4;
            CLOVER_HPGeCrystalNo = channelID%4;
            
            /*
            G4cout << "Here is the copyNo    "<< copyNo << G4endl;
            G4cout << "Here is the CLOVERNo    "<< CLOVERNo << G4endl;
            G4cout << "Here is the CLOVER_HPGeCrystalNo    "<< CLOVER_HPGeCrystalNo << G4endl;
            G4cout << " "<< G4endl;
            */
            
            CLOVER_HPGeCrystal_ITS = interactiontime/CLOVER_SamplingTime;
            edepCLOVER_HPGeCrystal = aStep->GetTotalEnergyDeposit()/keV;
            
            fEventAction->AddEnergyCLOVER_HPGeCrystal(CLOVERNo, CLOVER_HPGeCrystalNo, CLOVER_HPGeCrystal_ITS, edepCLOVER_HPGeCrystal);
        }
    }
    
    /*
    if (interactiontime < CLOVER_Shield_BGO_TotalSampledTime)
    {
        for(G4int i=0; i<8; i++)
        {
            for(G4int l=0; l<16; l++)
            {
                if ( volume == fDetConstruction->GetVolume_CLOVER_Shield_BGOCrystal(i, l) && interactiontime <     CLOVER_Shield_BGO_TotalSampledTime )
                {
                    CLOVER_BGO_ITS = interactiontime/CLOVER_Shield_BGO_SamplingTime;
                    edepCLOVER_BGOCrystal = aStep->GetTotalEnergyDeposit()/keV;
                    
                    fEventAction->AddEnergyBGODetectors(i, l, CLOVER_BGO_ITS, edepCLOVER_BGOCrystal);
                    //G4cout << "Here is the edepCLOVER_BGOCrystal    "<< edepBGO << G4endl;
                }
            }
        }
    }
    */
    
    
    //if (volumeName=="TIARA_Assembly" && !fEventAction->GA_GetLineOfSight() )   G4cout << "Here is the TIARA_Assembly Hit!" << G4endl;
    
    ////////////////////////////////////////////
    //              TIARA ARRAY
    ////////////////////////////////////////////
    
    if(GA_MODE)
    {
        
        if(((volumeName=="TIARA_AA_RS" || volumeName=="TIARA_SiliconWafer") && ((GA_LineOfSightMODE && fEventAction->GA_GetLineOfSight()==true) || !GA_LineOfSightMODE)) || (volumeName == "World" && GA_GenInputVar))
        {
            channelID = volume->GetCopyNo();
            worldPosition = preStepPoint->GetPosition();
            //worldPosition = worldPosition.unit();
 
            xPosW = worldPosition.x()/m;
            yPosW = worldPosition.y()/m;
            zPosW = worldPosition.z()/m;
            
            /*
            if(volumeName == "TIARA_SiliconWafer")
            {
                G4cout << " " << G4endl;
                G4cout << "Here is the TIARA_SiliconWafer HIT" << G4endl;
            }
            
            if(volumeName == "TIARA_AA_RS")
            {
                G4cout << " " << G4endl;
                G4cout << "Here is the TIARA_AA_RS HIT" << G4endl;
            }
            */
            
            if(volumeName == "TIARA_AA_RS")
            {
                fEventAction->FillGA_TIARAstor(channelID, 0, xPosW);
                fEventAction->FillGA_TIARAstor(channelID, 1, yPosW);
                fEventAction->FillGA_TIARAstor(channelID, 2, zPosW);
                fEventAction->FillGA_TIARAstor(channelID, 3, 1.);
            }
            
            if(GA_GenAngDist && fEventAction->GetGA_TIARA(channelID, 0)==0)
            {
                normVector = pow(pow(xPosW,2) + pow(yPosW,2) + pow(zPosW,2) , 0.5);
                theta = acos(zPosW/normVector)/deg;
                
                if(xPosW==0)
                {
                    if(yPosW==0) phi = 0;
                    if(yPosW>0) phi = 90;
                    if(yPosW<0) phi = 270;
                }
                else
                {
                    phi = atan(yPosW/xPosW)/deg;
                    
                    if(xPosW>0 && yPosW>0) phi = phi; // deg
                    if(xPosW<0 && yPosW>0) phi = phi + 180.; // deg
                    if(xPosW<0 && yPosW<0) phi = phi + 180.; // deg
                    if(xPosW>0 && yPosW<0) phi = phi + 360.; // deg
                }
                
                if(volumeName == "TIARA_AA_RS")
                {
                    fEventAction->SetGA_TIARA(channelID, 0, 1);
                    fEventAction->SetGA_TIARA(channelID, 1, theta);
                    fEventAction->SetGA_TIARA(channelID, 2, phi);
                }
            }
            
            if(GA_GenInputVar && volumeName == "World")
            {
                fEventAction->SetInputDist(0, theta);
                fEventAction->SetInputDist(1, phi);
            }
            //G4cout << "Here is the geantino Hit!     -->     " << G4endl;
        }
    }

    
    ////    Here, one declares the volumes that one considers will block the particles of interest and effectively mask the relevant volume of interest.
    if (GA_LineOfSightMODE && (volumeName == "TIARA_AA_RS" || volumeName=="TIARA_PCB" || volumeName=="TIARA_SiliconWafer"))
    {
        //G4cout << "Here is the volumeName    "<< volumeName << G4endl;
        fEventAction->GA_SetLineOfSight(false);
    }
    
    
    
    /*
// Collect energy and track length step by step

  // get volume of the current step
  G4VPhysicalVolume* volume
    = step->GetPreStepPoint()->GetTouchableHandle()->GetVolume();
  
  // energy deposit
  G4double edep = step->GetTotalEnergyDeposit();
  
    
  // step length
  G4double stepLength = 0.;
  if ( step->GetTrack()->GetDefinition()->GetPDGCharge() != 0. ) {
    stepLength = step->GetStepLength();
  }
      
  if ( volume == fDetConstruction->GetAbsorberPV() ) {
    fEventAction->AddAbs(edep,stepLength);
  }
  
  if ( volume == fDetConstruction->GetGapPV() ) {
    fEventAction->AddGap(edep,stepLength);
  }
    */
    
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
