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

#include "EventAction.hh"
#include "RunAction.hh"
#include "Analysis.hh"

#include "G4RunManager.hh"
#include "G4Event.hh"
#include "G4UnitsTable.hh"

#include "Randomize.hh"
#include <iomanip>

#include <fstream>
#include <string>

using namespace std;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

EventAction::EventAction()
 : G4UserEventAction(),
   fEnergyAbs(0.),
   fEnergyGap(0.),
   fTrackLAbs(0.),
   fTrackLGap(0.)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

EventAction::~EventAction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void EventAction::BeginOfEventAction(const G4Event* evt)
{  
  // initialisation per event
  //fEnergyAbs = 0.;
  //fEnergyGap = 0.;
  //fTrackLAbs = 0.;
  //fTrackLGap = 0.;
    
    evtNb = evt->GetEventID();
    
    GA_LineOfSight = true;
    
    if(GA_MODE)
    {
        if(evtNb==0)
        {
            for(G4int i=0; i<5; i++)
            {
                for(G4int j=0; j<16; j++)
                {
                    for(G4int k=0; k<8; k++)
                    {
                        GA_MMM_AngDist_counter[i][j][k]= 0;
                    }
                }
            }
            
            
            for(G4int i=0; i<640; i++)
            {
                for(G4int j=0; j<4; j++)
                {
                    GA_TIARA_AA_stor[i][j] = 0;
                }
            }
        }
        
        for(G4int i=0; i<640; i++)
        {
            for(G4int j=0; j<3; j++)
            {
                GA_TIARA_AA[i][j] = 0;
            }
        }
    }
    
    
    /////////////////////////////////////////////////////
    
    for(G4int i=0; i<9; i++)
    {
        for (G4int k=0; k<CLOVER_TotalTimeSamples; k++)
        {
            CLOVER_EDep[i][k] = 0;
            
            for(G4int j=0; j<4; j++)
            {
                CLOVER_HPGeCrystal_EDep[i][j][k] = 0;
                CLOVER_BGO_EDep[i][j][k] = 0;
                CLOVER_HPGeCrystal_EDepVETO[i][j][k] = false;
            }
        }
        
        for (G4int m=0; m<CLOVER_Shield_BGO_TotalTimeSamples; m++)
        {
            for(G4int l=0; l<16; l++)
            {
                CLOVER_BGO_EDep[i][l][m] = 0;
            }
        }
    }
    
    for(G4int i=0; i<5; i++)
    {
        for(G4int k = 0; k<TIARA_TotalTimeSamples; k++)
        {
            for(G4int j=0; j<16; j++)
            {
                for(G4int l=0; l<8; l++)
                {
                    for(G4int m=0; m<8; m++)
                    {
                        TIARA_AA[i][j][l][m][k] = 0;
                    }
                }
            }
        }
    }
    
    for(G4int i=0; i<3; i++)
    {
        for (G4int k=0; k<PADDLE_TotalTimeSamples; k++)
        {
            PADDLE_Trig[i] = 0;
            
            PADDLE_EDep[i][k] = 0;
            PADDLE_TOF[i][k] = 0;
            
            PADDLE_EWpositionX[i][k] = 0;
            PADDLE_EWpositionY[i][k] = 0;
            
            PADDLE_positionX[i][k] = 0;
            PADDLE_positionY[i][k] = 0;
        }
    }
    
    for(G4int j=0; j<4; j++)
    {
        for (G4int k=0; k<hit_buffersize; k++)
        {
            if(j==0) VDC_Observables[j][k] = -1;
            else{VDC_Observables[j][k] = 0;}
        }
    }
    
    ////    Input Variables
    InputDist[0] = 0;
    InputDist[1] = 0;

    
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void EventAction::EndOfEventAction(const G4Event* event)
{
  // Accumulate statistics
  //

  // get analysis manager
  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();

    
    ////////////////////////////////////////////////////////
    //
    //                TIARA DETECTOR ARRAY
    //
    ////////////////////////////////////////////////////////
    
    for(G4int i=0; i<5; i++)
    {
        shift_phi = (double)(i)*(10./8.)*(53.79-0.29);
        //shift_phi = (double)(i)*72.;
        
        for(G4int k = 0; k<TIARA_TotalTimeSamples; k++)
        {
            for(G4int j=0; j<16; j++)
            {
                for(G4int l=0; l<8; l++)
                {
                    ////    Processing the Energies of the TIARA Detectors
                    TIARA_AA[i][j][l][0][k] = G4RandGauss::shoot(TIARA_AA[i][j][l][0][k], 0.7);

                    if(TIARA_AA[i][j][l][0][k] >= G4RandGauss::shoot(TIARA_AA_ThresholdEnergy, 0.5))
                    {
                        ////      Counts versus Energy for each TIARA
                        analysisManager->FillH1(1+i, GainTIARA*TIARA_AA[i][j][l][0][k] + OffsetTIARA);
                        
                        ////      Counts versus Energy for the Entire TIARA Array
                        analysisManager->FillH1(6, GainTIARA*TIARA_AA[i][j][l][0][k] +  OffsetTIARA);
                        
                        ////////////////////////////////////////////////////////////
                        ////                Filling DataTreeSim
                        //analysisManager->FillNtupleIColumn(0, 0, i);
                        //analysisManager->FillNtupleIColumn(0, 1, j);
                        //analysisManager->FillNtupleIColumn(0, 2, k);
                        /*
                        //      Energy
                        analysisManager->FillNtupleDColumn(0, 12, TIARA_AA[i][j][l][0][k]);
                        //      Theta
                        analysisManager->FillNtupleDColumn(0, 13, TIARA_AA[i][j][l][1][k]);
                        //      Phi
                        analysisManager->FillNtupleDColumn(0, 14, TIARA_AA[i][j][l][2][k]);
                         */
                    }
                }
            }
        }
    }

    
    ////////////////////////////////////////////////////////
    //
    //          PADDLE DETECTORS - Plastic Scintillators
    //
    ////////////////////////////////////////////////////////
    
    GainPADDLE = 1.0;
    OffsetPADDLE = 0.0;

    
    for(G4int i=0; i<3; i++)
    {
        for (G4int k=0; k<PADDLE_TotalTimeSamples; k++)
        {
            ////              Calculating energy weighted positions
            PADDLE_positionX[i][k] = G4RandGauss::shoot(PADDLE_EWpositionX[i][k]/PADDLE_EDep[i][k], 4.8);
            PADDLE_positionY[i][k] = PADDLE_EWpositionY[i][k]/PADDLE_EDep[i][k];
            
            ////              Calculating a Gaussian Smeared Energy Deposition
            PADDLE_EDep[i][k] = G4RandGauss::shoot(PADDLE_EDep[i][k], 0.10*PADDLE_EDep[i][k]);
            
            if( PADDLE_EDep[i][k] >= G4RandGauss::shoot(PADDLE_ThresholdEnergy, 0.01*PADDLE_ThresholdEnergy))
            {
                ////////////////////////////////////////////////////////
                //      PADDLE DETECTORS - 1D, Counts versus Energy
                ////////////////////////////////////////////////////////
                
                analysisManager->FillH1(i+7, GainPADDLE*PADDLE_EDep[i][k] + OffsetPADDLE, 1);
                
                PADDLE_TOF[i][k] = G4RandGauss::shoot(PADDLE_TOF[i][k], 0.05*PADDLE_TOF[i][k]);
                
                
                ////////////////////////////////////////////////////////////////////
                //              PADDLE DETECTORS - 2D, Position versus Energy
                ////////////////////////////////////////////////////////////////////
                analysisManager->FillH2(i+1, PADDLE_positionX[i][k], PADDLE_positionY[i][k], PADDLE_EDep[i][k]);
                
                ////////////////////////////////////////////////////////////////////
                //              PADDLE DETECTORS - 2D, Energy versus T.O.F.
                ////////////////////////////////////////////////////////////////////
                
                analysisManager->FillH2(i+4, PADDLE_TOF[i][k], GainPADDLE*PADDLE_EDep[i][k] + OffsetPADDLE, 1);
                
            }
        }
    }
    
    
    
    
    ////////////////////////////////////////////////////////
    //
    //                CLOVER DETECTOR ARRAY
    //
    ////////////////////////////////////////////////////////
    
    for(G4int i=0; i<9; i++)
    {
        for(G4int k=0; k<CLOVER_TotalTimeSamples; k++)
        {
            for(G4int j=0; j<4; j++)
            {
                if(G4RandGauss::shoot(CLOVER_HPGeCrystal_EDep[i][j][k], 0.7) >= CLOVER_HPGeCrystal_ThresholdEnergy)
                {
                    CLOVER_HPGeCrystal_EDep[i][j][k] = G4RandGauss::shoot(CLOVER_HPGeCrystal_EDep[i][j][k], 1.7);
                    
                    if(Activate_CLOVER_ComptonSupression)
                    {
                        for(G4int l=0; l<CLOVER_ComptonSupression_TimeWindow; l++)
                        {
                            for(G4int m=0; m<16; m++)
                            {
                                //      COMPTON SUPRESSION - VETO CLOVER Energy Depositions in anti-coincidence with BGO Shield Energy Deposition
                                if (CLOVER_BGO_EDep[i][m][k+l] >= CLOVER_BGO_ThresholdEnergy)
                                {
                                    CLOVER_HPGeCrystal_EDepVETO[i][j][k] = true;
                                }
                            }
                        }
                        if (CLOVER_HPGeCrystal_EDepVETO[i][j][k]) CLOVER_HPGeCrystal_EDep[i][j][k] = 0;
                    }
                    
                    if(Activate_CLOVER_ADDBACK && CLOVER_HPGeCrystal_EDep[i][j][k] != 0)
                    {
                        //      ADDBACK
                        CLOVER_EDep[i][k] += CLOVER_HPGeCrystal_EDep[i][j][k];
                        
                        //      For each Clover
                        analysisManager->FillH1(i+10, GainCLOVER*CLOVER_EDep[i][k] + OffsetCLOVER);
                        
                        //      For the Entire Clover Array
                        analysisManager->FillH1(19, GainCLOVER*CLOVER_EDep[i][k] +  OffsetCLOVER);
                    }
                    
                    else if(CLOVER_HPGeCrystal_EDep[i][j][k] != 0)
                    {
                        //      For each Clover
                        analysisManager->FillH1(i+10, GainCLOVER*CLOVER_HPGeCrystal_EDep[i][j][k] + OffsetCLOVER);
                        
                        //      For the Entire Clover Array
                        analysisManager->FillH1(19, GainCLOVER*CLOVER_HPGeCrystal_EDep[i][j][k] +  OffsetCLOVER);
                    }
                }
            }
        }
    }
    
    
    
    ////////////////////////////////////////////////////////
    //
    //                VDC, VERTICAL DRIFT CHAMBER
    //
    ////////////////////////////////////////////////////////
    
    ////    VDC 1
    RayTrace(0, 0);     //RayTrace(VDCNo, XU_Wireplane)
    RayTrace(0, 1);
    CalcYFP(0);
    //G4cout << "Here is the Xpos[0] (VDC1)     -->     "<< Xpos[0] << G4endl;
    //G4cout << "Here is the ThetaFP[0] (VDC1)     -->     "<< ThetaFP[0] << G4endl;

    ////    VDC 2
    RayTrace(1, 0);     //RayTrace(VDCNo, XU_Wireplane)
    RayTrace(1, 1);
    CalcYFP(1);
    //G4cout << "Here is the Xpos[1] (VDC2)     -->     "<< Xpos[1] << G4endl;
    //G4cout << "Here is the ThetaFP[1] (VDC2)     -->     "<< ThetaFP[1] << G4endl;

    
    
    ////////////////////////////////////////////////////
    ////                                            ////
    ////                DataTreeSim                 ////
    ////                                            ////
    ////////////////////////////////////////////////////
    /*
    analysisManager->FillNtupleDColumn(0, 0, Xpos[0]);
    analysisManager->FillNtupleDColumn(0, 1, Y[0]);
    analysisManager->FillNtupleDColumn(0, 2, ThetaFP[0]);
    analysisManager->FillNtupleDColumn(0, 3, ThetaSCAT[0]);

    analysisManager->FillNtupleDColumn(0, 4, Xpos[1]);
    analysisManager->FillNtupleDColumn(0, 5, Y[1]);
    analysisManager->FillNtupleDColumn(0, 6, ThetaFP[1]);
    analysisManager->FillNtupleDColumn(0, 7, ThetaSCAT[1]);
    
    analysisManager->AddNtupleRow(0);

    */
    
    

    
    ////////////////////////////////////////////////////////////////
    ////                                                        ////
    ////                    GEOMETRY ANALYSIS                   ////
    ////                                                        ////
    ////////////////////////////////////////////////////////////////
    
    if(GA_MODE)
    {
        
        ////////////////////////////////
        ////    Input Variables
        analysisManager->FillNtupleDColumn(2, 0, InputDist[0]);
        analysisManager->FillNtupleDColumn(2, 1, InputDist[1]);
        
        analysisManager->AddNtupleRow(2);

        //G4cout << "Here is the value of InputDist[0]:    -->     " << InputDist[0] << G4endl;
        //G4cout << "Here is the value of InputDist[1]:    -->     " << InputDist[1] << G4endl;

        
        
        ////////////////////////////////////////////////////////////
        ////    Creating Distribution txt file for Data Sorting
        
        if(GA_GenAngDist)
        {
            if(evtNb==0)
            {
                sprintf(filenameV,"K600Veridical/MMM/AngDist.txt");
                fileNameHolder = string(filenameV);
                fileV_MMM.open(fileNameHolder, std::ios_base::app);
            }
            
            for(G4int i=0; i<512; i++)  //  i<640 for all 5 silicons
            {
                if(GA_TIARA_AA[i][0]!=0)
                {
                    TIARANo = i/128;
                    TIARA_RowNo = (i - (TIARANo*128))/8;
                    TIARA_SectorNo = (i - (TIARANo*128))%8;
                    
                    ////////////////////////////////////////////////////////////
                    ////            Filling GeometryAnalysisTree
                    
                    analysisManager->FillNtupleIColumn(1, 0, TIARANo);
                    analysisManager->FillNtupleIColumn(1, 1, TIARA_RowNo);
                    analysisManager->FillNtupleIColumn(1, 2, TIARA_SectorNo);
                    
                    //      Theta
                    analysisManager->FillNtupleDColumn(1, 3, GA_TIARA_AA[i][1]);
                    //      Phi
                    analysisManager->FillNtupleDColumn(1, 4, GA_TIARA_AA[i][2]);
                    
                    analysisManager->AddNtupleRow(1);
                    
                    fileV_MMM << TIARANo << "    " << TIARA_RowNo << "    " << TIARA_SectorNo << "    " << GA_TIARA_AA[i][1] << "    " << GA_TIARA_AA[i][2] << endl;
                    
                }
            }
            
            ////    Closing the Distribution File
            if((evtNb+1)%GA_numberOfEvents==0 )
            {
                fileV_MMM.close();
            }
            
            
            
        }
        
        
        
        
        /////////////////////////////////////////////////////////////////////
        ////    Writing the Output files for the GEOMETRY ANALYSIS mode
        if(evtNb == (GA_numberOfEvents-1))
        {
            G4double av_xPos, av_yPos, av_zPos;
            G4double normVector, theta, phi, solidAngle;
            
            // append text file
            G4String fileName1 = "K600SimOutput.txt";
            G4String fileName2 = "K600SimOutput.h";
            
            std::ofstream file1;
            file1.open(fileName1, std::ios_base::app);
            
            std::ofstream file2;
            file2.open(fileName2, std::ios_base::app);
            
            file1 << "(TIARA NUMBER)  (ROW NUMBER)  (SECTOR NUMBER)  (THETA)      (PHI)      (SOLID ANGLE)"<< endl;
            file2 << "                  "<< endl;
            file2 << "Double_t GA_TIARA[5][16][8][3];"<< endl;
            file2 << "                  "<< endl;
            file2 << "void initialize_GA()"<< endl;
            file2 << "{"<< endl;
            
            G4double GA_numberOfEvents_double = GA_numberOfEvents;
            
            ////    For the Silicon Array
            for(G4int i=0; i<512; i++)  //  i<640 for all 5 silicons
            {
                if( i%128 == 0) file1 << "   " << endl;
                
                TIARANo = i/128;
                TIARA_RowNo = (i - (TIARANo*128))/8;
                TIARA_SectorNo = (i - (TIARANo*128))%8;
                
                av_xPos = GA_TIARA_AA_stor[i][0]/GA_TIARA_AA_stor[i][3];
                av_yPos = GA_TIARA_AA_stor[i][1]/GA_TIARA_AA_stor[i][3];
                av_zPos = GA_TIARA_AA_stor[i][2]/GA_TIARA_AA_stor[i][3];
                
                normVector = pow(pow(av_xPos,2) + pow(av_yPos,2) + pow(av_zPos,2) , 0.5);
                theta = acos(av_zPos/normVector)/deg;
                solidAngle = (0.5)*(GA_TIARA_AA_stor[i][3]/GA_numberOfEvents_double);
                
                /*
                G4cout << "" << G4endl;
                G4cout << "Here is the GA_TIARA_AA_stor[3] value     -->     "<< GA_TIARA_AA_stor[i][3] << G4endl;
                G4cout << "Here is the GA_numberOfEvents_double value     -->     "<< GA_numberOfEvents_double << G4endl;
                G4cout << "Here is the solidAngle value     -->     "<< solidAngle << G4endl;
                G4cout << "" << G4endl;
                */
                
                if(av_xPos==0)
                {
                    if(av_yPos==0) phi = 0;
                    if(av_yPos>0) phi = 90;
                    if(av_yPos<0) phi = 270;
                }
                else
                {
                    phi = atan(av_yPos/av_xPos)/deg;
                    
                    if(av_xPos>0 && av_yPos>0) phi = phi; // deg
                    if(av_xPos<0 && av_yPos>0) phi = phi + 180.; // deg
                    if(av_xPos<0 && av_yPos<0) phi = phi + 180.; // deg
                    if(av_xPos>0 && av_yPos<0) phi = phi + 360.; // deg
                }
                
                ////////////////////////////////////////////////////////////////////////////////////////
                file1 << TIARANo << ",              " << TIARA_RowNo << ",            " << TIARA_SectorNo << ",               " << theta << ",     " << phi <<  ",   " << solidAngle << endl;
                
                file2 << "    GA_TIARA[" << TIARANo << "][" << TIARA_RowNo << "][" << TIARA_SectorNo << "][0]=" << theta << ";   GA_TIARA[" << TIARANo << "][" << TIARA_RowNo << "][" << TIARA_SectorNo << "][1]=" << phi << ";   GA_TIARA[" << TIARANo << "][" << TIARA_RowNo << "][" << TIARA_SectorNo << "][2]=" << solidAngle << ";" << endl;
                ////////////////////////////////////////////////////////////////////////////////////////
                
            }
            
            file2 << "}"<< endl;
            
            file1.close();
            file2.close();
            

        }
    }
    
    

    
    
    ///////////////////////////////////////////////////////
    //          PRINTING the VDC Values
    /*
     for (G4int k=0; k<hit_buffersize; k++)
     {
     if(VDC_Observables[0][k]>=0 && VDC_Observables[1][k]>20)
     {
     G4cout << "                         " << G4endl;
     G4cout << "Here is the k (index) value     -->     "<< k << G4endl;
     G4cout << "Here is the VDC_Observables[0][k], Channel ID     -->     "<< VDC_Observables[0][k] << G4endl;
     G4cout << "Here is the VDC_Observables[1][k], Edep     -->     "<< VDC_Observables[1][k] << G4endl;
     G4cout << "Here is the VDC_Observables[2][k], EW_zpos     -->     "<< VDC_Observables[2][k] << G4endl;
     G4cout << "Here is the VDC_Observables[3][k], EW_t     -->     "<< VDC_Observables[3][k] << G4endl;
     G4cout << "Here is the VDC_Observables[2][k], zpos     -->     "<< VDC_Observables[2][k]/VDC_Observables[1][k] << G4endl;
     G4cout << "Here is the VDC_Observables[3][k], t     -->     "<< VDC_Observables[3][k]/VDC_Observables[1][k] << G4endl;
     }
     }
     */
    
    /*
     for (G4int k=0; k<hit_buffersize; k++)
     {
     if(VDC_Observables[1][k]>20 && (VDC_Observables[0][k]>=0) && (VDC_Observables[0][k]<=142))
     {
     G4cout << "Here is the VDC_Observables[0][k] value     -->     "<< VDC_Observables[0][k] << G4endl;
     
     }
     }
     G4cout << "Here is the Upos[0] value     -->     "<< Upos[0] << G4endl;
     */
    
    
    
    /*
  // fill histograms
  analysisManager->FillH1(1, fEnergyAbs);
  analysisManager->FillH1(2, fEnergyGap);
  analysisManager->FillH1(3, fTrackLAbs);
  analysisManager->FillH1(4, fTrackLGap);
  */
    
  // fill ntuple
    /*
  analysisManager->FillNtupleDColumn(0, fEnergyAbs);
  analysisManager->FillNtupleDColumn(1, fEnergyGap);
  analysisManager->FillNtupleDColumn(2, fTrackLAbs);
  analysisManager->FillNtupleDColumn(3, fTrackLGap);
  analysisManager->AddNtupleRow();  
  */
    
    /*
  // Print per event (modulo n)
  //
  G4int eventID = event->GetEventID();
  G4int printModulo = G4RunManager::GetRunManager()->GetPrintProgress();
  if ( ( printModulo > 0 ) && ( eventID % printModulo == 0 ) ) {
    G4cout << "---> End of event: " << eventID << G4endl;     

    G4cout
       << "   Absorber: total energy: " << std::setw(7)
                                        << G4BestUnit(fEnergyAbs,"Energy")
       << "       total track length: " << std::setw(7)
                                        << G4BestUnit(fTrackLAbs,"Length")
       << G4endl
       << "        Gap: total energy: " << std::setw(7)
                                        << G4BestUnit(fEnergyGap,"Energy")
       << "       total track length: " << std::setw(7)
                                        << G4BestUnit(fTrackLGap,"Length")
       << G4endl;
  }
    */
    
}  

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
