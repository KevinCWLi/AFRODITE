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

#include "RunAction.hh"
#include "Analysis.hh"

#include "G4Run.hh"
#include "G4RunManager.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RunAction::RunAction()
 : G4UserRunAction()
{ 
  // set printing event number per each event
  G4RunManager::GetRunManager()->SetPrintProgress(1);     

  // Create analysis manager
  // The choice of analysis technology is done via selectin of a namespace
  // in Analysis.hh
  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
  G4cout << "Using " << analysisManager->GetType() << G4endl;

  // Create directories 
  //analysisManager->SetHistoDirectoryName("histograms");
  //analysisManager->SetNtupleDirectoryName("ntuple");
  analysisManager->SetVerboseLevel(1);
  analysisManager->SetFirstHistoId(1);

  // Book histograms, ntuple
  //
  
    // Creating 1D histograms
    
    ////    TIARA DETECTORS, Histograms 1-6
    analysisManager->CreateH1("CvsE_TIARA1","TIARA 1 - Counts versus Energy", 10000, 0., 10000.);
    analysisManager->CreateH1("CvsE_TIARA2","TIARA 2 - Counts versus Energy", 10000, 0., 10000.);
    analysisManager->CreateH1("CvsE_TIARA3","TIARA 3 - Counts versus Energy", 10000, 0., 10000.);
    analysisManager->CreateH1("CvsE_TIARA4","TIARA 4 - Counts versus Energy", 10000, 0., 10000.);
    analysisManager->CreateH1("CvsE_TIARA5","TIARA 5 - Counts versus Energy", 10000, 0., 10000.);
    analysisManager->CreateH1("CvsE_TIARA_EA","Entire TIARA Array - Counts versus Energy", 10000, 0., 10000.);
    
    ////    PADDLE DETECTORS, Histograms 7-9
    analysisManager->CreateH1("CvsE_PADDLE1","PADDLE 1 - Counts versus Energy", 1000, 0., 100.);
    analysisManager->CreateH1("CvsE_PADDLE2","PADDLE 2 - Counts versus Energy", 1000, 0., 100.);
    analysisManager->CreateH1("CvsE_PADDLE3","PADDLE 3 - Counts versus Energy", 1000, 0., 100.);
    
    ////    CLOVER DETECTORS, Histograms 10-19
    analysisManager->CreateH1("CvsE_CLOVER1","CLOVER 1 - Counts versus Energy", 3000, 0., 3000.);
    analysisManager->CreateH1("CvsE_CLOVER2","CLOVER 2 - Counts versus Energy", 3000, 0., 3000.);
    analysisManager->CreateH1("CvsE_CLOVER3","CLOVER 3 - Counts versus Energy", 3000, 0., 3000.);
    analysisManager->CreateH1("CvsE_CLOVER4","CLOVER 4 - Counts versus Energy", 3000, 0., 3000.);
    analysisManager->CreateH1("CvsE_CLOVER5","CLOVER 5 - Counts versus Energy", 3000, 0., 3000.);
    analysisManager->CreateH1("CvsE_CLOVER6","CLOVER 6 - Counts versus Energy", 3000, 0., 3000.);
    analysisManager->CreateH1("CvsE_CLOVER7","CLOVER 7 - Counts versus Energy", 3000, 0., 3000.);
    analysisManager->CreateH1("CvsE_CLOVER8","CLOVER 8 - Counts versus Energy", 3000, 0., 3000.);
    analysisManager->CreateH1("CvsE_CLOVER9","CLOVER 9 - Counts versus Energy", 3000, 0., 3000.);
    analysisManager->CreateH1("CvsE_CLOVER_EA","Entire CLOVER Array - Counts versus Energy", 3000, 0., 3000.);

    
    // Creating 2D histograms
    
    ////    PADDLE DETECTORS, Histograms 1-6
    analysisManager->CreateH2("PvsE_PADDLE1","PADDLE 1 - Position vs. Energy", 1400, -700., 700., 3, -3*(102/2), 3*(102/2));
    analysisManager->CreateH2("PvsE_PADDLE2","PADDLE 2 - Position vs. Energy", 1400, -700., 700., 3, -3*(102/2), 3*(102/2));
    analysisManager->CreateH2("PvsE_PADDLE3","PADDLE 3 - Position vs. Energy", 1400, -700., 700., 3, -3*(102/2), 3*(102/2));
    analysisManager->CreateH2("EvsTOF_PADDLE1","PADDLE 1 - Energy vs. Time of Flight", 300, 0., 500.*0.3, 40, 0., 40.);
    analysisManager->CreateH2("EvsTOF_PADDLE2","PADDLE 2 - Energy vs. Time of Flight", 300, 0., 500.*0.3, 40, 0., 40.);
    analysisManager->CreateH2("EvsTOF_PADDLE3","PADDLE 3 - Energy vs. Time of Flight", 300, 0., 500.*0.3, 40, 0., 40.);
    


    ////////////////////////////////////////////////////
    //                  DataTreeSim
    ////////////////////////////////////////////////////
    
    // Creating ntuple
    analysisManager->CreateNtuple("DataTreeSim", "K600 Spectrometer - Coincident Events");
    
    ////    VDC Detectors
    analysisManager->CreateNtupleDColumn(0, "VDC1_Xpos");
    analysisManager->CreateNtupleDColumn(0, "VDC1_Y");
    analysisManager->CreateNtupleDColumn(0, "VDC1_ThetaFP");
    analysisManager->CreateNtupleDColumn(0, "VDC1_ThetaSCAT");

    
    analysisManager->CreateNtupleDColumn(0, "VDC2_Xpos");
    analysisManager->CreateNtupleDColumn(0, "VDC2_Y");
    analysisManager->CreateNtupleDColumn(0, "VDC2_ThetaFP");
    analysisManager->CreateNtupleDColumn(0, "VDC2_ThetaSCAT");
    
    ////    TIARA Detectors
    analysisManager->CreateNtupleIColumn(0, "TIARANo[1]");
    analysisManager->CreateNtupleIColumn(0, "TIARA_RowNo[1]");
    analysisManager->CreateNtupleIColumn(0, "TIARA_SectorNo[1]");
    
    analysisManager->CreateNtupleDColumn(0, "TIARA_E[1]");
    analysisManager->CreateNtupleDColumn(0, "Theta[1]");
    analysisManager->CreateNtupleDColumn(0, "Phi[1]");

    //analysisManager->CreateNtupleDColumn("TIARANo[2]");
    //analysisManager->CreateNtupleDColumn("TIARA_RowNo[2]");
    //analysisManager->CreateNtupleDColumn("TIARA_SectorNo[2]");
    //analysisManager->CreateNtupleDColumn(0, "TIARA_t[1]");

    /*
    ////    PADDLE Detectors
    analysisManager->CreateNtupleDColumn(0, "pad1_E");
    analysisManager->CreateNtupleDColumn(0, "pad1_tof");
    analysisManager->CreateNtupleDColumn(0, "pad2_E");
    analysisManager->CreateNtupleDColumn(0, "pad2_tof");
    
    ////    HAGAR Detector
    analysisManager->CreateNtupleDColumn("HAGAR_E");
    //analysisManager->CreateNtupleDColumn("HAGAR_t");
    */
    
    analysisManager->FinishNtuple(0);
    
    
    ////////////////////////////////////////////////////////////
    //                  GeometryAnalysisTree
    ////////////////////////////////////////////////////////////
    
    // Creating ntuple
    analysisManager->CreateNtuple("GeometryAnalysisTree", "K600 Spectrometer - GeometryAnalysis");
    
    ////    TIARA Detectors
    analysisManager->CreateNtupleIColumn(1, "TIARANo");
    analysisManager->CreateNtupleIColumn(1, "TIARA_RowNo");
    analysisManager->CreateNtupleIColumn(1, "TIARA_SectorNo");
    analysisManager->CreateNtupleDColumn(1, "Theta");
    analysisManager->CreateNtupleDColumn(1, "Phi");

    analysisManager->FinishNtuple(1);
    
    
    ////////////////////////////////////////////////////////////
    //                  Input Variable Tree
    ////////////////////////////////////////////////////////////
    
    // Creating ntuple
    analysisManager->CreateNtuple("InputVariableTree", "K600 Spectrometer - InputVariable");
    
    ////    Initial Particle Angular Distributions
    analysisManager->CreateNtupleDColumn(2, "ThetaDist");
    analysisManager->CreateNtupleDColumn(2, "PhiDist");
    
    analysisManager->FinishNtuple(2);
    
    
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RunAction::~RunAction()
{
  delete G4AnalysisManager::Instance();  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::BeginOfRunAction(const G4Run* /*run*/)
{ 
  //inform the runManager to save random number seed
  //G4RunManager::GetRunManager()->SetRandomNumberStore(true);
  
  // Get analysis manager
  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();

  // Open an output file
  //
  G4String fileName = "K600Output";
  analysisManager->OpenFile(fileName);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::EndOfRunAction(const G4Run* /*run*/)
{
    
  // print histogram statistics
  //
  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
  if ( analysisManager->GetH1(1) ) {
      /*
    G4cout << "\n ----> print histograms statistic ";
    if(isMaster) {
      G4cout << "for the entire run \n" << G4endl; 
    }
    else {
      G4cout << "for the local thread \n" << G4endl; 
    }
    
    G4cout << " EAbs : mean = " 
       << G4BestUnit(analysisManager->GetH1(1)->mean(), "Energy") 
       << " rms = " 
       << G4BestUnit(analysisManager->GetH1(1)->rms(),  "Energy") << G4endl;
    
    G4cout << " EGap : mean = " 
       << G4BestUnit(analysisManager->GetH1(2)->mean(), "Energy") 
       << " rms = " 
       << G4BestUnit(analysisManager->GetH1(2)->rms(),  "Energy") << G4endl;
    
    G4cout << " LAbs : mean = " 
      << G4BestUnit(analysisManager->GetH1(3)->mean(), "Length") 
      << " rms = " 
      << G4BestUnit(analysisManager->GetH1(3)->rms(),  "Length") << G4endl;

    G4cout << " LGap : mean = " 
      << G4BestUnit(analysisManager->GetH1(4)->mean(), "Length") 
      << " rms = " 
      << G4BestUnit(analysisManager->GetH1(4)->rms(),  "Length") << G4endl;
       */
  }

  // save histograms & ntuple
  //
  analysisManager->Write();
  analysisManager->CloseFile();

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
