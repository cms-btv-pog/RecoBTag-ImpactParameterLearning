//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Thu Jul  4 14:09:09 2013 by ROOT version 5.32/00
// from TTree ttree/ttree
// found on file: JetTree_mc.root
//////////////////////////////////////////////////////////

#ifndef JetProbaValidation_h
#define JetProbaValidation_h

#include "CategoryDef.h"
#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <iostream>
#include <fstream>
#include <TH1.h>
#include "TList.h"
#include "CategoryDef.h"
#include "CategoriesDefinition.h"
#include "CategoryDefCollection.h"


// Header file for the classes stored in the TTree if any.

// Fixed size dimensions of array or collections stored in the TTree if any.

class JetProbaValidation {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
   Int_t           nBitTrigger;
   Int_t           BitTrigger[3];   //[nBitTrigger]
   Int_t           Run;
   Int_t           Evt;
   Int_t           LumiBlock;
   Float_t         pthat;
   Float_t         mcweight;
   Int_t           nPV;
   Float_t         PVz;
   Float_t         nPUtrue;
   Int_t           nPU;
   
   Int_t           nJet;
   Float_t         Jet_pt[200];   //[nJet]
   Float_t         Jet_genpt[200];   //[nJet]
   Float_t         Jet_residual[200];   //[nJet]
   Float_t         Jet_jes[200];   //[nJet]
   Float_t         Jet_eta[200];   //[nJet]
   Float_t         Jet_phi[200];   //[nJet]
   Float_t         Jet_mass[200];   //[nJet]
   Int_t           Jet_ntracks[200];   //[nJet]
   Int_t           Jet_nseltracks[200];   //[nJet]
   Int_t           Jet_flavour[200];   //[nJet]
   Int_t           Jet_nbHadrons[200];   //[nJet]
   Int_t           Jet_ncHadrons[200];   //[nJet]
   Float_t         Jet_Ip2N[200];   //[nJet]
   Float_t         Jet_Ip2P[200];   //[nJet]
   Float_t         Jet_Ip3N[200];   //[nJet]
   Float_t         Jet_Ip3P[200];   //[nJet]
   Float_t         Jet_ProbaN[200];   //[nJet]
   Float_t         Jet_ProbaP[200];   //[nJet]
   Float_t         Jet_Proba[200];   //[nJet]
   Float_t         Jet_BprobN[200];   //[nJet]
   Float_t         Jet_BprobP[200];   //[nJet]
   Float_t         Jet_Bprob[200];   //[nJet]
   Float_t         Jet_SvxN[200];   //[nJet]
   Float_t         Jet_Svx[200];   //[nJet]
   Float_t         Jet_SvxNHP[200];   //[nJet]
   Float_t         Jet_SvxHP[200];   //[nJet]
   Float_t         Jet_CombSvxN[200];   //[nJet]
   Float_t         Jet_CombSvxP[200];   //[nJet]
   Float_t         Jet_CombSvx[200];   //[nJet]
   Float_t         Jet_CombIVF[200];   //[nJet]
   Float_t         Jet_CombIVF_P[200];   //[nJet]
   Float_t         Jet_CombIVF_N[200];   //[nJet]
   Float_t         Jet_SoftMuN[200];   //[nJet]
   Float_t         Jet_SoftMuP[200];   //[nJet]
   Float_t         Jet_SoftMu[200];   //[nJet]
   Float_t         Jet_SoftElN[200];   //[nJet]
   Float_t         Jet_SoftElP[200];   //[nJet]
   Float_t         Jet_SoftEl[200];   //[nJet]
   Float_t         Jet_DoubleSV[200];   //[nJet]
   Int_t           Jet_hist1[200];   //[nJet]
   Int_t           Jet_hist2[200];   //[nJet]
   Int_t           Jet_hist3[200];   //[nJet]
   Int_t           Jet_histJet[200];   //[nJet]
   Int_t           Jet_histSvx[200];   //[nJet]
   Int_t           Jet_nFirstTrack[200];   //[nJet]
   Int_t           Jet_nLastTrack[200];   //[nJet]
   Int_t           Jet_nFirstSV[200];   //[nJet]
   Int_t           Jet_nLastSV[200];   //[nJet]
   Int_t           Jet_SV_multi[200];   //[nJet]
   Int_t           Jet_nFirstTrkInc[200];   //[nJet]
   Int_t           Jet_nLastTrkInc[200];   //[nJet]
   Int_t           Jet_nSM[200];   //[nJet]
   Int_t           Jet_nFirstSM[200];   //[nJet]
   Int_t           Jet_nLastSM[200];   //[nJet]
   Int_t           Jet_nSE[200];   //[nJet]
   Int_t           Jet_nFirstSE[200];   //[nJet]
   Int_t           Jet_nLastSE[200];   //[nJet]
   Int_t           Jet_looseID[200];   //[nJet]
   Int_t           Jet_tightID[200];   //[nJet]
   Int_t           nTrack;
   Float_t         Track_dxy[2000];   //[nTrack]
   Float_t         Track_dz[2000];   //[nTrack]
   Float_t         Track_length[2000];   //[nTrack]
   Float_t         Track_dist[2000];   //[nTrack]
   Float_t         Track_IP2D[2000];   //[nTrack]
   Float_t         Track_IP2Dsig[2000];   //[nTrack]
   Float_t         Track_IP2Derr[2000];   //[nTrack]
   Float_t         Track_IP[2000];   //[nTrack]
   Float_t         Track_IPsig[2000];   //[nTrack]
   Float_t         Track_IPerr[2000];   //[nTrack]
   Float_t         Track_Proba[2000];   //[nTrack]
   Float_t         Track_p[2000];   //[nTrack]
   Float_t         Track_pt[2000];   //[nTrack]
   Float_t         Track_eta[2000];   //[nTrack]
   Float_t         Track_phi[2000];   //[nTrack]
   Float_t         Track_chi2[2000];   //[nTrack]
   Int_t           Track_charge[2000];   //[nTrack]
   Int_t           Track_history[2000];   //[nTrack]
   Int_t           Track_nHitStrip[2000];   //[nTrack]
   Int_t           Track_nHitPixel[2000];   //[nTrack]
   Int_t           Track_nHitAll[2000];   //[nTrack]
   Int_t           Track_nHitTIB[2000];   //[nTrack]
   Int_t           Track_nHitTID[2000];   //[nTrack]
   Int_t           Track_nHitTOB[2000];   //[nTrack]
   Int_t           Track_nHitTEC[2000];   //[nTrack]
   Int_t           Track_nHitPXB[2000];   //[nTrack]
   Int_t           Track_nHitPXF[2000];   //[nTrack]
   Int_t           Track_isHitL1[2000];   //[nTrack]
   Int_t           Track_PV[2000];   //[nTrack]
   Int_t           Track_SV[2000];   //[nTrack]
   Float_t         Track_PVweight[2000];   //[nTrack]
   Float_t         Track_SVweight[2000];   //[nTrack]
   Int_t           Track_isfromSV[2000];   //[nTrack]
   Int_t           Track_category[2000];   //[nTrack]
   
   
   int number_of_category;
   int nbins;
   float xmin;
   float xmax;
      
   // List of branches
   TBranch        *b_nBitTrigger;   //!
   TBranch        *b_BitTrigger;   //!
   TBranch        *b_Run;   //!
   TBranch        *b_Evt;   //!
   TBranch        *b_LumiBlock;   //!
   TBranch        *b_pthat;   //!
   TBranch        *b_mcweight;   //!
   TBranch        *b_nPV;   //!
   TBranch        *b_PVz;   //!
   TBranch        *b_nPUtrue;   //!
   TBranch        *b_nPU;   //!
   TBranch        *b_nJet;   //!
   TBranch        *b_Jet_pt;   //!
   TBranch        *b_Jet_genpt;   //!
   TBranch        *b_Jet_residual;   //!
   TBranch        *b_Jet_jes;   //!
   TBranch        *b_Jet_eta;   //!
   TBranch        *b_Jet_phi;   //!
   TBranch        *b_Jet_mass;   //!
   TBranch        *b_Jet_ntracks;   //!
   TBranch        *b_Jet_nseltracks;   //!
   TBranch        *b_Jet_flavour;   //!
   TBranch        *b_Jet_nbHadrons;   //!
   TBranch        *b_Jet_ncHadrons;   //!
   TBranch        *b_Jet_Ip2N;   //!
   TBranch        *b_Jet_Ip2P;   //!
   TBranch        *b_Jet_Ip3N;   //!
   TBranch        *b_Jet_Ip3P;   //!
   TBranch        *b_Jet_ProbaN;   //!
   TBranch        *b_Jet_ProbaP;   //!
   TBranch        *b_Jet_Proba;   //!
   TBranch        *b_Jet_BprobN;   //!
   TBranch        *b_Jet_BprobP;   //!
   TBranch        *b_Jet_Bprob;   //!
   TBranch        *b_Jet_SvxN;   //!
   TBranch        *b_Jet_Svx;   //!
   TBranch        *b_Jet_SvxNHP;   //!
   TBranch        *b_Jet_SvxHP;   //!
   TBranch        *b_Jet_CombSvxN;   //!
   TBranch        *b_Jet_CombSvxP;   //!
   TBranch        *b_Jet_CombSvx;   //!
   TBranch        *b_Jet_CombIVF;   //!
   TBranch        *b_Jet_CombIVF_P;   //!
   TBranch        *b_Jet_CombIVF_N;   //!
   TBranch        *b_Jet_SoftMuN;   //!
   TBranch        *b_Jet_SoftMuP;   //!
   TBranch        *b_Jet_SoftMu;   //!
   TBranch        *b_Jet_SoftElN;   //!
   TBranch        *b_Jet_SoftElP;   //!
   TBranch        *b_Jet_SoftEl;   //!
   TBranch        *b_Jet_DoubleSV;   //!
   TBranch        *b_Jet_hist1;   //!
   TBranch        *b_Jet_hist2;   //!
   TBranch        *b_Jet_hist3;   //!
   TBranch        *b_Jet_histJet;   //!
   TBranch        *b_Jet_histSvx;   //!
   TBranch        *b_Jet_nFirstTrack;   //!
   TBranch        *b_Jet_nLastTrack;   //!
   TBranch        *b_Jet_nFirstSV;   //!
   TBranch        *b_Jet_nLastSV;   //!
   TBranch        *b_Jet_SV_multi;   //!
   TBranch        *b_Jet_nFirstTrkInc;   //!
   TBranch        *b_Jet_nLastTrkInc;   //!
   TBranch        *b_Jet_nSM;   //!
   TBranch        *b_Jet_nFirstSM;   //!
   TBranch        *b_Jet_nLastSM;   //!
   TBranch        *b_Jet_nSE;   //!
   TBranch        *b_Jet_nFirstSE;   //!
   TBranch        *b_Jet_nLastSE;   //!
   TBranch        *b_Jet_looseID;   //!
   TBranch        *b_Jet_tightID;   //!
   TBranch        *b_nTrack;   //!
   TBranch        *b_Track_dxy;   //!
   TBranch        *b_Track_dz;   //!
   TBranch        *b_Track_length;   //!
   TBranch        *b_Track_dist;   //!
   TBranch        *b_Track_IP2D;   //!
   TBranch        *b_Track_IP2Dsig;   //!
   TBranch        *b_Track_IP2Derr;   //!
   TBranch        *b_Track_IP;   //!
   TBranch        *b_Track_IPsig;   //!
   TBranch        *b_Track_IPerr;   //!
   TBranch        *b_Track_Proba;   //!
   TBranch        *b_Track_p;   //!
   TBranch        *b_Track_pt;   //!
   TBranch        *b_Track_eta;   //!
   TBranch        *b_Track_phi;   //!
   TBranch        *b_Track_chi2;   //!
   TBranch        *b_Track_charge;   //!
   TBranch        *b_Track_history;   //!
   TBranch        *b_Track_nHitStrip;   //!
   TBranch        *b_Track_nHitPixel;   //!
   TBranch        *b_Track_nHitAll;   //!
   TBranch        *b_Track_nHitTIB;   //!
   TBranch        *b_Track_nHitTID;   //!
   TBranch        *b_Track_nHitTOB;   //!
   TBranch        *b_Track_nHitTEC;   //!
   TBranch        *b_Track_nHitPXB;   //!
   TBranch        *b_Track_nHitPXF;   //!
   TBranch        *b_Track_isHitL1;   //!
   TBranch        *b_Track_PV;   //!
   TBranch        *b_Track_SV;   //!
   TBranch        *b_Track_PVweight;   //!
   TBranch        *b_Track_SVweight;   //!
   TBranch        *b_Track_isfromSV;   //!
   TBranch        *b_Track_category;   //!

   JetProbaValidation(TTree *tree=0);
   virtual ~JetProbaValidation();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
   
   
   
   virtual double jetProbability( std::vector<double>  v) ;
   virtual int    IsInCategory(float trkEta, float trkHTrk, float trkHPix, float trkp, 
   float trkChi2, std::vector<CategoryDef > d);
   virtual void    ComputeProba(const char*fileInPutName);
   virtual double  calculTrackProba( float ipsig, CategoryDef vectCat  );
   virtual bool    passTrackSel(int trk, int pix, float ip2d, float pt, float len, float chi2, float zip, float dist, float eta, float etaJet, float phi, float phiJet, float ptJet);
   virtual void    plotPerf();



};

#endif

#ifdef JetProbaValidation_cxx
JetProbaValidation::JetProbaValidation(TTree *tree) : fChain(0) 
{  
   
   TChain *superTree = new TChain("btagana/ttree");
   if (tree == 0) {
  
   //  superTree->Add("/opt/sbg/data/data2/cms/cbeluffi/BTag_2013_06_24/DataLegacy/BTagAna/CMSSW_5_3_11_patch5/src/RecoBTag/PerformanceMeasurements/test/mc_calib_confirm2/JetTree_*.root");
     
   Init(superTree);
   } else Init(tree);
   
}

JetProbaValidation::~JetProbaValidation()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t JetProbaValidation::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t JetProbaValidation::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void JetProbaValidation::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("nBitTrigger", &nBitTrigger, &b_nBitTrigger);
   fChain->SetBranchAddress("BitTrigger", BitTrigger, &b_BitTrigger);
   fChain->SetBranchAddress("Run", &Run, &b_Run);
   fChain->SetBranchAddress("Evt", &Evt, &b_Evt);
   fChain->SetBranchAddress("LumiBlock", &LumiBlock, &b_LumiBlock);
   fChain->SetBranchAddress("pthat", &pthat, &b_pthat);
   fChain->SetBranchAddress("mcweight", &mcweight, &b_mcweight);
   fChain->SetBranchAddress("nPV", &nPV, &b_nPV);
   fChain->SetBranchAddress("PVz", &PVz, &b_PVz);
   fChain->SetBranchAddress("nPUtrue", &nPUtrue, &b_nPUtrue);
   fChain->SetBranchAddress("nPU", &nPU, &b_nPU);
   fChain->SetBranchAddress("nJet", &nJet, &b_nJet);
   fChain->SetBranchAddress("Jet_pt", Jet_pt, &b_Jet_pt);
   fChain->SetBranchAddress("Jet_genpt", Jet_genpt, &b_Jet_genpt);
   fChain->SetBranchAddress("Jet_residual", Jet_residual, &b_Jet_residual);
   fChain->SetBranchAddress("Jet_jes", Jet_jes, &b_Jet_jes);
   fChain->SetBranchAddress("Jet_eta", Jet_eta, &b_Jet_eta);
   fChain->SetBranchAddress("Jet_phi", Jet_phi, &b_Jet_phi);
   fChain->SetBranchAddress("Jet_mass", Jet_mass, &b_Jet_mass);
   fChain->SetBranchAddress("Jet_ntracks", Jet_ntracks, &b_Jet_ntracks);
   fChain->SetBranchAddress("Jet_nseltracks", Jet_nseltracks, &b_Jet_nseltracks);
   fChain->SetBranchAddress("Jet_flavour", Jet_flavour, &b_Jet_flavour);
   fChain->SetBranchAddress("Jet_nbHadrons", Jet_nbHadrons, &b_Jet_nbHadrons);
   fChain->SetBranchAddress("Jet_ncHadrons", Jet_ncHadrons, &b_Jet_ncHadrons);
   fChain->SetBranchAddress("Jet_Ip2N", Jet_Ip2N, &b_Jet_Ip2N);
   fChain->SetBranchAddress("Jet_Ip2P", Jet_Ip2P, &b_Jet_Ip2P);
   fChain->SetBranchAddress("Jet_Ip3N", Jet_Ip3N, &b_Jet_Ip3N);
   fChain->SetBranchAddress("Jet_Ip3P", Jet_Ip3P, &b_Jet_Ip3P);
   fChain->SetBranchAddress("Jet_ProbaN", Jet_ProbaN, &b_Jet_ProbaN);
   fChain->SetBranchAddress("Jet_ProbaP", Jet_ProbaP, &b_Jet_ProbaP);
   fChain->SetBranchAddress("Jet_Proba", Jet_Proba, &b_Jet_Proba);
   fChain->SetBranchAddress("Jet_BprobN", Jet_BprobN, &b_Jet_BprobN);
   fChain->SetBranchAddress("Jet_BprobP", Jet_BprobP, &b_Jet_BprobP);
   fChain->SetBranchAddress("Jet_Bprob", Jet_Bprob, &b_Jet_Bprob);
   fChain->SetBranchAddress("Jet_SvxN", Jet_SvxN, &b_Jet_SvxN);
   fChain->SetBranchAddress("Jet_Svx", Jet_Svx, &b_Jet_Svx);
   fChain->SetBranchAddress("Jet_SvxNHP", Jet_SvxNHP, &b_Jet_SvxNHP);
   fChain->SetBranchAddress("Jet_SvxHP", Jet_SvxHP, &b_Jet_SvxHP);
   fChain->SetBranchAddress("Jet_CombSvxN", Jet_CombSvxN, &b_Jet_CombSvxN);
   fChain->SetBranchAddress("Jet_CombSvxP", Jet_CombSvxP, &b_Jet_CombSvxP);
   fChain->SetBranchAddress("Jet_CombSvx", Jet_CombSvx, &b_Jet_CombSvx);
   fChain->SetBranchAddress("Jet_CombIVF", Jet_CombIVF, &b_Jet_CombIVF);
   fChain->SetBranchAddress("Jet_CombIVF_P", Jet_CombIVF_P, &b_Jet_CombIVF_P);
   fChain->SetBranchAddress("Jet_CombIVF_N", Jet_CombIVF_N, &b_Jet_CombIVF_N);
   fChain->SetBranchAddress("Jet_SoftMuN", Jet_SoftMuN, &b_Jet_SoftMuN);
   fChain->SetBranchAddress("Jet_SoftMuP", Jet_SoftMuP, &b_Jet_SoftMuP);
   fChain->SetBranchAddress("Jet_SoftMu", Jet_SoftMu, &b_Jet_SoftMu);
   fChain->SetBranchAddress("Jet_SoftElN", Jet_SoftElN, &b_Jet_SoftElN);
   fChain->SetBranchAddress("Jet_SoftElP", Jet_SoftElP, &b_Jet_SoftElP);
   fChain->SetBranchAddress("Jet_SoftEl", Jet_SoftEl, &b_Jet_SoftEl);
   fChain->SetBranchAddress("Jet_DoubleSV", Jet_DoubleSV, &b_Jet_DoubleSV);
   fChain->SetBranchAddress("Jet_hist1", Jet_hist1, &b_Jet_hist1);
   fChain->SetBranchAddress("Jet_hist2", Jet_hist2, &b_Jet_hist2);
   fChain->SetBranchAddress("Jet_hist3", Jet_hist3, &b_Jet_hist3);
   fChain->SetBranchAddress("Jet_histJet", Jet_histJet, &b_Jet_histJet);
   fChain->SetBranchAddress("Jet_histSvx", Jet_histSvx, &b_Jet_histSvx);
   fChain->SetBranchAddress("Jet_nFirstTrack", Jet_nFirstTrack, &b_Jet_nFirstTrack);
   fChain->SetBranchAddress("Jet_nLastTrack", Jet_nLastTrack, &b_Jet_nLastTrack);
   fChain->SetBranchAddress("Jet_nFirstSV", Jet_nFirstSV, &b_Jet_nFirstSV);
   fChain->SetBranchAddress("Jet_nLastSV", Jet_nLastSV, &b_Jet_nLastSV);
   fChain->SetBranchAddress("Jet_SV_multi", Jet_SV_multi, &b_Jet_SV_multi);
   fChain->SetBranchAddress("Jet_nFirstTrkInc", Jet_nFirstTrkInc, &b_Jet_nFirstTrkInc);
   fChain->SetBranchAddress("Jet_nLastTrkInc", Jet_nLastTrkInc, &b_Jet_nLastTrkInc);
   fChain->SetBranchAddress("Jet_nSM", Jet_nSM, &b_Jet_nSM);
   fChain->SetBranchAddress("Jet_nFirstSM", Jet_nFirstSM, &b_Jet_nFirstSM);
   fChain->SetBranchAddress("Jet_nLastSM", Jet_nLastSM, &b_Jet_nLastSM);
   fChain->SetBranchAddress("Jet_nSE", Jet_nSE, &b_Jet_nSE);
   fChain->SetBranchAddress("Jet_nFirstSE", Jet_nFirstSE, &b_Jet_nFirstSE);
   fChain->SetBranchAddress("Jet_nLastSE", Jet_nLastSE, &b_Jet_nLastSE);
   fChain->SetBranchAddress("Jet_looseID", Jet_looseID, &b_Jet_looseID);
   fChain->SetBranchAddress("Jet_tightID", Jet_tightID, &b_Jet_tightID);
   fChain->SetBranchAddress("nTrack", &nTrack, &b_nTrack);
   fChain->SetBranchAddress("Track_dxy", Track_dxy, &b_Track_dxy);
   fChain->SetBranchAddress("Track_dz", Track_dz, &b_Track_dz);
   fChain->SetBranchAddress("Track_length", Track_length, &b_Track_length);
   fChain->SetBranchAddress("Track_dist", Track_dist, &b_Track_dist);
   fChain->SetBranchAddress("Track_IP2D", Track_IP2D, &b_Track_IP2D);
   fChain->SetBranchAddress("Track_IP2Dsig", Track_IP2Dsig, &b_Track_IP2Dsig);
   fChain->SetBranchAddress("Track_IP2Derr", Track_IP2Derr, &b_Track_IP2Derr);
   fChain->SetBranchAddress("Track_IP", Track_IP, &b_Track_IP);
   fChain->SetBranchAddress("Track_IPsig", Track_IPsig, &b_Track_IPsig);
   fChain->SetBranchAddress("Track_IPerr", Track_IPerr, &b_Track_IPerr);
   fChain->SetBranchAddress("Track_Proba", Track_Proba, &b_Track_Proba);
   fChain->SetBranchAddress("Track_p", Track_p, &b_Track_p);
   fChain->SetBranchAddress("Track_pt", Track_pt, &b_Track_pt);
   fChain->SetBranchAddress("Track_eta", Track_eta, &b_Track_eta);
   fChain->SetBranchAddress("Track_phi", Track_phi, &b_Track_phi);
   fChain->SetBranchAddress("Track_chi2", Track_chi2, &b_Track_chi2);
   fChain->SetBranchAddress("Track_charge", Track_charge, &b_Track_charge);
   fChain->SetBranchAddress("Track_history", Track_history, &b_Track_history);
   fChain->SetBranchAddress("Track_nHitStrip", Track_nHitStrip, &b_Track_nHitStrip);
   fChain->SetBranchAddress("Track_nHitPixel", Track_nHitPixel, &b_Track_nHitPixel);
   fChain->SetBranchAddress("Track_nHitAll", Track_nHitAll, &b_Track_nHitAll);
   fChain->SetBranchAddress("Track_nHitTIB", Track_nHitTIB, &b_Track_nHitTIB);
   fChain->SetBranchAddress("Track_nHitTID", Track_nHitTID, &b_Track_nHitTID);
   fChain->SetBranchAddress("Track_nHitTOB", Track_nHitTOB, &b_Track_nHitTOB);
   fChain->SetBranchAddress("Track_nHitTEC", Track_nHitTEC, &b_Track_nHitTEC);
   fChain->SetBranchAddress("Track_nHitPXB", Track_nHitPXB, &b_Track_nHitPXB);
   fChain->SetBranchAddress("Track_nHitPXF", Track_nHitPXF, &b_Track_nHitPXF);
   fChain->SetBranchAddress("Track_isHitL1", Track_isHitL1, &b_Track_isHitL1);
   fChain->SetBranchAddress("Track_PV", Track_PV, &b_Track_PV);
   fChain->SetBranchAddress("Track_SV", Track_SV, &b_Track_SV);
   fChain->SetBranchAddress("Track_PVweight", Track_PVweight, &b_Track_PVweight);
   fChain->SetBranchAddress("Track_SVweight", Track_SVweight, &b_Track_SVweight);
   fChain->SetBranchAddress("Track_isfromSV", Track_isfromSV, &b_Track_isfromSV);
   fChain->SetBranchAddress("Track_category", Track_category, &b_Track_category);
   Notify();
   
}

Bool_t JetProbaValidation::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void JetProbaValidation::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t JetProbaValidation::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef JetProbaValidation_cxx
