#define JetProbaValidation_cxx
#include "JetProbaValidation.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TH1.h>
#include "CategoryDef.h"
#include <vector>
#include "TClass.h"
#include "TMath.h"
#include "TLorentzVector.h"
#include <iomanip>
#include "TGraph.h"
#include "TPaveLabel.h"
#include "TPad.h"
#include "TGraphErrors.h"
#include "TLegend.h"

using namespace std;


//-----------------------------------------------------Find the category of a given track----------------------------------------------------------------------------------------//
int JetProbaValidation::IsInCategory(float trkEta, float trkHTrk, float trkHPix, float trkp, float trkChi2, std::vector<CategoryDef > d){
  
  int isIn=-1;
  
  for(unsigned int i=0; i<d.size(); i++){
   if (fabs(trkEta)>d[i].etaMin&&fabs(trkEta)<d[i].etaMax&&trkHTrk<=d[i].nHitsMax&&trkHTrk>=d[i].nHitsMin&&
     trkHPix<=d[i].nPixelHitsMax&&trkHPix>=d[i].nPixelHitsMin&&trkp<d[i].pMax&&trkp>d[i].pMin&&
     trkChi2>=d[i].chiMin&&trkChi2<=d[i].chiMax) isIn=i;
  }
  
  return isIn;
}


//-----------------------------------------------------Return true if track passes btag selection ---------------------------------------------------------------------------------//
bool JetProbaValidation::passTrackSel(int trk, int pix, float ip2d, float pt, float len, float chi2, float zip, float dist, float eta, float etaJet, float phi, float phiJet, float ptJet)
{
  TLorentzVector track4P, jet4P;
  track4P.SetPtEtaPhiM(pt, eta,phi , 0. );
  jet4P.SetPtEtaPhiM(ptJet ,etaJet ,phiJet , 0 );

  float deltaR = jet4P.DeltaR(track4P);
  bool passSel=false;
  
  bool passTrk=false;
  bool passIP=false;
  bool passPix=false;
  bool passPt=false;
  bool passLen=false;
  bool passChi2=false;
  bool passDist=false;
  bool passZip=false;
  bool passDR=false;
  
  //if ((trk+pix)>=8)   passTrk=true;
  if ((trk)>=8)   passTrk=true;
  
  if (fabs(ip2d)<0.2) passIP=true;
  if (pix>=2)         passPix=true;
  if (pt>1.0)         passPt=true;
  if (len<5.0)        passLen=true;
  if (chi2<5.0)       passChi2=true;
  if (fabs(dist)<0.07)passDist=true;
  if (fabs(zip)<17.0) passZip=true;
  if (deltaR<0.3)     passDR=true;
  
  if (passTrk&&passIP&&passPix&&passPt&&passLen&&passChi2&&passDist&&passZip&&passDR) {
  
    passSel=true;
  }
  
  return passSel;
}
  
  
//-----------------------------------------------------Compute the probability that a track comes from PV ---------------------------------------------------------------------------//
double JetProbaValidation::calculTrackProba( float ipsig, CategoryDef theCat){
  
  //int catnum = IsInCategory(eta,nhit ,npix ,p, chi, vectCat);
  
  double proba= -1000.0;
  double histoIntegral = 1.0;
  
  TH1D * thecathisto = theCat.histo;
  
  double nbins=thecathisto->GetNbinsX();
//   double xmin =thecathisto->GetXaxis()->GetXmin();
//   double xmax =thecathisto->GetXaxis()->GetXmax();
  
  int nbin = int(fabs(ipsig)*(nbins+1)/xmax);
  
  double sum=0.0;
  for(int i=nbin; i<=(nbins) ; i++){
    double x=thecathisto->GetBinContent(i);
    sum=sum+x;
  }
  
  double sum_tot=0.0;
  for(int i=1; i<=(nbins) ; i++){
    //double x=myHisto->GetBinContent(i);
    double x=thecathisto->GetBinContent(i);
    sum_tot=sum_tot+x;
  }    
  proba=sum/sum_tot;
  
  //if(proba < 0.00001)   cout << "found a proba of " << proba << "  ip/sigma "  << ipsig << "  xmax " << xmax << endl;
  return proba;
}


double JetProbaValidation::jetProbability( std::vector<double>  v) 
{
  std::vector<double> probabilities;
  int i=0;

  for(std::vector<double>::const_iterator it = v.begin(); it!=v.end(); ++it, i++)
  {
      float p;
      //if(*it >=0 ) p=*it; else continue; 
      if (*it >=0){p=*it/2.;}else{p=1.-*it/2.;}
      //p=*it ; 
      probabilities.push_back(p);
  } 
 
  double m_minTrackProb= 0.005;
  int ngoodtracks=v.size();
  double SumJet=0.;

  for(std::vector<double>::const_iterator q = probabilities.begin(); q != probabilities.end(); q++){
    SumJet+=(*q>m_minTrackProb)?log(*q):log(m_minTrackProb);
  }

  double ProbJet;
  double Loginvlog=0;

  if(SumJet<0.){
    if(ngoodtracks>=2){
      Loginvlog=log(-SumJet);
    }
    double Prob=1.;
    double lfact=1.;
    for(int l=1; l!=ngoodtracks; l++){
       lfact*=l;
      Prob+=exp(l*Loginvlog-log(1.*lfact));
    }
    double LogProb=log(Prob);
    ProbJet=
      std::min(exp(std::max(LogProb+SumJet,-30.)),1.);
  }else{
    ProbJet=1.;
  }
  if(ProbJet>1)
   std::cout << "ProbJet too high: "  << ProbJet << std::endl;

   return -log10(ProbJet)/4.;
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------//
//--------------------------------------------------SECOND LOOP TO FILL COMPUTE TRACK PROBABILTY-----------------------------------------------------------------------------------//
//------------------------------------------------------------------------------------------------------------------------------------------------------------------------------//

void JetProbaValidation::ComputeProba(const char*fileInPutName)
{ 
  cout << "initialize category files " << endl;
  TFile* myinputfile = new TFile(fileInPutName );
  CategoryDefCollection * pcatDefVector = (CategoryDefCollection*) myinputfile->Get("CategoryDefCollection");
  
  cout << "category files initialized" << endl;
  if(pcatDefVector == 0) cout << "null pointer " << endl;
  std::vector<CategoryDef > catDefList = pcatDefVector->ListOfCategory;

  for(unsigned int i=0; i<catDefList.size(); i++){
  
    catDefList[i].GetDistribFromFile(myinputfile );
  }

  TFile* myfileout = new TFile("JP_myCalib.root", "recreate");
  myfileout->cd();
  
//$$
  int nbin_trkH=52;
  float xmin_trkH=-0.02;
  float xmax_trkH=1.02;  
//$$

  cout << "initialize control histograms " << endl;
  TH1D * JetProba = new TH1D("JetProba", "JetProba", 100, 0.0, 2.5);
  
  TH1D * trackPCat_all = new TH1D("trackPCat_all", "trackPCat_all", nbin_trkH, xmin_trkH, xmax_trkH);

  std::vector<TH1D *> vectTrkProba; 
  
  for(unsigned int i = 0; i<catDefList.size(); i++){
    TH1D* tmphisto = new TH1D( (catDefList[i].histoname+"proba").Data(), (catDefList[i].histoname+"proba").Data()  , nbin_trkH, xmin_trkH, xmax_trkH);
    vectTrkProba.push_back(tmphisto);  
  }
  
  TH1F * jetProbaDiscriB = new TH1F("jetProbaDiscriB", "jetProbaDiscriB", 101, 0, 1); 
  TH1F * jetProbaDiscriUDSG = new TH1F("jetProbaDiscriUDSG", "jetProbaDiscriUDSG", 101, 0, 1);
  
  std::vector<double > vectProba;
  //std::vector<double>  vectProba_trueMC ;
  
  Long64_t nentries = fChain->GetEntriesFast();
  
//$$
//  int itest50 = 0, itest80 = 0, itest120 = 0, itest170 = 0;
// // // 50ns:
// //  int test50  = int(1920430./1938260. * 337747.); 
// //  int test80  = int(276253./278517. * 222099.); 
// //  int test120 = int(47110./49997. * 49309.); 
// //  int test170 = int(11727./53976. * 53936.);
// // 25ns:
//  int test50  = int(1920430./2359910. * 405733.); 
//  int test80  = int(276253./408702. * 314011.); 
//  int test120 = int(47110./99196. * 97565.); 
//  int test170 = int(11727./99648. * 99577.);
//  std::cout << std::endl; 
//  std::cout << " test50-80-120-170 " 
//  << test50 << " " << test80 << " " << test120 << " " << test170 << " " << std::endl;
//  std::cout << std::endl; 
//$$

//*********************************
//*********************************
// Event loop

  Long64_t nbytes = 0, nb = 0;  
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0 ) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;
    
//*********************************
// pthat reweighting
//$$
//     if ( pthat >= 50 && pthat < 80 ) itest50++; 
//     if ( pthat >= 80 && pthat <120 ) itest80++; 
//     if ( pthat >=120 && pthat <170 ) itest120++; 
//     if ( pthat >=170 && pthat <300 ) itest170++;
// 
//     if ( pthat >= 50 && pthat < 80 && itest50 > test50 ) continue; 
//     if ( pthat >= 80 && pthat <120 && itest80 >  test80 ) continue; 
//     if ( pthat >=120 && pthat <170 && itest120 >  test120 ) continue; 
//     if ( pthat >=170 && pthat <300 && itest170 >  test170 ) continue;
//$$
    
//*********************************
// Run certification:
//$$$$
//     if (pthat <= 0 && 
//      !( ( Run == 251244 && ( (LumiBlock >= 85 && LumiBlock <= 86) ||
//         		     (LumiBlock >= 88 && LumiBlock <= 93) ||
//         		     (LumiBlock >= 96 && LumiBlock <=121) ||
//         		     (LumiBlock >=123 && LumiBlock <=156) ||
//         		     (LumiBlock >=158 && LumiBlock <=428) ||
//         		     (LumiBlock >=430 && LumiBlock <=442) ) ) ||
//         ( Run == 251251 && ( (LumiBlock >=  1 && LumiBlock <= 31) ||
//         		     (LumiBlock >= 33 && LumiBlock <= 97) ||
//         		     (LumiBlock >= 99 && LumiBlock <=167) ) ) ||
//         ( Run == 251252 && ( (LumiBlock >=  1 && LumiBlock <=283) ||
//         		     (LumiBlock >=285 && LumiBlock <=505) ||
//         		     (LumiBlock >=507 && LumiBlock <=554) ) ) ||
//         ( Run == 251561 && ( (LumiBlock >=  1 && LumiBlock <= 94) ) ) ||
//         ( Run == 251562 && ( (LumiBlock >=  1 && LumiBlock <=439) ||
//         		     (LumiBlock >=443 && LumiBlock <=691) ) ) ||
//         ( Run == 251643 && ( (LumiBlock >=  1 && LumiBlock <=216) ||
//         		     (LumiBlock >=222 && LumiBlock <=606) ) ) ||
//         ( Run == 251721 && ( (LumiBlock >= 21 && LumiBlock <= 36) ) ) ||
//         ( Run == 251883 && ( (LumiBlock >= 56 && LumiBlock <= 56) ||
//         		     (LumiBlock >= 58 && LumiBlock <= 60) ||
//         		     (LumiBlock >= 62 && LumiBlock <=144) ||
//         		     (LumiBlock >=156 && LumiBlock <=437) ) )
//     ) ) continue;
//$$$$
    
    if(jentry%10000 == 0) std::cout << "the processed event number is " << jentry << " over " << nentries 
//$$
//     << "   count 50-80-120-170 " 
//     << itest50 << " " << itest80 << " " << itest120 << " " << itest170
//$$
    << std::endl;
   

//*********************************
// Leading Jet
//$$
//     int njet80 = 0;
//     for (int ijet = 0; ijet < nJet; ijet++) {
//       float residual = 1.;
//       float ptjet = Jet_pt[ijet] * residual;
//       float etajet = fabs(Jet_eta[ijet]);
//     if ( !(etajet < 2.4) ) continue;
//       if ( ptjet >  80. ) njet80++;
//     }
//     if ( njet80  < 1 ) continue;
//$$


//*********************************
// Loop on jets
   
    for (int ijet=0; ijet<nJet ; ijet++){
    
//$$
      if (Jet_pt[ijet]>30 && fabs(Jet_eta[ijet]) < 2.4 ){
//$$
      
        double jetP=0.0;
        vectProba.clear();
	
	for( int itrack=Jet_nFirstTrack[ijet]; itrack< Jet_nLastTrack[ijet]; itrack++){
	
	  double proba=-1000;
	   	  
	  bool passSelection=passTrackSel(Track_nHitAll[itrack],
	   Track_nHitPixel[itrack],Track_IP2D[itrack],
	   Track_pt[itrack],Track_length[itrack],Track_chi2[itrack],
//$$	   Track_zIP[itrack],Track_dist[itrack],Track_eta[itrack],
	   Track_dz[itrack],Track_dist[itrack],Track_eta[itrack],
	   Jet_eta[ijet],Track_phi[itrack],Jet_phi[ijet],Jet_pt[ijet]);   
	  
	  if(passSelection && Track_IPsig[itrack]<0) {          
	    
	    int catnum =-1;
	    
	    catnum = IsInCategory(Track_eta[itrack],  Track_nHitAll[itrack] , Track_nHitPixel[itrack], Track_p[itrack] 
	    , Track_chi2[itrack], catDefList);

	     if(catnum >=0) {
	      proba=calculTrackProba( Track_IPsig[itrack], catDefList[catnum]); 
	    
	      vectProba.push_back(proba);
	      vectTrkProba[catnum]->Fill(proba);
	     }else{
//$$	      
// 	      cout << "no category found for this track " <<  endl;
// 	      cout << "Track_p[itrack]              " << Track_p[itrack]    << endl;
// 	      cout << "Track_eta[itrack]            " << Track_eta[itrack]  << endl;
// 	      cout << "Track_chi2[itrack]           " << Track_chi2[itrack] << endl;
// 	      cout << "Track_nHitPixel[itrack]      " << Track_nHitPixel[itrack]<< endl;
// 	      cout << "Track_nHitAll[itrack]        " << Track_nHitAll[itrack]  << endl;
//$$	      
	    }
	  }	    
	}
	
        if (vectProba.size()!=0) jetP=jetProbability( vectProba);
	
	JetProba->Fill(jetP);
	
	if(Jet_flavour[ijet]==5)jetProbaDiscriB->Fill(Jet_Proba[ijet]);
	
	if(Jet_flavour[ijet]== 1 || Jet_flavour[ijet]== 2 || Jet_flavour[ijet]==3 || Jet_flavour[ijet]==21) jetProbaDiscriUDSG->Fill(Jet_Proba[ijet]);
      }
    }
  }  
  
  myfileout->cd();
  
  for(unsigned int i = 0; i<catDefList.size(); i++) vectTrkProba[i]->Write();
  
  trackPCat_all->Write();  
  jetProbaDiscriB ->Write(); 
  jetProbaDiscriUDSG->Write(); 
}


void JetProbaValidation::plotPerf()
{
}
