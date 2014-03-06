#define JetProbaCalib_cxx
#include "JetProbaCalib.h"
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



std::ofstream out_txt("CalibrationFiles/Calib3D.xml", ios::out);

//-----------------------------------------------------Find the category of a given track----------------------------------------------------------------------------------------//
int JetProbaCalib::IsInCategory(float trkEta, float trkHTrk, float trkHPix, float trkp, float trkChi2, std::vector<CategoryDef >  d){
  
  int isIn=-1;
  
  
  for(unsigned int i=0; i<d.size(); i++){
   if (fabs(trkEta)>d[i].etaMin&&fabs(trkEta)<d[i].etaMax&&trkHTrk<=d[i].nHitsMax&&trkHTrk>=d[i].nHitsMin&&
     trkHPix<=d[i].nPixelHitsMax&&trkHPix>=d[i].nPixelHitsMin&&trkp<d[i].pMax&&trkp>d[i].pMin&&
     trkChi2>=d[i].chiMin&&trkChi2<=d[i].chiMax) isIn=i;
  }
  
  return isIn;
  
}


//-----------------------------------------------------Return true if track passes btag selection ---------------------------------------------------------------------------------//
bool JetProbaCalib::passTrackSel(int trk, int pix, float ip2d, float pt, float len, float chi2, float zip, float dist, float eta, float etaJet, float phi, float phiJet, float ptJet)
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


  
 
//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------//
//--------------------------------------------------FIRST LOOP TO FILL CATEGORIES------------------------------------------------------------------------------------------------//
//------------------------------------------------------------------------------------------------------------------------------------------------------------------------------//

void JetProbaCalib::Loop()
{
  double ptmax=10;
  
  TFile* myfile = new TFile("CalibrationFiles/calibeHistoWrite_std.root", "recreate");
  myfile->cd();
  
  if (fChain == 0) return;
  
  Long64_t nentries = fChain->GetEntriesFast();
  
  Long64_t nbytes = 0, nb = 0;
  
  //----------------------------------
  //Setup categories
  //----------------------------------
  CategoriesDefinition catDef;
  std::vector<CategoryDef >  vectCategories =  catDef.CategoryList();
  
  
  
  
  //-------------------------------------// EVENT LOOP //---------------------------------------------------------------------------------------------//
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0 ) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;
    // if (Cut(ientry) < 0) continue;
    
    if(jentry%10000 == 0) std::cout << "the processed event number is " << jentry << " over " << nentries << endl;
   
    //if(jentry>=1000000) break; 
    
    
    for (int ijet=0; ijet<nJet ; ijet++){
      
      if (fabs(Jet_eta[ijet]) > 2.5 || Jet_pt[ijet]<10 ) continue;
      for( int itrack=Jet_nFirstTrack[ijet]; itrack< Jet_nLastTrack[ijet]; itrack++){
        
        bool passSelection=passTrackSel(Track_nHitAll[itrack],Track_nHitPixel[itrack],Track_dxy[itrack],Track_pt[itrack],Track_length[itrack],Track_chi2[itrack],Track_dz[itrack],Track_dist[itrack],Track_eta[itrack],Jet_eta[ijet],Track_phi[itrack],Jet_phi[ijet],Jet_pt[ijet]);		
        
	if(Track_IP[itrack]>0 || !passSelection) continue;

        		
        int catnumb = -1;   
        catnumb = IsInCategory(Track_eta[itrack],Track_nHitAll[itrack] ,Track_nHitPixel[itrack] ,Track_p[itrack] ,
	Track_chi2[itrack], vectCategories);
        
        
        if(catnumb >= 0) vectCategories[catnumb].histo->Fill(-1*Track_IPsig[itrack]) ;
	else{
	  //cout << "no category found for this track " << endl;
	  //cout << "Track_p[itrack]          " << Track_p[itrack] << endl;
	  //cout << "Track_eta[itrack]        " << Track_eta[itrack] << endl;
	  //cout << "Track_nHitAll[itrack]    " << Track_nHitAll[itrack] << endl;
	  //cout << "Track_nHitPixel[itrack]  " << Track_nHitPixel[itrack] << endl;
	  //cout << "Track_chi2[itrack]       " << Track_chi2[itrack] << endl;
	}
	
        
      } //end loop on tracks
    } //end lopp on jets
  }//end lopp on events
  

	
	  
   //Store histograms into xml for the DB

   xmlFirst();
   
   for(unsigned int i=0; i< vectCategories.size(); i++){
    producexml( vectCategories[i]);
   }  
   
   xmlLast();


  
  //store histograms into root format
  CategoryDefCollection  catDefVector;
  for(int i=0; i< vectCategories.size(); i++){
    (catDefVector.ListOfCategory).push_back((vectCategories[i]));
  }   
  
  CategoryDefCollection * pcatDefVector = & catDefVector;
  //histoVector             ClassvectHisto;
  //ClassvectHisto.ListOfHisto = vectHisto;
  
  //ClassvectHisto.Write();
  pcatDefVector->Write();
  
  for(int i=0; i< vectCategories.size(); i++){
    ((*pcatDefVector).ListOfCategory)[i].histo->Write();
  }   




    
  
}




//-----------------------------------------------------------//
//-----------For output into XML format----------------------//
//-----------------------------------------------------------//



void  JetProbaCalib::xmlFirst()
{
  out_txt << "<Object class=\"TrackProbabilityCalibration\">" << endl;
  out_txt << "  <Class name=\"TrackProbabilityCalibration\" version=\"1\">" << endl;
  out_txt << "    <Member name=\"data\">" << endl;
  out_txt << "      <Item name=\"Version\" v=\"6\"/>" << endl;
  out_txt << "      <Item name=\"Int_t\" v=\"10\"/>" << endl;
 
}


void  JetProbaCalib::producexml( CategoryDef cat )
{ 

  // TH1D* hcat = (TH1D*)gROOT->FindObject((catname.c_str()));  
  
  
  //TH1F *hcat = &cat.histo;
  
  TH1D* hcat = cat.histo;  
  
  out_txt << "       <Object class=\"TrackProbabilityCalibration::Entry\">" << endl; 
  out_txt << "        <Class name=\"TrackProbabilityCalibration::Entry\" version=\"1\">" << endl;
  out_txt << "  	<Member name=\"category\">"  << endl;
  out_txt << "  	  <Class name=\"TrackProbabilityCategoryData\" version=\"1\"> " << endl;
  

  out_txt << "  	       <Member name=\"pMin\" v=\""<<	       cat.pMin  << "\"/>" << endl;
  out_txt << "  	       <Member name=\"pMax\" v=\""<<	       cat.pMax  << "\"/>" << endl;
  out_txt << "  	       <Member name=\"etaMin\" v=\""<<         cat.etaMin << "\"/>" << endl;
  out_txt << "  	       <Member name=\"etaMax\" v=\""<<         cat.etaMax << "\"/>" << endl;
  out_txt << "  	       <Member name=\"nHitsMin\" v=\""<<       cat.nHitsMin << "\"/>" << endl;
  out_txt << "  	       <Member name=\"nHitsMax\" v=\""<<       cat.nHitsMax << "\"/>" << endl;
  out_txt << "  	       <Member name=\"nPixelHitsMin\" v=\""<<  cat.nPixelHitsMin  << "\"/>" << endl;
  out_txt << "  	       <Member name=\"nPixelHitsMax\" v=\""<<  cat.nPixelHitsMax << "\"/>" << endl;
  out_txt << "  	       <Member name=\"chiMin\" v=\""<<         cat.chiMin << "\"/>" << endl;
  out_txt << "  	       <Member name=\"chiMax\" v=\""<<         cat.chiMax << "\"/>" << endl;
  out_txt << "  	       <Member name=\"withFirstPixel\" v=\""<< cat.withFirstPixelHit << "\"/>" << endl;
  out_txt << "  	       <Member name=\"trackQuality\" v=\""<< 1000 << "\"/>" << endl;
  


  out_txt << "  	       </Class>"									 << endl;
  out_txt << "  	      </Member>"									 << endl;
  out_txt << "  	      <Member name=\"histogram\">"							 << endl;
  out_txt << "  	       <Class name=\"PhysicsTools::Calibration::Histogram_float_float_\" version=\"1\">" << endl;
  out_txt << "  		 <Member name=\"binULimits\">"  						 << endl;
  out_txt << "  		  <Item name=\"Version\" v=\"6\"/>"						 << endl;
  out_txt << "  		  <Item name=\"Int_t\" v=\"0\"/>"						 << endl;
  out_txt << "  		 </Member>"									 << endl;
  out_txt << "  		 <Member name=\"binValues\">"							 << endl;
  out_txt << "  		  <Item name=\"Version\" v=\"6\"/>"						 << endl;
  out_txt << "  		  <Item name=\"Int_t\" v=\""<<hcat->GetNbinsX() <<"\"/>"		       << endl;
  out_txt << "  		  <Item name=\"Array\">  "							 << endl;
  
       
  nbins = hcat->GetNbinsX();

  for(int i=1;i <= nbins;i++)
    {   double entry = hcat->GetBinContent(i);
        out_txt << "                      <Item name=\"Float_t\" v=\""<<
	fixed<<setprecision(6)<<entry <<  "\"/>" << endl; 
      //cout << " get bin x " << hcat->GetBinContent(i) <<  endl; 
    }

 out_txt << "		    </Item> "<< endl;
 out_txt << "		   </Member> "<< endl;
 out_txt << "		   <Member name=\"limits\"> "<< endl;
 out_txt << "		     <Class name=\"PhysicsTools::Calibration::Range_float_\" version=\"1\"> "<< endl;
 out_txt << "		       <Member name=\"min\" v=\"0.000000\"/> "<< endl;
 out_txt << "		       <Member name=\"max\" v=\""<< hcat->GetXaxis()->GetXmax() <<"\"/> "<< endl;
 out_txt << "		     </Class> "<< endl;
 out_txt << "		   </Member> "<< endl;
 out_txt << "		 </Class> "<< endl;
 out_txt << "	       </Member> "<< endl;
 out_txt << "	     </Class> "<< endl;
 out_txt << "	   </Object> "<< endl;




}


void  JetProbaCalib::xmlLast(){

 out_txt << "	 </Member>"<< endl;
 out_txt << "  </Class>"<< endl;
 out_txt << "</Object>"<< endl;

}






