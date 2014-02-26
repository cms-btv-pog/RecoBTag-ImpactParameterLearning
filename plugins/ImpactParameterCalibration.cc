// -*- C++ -*-
//
// Package:    ImpactParameterCalibration
// Class:      ImpactParameterCalibration
// 
/**\class ImpactParameterCalibration ImpactParameterCalibration.cc RecoBTag/ImpactParameterCalibration/src/ImpactParameterCalibration.cc
   
Description: <one line class summary>

Implementation:
<Notes on implementation>
*/
//
// Original Author:  Jeremy Andrea/Andrea Rizzi
//         Created:  Mon Aug  6 16:10:38 CEST 2007
// $Id: ImpactParameterCalibration.cc,v 1.14 2010/02/11 00:13:30 wmtan Exp $
//
//
// system include files

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "Utilities/General/interface/FileInPath.h"
#include "RecoVertex/VertexPrimitives/interface/TransientVertex.h"

#include "FWCore/Framework/interface/IOVSyncValue.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CondCore/DBOutputService/interface/PoolDBOutputService.h"

#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "DataFormats/BTauReco/interface/TrackIPTagInfo.h"
#include "CondFormats/BTauObjects/interface/TrackProbabilityCalibration.h"
#include "CondFormats/BTauObjects/interface/CalibratedHistogram.h"
#include "RecoVertex/VertexPrimitives/interface/ConvertToFromReco.h"

#include "CondFormats/DataRecord/interface/BTagTrackProbability2DRcd.h"
#include "CondFormats/DataRecord/interface/BTagTrackProbability3DRcd.h"

#include "DataFormats/PatCandidates/interface/Jet.h"
#include "TClass.h"
#include "TrackingTools/IPTools/interface/IPTools.h"
#include "TBufferFile.h"

#include "TLorentzVector.h"
#include "TBufferXML.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <memory>
#include "TMath.h"
#include "TrackClassMatch.h"


using namespace reco;
using namespace std;
//
// class decleration
//

typedef std::vector<pat::Jet> PatJetCollection;

class ImpactParameterCalibration : public edm::EDAnalyzer {
   public:
      explicit ImpactParameterCalibration(const edm::ParameterSet&);
      ~ImpactParameterCalibration();
  
  
private:
  virtual void beginJob() ;
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;
  
  virtual void initFromFirstES(const edm::EventSetup&);
  virtual void processJets(const edm::Handle<PatJetCollection>& jetsColl, const  reco::Vertex  *pv, const edm::Event&, const edm::EventSetup&);
  edm::ParameterSet config;
  bool m_needInitFromES;
  TrackProbabilityCalibration * fromXml(edm::FileInPath xmlCalibration);
  
  static TrackProbabilityCategoryData createCategory(double  pmin,double  pmax,
						     double  etamin,  double  etamax,
						     int  nhitmin,  int  nhitmax,
						     int  npixelhitsmin, int  npixelhitsmax,
						     double cmin, double cmax, int withFirst)
  {
    TrackProbabilityCategoryData c;
    c.pMin=pmin;
    c.pMax=pmax;
    c.etaMin=etamin;
    c.etaMax=etamax;
    c.nHitsMin=nhitmin;
    c.nHitsMax=nhitmax;
    c.nPixelHitsMin=npixelhitsmin;
    c.nPixelHitsMax=npixelhitsmax;
    c.chiMin=cmin;
    c.chiMax=cmax;
    c.withFirstPixel=withFirst;
    return c;
  }
  
  TrackProbabilityCalibration * m_calibration[2];
  //edm::InputTag m_iptaginfo;
  edm::InputTag m_pv;
  edm::InputTag JetCollectionTag_;
  unsigned int minLoop, maxLoop;
  edm::ESHandle<TransientTrackBuilder> trackBuilder ;
};

ImpactParameterCalibration::ImpactParameterCalibration(const edm::ParameterSet& iConfig):config(iConfig)
{
  m_needInitFromES = false;
  //m_iptaginfo = iConfig.getParameter<edm::InputTag>("tagInfoSrc");
  m_pv = iConfig.getParameter<edm::InputTag>("primaryVertexColl");
  
  JetCollectionTag_ = iConfig.getParameter<edm::InputTag>("Jets");
  
  bool createOnlyOne = iConfig.getUntrackedParameter<bool>("createOnlyOneCalibration", true);
  minLoop=0;
  maxLoop=1;
  if (createOnlyOne == true){
    int whichCalib = iConfig.getUntrackedParameter<int>("dimension", 3);
    if (whichCalib==2){
      std::cout <<" Writing only 2D calibrations"<<std::endl;
      minLoop=1;
      maxLoop=1;
    }else if (whichCalib==3){
      std::cout <<" Writing only 3D calibrations"<<std::endl;
      minLoop=0;
      maxLoop=0;
    }else {
      std::cout <<" Dimension not found: "<<whichCalib<<"; it must be either 2 or 3"<<std::endl;
    }
  }
  
}


ImpactParameterCalibration::~ImpactParameterCalibration()
{
}


void ImpactParameterCalibration::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  
  if(m_needInitFromES) initFromFirstES(iSetup);
  using namespace edm;
  using namespace reco;
  
  
  iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder", trackBuilder);   
  
  
  
  edm::Handle <PatJetCollection> jetsColl;
  iEvent.getByLabel (JetCollectionTag_, jetsColl);
  
  Handle<reco::VertexCollection> primaryVertex;
  iEvent.getByLabel(m_pv,primaryVertex);
  
  const  reco::Vertex  *pv;
  bool pvFound = (primaryVertex->size() != 0);
  if ( pvFound ) {
    pv = &(*primaryVertex->begin());
  }
  else {
    reco::Vertex::Error e;
    e(0,0)=0.0015*0.0015;
    e(1,1)=0.0015*0.0015;
    e(2,2)=15.*15.;
    reco::Vertex::Point p(0,0,0);
    pv=  new reco::Vertex(p,e,1,1,1);
    //newvertex = true;
  }
  
  
  processJets(jetsColl, pv, iEvent, iSetup) ;
  
  
  
}





void ImpactParameterCalibration::initFromFirstES(const edm::EventSetup& iSetup)
{
  using namespace edm;
  CalibratedHistogram hist(config.getParameter<int>("nBins"),0,config.getParameter<double>("maxSignificance"));
  
  bool resetHistogram = config.getParameter<bool>("resetHistograms");
  ESHandle<TrackProbabilityCalibration> calib2DHandle;
  iSetup.get<BTagTrackProbability2DRcd>().get(calib2DHandle);
  ESHandle<TrackProbabilityCalibration> calib3DHandle;
  iSetup.get<BTagTrackProbability3DRcd>().get(calib3DHandle);
  const TrackProbabilityCalibration * ca[2];
  ca[0]  = calib3DHandle.product();
  ca[1]  = calib2DHandle.product();
  for(unsigned int i=minLoop;i <=maxLoop ;i++)
    for(unsigned int j=0;j<ca[i]->data.size() ; j++)
      {
	TrackProbabilityCalibration::Entry e;
	e.category=ca[i]->data[j].category;
	
	if(resetHistogram)
	  e.histogram=hist;
	else
	  e.histogram=ca[i]->data[j].histogram;
	
	m_calibration[i]->data.push_back(e);
      }
  
  
}


// ------------ method called once each job just before starting event loop  ------------
void 
ImpactParameterCalibration::beginJob()
{
  using namespace edm;
  m_calibration[0] =   new TrackProbabilityCalibration();
  m_calibration[1] =   new TrackProbabilityCalibration();
  
  CalibratedHistogram hist(config.getParameter<int>("nBins"),0,config.getParameter<double>("maxSignificance"));
  
  
  std::string categories = config.getParameter<std::string>("inputCategories");
  
  if(categories == "HardCoded")
    {
      vector<TrackProbabilityCategoryData> v;
      //TrackProbabilityCategoryData {pMin, pMax, etaMin, etaMax,
      //nHitsMin, nHitsMax, nPixelHitsMin, nPixelHitsMax, chiMin,chiMax, withFirstPixel;
      //trackQuality;
      v.push_back(createCategory(0, 5000, 0  , 2.5, 8 , 50, 1, 1, 0  , 5  , 0));
      v.push_back(createCategory(0, 5000, 0  , 2.5, 8 , 50, 2, 8, 2.5, 5  , 0));
      v.push_back(createCategory(0, 8   , 0  , 0.8, 8 , 50, 3, 8, 0  , 2.5, 0));
      v.push_back(createCategory(0, 8   , 0.8, 1.6, 8 , 50, 3, 8, 0  , 2.5, 0));
      v.push_back(createCategory(0, 8   , 1.6, 2.5, 8 , 50, 3, 8, 0  , 2.5, 0));
      v.push_back(createCategory(0, 8   , 0  , 2.5, 8 , 50, 2, 2, 0  , 2.5, 0));
      v.push_back(createCategory(8, 5000, 0  , 0.8, 8 , 50, 3, 8, 0  , 2.5, 0));
      v.push_back(createCategory(8, 5000, 0.8, 1.6, 8 , 50, 3, 8, 0  , 2.5, 0));
      v.push_back(createCategory(8, 5000, 1.6, 2.5, 8 , 50, 3, 8, 0  , 2.5, 0));
      v.push_back(createCategory(8, 5000, 0  , 2.5, 8 , 50, 2 ,2, 0  , 2.5, 0));
      for(unsigned int i=minLoop;i <=maxLoop ;i++)
	for(unsigned int j=0;j<v.size() ; j++)
	  {
	    TrackProbabilityCalibration::Entry e;
	    e.category=v[j];
	    e.histogram=hist;
	    m_calibration[i]->data.push_back(e);
	  }
      
    }
  if(categories == "RootXML")
    {
      bool resetHistogram = config.getParameter<bool>("resetHistograms");
      const TrackProbabilityCalibration * ca[2];
      ca[0]  = fromXml(config.getParameter<edm::FileInPath>("calibFile3d"));
      ca[1]  = fromXml(config.getParameter<edm::FileInPath>("calibFile2d"));
      
      for(unsigned int i=minLoop;i <=maxLoop ;i++)
	for(unsigned int j=0;j<ca[i]->data.size() ; j++)
	  {
	    TrackProbabilityCalibration::Entry e;
	    e.category=ca[i]->data[j].category;
	    
	    if(resetHistogram)
	      e.histogram=hist;
	    else
	      e.histogram=ca[i]->data[j].histogram;
	    
	    m_calibration[i]->data.push_back(e);
	  }
      
      delete ca[0];
      delete ca[1];
      
    }
  if(categories == "EventSetup")
    {
      m_needInitFromES=true;
    }

  
  
}

TrackProbabilityCalibration * ImpactParameterCalibration::fromXml(edm::FileInPath xmlCalibration)   
{
  std::ifstream xmlFile(xmlCalibration.fullPath().c_str());
  if (!xmlFile.good())
    throw cms::Exception("BTauFakeMVAJetTagConditions")
      << "File \"" << xmlCalibration.fullPath()
      << "\" could not be opened for reading."
      << std::endl;
  std::ostringstream ss;
  ss << xmlFile.rdbuf();
  xmlFile.close();
  TClass *classType = 0;
  void *ptr = TBufferXML(TBuffer::kRead).ConvertFromXMLAny(
							   ss.str().c_str(), &classType, kTRUE, kFALSE);
  if (!ptr)
    throw cms::Exception("ImpactParameterCalibration")
      << "Unknown error parsing XML serialization"
      << std::endl;
  
  if (std::strcmp(classType->GetName(),
		  "TrackProbabilityCalibration")) {
    classType->Destructor(ptr);
    throw cms::Exception("ImpactParameterCalibration")
      << "Serialized object has wrong C++ type."
      << std::endl;
  }
  
  return static_cast<TrackProbabilityCalibration*>(ptr);
}




// ------------ method called once each job just after ending the event loop  ------------
void 
ImpactParameterCalibration::endJob() {
  
  if(config.getParameter<bool>("writeToDB"))
    {
      edm::Service<cond::service::PoolDBOutputService> mydbservice;
      if( !mydbservice.isAvailable() ) return;
      if(minLoop == 0 )  mydbservice->createNewIOV<TrackProbabilityCalibration>(m_calibration[0], mydbservice->beginOfTime(), mydbservice->endOfTime(),"BTagTrackProbability3DRcd");
      if(maxLoop == 1)   mydbservice->createNewIOV<TrackProbabilityCalibration>(m_calibration[1],  mydbservice->beginOfTime(), mydbservice->endOfTime(),"BTagTrackProbability2DRcd");
    } 
  
  
  if(config.getParameter<bool>("writeToRootXML"))
    {
      if(maxLoop == 1 ){
	std::ofstream of2("2d.xml");
	TBufferXML b2(TBuffer::kWrite);
	of2 << b2.ConvertToXML(const_cast<void*>(static_cast<const void*>(m_calibration[1])),
			       TClass::GetClass("TrackProbabilityCalibration"),
			       kTRUE, kFALSE);
	of2.close();
      }
      if(minLoop == 0 ){
	std::ofstream of3("3d.xml");
	TBufferXML b3(TBuffer::kWrite);
	of3 << b3.ConvertToXML(const_cast<void*>(static_cast<const void*>(m_calibration[0])),
			       TClass::GetClass("TrackProbabilityCalibration"),
			       kTRUE, kFALSE);
	of3.close();
      }
    }
  
  
  if(config.getParameter<bool>("writeToBinary"))
    {
      if(maxLoop == 1 ){
	
	std::ofstream ofile("2d.dat");
	TBufferFile buffer(TBuffer::kWrite);
	buffer.StreamObject(const_cast<void*>(static_cast<const void*>(m_calibration[1])),
			    TClass::GetClass("TrackProbabilityCalibration"));
	Int_t size = buffer.Length();
	ofile.write(buffer.Buffer(),size);
	ofile.close();
      }
      if(minLoop == 0 ){
	std::ofstream ofile3("3d.dat");
	TBufferFile buffer3(TBuffer::kWrite);
	buffer3.StreamObject(const_cast<void*>(static_cast<const void*>(m_calibration[0])),
			     TClass::GetClass("TrackProbabilityCalibration"));
	Int_t size3 = buffer3.Length();
	ofile3.write(buffer3.Buffer(),size3);
	ofile3.close();
      }
    }
  
  
  
  
  
}
void ImpactParameterCalibration::processJets(const edm::Handle<PatJetCollection>& jetsColl, const  reco::Vertex  *pv, const edm::Event& iEvent, const edm::EventSetup& iSetup)
{  
  
  vector<TrackProbabilityCalibration::Entry>::iterator found;
  vector<TrackProbabilityCalibration::Entry>::iterator it_begin;
  vector<TrackProbabilityCalibration::Entry>::iterator it_end; 
  
  for ( PatJetCollection::const_iterator pjet = jetsColl->begin(); pjet != jetsColl->end(); ++pjet ) {
    
    double ptjet  = pjet->pt()  ;
    double jeteta = pjet->eta() ;
    double jetphi = pjet->phi() ;
    
    
    if (ptjet>10 && fabs(jeteta)<2.5) {
    
      float etajet = TMath::Abs( pjet->eta() );
      float phijet = pjet->phi();
      
      if (phijet < 0.) phijet += 2*TMath::Pi();
      
      const edm::RefVector<reco::TrackCollection>  &selected_tracks( pjet->tagInfoTrackIP("impactParameter")->selectedTracks() );
      
      int ntagtracks = pjet->tagInfoTrackIP("impactParameter")->probabilities(0).size();
      
      unsigned int trackSize = selected_tracks.size();
      
      for(unsigned int i=minLoop; i <= maxLoop;i++){ 
      
        int k=0;
	it_begin=m_calibration[i]->data.begin();
	it_end=m_calibration[i]->data.end();
	
	for (unsigned int itt=0; itt < trackSize; ++itt)
	  {
	    reco::Track  ptrack;
	    ptrack = *selected_tracks[itt];
	    
	    TransientTrack transientTrack = trackBuilder->build(ptrack);
	    GlobalVector direction(pjet->px(), pjet->py(), pjet->pz());
	    
	    //--------------------------------
	    Double_t decayLength=-1;
	    TrajectoryStateOnSurface closest = IPTools::closestApproachToJet(transientTrack.impactPointState(), *pv, direction, transientTrack.field());
	    if (closest.isValid())
	      decayLength =  (closest.globalPosition()-   RecoVertex::convertPos(pv->position())).mag();
	    else
	      decayLength = -1;
	    
	    Double_t distJetAxis =  IPTools::jetTrackDistance(transientTrack, direction, *pv).second.value();
	    
	    double Track_dist     = distJetAxis;
	    double Track_length   = decayLength;
	    
	    double Track_dxy      = ptrack.dxy(pv->position());
	    double Track_dz       = ptrack.dz(pv->position());
	    double Track_zIP      = ptrack.dz()-(*pv).z();
	    
	    float deta = ptrack.eta() - jeteta;
	    float dphi = ptrack.phi() - jetphi;
	    
	    if ( dphi > TMath::Pi() ) dphi = 2.*TMath::Pi() - dphi;
	    float deltaR = TMath::Sqrt(deta*deta + dphi*dphi);
	    
	    bool pass_cut_trk =false;
	    
	    TLorentzVector track4P, jet4P;
	    track4P.SetPtEtaPhiM(ptrack.pt(), ptrack.eta(), ptrack.phi(), 0. );
	    jet4P.SetPtEtaPhiM( pjet->pt(), pjet->eta(), pjet->phi(), 0 );
	    
	    if( jet4P.DeltaR(track4P) < 0.3 &&std::fabs(distJetAxis) < 0.07 && decayLength < 5.0 )pass_cut_trk=true;
	    
	    TrackClassMatch::Input input(ptrack,*pjet,*pv);
	    
	    if (pass_cut_trk ){
	      double Track_IP2D	 = pjet->tagInfoTrackIP("impactParameter")->impactParameterData()[k].ip2d.value();
	      double Track_IP2Dsig   = pjet->tagInfoTrackIP("impactParameter")->impactParameterData()[k].ip2d.significance();
	      double Track_IP	 = pjet->tagInfoTrackIP("impactParameter")->impactParameterData()[k].ip3d.value();
	      double Track_IPsig	 = pjet->tagInfoTrackIP("impactParameter")->impactParameterData()[k].ip3d.significance();
	      double Track_IP2Derr   = pjet->tagInfoTrackIP("impactParameter")->impactParameterData()[k].ip2d.error();
	      double Track_IPerr	 = pjet->tagInfoTrackIP("impactParameter")->impactParameterData()[k].ip3d.error();
	      double Track_Proba	 = pjet->tagInfoTrackIP("impactParameter")->probabilities(0)[k];	
	      double Track_p	 = ptrack.p();
	      double Track_pt	 = ptrack.pt();
	      double Track_eta	 = ptrack.eta();
	      double Track_phi	 = ptrack.phi();
	      double Track_chi2	 = ptrack.normalizedChi2();
	      int    Track_nHitAl   = ptrack.numberOfValidHits();
	      int    Track_nHitPixel= ptrack.hitPattern().numberOfValidPixelHits();
	      int    Track_nHitStrip= ptrack.hitPattern().numberOfValidStripHits();
	      int    Track_nHitTIB  = ptrack.hitPattern().numberOfValidStripTIBHits();
	      int    Track_nHitTID  = ptrack.hitPattern().numberOfValidStripTIDHits();
	      int    Track_nHitTOB  = ptrack.hitPattern().numberOfValidStripTOBHits();
	      int    Track_nHitTEC  = ptrack.hitPattern().numberOfValidStripTECHits();
	      int    Track_nHitPXB  = ptrack.hitPattern().numberOfValidPixelBarrelHits();
	      int    Track_nHitPXF  = ptrack.hitPattern().numberOfValidPixelEndcapHits(); 	    
	
	      
//  	    cout << "-------------------------- " << endl;
// // 	    
//  	    cout << "jet pt " <<   pjet->pt()<<endl;
//  	    cout << "jet eta " <<   pjet->eta()<<endl;
// 
// 	  
// // 	      
//  	    cout << "track pt " << ptrack.pt() <<endl;
//  	    cout << "track eta " << ptrack.eta() <<endl;	
//  	    cout << "track zIP " << ptrack.dz()-(*pv).z()  <<endl; 
//  	    cout << "track IP2D " << pjet->tagInfoTrackIP("impactParameter")->impactParameterData()[k].ip2d.value()  <<endl;	    
//  	    cout << "track dist " <<  distJetAxis <<endl;	    
//  	    cout << "track lenght " << decayLength  <<endl;
// // 	    
//  	    cout << "pix hits " << ptrack.hitPattern().numberOfValidPixelHits() <<endl;
//  	    cout << "pix+sistrip hits " << ptrack.hitPattern().numberOfValidStripHits()+ptrack.hitPattern().numberOfValidPixelHits() <<endl;
//  	    cout << "norm. chi2 " << ptrack.normalizedChi2()<<endl;
	      
	      found = std::find_if(it_begin,it_end,bind1st(TrackClassMatch(),input));
	      
	      if(found!=it_end && Track_IPsig<0) found->histogram.fill(-Track_IPsig);
	      
 	      //else                               std::cout << "No category for this track!!" << std::endl;
	      
	    }
	    
	    ++k;	
	  }//end tracks
	
      }//end loop categories  
      
    }//end jet conditions
    
    
  }//end loop jets
  
  
}




//define this as a plug-in
DEFINE_FWK_MODULE(ImpactParameterCalibration);
