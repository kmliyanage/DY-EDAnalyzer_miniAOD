// -*- C++ -*-
//
// Package:    Demo/DemoAnalyzer
// Class:      DemoAnalyzer
//
/**\class DemoAnalyzer DemoAnalyzer.cc Demo/DemoAnalyzer/plugins/DemoAnalyzer.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Kalpanie Madara Liyanage
//         Created:  Sat, 30 May 2020 06:29:45 GMT
//
//


// system include files
#include <memory>
#include <iomanip>

#include <TLorentzVector.h>
#include <TVector3.h>
#include <algorithm>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
 #include "FWCore/Utilities/interface/InputTag.h"
 #include "DataFormats/TrackReco/interface/Track.h"
 #include "DataFormats/TrackReco/interface/TrackFwd.h"



//////////////////////////////////////

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonSelectors.h"
#include "DataFormats/MuonReco/interface/MuonPFIsolation.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/PatCandidates/interface/PackedGenParticle.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/HLTReco/interface/TriggerObject.h"
#include "DataFormats/PatCandidates/interface/TriggerObjectStandAlone.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "Math/GenVector/VectorUtil.h"
#include "TTree.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/Math/interface/deltaPhi.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"


#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"

#include "TrackingTools/Records/interface/TrackingComponentsRecord.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "RecoVertex/VertexPrimitives/interface/TransientVertex.h"
#include "DataFormats/MuonReco/interface/MuonCocktails.h"
#include "MagneticField/Engine/interface/MagneticField.h"
#include "CommonTools/Statistics/interface/ChiSquaredProbability.h"


#include "CommonTools/Utils/interface/StringCutObjectSelector.h"
#include "DataFormats/PatCandidates/interface/CompositeCandidate.h"

#include "RecoVertex/KalmanVertexFit/interface/KalmanVertexFitter.h"
#include "RecoVertex/VertexTools/interface/InvariantMassFromVertex.h"


#include "RecoVertex/VertexTools/interface/SequentialVertexFitter.h"
#include "RecoVertex/KinematicFit/interface/KinematicParticleVertexFitter.h"
#include "RecoVertex/KinematicFitPrimitives/interface/KinematicParticleFactoryFromTransientTrack.h"
#include "RecoVertex/KinematicFit/interface/MassKinematicConstraint.h"
#include "RecoVertex/KinematicFit/interface/KinematicParticleFitter.h"
#include "RecoVertex/KinematicFitPrimitives/interface/MultiTrackKinematicConstraint.h"
#include "RecoVertex/KinematicFit/interface/KinematicConstrainedVertexFitter.h"
#include "RecoVertex/KinematicFit/interface/TwoTrackMassKinematicConstraint.h"
#include "RecoVertex/KinematicFitPrimitives/interface/KinematicParticle.h"
#include "RecoVertex/KinematicFitPrimitives/interface/RefCountedKinematicParticle.h"
#include "RecoVertex/KinematicFitPrimitives/interface/TransientTrackKinematicParticle.h"
#include "DataFormats/Candidate/interface/Particle.h"

//
// class declaration
//

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<>
// This will improve performance in multithreaded jobs.


//using reco::TrackCollection;

class DemoAnalyzer : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
   public:
      explicit DemoAnalyzer(const edm::ParameterSet&);
      ~DemoAnalyzer();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
      virtual void beginJob() override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;

      // ----------member data ---------------------------
    //edm::EDGetTokenT<TrackCollection> tracksToken_;  //used to select what tracks to read from configuration file
    edm::EDGetTokenT<pat::MuonCollection> muonCollToken;
    edm::EDGetTokenT<pat::PackedGenParticleCollection> genCollToken;
    edm::EDGetTokenT<reco::VertexCollection> vertexCollToken;
    //edm::ESHandle<TransientTrackBuilder> ttkb;
    
    //variables
    
    float HLT_pt,HLT_eta,HLT_phi;
    float ptEffCut;
    float PtDYTRecMu1,PtDYTRecMu2,PtRecTunePMu1,PtRecTunePMu2,PtRecTunePMu3,
    PtRecMuBestTrack1,PtRecMuBestTrack2,PtRecMuBestTrack3;
    float RecoHLTMatchingDeltaRcut,deltaRcut,minMassCut,maxMassCut;
    float vtxChi2Mu,vtxMassMu;
    float mPtGen1,mPhiGen1,mEtaGen1,mEnGen1;
    unsigned mGenFlag1;
    float mPtGen2,mPhiGen2,mEtaGen2,mEnGen2;
    int ChargeRecMu1,ChargeRecMu2,ChargeRecMu3;
    unsigned flagmu1;
    unsigned flag1;
    float PtRecTunePMuBestTrack1,EnRecMu1,EtaRecMu1,PhiRecMu1;
    float PtRecTunePMuBestTrack2,EnRecMu2,EtaRecMu2,PhiRecMu2;
    float PtRecTunePMuBestTrack3,EnRecMu3,EtaRecMu3,PhiRecMu3;
    float pxRecMu1,pyRecMu1,pzRecMu1,pRecMu1,dxyRecMu1;
    float pxRecMu2,pyRecMu2,pzRecMu2,pRecMu2,dxyRecMu2;
    float pxRecMu3,pyRecMu3,pzRecMu3,pRecMu3,dxyRecMu3;
    float genET1,genPhi1,genEta1,genEn1;
    int genID1,genStat1;
    float genET2,genPhi2,genEta2,genEn2;
    int genID2,genStat2;
    float MassGen,RecoMass;
    int NbGen,NbReco;
    int nbTP,nbTT,nbTF;
    float TagProbeEtaCut;
    float Eff;
    float MassCutMin,MassCutMax;
    float MassResolution;
    float EtaCut;
    
    
    
    ///make high_pt muon pairs.....////////////////////////
    
    int index1 = -1;
    int index2 = -1;
    int temp = -1;
    int dil_charge = 0;
    double cos_angle = -1.0;
    double vertex_m = -1.0;
    double dil_m = -1.0;
    double vertex_chi2 = -1.0;
    bool GoodDataRan = false;
    bool GoodVtx = false;
    double sum_pt = -1.0;
    double max_sum_pt = -1.0;
    const double Z_mass = 91.188;
    float close_md = -1.0;
    float md = -1.0;
    
    
 
    
    bool choose = false;
    
    
    
    UInt_t  run;
    UInt_t  lumi;
    UInt_t event;
    
    
    //histos
    TH1F* h_genpt;
    TH1F* h_pt;
    TH1F* high_pt_muon_pt;

};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
DemoAnalyzer::DemoAnalyzer(const edm::ParameterSet& iConfig)

{
   //now do what ever initialization is needed
    
    
    edm::InputTag theMuonLabel("slimmedMuons");
    muonCollToken = consumes<pat::MuonCollection>(theMuonLabel);
    
    
    edm::InputTag theGenMuonLabel("packedGenParticles");
    genCollToken = consumes<pat::PackedGenParticleCollection>(theGenMuonLabel);
    
    edm::InputTag theVertexLabel("offlineSlimmedPrimaryVertices");
    vertexCollToken = consumes<reco::VertexCollection>(theVertexLabel);
    
    edm::InputTag theTrackLabel("offlineSlimmedPrimaryVertices");
    vertexCollToken = consumes<reco::VertexCollection>(theVertexLabel);
    
    usesResource("TFileService");
    edm::Service<TFileService> fs;
    
    h_pt = fs->make<TH1F>("pt", "RECO pt", 36, 50.0, 1850.0);
    h_genpt = fs->make<TH1F>("genpt", "GEN pt", 36, 50.0, 1850.0);
    high_pt_muon_pt = fs->make<TH1F>("high_pt_muon_pt", "high_pt_muon_pt", 36, 50.0, 1850.0);
}


DemoAnalyzer::~DemoAnalyzer()
{

   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void DemoAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
    using namespace edm;
    using namespace std;
    using namespace reco;
    using namespace pat;
    using namespace l1extra;
    
    
    /////////////////////////////////////////////////////////////////////////////
    // Event Information
    //
    run = iEvent.id().run();
    lumi = iEvent.luminosityBlock();
    event = iEvent.id().event();
    
    cout << " Run : " << run << " Lumi : " << lumi << " Event : " << event <<endl;
    
    /////////////////////////////////////////////////////////////////////////////
    
    
    
   /*
    edm::Handle<pat::PackedGenParticleCollection> genColl;
    iEvent.getByToken(genCollToken, genColl);
    int n = 0;
    for (auto it = genColl->cbegin(); it != genColl->cend(); ++it){
        if ( abs((*it).pdgId()) == 13 && fabs((*it).eta()) < 2.4 && (*it).pt() > 1.5 )
            n++;
    }
    cout << "Number of GEN muons: " << n << endl;

    for (auto it = genColl->cbegin(); it != genColl->cend(); ++it) {
        
        const pat::PackedGenParticle& mcParticle = (*it);
        if ( abs(mcParticle.pdgId()) != 13 ) continue; // skip this particle if it isn't a muon
        h_genpt->SetFillColor(38);
        h_genpt->Fill(mcParticle.pt());
        
    }
    
  */
    
    
    //////////////////////////////////////////////////////////////////////////////////
    // Primary Vertex(PV) selection
    //
    edm::Handle<reco::VertexCollection>vertexColl;
    iEvent.getByToken(vertexCollToken, vertexColl);
    
    ///To find a good pv
    bool goodVertex = false;
    //const reco::Vertex* bestVtx;
    auto bestVtx = vertexColl->cbegin();
    std::cout<<"Number of PVs (pile-up) : "<<vertexColl->size()<<std::endl;
    
    for(auto pv = vertexColl->cbegin(); pv != vertexColl->cend(); ++pv){
        if (!pv->isValid()) continue;
        if (fabs(pv->position().z()) > 24.0) continue;
        if (sqrt((pv->position().x()*pv->position().x()) + (pv->position().y()*pv->position().y())) > 2.0) continue;
        if (pv->nTracks() < 4) continue;
        
        goodVertex = true;
        bestVtx = pv;
        break;
    }
    
    if(goodVertex)
        std::cout<< "Vertex info : " << "Z : " << bestVtx->position().z() << "    " << "rho : " << sqrt((bestVtx->position().x()*bestVtx->position().x()) + (bestVtx->position().y()*bestVtx->position().y())) << "    " << "No.tracks : " << bestVtx->nTracks() << std::endl;
    
   
    
    
        //////////////////////////////////////////////////////////////////////////////////
        // RECO Muons
        //
        edm::Handle<vector<pat::Muon>> muonColl;
        iEvent.getByToken(muonCollToken, muonColl);
        // cout << "Number of RECO muons: " << muonColl->size() << endl;
        
        
        //get the builder:
        //const edm::InputTag src;
        edm::ESHandle<TransientTrackBuilder> ttkb;
        //iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder" ,ttkb);
    
    
        for(auto it1 = muonColl->cbegin(); it1 != muonColl->cend(); ++it1){
        
            h_pt->SetFillColor(44);
            h_pt->Fill(it1->pt());
        
            if( fabs(it1->eta()) < 2.4 &&
               //it1->isGoodMuon() &&
               it1->isGlobalMuon() &&
               it1->isTrackerMuon() &&
               it1->tunePMuonBestTrack()->pt() > 53. &&
               it1->tunePMuonBestTrack()->ptError() / it1->tunePMuonBestTrack()->pt() < 0.3 &&
               fabs(it1->dB()) < 0.2 &&
               it1->isolationR03().sumPt / it1->pt() &&
               it1->innerTrack()->hitPattern().trackerLayersWithMeasurement() > 5 &&
               it1->innerTrack()->hitPattern().numberOfValidPixelHits() > 0 &&
               (it1->globalTrack()->hitPattern().numberOfValidMuonHits() > 0 || it1->tunePMuonBestTrack()->hitPattern().numberOfValidMuonHits() > 0 ) &&
               (it1->numberOfMatchedStations() > 1 || (it1->numberOfMatchedStations() == 1 && (it1->expectedNnumberOfMatchedStations()< 2 || !(it1->stationMask() == 1 || it1->stationMask() == 16) || it1->numberOfMatchedRPCLayers() > 2)))
           
               ){
            
                    index1 = std::distance(muonColl->cbegin(), it1);
                    index2 = index1+1;
                    //int ddd = 0;
            
                    temp = index1;
            
                    high_pt_muon_pt->SetFillColor(44);
                    high_pt_muon_pt->Fill(it1->pt());
            
            
                    for(auto it2 = it1+1;  it2 != muonColl->cend(); ++it2){
                    
                        if(it1==it2) continue;
                    
                        //opposite charge
                        if( (it1->charge())*(it2->charge()) == 1) continue;

                    
                        if( fabs(it2->eta()) < 2.4 &&
                           //it2->isGoodMuon() &&
                           it2->isGlobalMuon() &&
                           it2->isTrackerMuon() &&
                           it2->tunePMuonBestTrack()->pt() > 53. &&
                           it2->tunePMuonBestTrack()->ptError() / it2->tunePMuonBestTrack()->pt() < 0.3 &&
                           fabs(it1->dB()) < 0.2 && //for PAT Muons
                           //fabs(it1->muonBestTrack()->dxy(firstGoodVertex->position())) < 0.2 && //for Reco Muons
                           it2->isolationR03().sumPt / it2->pt() &&
                           it2->innerTrack()->hitPattern().trackerLayersWithMeasurement() > 5 &&
                           it2->innerTrack()->hitPattern().numberOfValidPixelHits() > 0 &&
                           (it2->globalTrack()->hitPattern().numberOfValidMuonHits() > 0 || it2->tunePMuonBestTrack()->hitPattern().numberOfValidMuonHits() > 0 ) &&
                           (it2->numberOfMatchedStations() > 1 || (it2->numberOfMatchedStations() == 1 && (it2->expectedNnumberOfMatchedStations()< 2 || !(it2->stationMask() == 1 || it2->stationMask() == 16) || it2->numberOfMatchedRPCLayers() > 2)))
                                                              
                           ){
                    
                        
                        
                                ///////////////////////coditions for dimuon////////////////////////////////
                        
                                cos_angle = it1->momentum().Dot(it2->momentum()) / it1->p() / it2->p();
                                dil_charge = it1->charge()*it2->charge();
                                max_sum_pt = it1->pt() + it2->pt();
                                dil_m = (it1->p4() + it2->p4()).M();
                                //close_md = fabs(dil_m-Z_mass);
                        
                                //For the vertex fits
                            
                                //Get TransientTrack for each selected muon
                                std::vector<reco::TransientTrack> ttv; //vector to store TransientTracks
                            
                            
                              
                                const reco::TrackRef& tk1 = it1->tunePMuonBestTrack();
                                ttv.push_back(ttkb->build(tk1));
                                const reco::TrackRef tk2 = it2->tunePMuonBestTrack();
                                //ttv.push_back(tk2);
                                ttv.push_back(ttkb->build(tk2));
                            
                                //KalmanVertexFitter kvf(true);
                                //CachingVertex<5> v = kvf.vertex(ttv);
                            
                                //vertex_chi2 = v.totalChiSquared()/v.degreesOfFreedom();
                                //static const double muon_mass = 0.1056583;
                                //InvariantMassFromVertex imfv;
                                //InvariantMassFromVertex::LorentzVector p4 = imfv.p4(v, muon_mass);
                                //Measurement1D mass = imfv.invariantMass(v, muon_mass);
                                //vertex_m = mass.value();
                        
                            
                                //std::cout<< "Found: " << ttv.size() << std::endl;
                            
                               // std::cout<<index1<<"    "<<index2<<"    "<<cos_angle<<" "<<dil_charge<<"    "<<max_sum_pt<<"    "<<dil_m<<" "<<vertex_m<<"   "<<vertex_chi2<<std::endl;
                        
                    
                            } //end conditions on it2
                    
                    } //end it2
    
            } //end conditions on it1
    
        } //end it1
        
    //} //end goodVertex

}
// ------------ method called once each job just before starting event loop  ------------
void
DemoAnalyzer::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void
DemoAnalyzer::endJob()
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
DemoAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);

  //Specify that only 'tracks' is allowed
  //To use, remove the default given above and uncomment below
  //ParameterSetDescription desc;
  //desc.addUntracked<edm::InputTag>("tracks","ctfWithMaterialTracks");
  //descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(DemoAnalyzer);

