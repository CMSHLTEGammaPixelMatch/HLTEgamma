/** \class HLTElectronPixelMatchFilter
 *
 * $Id: HLTElectronPixelMatchFilter.cc,v 1.15 2012/03/06 10:13:59 sharper Exp $
 *
 *  \author Aidan Randle-Conde (ULB)
 *
 */

#include "HLTrigger/Egamma/interface/HLTElectronPixelMatchFilterAnalysis.h"

#include "DataFormats/Common/interface/Handle.h"

#include "DataFormats/HLTReco/interface/TriggerFilterObjectWithRefs.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "DataFormats/Common/interface/AssociationMap.h"

#include "DataFormats/CaloRecHit/interface/CaloCluster.h"
#include "DataFormats/CaloRecHit/interface/CaloClusterFwd.h"

#include "DataFormats/EgammaReco/interface/SuperCluster.h"
#include "DataFormats/EgammaReco/interface/SuperClusterFwd.h"
#include "DataFormats/EgammaReco/interface/ElectronSeed.h"
#include "DataFormats/EgammaReco/interface/ElectronSeedFwd.h"

#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackBase.h"

#include "DataFormats/RecoCandidate/interface/RecoEcalCandidate.h"
#include "DataFormats/RecoCandidate/interface/RecoEcalCandidateFwd.h"
#include "DataFormats/RecoCandidate/interface/RecoCandidate.h"

#include "DataFormats/TrajectorySeed/interface/TrajectorySeedCollection.h"

#include <string>
#include <math.h>
#include <iostream>

//
// constructors and destructor
//
HLTElectronPixelMatchFilterAnalysis::HLTElectronPixelMatchFilterAnalysis(const edm::ParameterSet& iConfig) : HLTFilter(iConfig)
,file_out_(0)
,tree_out_(0)
,tree_par_(0)
,filename_("test.root")
,nmatch_store_limit_(-1)
,el_charge_(0)
,el_pt_(0)
,el_eta_(0)
,el_phi_(0)
,el_E_(0)
,el_nPixelMatch_(0)
,el_matchSuccess_(0)
,el_n_(0)
,el_nPass_(0)
,el_dRz1_neg_(0)
,el_dRz2_neg_(0)
,el_dPhi1_neg_(0)
,el_dPhi2_neg_(0)
,el_dRz1_pos_(0)
,el_dRz2_pos_(0)
,el_dPhi1_pos_(0)
,el_dPhi2_pos_(0)
,el_s2_neg_(0)
,el_s2_pos_(0)
,el_subDet1_(0)
,el_subDet2_(0)
,el_matchCharge_(0)
,el_n_neg_(0)
,el_n_pos_(0)
,el_chainNumber_(0)
,ePhiMin1_(0)
,ePhiMax1_(0)
,pPhiMin1_(0)
,pPhiMax1_(0)
,PhiMin2_(0)
,PhiMax2_(0)
,r2MinF_(0)
,r2MaxF_(0)
,rMinI_(0)
,rMaxI_(0)
,z2MinB_(0)
,z2MaxB_(0)
,s_a_phi1B_(-999)
,s_a_phi1I_(-999)
,s_a_phi1F_(-999)
,s_a_phi2B_(-999)
,s_a_phi2I_(-999)
,s_a_phi2F_(-999)
,s_a_zB_(-999)
,s_a_rI_(-999)
,s_a_rF_(-999)
,s_b_phi1B_(0)
,s_b_phi1I_(0)
,s_b_phi1F_(0)
,s_b_phi2B_(0)
,s_b_phi2I_(0)
,s_b_phi2F_(0)
,s_b_zB_(0)
,s_b_rI_(0)
,s_b_rF_(0)
,chainNumber_(0)
,el_best_s2_(0)
,el_best_index_(0)
,el_best_charge_(0)
{
  candTag_               = iConfig.getParameter< edm::InputTag > ("candTag");
  L1IsoPixelSeedsTag_    = iConfig.getParameter< edm::InputTag > ("L1IsoPixelSeedsTag");
  L1NonIsoPixelSeedsTag_ = iConfig.getParameter< edm::InputTag > ("L1NonIsoPixelSeedsTag");
  npixelmatchcut_        = iConfig.getParameter< double >        ("npixelmatchcut");
  ncandcut_              = iConfig.getParameter< int >           ("ncandcut");
  doIsolated_            = iConfig.getParameter< bool >          ("doIsolated");
  L1IsoCollTag_          = iConfig.getParameter< edm::InputTag > ("L1IsoCand");
  L1NonIsoCollTag_       = iConfig.getParameter< edm::InputTag > ("L1NonIsoCand");
  filename_              = iConfig.getParameter< std::string   > ("Ntuple_filename");
  
  ePhiMin1_  = iConfig.getParameter< double >("ePhiMin1" ) ;
  ePhiMax1_  = iConfig.getParameter< double >("ePhiMax1" ) ;
  pPhiMin1_  = iConfig.getParameter< double >("pPhiMin1" ) ;
  pPhiMax1_  = iConfig.getParameter< double >("pPhiMax1" ) ;
  PhiMin2_   = iConfig.getParameter< double >("PhiMin2"  ) ;
  PhiMax2_   = iConfig.getParameter< double >("PhiMax2"  ) ;
  r2MinF_    = iConfig.getParameter< double >("r2MinF"   ) ;
  r2MaxF_    = iConfig.getParameter< double >("r2MaxF"   ) ;
  rMinI_     = iConfig.getParameter< double >("rMinI"    ) ;
  rMaxI_     = iConfig.getParameter< double >("rMaxI"    ) ;
  z2MinB_    = iConfig.getParameter< double >("z2MinB"   ) ;
  z2MaxB_    = iConfig.getParameter< double >("z2MaxB"   ) ;
  s_a_phi1B_ = iConfig.getParameter< double >("s_a_phi1B") ;
  s_a_phi1I_ = iConfig.getParameter< double >("s_a_phi1I") ;
  s_a_phi1F_ = iConfig.getParameter< double >("s_a_phi1F") ;
  s_a_phi2B_ = iConfig.getParameter< double >("s_a_phi2B") ;
  s_a_phi2I_ = iConfig.getParameter< double >("s_a_phi2I") ;
  s_a_phi2F_ = iConfig.getParameter< double >("s_a_phi2F") ;
  s_a_zB_    = iConfig.getParameter< double >("s_a_zB"   ) ;
  s_a_rI_    = iConfig.getParameter< double >("s_a_rI"   ) ;
  s_a_rF_    = iConfig.getParameter< double >("s_a_rF"   ) ;
  int use_s_tmp = iConfig.getParameter< int >("use_s"    ) ;
  chainNumber_  = iConfig.getParameter< int >("chainNumber") ;
  nmatch_store_limit_ = iConfig.getParameter< int >("nmatch_store_limit") ;
  
  use_s_ = (use_s_tmp==1) ? true : false ;
  
  s_b_phi1B_ = 1.0/s_a_phi1B_ ;
  s_b_phi1I_ = 1.0/s_a_phi1I_ ;
  s_b_phi1F_ = 1.0/s_a_phi1F_ ;
  s_b_phi2B_ = 1.0/s_a_phi2B_ ;
  s_b_phi2I_ = 1.0/s_a_phi2I_ ;
  s_b_phi2F_ = 1.0/s_a_phi2F_ ;
  s_b_zB_    = 1.0/s_a_zB_    ;
  s_b_rI_    = 1.0/s_a_rI_    ;
  s_b_rF_    = 1.0/s_a_rF_    ;
  
  file_out_ = 0 ;
  file_out_ = new TFile(filename_.c_str(),"RECREATE") ;
  if(0!=file_out_){
    file_out_->cd() ;
    tree_out_ = new TTree("electrons" , "") ;
    tree_par_ = new TTree("parameters", "") ;
    
    // Electron variables
    tree_out_->Branch("el_charge"      , &el_charge_      ) ;
    tree_out_->Branch("el_pt"          , &el_pt_          ) ;
    tree_out_->Branch("el_eta"         , &el_eta_         ) ;
    tree_out_->Branch("el_phi"         , &el_phi_         ) ;
    tree_out_->Branch("el_E"           , &el_E_           ) ;
    tree_out_->Branch("el_nPixelMatch" , &el_nPixelMatch_ ) ;
    tree_out_->Branch("el_matchSuccess", &el_matchSuccess_) ;
    tree_out_->Branch("el_n"           , &el_n_    , "el_n/I"    ) ;
    tree_out_->Branch("el_nPass"       , &el_nPass_, "el_nPass/I") ;
    
    tree_out_->Branch("el_dRz1_neg"    , &el_dRz1_neg_   ) ;
    tree_out_->Branch("el_dRz2_neg"    , &el_dRz2_neg_   ) ;
    tree_out_->Branch("el_dPhi1_neg"   , &el_dPhi1_neg_  ) ;
    tree_out_->Branch("el_dPhi2_neg"   , &el_dPhi2_neg_  ) ;
    tree_out_->Branch("el_dRz1_pos"    , &el_dRz1_pos_   ) ;
    tree_out_->Branch("el_dRz2_pos"    , &el_dRz2_pos_   ) ;
    tree_out_->Branch("el_dPhi1_pos"   , &el_dPhi1_pos_  ) ;
    tree_out_->Branch("el_dPhi2_pos"   , &el_dPhi2_pos_  ) ;
    
    tree_out_->Branch("el_subDet1"     , &el_subDet1_    ) ;
    tree_out_->Branch("el_subDet2"     , &el_subDet2_    ) ;
    tree_out_->Branch("el_matchCharge" , &el_matchCharge_) ;
    tree_out_->Branch("el_n_neg"       , &el_n_neg_      ) ;
    tree_out_->Branch("el_n_pos"       , &el_n_pos_      ) ;
    tree_out_->Branch("el_chainNumber" , &el_chainNumber_) ;
    
    tree_out_->Branch("el_s2_neg"      , &el_s2_neg_     ) ;
    tree_out_->Branch("el_s2_pos"      , &el_s2_pos_     ) ;
    tree_out_->Branch("el_best_s2"     , &el_best_s2_     ) ;
    tree_out_->Branch("el_best_index"  , &el_best_index_  ) ;
    tree_out_->Branch("el_best_charge" , &el_best_charge_ ) ;
    
    // Parameters
    tree_par_->Branch("ePhiMin1", &ePhiMin1_, "ePhiMin1/F") ;
    tree_par_->Branch("ePhiMax1", &ePhiMax1_, "ePhiMax1/F") ;
    tree_par_->Branch("pPhiMin1", &pPhiMin1_, "pPhiMin1/F") ;
    tree_par_->Branch("pPhiMax1", &pPhiMax1_, "pPhiMax1/F") ;
    tree_par_->Branch("PhiMin2" , &PhiMin2_ , "PhiMin2/F" ) ;
    tree_par_->Branch("PhiMax2" , &PhiMax2_ , "PhiMax2/F" ) ;
    tree_par_->Branch("r2MinF"  , &r2MinF_  , "r2MinF/F"  ) ;
    tree_par_->Branch("r2MaxF"  , &r2MaxF_  , "r2MaxF/F"  ) ;
    tree_par_->Branch("rMinI"   , &rMinI_   , "rMinI/F"   ) ;
    tree_par_->Branch("rMaxI"   , &rMaxI_   , "rMaxI/F"   ) ;
    tree_par_->Branch("z2MinB"  , &z2MinB_  , "z2MinB/F"  ) ;
    tree_par_->Branch("z2MaxB"  , &z2MaxB_  , "z2MaxB/F"  ) ;
    
    tree_par_->Branch("npixelmatchcut", &npixelmatchcut_ , "npixelmatchcut/F") ;
    tree_par_->Branch("ncandcut"      , &ncandcut_       , "ncandcut/I"      ) ;
    tree_par_->Branch("chainNumber"   , &chainNumber_    , "chainNumber/I"   ) ;
  }
  if(0!=tree_par_){
    // Store parameters to file
    // Including chainNumber (so store once per constructor)
    tree_par_->Fill() ;
  }
};


HLTElectronPixelMatchFilterAnalysis::~HLTElectronPixelMatchFilterAnalysis(){
  if(file_out_){
    file_out_->Write() ;
    if(tree_out_) delete tree_out_ ;
    if(tree_par_) delete tree_par_ ;
    delete file_out_ ;
  }
  if(el_charge_      ) delete el_charge_ ;
  if(el_pt_          ) delete el_pt_ ;
  if(el_eta_         ) delete el_eta_ ;
  if(el_phi_         ) delete el_phi_ ;
  if(el_E_           ) delete el_E_ ;
  if(el_nPixelMatch_ ) delete el_nPixelMatch_ ;
  if(el_matchSuccess_) delete el_matchSuccess_ ;
  if(el_chainNumber_ ) delete el_chainNumber_ ;
  
  if(el_best_s2_     ) delete el_best_s2_ ;
  if(el_best_charge_ ) delete el_best_charge_ ;
  if(el_best_index_  ) delete el_best_index_ ;
  
  if(el_dRz1_neg_ ) delete el_dRz1_neg_  ;
  if(el_dRz2_neg_ ) delete el_dRz2_neg_  ;
  if(el_dPhi1_neg_) delete el_dPhi1_neg_ ;
  if(el_dPhi2_neg_) delete el_dPhi2_neg_ ;
  if(el_dRz1_pos_ ) delete el_dRz1_pos_  ;
  if(el_dRz2_pos_ ) delete el_dRz2_pos_  ;
  if(el_dPhi1_pos_) delete el_dPhi1_pos_ ;
  if(el_dPhi2_pos_) delete el_dPhi2_pos_ ;
  if(el_subDet1_  ) delete el_subDet1_   ;
  if(el_subDet2_  ) delete el_subDet1_   ;
  if(el_s2_neg_   ) delete el_s2_neg_    ;
  if(el_s2_pos_   ) delete el_s2_pos_    ;
  
  if(el_matchCharge_) delete el_matchCharge_ ;
  
  if(el_n_neg_) delete el_n_neg_ ;
  if(el_n_pos_) delete el_n_pos_ ;
}

float HLTElectronPixelMatchFilterAnalysis::calculate_s2(reco::ElectronSeedCollection::const_iterator it, int charge)
{
  int subDet1 = it->subDet1() ;
  int subDet2 = it->subDet2() ;
  if(charge<0){ // Negative
    if(subDet1==1 && subDet2==1){ // Barrel
      return pow(s_b_phi1B_*it->dPhi1(),2) + pow(s_b_phi2B_*it->dPhi2(),2) +pow(s_b_zB_*it->dRz1(),2) ;
    }
    else if(subDet1==1 && subDet2!=1){ // Intermediate
      return pow(s_b_phi1I_*it->dPhi1(),2) + pow(s_b_phi1I_*it->dPhi2(),2) +pow(s_b_rI_*it->dRz1(),2) ;
    }
    else if(subDet1!=1 && subDet2!=1){ // Forward
      return pow(s_b_phi1F_*it->dPhi1(),2) + pow(s_b_phi1F_*it->dPhi2(),2) +pow(s_b_rF_*it->dRz1(),2) ;
    }
  }
  else{ // Positive
    if(subDet1==1 && subDet2==1){ // Barrel
      return pow(s_b_phi1B_*it->dPhi1(),2) + pow(s_b_phi1B_*it->dPhi2(),2) +pow(s_b_zB_*it->dRz1Pos(),2) ;
    }
    else if(subDet1==1 && subDet2!=1){ // Intermediate
      return pow(s_b_phi1I_*it->dPhi1(),2) + pow(s_b_phi1I_*it->dPhi2(),2) +pow(s_b_rI_*it->dRz1Pos(),2) ;
    }
    else if((subDet1=!1) && (subDet2!=1)){ // Forward
      return pow(s_b_phi1F_*it->dPhi1(),2) + pow(s_b_phi1F_*it->dPhi2(),2) +pow(s_b_rF_*it->dRz1Pos(),2) ;
    }
  }
  return 999 ;
}

// ------------ method called to produce the data  ------------
bool
HLTElectronPixelMatchFilterAnalysis::hltFilter(edm::Event& iEvent, const edm::EventSetup& iSetup, trigger::TriggerFilterObjectWithRefs & filterproduct)
{
  // The filter object
  using namespace trigger;
  if (saveTags()) {
    filterproduct.addCollectionTag(L1IsoCollTag_);
    if (not doIsolated_) filterproduct.addCollectionTag(L1NonIsoCollTag_);
  }

  // Ref to Candidate object to be recorded in filter object
  edm::Ref<reco::RecoEcalCandidateCollection> ref;

  edm::Handle<trigger::TriggerFilterObjectWithRefs> PrevFilterOutput;

  iEvent.getByLabel (candTag_,PrevFilterOutput);

  std::vector<edm::Ref<reco::RecoEcalCandidateCollection> > recoEcalCands;
  PrevFilterOutput->getObjects(TriggerCluster, recoEcalCands);
  if(recoEcalCands.empty()) PrevFilterOutput->getObjects(TriggerPhoton,recoEcalCands);  //we dont know if its type trigger cluster or trigger photon
  
  //get hold of the pixel seed - supercluster association map
  edm::Handle<reco::ElectronSeedCollection> L1IsoSeeds;
  iEvent.getByLabel (L1IsoPixelSeedsTag_,L1IsoSeeds);

  edm::Handle<reco::ElectronSeedCollection> L1NonIsoSeeds;
  if(!doIsolated_){
    iEvent.getByLabel (L1NonIsoPixelSeedsTag_,L1NonIsoSeeds);
  }
  
  el_charge_->clear() ;
  el_pt_ ->clear() ;
  el_eta_->clear() ;
  el_phi_->clear() ;
  el_E_  ->clear() ;
  
  el_nPixelMatch_ ->clear() ;
  el_matchSuccess_->clear() ;
  
  el_n_     = 0 ;
  el_nPass_ = 0 ;
  
  el_dRz1_neg_   ->clear() ;
  el_dRz2_neg_   ->clear() ;
  el_dPhi1_neg_  ->clear() ;
  el_dPhi2_neg_  ->clear() ;
  el_dRz1_pos_   ->clear() ;
  el_dRz2_pos_   ->clear() ;
  el_dPhi1_pos_  ->clear() ;
  el_dPhi2_pos_  ->clear() ;
  el_subDet1_    ->clear() ;
  el_subDet2_    ->clear() ;
  el_matchCharge_->clear() ;
  el_chainNumber_->clear() ;

  el_s2_neg_->clear() ;
  el_s2_pos_->clear() ;
  
  el_best_s2_    ->clear() ;
  el_best_index_ ->clear() ;
  el_best_charge_->clear() ;
  
  // look at all egammas,  check cuts and add to filter object
  int n = 0 ;
  for (unsigned int i=0; i<recoEcalCands.size(); i++) {

    ref = recoEcalCands[i];
    reco::SuperClusterRef recr2 = ref->superCluster();

    int nmatch = 0;
    int nfail  = 0;
    
    // Update TTree contents
    math::XYZTLorentzVector p4 = ref->p4() ;
    el_pt_ ->push_back(p4.pt() ) ;
    el_eta_->push_back(p4.eta()) ;
    el_phi_->push_back(p4.phi()) ;
    el_E_  ->push_back(p4.e()  ) ;
    el_n_++ ;
    
    int n_positive = 0 ;
    int n_negative = 0 ;
    
    el_n_neg_->clear() ;
    el_n_pos_->clear() ;

    std::vector<float> el_dRz1_neg_tmp    ;
    std::vector<float> el_dRz2_neg_tmp    ;
    std::vector<float> el_dPhi1_neg_tmp   ;
    std::vector<float> el_dPhi2_neg_tmp   ;
    std::vector<float> el_s2_neg_tmp      ;
    std::vector<float> el_dRz1_pos_tmp    ;
    std::vector<float> el_dRz2_pos_tmp    ;
    std::vector<float> el_dPhi1_pos_tmp   ;
    std::vector<float> el_dPhi2_pos_tmp   ;
    std::vector<float> el_s2_pos_tmp      ;
    std::vector<int  > el_subDet1_tmp     ;
    std::vector<int  > el_subDet2_tmp     ;
    std::vector<int  > el_matchCharge_tmp ;
    
    int el_best_index_tmp  = -999 ;
    int el_best_charge_tmp =    0 ;
    float el_best_s2_tmp   =  999 ;
    
    for(reco::ElectronSeedCollection::const_iterator it = L1IsoSeeds->begin(); it != L1IsoSeeds->end(); it++){
      edm::RefToBase<reco::CaloCluster> caloCluster = it->caloCluster() ;
      
      reco::SuperClusterRef scRef = caloCluster.castTo<reco::SuperClusterRef>() ;
      if(&(*recr2) ==  &(*scRef)){
        nmatch++;
        int charge_tmp = it->getCharge() ;
        if(charge_tmp>0) n_positive++ ;
        if(charge_tmp<0) n_negative++ ;
        el_matchCharge_tmp.push_back(charge_tmp) ;
        float el_s2_neg = -999 ; 
        float el_s2_pos = -999 ;
        if(nmatch<nmatch_store_limit_ && nmatch_store_limit_>0){
          el_dRz1_neg_tmp .push_back( it->dRz1()     ) ;
          el_dRz2_neg_tmp .push_back( it->dRz2()     ) ;
          el_dPhi1_neg_tmp.push_back( it->dPhi1()    ) ;
          el_dPhi2_neg_tmp.push_back( it->dPhi2()    ) ;
          el_dRz1_pos_tmp .push_back( it->dRz1Pos()  ) ;
          el_dRz2_pos_tmp .push_back( it->dRz2Pos()  ) ;
          el_dPhi1_pos_tmp.push_back( it->dPhi1Pos() ) ;
          el_dPhi2_pos_tmp.push_back( it->dPhi2Pos() ) ;
          el_subDet1_tmp  .push_back( it->subDet1()  ) ;
          el_subDet2_tmp  .push_back( it->subDet2()  ) ;
        }
        if(use_s_){
          el_s2_neg = calculate_s2(it,-1) ;
          el_s2_pos = calculate_s2(it, 1) ;
          if(el_s2_neg<el_best_s2_tmp){
            el_best_s2_tmp = el_s2_neg ;
            el_best_index_tmp = el_dRz1_neg_tmp.size() ;
            el_best_charge_tmp = -1 ;
          }
          if(el_s2_pos<el_best_s2_tmp){
            el_best_s2_tmp = el_s2_pos ;
            el_best_index_tmp = el_dRz1_pos_tmp.size() ;
            el_best_charge_tmp =  1 ;
          }
        }
        if(nmatch<nmatch_store_limit_ && nmatch_store_limit_>0){
          el_s2_neg_tmp   .push_back( el_s2_neg      ) ;
          el_s2_pos_tmp   .push_back( el_s2_pos      ) ;
        }
      }
      else{
        nfail++ ;
      }
    }
    el_dRz1_neg_   ->push_back(el_dRz1_neg_tmp   ) ;
    el_dRz2_neg_   ->push_back(el_dRz2_neg_tmp   ) ;
    el_dPhi1_neg_  ->push_back(el_dPhi1_neg_tmp  ) ;
    el_dPhi2_neg_  ->push_back(el_dPhi2_neg_tmp  ) ;
    el_dRz1_pos_   ->push_back(el_dRz1_pos_tmp   ) ;
    el_dRz2_pos_   ->push_back(el_dRz2_pos_tmp   ) ;
    el_dPhi1_pos_  ->push_back(el_dPhi1_pos_tmp  ) ;
    el_dPhi2_pos_  ->push_back(el_dPhi2_pos_tmp  ) ;
    el_subDet1_    ->push_back(el_subDet1_tmp    ) ;
    el_subDet2_    ->push_back(el_subDet2_tmp    ) ;
    el_s2_neg_     ->push_back(el_s2_neg_tmp     ) ;
    el_s2_pos_     ->push_back(el_s2_pos_tmp     ) ;
    el_matchCharge_->push_back(el_matchCharge_tmp) ;
    el_n_neg_      ->push_back(el_dRz1_neg_tmp.size()) ;
    el_n_pos_      ->push_back(el_dRz1_pos_tmp.size()) ;
    el_best_s2_    ->push_back(el_best_s2_tmp    ) ;
    el_best_index_ ->push_back(el_best_index_tmp ) ;
    el_best_charge_->push_back(el_best_charge_tmp) ;

    if(!doIsolated_){
      for(reco::ElectronSeedCollection::const_iterator it = L1NonIsoSeeds->begin(); it != L1NonIsoSeeds->end(); it++){
        edm::RefToBase<reco::CaloCluster> caloCluster = it->caloCluster() ;
        reco::SuperClusterRef scRef = caloCluster.castTo<reco::SuperClusterRef>() ;
        if(&(*recr2) ==  &(*scRef)) {
          nmatch++;
        }
      }

    }//end if(!doIsolated_)
    bool accept_electron = false ;
    if ( nmatch >= npixelmatchcut_) {
      n++;
      filterproduct.addObject(TriggerCluster, ref);
      accept_electron = true ;
    }
    el_nPixelMatch_->push_back(nmatch) ;
    el_matchSuccess_->push_back(accept_electron) ;
    el_chainNumber_->push_back(chainNumber_) ;
    
    float charge = 0 ;
    if(n_positive+n_negative>0){
      if(n_positive>n_negative){
        charge =  1.0*n_positive/(n_positive+n_negative) ;
      }
      else if(n_positive<n_negative){
        charge = -1.0*n_negative/(n_positive+n_negative) ;
      }
    }
    el_charge_->push_back(charge) ;
    
  }//end of loop over candidates
  el_nPass_ = n ;
  
  // filter decision
  bool accept(n>=ncandcut_);
  accept = true ;
  
  if(tree_out_){
    tree_out_->Fill() ;
  }
  
  return accept;
}


