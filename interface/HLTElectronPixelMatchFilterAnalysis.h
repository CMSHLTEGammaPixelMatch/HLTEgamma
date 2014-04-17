#ifndef HLTElectronPixelMatchFilterAnalysis_h
#define HLTElectronPixelMatchFilterAnalysis_h

/** \class HLTElectronPixelMatchFilterAnalysis
 *
 *  \author Aidan Randle-Conde (ULB)
 *
 */

#include "HLTrigger/HLTcore/interface/HLTFilter.h"

#include "DataFormats/EgammaReco/interface/ElectronSeed.h"

#include <TFile.h>
#include <TTree.h>

#include <vector>

//
// class decleration
//

class HLTElectronPixelMatchFilterAnalysis : public HLTFilter {

   public:
      explicit HLTElectronPixelMatchFilterAnalysis(const edm::ParameterSet&);
      ~HLTElectronPixelMatchFilterAnalysis();
      virtual bool hltFilter(edm::Event&, const edm::EventSetup&, trigger::TriggerFilterObjectWithRefs & filterproduct);

   private:
      edm::InputTag candTag_;     // input tag identifying product contains filtered egammas

      edm::InputTag L1IsoPixelSeedsTag_; // input tag for the pixel seed - supercluster map
      //edm::InputTag L1IsoPixelmapendcapTag_; // input tag for the pixel seed - supercluster map

      edm::InputTag L1NonIsoPixelSeedsTag_; // input tag for the pixel seed - supercluster map
      //edm::InputTag L1NonIsoPixelmapendcapTag_; // input tag for the pixel seed - supercluster map

      double npixelmatchcut_;     // number of pixelmatch hits
      int    ncandcut_;           // number of electrons required
      
      bool doIsolated_;
      edm::InputTag L1IsoCollTag_; 
      edm::InputTag L1NonIsoCollTag_;
      
      float calculate_s2(reco::ElectronSeedCollection::const_iterator, int) ;
      
      // Information to be saved to file
      TFile* file_out_ ;
      TTree* tree_out_ ;
      TTree* tree_par_ ;
      std::string filename_ ;
      
      // Limit on nmatch to store to protect from excessive RAM usage
      int nmatch_store_limit_;
      
      // Charge information
      std::vector<float>* el_charge_ ;
      std::vector<float>* el_pt_ ;
      std::vector<float>* el_eta_ ;
      std::vector<float>* el_phi_ ;
      std::vector<float>* el_E_ ;
      std::vector<int  >* el_nPixelMatch_ ;
      std::vector<int  >* el_matchSuccess_ ;
      int el_n_ ;
      int el_nPass_ ;
      
      // Matching variables
      std::vector< std::vector<float> >* el_dRz1_neg_  ;
      std::vector< std::vector<float> >* el_dRz2_neg_  ;
      std::vector< std::vector<float> >* el_dPhi1_neg_ ;
      std::vector< std::vector<float> >* el_dPhi2_neg_ ;
      std::vector< std::vector<float> >* el_dRz1_pos_  ;
      std::vector< std::vector<float> >* el_dRz2_pos_  ;
      std::vector< std::vector<float> >* el_dPhi1_pos_ ;
      std::vector< std::vector<float> >* el_dPhi2_pos_ ;
      std::vector< std::vector<float> >* el_s2_neg_    ;
      std::vector< std::vector<float> >* el_s2_pos_    ;
      std::vector< std::vector<int  > >* el_subDet1_   ;
      std::vector< std::vector<int  > >* el_subDet2_   ;
      std::vector< std::vector<int  > >* el_matchCharge_ ;
      
      std::vector< int >* el_n_neg_ ;
      std::vector< int >* el_n_pos_ ;
      std::vector< int >* el_chainNumber_ ;
      
      // Parameters
      float ePhiMin1_ ;
      float ePhiMax1_ ;
      float pPhiMin1_ ;
      float pPhiMax1_ ;
      float PhiMin2_  ;
      float PhiMax2_  ;
      float r2MinF_   ;
      float r2MaxF_   ;
      float rMinI_    ;
      float rMaxI_    ;
      float z2MinB_   ;
      float z2MaxB_   ;
      
      // S parameter values
      // Divide by s_a_
      float s_a_phi1B_ ;
      float s_a_phi1I_ ;
      float s_a_phi1F_ ;
      float s_a_phi2B_ ;
      float s_a_phi2I_ ;
      float s_a_phi2F_ ;
      float s_a_zB_ ;
      float s_a_rI_ ;
      float s_a_rF_ ;
      
      // Multiply by s_b_ (marginally quicker)
      float s_b_phi1B_ ;
      float s_b_phi1I_ ;
      float s_b_phi1F_ ;
      float s_b_phi2B_ ;
      float s_b_phi2I_ ;
      float s_b_phi2F_ ;
      float s_b_zB_ ;
      float s_b_rI_ ;
      float s_b_rF_ ;
      
      bool use_s_ ;
      int chainNumber_ ;
      
      std::vector< float >* el_best_s2_     ;
      std::vector< int   >* el_best_index_  ;
      std::vector< int   >* el_best_charge_ ;
};

#endif //HLTElectronPixelMatchFilter_h


