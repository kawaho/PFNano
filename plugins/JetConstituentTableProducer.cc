#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/StreamID.h"

#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "RecoVertex/VertexTools/interface/VertexDistance3D.h"
#include "RecoVertex/VertexTools/interface/VertexDistanceXY.h"
#include "RecoVertex/VertexPrimitives/interface/ConvertToFromReco.h"
#include "RecoVertex/VertexPrimitives/interface/VertexState.h"

#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"

#include "DataFormats/Candidate/interface/CandidateFwd.h"

#include "RecoBTag/FeatureTools/interface/TrackInfoBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"

#include "DataFormats/BTauReco/interface/TrackIPTagInfo.h"
#include "DataFormats/BTauReco/interface/SecondaryVertexTagInfo.h"
#include "RecoBTag/FeatureTools/interface/deep_helpers.h"
#include "DataFormats/Candidate/interface/VertexCompositePtrCandidate.h"
using namespace btagbtvdeep;

#include "CommonTools/Utils/interface/StringCutObjectSelector.h"
#include "DataFormats/NanoAOD/interface/FlatTable.h"

template<typename T>
class JetConstituentTableProducer : public edm::stream::EDProducer<> {
public:
  explicit JetConstituentTableProducer(const edm::ParameterSet &);
  ~JetConstituentTableProducer() override;

  static void fillDescriptions(edm::ConfigurationDescriptions &descriptions);

private:
  void produce(edm::Event &, const edm::EventSetup &) override;

  bool isBAncestor(const reco::Candidate &particle);

  typedef reco::VertexCollection VertexCollection;
  //=====
  typedef reco::VertexCompositePtrCandidateCollection SVCollection;

  const std::string name_;
  const std::string name_in_evt;
  const std::string nameSV_;
  const std::string idx_name_;
  const std::string sv_idx_name_;
  const std::string idx_nameSV_;
  const bool readBtag_;
  const double jet_radius_;

  edm::EDGetTokenT<edm::View<T>> jet_token_;
  edm::EDGetTokenT<VertexCollection> vtx_token_;
  edm::EDGetTokenT<reco::CandidateView> cand_token_;
//  edm::EDGetTokenT<reco::CandidateView> pfcand_token_;
  edm::EDGetTokenT<SVCollection> sv_token_;

  edm::Handle<VertexCollection> vtxs_;
  edm::Handle<reco::CandidateView> cands_;
//  edm::Handle<reco::CandidateView> pfcands_;
  edm::Handle<SVCollection> svs_;
  edm::ESHandle<TransientTrackBuilder> track_builder_;

  const reco::Vertex *pv_ = nullptr;
  
};

//
// constructors and destructor
//
template< typename T>
JetConstituentTableProducer<T>::JetConstituentTableProducer(const edm::ParameterSet &iConfig)
    : name_(iConfig.getParameter<std::string>("name")),
      name_in_evt(iConfig.getParameter<std::string>("name_in_evt")),
      nameSV_(iConfig.getParameter<std::string>("nameSV")),
      idx_name_(iConfig.getParameter<std::string>("idx_name")),
      sv_idx_name_(iConfig.getParameter<std::string>("sv_idx_name")),
      idx_nameSV_(iConfig.getParameter<std::string>("idx_nameSV")),
      readBtag_(iConfig.getParameter<bool>("readBtag")),
      jet_radius_(iConfig.getParameter<double>("jet_radius")),
      jet_token_(consumes<edm::View<T>>(iConfig.getParameter<edm::InputTag>("jets"))),
      vtx_token_(consumes<VertexCollection>(iConfig.getParameter<edm::InputTag>("vertices"))),
      cand_token_(consumes<reco::CandidateView>(iConfig.getParameter<edm::InputTag>("candidates"))),
//      pfcand_token_(consumes<reco::CandidateView>(iConfig.getParameter<edm::InputTag>("pfcandidates"))),
      sv_token_(consumes<SVCollection>(iConfig.getParameter<edm::InputTag>("secondary_vertices"))){

  produces<nanoaod::FlatTable>(name_);
  produces<nanoaod::FlatTable>(nameSV_);
  produces<std::vector<reco::CandidatePtr>>(name_in_evt);
}

template< typename T>
JetConstituentTableProducer<T>::~JetConstituentTableProducer() {}

template< typename T>
bool JetConstituentTableProducer<T>::isBAncestor(const reco::Candidate &particle) {
   //particle is already the B
   if (((int)((abs(particle.pdgId()) / 100) % 10) == 5) || ((int)((abs(particle.pdgId()) / 1000) % 10) == 5)) return true;
   //if (abs(particle.pdgId()) == 5) return true;

   //otherwise loop on mothers, if any and return true if the B ancestor is found
   for (size_t i = 0; i < particle.numberOfMothers(); i++) {
      if (isBAncestor(*particle.mother(i))) return true;
                                        
   }

   return false;
}


template< typename T>
void JetConstituentTableProducer<T>::produce(edm::Event &iEvent, const edm::EventSetup &iSetup) {
  // elements in all these collections must have the same order!
  auto outCands = std::make_unique<std::vector<reco::CandidatePtr>>();
  auto outSVs = std::make_unique<std::vector<const reco::VertexCompositePtrCandidate *>> ();
  std::vector<int> jetIdx_pf, jetIdx_sv, pfcandIdx, svcandIdx, svIdx;
  std::vector<int> from_b;
  // PF Cands
  std::vector<float> btagEtaRel, btagPtRatio, btagPParRatio, btagSip3dVal, btagSip3dSig, btagJetDistVal, btagTrackDeltaR, btagSip2dVal, btagSip2dSig, btagTrackJetDecayLen, cand_pt;

  // Secondary vertices
  std::vector<float> sv_mass, sv_pt, sv_ntracks, sv_chi2, sv_normchi2, sv_dxy, sv_dxysig, sv_d3d, sv_d3dsig, sv_costhetasvpv;
  std::vector<float> sv_ptrel, sv_phirel, sv_deltaR, sv_enratio;

  auto jets = iEvent.getHandle(jet_token_);
  iEvent.getByToken(vtx_token_, vtxs_);
  iEvent.getByToken(cand_token_, cands_);
//  iEvent.getByToken(pfcand_token_, pfcands_);
  iEvent.getByToken(sv_token_, svs_);

  if(readBtag_){
    iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder", track_builder_);
  }

  for (unsigned i_jet = 0; i_jet < jets->size(); ++i_jet) {
    const auto &jet = jets->at(i_jet);
    math::XYZVector jet_dir = jet.momentum().Unit();
    GlobalVector jet_ref_track_dir(jet.px(), jet.py(), jet.pz());
    VertexDistance3D vdist;

    pv_ = &vtxs_->at(0);
    
    //////////////////////
    // Secondary Vertices
    std::vector<const reco::VertexCompositePtrCandidate *> jetSVs;
    std::vector<const reco::VertexCompositePtrCandidate *> allSVs;  
    for (const auto &sv : *svs_) {
      // Factor in cuts in NanoAOD for indexing
      Measurement1D dl= vdist.distance(vtxs_->front(), VertexState(RecoVertex::convertPos(sv.position()),RecoVertex::convertError(sv.error())));
      if(dl.significance() > 3){
        allSVs.push_back(&sv);
      }
      if (reco::deltaR2(sv, jet) < jet_radius_ * jet_radius_) {
        jetSVs.push_back(&sv);
      }
    }
    // sort by dxy significance
    std::sort(jetSVs.begin(),
              jetSVs.end(),
              [&](const reco::VertexCompositePtrCandidate *sva, const reco::VertexCompositePtrCandidate *svb) {
                return sv_vertex_comparator(*sva, *svb, *pv_);
              });

    for (const auto &sv : jetSVs) {
      // auto svPtrs = svs_->ptrs();
      auto svInNewList = std::find(allSVs.begin(), allSVs.end(), sv );
      if (svInNewList == allSVs.end()) {
        // continue;
        svIdx.push_back(-1);
      } else {
        svIdx.push_back(svInNewList - allSVs.begin());
      }
      outSVs->push_back(sv);
      jetIdx_sv.push_back(i_jet);
      if (readBtag_ && !vtxs_->empty()) {

	// Jet independent
        sv_mass.push_back(sv->mass());
        sv_pt.push_back(sv->pt());

        sv_ntracks.push_back(sv->numberOfDaughters());
        sv_chi2.push_back(sv->vertexChi2());
        sv_normchi2.push_back(catch_infs_and_bound(sv->vertexChi2() / sv->vertexNdof(), 1000, -1000, 1000));
        const auto& dxy_meas = vertexDxy(*sv, *pv_);
        sv_dxy.push_back(dxy_meas.value());
        sv_dxysig.push_back(catch_infs_and_bound(dxy_meas.value() / dxy_meas.error(), 0, -1, 800));
        const auto& d3d_meas = vertexD3d(*sv, *pv_);
        sv_d3d.push_back(d3d_meas.value());
        sv_d3dsig.push_back(catch_infs_and_bound(d3d_meas.value() / d3d_meas.error(), 0, -1, 800));
        sv_costhetasvpv.push_back(vertexDdotP(*sv, *pv_));
        // Jet related
        sv_ptrel.push_back(sv->pt() / jet.pt());
        sv_phirel.push_back(reco::deltaPhi(*sv, jet));
        sv_deltaR.push_back(catch_infs_and_bound(std::fabs(reco::deltaR(*sv, jet_dir)) - 0.5, 0, -2, 0));
        sv_enratio.push_back(sv->energy() / jet.energy());
      }
    }

    // PF Cands    
    std::vector<reco::CandidatePtr> const & daughters = jet.daughterPtrVector();

    for (const auto &cand : daughters) {
      auto candPtrs = cands_->ptrs();
      auto candInNewList = std::find( candPtrs.begin(), candPtrs.end(), cand );
      if ( candInNewList == candPtrs.end() ) {
        //std::cout << "Cannot find candidate : " << cand.id() << ", " << cand.key() << ", pt = " << cand->pt() << std::endl;
        continue;
      }
      outCands->push_back(cand);
      jetIdx_pf.push_back(i_jet);
      pfcandIdx.push_back(candInNewList - candPtrs.begin());
      cand_pt.push_back(cand->pt());

      if (!readBtag_) {
         if (isBAncestor(*cand)) {
           from_b.push_back(1);
         } else {
           from_b.push_back(0);
         }
      }

      if (readBtag_ && !vtxs_->empty()) {
         svcandIdx.push_back(-1);
         if (cand->charge() != 0) {
           for (size_t jetSVs_index = 0; jetSVs_index < jetSVs.size(); ++jetSVs_index) {
             std::vector<reco::CandidatePtr> const & sv_daughters = jetSVs[jetSVs_index]->daughterPtrVector();
             auto candInNewList = std::find( sv_daughters.begin(), sv_daughters.end(), cand );
             if ( candInNewList == sv_daughters.end() ) {
               continue;
             } else {
                // std::cout << "Found SV candidate : " << cand.id() << ", " << cand.key() << ", pt = " << cand->pt() << std::endl;
                svcandIdx.back() = jetSVs_index;
             }
           }
        }
        if ( cand.isNull() ) continue;
        auto const *packedCand = dynamic_cast <pat::PackedCandidate const *>(cand.get());
        if ( packedCand == nullptr ) continue;
        if ( packedCand && packedCand->hasTrackDetails()){
          btagbtvdeep::TrackInfoBuilder trkinfo(track_builder_);
          trkinfo.buildTrackInfo(&(*packedCand), jet_dir, jet_ref_track_dir, vtxs_->at(0));
          btagEtaRel.push_back(trkinfo.getTrackEtaRel());
          btagPtRatio.push_back(trkinfo.getTrackPtRatio());
          btagPParRatio.push_back(trkinfo.getTrackPParRatio());
          btagSip2dVal.push_back(trkinfo.getTrackSip2dVal());
          btagSip2dSig.push_back(trkinfo.getTrackSip2dSig());
          btagSip3dVal.push_back(trkinfo.getTrackSip3dVal());
          btagSip3dSig.push_back(trkinfo.getTrackSip3dSig());
          btagJetDistVal.push_back(trkinfo.getTrackJetDistVal());
          btagTrackDeltaR.push_back(trkinfo.getTrackDeltaR());
          btagTrackJetDecayLen.push_back(trkinfo.getTrackJetDecayLen());

        } else {
                btagEtaRel.push_back(0);
                btagPtRatio.push_back(0);
                btagPParRatio.push_back(0);
                btagSip3dVal.push_back(0);
                btagSip3dSig.push_back(0);
                btagSip2dVal.push_back(0);
                btagSip2dSig.push_back(0);
                btagJetDistVal.push_back(0);
                btagTrackDeltaR.push_back(0);
                btagTrackJetDecayLen.push_back(0);
        }
      }
    }  // end jet loop
  }

  auto candTable = std::make_unique<nanoaod::FlatTable>(outCands->size(), name_, false);
  // We fill from here only stuff that cannot be created with the SimpleFlatTableProducer
  candTable->addColumn<int>(idx_name_, pfcandIdx, "Index in the candidate list", nanoaod::FlatTable::IntColumn);
  candTable->addColumn<int>("jetIdx", jetIdx_pf, "Index of the parent jet", nanoaod::FlatTable::IntColumn);
  if (readBtag_) {
    candTable->addColumn<int>(sv_idx_name_, svcandIdx, "Index of candidate in the jetSVs list", nanoaod::FlatTable::IntColumn);
    candTable->addColumn<float>("pt", cand_pt, "pt", nanoaod::FlatTable::FloatColumn, 10);  // to check matchind down the line
    candTable->addColumn<float>("btagEtaRel", btagEtaRel, "btagEtaRel", nanoaod::FlatTable::FloatColumn, 10);
    candTable->addColumn<float>("btagPtRatio", btagPtRatio, "btagPtRatio", nanoaod::FlatTable::FloatColumn, 10);
    candTable->addColumn<float>("btagPParRatio", btagPParRatio, "btagPParRatio", nanoaod::FlatTable::FloatColumn, 10);
    candTable->addColumn<float>("btagSip3dVal", btagSip3dVal, "btagSip3dVal", nanoaod::FlatTable::FloatColumn, 10);
    candTable->addColumn<float>("btagSip3dSig", btagSip3dSig, "btagSip3dSig", nanoaod::FlatTable::FloatColumn, 10);
    candTable->addColumn<float>("btagSip2dVal", btagSip2dVal, "btagSip2dVal", nanoaod::FlatTable::FloatColumn, 10);
    candTable->addColumn<float>("btagSip2dSig", btagSip2dSig, "btagSip2dSig", nanoaod::FlatTable::FloatColumn, 10);
    candTable->addColumn<float>("btagJetDistVal", btagJetDistVal, "btagJetDistVal", nanoaod::FlatTable::FloatColumn, 10);
    candTable->addColumn<float>("btagTrackDeltaR", btagTrackDeltaR, "btagTrackDeltaR", nanoaod::FlatTable::FloatColumn, 10);
    candTable->addColumn<float>("btagTrackJetDecayLen", btagTrackJetDecayLen, "btagTrackJetDecayLen", nanoaod::FlatTable::FloatColumn, 10);
  } else {
    candTable->addColumn<int>("isbProduct", from_b, "True if ancestor is a b hadron", nanoaod::FlatTable::IntColumn);
  }
  iEvent.put(std::move(candTable), name_);

   // SV table
  auto svTable = std::make_unique<nanoaod::FlatTable>(outSVs->size(), nameSV_, false);
  // We fill from here only stuff that cannot be created with the SimpleFlatTnameableProducer
  svTable->addColumn<int>("jetIdx", jetIdx_sv, "Index of the parent jet", nanoaod::FlatTable::IntColumn);
  svTable->addColumn<int>(idx_nameSV_, svIdx, "Index in the SV list", nanoaod::FlatTable::IntColumn);
  if (readBtag_) {
    svTable->addColumn<float>("mass", sv_mass, "SV mass", nanoaod::FlatTable::FloatColumn, 10);
    svTable->addColumn<float>("pt", sv_pt, "SV pt", nanoaod::FlatTable::FloatColumn, 10);
    svTable->addColumn<float>("ntracks", sv_ntracks, "Number of tracks associated to SV", nanoaod::FlatTable::FloatColumn, 10);
    svTable->addColumn<float>("chi2", sv_chi2, "chi2", nanoaod::FlatTable::FloatColumn, 10);
    svTable->addColumn<float>("normchi2", sv_normchi2, "chi2/ndof", nanoaod::FlatTable::FloatColumn, 10);
    svTable->addColumn<float>("dxy", sv_dxy, "", nanoaod::FlatTable::FloatColumn, 10);
    svTable->addColumn<float>("dxysig", sv_dxysig, "", nanoaod::FlatTable::FloatColumn, 10);
    svTable->addColumn<float>("d3d", sv_d3d, "", nanoaod::FlatTable::FloatColumn, 10);
    svTable->addColumn<float>("d3dsig", sv_d3dsig, "", nanoaod::FlatTable::FloatColumn, 10);
    svTable->addColumn<float>("costhetasvpv", sv_costhetasvpv, "", nanoaod::FlatTable::FloatColumn, 10);
    // Jet related
    svTable->addColumn<float>("phirel", sv_phirel, "DeltaPhi(sv, jet)", nanoaod::FlatTable::FloatColumn, 10);
    svTable->addColumn<float>("ptrel", sv_ptrel, "pT relative to parent jet", nanoaod::FlatTable::FloatColumn, 10);
    svTable->addColumn<float>("deltaR", sv_deltaR, "dR from parent jet", nanoaod::FlatTable::FloatColumn, 10);
    svTable->addColumn<float>("enration", sv_enratio, "energy relative to parent jet", nanoaod::FlatTable::FloatColumn, 10);
  }
  iEvent.put(std::move(svTable), nameSV_);

  iEvent.put(std::move(outCands), name_in_evt);

}

template< typename T>
void JetConstituentTableProducer<T>::fillDescriptions(edm::ConfigurationDescriptions &descriptions) {
  edm::ParameterSetDescription desc;
  desc.add<std::string>("name", "JetPFCands");
  desc.add<std::string>("name_in_evt", "JetPFCands");
  desc.add<std::string>("nameSV", "JetSV");
  desc.add<std::string>("idx_name", "candIdx");
  desc.add<std::string>("sv_idx_name", "svcandIdx");
  desc.add<std::string>("idx_nameSV", "svIdx");
  desc.add<double>("jet_radius", true);
  desc.add<bool>("readBtag", true);
  desc.add<edm::InputTag>("jets", edm::InputTag("slimmedJetsAK8"));
  desc.add<edm::InputTag>("vertices", edm::InputTag("offlineSlimmedPrimaryVertices"));
  desc.add<edm::InputTag>("candidates", edm::InputTag("packedPFCandidates"));
//  desc.add<edm::InputTag>("pfcandidates", edm::InputTag("packedPFCandidates"));
  desc.add<edm::InputTag>("secondary_vertices", edm::InputTag("slimmedSecondaryVertices"));
  descriptions.addWithDefaultLabel(desc);
}

typedef JetConstituentTableProducer<pat::Jet> PatJetConstituentTableProducer;
typedef JetConstituentTableProducer<reco::GenJet> GenJetConstituentTableProducer;

DEFINE_FWK_MODULE(PatJetConstituentTableProducer);
DEFINE_FWK_MODULE(GenJetConstituentTableProducer);
