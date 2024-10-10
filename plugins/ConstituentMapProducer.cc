#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/StreamID.h"

#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"

#include "DataFormats/Candidate/interface/CandidateFwd.h"

#include "CommonTools/Utils/interface/StringCutObjectSelector.h"
#include "DataFormats/NanoAOD/interface/FlatTable.h"

class ConstituentMapProducer : public edm::stream::EDProducer<> {
public:
  explicit ConstituentMapProducer(const edm::ParameterSet &);
  ~ConstituentMapProducer() override;

  static void fillDescriptions(edm::ConfigurationDescriptions &descriptions);

private:
  void produce(edm::Event &, const edm::EventSetup &) override;

  edm::EDGetTokenT<std::vector<reco::CandidatePtr>> cand_token_;
  edm::EDGetTokenT<std::vector<reco::CandidatePtr>> gencand_token_;


};

ConstituentMapProducer::ConstituentMapProducer(const edm::ParameterSet &iConfig)
    : 
      cand_token_(consumes< std::vector <reco::CandidatePtr> >(iConfig.getParameter<edm::InputTag>("candidates"))),
      gencand_token_(consumes< std::vector <reco::CandidatePtr> >(iConfig.getParameter<edm::InputTag>("gencandidates")))
{
  produces<nanoaod::FlatTable>("PFCandsMatch");
}

ConstituentMapProducer::~ConstituentMapProducer() {}

void ConstituentMapProducer::produce(edm::Event &iEvent, const edm::EventSetup &iSetup) {

  std::vector<int> genIdx;

  auto cands_ = iEvent.getHandle(cand_token_);
  auto gencands_ = iEvent.getHandle(gencand_token_);

  for (const auto& cand : *cands_) { 
     
     if (cand->charge()==0) {
         genIdx.push_back(-1); 
         continue;
     }
   
     int matchGenIdx = -1;
     double minDR2 = std::numeric_limits<double>::max();
     for (unsigned i = 0; i < gencands_->size(); ++i) {
       
         if (std::find( genIdx.begin(), genIdx.end(), i ) != genIdx.end()) continue;

         auto gencand = gencands_->at(i);

         double DR2 = reco::deltaR2(*cand, *gencand);
         if (DR2 > 0.0004) continue;

         double relPt = cand->pt() / gencand->pt();
         if (relPt < 0.8 || relPt > 1.2) continue;

         if (DR2 < minDR2) {
            minDR2 = DR2;
            matchGenIdx = i;
         }
          
     }

     genIdx.push_back(matchGenIdx);
  }

  auto candTable = std::make_unique<nanoaod::FlatTable>(cands_->size(), "PFCandsMatch", false);
  candTable->addColumn<int>("genIdx", genIdx, "Index of the matched (charged only) gen candidates", nanoaod::FlatTable::IntColumn);
  iEvent.put(std::move(candTable), "PFCandsMatch");
}


void ConstituentMapProducer::fillDescriptions(edm::ConfigurationDescriptions &descriptions) {
  edm::ParameterSetDescription desc;
  desc.add<edm::InputTag>("candidates", edm::InputTag("packedPFCandidates"));
  desc.add<edm::InputTag>("gencandidates", edm::InputTag("packedPFCandidates"));
  descriptions.addWithDefaultLabel(desc);
}

typedef ConstituentMapProducer GenConstituentMapProducer;

DEFINE_FWK_MODULE(GenConstituentMapProducer);
