// Writes a NanoAOD FlatTable for the IFN jet collection produced by
// PhysicsTools/IFNFlavour/JetFlavourClusteringIFN. Columns: pt, eta, phi, mass,
// partonFlavourIFN, genJetIdx.
//
// Inputs:
//   src         -- reco::BasicJetCollection (e.g. genJetFlavourAssociationIFN:ifnJets)
//   flavour     -- std::vector<int> parallel to src (e.g.
//                  genJetFlavourAssociationIFN:ifnJetNetFlavourOverall)
//   genIndex    -- std::vector<int> parallel to src, index of the matched gen jet
//                  (e.g. genJetFlavourAssociationIFN:IFNJetgenIndex); -1 if none.

#include <memory>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/Utilities/interface/Exception.h"

#include "DataFormats/JetReco/interface/BasicJet.h"
#include "DataFormats/JetReco/interface/BasicJetCollection.h"
#include "DataFormats/NanoAOD/interface/FlatTable.h"

class IFNJetTableProducer : public edm::stream::EDProducer<> {
public:
  explicit IFNJetTableProducer(const edm::ParameterSet& iConfig)
      : name_(iConfig.getParameter<std::string>("name")),
        doc_(iConfig.getParameter<std::string>("doc")),
        src_(consumes<reco::BasicJetCollection>(iConfig.getParameter<edm::InputTag>("src"))),
        flavour_(consumes<std::vector<int>>(iConfig.getParameter<edm::InputTag>("flavour"))),
        genIndex_(consumes<std::vector<int>>(iConfig.getParameter<edm::InputTag>("genIndex"))),
        precision_(iConfig.getParameter<int>("precision")) {
    produces<nanoaod::FlatTable>();
  }

  ~IFNJetTableProducer() override = default;

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
    edm::ParameterSetDescription desc;
    desc.add<std::string>("name", "IFNJet")->setComment("name of the FlatTable / branch prefix");
    desc.add<std::string>("doc", "IFN clustered jets (parton or hadron level, see IFNFlavour)");
    desc.add<edm::InputTag>("src", edm::InputTag("genJetFlavourAssociationIFN", "ifnJets"))
        ->setComment("IFN jet BasicJet collection");
    desc.add<edm::InputTag>("flavour", edm::InputTag("genJetFlavourAssociationIFN", "ifnJetNetFlavourOverall"))
        ->setComment("per-IFN-jet parton flavour, parallel to src");
    desc.add<edm::InputTag>("genIndex", edm::InputTag("genJetFlavourAssociationIFN", "IFNJetgenIndex"))
        ->setComment("per-IFN-jet matched gen-jet index, parallel to src; -1 if unmatched");
    desc.add<int>("precision", 10)->setComment("mantissa precision for float columns");
    descriptions.addWithDefaultLabel(desc);
  }

private:
  void produce(edm::Event& iEvent, edm::EventSetup const&) override {
    edm::Handle<reco::BasicJetCollection> jets;
    iEvent.getByToken(src_, jets);

    edm::Handle<std::vector<int>> flavour;
    iEvent.getByToken(flavour_, flavour);

    edm::Handle<std::vector<int>> genIndex;
    iEvent.getByToken(genIndex_, genIndex);

    if (jets->size() != flavour->size() || jets->size() != genIndex->size())
      throw cms::Exception("IFNJetTableProducer")
          << "IFN jet collection (" << jets->size() << "), flavour (" << flavour->size()
          << ") and genIndex (" << genIndex->size() << ") have different sizes";

    const unsigned n = jets->size();
    std::vector<float> pt, eta, phi, mass;
    pt.reserve(n);
    eta.reserve(n);
    phi.reserve(n);
    mass.reserve(n);
    for (const auto& j : *jets) {
      pt.push_back(j.pt());
      eta.push_back(j.eta());
      phi.push_back(j.phi());
      mass.push_back(j.mass());
    }

    auto tab = std::make_unique<nanoaod::FlatTable>(n, name_, false, false);
    tab->setDoc(doc_);
    tab->addColumn<float>("pt", pt, "pt", nanoaod::FlatTable::FloatColumn, precision_);
    tab->addColumn<float>("eta", eta, "eta", nanoaod::FlatTable::FloatColumn, precision_);
    tab->addColumn<float>("phi", phi, "phi", nanoaod::FlatTable::FloatColumn, precision_);
    tab->addColumn<float>("mass", mass, "mass", nanoaod::FlatTable::FloatColumn, precision_);
    tab->addColumn<int>("partonFlavourIFN", *flavour, "IFN parton flavour (net)", nanoaod::FlatTable::IntColumn);
    tab->addColumn<int>("genJetIdx", *genIndex, "index of matched gen jet (-1 if unmatched)", nanoaod::FlatTable::IntColumn);
    iEvent.put(std::move(tab));
  }

  const std::string name_;
  const std::string doc_;
  const edm::EDGetTokenT<reco::BasicJetCollection> src_;
  const edm::EDGetTokenT<std::vector<int>> flavour_;
  const edm::EDGetTokenT<std::vector<int>> genIndex_;
  const int precision_;
};

DEFINE_FWK_MODULE(IFNJetTableProducer);
