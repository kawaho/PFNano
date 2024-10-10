import FWCore.ParameterSet.Config as cms
from  PhysicsTools.NanoAOD.common_cff import *

def addPFCands(process, runOnMC=False, allPF = False, onlyAK4=False, onlyAK8=False):
    process.customizedPFCandsTask = cms.Task( )
    process.schedule.associate(process.customizedPFCandsTask)

    process.finalJetsAK8Constituents = cms.EDProducer("PatJetConstituentPtrSelector",
                                            src = cms.InputTag("finalJetsAK8"),
                                            cut = cms.string("")
                                            )
    process.finalJetsAK4Constituents = cms.EDProducer("PatJetConstituentPtrSelector",
                                            src = cms.InputTag("finalJets"),
                                            cut = cms.string("")
                                            )
    if allPF:
        candInput = cms.InputTag("packedPFCandidates")
    elif onlyAK4:
        candList = cms.VInputTag(cms.InputTag("finalJetsAK4Constituents", "constituents"))
        process.customizedPFCandsTask.add(process.finalJetsAK4Constituents)
        process.finalJetsConstituents = cms.EDProducer("PackedCandidatePtrMerger", src = candList, skipNulls = cms.bool(True), warnOnSkip = cms.bool(True))
        candInput = cms.InputTag("finalJetsConstituents")
    elif onlyAK8:
        candList = cms.VInputTag(cms.InputTag("finalJetsAK8Constituents", "constituents"))
        process.customizedPFCandsTask.add(process.finalJetsAK8Constituents)
        process.finalJetsConstituents = cms.EDProducer("PackedCandidatePtrMerger", src = candList, skipNulls = cms.bool(True), warnOnSkip = cms.bool(True))
        candInput = cms.InputTag("finalJetsConstituents")
    else:
        candList = cms.VInputTag(cms.InputTag("finalJetsAK4Constituents", "constituents"), cms.InputTag("finalJetsAK8Constituents", "constituents"))
        process.customizedPFCandsTask.add(process.finalJetsAK4Constituents)
        process.customizedPFCandsTask.add(process.finalJetsAK8Constituents)
        process.finalJetsConstituents = cms.EDProducer("PackedCandidatePtrMerger", src = candList, skipNulls = cms.bool(True), warnOnSkip = cms.bool(True))
        candInput = cms.InputTag("finalJetsConstituents")
    process.customConstituentsExtTable = cms.EDProducer("SimpleCandidateFlatTableProducer",
                                                        src = candInput,
                                                        cut = cms.string(""), #we should not filter after pruning
                                                        name = cms.string("PFCands"),
                                                        doc = cms.string("interesting particles from AK4 and AK8 jets"),
                                                        singleton = cms.bool(False), # the number of entries is variable
                                                        extension = cms.bool(False), # this is the extension table for the AK8 constituents
                                                        variables = cms.PSet(CandVars,
                                                            puppiWeight = Var("puppiWeight()", float, doc="Puppi weight",precision=10),
                                                            puppiWeightNoLep = Var("puppiWeightNoLep()", float, doc="Puppi weight removing leptons",precision=10),
                                                            vtxChi2 = Var("?hasTrackDetails()?vertexChi2():-1", float, doc="vertex chi2",precision=10),
                                                            trkChi2 = Var("?hasTrackDetails()?pseudoTrack().normalizedChi2():-1", float, doc="normalized trk chi2", precision=10),
                                                            dz = Var("?hasTrackDetails()?dz():-1", float, doc="pf dz", precision=10),
                                                            dzErr = Var("?hasTrackDetails()?dzError():-1", float, doc="pf dz err", precision=10),
                                                            d0 = Var("?hasTrackDetails()?dxy():-1", float, doc="pf d0", precision=10),
                                                            d0Err = Var("?hasTrackDetails()?dxyError():-1", float, doc="pf d0 err", precision=10),
                                                            pvAssocQuality = Var("pvAssociationQuality()", int, doc="primary vertex association quality"),
                                                            lostInnerHits = Var("lostInnerHits()", int, doc="lost inner hits"),
                                                            trkQuality = Var("?hasTrackDetails()?pseudoTrack().qualityMask():0", int, doc="track quality mask"),
                                                         )
                                    )
    process.customAK8ConstituentsTable = cms.EDProducer("PatJetConstituentTableProducer",
                                                        candidates = candInput,
                                                        jets = cms.InputTag("finalJetsAK8"),
                                                        jet_radius = cms.double(0.8),
                                                        name = cms.string("FatJetPFCands"),
                                                        name_in_evt = cms.string("evtFatJetPFCands"),
                                                        idx_name = cms.string("pFCandsIdx"),
                                                        sv_idx_name = cms.string("svCandsIdx"),
                                                        nameSV = cms.string("FatJetSVs"),
                                                        idx_nameSV = cms.string("sVIdx"),
                                                        )
    process.customAK4ConstituentsTable = cms.EDProducer("PatJetConstituentTableProducer",
                                                        candidates = candInput,
                                                        jets = cms.InputTag("finalJets"),
                                                        jet_radius = cms.double(0.4),
                                                        name = cms.string("JetPFCands"),
                                                        name_in_evt = cms.string("evtJetPFCands"),
                                                        idx_name = cms.string("pFCandsIdx"),
                                                        sv_idx_name = cms.string("svCandsIdx"),
                                                        nameSV = cms.string("JetSVs"),
                                                        idx_nameSV = cms.string("sVIdx"),
                                                        )
    if not allPF:
        process.customizedPFCandsTask.add(process.finalJetsConstituents)
    process.customizedPFCandsTask.add(process.customConstituentsExtTable)
    process.customizedPFCandsTask.add(process.customAK8ConstituentsTable)
    process.customizedPFCandsTask.add(process.customAK4ConstituentsTable)
    
    if runOnMC:

        process.genJetsAK8Constituents = cms.EDProducer("GenJetPackedConstituentPtrSelector",
                                                    src = cms.InputTag("slimmedGenJetsAK8"),
                                                    cut = cms.string("pt > 100.")
                                                    )

      
        process.genJetsAK4Constituents = process.genJetsAK8Constituents.clone(
                                                    src = cms.InputTag("slimmedGenJets"),
                                                    cut = cms.string("pt > 20")
                                                    )

        if allPF:
            genCandInput = cms.InputTag("packedGenParticles") 
        elif onlyAK4:
            genCandList = cms.VInputTag(cms.InputTag("genJetsAK4Constituents", "constituents"))
            genCandInput =  cms.InputTag("genJetsConstituents")
            process.genJetsConstituents = cms.EDProducer("PackedGenParticlePtrMerger", src = genCandList, skipNulls = cms.bool(True), warnOnSkip = cms.bool(True))
        elif onlyAK8:
            genCandList = cms.VInputTag(cms.InputTag("genJetsAK8Constituents", "constituents"))
            genCandInput =  cms.InputTag("genJetsConstituents")
            process.genJetsConstituents = cms.EDProducer("PackedGenParticlePtrMerger", src = genCandList, skipNulls = cms.bool(True), warnOnSkip = cms.bool(True))
        else:
            genCandList = cms.VInputTag(cms.InputTag("genJetsAK4Constituents", "constituents"), cms.InputTag("genJetsAK8Constituents", "constituents"))
            genCandInput =  cms.InputTag("genJetsConstituents")
            process.genJetsConstituents = cms.EDProducer("PackedGenParticlePtrMerger", src = genCandList, skipNulls = cms.bool(True), warnOnSkip = cms.bool(True))
        process.genJetsParticleTable = cms.EDProducer("SimpleCandidateFlatTableProducer",
                                                         src = genCandInput,
                                                         cut = cms.string(""), #we should not filter after pruning
                                                         name = cms.string("GenCands"),
                                                         doc = cms.string("interesting gen particles from AK4 and AK8 jets"),
                                                         singleton = cms.bool(False), # the number of entries is variable
                                                         extension = cms.bool(False), # this is the main table for the AK8 constituents
                                                         variables = cms.PSet(
                                                                               CandVars,
#                                                                               genPartIdxMother = Var("?numberOfMothers>0?motherRef().key():-1", int, doc="index of the mother particle"),
#                                                                               status  = Var("status", int, doc="Particle status. 1=stable"),
#                                                                               statusFlags = (Var(
#            "statusFlags().isLastCopyBeforeFSR()                  * 16384 +"
#            "statusFlags().isLastCopy()                           * 8192  +"
#            "statusFlags().isFirstCopy()                          * 4096  +"
#            "statusFlags().fromHardProcessBeforeFSR()             * 2048  +"
#            "statusFlags().isDirectHardProcessTauDecayProduct()   * 1024  +"
#            "statusFlags().isHardProcessTauDecayProduct()         * 512   +"
#            "statusFlags().fromHardProcess()                      * 256   +"
#            "statusFlags().isHardProcess()                        * 128   +"
#            "statusFlags().isDirectHadronDecayProduct()           * 64    +"
#            "statusFlags().isDirectPromptTauDecayProduct()        * 32    +"
#            "statusFlags().isDirectTauDecayProduct()              * 16    +"
#            "statusFlags().isPromptTauDecayProduct()              * 8     +"
#            "statusFlags().isTauDecayProduct()                    * 4     +"
#            "statusFlags().isDecayedLeptonHadron()                * 2     +"
#            "statusFlags().isPrompt()                             * 1      ",
#            int, doc=("gen status flags stored bitwise, bits are: "
#                "0 : isPrompt, "
#                "1 : isDecayedLeptonHadron, "
#                "2 : isTauDecayProduct, "
#                "3 : isPromptTauDecayProduct, "
#                "4 : isDirectTauDecayProduct, "
#                "5 : isDirectPromptTauDecayProduct, "
#                "6 : isDirectHadronDecayProduct, "
#                "7 : isHardProcess, "
#                "8 : fromHardProcess, "
#                "9 : isHardProcessTauDecayProduct, "
#                "10 : isDirectHardProcessTauDecayProduct, "
#                "11 : fromHardProcessBeforeFSR, "
#                "12 : isFirstCopy, "
#                "13 : isLastCopy, "
#                "14 : isLastCopyBeforeFSR, ")
#            )),



                                                                          )
                                                     )
        process.genAK8ConstituentsTable = cms.EDProducer("GenJetConstituentTableProducer",
                                                         candidates = genCandInput,
                                                         jets = cms.InputTag("genJetsAK8Constituents"), # Note: The name has "Constituents" in it, but these are the jets
                                                         name = cms.string("GenFatJetCands"),
                                                         name_in_evt = cms.string("evtGenFatJetCands"),
                                                         nameSV = cms.string("GenFatJetSVs"),
                                                         idx_name = cms.string("pFCandsIdx"),
                                                         sv_idx_name = cms.string("svCandsIdx"),
                                                         idx_nameSV = cms.string("sVIdx"),
                                                         readBtag = cms.bool(False))
        process.genAK4ConstituentsTable = cms.EDProducer("GenJetConstituentTableProducer",
                                                         candidates = genCandInput,
                                                         jets = cms.InputTag("genJetsAK4Constituents"), # Note: The name has "Constituents" in it, but these are the jets
                                                         name = cms.string("GenJetCands"),
                                                         name_in_evt = cms.string("evtGenJetCands"),
                                                         nameSV = cms.string("GenJetSVs"),
                                                         idx_name = cms.string("pFCandsIdx"),
                                                         sv_idx_name = cms.string("svCandsIdx"),
                                                         idx_nameSV = cms.string("sVIdx"),
                                                         readBtag = cms.bool(False))
        process.customizedPFCandsTask.add(process.genJetsAK8Constituents)
        process.customizedPFCandsTask.add(process.genJetsAK4Constituents) #Note: For gen need to add jets to the process to keep pt cuts.
        if not allPF:
            process.customizedPFCandsTask.add(process.genJetsConstituents)
        process.customizedPFCandsTask.add(process.genJetsParticleTable)
        process.customizedPFCandsTask.add(process.genAK8ConstituentsTable)
        process.customizedPFCandsTask.add(process.genAK4ConstituentsTable)
    
        process.AK4ConstituentsMatchTable = cms.EDProducer("GenConstituentMapProducer",
                                                         candidates = cms.InputTag("customAK4ConstituentsTable:evtJetPFCands"),
                                                         gencandidates = cms.InputTag("genAK4ConstituentsTable:evtGenJetCands"),
                                                         )
        
        process.customizedPFCandsTask.add(process.AK4ConstituentsMatchTable)
    return process
