import FWCore.ParameterSet.Config as cms
from PhysicsTools.PFNano.addPFCands_cff import addPFCands
from PhysicsTools.PFNano.addBTV import add_BTV
from PhysicsTools.NanoAOD.common_cff import Var
from PhysicsTools.IFNFlavour.ifnFlavour_cff import addIFNFlavour
from PhysicsTools.IFNFlavour.ifnFlavourValidation_cff import addIFNFlavourValidation


# Opt-in switch for the Interleaved Flavour Neutralisation (IFN) parton-flavour
# branches. Enable by ADDING this customiser on the cmsDriver command line, e.g.
#   --customise PhysicsTools/PFNano/pfnano_cff.PFnano_addIFNFlavour
# (combine with a PFnano_customizeMC* customiser). It adds, in parallel with the
# existing ghost-based flavour, the columns:
#   Jet_partonFlavourIFN, GenJet_partonFlavourIFN, GenJetAK8_partonFlavourIFN
# leaving the stock Jet_partonFlavour / GenJet_partonFlavour untouched for
# per-jet comparison. MC only.
def PFnano_addIFNFlavour(process):
    addIFNFlavour(process, addReco=True, addGen=True, addGenAK8=True)
    return process

# Side-by-side validation of IFN vs ghost parton/hadron flavour. Writes
# ifn_validation.root via TFileService alongside the NanoAOD output. Enable on
# the cmsDriver line by ADDING:
#   --customise PhysicsTools/PFNano/pfnano_cff.PFnano_addIFNFlavourValidation
# (combine with a PFnano_customizeMC* customiser -- the latter already adds the
# IFN flavour producers this analyzer reads).
def PFnano_addIFNFlavourValidation(process):
    addIFNFlavourValidation(process)
    return process

# keepInputs can take DeepCSV, DeepJet and DDX (any combination, or use empty placeholder list if no inputs are required)
def PFnano_customizeMC(process):
    addPFCands(process, True)
    add_BTV(process, True, keepInputs=['DeepCSV','DDX'])
    addIFNFlavour(process, addReco=True, addGen=True, addGenAK8=True, useHadrons=True)
#    addIFNFlavourValidation(process)
    process.NANOAODSIMoutput.fakeNameForCrab = cms.untracked.bool(True)  # needed for crab publication
    return process

def PFnano_customizeMC_add_DeepJet(process):
    addPFCands(process, True)
    add_BTV(process, True, keepInputs=['DeepCSV','DeepJet','DDX'])
    process.NANOAODSIMoutput.fakeNameForCrab = cms.untracked.bool(True)  # needed for crab publication
    return process

def PFnano_customizeMC_add_DeepJet_and_Truth(process):
    addPFCands(process, True)
    add_BTV(process, True, keepInputs=['DeepCSV','DeepJet','DDX'], storeAK4Truth="yes")
    process.NANOAODSIMoutput.fakeNameForCrab = cms.untracked.bool(True)  # needed for crab publication
    return process

def PFnano_customizeMC_allPF(process):
    addPFCands(process, True, True)
    add_BTV(process, True, keepInputs=['DeepCSV','DDX'])
    process.NANOAODSIMoutput.fakeNameForCrab = cms.untracked.bool(True)  # needed for crab publication
    return process

def PFnano_customizeMC_allPF_add_DeepJet(process):
    addPFCands(process, True, True)
    add_BTV(process, True, keepInputs=['DeepCSV','DeepJet','DDX'])
    process.NANOAODSIMoutput.fakeNameForCrab = cms.untracked.bool(True)  # needed for crab publication
    return process

def PFnano_customizeMC_AK4JetsOnly(process):
    addPFCands(process, True, False, True)
    add_BTV(process, True, True, keepInputs=['DeepCSV'])
    process.NANOAODSIMoutput.fakeNameForCrab = cms.untracked.bool(True)  # needed for crab publication
    return process

def PFnano_customizeMC_AK4JetsOnly_add_DeepJet(process):
    addPFCands(process, True, False, True)
    add_BTV(process, True, True, keepInputs=['DeepCSV','DeepJet'])
    process.NANOAODSIMoutput.fakeNameForCrab = cms.untracked.bool(True)  # needed for crab publication
    return process

def PFnano_customizeMC_AK8JetsOnly(process):
    addPFCands(process, True, False, False, True)
    add_BTV(process, True, False, True, keepInputs=['DDX'])
    process.NANOAODSIMoutput.fakeNameForCrab = cms.untracked.bool(True)  # needed for crab publication
    return process

def PFnano_customizeMC_noInputs(process):
    add_BTV(process, True, keepInputs=[])
    process.NANOAODSIMoutput.fakeNameForCrab = cms.untracked.bool(True)  # needed for crab publication
    return process


#### DATA customization
def PFnano_customizeData(process):
    addPFCands(process, False)
    add_BTV(process, False, keepInputs=['DeepCSV','DDX'])
    process.NANOAODSIMoutput.fakeNameForCrab = cms.untracked.bool(True)  # needed for crab publication
    return process

def PFnano_customizeData_add_DeepJet(process):
    addPFCands(process, False)
    add_BTV(process, False, keepInputs=['DeepCSV','DeepJet','DDX'])
    process.NANOAODSIMoutput.fakeNameForCrab = cms.untracked.bool(True)  # needed for crab publication
    return process

def PFnano_customizeData_allPF(process):
    addPFCands(process, False, True)
    add_BTV(process, False, keepInputs=['DeepCSV','DDX'])
    process.NANOAODSIMoutput.fakeNameForCrab = cms.untracked.bool(True)  # needed for crab publication
    return process

def PFnano_customizeData_allPF_add_DeepJet(process):
    addPFCands(process, False, True)
    add_BTV(process, False, keepInputs=['DeepCSV','DeepJet','DDX'])
    process.NANOAODSIMoutput.fakeNameForCrab = cms.untracked.bool(True)  # needed for crab publication
    return process

def PFnano_customizeData_AK4JetsOnly(process):
    addPFCands(process, False, False, True)
    add_BTV(process, False, True, keepInputs=['DeepCSV'])
    process.NANOAODSIMoutput.fakeNameForCrab = cms.untracked.bool(True)  # needed for crab publication
    return process

def PFnano_customizeData_AK4JetsOnly_add_DeepJet(process):
    addPFCands(process, False, False, True)
    add_BTV(process, False, True, keepInputs=['DeepCSV','DeepJet'])
    process.NANOAODSIMoutput.fakeNameForCrab = cms.untracked.bool(True)  # needed for crab publication
    return process

def PFnano_customizeData_AK8JetsOnly(process):
    addPFCands(process, False, False, False, True)
    add_BTV(process, False, False, True, keepInputs=['DDX'])
    process.NANOAODSIMoutput.fakeNameForCrab = cms.untracked.bool(True)  # needed for crab publication
    return process

def PFnano_customizeData_noInputs(process):
    add_BTV(process, False, keepInputs=[])
    process.NANOAODSIMoutput.fakeNameForCrab = cms.untracked.bool(True)  # needed for crab publication
    return process
