
##*****************************************************
##*****************************************************
##******** for DATA
##*****************************************************
##*****************************************************
import FWCore.ParameterSet.Config as cms

from SimTracker.TrackAssociation.TrackAssociatorByChi2_cfi import *
from SimTracker.TrackAssociation.TrackAssociatorByHits_cfi import *
from RecoBTag.ImpactParameter.impactParameter_cfi import *

process = cms.Process("Calib")

## b-tag infos
bTagInfos = ['impactParameterTagInfos'
]
## b-tag discriminators
bTagDiscriminators = ['jetBProbabilityBJetTags','jetProbabilityBJetTags','trackCountingHighPurBJetTags','trackCountingHighEffBJetTags',
]


process.source = cms.Source(
  "PoolSource",
  # replace 'myfile.root' with the source file you want to use
  fileNames = cms.untracked.vstring(   
  # data from /Jet/Run2012A-PromptReco-v1/AOD
  #        'file:FE2E21F0-0583-E111-8DA9-001D09F241F0.root'
  # data from /BTag/Run2012A-PromptReco-v1/AOD
#   'rfio:/dpm/in2p3.fr/home/cms/phedex/store/mc/Summer12_DR53X/QCD_Pt-120to170_TuneZ2star_8TeV_pythia6/GEN-SIM-RECO/PU_S10_START53_V7A-v1/00000/F88F3DF3-53E7-E111-91BE-003048D4DCD8.root',
# /QCD_Pt-80to120_TuneZ2star_8TeV_pythia6/Summer12_DR53X-PU_S10_START53_V7A-v1/GEN-SIM-RECO
'/store/mc/Summer12_DR53X/QCD_Pt-470to600_TuneZ2star_8TeV_pythia6/AODSIM/PU_S10_START53_V7A-v2/00000/FADB0913-1708-E211-BBB1-00261894383C.root'
 )
  )

process.maxEvents = cms.untracked.PSet(
  input = cms.untracked.int32(10)
  )


process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
#process.GlobalTag.globaltag = "GR_R_52_V7::All"
process.GlobalTag.globaltag = "START53_V7A::All"
#process.GlobalTag.globaltag = "START53_V7A::All"
#process.GlobalTag.globaltag = "GR_R_53_V2B::All"


###-------------------- Import the JEC services -----------------------
#process.load('JetMETCorrections.Configuration.DefaultJEC_cff')
##$$
###-------------------- Import the Jet RECO modules -----------------------
#process.load('RecoJets.Configuration.RecoPFJets_cff')
###-------------------- Turn-on the FastJet density calculation -----------------------
#process.kt6PFJets.doRhoFastjet = True
###-------------------- Turn-on the FastJet jet area calculation for your favorite algorithm -----------------------
#process.ak5PFJets.doAreaFastjet = True
#$$

process.load("Configuration.StandardSequences.Reconstruction_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.StandardSequences.Geometry_cff")

process.load("SimTracker.TrackAssociation.TrackAssociatorByChi2_cfi")
process.load("SimTracker.TrackAssociation.TrackAssociatorByHits_cfi")

process.load("RecoBTau.JetTagComputer.jetTagRecord_cfi")

process.load("SimTracker.TrackHistory.TrackHistory_cff")

# for Impact Parameter based taggers
process.load("RecoBTag.ImpactParameter.negativeOnlyJetBProbabilityComputer_cfi")
process.load("RecoBTag.ImpactParameter.negativeOnlyJetProbabilityComputer_cfi")
process.load("RecoBTag.ImpactParameter.positiveOnlyJetProbabilityComputer_cfi")
process.load("RecoBTag.ImpactParameter.positiveOnlyJetBProbabilityComputer_cfi")
process.load("RecoBTag.ImpactParameter.negativeTrackCounting3D2ndComputer_cfi")
process.load("RecoBTag.ImpactParameter.negativeTrackCounting3D3rdComputer_cfi")
process.load("RecoBTag.Configuration.RecoBTag_cff")
process.load("RecoJets.JetAssociationProducers.ak5JTA_cff")
process.load("RecoBTag.ImpactParameter.negativeOnlyJetBProbabilityJetTags_cfi")
process.load("RecoBTag.ImpactParameter.negativeOnlyJetProbabilityJetTags_cfi")
process.load("RecoBTag.ImpactParameter.positiveOnlyJetProbabilityJetTags_cfi")
process.load("RecoBTag.ImpactParameter.positiveOnlyJetBProbabilityJetTags_cfi")
process.load("RecoBTag.ImpactParameter.negativeTrackCountingHighPur_cfi")
process.load("RecoBTag.ImpactParameter.negativeTrackCountingHighEffJetTags_cfi")
process.load("RecoBTag.ImpactParameter.jetProbabilityBJetTags_cfi")
process.load("RecoBTag.ImpactParameter.jetBProbabilityBJetTags_cfi")

process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")

# for Secondary Vertex taggers
process.load("RecoBTag.SecondaryVertex.secondaryVertexTagInfos_cfi")
process.load("RecoBTag.SecondaryVertex.secondaryVertexNegativeTagInfos_cfi")
process.load("RecoBTag.SecondaryVertex.simpleSecondaryVertexHighEffBJetTags_cfi")
process.load("RecoBTag.SecondaryVertex.simpleSecondaryVertexHighPurBJetTags_cfi")
process.load("RecoBTag.SecondaryVertex.simpleSecondaryVertexNegativeHighEffBJetTags_cfi")
process.load("RecoBTag.SecondaryVertex.simpleSecondaryVertexNegativeHighPurBJetTags_cfi")
process.load("RecoBTag.SecondaryVertex.combinedSecondaryVertexNegativeBJetTags_cfi")
process.load("RecoBTag.SecondaryVertex.combinedSecondaryVertexNegativeES_cfi")
process.load("RecoBTag.SecondaryVertex.combinedSecondaryVertexPositiveBJetTags_cfi")
process.load("RecoBTag.SecondaryVertex.combinedSecondaryVertexPositiveES_cfi")

# for Soft Muon tagger
process.load("RecoBTag.SoftLepton.negativeSoftMuonES_cfi")
process.load("RecoBTag.SoftLepton.positiveSoftMuonES_cfi")
#process.load("RecoBTag.SoftLepton.negativeSoftMuonBJetTags_cfi")
#process.load("RecoBTag.SoftLepton.positiveSoftMuonBJetTags_cfi")

process.load("RecoBTag.SoftLepton.negativeSoftLeptonByPtES_cfi")
process.load("RecoBTag.SoftLepton.positiveSoftLeptonByPtES_cfi")
#process.load("RecoBTag.SoftLepton.negativeSoftMuonByPtBJetTags_cfi")
#process.load("RecoBTag.SoftLepton.positiveSoftMuonByPtBJetTags_cfi")

process.load("PhysicsTools.JetMCAlgos.CaloJetsMCFlavour_cfi")  

process.load("SimTracker.TrackHistory.TrackClassifier_cff")
process.load("RecoBTag.PerformanceMeasurements.MistagAnalyzer_cff")


process.out = cms.OutputModule("PoolOutputModule",
                               outputCommands = cms.untracked.vstring('keep *', 
                                                                      'keep recoJetTags_*_*_*'),
                               fileName = cms.untracked.string('testtt.root')
                               )


#-------------------------------------
#Produce PFJets from PF2PAT

# PAT Layer 0+1
process.load("PhysicsTools.PatAlgos.patSequences_cff")

#process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(False))
#from PhysicsTools.PatAlgos.tools.pfTools import *
#
##postfix = "PF2PAT"  # to have only PF2PAT
#postfix = "PFlow"
#
#usePF2PAT(process,runPF2PAT=True,
#          jetAlgo='AK5', runOnMC=True, postfix=postfix,
#          jetCorrections=('AK5PFchs', ['L1FastJet', 'L2Relative', 'L3Absolute'])#,
#          #pvCollection=cms.InputTag('goodOfflinePrimaryVertices')
#          )
#process.pfPileUpPF2PAT.checkClosestZVertex = False
##applyPostfix(process,"pfPileUp",postfix).checkClosestZVertex = cms.bool(False) 
#
#process.pfPileUpPF2PAT.Enable = True
#process.pfPileUpPF2PAT.Vertices = cms.InputTag('goodOfflinePrimaryVertices')

## Standard PAT Configuration File
process.load("PhysicsTools.PatAlgos.patSequences_cff")
postfix = "PFlow"

## Jet energy corrections
jetCorrectionsAK5 = ('AK5PFchs', ['L1FastJet', 'L2Relative', 'L3Absolute'])

from PhysicsTools.PatAlgos.tools.pfTools import *
usePF2PAT(process,runPF2PAT=True, jetAlgo='AK5', runOnMC=True, postfix=postfix,
          jetCorrections=jetCorrectionsAK5, pvCollection=cms.InputTag('goodOfflinePrimaryVertices'))

## Top projections in PF2PAT
getattr(process,"pfNoPileUp"+postfix).enable = True

from PhysicsTools.PatAlgos.tools.jetTools import *
## Switch the default jet collection (done in order to use the above specified b-tag infos and discriminators)
switchJetCollection(process,
    cms.InputTag('pfNoTau'+postfix),
    doJTA        = True,
    doBTagging   = True,
    btagInfo     = bTagInfos,
    btagdiscriminators = bTagDiscriminators,
    jetCorrLabel = jetCorrectionsAK5,
    doType1MET   = False,
    genJetCollection = cms.InputTag('ak5GenJetsNoNu'),
    doJetID      = False,
    postfix      = postfix
)

## Add TagInfos to PAT jets
patJets = ['patJets'+postfix]

for m in patJets:
    if hasattr(process,m):
        print "Switching 'addTagInfos' for " + m + " to 'True'"
        setattr( getattr(process,m), 'addTagInfos', cms.bool(True) )



#-------------------------------------
#Redo the primary vertex

from PhysicsTools.SelectorUtils.pvSelector_cfi import pvSelector

process.goodOfflinePrimaryVertices = cms.EDFilter(
  "PrimaryVertexObjectFilter",
  filterParams = pvSelector.clone( minNdof = cms.double(4.0), maxZ = cms.double(24.0) ),
  src=cms.InputTag('offlinePrimaryVertices')
  )

#-------------------------------------
#JetID
process.load('PhysicsTools.SelectorUtils.pfJetIDSelector_cfi')
process.load('PhysicsTools.SelectorUtils.jetIDSelector_cfi')
process.jetIDSelector.version = cms.string('PURE09')


#-------------------------------------
#Filter for PFJets
process.PATJetsFilter = cms.EDFilter("PATJetSelector",    
                                    src = cms.InputTag('selectedPatJets'+postfix),
                                    cut = cms.string("pt > 20.0 && abs(eta) < 2.5 && neutralHadronEnergyFraction < 0.99 && neutralEmEnergyFraction < 0.99 && nConstituents > 1 && chargedHadronEnergyFraction > 0.0 && chargedMultiplicity > 0.0 && chargedEmEnergyFraction < 0.99"),
                                    #filter = cms.bool(True)
                                    )
#---------------------------------------
process.load("bTag.CommissioningCommonSetup.caloJetIDFilter_cff")


#Filter for removing scraping events
process.noscraping = cms.EDFilter("FilterOutScraping",
                                  applyfilter = cms.untracked.bool(True),
                                  debugOn = cms.untracked.bool(False),
                                  numtrack = cms.untracked.uint32(10),
                                  thresh = cms.untracked.double(0.25)
                                  )


#Filter for good primary vertex
process.load("RecoVertex.PrimaryVertexProducer.OfflinePrimaryVertices_cfi")
process.primaryVertexFilter = cms.EDFilter("GoodVertexFilter",
                                           vertexCollection = cms.InputTag('offlinePrimaryVertices'),
                                           minimumNDOF = cms.uint32(4) ,
                                           maxAbsZ = cms.double(24), 
                                           maxd0 = cms.double(2)	
                                           )


#---------------------------------------
# process.TFileService = cms.Service("TFileService", fileName = cms.string("TrackTree.root") )
#process.TFileService = cms.Service("TFileService", fileName = cms.string("JetTree.root") )



process.load("RecoBTag.ImpactParameterLearning.ImpactParameterCalibration_cfi")
process.ipCalib.Jets   = cms.InputTag('PATJetsFilter')
process.ipCalib.primaryVertexColl     = cms.InputTag('goodOfflinePrimaryVertices')

#process.ipCalib.jetCorrector   = cms.string('ak5PFL1FastL2L3')

  
    
#---------------------------------------
# trigger selection !
import HLTrigger.HLTfilters.triggerResultsFilter_cfi as hlt
process.JetHLTFilter = hlt.triggerResultsFilter.clone(
     triggerConditions = cms.vstring(
        "HLT_PFJet80_v5"
     ),
     hltResults = cms.InputTag("TriggerResults","","HLT"),
     l1tResults = cms.InputTag( "" ),
     throw = cms.bool( False) #set to false to deal with missing triggers while running over different trigger menus
     )
#---------------------------------------


process.p = cms.Path(
  #$$
  #process.JetHLTFilter*
  process.noscraping
  *process.primaryVertexFilter
  *process.goodOfflinePrimaryVertices
  *getattr(process,"patPF2PATSequence"+postfix)
  *process.PATJetsFilter
#  *process.myPartons
#  *process.AK5Flavour
#  *process.noscraping
#  *process.primaryVertexFilter
#  *process.ak5JetTracksAssociatorAtVertex
#  *process.btagging
#  *process.positiveOnlyJetProbabilityJetTags*process.negativeOnlyJetProbabilityJetTags
#  *process.negativeTrackCountingHighEffJetTags*process.negativeTrackCountingHighPur
#  *process.secondaryVertexTagInfos*process.simpleSecondaryVertexHighEffBJetTags*process.simpleSecondaryVertexHighPurBJetTags
#  *process.secondaryVertexNegativeTagInfos*process.simpleSecondaryVertexNegativeHighEffBJetTags*process.simpleSecondaryVertexNegativeHighPurBJetTags
#  *process.combinedSecondaryVertexNegativeBJetTags*process.combinedSecondaryVertexPositiveBJetTags
  #*process.negativeSoftMuonBJetTags*process.positiveSoftMuonBJetTags	
  #*process.negativeSoftLeptonByPtBJetTags*process.positiveSoftLeptonByPtBJetTags
  #*process.positiveSoftLeptonByPtBJetTags	
  #*process.negativeOnlyJetBProbabilityJetTags*process.positiveOnlyJetBProbabilityJetTags
  *process.ipCalib
  
  )

from PhysicsTools.PatAlgos.patEventContent_cff import patEventContent
process.out = cms.OutputModule("PoolOutputModule",
                               fileName = cms.untracked.string('patTuple.root'),
                               # save only events passing the full path
                               SelectEvents   = cms.untracked.PSet( SelectEvents = cms.vstring('p') ),
                               # save PAT Layer 1 output; you need a '*' to
                               # unpack the list of commands 'patEventContent'
                               outputCommands = cms.untracked.vstring('drop *',
                                                'keep *_*_*_Mistag'  )
                               )

#process.outpath = cms.EndPath(process.out)
## Output Module Configuration (expects a path 'p')




