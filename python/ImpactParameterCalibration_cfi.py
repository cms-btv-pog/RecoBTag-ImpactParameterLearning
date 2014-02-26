import FWCore.ParameterSet.Config as cms

ipCalib = cms.EDAnalyzer("ImpactParameterCalibration",
                                 
                                 writeToDB                = cms.bool(False),
                                 writeToBinary            = cms.bool(False),
                                 nBins                    = cms.int32(50),
                                 maxSignificance          = cms.double(50.0),
                                 writeToRootXML           = cms.bool(True),
                                 inputCategories          = cms.string('HardCoded'),
                                 primaryVertexColl        = cms.InputTag("goodOfflinePrimaryVertices"),
                                 Jets                     = cms.InputTag('ak5PFJets'),
                                 MinPt                    = cms.double(10.0),
                                 MaxEta                   = cms.double(2.5)    
                                 
                                 )
