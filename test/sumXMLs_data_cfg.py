import FWCore.ParameterSet.Config as cms

process = cms.Process("sumXMLs")
process.load("FWCore.MessageLogger.MessageLogger_cfi")

process.source = cms.Source("EmptySource")

process.maxEvents = cms.untracked.PSet(
     input = cms.untracked.int32(1)
)
process.sumCalib = cms.EDAnalyzer("SumHistoCalibration",
	
	 xmlfiles2d = cms.vstring(
"RecoBTag/ImpactParameterLearning/test/Histo_Data80X_2016B_v1_274954-275125.xml"
	 ),

	 xmlfiles3d = cms.vstring(
"RecoBTag/ImpactParameterLearning/test/Histo_Data80X_2016B_v1_274954-275125.xml"
	 ),

	 sum2D = cms.bool(True),
	 sum3D = cms.bool(True),
         writeToDB       = cms.bool(True),
         writeToRootXML  = cms.bool(False),
         writeToBinary   = cms.bool(False)
)

process.PoolDBOutputService = cms.Service("PoolDBOutputService",
    authenticationMethod = cms.untracked.uint32(1),
    loadBlobStreamer = cms.untracked.bool(True),
    catalog = cms.untracked.string('file:mycatalog_new.xml'),
    DBParameters = cms.PSet(
        messageLevel = cms.untracked.int32(0),
        authenticationPath = cms.untracked.string('.')
    ),
    timetype = cms.untracked.string('runnumber'),
    connect = cms.string('sqlite_file:JPcalib_Data80X_2016B_v1_274954-275125.db'),
    toPut = cms.VPSet(cms.PSet(
        record = cms.string('BTagTrackProbability2DRcd'), 
        tag = cms.string('probBTagPDF2D_tag_mc')
    ), 
        cms.PSet(
            record = cms.string('BTagTrackProbability3DRcd'),
            tag = cms.string('probBTagPDF3D_tag_mc')
        ))
)

process.p = cms.Path(process.sumCalib)
