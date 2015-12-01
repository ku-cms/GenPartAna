import FWCore.ParameterSet.Config as cms

from Analysis.VLQAna.JetSelector_cfi import *
from FWCore.ParameterSet.VarParsing import VarParsing

options = VarParsing ()

options.register ('dataset',
                                  'none',
                                  VarParsing.multiplicity.singleton,
                                  VarParsing.varType.string,
                                  "Dataset to process")
options.parseArguments()

process = cms.Process("Analyze")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.MessageLogger.cerr.FwkReport.reportEvery = 1000

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.load("UserCode.TprimeAna."+options.dataset+"_cfi")
process.load("Analysis.VLQAna.HbbCandidateProducer_cfi")

process.genana = cms.EDAnalyzer('GenPartAna'
)

process.analyze = cms.EDAnalyzer('TprimeAna'
)


process.TFileService = cms.Service("TFileService", fileName = cms.string(options.dataset+'_gen.root') )

#if options.dataset[0:5] == "Tprime":
process.p = cms.Path(process.genana)
#else:
#	process.p = cms.Path(process.analyze)

#process.p = cms.Schedule(p2)

