// -*- C++ -*-
//
// Package:    UserCode/GenPartAna
// Class:      GenPartAna
// 
/**\class GenPartAna GenPartAna.cc UserCode/GenPartAna/plugins/GenPartAna.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Robert Stringer
//         Created:  Fri, 16 Oct 2015 18:59:00 GMT
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/Math/interface/deltaPhi.h"
#include <iostream>
#include <TH2F.h>
#include <math.h>

//
//

// class declaration
//

class GenPartAna : public edm::EDAnalyzer {
   public:
      explicit GenPartAna(const edm::ParameterSet&);
      ~GenPartAna();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
      virtual void beginJob() override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;

      //virtual void beginRun(edm::Run const&, edm::EventSetup const&) override;
      //virtual void endRun(edm::Run const&, edm::EventSetup const&) override;
      //virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
      //virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;

      // ----------member data ---------------------------

      edm::Service<TFileService> fs;
      std::map< std::string, TH1F* > histos;
      std::map< std::string, TH2F* > histos2;


};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
GenPartAna::GenPartAna(const edm::ParameterSet& iConfig)

{
   //now do what ever initialization is needed

}


GenPartAna::~GenPartAna()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
GenPartAna::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;
   using namespace std;


#ifdef THIS_IS_AN_EVENT_EXAMPLE
   Handle<ExampleData> pIn;
   iEvent.getByLabel("example",pIn);
#endif
   
#ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
   ESHandle<SetupData> pSetup;
   iSetup.get<SetupRecord>().get(pSetup);
#endif
      edm::Handle<std::vector<float> > mh_genPartPt;
      edm::Handle<std::vector<float> > mh_genPartEta;
      edm::Handle<std::vector<float> > mh_genPartPhi;
      edm::Handle<std::vector<float> > mh_genPartE;
      edm::Handle<std::vector<float> > mh_genPartID;
      edm::Handle<std::vector<float> > mh_genPartStatus;
      edm::Handle<std::vector<float> > mh_genPartMomID;
      edm::Handle<std::vector<float> > mh_genPartMom0ID;
      edm::Handle<std::vector<float> > mh_genPartMom1ID;
      edm::Handle<std::vector<float> > mh_genPartDau0ID;
      edm::Handle<std::vector<float> > mh_genPartDau1ID;

      edm::Handle<std::vector<unsigned> > mh_bjetIdxs;
      edm::Handle<std::vector<unsigned> > mh_tjetIdxs;
      edm::Handle<std::vector<unsigned> > mh_hjetIdxs;

      edm::Handle<std::vector<float> > mh_jetAK4Pt;
      edm::Handle<std::vector<float> > mh_jetAK4GenPartonPt;
      edm::Handle<std::vector<float> > mh_jetAK4Eta;
      edm::Handle<std::vector<float> > mh_jetAK4GenPartonY;
      edm::Handle<std::vector<float> > mh_jetAK4Phi;
      edm::Handle<std::vector<float> > mh_jetAK4E;
      edm::Handle<std::vector<float> > mh_jetAK4HadronFlavour;
      edm::Handle<std::vector<float> > mh_jetAK4CSV;
      edm::Handle<std::vector<float> > mh_jetAK4Mass;
      edm::Handle<std::vector<unsigned> > mh_jetAK4Good;

      edm::Handle<std::vector<float> > mh_jetAK8Pt;
      edm::Handle<std::vector<float> > mh_jetAK8Eta;
      edm::Handle<std::vector<float> > mh_jetAK8Phi;
      edm::Handle<std::vector<float> > mh_jetAK8E;
      edm::Handle<std::vector<float> > mh_jetAK8Mass;
      edm::Handle<std::vector<float> > mh_jetAK8subjetIndex0;
      edm::Handle<std::vector<float> > mh_jetAK8subjetIndex1;
      edm::Handle<std::vector<float> > mh_jetAK8HadronFlavour;
      edm::Handle<std::vector<float> > mh_jetAK8tau1;
      edm::Handle<std::vector<float> > mh_jetAK8tau2;
      edm::Handle<std::vector<float> > mh_jetAK8tau3;
      edm::Handle<std::vector<unsigned> > mh_jetAK8Good;
      edm::Handle<std::vector<float> > mh_jetAK8nSubjets;

      edm::Handle<double> mh_htak4jets;
      edm::Handle<double> mh_htak4trigjets;

      edm::Handle<std::vector<unsigned> > mh_tjetsIdx_lsubj;
      edm::Handle<std::vector<unsigned> > mh_hjetsIdx_lsubj;

//      iEvent.getByLabel("jetsAK4", "jetAK4Pt", mh_jetAK4Pt);
//      iEvent.getByLabel("jetsAK4", "jetAK4GenPartonPt", mh_jetAK4GenPartonPt);
      iEvent.getByLabel("jetsAK4", "jetAK4Eta", mh_jetAK4Eta);
//      iEvent.getByLabel("jetsAK4", "jetAK4GenPartonY", mh_jetAK4GenPartonY);
      iEvent.getByLabel("jetsAK4", "jetAK4Phi", mh_jetAK4Phi);
      iEvent.getByLabel("jetsAK4", "jetAK4E", mh_jetAK4E);
//      iEvent.getByLabel("jetsAK4", "jetAK4HadronFlavour", mh_jetAK4HadronFlavour);
//      iEvent.getByLabel("jetsAK4", "jetAK4CSV", mh_jetAK4CSV);
      iEvent.getByLabel("jetsAK4", "jetAK4Mass", mh_jetAK4Mass);

//      iEvent.getByLabel("jetsAK8", "jetAK8Pt", mh_jetAK8Pt);
      iEvent.getByLabel("jetsAK8", "jetAK8Eta", mh_jetAK8Eta);
      iEvent.getByLabel("jetsAK8", "jetAK8Phi", mh_jetAK8Phi);
//      iEvent.getByLabel("jetsAK8", "jetAK8E", mh_jetAK8E);
      iEvent.getByLabel("jetsAK8", "jetAK8Mass", mh_jetAK8Mass);
//      iEvent.getByLabel("jetsAK8", "jetAK8subjetIndex0", mh_jetAK8subjetIndex0);
//      iEvent.getByLabel("jetsAK8", "jetAK8subjetIndex1", mh_jetAK8subjetIndex1);
//      iEvent.getByLabel("jetsAK8", "jetAK8HadronFlavour", mh_jetAK8HadronFlavour);
//      iEvent.getByLabel("jetsAK8", "jetAK8nSubJets", mh_jetAK8nSubjets );

      iEvent.getByLabel("presel", "htak4jets", mh_htak4jets);
      iEvent.getByLabel("presel", "htak4trigjets", mh_htak4trigjets);

      iEvent.getByLabel("genPart", "genPartPt", mh_genPartPt);
      iEvent.getByLabel("genPart", "genPartEta", mh_genPartEta);
      iEvent.getByLabel("genPart", "genPartPhi", mh_genPartPhi);
//      iEvent.getByLabel("genPart", "genPartE", mh_genPartE);
      iEvent.getByLabel("genPart", "genPartID", mh_genPartID);
      iEvent.getByLabel("genPart", "genPartStatus", mh_genPartStatus);
//      iEvent.getByLabel("genPart", "genPartMomID", mh_genPartMomID);
      iEvent.getByLabel("genPart", "genPartMom0ID", mh_genPartMom0ID);
//      iEvent.getByLabel("genPart", "genPartMom1ID", mh_genPartMom1ID);
//      iEvent.getByLabel("genPart", "genPartDau0ID", mh_genPartDau0ID);
//      iEvent.getByLabel("genPart", "genPartDau1ID", mh_genPartDau1ID);


      iEvent.getByLabel("presel", "bjetIdxs", mh_bjetIdxs);
      iEvent.getByLabel("presel", "tjetIdxs", mh_tjetIdxs);
      iEvent.getByLabel("presel", "hjetIdxs", mh_hjetIdxs);

      iEvent.getByLabel("presel","ak4goodjets", mh_jetAK4Good);
      iEvent.getByLabel("presel","ak8goodjets", mh_jetAK8Good);


      iEvent.getByLabel("anavars", "tjetsIdx-lsubj",mh_tjetsIdx_lsubj);
      iEvent.getByLabel("anavars", "hjetsIdx-lsubj",mh_hjetsIdx_lsubj);

      const std::vector<float> * genPartPt = mh_genPartPt.product();
      const std::vector<float> * genPartID = mh_genPartID.product();
      const std::vector<float> * genPartMomID = mh_genPartMom0ID.product();
      const std::vector<float> * genPartStatus= mh_genPartStatus.product();
      const std::vector<float> * genPartEta = mh_genPartEta.product();
      const std::vector<float> * genPartPhi = mh_genPartPhi.product();

      const std::vector<unsigned> * bjetIdxs = mh_bjetIdxs.product();
      const std::vector<unsigned> * tjetIdxs = mh_tjetIdxs.product();
      const std::vector<unsigned> * hjetIdxs = mh_hjetIdxs.product();

      const std::vector<unsigned> * tjetIdxs_lsubj = mh_tjetIdxs.product();
      const std::vector<unsigned> * hjetIdxs_lsubj = mh_hjetIdxs.product();

//      const std::vector<float> * jetcsv = mh_jetAK4CSV.product();
//      const std::vector<float> * jetpt = mh_jetAK4Pt.product();
     const std::vector<float> * jeteta = mh_jetAK4Eta.product();
      const std::vector<float> * jetphi = mh_jetAK4Phi.product();
//      const std::vector<float> * jetM = mh_jetAK4Mass.product();
//      const std::vector<float> * jetflavor = mh_jetAK4HadronFlavour.product();
      const std::vector<unsigned> * jetgood = mh_jetAK4Good.product();

//      const std::vector<float> * jet8pt = mh_jetAK8Pt.product();
      const std::vector<float> * jet8eta = mh_jetAK8Eta.product();
      const std::vector<float> * jet8phi = mh_jetAK8Phi.product();
      const std::vector<float> * jet8M = mh_jetAK8Mass.product();
      const std::vector<unsigned> * jet8good = mh_jetAK8Good.product();

   //   const double * htak4jets = mh_htak4jets.product();
      const double * htak4trigjets = mh_htak4trigjets.product();


	if (*htak4trigjets < 800) return;

        unsigned bfromtidx = 99;
        unsigned tidx = 99;
        unsigned h1idx = 99;
        unsigned h2idx = 99;
	unsigned hidx = 99;
	unsigned bAssoc = 99;

	cout << "Find t h " << endl;

        for(unsigned int i = 0 ; i < genPartID->size() ; i++)
        {
	   int status = (*genPartStatus)[i];
	   if(status != 23 and status != 22) continue;

                if(fabs((*genPartID)[i]) == 6)
                {
                   tidx = i;
                }
                if(fabs((*genPartID)[i]) == 25)
                {
                   hidx = i;
                } 
        }
	cout << "H: " << hidx << " t: " << tidx << endl;	

        for(unsigned int i = 0 ; i < genPartID->size() ; i++)
        {
	   int status = (*genPartStatus)[i];
	   if(status != 23 and status != 22) continue;

           if(fabs((*genPartID)[i]) == 5)
           {
                if(fabs((*genPartMomID)[i]) == 6)
                {
                   bfromtidx = i;    
                }       
                if((*genPartMomID)[i] == 25 && h1idx == 99)
                {
                   h1idx = i;   
                } else if ((*genPartMomID)[i] == 25)
                  {     
                        h2idx = i;
                  }
		if(fabs((*genPartMomID)[i]) == 1 or fabs((*genPartMomID)[i]) == 2 or fabs((*genPartMomID)[i]) == 21 )
		{
		   bAssoc = i;
		}
           }

        }
	

	histos[ "ptbAssoc" ] ->Fill((*genPartPt)[bAssoc]);
	histos[ "etabAssoc" ] ->Fill((*genPartEta)[bAssoc]);
	histos2[ "ptetabAssoc" ] ->Fill((*genPartPt)[bAssoc],(*genPartEta)[bAssoc]); 

//	if(h1idx == 99)
//	for(unsigned int i = 0 ; i < genPartID->size() ; i++) {
	//		if((*genPartStatus)[i] != 23 and (*genPartStatus)[i] != 22) continue;
		//	std::cout << (*genPartID)[i]  << " " << (*genPartEta)[i] << "  " << (*genPartPhi)[i] << " " << (*genPartStatus)[i]<<" "<< (*genPartMomID)[i] << endl;
//	}
        cout << "h1: "<< h1idx << " h2: " << h2idx << " t->b: " << bfromtidx << endl;  

	unsigned td1idx,td2idx, td3idx;
	td1idx=td2idx=td3idx = 99;


	    for(unsigned int i = 0 ; i < genPartMomID->size() ; i++) 
	    {
		if( fabs((*genPartMomID)[i]) == 6 and ((*genPartStatus)[i] == 23 or (*genPartStatus)[i] == 22) ) 
		{
		    td1idx = i;
		} 
		if( fabs((*genPartMomID)[i]) ==  24 and ((*genPartStatus)[i] == 23 or (*genPartStatus)[i] == 22))
		{

		   td2idx = i;

		}
               if( fabs((*genPartMomID)[i]) ==  24 and ((*genPartStatus)[i] == 23 or (*genPartStatus)[i] == 22) and td2idx != 99)
                {

                   td3idx = i;

                }
	    }

	    if ( h1idx != 99 and h2idx != 99 and td2idx != 99 and td3idx != 99 and fabs((*genPartID)[td2idx]) >= 1 and fabs((*genPartID)[td3idx]) <=5)  // W->qq ; h->bb
	    {
		histos["isHbbWqq"]->Fill(1);
		cout << "isHbbWqq" << endl;
	    }
	    else
		return;
	//	cout << "notHbbWqq" << endl;

	if( h1idx != 99 and h2idx != 99) 
	{ 
           histos["drHiggs"]->Fill(deltaR((*genPartEta)[h1idx],(*genPartPhi)[h1idx],(*genPartEta)[h2idx],(*genPartPhi)[h2idx]));
           histos2["drHiggsPt"]->Fill((*genPartPt)[hidx],deltaR((*genPartEta)[h1idx],(*genPartPhi)[h1idx],(*genPartEta)[h2idx],(*genPartPhi)[h2idx])) ;

	    if( td1idx ==99 or td2idx ==99 or td3idx ==99) {
		cout << "Couldn't find all decay products. "<< td1idx <<" "<< td2idx << " " << td3idx << endl;
	//		cout << "--------" << endl;	
	//	        for(unsigned int i = 0 ; i < genPartID->size() ; i++) {
          //              std::cout << (*genPartID)[i]  << " " << (*genPartEta)[i] << "  " << (*genPartPhi)[i] << " " << (*genPartStatus)[i]<<" "<< (*genPartMomID)[i] << endl;
        //		}
	//		cout << "--------" << endl;	
	    }
	    else {
  
	    	float dr12 = deltaR((*genPartEta)[td1idx],(*genPartPhi)[td1idx],(*genPartEta)[td2idx],(*genPartPhi)[td2idx]) ;
	    	float dr23 = deltaR((*genPartEta)[td2idx],(*genPartPhi)[td2idx],(*genPartEta)[td3idx],(*genPartPhi)[td3idx]) ;
	    	float dr13= deltaR((*genPartEta)[td1idx],(*genPartPhi)[td1idx],(*genPartEta)[td3idx],(*genPartPhi)[td3idx]) ;

	        histos2["drt-bjj"]->Fill((*genPartPt)[tidx], fmax(fmax(dr12,dr23),dr13) );
	    }	
	}
	
	
	
       
        bool bh1,bh2,bt;
        bh1=bh2=bt=false;       
	bool btag1,btag2;
	btag1 = btag2 = false;
 
        for(unsigned int i = 0 ; i < jetgood->size() ; i++)
        {
            if(deltaR((*genPartEta)[bfromtidx],(*genPartPhi)[bfromtidx],(*jeteta)[(*jetgood)[i]],(*jetphi)[(*jetgood)[i]]) < 0.5)
            {
            	bt = true;
            }   
            if(deltaR((*genPartEta)[h1idx],(*genPartPhi)[h1idx],(*jeteta)[(*jetgood)[i]],(*jetphi)[(*jetgood)[i]]) < 0.3)
            {
                bh1 = true;
		btag1 = (bjetIdxs->end() != std::find(bjetIdxs->begin(), bjetIdxs->end(), (*jetgood)[i]));
		cout << "matched idx: " <<(*jetgood)[i] << "  Bjets: ";
		for(unsigned j : *bjetIdxs) cout << j << ",";
		cout <<endl;    
            }           
            if(deltaR((*genPartEta)[h2idx],(*genPartPhi)[h2idx],(*jeteta)[(*jetgood)[i]],(*jetphi)[(*jetgood)[i]]) < 0.3)
            {
               bh2 = true;
	       btag2 = (bjetIdxs->end() != std::find(bjetIdxs->begin(), bjetIdxs->end(), (*jetgood)[i]));
            }           
        }
        if(bt)
        	histos["ak4tMatched"]->Fill(1);
        if(bh1 and bh2)
        	histos["ak4hMatched"]->Fill(1);

	if((bh1 and btag1) and (bh2 and btag2))
		histos["ak4hMatchedTagged"]->Fill(2);
	else
	  if((bh1 and btag1) or (bh2 and btag2))
		histos["ak4hMatchedTagged"]->Fill(1);
	

        //tidx = -1;
        //h1idx = -1;
       // h2idx = -1;

	//cout << "Match t h" << endl;

	bt = false;
	bh1 = false;

        for(unsigned int i = 0 ; i < jet8good->size() ; i++)
        {

            if(deltaR((*genPartEta)[tidx],(*genPartPhi)[tidx],(*jet8eta)[(*jet8good)[i]],(*jet8phi)[(*jet8good)[i]]) < 0.5)
            {
                bt = true;
                histos["ak8tMatchedMass"]->Fill((*jet8M)[(*jet8good)[i]]);
		histos["ak8tMatchedTagged"]->Fill(tjetIdxs->end() != std::find(tjetIdxs->begin(), tjetIdxs->end(), i));
		histos["ak8tMatchedTagged_lsubj"]->Fill(tjetIdxs_lsubj->end() != std::find(tjetIdxs_lsubj->begin(), tjetIdxs_lsubj->end(), i));
		cout << "ttags: ";
		for (auto k: *tjetIdxs) std::cout << k << ' '; 
		cout << " jet: " << i <<endl;
            }
            if(deltaR((*genPartEta)[hidx],(*genPartPhi)[hidx],(*jet8eta)[(*jet8good)[i]],(*jet8phi)[(*jet8good)[i]]) < 0.5)
            {
                bh1 = true;
                histos["ak8hMatchedMass"]->Fill((*jet8M)[(*jet8good)[i]]);
		histos["ak8hMatchedTagged"]->Fill (hjetIdxs->end() != std::find(hjetIdxs->begin(), hjetIdxs->end(), i));
		histos["ak8hMatchedTagged_lsubj"]->Fill (hjetIdxs_lsubj->end() != std::find(hjetIdxs_lsubj->begin(), hjetIdxs_lsubj->end(), i));
            }
        }
         
	if (!bt and !bh1) {
 		for (unsigned int i = 0 ; i < genPartID->size() ; i++) {
	//		if((*genPartStatus)[i] != 23 ) continue;
			std::cout << (*genPartID)[i]  << " " << (*genPartEta)[i] << "  " << (*genPartPhi)[i] << " " << (*genPartStatus)[i]<<" "<< (*genPartMomID)[i] << endl;
		}
		for(unsigned int i = 0 ; i < jet8good->size() ; i++)
		     	cout << "Jet: " << (*jet8eta)[(*jet8good)[i]] << " " << (*jet8phi)[(*jet8good)[i]]	<< endl;
		
	}
        if(bt) {
                histos["ak8tMatched"]->Fill(1);
	}
        if(bh1)
	{ 
                histos["ak8hMatched"]->Fill(1);         
	}


//      histos["ak4hMatched"]->Fill(
//
//      //      histos["ak8tMatched"]->Fill
//      //      histos["ak8hMatched"]


}


// ------------ method called once each job just before starting event loop  ------------
void 
GenPartAna::beginJob()
{

  histos[ "ak4tMatched" ]= fs->make<TH1F>( "ak4tMatched","ak4tMatched",3,0,3);
  histos[ "ak4hMatched" ]= fs->make<TH1F>( "ak4hMatched","ak4hMatched",3,0,3);
  histos[ "ak8tMatched" ]= fs->make<TH1F>( "ak8tMatched","ak8tMatched",3,0,3);
  histos[ "ak8hMatched" ]= fs->make<TH1F>( "ak8hMatched","ak8hMatched",3,0,3);
  histos[ "ak8tMatchedTagged" ]= fs->make<TH1F>( "ak8tMatchedTagged","ak8tMatchedTagged",3,0,3);
  histos[ "ak8hMatchedTagged" ]= fs->make<TH1F>( "ak8hMatchedTagged","ak8hMatchedTagged",3,0,3);
  histos[ "ak8tMatchedTagged_lsubj" ]= fs->make<TH1F>( "ak8tMatchedTagged_lsubj","ak8tMatchedTagged_lsubj",3,0,3);
  histos[ "ak8hMatchedTagged_lsubj" ]= fs->make<TH1F>( "ak8hMatchedTagged_lsubj","ak8hMatchedTagged_lsubj",3,0,3);
  histos[ "ak4hMatchedTagged" ]= fs->make<TH1F>( "ak4hMatchedTagged","ak4hMatchedTagged",3,0,3);

  histos[ "ak8tMatchedMass" ]= fs->make<TH1F>( "ak8tMatchedMass","ak8tMatchedMass",300,0,300);
  histos[ "ak8hMatchedMass" ]= fs->make<TH1F>( "ak8hMatchedMass","ak8hMatchedMass",300,0,300);
  histos[ "drHiggs" ] = fs->make<TH1F> ( "drHiggs", "drHiggs", 50 ,0, 5);
  histos["isHbbWqq"] = fs->make<TH1F>( "isHbbWqq","isHbbWqq",3,0,3);

  histos["ptbAssoc"] = fs->make<TH1F>( "ptbAssoc","p_{T} b assoc",1000,0,1000);
  histos["etabAssoc"] = fs->make<TH1F>( "etabAssoc","eta b assoc",50,-5,5);

  histos2[ "ptetabAssoc" ] = fs->make<TH2F> ( "ptetabAssoc" , "pt vs eta b assoc",1000,0,1000,50,-5,5);
  histos2[ "drHiggsPt" ] = fs->make<TH2F> ( "drHiggsPt", "drHiggsPt", 500,0,1500,50 ,0 ,5 );
  histos2[ "drt-bjj" ] = fs->make<TH2F> ( "drt-bjj", "drt-bjj", 500,0,1500,50 ,0 ,5 );

  

}

// ------------ method called once each job just after ending the event loop  ------------
void 
GenPartAna::endJob() 
{
}

// ------------ method called when starting to processes a run  ------------
/*
void 
GenPartAna::beginRun(edm::Run const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when ending the processing of a run  ------------
/*
void 
GenPartAna::endRun(edm::Run const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when starting to processes a luminosity block  ------------
/*
void 
GenPartAna::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when ending the processing of a luminosity block  ------------
/*
void 
GenPartAna::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
GenPartAna::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(GenPartAna);
