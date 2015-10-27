#include "specialAna.hh"
#include "HistClass.hh"
#include "Tools/Tools.hh"
#include <csignal>
#include <algorithm>
#include "TVector3.h"
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-local-typedefs"
#include "boost/format.hpp"
#pragma GCC diagnostic pop

specialAna::specialAna( const Tools::MConfig &cfg,
			Systematics *syst_shifter,
			EventSelector *selector
			) :
  runOnData(       cfg.GetItem< bool >( "General.RunOnData" ) ),
  debug ( cfg.GetItem< bool > ("ZH.DebugSwitch")  ),
  //   m_TauType(       cfg.GetItem< string >( "Tau.Type.Rec" ) ),
  doSampleWeighting(cfg.GetItem< bool >("General.DoSampleWeighting")),
  m_trigger_string(Tools::splitString< std::string >(cfg.GetItem< std::string >("ZH.trigger_list"))),
  m_eventSelector( selector )
{
  string safeFileName = "SpecialHistos.root";
  events_=0;
  file1 = new TFile(safeFileName.c_str(), "RECREATE");
  file1->cd();
  // number of events, saved in a histogram
  HistClass::CreateHistoUnchangedName("h_counters", 10, 0, 11, "N_{events}");

  sample_weight = 1.;

  HistClass::CreateHisto("nTau",   8, -0.5, 7.5, "nTau  (Reco)");
  HistClass::CreateHisto("nMuon",  8, -0.5, 7.5, "nMuon (Reco)");
  HistClass::CreateHisto("nEle",   8, -0.5, 7.5, "nEle  (Reco)");
  HistClass::CreateHisto("nLep",   8, -0.5, 7.5, "nLep  (Reco)");
  
  HistClass::CreateHisto("Muon_Pt_Lead",     250, 0, 1000, "p_{T}^{Mu}  lead [GeV]");
  HistClass::CreateHisto("Muon_Pt_SubLead",  250, 0, 1000, "p_{T}^{Mu}  sublead [GeV]");
  HistClass::CreateHisto("Ele_Pt_Lead",      250, 0, 1000, "p_{T}^{Ele} lead [GeV]");
  HistClass::CreateHisto("Ele_Pt_SubLead",   250, 0, 1000, "p_{T}^{Ele} sublead [GeV]");

  HistClass::CreateHisto("Zmass_mu",  400, 0, 200, "M_{Z} [GeV]");
  HistClass::CreateHisto("Zmass_ele", 400, 0, 200, "M_{Z} [GeV]");
  HistClass::CreateHisto("Zmass_lep", 400, 0, 200, "M_{Z} [GeV]");
    
  HistClass::CreateHisto("Zmass_mu_selEvt",  400, 0, 200, "M_{Z} [GeV]");
  HistClass::CreateHisto("Zmass_ele_selEvt", 400, 0, 200, "M_{Z} [GeV]");
  HistClass::CreateHisto("Zmass_lep_selEvt", 400, 0, 200, "M_{Z} [GeV]");

  HistClass::CreateHisto("Z_pt_selEvt", 500, 0, 2000, "Dilepton p_{T} [GeV]");
  HistClass::CreateHisto("H_pt_selEvt", 500, 0, 2000, "DiTau p_{T} [GeV]");

  HistClass::CreateHisto("deltaRmumuReco_selEvt",  200, 0, 2.0, "#Delta R^{Mu}");
  HistClass::CreateHisto("deltaReeReco_selEvt",    200, 0, 2.0, "#Delta R^{Ele}");

  HistClass::CreateHisto("Hmass_tau",  250, 0, 1000,  "M_{H} [GeV]");
  HistClass::CreateHisto("Hmass_tau_vis",  250, 0, 1000,  "M_{H} [GeV]");
  HistClass::CreateHisto("Hmass_tau_vis_selEvt",  250, 0, 1000,  "M_{H} [GeV]");
  HistClass::CreateHisto("deltaRTauTauReco_selEvt",  200, 0, 2.0, "#Delta R^{Tau_{h}}");

  HistClass::CreateHisto("ZPrimemass", 100, 0, 4000, "M_{X} [GeV]");
  HistClass::CreateHisto("ZPrimemass_vis", 100, 0, 4000, "M_{X} [GeV]");

  
  ///// GEN LEVEL HIST ///////
  HistClass::CreateHisto("Gen_MuonFromZ_Pt_Lead",    250, 0, 1000, "p_{T}^{GenMu} lead [GeV]");
  HistClass::CreateHisto("Gen_MuonFromZ_Pt_SubLead", 250, 0, 1000, "p_{T}^{GenMu} sublead [GeV]");
  HistClass::CreateHisto("Gen_MuonFromZ_Eta_Lead",    250, -2.5, 2.5, "#eta^{GenMu} lead");
  HistClass::CreateHisto("Gen_MuonFromZ_Eta_SubLead", 250, -2.5, 2.5, "#eta^{GenMu} sublead");
  HistClass::CreateHisto("Gen_MuonFromZ_Phi_Lead",    350, -3.5, 3.5, "#phi^{GenMu} lead");
  HistClass::CreateHisto("Gen_MuonFromZ_Phi_SubLead", 350, -3.5, 3.5, "#phi^{GenMu} sublead");
  HistClass::CreateHisto("Gen_DiMuonFromZ_dR",   200, 0, 2.0, "#Delta R^{GenMu}");

  HistClass::CreateHisto("Gen_EleFromZ_Pt_Lead",    250, 0, 1000, "p_{T}^{GenEle} lead [GeV]");
  HistClass::CreateHisto("Gen_EleFromZ_Pt_SubLead", 250, 0, 1000, "p_{T}^{GenEle} sublead [GeV]");
  HistClass::CreateHisto("Gen_EleFromZ_Eta_Lead",    250, -2.5, 2.5, "#eta^{GenEle} lead");
  HistClass::CreateHisto("Gen_EleFromZ_Eta_SubLead", 250, -2.5, 2.5, "#eta^{GenEle} sublead");
  HistClass::CreateHisto("Gen_EleFromZ_Phi_Lead",    350, -3.5, 3.5, "#phi^{GenEle} lead");
  HistClass::CreateHisto("Gen_EleFromZ_Phi_SubLead", 350, -3.5, 3.5, "#phi^{GenEle} sublead");
  HistClass::CreateHisto("Gen_DiEleFromZ_dR",   200, 0, 2.0, "#Delta R^{GenEle}");

  HistClass::CreateHisto("Gen_TauAllFromH_Pt_Lead",     250, 0, 1000, "p_{T}^{GenTauAll} lead [GeV]");
  HistClass::CreateHisto("Gen_TauAllFromH_Pt_SubLead",  250, 0, 1000, "p_{T}^{GenTauAll} sublead [GeV]");
  HistClass::CreateHisto("Gen_TauAllFromH_Eta_Lead",    250, -2.5, 2.5, "#eta^{GenTauAll} lead");
  HistClass::CreateHisto("Gen_TauAllFromH_Eta_SubLead", 250, -2.5, 2.5, "#eta^{GenTaauAll} sublead");
  HistClass::CreateHisto("Gen_TauAllFromH_Phi_Lead",    350, -3.5, 3.5, "#phi^{GenTauAll} lead");
  HistClass::CreateHisto("Gen_TauAllFromH_Phi_SubLead", 350, -3.5, 3.5, "#phi^{GenTauAll} sublead");
  HistClass::CreateHisto("Gen_DiTauAllFromH_dR",        200, 0, 2.0, "#Delta R^{GenTauAll}");


  for(auto syst : syst_shifter->m_activeSystematics){
    if(syst->m_particleType=="met"){
      syst->m_particleType=m_METType;
    }
  }
}

specialAna::~specialAna() {
}

void specialAna::analyseEvent( const pxl::Event* event ) {
  initEvent( event );
  //std::cout << "\nEVT" << std::endl;

  if(not runOnData){
    if (debug) std::cout << "runOnData=" << runOnData << "   Will call Fill_Gen_Controll_histo()" << std::endl;
    Fill_Gen_Controll_histo();
    if (debug) std::cout << "Called it successfully" << std::endl;
    ////    GeneratorSelection();
  }
  

  pxl::EventView *RecEvtView = event->getObjectOwner().findObject< pxl::EventView >( "Rec" );

  bool isRec = (RecEvtView->getUserRecord("Type").asString() == "Rec");
  
  //m_eventSelector;
  std::map< std::string, std::vector< pxl::Particle* > > particleMap = m_eventSelector->getParticleLists( RecEvtView, isRec );
  std::vector< pxl::Particle* > muoList = particleMap["Muon"];
  std::vector< pxl::Particle* > eleList = particleMap["Ele"];
  std::vector< pxl::Particle* > tauList = particleMap["Tau"];
  std::vector< pxl::Particle* > metList = particleMap["MET"];
  
  if (debug) std::cout << "metList.size() " << metList.size() << "  tauList.size() " << tauList.size() << std::endl;
  
  HistClass::Fill("nTau",  tauList.size(), weight);
  HistClass::Fill("nMuon", muoList.size(), weight);
  HistClass::Fill("nEle",  eleList.size(), weight);
  HistClass::Fill("nLep",  (eleList.size()+muoList.size()), weight);

  //   std::vector< pxl::Particle* > zCandidateList;
  //  pxl::Particle*  RecoMu = 0;
  pxl::Particle*  RecoMu1 = 0;
  pxl::Particle*  RecoMu2 = 0;
  pxl::Particle*  RecoEle1 = 0;
  pxl::Particle*  RecoEle2 = 0;
  pxl::Particle*  RecoTau1 = 0;
  pxl::Particle*  RecoTau2 = 0;

  pxl::Particle*  zCandidate_Mu = 0;
  pxl::Particle*  zCandidate_Ele = 0;
  pxl::Particle*  zCandidate = 0;
  pxl::Particle*  HiggsCandidate_Tau = 0;
  pxl::Particle*  HiggsCandidate_Tau_vis = 0;
  pxl::Particle*  ZPrime_Candidate = 0;
  pxl::Particle*  ZPrime_Candidate_vis = 0;

  double zmass = 91.2;
  double bestdiff_mu = 10000000.0;
  double bestdiff_ele = 10000000.0;
  double bestdiff_tau = 10000000.0;

  double Hmass = 125.0;
  int ZeeChannel = 0;
  int ZmumuChannel = 0;

  
  if (TriggerSelector(event)) {
    //std::cout << "Trigger passed " << std::endl;
    /////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////
    if (muoList.size() > 1){
      //~ for( auto& muo : muoList ){
      for( vector< pxl::Particle* >::const_iterator muon1 = muoList.begin(); muon1 != muoList.end(); muon1++ ){
	pxl::Particle* mu1 = RecEvtView->getObjectOwner().create< pxl::Particle >( *muon1 );
	for( vector< pxl::Particle* >::const_iterator muon2 = muon1+1; muon2 != muoList.end(); ++muon2 ){
	  pxl::Particle* mu2 = RecEvtView->getObjectOwner().create< pxl::Particle >( *muon2 );
	  if ( ((Check_Muo_ID_Strict(mu1)) &&  (Check_Muo_ID_Mild(mu2))) || ((Check_Muo_ID_Strict(mu2)) &&  (Check_Muo_ID_Mild(mu1)))  ) {
	    //std::cout << "mu1->getCharge() = " << mu1->getCharge() << std::endl;
	    if ( (mu1->getCharge()+mu2->getCharge())==0 ) {
	      //pxl::Particle* zcand_mu = RecEvtView->getObjectOwner().create< pxl::Particle >( *muon1 );
	      pxl::Particle* zcand_mu = (pxl::Particle*) mu1->clone();
	      zcand_mu->addP4( mu2->getVector() );
	      // zCandidate_MuList.push_back( zcand );
	      //std::cout << "Mu : Will use bestdiff = " << bestdiff_mu << std::endl;
	      if( fabs( zcand_mu->getMass() - zmass) < bestdiff_mu) {
		zCandidate_Mu = zcand_mu;
		RecoMu1 = mu1;
		RecoMu2 = mu2;
		bestdiff_mu = fabs( zcand_mu->getMass() - zmass);
	      }
	    }
	  }		
	}
      }
      if (zCandidate_Mu != 0) {
	
	//std::cout << "**Found z mass Muon** " << zCandidate_Mu->getMass() << std::endl;
	MyRecoMu->push_back(RecoMu1);
	MyRecoMu->push_back(RecoMu2);
	std::sort (MyRecoMu->begin(), MyRecoMu->end(), ptsort_pxl);
	
	HistClass::Fill("Zmass_mu",         zCandidate_Mu->getMass(), weight);
	HistClass::Fill("Muon_Pt_Lead",     MyRecoMu->at(0)->getPt(), weight);
	HistClass::Fill("Muon_Pt_SubLead",  MyRecoMu->at(1)->getPt(), weight);
      }  
    }
    
    // Electron //
    if (eleList.size() > 1){
      for( vector< pxl::Particle* >::const_iterator ele1 = eleList.begin(); ele1 != eleList.end(); ele1++ ){
	pxl::Particle* e1 = RecEvtView->getObjectOwner().create< pxl::Particle >( *ele1 );
	//if ( ((e1->getUserRecord("ISOfailed").asBool()) == true) ||  ((e1->getUserRecord("IDpassed").asBool()) == true)   ) {
	for( vector< pxl::Particle* >::const_iterator ele2 = ele1+1; ele2 != eleList.end(); ++ele2 ){
	  pxl::Particle* e2 = RecEvtView->getObjectOwner().create< pxl::Particle >( *ele2 );
	  // if ( ((e2->getUserRecord("ISOfailed").asBool()) == true) || ((e2->getUserRecord("IDpassed").asBool()) == true) ) {
	  if (debug) std::cout << "e1 " <<  e1->getUserRecord("ISOfailed").asBool() << "  " <<  e1->getUserRecord("IDpassed").asBool() << std::endl;
	  if (debug) std::cout << "e2 " <<  e2->getUserRecord("ISOfailed").asBool() << "  " <<  e2->getUserRecord("IDpassed").asBool() << std::endl;

	  if ( ((Check_Ele_ID(e1)) &&  (Check_Ele_ID(e2))) ) {
	    // std::cout << "PASSED ele ID for both ele" << std::endl;
	    if ( (e1->getCharge()+e2->getCharge())==0 ){
	      pxl::Particle* zcand_ele = (pxl::Particle*) e1->clone();
	      zcand_ele->addP4( e2->getVector() );
	      //std::cout<< "Ele : Will use bestdiff = " << bestdiff_ele << std::endl;
	      if( fabs( zcand_ele->getMass() - zmass) < bestdiff_ele) {
		zCandidate_Ele = zcand_ele;
		RecoEle1 = e1;
		RecoEle2 = e2;
		bestdiff_ele = fabs( zcand_ele->getMass() - zmass);
	      }
	    }
	  }
	}
      }
      if (zCandidate_Ele != 0) {
	//std::cout << "**Found z mass ele** " << zCandidate_Ele->getMass() << std::endl;
	MyRecoEle->push_back(RecoEle1);
        MyRecoEle->push_back(RecoEle2);
	std::sort (MyRecoEle->begin(), MyRecoEle->end(), ptsort_pxl);

	HistClass::Fill("Zmass_ele",        zCandidate_Ele->getMass(), weight);
	HistClass::Fill("Ele_Pt_Lead",      MyRecoEle->at(0)->getPt(), weight);
	HistClass::Fill("Ele_Pt_SubLead",   MyRecoEle->at(1)->getPt(), weight);
      }
    }
    
    //  std::cout << "zCandidate_Mu=" << zCandidate_Mu << "   zCandidate_Ele=" << zCandidate_Ele << std::endl; 
    if (zCandidate_Mu && zCandidate_Ele) {
      //  std::cout << "z cand found from both mu and ele collection. Have to choose one." << std::endl;
      if  ( (fabs(zmass-zCandidate_Mu->getMass())) < (fabs(zmass-zCandidate_Ele->getMass())) ) {  
	zCandidate=zCandidate_Mu;
	ZmumuChannel=1;
      }
      else if ((fabs(zmass-zCandidate_Mu->getMass())) > (fabs(zmass-zCandidate_Ele->getMass())) ) {
	zCandidate=zCandidate_Ele;
	ZeeChannel=1;

      }
    }
    
    else if (zCandidate_Mu && !zCandidate_Ele) {
      zCandidate=zCandidate_Mu;
      ZmumuChannel=1;

    }
    else if (zCandidate_Ele && !zCandidate_Mu) {
      zCandidate=zCandidate_Ele;
      ZeeChannel=1;

    }
    
    if (zCandidate != 0){
      HistClass::Fill("Zmass_lep", zCandidate->getMass(), weight);
    }


    ////////////////////////////////////////////////////////////////////////////
    ///        USE ONLY VISIBLE MASS : DO NOT TAKE ET_MISS INTO ACCOUNT      ///
    ////////////////////////////////////////////////////////////////////////////
    if ( tauList.size() > 1 ) {
      for( vector< pxl::Particle* >::const_iterator tau1 = tauList.begin(); tau1 != tauList.end(); tau1++ ){
	pxl::Particle* t1 = RecEvtView->getObjectOwner().create< pxl::Particle >( *tau1 );
	//if ( Check_Tau_ID(t1) ) {
	for( vector< pxl::Particle* >::const_iterator tau2 = tau1+1; tau2 != tauList.end(); ++tau2 ){
	  pxl::Particle* t2 = RecEvtView->getObjectOwner().create< pxl::Particle >( *tau2 );
	  if ( (Check_Tau_ID(t1)) && (Check_Tau_ID(t2)) ) {
	    if ( (t1->getCharge()+t2->getCharge())==0 ) {
	      
	      //pxl::Particle* Hcand_tau_vis = RecEvtView->getObjectOwner().create< pxl::Particle >( *tau1 );
	      pxl::Particle* Hcand_tau_vis = (pxl::Particle*) t1->clone();
	      //Hcand_tau_vis->addP4( (*tau2)->getVector() );
	      Hcand_tau_vis->addP4( (t2)->getVector() );
	      //std::cout << "Tau : Will use bestdiff = " << bestdiff_tau << std::endl;
	      if( fabs( Hcand_tau_vis->getMass() - Hmass) < bestdiff_tau) {
		HiggsCandidate_Tau_vis = Hcand_tau_vis;
		RecoTau1 = t1;
                RecoTau2 = t2;
		bestdiff_tau = fabs( Hcand_tau_vis->getMass() - Hmass);
	      } 
	    }
	    //  }
	  }		
	}
      }
    }	     	      
    
    /*
    ////////////////////////////////////////////////////////////////////////////
    ////    Collinear Approximation ////  
    if ( (tauList.size() > 1) && (metList.size() > 0) ) {
      if (debug) std::cout << ">1 Taus and >0 Met Found" << std::endl;
      for( vector< pxl::Particle* >::const_iterator tau1 = tauList.begin(); tau1 != tauList.end(); tau1++ ){
	pxl::Particle* t1 = RecEvtView->getObjectOwner().create< pxl::Particle >( *tau1 );
	if ( Check_Tau_ID(t1) ) {
	  for( vector< pxl::Particle* >::const_iterator tau2 = tau1+1; tau2 != tauList.end(); ++tau2 ){
	    pxl::Particle* t2 = RecEvtView->getObjectOwner().create< pxl::Particle >( *tau2 );
	    if ( Check_Tau_ID(t2) ) {
	      if ( (t1->getCharge()+t2->getCharge())==0 ) {
		double EtMiss = metList.at(0)->getPt();
		double EtL1 = t1->getPt();
		double EtL2 = t2->getPt();
		double cosPhiMiss = TMath::Cos(metList.at(0)->getPhi());
		double sinPhiMiss = TMath::Sin(metList.at(0)->getPhi());
		double cosPhiL1 = TMath::Cos(t1->getPhi());
		double sinPhiL1 = TMath::Sin(t1->getPhi());
		double cosPhiL2 = TMath::Cos(t2->getPhi());
		double sinPhiL2 = TMath::Sin(t2->getPhi());
		
		if (debug) {
		  std::cout << "EtMiss=" << EtMiss << " EtL1=" << EtL1 << " EtL2=" << EtL2 << " cosPhiMiss=" << cosPhiMiss
			    << " sinPhiMiss=" << sinPhiMiss << " cosPhiL1=" << cosPhiL1 << " sinPhiL1=" << sinPhiL1
			    << " cosPhiL2=" << cosPhiL2 << " sinPhiL2=" << sinPhiL2 << std::endl;
		}
		double a = (EtMiss/EtL1) * ( (cosPhiMiss*sinPhiL2 - sinPhiMiss*cosPhiL2)/ (cosPhiL1*sinPhiL2 - sinPhiL1*cosPhiL2) ) ;
		double b = (EtMiss/EtL2) * ( (cosPhiMiss*sinPhiL1 - sinPhiMiss*cosPhiL1)/ (cosPhiL2*sinPhiL1 - sinPhiL2*cosPhiL1) ) ;
		if (debug) std::cout << "a=" << a << "  b=" << b << std::endl;
		if (a<0.0 || b<0.0) continue;
		pxl::Particle* Hcand_tau = RecEvtView->getObjectOwner().create< pxl::Particle >( *tau1 );
		
		Hcand_tau->addP4( (*tau2)->getVector() );
		
		if (debug) std::cout << "Before adding a, pt=" << Hcand_tau->getPt() << " et=" << Hcand_tau->getEt() << " m=" << Hcand_tau->getMass() << std::endl;
	        Hcand_tau->addP4( ((*tau1)->getVector())*a );
		if (debug) std::cout << "After adding a, pt=" << Hcand_tau->getPt() << " et=" << Hcand_tau->getEt() << " m=" << Hcand_tau->getMass() << std::endl;
		Hcand_tau->addP4( ((*tau2)->getVector())*b );
		if (debug) std::cout << "After adding b, pt=" << Hcand_tau->getPt() << " et=" << Hcand_tau->getEt() << " m=" << Hcand_tau->getMass() << std::endl;

		if( fabs( Hcand_tau->getMass() - Hmass) < bestdiff) {
		  HiggsCandidate_Tau = Hcand_tau;
		  bestdiff = fabs( Hcand_tau->getMass() - Hmass);
		}
	      }
	    }
	  }		
	}
      }
    }
    
    if (HiggsCandidate_Tau != 0) {
      if (debug) {
	std::cout << "**Found H mass Tau** " << HiggsCandidate_Tau->getMass() << std::endl;
	std::cout << "Higgs-> Px=" << HiggsCandidate_Tau->getPx() << " Py=" << HiggsCandidate_Tau->getPy() << 
	  " Pz=" << HiggsCandidate_Tau->getPz() << " E=" << HiggsCandidate_Tau->getE() << 
	  " pT=" << HiggsCandidate_Tau->getPt() << " Et=" << HiggsCandidate_Tau->getEt() << std::endl;
      }
      HistClass::Fill("Hmass_tau", HiggsCandidate_Tau->getMass(), weight);
    }
    */
    if (HiggsCandidate_Tau_vis != 0) {
      if (debug) {
      std::cout << "**Found H mass Tau visible** " << HiggsCandidate_Tau_vis->getMass() << std::endl;
      std::cout << "HiggsVis-> Px=" << HiggsCandidate_Tau_vis->getPx() << " Py=" << 
	HiggsCandidate_Tau_vis->getPy() << " Pz=" << HiggsCandidate_Tau_vis->getPz() 
		<< " E=" << HiggsCandidate_Tau_vis->getE() <<
	" pT=" << HiggsCandidate_Tau_vis->getPt() << " Et=" << HiggsCandidate_Tau_vis->getEt()<< std::endl;
      }
      HistClass::Fill("Hmass_tau_vis", HiggsCandidate_Tau_vis->getMass(), weight);
      MyRecoTau->push_back(RecoTau1);
      MyRecoTau->push_back(RecoTau2);
      std::sort (MyRecoTau->begin(), MyRecoTau->end(), ptsort_pxl);

    }  


    /*    
    if ( zCandidate  &&  HiggsCandidate_Tau ) {
      if (zCandidate->getMass()>70.0 && zCandidate->getMass()<110.0) {
	if (debug) std::cout << "SIGNAL-LIKE EVENT FOUND !!!!!" << std::endl ;
	ZPrime_Candidate = zCandidate;
	ZPrime_Candidate->addP4(HiggsCandidate_Tau->getVector() );
	if (debug) std::cout << "Zprime mass=" << ZPrime_Candidate->getMass() << std::endl;
	HistClass::Fill("ZPrimemass", ZPrime_Candidate->getMass(), weight);
      }
    }
    */

    if ( zCandidate  &&  HiggsCandidate_Tau_vis ) {
      if (zCandidate->getMass()>60.0 && zCandidate->getMass()<120.0) {
	if (debug) std::cout << "SIGNAL-LIKE EVENT FOUND (vis) !!!!!" << std::endl ;
	//ZPrime_Candidate_vis->addP4(zCandidate->getVector() );
	ZPrime_Candidate_vis = (pxl::Particle*) zCandidate->clone();
	ZPrime_Candidate_vis->addP4(HiggsCandidate_Tau_vis->getVector() );
	if (debug) std::cout << "Zprime mass vis=" << ZPrime_Candidate_vis->getMass() << std::endl;
	HistClass::Fill("ZPrimemass_vis", ZPrime_Candidate_vis->getMass(), weight);
	HistClass::Fill("Hmass_tau_vis_selEvt", HiggsCandidate_Tau_vis->getMass(), weight);
	if (ZmumuChannel>0) { 
	  //std::cout << "ZmumuChannel = " << ZmumuChannel <<  "  mass = " << zCandidate->getMass()  << " , " << ZPrime_Candidate_vis->getMass()  << std::endl;
	  HistClass::Fill("Zmass_mu_selEvt", zCandidate->getMass(), weight);
	  double deltaRmumuReco = DeltaR(MyRecoMu->at(0),MyRecoMu->at(1));
	  HistClass::Fill("deltaRmumuReco_selEvt", deltaRmumuReco, weight);

	}
        if (ZeeChannel>0) {
	  //std::cout << "ZeeChannel = " << ZeeChannel <<  "  mass = " << zCandidate->getMass() << " , " << ZPrime_Candidate_vis->getMass()  << std::endl;
	  HistClass::Fill("Zmass_ele_selEvt", zCandidate->getMass(), weight);
	  double deltaReeReco= DeltaR(MyRecoEle->at(0),MyRecoEle->at(1));
	  HistClass::Fill("deltaReeReco_selEvt", deltaReeReco, weight);
	}

	HistClass::Fill("Zmass_lep_selEvt", zCandidate->getMass(), weight);
	HistClass::Fill("Z_pt_selEvt", zCandidate->getPt(), weight);
	HistClass::Fill("H_pt_selEvt", HiggsCandidate_Tau_vis->getPt(), weight);

	double deltaRTauTauReco = DeltaR(MyRecoTau->at(0),MyRecoTau->at(1));
	HistClass::Fill("deltaRTauTauReco_selEvt", deltaRTauTauReco, weight);

      }
    }
 }

 if (debug) std::cout << "Will call endEvent()" << std::endl;
 endEvent( event );
 
}

void specialAna::endEvent( const pxl::Event* event ){
  if( not runOnData ){

    delete EleListGen;
    delete MuonListGen;
    delete GammaListGen;
    delete METListGen;
    delete JetListGen;
    delete TauListGen;
    delete S3ListGen;        
    delete SelectedGen;     
    delete MyGenMu;
    delete MyGenEle;
    delete MyGenTauAll;
    delete MyRecoMu;
    delete MyRecoEle;
    delete MyRecoTau;

    EleListGen = 0;
    MuonListGen = 0;
    GammaListGen = 0;
    METListGen = 0;
    JetListGen = 0;
    TauListGen = 0;
    S3ListGen = 0;
    SelectedGen = 0;
    MyGenMu = 0;
    MyGenEle = 0;
    MyGenTauAll = 0;
    MyRecoMu = 0;
    MyRecoEle = 0;
    MyRecoTau = 0;

  }

}


void specialAna::endJob( const Serializable* ) { 
  //std::cout << "Inside endJob" << std::endl;		
  file1->cd();
  HistClass::WriteAll("counters");
  file1->mkdir("RECO");
  file1->cd("RECO/");
  HistClass::Write("nTau");
  HistClass::Write("nMuon");
  HistClass::Write("nEle");
  HistClass::Write("nLep");

  HistClass::Write("Muon_Pt_Lead");
  HistClass::Write("Muon_Pt_SubLead");
  HistClass::Write("Ele_Pt_Lead");
  HistClass::Write("Ele_Pt_SubLead");

  HistClass::Write("Zmass_mu");
  HistClass::Write("Zmass_ele");
  HistClass::Write("Zmass_lep");

  HistClass::Write("Zmass_mu_selEvt");
  HistClass::Write("Zmass_ele_selEvt");
  HistClass::Write("Zmass_lep_selEvt");

  HistClass::Write("Z_pt_selEvt");
  HistClass::Write("H_pt_selEvt");

  HistClass::Write("deltaRmumuReco_selEvt");
  HistClass::Write("deltaReeReco_selEvt");
  HistClass::Write("deltaRTauTauReco_selEvt");

  HistClass::Write("Hmass_tau");
  HistClass::Write("Hmass_tau_vis");
  HistClass::Write("Hmass_tau_vis_selEvt");
  HistClass::Write("ZPrimemass");
  HistClass::Write("ZPrimemass_vis");


  if (debug) std::cout << "will create GEN folder" << std::endl;
  file1->mkdir("GEN");
  file1->cd("GEN/");
  HistClass::Write("Gen_MuonFromZ_Pt_Lead");
  HistClass::Write("Gen_MuonFromZ_Pt_SubLead");
  HistClass::Write("Gen_MuonFromZ_Eta_Lead");
  HistClass::Write("Gen_MuonFromZ_Eta_SubLead");
  HistClass::Write("Gen_MuonFromZ_Phi_Lead");
  HistClass::Write("Gen_MuonFromZ_Phi_SubLead");

  HistClass::Write("Gen_EleFromZ_Pt_Lead");
  HistClass::Write("Gen_EleFromZ_Pt_SubLead");
  HistClass::Write("Gen_EleFromZ_Eta_Lead");
  HistClass::Write("Gen_EleFromZ_Eta_SubLead");
  HistClass::Write("Gen_EleFromZ_Phi_Lead");
  HistClass::Write("Gen_EleFromZ_Phi_SubLead");

  HistClass::Write("Gen_TauAllFromH_Pt_Lead");
  HistClass::Write("Gen_TauAllFromH_Pt_SubLead");
  HistClass::Write("Gen_TauAllFromH_Eta_Lead");
  HistClass::Write("Gen_TauAllFromH_Eta_SubLead");
  HistClass::Write("Gen_TauAllFromH_Phi_Lead");
  HistClass::Write("Gen_TauAllFromH_Phi_SubLead");

  HistClass::Write("Gen_DiMuonFromZ_dR");
  HistClass::Write("Gen_DiEleFromZ_dR");
  HistClass::Write("Gen_DiTauAllFromH_dR");


  if (debug) std::cout << "DONE" << std::endl;

  file1->Close();
  
  delete file1;
  
}


void specialAna::initEvent( const pxl::Event* event ){
  HistClass::Fill("h_counters", 1, 1); // increment number of events
  events_++;
  m_GenEvtView = event->getObjectOwner().findObject< pxl::EventView >( "Gen" );
  EleListGen     = new std::vector< pxl::Particle* >;
  MuonListGen    = new std::vector< pxl::Particle* >;
  GammaListGen   = new std::vector< pxl::Particle* >;
  METListGen     = new std::vector< pxl::Particle* >;
  JetListGen     = new std::vector< pxl::Particle* >;
  TauListGen     = new std::vector< pxl::Particle* >;
  TauVisListGen  = new std::vector< pxl::Particle* >;
  S3ListGen      = new std::vector< pxl::Particle* >;
    
  SelectedGen    = new std::vector< pxl::Particle* >;
  MyGenMu        = new std::vector< pxl::Particle* >;
  MyGenEle       = new std::vector< pxl::Particle* >;
  MyGenTauAll    = new std::vector< pxl::Particle* >;

  MyRecoMu        = new std::vector< pxl::Particle* >;
  MyRecoEle       = new std::vector< pxl::Particle* >;
  MyRecoTau       = new std::vector< pxl::Particle* >;

  if (event->getObjectOwner().findObject< pxl::EventView >("Trig")) {
    m_TrigEvtView = event->getObjectOwner().findObject< pxl::EventView >("Trig");
  } else {
    m_TrigEvtView = event->getObjectOwner().findObject< pxl::EventView >("Rec");
  }

  vector< pxl::Particle* > AllParticlesGen;
  m_GenEvtView->getObjectsOfType< pxl::Particle >( AllParticlesGen );
  pxl::sortParticles( AllParticlesGen );
  // push them into the corresponding vectors
  string genCollection="gen";
        
  for( vector< pxl::Particle* >::const_iterator part_it = AllParticlesGen.begin(); part_it != AllParticlesGen.end(); ++part_it ) {
    pxl::Particle *part = *part_it;
    string Name = part->getName();
    // Only fill the collection if we want to use the particle!
    if(      Name == "Muon"    ) MuonListGen->push_back( part );
    else if( Name == "Ele"     ) EleListGen->push_back( part );
    else if( Name == "Gamma"   ) GammaListGen->push_back( part );
    else if( Name == "Tau"     ) TauListGen->push_back( part );
    else if( Name == (m_METType+"_gen") ) METListGen->push_back( part );
    else if( Name == m_JetAlgo ) JetListGen->push_back( part );
    else if( Name == genCollection) {
      S3ListGen->push_back( part );
    }
  }
  

  weight = 1.;

  event_weight = 1;
  pileup_weight = 1;

  if (not runOnData) {

    event_weight = m_GenEvtView->getUserRecord("Weight");
    if (debug) std::cout << "event_weight=" << event_weight << std::endl;
    // event_weight = 1;
    // double varKfactor_weight = m_GenEvtView->getUserRecord_def( "kfacWeight",1. );
    pileup_weight = m_GenEvtView->getUserRecord_def("PUWeight", 1.);
    if (debug) std::cout << "pileup_weight=" << pileup_weight << std::endl;

    weight = event_weight * pileup_weight;
    if (debug) std::cout << "weight=" << weight << std::endl;

    if (doSampleWeighting) {
      std::string datastream = event->getUserRecord("Dataset").asString();
      TString Datastream = datastream;
      //      std::cout << "Datastream = " << Datastream << std::endl;
      sample_weight = 1.;   
    }
  }
  
  weight = weight * sample_weight;

} 
 
  
bool specialAna::ptsort_pxl (pxl::Particle *i, pxl::Particle *j) { return ( (i->getPt()) > (j->getPt()) ); }
//bool specialAna::ptsort (double i, double j) { return (i<j); }


bool specialAna::Check_Muo_ID_Strict(pxl::Particle* muon, bool do_pt_cut, bool do_eta_cut) {
  bool muon_ID = muon->getUserRecord("isHighPtMuon").asBool() ? true : false;
  //bool muon_ISO = muon -> getUserRecord("IsoR3SumPt").asFloat() / muon -> getPt() < 0.1 ? true : false;
  bool muon_eta = TMath::Abs(muon -> getEta()) < 2.4 ? true : false;
  bool muon_pt = muon -> getPt() > 25. ? true : false;
  if (not do_pt_cut) {
    muon_pt = true;
  }
  if (not do_eta_cut) {
    muon_eta = true;
  }
  if (muon_ID && muon_eta && muon_pt) return true;
  return false;
}

bool specialAna::Check_Muo_ID_Mild(pxl::Particle* muon, bool do_pt_cut, bool do_eta_cut) {

  bool muon_ID = muon->getUserRecord("isTrackerMuon").asBool() ? true : false;
    //bool muon_ISO = muon -> getUserRecord("IsoR3SumPt").asFloat() / muon -> getPt() < 0.1 ? true : false;
  bool muon_eta = TMath::Abs(muon -> getEta()) < 2.4 ? true : false;
  bool muon_pt = muon -> getPt() > 25.0 ? true : false;
  if (not do_pt_cut) {
    muon_pt = true;
  }
  if (not do_eta_cut) {
    muon_eta = true;
  }
  if (muon_ID && muon_eta && muon_pt) return true;
  return false;
}

bool specialAna::Check_Ele_ID(pxl::Particle* ele, bool do_pt_cut, bool do_eta_cut) {
  bool ele_ID = false;
  if ( ((ele->getUserRecord("ISOfailed").asBool()) == true) ||  ((ele->getUserRecord("IDpassed").asBool()) == true) ) ele_ID=true; 
  bool ele_eta = TMath::Abs(ele-> getEta()) < 2.5 ? true : false;
  bool ele_pt = ele-> getPt() > 35.0 ? true : false;
  if (not do_pt_cut) {
    ele_pt = true;
  }
  if (not do_eta_cut) {
    ele_eta = true;
  }
  if (ele_ID && ele_eta && ele_pt) return true;
  return false;

}

bool specialAna::Check_Tau_ID(pxl::Particle* tau, bool do_pt_cut, bool do_eta_cut) {
  bool tau_ID = tau->getUserRecord("decayModeFindingNewDMs").asFloat() >= 0.5  ?  true : false; 
  bool tau_ELE = tau->getUserRecord("againstElectronLooseMVA5").asFloat() >= 0.5  ?  true : false;
  bool tau_MUO = tau->getUserRecord("againstMuonLoose3").asFloat() >= 0.5  ?  true : false;
  bool tau_eta = TMath::Abs(tau -> getEta()) < 2.3 ? true : false;
  bool tau_pt = tau -> getPt() > 20. ? true : false;
  
  if (tau_ID && tau_ELE && tau_MUO && tau_eta && tau_pt) return true;
  return false;
}


bool specialAna::TriggerSelector(const pxl::Event* event) {
  bool triggered = false;
  for (auto const it : m_trigger_string) {
    // triggered = m_TrigEvtView->hasUserRecord(it) ? true : false;
    for (auto us : m_TrigEvtView->getUserRecords()) {
       if (std::string::npos != us.first.find(it)) {
	triggered = true;
	triggers.insert(us.first);
      }
    }
  }
  return (triggered);
}

void specialAna::Fill_Gen_Controll_histo() {
  if (debug) std::cout << "Inside Fill_Gen_Controll_histo" << std::endl;
  if (debug) std::cout << "S3ListGen->size() = " << S3ListGen->size() << std::endl;
  for(uint i = 0; i < S3ListGen->size(); i++){
    if (debug) std::cout << "i=" << i << " S3ListGen->at(i)->getPdgNumber()=" << S3ListGen->at(i)->getPdgNumber() <<  std::endl; 
    
    // PDG ID //
    // Muon : 13
    // Electron : 11
    // Tau : 15
    // Z : 23
    // Higgs : 25

    if( (TMath::Abs(S3ListGen->at(i)->getPdgNumber()) == 13)  &&  (TMath::Abs( (( pxl::Particle*)  S3ListGen->at(i)->getMother())->getPdgNumber()) == 23) ) {
      if (debug) std::cout << "Found muon with pt " << S3ListGen->at(i)->getPt()  << std::endl ;
      MyGenMu->push_back(S3ListGen->at(i));
    }
    
    if( (TMath::Abs(S3ListGen->at(i)->getPdgNumber()) == 11)  &&  (TMath::Abs( (( pxl::Particle*)  S3ListGen->at(i)->getMother())->getPdgNumber()) == 23) ) {
      MyGenEle->push_back(S3ListGen->at(i));
    }

    /////
    if( (TMath::Abs(S3ListGen->at(i)->getPdgNumber()) == 15)  &&  (TMath::Abs( (( pxl::Particle*)  S3ListGen->at(i)->getMother())->getPdgNumber()) == 25) ) {
      //std::cout << "\n Found a tau from Higgs!!  No. of tau daughters : " << S3ListGen->at(i)->getDaughters().size() << std::endl;
      MyGenTauAll->push_back(S3ListGen->at(i));
      
      for (std::set< pxl::Relative* >::const_iterator part_it = S3ListGen->at(i)->getDaughters().begin(); part_it != S3ListGen->at(i)->getDaughters().end(); ++part_it) {
	pxl::Relative *part_i = *part_it;
	pxl::Particle *part_tmp = (pxl::Particle*)part_i;
	//std::cout << "Tau daughter is -> " << part_tmp->getPdgNumber() << std::endl; 
	//if ( (TMath::Abs(part_tmp->getPdgNumber()) == 13) ||  (TMath::Abs(part_tmp->getPdgNumber()) == 11) ) continue;
	//HistClass::Fill("GenTauHadronicPt",S3ListGen->at(i)->getPt(),1);
      }
    }    
    /////
  }

  //  std::sort (mu_pt->begin(), mu_pt->end(), ptsort);
  std::sort (MyGenMu->begin(), MyGenMu->end(), ptsort_pxl);
  std::sort (MyGenEle->begin(), MyGenEle->end(), ptsort_pxl);
  std::sort (MyGenTauAll->begin(), MyGenTauAll->end(), ptsort_pxl);


  if (debug) std::cout << "size of MyGenMu = " << MyGenMu->size() << std::endl ;
  if (debug) std::cout << "size of MyGenTauAll = " << MyGenTauAll->size() << std::endl ;

  if (MyGenMu->size()==2) {
    if (debug) std::cout << "mu pt1 = " << MyGenMu->at(0)->getPt() << " and mu pt2 = " << MyGenMu->at(1)->getPt() << std::endl ;
    HistClass::Fill("Gen_MuonFromZ_Pt_Lead",    MyGenMu->at(0)->getPt(),1);
    HistClass::Fill("Gen_MuonFromZ_Pt_SubLead", MyGenMu->at(1)->getPt(),1);
    
    HistClass::Fill("Gen_MuonFromZ_Eta_Lead",    MyGenMu->at(0)->getEta(),1);
    HistClass::Fill("Gen_MuonFromZ_Eta_SubLead", MyGenMu->at(1)->getEta(),1);

    HistClass::Fill("Gen_MuonFromZ_Phi_Lead",    MyGenMu->at(0)->getPhi(),1);
    HistClass::Fill("Gen_MuonFromZ_Phi_SubLead", MyGenMu->at(1)->getPhi(),1);
    
    double mu_dR = DeltaR(MyGenMu->at(0),MyGenMu->at(1));
    HistClass::Fill("Gen_DiMuonFromZ_dR", mu_dR, 1);
    
  }

  if (MyGenEle->size()==2) {
    HistClass::Fill("Gen_EleFromZ_Pt_Lead",     MyGenEle->at(0)->getPt(),1);
    HistClass::Fill("Gen_EleFromZ_Pt_SubLead",  MyGenEle->at(1)->getPt(),1);

    HistClass::Fill("Gen_EleFromZ_Eta_Lead",    MyGenEle->at(0)->getEta(),1);
    HistClass::Fill("Gen_EleFromZ_Eta_SubLead", MyGenEle->at(1)->getEta(),1);

    HistClass::Fill("Gen_EleFromZ_Phi_Lead",    MyGenEle->at(0)->getPhi(),1);
    HistClass::Fill("Gen_EleFromZ_Phi_SubLead", MyGenEle->at(1)->getPhi(),1);

    double ele_dR = DeltaR(MyGenEle->at(0),MyGenEle->at(1));
    HistClass::Fill("Gen_DiEleFromZ_dR", ele_dR, 1);

  }

  if (MyGenTauAll->size()==2) {
    HistClass::Fill("Gen_TauAllFromH_Pt_Lead",     MyGenTauAll->at(0)->getPt(),1);
    HistClass::Fill("Gen_TauAllFromH_Pt_SubLead",  MyGenTauAll->at(1)->getPt(),1);

    HistClass::Fill("Gen_TauAllFromH_Eta_Lead",    MyGenTauAll->at(0)->getEta(),1);
    HistClass::Fill("Gen_TauAllFromH_Eta_SubLead", MyGenTauAll->at(1)->getEta(),1);

    HistClass::Fill("Gen_TauAllFromH_Phi_Lead",    MyGenTauAll->at(0)->getPhi(),1);
    HistClass::Fill("Gen_TauAllFromH_Phi_SubLead", MyGenTauAll->at(1)->getPhi(),1);

    double tauAll_dR = DeltaR(MyGenTauAll->at(0),MyGenTauAll->at(1));
    HistClass::Fill("Gen_DiTauAllFromH_dR", tauAll_dR, 1);

  }
}


double specialAna::DeltaPhi(pxl::Particle* p1, pxl::Particle* p2) {
  double a = p1->getPhi();
  double b = p2->getPhi();
  double temp = fabs(a-b);
  if (temp <= TMath::Pi()) {
    return temp;
  } else {
    return  2.*TMath::Pi() - temp;
  }
}

double specialAna::DeltaR(pxl::Particle* p1, pxl::Particle* p2) {
  double d_eta = TMath::Abs(p1->getEta() - p2->getEta());
  double d_phi = DeltaPhi(p1, p2);
  return sqrt(pow(d_eta, 2) + pow(d_phi, 2));
}
