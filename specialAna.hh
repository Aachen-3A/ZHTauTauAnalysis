#ifndef specialAna_hh
#define specialAna_hh

#include <iostream>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <fstream>

///clean up the header!!!
#include "Pxl/Pxl/interface/pxl/core.hh"
#include "Pxl/Pxl/interface/pxl/hep.hh"
#include "Tools/PXL/Sort.hh"
//#include "Tools/Tools.hh"
#include "Tools/MConfig.hh"
#include "TF1.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TString.h"
#include "TLorentzVector.h"
#include <TFile.h>

#include "Main/Systematics.hh"
#include "Main/EventSelector.hh"


//----------------------------------------------------------------------
using namespace std;

//class Systematics;

class specialAna : public pxl::AnalysisProcess  {
public:
  specialAna( const Tools::MConfig &config,
	      Systematics *syst_shifter,
	      EventSelector *selector);
  virtual ~specialAna();
  
  virtual void endJob(const Serializable*);
  virtual void analyseEvent( const pxl::Event* event );
  bool Check_Muo_ID_Strict(pxl::Particle* muon, bool do_pt_cut = true, bool do_eta_cut = true);
  bool Check_Muo_ID_Mild(pxl::Particle* muon, bool do_pt_cut = true, bool do_eta_cut = true);
  bool Check_Ele_ID(pxl::Particle* ele, bool do_pt_cut = true, bool do_eta_cut = true);

  bool Check_Tau_ID(pxl::Particle* tau, bool do_pt_cut = true, bool do_eta_cut = true);
  bool TriggerSelector(const pxl::Event* event);
  double DeltaPhi(pxl::Particle* p1, pxl::Particle* p2);
  double DeltaR(pxl::Particle* p1, pxl::Particle* p2);
  TFile* file1;
   
  ofstream eventDisplayFile;
  std::stringstream eventsAfterCuts;
  std::stringstream eventsAfterCutsEvents;
  
  void Fill_Gen_Controll_histo( );
  
  void initEvent( const pxl::Event* event );
  void endEvent( const pxl::Event* event );
  //  static bool ptsort (double, double) ;
  static bool ptsort_pxl (pxl::Particle *i, pxl::Particle *j) ;


  EventSelector *m_eventSelector;
  Systematics *m_syst_shifter;
  
  pxl::EventView *m_RecEvtView;
  pxl::EventView *m_GenEvtView;
  pxl::EventView *m_TrigEvtView;
  pxl::EventView *m_FilterEvtView;

  int events_;
  bool doSampleWeighting;
  bool runOnData;
  bool debug;
  bool useSyst;
  std::string m_dataPeriod;
  string const m_JetAlgo, m_BJets_algo, m_METType, m_TauType;
  const std::vector< std::string > m_trigger_string;
  std::unordered_set< std::string > triggers;
  double weight;
  double sample_weight;
  double event_weight;
  double pileup_weight;

  //  std::vector< pxl::Particle* > * EleList;
  // std::vector< pxl::Particle* > * MuonList;
  // std::vector< pxl::Particle* > * TauList;
	

  std::vector< pxl::Particle* > * EleListGen;
  std::vector< pxl::Particle* > * MuonListGen;
  std::vector< pxl::Particle* > * TauListGen;
  std::vector< pxl::Particle* > * TauVisListGen;
  std::vector< pxl::Particle* > * GammaListGen;
  std::vector< pxl::Particle* > * METListGen;
  std::vector< pxl::Particle* > * JetListGen;
  std::vector< pxl::Particle* > * S3ListGen;
  std::vector< pxl::Particle* > * SelectedGen;

  std::vector< pxl::Particle* > * MyGenMu;
  std::vector< pxl::Particle* > * MyGenEle;
  std::vector< pxl::Particle* > * MyGenTauAll;

  std::vector< pxl::Particle* > * MyRecoMu;
  std::vector< pxl::Particle* > * MyRecoEle;
  std::vector< pxl::Particle* > * MyRecoTau;


};

#endif
