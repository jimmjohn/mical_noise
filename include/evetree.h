//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Sat Mar  9 03:36:31 2019 by ROOT version 5.34/36
// from TTree evetree/event tree
// found on file: RPCv4t_evtraw-24082017-094530.rre
//////////////////////////////////////////////////////////

#ifndef evetree_h
#define evetree_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

#include "TString.h"
#include "TTimeStamp.h"
#include "TBits.h"
#include <vector>

//#define DATA_ANALYSIS

class evetree {


 public :
  TTree          *fChain;   //!pointer to the analyzed TTree or TChain
  Int_t           fCurrent; //!current Tree number in a TChain

  static const int nlayer =12;
  static const int  nchannel=8;  //Number of TDC Channels

  // Declaration of leaf types from Geant4 Output
  static const unsigned int ndigihtmx=5000;
  UInt_t          irun;
  UInt_t          ievt;
  UInt_t          ngent;
  Int_t           pidin[1];   //[ngent]
  Float_t         ievt_wt;
  Int_t           intxn_id;
  Float_t         momin[1];   //[ngent]
  Float_t         thein[1];   //[ngent]
  Float_t         phiin[1];   //[ngent]
  Float_t         posxin[1];   //[ngent]
  Float_t         posyin[1];   //[ngent]
  Float_t         poszin[1];   //[ngent]
  UInt_t          ngenerated;
  UInt_t          naperture;
  UInt_t          ndigiht;
  Int_t           trigx;
  Int_t           trigy;
  UInt_t          stripid[ndigihtmx];   //[ndigiht]
  Int_t           digipdgid[ndigihtmx];   //[ndigiht]
  Int_t           digitime[ndigihtmx];   //[ndigiht]
  Int_t           digitruetime[ndigihtmx];   //[ndigiht]
  Float_t         digienr[ndigihtmx];   //[ndigiht]
  Float_t         digivx[ndigihtmx];   //[ndigiht]
  Float_t         digivy[ndigihtmx];   //[ndigiht]
  Float_t         digivz[ndigihtmx];   //[ndigiht]
  Float_t         digipx[ndigihtmx];   //[ndigiht]
  Float_t         digipy[ndigihtmx];   //[ndigiht]
  Float_t         digipz[ndigihtmx];   //[ndigiht]

  // Declaration of input data format
  Int_t           ENum[nlayer];
  Int_t           REnum[nlayer];
  ULong64_t       CEnum;
  TTimeStamp     *EveTS[nlayer];
  TBits          *xLayer[nlayer];
  TBits          *yLayer[nlayer];

  Int_t           tdc_ref_l[nlayer];
  Int_t           tdc_ref_t[nlayer];
  Int_t           trigCntDiff[nlayer];

  std::vector<unsigned int> *vxtdc_l[nlayer][nchannel];
  std::vector<unsigned int> *vytdc_l[nlayer][nchannel];

  std::vector<unsigned int> *vxtdc_t[nlayer][nchannel];
  std::vector<unsigned int> *vytdc_t[nlayer][nchannel];

#if !defined(DATA_ANALYSIS)
     // List of branches
     TBranch        *b_irun;   //!
     TBranch        *b_ievt;   //!
     TBranch        *b_ngent;   //!
     TBranch        *b_pidin;   //!
     TBranch        *b_ievt_wt;   //!
     TBranch        *b_intxn_id;   //!
     TBranch        *b_momin;   //!
     TBranch        *b_thein;   //!
     TBranch        *b_phiin;   //!
     TBranch        *b_posxin;   //!
     TBranch        *b_posyin;   //!
     TBranch        *b_poszin;   //!
     TBranch        *b_ngenerated;   //!
     TBranch        *b_naperture;   //!
     TBranch        *b_ndigiht;   //!
     TBranch        *b_trigx;   //!
     TBranch        *b_trigy;   //!
     TBranch        *b_stripid;   //!
     TBranch        *b_digipdgid;   //!
     TBranch        *b_digitime;   //!
     TBranch        *b_digitruetime;   //!
     TBranch        *b_digienr;   //!
     TBranch        *b_digivx;   //!
     TBranch        *b_digivy;   //!
     TBranch        *b_digivz;   //!
     TBranch        *b_digipx;   //!
     TBranch        *b_digipy;   //!
     TBranch        *b_digipz;   //!
#endif

#if defined(DATA_ANALYSIS)
  // List of branches
  TBranch        *b_ENum;   //!
  TBranch        *b_REnum;   //!
  TBranch        *b_CEnum;   //!
  TBranch        *b_Evetime[nlayer];   //!
  TBranch        *b_xstriphitsL[nlayer];   //!
  TBranch        *b_ystriphitsL[nlayer];   //!
  TBranch        *b_tdc_ref_l;   //!
  TBranch        *b_tdc_ref_t;   //!
  TBranch        *b_trigCntDiff;   //!
  TBranch        *b_xtdc_l[nlayer][nchannel];   //!
  TBranch        *b_ytdc_l[nlayer][nchannel];   //!
  TBranch        *b_xtdc_t[nlayer][nchannel];   //!
  TBranch        *b_ytdc_t[nlayer][nchannel];
#endif

  evetree(TTree *tree=0);
  virtual ~evetree();
  virtual Int_t    Cut(Long64_t entry);
  virtual Int_t    GetEntry(Long64_t entry);
  virtual Long64_t LoadTree(Long64_t entry);
  virtual void     Init(TTree *tree);
  virtual void     Loop();
  virtual Bool_t   Notify();
  virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef evetree_cxx
evetree::evetree(TTree *tree) : fChain(0)
{
  // if parameter tree is not specified (or zero), connect the file
  // used to generate this class and read the Tree.
  if (tree == 0) {
    TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("RPCv4t_evtraw-24082017-094530.rre");
    if (!f || !f->IsOpen()) {
      f = new TFile("RPCv4t_evtraw-24082017-094530.rre");
    }
    f->GetObject("evetree",tree);

  }
  Init(tree);
}

evetree::~evetree()
{
  if (!fChain) return;
  delete fChain->GetCurrentFile();
}

Int_t evetree::GetEntry(Long64_t entry)
{
  // Read contents of entry.
  if (!fChain) return 0;
  return fChain->GetEntry(entry);
}
Long64_t evetree::LoadTree(Long64_t entry)
{
  // Set the environment to read one entry
  if (!fChain) return -5;
  Long64_t centry = fChain->LoadTree(entry);
  if (centry < 0) return centry;
  if (fChain->GetTreeNumber() != fCurrent) {
    fCurrent = fChain->GetTreeNumber();
    Notify();
  }
  return centry;
}

void evetree::Init(TTree *tree)
{
  // The Init() function is called when the selector needs to initialize
  // a new tree or chain. Typically here the branch addresses and branch
  // pointers of the tree will be set.
  // It is normally not necessary to make changes to the generated
  // code, but the routine can be extended by the user if needed.
  // Init() will be called many times when running on PROOF
  // (once per file to be processed).

  // Set object pointer

  // Set branch addresses and branch pointers
  if (!tree) return;
  fChain = tree;
  fCurrent = -1;
  fChain->SetMakeClass(1);


  #if !defined(DATA_ANALYSIS)
     fChain->SetBranchAddress("irun", &irun, &b_irun);
     fChain->SetBranchAddress("ievt", &ievt, &b_ievt);
     fChain->SetBranchAddress("ngent", &ngent, &b_ngent);
     fChain->SetBranchAddress("pidin", pidin, &b_pidin);
     fChain->SetBranchAddress("ievt_wt", &ievt_wt, &b_ievt_wt);
     fChain->SetBranchAddress("intxn_id", &intxn_id, &b_intxn_id);
     fChain->SetBranchAddress("momin", momin, &b_momin);
     fChain->SetBranchAddress("thein", thein, &b_thein);
     fChain->SetBranchAddress("phiin", phiin, &b_phiin);
     fChain->SetBranchAddress("posxin", posxin, &b_posxin);
     fChain->SetBranchAddress("posyin", posyin, &b_posyin);
     fChain->SetBranchAddress("poszin", poszin, &b_poszin);
     fChain->SetBranchAddress("ngenerated", &ngenerated, &b_ngenerated);
     fChain->SetBranchAddress("naperture", &naperture, &b_naperture);
     fChain->SetBranchAddress("ndigiht", &ndigiht, &b_ndigiht);
     fChain->SetBranchAddress("trigx", &trigx, &b_trigx);
     fChain->SetBranchAddress("trigy", &trigy, &b_trigy);
     fChain->SetBranchAddress("stripid", stripid, &b_stripid);
     fChain->SetBranchAddress("digipdgid", digipdgid, &b_digipdgid);
     fChain->SetBranchAddress("digitime", digitime, &b_digitime);
     fChain->SetBranchAddress("digitruetime", digitruetime, &b_digitruetime);
     fChain->SetBranchAddress("digienr", digienr, &b_digienr);
     fChain->SetBranchAddress("digivx", digivx, &b_digivx);
     fChain->SetBranchAddress("digivy", digivy, &b_digivy);
     fChain->SetBranchAddress("digivz", digivz, &b_digivz);
     fChain->SetBranchAddress("digipx", digipx, &b_digipx);
     fChain->SetBranchAddress("digipy", digipy, &b_digipy);
     fChain->SetBranchAddress("digipz", digipz, &b_digipz);
  #endif

#if defined(DATA_ANALYSIS)
  for(int i=0; i<nlayer; i++) {
    ENum[i]=0;
    REnum[i]=0;
    tdc_ref_l[i]=0;
    tdc_ref_t[i]=0;
    trigCntDiff[i]=0;
    for(int j=0; j<nchannel; j++) {
      vxtdc_l[i][j]=0;
      vytdc_l[i][j]=0;
      vxtdc_t[i][j]=0;
      vytdc_t[i][j]=0;
    }
  }

  fChain->SetBranchAddress("ENum", ENum, &b_ENum);
  fChain->SetBranchAddress("REnum", REnum, &b_REnum);
  fChain->SetBranchAddress("CEnum", &CEnum, &b_CEnum);

  fChain->SetBranchAddress("tdc_ref_l", tdc_ref_l, &b_tdc_ref_l);
  fChain->SetBranchAddress("tdc_ref_t", tdc_ref_t, &b_tdc_ref_t);
  fChain->SetBranchAddress("trigCntDiff", trigCntDiff, &b_trigCntDiff);
  for(int ij=0;ij<nlayer;ij++) {
    EveTS[ij] = 0;
    xLayer[ij] = 0;
    yLayer[ij] = 0;
    fChain->SetBranchAddress(TString::Format("Evetime_%i",ij), &EveTS[ij], &b_Evetime[ij]);
    fChain->SetBranchAddress(TString::Format("xstriphitsL%i",ij), &xLayer[ij], &b_xstriphitsL[ij]);
    fChain->SetBranchAddress(TString::Format("ystriphitsL%i",ij), &yLayer[ij], &b_ystriphitsL[ij]);
    for(int jk=0;jk<nchannel;jk++) {
      fChain->SetBranchAddress(TString::Format("xtdc_l_%i_%i",ij,jk), &vxtdc_l[ij][jk], &b_xtdc_l[ij][jk]);
      fChain->SetBranchAddress(TString::Format("ytdc_l_%i_%i",ij,jk), &vytdc_l[ij][jk], &b_ytdc_l[ij][jk]);

      fChain->SetBranchAddress(TString::Format("xtdc_t_%i_%i",ij,jk), &vxtdc_t[ij][jk], &b_xtdc_t[ij][jk]);
      fChain->SetBranchAddress(TString::Format("ytdc_t_%i_%i",ij,jk), &vytdc_t[ij][jk], &b_ytdc_t[ij][jk]);
    }
  }
#endif


  Notify();
}

Bool_t evetree::Notify()
{
  // The Notify() function is called when a new file is opened. This
  // can be either for a new TTree in a TChain or when when a new TTree
  // is started when using PROOF. It is normally not necessary to make changes
  // to the generated code, but the routine can be extended by the
  // user if needed. The return value is currently not used.

  return kTRUE;
}

void evetree::Show(Long64_t entry)
{
  // Print contents of entry.
  // If entry is not specified, print current entry
  if (!fChain) return;
  fChain->Show(entry);
}
Int_t evetree::Cut(Long64_t entry)
{
  // This function may be called from Loop.
  // returns  1 if entry is accepted.
  // returns -1 otherwise.
  return 1;
}
#endif // #ifdef evetree_cxx
