
#include <iostream>
#include <functional>
#include <iomanip>
#include <cmath>
#include <fstream>
#include <new>
#include <climits>
#include <vector>


#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
#include "TH1F.h"
#include "TH2F.h"
#include "evetree.h"          //Contains root file format
#include "Constants.h"

#include "StraightLineFit.h"

using namespace std;

//#define DATA_ANALYSIS

Double_t cal_slope2(Double_t x, Double_t* par) {
  if (x<nstrip/2. -0.5) {
    return par[0] + par[1]*(x - nstrip/4. +0.5);
  } else {
    double par3 = (par[2]-par[0]-par[1]*nstrip/4.)/(nstrip/4.);
    return par[2] + par3*(x - 3*nstrip/4. +0.5);
  }
}

//After second fitting if we called GetPosInStrip two times alignment correction will be done. Is this correct?
void GetPosInStrip(int ixy, double* ext, double* otherext, double* off, double* pos, double* local) {
  for (int ij=0; ij<nlayer; ij++) {
    local[ij] = ext[ij]; //+off[ij];
    #if defined(DATA_ANALYSIS)
    if (ixy==0) {
      local[ij] +=cal_slope2(otherext[ij], &align_ystr_xdev[ij][0]);
      local[ij] +=cal_slope2(local[ij], &align_xstr_xdev[ij][0]);
    } else {
      local[ij] +=cal_slope2(otherext[ij], &align_xstr_ydev[ij][0]);
      local[ij] +=cal_slope2(local[ij], &align_ystr_ydev[ij][0]);
    }
    local[ij] +=off[ij];
    #endif
    int istr = int(local[ij]);
    if (local[ij]<0.0) {
      pos[ij] = local[ij] - istr + 0.5;
    } else {
      pos[ij] = local[ij] - istr - 0.5;
    }

    local[ij] -= 0.5;
  }
}


int main() {



  //TFile *fileIn_Mag = new TFile("gen4_inout1_SingleStack_with_noise.root","read");
#if !defined(DATA_ANALYSIS)
  TFile *fileIn_NonMag = new TFile("gen4_inout1_SingleStack_with_noise_30_percent.root","read");
#endif
#if defined(DATA_ANALYSIS)
  TFile *fileIn_NonMag = new TFile("RPCv4t_evtraw-20181226_192924.rre","read");
#endif

  //////////////////////////////////////////////
  //                                          //
  //  Saving Positions and time if the time   //
  //   of hit is between 310 ns to 360 ns     //
  //                                          //
  //////////////////////////////////////////////

  double xFirstPosition[nlayer];
  double yFirstPosition[nlayer];
  double xFirstTime[nlayer];
  double yFirstTime[nlayer];
  double xFirstPID[nlayer];
  double yFirstPID[nlayer];

  //////////////////////////////////////////////
  //                                          //
  //           Position Analysis              //
  //                                          //
  //////////////////////////////////////////////

  //X-side
  int    xhits[nlayer],xfitfailed,xndof;
  double xpos[nlayer],xxerr[nlayer],xintersect,xslope,xerrinter,xerrslope,xerrcovar,xchisquare;
  double xext[nlayer],xdev[nlayer],xexter[nlayer],xextloc[nlayer],xposinstr[nlayer];
  bool   xusedpos[nlayer];

  //Y-side
  int    yhits[nlayer],yfitfailed,yndof;
  double ypos[nlayer],yyerr[nlayer],yintersect,yslope,yerrinter,yerrslope,yerrcovar,ychisquare;
  double yext[nlayer],ydev[nlayer],yexter[nlayer],yextloc[nlayer],yposinstr[nlayer];
  bool   yusedpos[nlayer];


  double dist[nlayer], xval, yval;
  //X-side
  int    xthits[nlayer],xtfitfailed,xtndof;
  double xtime[nlayer],xterr[nlayer],xtintersect,xtslope,xterrinter,xterrslope,xterrcovar,xtchisquare;
  double xtext[nlayer],xtdev[nlayer],xtexter[nlayer],xtextloc[nlayer],xtposinstr[nlayer];
  bool   xusedtime[nlayer];

  //Y-side
  int    ythits[nlayer],ytfitfailed,ytndof;
  double ytime[nlayer],yterr[nlayer],ytintersect,ytslope,yterrinter,yterrslope,yterrcovar,ytchisquare;
  double ytext[nlayer],ytdev[nlayer],ytexter[nlayer],ytextloc[nlayer],ytposinstr[nlayer];
  bool   yusedtime[nlayer];


  //Additional variables
  std::vector<int>    xpts[nlayer];                   //Where is each hit
  std::vector<int>    ypts[nlayer];
  std::vector<int>    xyptsfull[nlayer];              //For filling correlated hits

  std::vector<double> xptsall[nlayer];                //Stores all clusters
  std::vector<double> yptsall[nlayer];

  //Variables for first hits
  std::vector<int>    xPoints_all[nlayer];      //Where is each hit (first hit)
  std::vector<int>    yPoints_all[nlayer];
  std::vector<double> xTime_all[nlayer];        //Corrected timing of all clusters
  std::vector<double> yTime_all[nlayer];
  std::vector<int>    xPIDall[nlayer];
  std::vector<int>    yPIDall[nlayer];

  std::vector<double> xPoints[nlayer];          //Stores all clusters---combine hits if they are near
  std::vector<double> yPoints[nlayer];
  std::vector<double> xTime[nlayer];            //Stores all clusters---if they are near find minimum
  std::vector<double> yTime[nlayer];
  std::vector<int>    xPID[nlayer];
  std::vector<int>    yPID[nlayer];


  TH1F* xlayer_occu[nlayer][nset];
  TH1F* ylayer_occu[nlayer][nset];
  TH1F* xPID_layer[nlayer][nset];
  TH1F* yPID_layer[nlayer][nset];
  TH1F* xlayer_mult[nlayer][nset];
  TH1F* ylayer_mult[nlayer][nset];

  TH1F* xlayer_mult_noise[nlayer][nset];
  TH1F* ylayer_mult_noise[nlayer][nset];

  TH1F* x_total_hits[nset];
  TH1F* y_total_hits[nset];
  TH1F* orighits[nset];
  TH2F* rawhits_corr_xymul[nlayer][nset];
  TH2F* raw_occu[nlayer][nset];
  TH2F* rawhits_xlay_corr_mul[nlayer][nlayer][nset];
  TH2F* rawhits_ylay_corr_mul[nlayer][nlayer][nset];
  TH2F* rawhits_x_vs_y_corr_mul[nlayer][nset];
  TH2F* raw_xmul_vs_time[nlayer][nset];
  TH2F* raw_ymul_vs_time[nlayer][nset];
  TH2F* xtotalhits_vs_layer_corr[nset];
  TH2F* ytotalhits_vs_layer_corr[nset];

  //Noise information are saved in these histograms
  TH1F* n_multall[nset];
  TH1F* layer_hitsall[nlayer][nset];
  TH1F* layer_multall[nlayer][nset];
  TH1F* layer_timeall[nlayer][nset];
  TH2F* total_vs_indmulall[nlayer][nset];
  TH2F* mult_correlall[nlayer][nset];
  TH2F* hist_correlall[nlayer][nset];
  //After removal of muon hits
  TH1F* n_mult[nset];
  TH1F* layer_hits[nlayer][nset];
  TH1F* layer_mult[nlayer][nset];
  TH1F* layer_time[nlayer][nset];
  TH2F* total_vs_indmul[nlayer][nset];
  TH2F* mult_correl[nlayer][nset];
  TH2F* hist_correl[nlayer][nset];

  Float_t threco, phreco;
  Int_t nXStripsRaw[nlayer],nYStripsRaw[nlayer];



  //TFile* outputFile = new TFile("output.root","recreate");
  TFile *fileOut = new TFile("multiplicity.root","RECREATE");
  TTree* T1 = new TTree("T1","T1");
  T1->Branch("xndof",&xndof,"xndof/i");
  T1->Branch("yndof",&yndof,"yndof/i");
  T1->Branch("threco",&threco,"threco/F");
  T1->Branch("phreco",&phreco,"phreco/F");
  T1->Branch("xfitfailed",&xfitfailed,"xfitfailed/F");
  T1->Branch("yfitfailed",&yfitfailed,"yfitfailed/F");
  T1->Branch("nXStripsRaw",nXStripsRaw,"nXStripsRaw[12]/I");
  T1->Branch("nYStripsRaw",nYStripsRaw,"nYStripsRaw[12]/I");
  T1->Branch("xchisquare",&xchisquare,"xchisquare/F");
  T1->Branch("ychisquare",&ychisquare,"ychisquare/F");



  char title[200];

  for(int set=0; set<nset; set++) {
    sprintf(title, "n_mult_m%i",set);
    n_mult[set] = new TH1F(title, title, 2*nlayer*nstrip+1, -0.5, 2*nlayer*nstrip+0.5);
    sprintf(title, "n_multall_m%i",set);
    n_multall[set] = new TH1F(title, title, 2*nlayer*nstrip+1, -0.5, 2*nlayer*nstrip+0.5);

    sprintf(title, "orighits_m%i",set);
    orighits[set] = new TH1F(title, title, 2*nstrip, -0.5, 2*nstrip-0.5);
    sprintf(title, "x_total_hits_m%i",set);
    x_total_hits[set] = new TH1F(title, title, 2*nstrip, -0.5, 2*nstrip-0.5);
    sprintf(title, "y_total_hits_m%i",set);
    y_total_hits[set] = new TH1F(title, title, 2*nstrip, -0.5, 2*nstrip-0.5);
    sprintf(title, "xtotalhits_vs_layer_corr_m%i",set);
    xtotalhits_vs_layer_corr[set] = new TH2F(title, title, 2*nstrip, -0.5, 2*nstrip-0.5, nlayer-2, -0.5, nlayer-2.5);
    sprintf(title, "ytotalhits_vs_layer_corr_m%i",set);
    ytotalhits_vs_layer_corr[set] = new TH2F(title, title, 2*nstrip, -0.5, 2*nstrip-0.5, nlayer-2, -0.5, nlayer-2.5);
    for (int ij=0; ij<nlayer; ij++) {
      sprintf(title, "xlayer_occu_l%i_m%i", ij,set);
      xlayer_occu[ij][set] = new TH1F(title, title, nstrip, -0.5, nstrip-0.5);
      //xlayer_occu[ij]->SetDirectory(fFile);
      sprintf(title, "ylayer_occu_l%i_m%i", ij,set);
      ylayer_occu[ij][set] = new TH1F(title, title, nstrip, -0.5, nstrip-0.5);
      sprintf(title, "xPID_layer_l%i_m%i", ij,set);
      xPID_layer[ij][set] = new TH1F(title, title, 2000, -1000, 1000);
      sprintf(title, "yPID_layer_l%i_m%i", ij,set);
      yPID_layer[ij][set] = new TH1F(title, title, 2000, -1000, 1000);

      sprintf(title, "xlayer_mult_l%i_m%i", ij,set);
      xlayer_mult[ij][set] = new TH1F(title, title, nstrip, -0.5, nstrip-0.5);
      sprintf(title, "ylayer_mult_l%i_m%i", ij,set);
      ylayer_mult[ij][set] = new TH1F(title, title, nstrip, -0.5, nstrip-0.5);

      sprintf(title, "xlayer_mult_noise_l%i_m%i", ij,set);
      xlayer_mult_noise[ij][set] = new TH1F(title, title, nstrip, -0.5, nstrip-0.5);
      sprintf(title, "ylayer_mult_noise_l%i_m%i", ij,set);
      ylayer_mult_noise[ij][set] = new TH1F(title, title, nstrip, -0.5, nstrip-0.5);

      sprintf(title, "rawhits_corr_xymul_l%i_m%i", ij,set);
      rawhits_corr_xymul[ij][set] = new TH2F(title, title, 8*(nstrip+1), -0.5, 8*(nstrip+0.5),  8*(nstrip+1), -0.5, 8*(nstrip+0.5));
      sprintf(title, "raw_occu_l%i_m%i", ij,set);
      raw_occu[ij][set] = new TH2F(title, title, nstrip, -0.5, nstrip-0.5, nstrip, -0.5, nstrip-0.5);
      for (int jk=ij+1; jk<nlayer; jk++) {
        sprintf(title, "rawhits_xlay_corr_mul_l%i_l%i_m%i", ij, jk, set);
        rawhits_xlay_corr_mul[ij][jk][set] = new TH2F(title, title, 8*(nstrip+1), -0.5, 8*(nstrip+0.5),  8*(nstrip+1), -0.5, 8*(nstrip+0.5));
        sprintf(title, "rawhits_ylay_corr_mul_l%i_l%i_m%i", ij, jk, set);
        rawhits_ylay_corr_mul[ij][jk][set] = new TH2F(title, title, 8*(nstrip+1), -0.5, 8*(nstrip+0.5),  8*(nstrip+1), -0.5, 8*(nstrip+0.5));
      }
      sprintf(title, "rawhits_x_vs_y_corr_mul_l%i_m%i", ij, set);
      rawhits_x_vs_y_corr_mul[ij][set] = new TH2F(title, title, 8*(nstrip+1), -0.5, 8*(nstrip+0.5),8*(nstrip+1), -0.5, 8*(nstrip+0.5));
      sprintf(title, "raw_xmul_vs_time_l%i_m%i", ij, set);
      raw_xmul_vs_time[ij][set] = new TH2F(title, title, 8*(nstrip+1), -0.5, 8*(nstrip+0.5),1000, -0.5, 1000);
      sprintf(title, "raw_ymul_vs_time_l%i_m%i", ij, set);
      raw_ymul_vs_time[ij][set] = new TH2F(title, title, 8*(nstrip+1), -0.5, 8*(nstrip+0.5),1000, -0.5, 1000);

      sprintf(title, "layer_hitsall_l%i_m%i", ij, set);
      layer_hitsall[ij][set] = new TH1F(title, title, 2*nstrip, -0.5, 2*nstrip-0.5);
      sprintf(title, "layer_multall_l%i_m%i", ij, set);
      layer_multall[ij][set] = new TH1F(title, title, 2*nstrip+1, -0.5, 2*nstrip+0.5);
      sprintf(title, "layer_timeall_l%i_m%i", ij, set);
      layer_timeall[ij][set] = new TH1F(title, title, 1000, -10., 10.);
      sprintf(title, "hist_correlall_l%i_m%i", ij, set);
      hist_correlall[ij][set] = new TH2F(title, title, 2*nstrip, -0.5, 2*nstrip-0.5, 2*nstrip, -0.5, 2*nstrip-0.5);
      sprintf(title, "total_vs_indmulall_l%i_m%i", ij, set);
      total_vs_indmulall[ij][set] = new TH2F(title, title, 2*nlayer*nstrip+1, -0.5, 2*nlayer*nstrip+0.5, 2*nstrip+1, -0.5, 2*nstrip+0.5);
      sprintf(title, "mult_correlall_l%i_m%i", ij, set);
      mult_correlall[ij][set] = new TH2F(title, title, 2*nstrip+1, -0.5, 2*nstrip+0.5, 2*nstrip+1, -0.5, 2*nstrip+0.5);

      sprintf(title, "layer_hits_l%i_m%i", ij, set);
      layer_hits[ij][set] = new TH1F(title, title, 2*nstrip, -0.5, 2*nstrip-0.5);
      sprintf(title, "layer_mult_l%i_m%i", ij, set);
      layer_mult[ij][set] = new TH1F(title, title, 2*nstrip+1, -0.5, 2*nstrip+0.5);
      sprintf(title, "layer_time_l%i_m%i", ij, set);
      layer_time[ij][set] = new TH1F(title, title, 1000, -10., 10.);
      sprintf(title, "hist_correl_l%i_m%i", ij, set);
      hist_correl[ij][set] = new TH2F(title, title, 2*nstrip, -0.5, 2*nstrip-0.5, 2*nstrip, -0.5, 2*nstrip-0.5);
      sprintf(title, "total_vs_indmul_l%i_m%i", ij, set);
      total_vs_indmul[ij][set] = new TH2F(title, title, 2*nlayer*nstrip+1, -0.5, 2*nlayer*nstrip+0.5, 2*nstrip+1, -0.5, 2*nstrip+0.5);
      sprintf(title, "mult_correl_l%i_m%i", ij, set);
      mult_correl[ij][set] = new TH2F(title, title, 2*nstrip+1, -0.5, 2*nstrip+0.5, 2*nstrip+1, -0.5, 2*nstrip+0.5);


    }
  }

  for(int set=0; set<1; set++) {

    TTree *event_tree;
#if !defined(DATA_ANALYSIS)
    if(set==0){event_tree=(TTree*)fileIn_NonMag->Get("T2");}
#endif
#if defined(DATA_ANALYSIS)
    if(set==0){event_tree=(TTree*)fileIn_NonMag->Get("evetree");}
#endif
  //  if(set==1){event_tree=(TTree*)fileIn_Mag->Get("T2");}
    evetree *event = new evetree(event_tree);

    //__________________________________________________________________________
    event->Loop();
    Int_t nentry = event_tree->GetEntries();
    int nTrigentry =0;
    nentry = 3000000;
    for(int iev=0;iev<nentry;iev++) {
      //cout<<iev<<endl;
      event_tree->GetEntry(iev);

      #if !defined(DATA_ANALYSIS)
         // Set object pointer
         //  const char *sideMark[2] = {"x","y"};
         for(int i=0;i<nlayer;i++){
           for(int j=0;j<nchannel;j++){
             event->vxtdc_l[i][j]=0;
             event->vytdc_l[i][j]=0;
             event->vxtdc_t[i][j]=0;
             event->vytdc_t[i][j]=0;
           }
         }
      #endif

      //Initialize vectors
      for(int jk=0;jk<nlayer;jk++) {
        xpts[jk].clear(); ypts[jk].clear();
        xptsall[jk].clear(); yptsall[jk].clear();
        xPoints_all[jk].clear(); yPoints_all[jk].clear();
        xTime_all[jk].clear(); yTime_all[jk].clear();
        xPoints[jk].clear(); yPoints[jk].clear();
        xTime[jk].clear(); yTime[jk].clear();
        xPIDall[jk].clear(); yPIDall[jk].clear();
        xyptsfull[jk].clear();
      }



      /////////////////////////////////////
      //                                 //
      //  Geant Digitization to data     //
      //                                 //
      /////////////////////////////////////


   #if !defined(DATA_ANALYSIS)
      //if(event->ndigiht>500) {continue;}// hard-coded max-particle in T2.h
      for(int ij=0; ij<nlayer; ij++) {
        event->ENum[ij] = event->ievt;
        event->tdc_ref_l[ij] = 0;
        event->xLayer[ij] = new TBits(0);
        event->yLayer[ij] = new TBits(0);
        for(int jk=0; jk<nchannel; jk++) {
          event->vxtdc_l[ij][jk] = new std::vector<unsigned int>(0);
          event->vytdc_l[ij][jk] = new std::vector<unsigned int>(0);
        }
      }

      // if(ievt==481249) {
      //  TString msg1 = TString::Format("Ngenerated:%i\n", ngenerated);
      //  if (gProofServ) gProofServ->SendAsynMessage(msg1);
      // }
      //
      int TrgLayer[4];
      int nTriggerX[nlayer];
      int nTriggerY[nlayer];

      TrgLayer[0] = 6;
      TrgLayer[1] = 7;
      TrgLayer[2] = 8;
      TrgLayer[3] = 9;
      int triglays = 4;
      for(int ij=0; ij<nlayer; ij++) {
        nTriggerX[ij] = 0; nTriggerY[ij]=0;
      }


      //cout<<"\nNew Event, Ndigihit="<<event->ndigiht<<endl;
      if(event->ndigiht==0){
        for(int ij=0; ij<nlayer; ij++) {
          if(event->xLayer[ij]){delete event->xLayer[ij];}
          if(event->yLayer[ij]){delete event->yLayer[ij];}
          for(int jk=0; jk<nchannel; jk++) {
            if(event->vxtdc_l[ij][jk]){delete event->vxtdc_l[ij][jk];}
            if(event->vytdc_l[ij][jk]){delete event->vytdc_l[ij][jk];}
          }
        }
        continue;
      }

      for(unsigned int ki=0;ki<event->ndigiht;ki++) {

        UInt_t strpID = event->stripid[ki];

        int nInSide = (strpID>>31)==0?0:1; // x or y

        strpID>>=8;
        int nInX = strpID%128;
        strpID>>=7;
        int nInCH = strpID%8;
        strpID>>=3;
        int nInMO = strpID%8;
        strpID>>=3;
        int nInLA = strpID%256;
        strpID>>=8;
        int nInDT = strpID%4;

        if(nInSide == 0 && nInX<64 && nInX>=0 && nInLA<10 && nInLA>=0) {
          event->xLayer[nInLA]->SetBitNumber(nInX);
          event->vxtdc_l[nInLA][nInX%8]->push_back(event->digitime[ki]-3010.0);
          //timex_reso_Chk[nInLA]->Fill(600.0+0.1*(event->digitime[ki]-3010.0));
          //vxtdc_l[nInLA][nInX%8]->push_back(digitime[ki]*2.0-2950.0);
          //timex_reso_Chk[nInLA]->Fill(600.0+0.1*(digitime[ki]*2.0-2950.0));
          //TString msg1 = TString::Format("digitime=%d\n", digitime[ki]);
          //if (gProofServ) gProofServ->SendAsynMessage(msg1);
          if(TrgLayer[0]==nInLA || TrgLayer[1]==nInLA || TrgLayer[2]==nInLA || TrgLayer[3]==nInLA) {
            nTriggerX[nInLA]++;
            //cout<<"Trigger in X L"<<nInLA<<"\tstrip"<<nInX<<endl;
          }
        }

        if(nInSide == 1 && nInX<64 && nInX>=0 && nInLA<10 && nInLA>=0) {
          event->yLayer[nInLA]->SetBitNumber(nInX);
          event->vytdc_l[nInLA][nInX%8]->push_back(event->digitime[ki]-3010.0);
          //timey_reso_Chk[nInLA]->Fill(600.0+0.1*(digitime[ki]-3010.0));
          //vytdc_l[nInLA][nInX%8]->push_back(digitime[ki]*2.0-2950.0);
          //timey_reso_Chk[nInLA]->Fill(600.0+0.1*(digitime[ki]*2.0-2950.0));
          if(TrgLayer[0]==nInLA || TrgLayer[1]==nInLA || TrgLayer[2]==nInLA || TrgLayer[3]==nInLA) {
            nTriggerY[nInLA]++;
            //cout<<"Trigger in Y L"<<nInLA<<"\tstrip"<<nInX<<endl;
          }
        }

      }  // digihit ending

      int sw_trigx = 0;
      int sw_trigy = 0;

      for(int tr1=0; tr1<nlayer; tr1++) {
        for(int tr2=0; tr2<triglays; tr2++) {
          if(tr1==TrgLayer[tr2]) {
    	       if(nTriggerX[tr1]>0) {sw_trigx++;}
    	       if(nTriggerY[tr1]>0) {sw_trigy++;}
          }
        }
      }

      if(sw_trigx<=3 && sw_trigy<=3) {
        for(int ij=0; ij<nlayer; ij++) {
          if(event->xLayer[ij]){delete event->xLayer[ij];}
          if(event->yLayer[ij]){delete event->yLayer[ij];}
          for(int jk=0; jk<nchannel; jk++) {
            if(event->vxtdc_l[ij][jk]){delete event->vxtdc_l[ij][jk];}
            if(event->vytdc_l[ij][jk]){delete event->vytdc_l[ij][jk];}
          }
        }
        continue;
      }

   #endif

      nTrigentry++;

      //Fill position hits in xptsall vector. Also corrosponding timing
      //  in xTime_all and yTime_all vectors provided the time is between 310 to 360
      //cout<<"\n"<<iev;

      int xxnhits=0;
      int yynhits=0;
      int Orighits=0;

      for(int jk=0;jk<nlayer-2;jk++) {
        nXStripsRaw[jk] = 100;
        nYStripsRaw[jk] = 100;
      }
      for(int jk=0;jk<nlayer-2;jk++) {
        //cout<<"\nxLayer"<<jk<<"\t";
        for(int kl=0; kl<nstrip; kl++) {
          //X-Side
          if(event->xLayer[jk]->TestBitNumber(kl)) {
            //cout<<"\t"<<kl;
            xpts[jk].push_back(kl);
            xyptsfull[jk].push_back(kl);
            //cout <<kl <<"\t";
            //if(vxtdc_l[jk][kl%8]->size()>0 && pdgidx[jk][kl]->size()>0) { //DATA_ANALYSIS
            //for(int ix=0; ix<event->vxtdc_l[jk][kl%8]->size();ix++) {
            if(event->vxtdc_l[jk][kl%8]->size()){
              double tmpxtime = 600+0.1*int(event->vxtdc_l[jk][kl%8]->at(0)-event->tdc_ref_l[jk]);
              //double tmpxtime = 600+0.1*int(event->vxtdc_l[jk][kl%8]->at(ix)-event->tdc_ref_l[jk]);
              //TString msg1 = TString::Format("%f\n", tmpxtime);
              //if (gProofServ) gProofServ->SendAsynMessage(msg1);
              //if(tmpxtime > 310. && tmpxtime < 350){
              xPoints_all[jk].push_back(kl);
              xTime_all[jk].push_back(tmpxtime);
              //xPIDall[jk].push_back(pdgidx[jk][kl]->at(0));
              //}
            }
          }
        }
        //cout<<"\nyLayer"<<jk<<"\t";
        for(int kl=0; kl<nstrip; kl++) {
          //Y-Side
          if(event->yLayer[jk]->TestBitNumber(kl)) {
            //cout<<"\t"<<kl;
            ypts[jk].push_back(kl);
            xyptsfull[jk].push_back(numberInX+kl);
            //if(vytdc_l[jk][kl%8]->size()>0 && pdgidy[jk][kl]->size()>0) {
            //for(int iy=0;iy<event->vytdc_l[jk][kl%8]->size();iy++) {
            if(event->vytdc_l[jk][kl%8]->size()){
              double tmpytime = 600.+ 0.1*int(event->vytdc_l[jk][kl%8]->at(0)-event->tdc_ref_l[jk]);
              //double tmpytime = 600.+ 0.1*int(event->vytdc_l[jk][kl%8]->at(iy)-event->tdc_ref_l[jk]);
              //TString msg1 = TString::Format("%f\n", tmpytime);
              //if (gProofServ) gProofServ->SendAsynMessage(msg1);
              //if(tmpytime > 310. && tmpytime < 350){
              yPoints_all[jk].push_back(kl);
              yTime_all[jk].push_back(tmpytime);
              //yPIDall[jk].push_back(pdgidy[jk][kl]->at(0));
              //}
            }
          }
        }

        //cout<<endl;




        for (unsigned int ix=0; ix<xpts[jk].size(); ix++) {
          xlayer_occu[jk][set]->Fill(xpts[jk][ix]);
          //xPID_layer[jk]->Fill(xPIDall[jk][ix]);
        }
        for (unsigned int iy=0; iy<ypts[jk].size(); iy++) {
          ylayer_occu[jk][set]->Fill(ypts[jk][iy]);
          //yPID_layer[jk]->Fill(yPIDall[jk][iy]);
        }
        nXStripsRaw[jk] = xpts[jk].size();
        nYStripsRaw[jk] = ypts[jk].size();
        xlayer_mult[jk][set]->Fill(xpts[jk].size());
        ylayer_mult[jk][set]->Fill(ypts[jk].size());
        rawhits_corr_xymul[jk][set]->Fill(xpts[jk].size(), ypts[jk].size());

        for (unsigned int ix=0; ix<xpts[jk].size(); ix++) {
          for (unsigned int iy=0; iy<ypts[jk].size(); iy++) {
            raw_occu[jk][set]->Fill(xpts[jk][ix], ypts[jk][iy]);
          }
        }

        rawhits_x_vs_y_corr_mul[jk][set]->Fill(xPoints_all[jk].size(), yPoints_all[jk].size());

        for (unsigned int ix=0; ix<xTime_all[jk].size(); ix++) {
          raw_xmul_vs_time[jk][set]->Fill(xPoints_all[jk].size(),xTime_all[jk][ix]);
        }

        for (unsigned int ix=0; ix<yTime_all[jk].size(); ix++) {
          raw_ymul_vs_time[jk][set]->Fill(yPoints_all[jk].size(),yTime_all[jk][ix]);
        }

        //cout<<"-------------------------------------------------"<<endl;

        xxnhits+=xPoints_all[jk].size();
        yynhits+=yPoints_all[jk].size();
        Orighits+=max(xxnhits,yynhits);

      }



      for(int jk=0; jk<nlayer; jk++) {
        for (int ij=jk+1; ij<nlayer; ij++) {
          rawhits_xlay_corr_mul[jk][ij][set]->Fill(xPoints_all[jk].size(), xPoints_all[ij].size());
          rawhits_ylay_corr_mul[jk][ij][set]->Fill(yPoints_all[jk].size(), yPoints_all[ij].size());
        }
      }

      //cout<<"===================================================="<<endl;

      x_total_hits[set]->Fill(xxnhits);
      y_total_hits[set]->Fill(yynhits);
      orighits[set]->Fill(Orighits);
      for(int jk=0; jk<nlayer-2; jk++) {
        xtotalhits_vs_layer_corr[set]->Fill(xxnhits, jk, xPoints_all[jk].size());
        ytotalhits_vs_layer_corr[set]->Fill(yynhits, jk, yPoints_all[jk].size());
      }
      for(int jk=0; jk<nlayer-2; jk++) {
        xhits[jk] = xpts[jk].size();
        yhits[jk] = ypts[jk].size();
        xdev[jk] = 100; xpos[jk]=0.0;
        ydev[jk] = 100; ypos[jk]=0.0;
        //X-Side
        if (xhits[jk]<=0 || xhits[jk] > nmxhits) {
          xpos[jk]= -100;
        } else {
          for (int ix=0; ix<xhits[jk]; ix++) {
            // Only layers with one hit or one cluster is used
            if (ix<xhits[jk]-1 && abs(xpts[jk][ix]-xpts[jk][ix+1])>1)
            {xpos[jk]=-100; break;}
            xpos[jk] += xpts[jk][ix];
          }
          int tempx=0;
          int mul=1;
          for(int ix=0; ix<xhits[jk]; ix++) {
            //Find clusters and save to xptsall[jk]
            if(mul==1){tempx += xpts[jk][ix];}
            if(ix<xhits[jk]-1 &&  abs(xpts[jk][ix]-xpts[jk][ix+1])==1)
            {tempx += xpts[jk][ix+1]; mul++;}
            else {
              double val = (double)tempx/(double)mul + 0.5 - xoff[jk];
              #if defined(DATA_ANALYSIS)
              val -= cal_slope2(val, &align_xstr_xdev[jk][0]);
              #endif
              xptsall[jk].push_back(val);
              mul=1;
              tempx=0;
            }
          }
        }
        xxerr[jk] = errxco[jk]*errxco[jk];
        if (xpos[jk]>=0.0) {
          xpos[jk]  = xpos[jk]/xhits[jk] + 0.5 - xoff[jk];
          #if defined(DATA_ANALYSIS)
          //Aligmnent Correction
          xpos[jk] -= cal_slope2(xpos[jk], &align_xstr_xdev[jk][0]);
          #endif
          xxerr[jk] = xposerrsq[xhits[jk]-1][jk];
        }
        //Y-Side
        if (yhits[jk]<=0 || yhits[jk] > nmxhits) {
            ypos[jk]= -100;
        } else {
          for (int ix=0; ix<yhits[jk]; ix++) {
              //Only layers with one hit or one cluster is used.
              if (ix<yhits[jk]-1 && abs(ypts[jk][ix]-ypts[jk][ix+1])>1)
              {ypos[jk]=-100; break;}
              // cout<<"No of hits="<<yhits[jk]<<endl;
              // cout<<"jk="<<jk<<"\tix="<<ix<<endl;
              // cout<<"ypts="<<ypts[jk][ix]<<endl;
              ypos[jk] += ypts[jk][ix];
          }
          int tempy=0;
          int mul=1;
          for(int ix=0; ix<yhits[jk]; ix++) {
            //Find clusters and save to yptsall[jk]
            if(mul==1){tempy += ypts[jk][ix];}
            if(ix<yhits[jk]-1 &&  abs(ypts[jk][ix]-ypts[jk][ix+1])==1)
             {tempy += ypts[jk][ix+1]; mul++;}
            else {
              double val = (double)tempy/(double)mul + 0.5 - yoff[jk];
              #if defined(DATA_ANALYSIS)
              val -= cal_slope2(val, &align_ystr_ydev[jk][0]);
              #endif
              yptsall[jk].push_back(val);
              mul=1;
              tempy=0;
            }
          }
        }
        yyerr[jk] = errxco[jk]*errxco[jk];
        if (ypos[jk]>=0.0) {
          ypos[jk]  = ypos[jk]/yhits[jk] + 0.5 - yoff[jk];
          //Aligmnent Correction
          #if defined(DATA_ANALYSIS)
          ypos[jk] -= cal_slope2(ypos[jk], &align_ystr_ydev[jk][0]);
          #endif
          yyerr[jk] = yposerrsq[yhits[jk]-1][jk];
        }
        //Sort out hits, which can be used for fit
        xusedpos[jk] = (xpos[jk]>-99 && xhits[jk]<=nmxhits) ? true : false;
        yusedpos[jk] = (ypos[jk]>-99 && yhits[jk]<=nmxhits) ? true : false;
      }//NLAYER Loop

      //Fitting
      //Combined Position Fit

      //Temporary storing of position is required. Otherwise after each fitting alignment correction is done
      double tempxpos[nlayer], tempypos[nlayer];
      copy(begin(xpos), end(xpos), begin(tempxpos));
      copy(begin(ypos), end(ypos), begin(tempypos));

      double zvalue[nlayer];
      for (unsigned int ix=0; ix<nlayer; ix++) {
        zvalue[ix]=layerzpos[ix];
        xext[ix]= xextloc[ix] = xexter[ix] =xposinstr[ix] =  100000;
        yext[ix]= yextloc[ix] = yexter[ix] =yposinstr[ix] =  100000;
      }

      //X-Side fit
      StraightLineFit xposfit(1, zvalue, xpos,  xxerr, xusedpos, 10, 11, layfirst, laylast, xyPosDev);
      xposfit.GetParameters(xfitfailed, xintersect, xslope);
      xposfit.GetFitValues(xext, xdev, xexter);
      //Y-Side fit
      StraightLineFit yposfit(1, zvalue, ypos,  yyerr, yusedpos, 10, 11, layfirst, laylast, xyPosDev);
      yposfit.GetParameters(yfitfailed, yintersect, yslope);
      yposfit.GetFitValues(yext, ydev, yexter);

      #if defined(DATA_ANALYSIS)
      for (int ix=0; ix<nlayer; ix++) {
        xpos[ix] -=cal_slope2(yext[ix], &align_ystr_xdev[ix][0]);
        ypos[ix] -=cal_slope2(xext[ix], &align_xstr_ydev[ix][0]);
      }
      #endif

      //X-Side fit
      xposfit = StraightLineFit(1, zvalue, xpos,  xxerr, xusedpos, 10, 11, layfirst, laylast, xyPosDev);
      xposfit.GetParameters(xfitfailed, xintersect, xslope);
      xposfit.GetFitValues(xext, xdev, xexter);
      xposfit.GetChisqure(xndof,xchisquare);
      GetPosInStrip(0, xext, yext, xoff, xposinstr, xextloc);
      //Y-Side fit
      yposfit = StraightLineFit(1, zvalue, ypos,  yyerr, yusedpos, 10, 11, layfirst, laylast, xyPosDev);
      yposfit.GetParameters(yfitfailed, yintersect, yslope);
      yposfit.GetFitValues(yext, ydev, yexter);
      yposfit.GetChisqure(yndof,ychisquare);
      GetPosInStrip(1, yext, xext, yoff, yposinstr, yextloc);


      copy(begin(tempxpos), end(tempxpos), begin(xpos));
      copy(begin(tempypos), end(tempypos), begin(ypos));

      double stripwidth = 3.0; //cm
      const double pival=acos(-1);

      threco = acos(sqrt(1./(1+pow(stripwidth*xslope,2.)+pow(stripwidth*yslope,2.)))); //10Nov//S= sqrt(dx^2+dy^2+dz^2)) theta= acos(height / S) height = dz;
      threco = (180./pival)*threco; //acos(sqrt(1./(1+pow(stripwidth*xslope,2.)+pow(stripwidth*yslope,2.))));
      phreco = atan2(yslope, xslope);  // What is the direction



      for(int jk=0;jk<nlayer; jk++) {
        xlayer_mult[jk][set]->Fill(xpts[jk].size());
        ylayer_mult[jk][set]->Fill(ypts[jk].size());
      }


    //  if(set==0 && phreco>0.6 && phreco<0.65 &&xndof==10) {cout<<phreco<<endl;}


      //cout<<"0-"<<xext[0]<<"\t1-"<<xext[1]<<"\t2-"<<xext[2]<<"\t3-"<<xext[3]<<"\t4-"<<xext[4]<<"\t5-"<<xext[5]<<"\t6-"<<xext[6]<<"\t7-"<<xext[7]<<"\t8-"<<xext[8]<<"\t9-"<<xext[9]<<endl;

      if(set==0) {T1->Fill();}


      int init=-1;
      double initzpos=0;
      //Timing Correction   ----- We need ortoganal position for correct for timing delay
      for(int jk=0; jk<nlayer-2; jk++) {
          double tshft=1000.0;
          xtdev[jk] = -10000.0; xtime[jk]=-1000.0;
          ytdev[jk] = -10000.0; ytime[jk]=-1000.0;
          //X-Side
          double inittime = 11000;
          for (int xix=0; xix<int(xpts[jk].size()); xix++) {
            int chNum = xpts[jk][xix]%8;
            if(event->vxtdc_l[jk][chNum]->size()) {
               double tmpxtime = 600+0.1*int(event->vxtdc_l[jk][chNum]->at(0)-event->tdc_ref_l[jk]);
               if(tmpxtime<inittime) {
                 inittime = tmpxtime; tshft = xtoffset[jk][xpts[jk][xix]];
               }
            }
          }
          if(inittime<10000) {
            xtime[jk] = inittime;
            xtime[jk] -=timeoffsetx[jk];
            xtime[jk] -=tshft;
            xtime[jk] -=slope_path*yext[jk];
          }
          //Y-Side
          inittime = 11000;
          for (int yiy=0; yiy<int(ypts[jk].size()); yiy++) {
            int chNum = ypts[jk][yiy]%8;
            if(event->vytdc_l[jk][chNum]->size()) {
               double tmpytime = 600+0.1*int(event->vytdc_l[jk][chNum]->at(0)-event->tdc_ref_l[jk]);
               if(tmpytime<inittime) {
                 inittime = tmpytime; tshft = ytoffset[jk][ypts[jk][yiy]];
               }
            }
          }
          if(inittime<10000) {
            ytime[jk] = inittime;
            ytime[jk] -=timeoffsety[jk];
            ytime[jk] -=tshft;
            ytime[jk] -=slope_path*xext[jk];
          }
          //Sort out hits, which can be used for fit
          xusedtime[jk] = (xtime[jk]>-99) ? true : false;
          yusedtime[jk] = (ytime[jk]>-99) ? true : false;
          dist[jk] =-100;
          xval=-100; yval=-100;
          if (init<0) {
            xval = xextloc[jk];
            yval = yextloc[jk];
            dist[jk] = 0.0;
            init = jk;
            initzpos = layerzpos[jk];
          } else {
            dist[jk] = sqrt( pow((xextloc[jk] - xval)*stripwidth, 2.) +
            pow((yextloc[jk] - yval)*stripwidth, 2.) +
            pow(layerzpos[jk] - initzpos, 2.));
            //                                     pow((ij - init)*layergap, 2.));
          }

      }

      // // Xtime fit
	    // int iTimeSlopeConst = 1;
	    // StraightLineFit xtimefit(iTimeSlopeConst, dist, xtime,  timeserrx2, xusedtime, 10, 11, layfirst, laylast, float(7.0));
	    // xtimefit.GetParameters(xtfitfailed, xtintersect, xtslope);
	    // xtimefit.GetFitValues(xtext, xtdev, xtexter);
	    // xtimefit.GetChisqure(xtndof, xtchisquare);
      //
      // // Ytime fit
	    // StraightLineFit ytimefit(iTimeSlopeConst, dist, ytime,  timeserry2, yusedtime, 10, 11, layfirst, laylast, float(7.0));
	    // ytimefit.GetParameters(ytfitfailed, ytintersect, ytslope);
	    // ytimefit.GetFitValues(ytext, ytdev, ytexter);
	    // ytimefit.GetChisqure(ytndof, ytchisquare);


      //Before removing noise hits
      int nmulall=0;
      int xhitsall[nlayer];
      for (int ij=0; ij<nlayer; ij++) {
	       xhits[ij] = xpts[ij].size();
	       yhits[ij] = ypts[ij].size();
	       nmulall += xhitsall[ij] = xhits[ij] + yhits[ij]; // = xptsall[ij].size();  + yptsall[ij].size();
      }

      n_multall[set]->Fill(nmulall);

      //cout<<"\nnmulall="<<nmulall<<endl;
      //if(nmulall==0){return 0;}


      for (int ij=0; ij<nlayer; ij++) {
	       for (int jk=0; jk<=xhitsall[ij]; jk++) {
	          total_vs_indmulall[ij][set]->Fill(nmulall, jk); //check this
	       }
      }
      for (int ij=0; ij<nlayer; ij++) {
        for(int jk=0; jk<xyptsfull[ij].size();jk++) {
          layer_hitsall[ij][set]->Fill(xyptsfull[ij][jk]);
          hist_correlall[ij][set]->Fill(xyptsfull[ij][jk], xyptsfull[ij][jk]);
          for(int kl=jk+1; kl<xyptsfull[ij].size(); kl++) {
            hist_correlall[ij][set]->Fill(xyptsfull[ij][jk], xyptsfull[ij][kl]);
          }
        }
        for (int jk=ij+1; jk<nlayer; jk++) {
      	  mult_correlall[jk][set]->Fill(xyptsfull[ij].size(), xyptsfull[jk].size());
      	}
        layer_timeall[ij][set]->Fill(xtime[ij]-xtext[ij]);
      	layer_timeall[ij][set]->Fill(ytime[ij]-ytext[ij]);
      }





      //Remove the hits near the muon.... select only noise
      if (xndof >4 && xchisquare<50) {
	       for (int ij=0; ij<nlayer; ij++) {
	          for (int ix=0; ix<xpts[ij].size(); ix++) {
	             if (abs(xext[ij] - xpts[ij][ix])<2.6) {
	                xpts[ij].erase(xpts[ij].begin()+ix); ix--;
	             }
	          }
	       }
      }

      if (yndof >4 && ychisquare<50) {
	       for (int ij=0; ij<nlayer; ij++) {
	          for (int iy=0; iy<ypts[ij].size(); iy++) {
	             if (abs(yext[ij] - ypts[ij][iy])<2.6) {
	                ypts[ij].erase(ypts[ij].begin()+iy); iy--;
	             }
	          }
	       }
      }



      //After removing noise hits
      //copying all yhits to xpts
      int nmul = 0;
      for (int ij=0; ij<nlayer; ij++) {
      	xhits[ij] = 0;
        xlayer_mult_noise[ij][set]->Fill(xpts[ij].size());
        ylayer_mult_noise[ij][set]->Fill(ypts[ij].size());
      	for (int iy=0; iy<ypts[ij].size(); iy++) {
      	  xpts[ij].push_back(numberInX+ypts[ij][iy]);
      	}
      	nmul += xhits[ij] = xpts[ij].size();
      }

      n_mult[set]->Fill(nmul);
      for (int ij=0; ij<nlayer; ij++) {
      	for (int jk=0; jk<xhits[ij]; jk++) {
          layer_hits[ij][set]->Fill(xpts[ij][jk]);
          hist_correl[ij][set]->Fill(xpts[ij][jk], xpts[ij][jk]);  ///  ?
          for (int kl=jk+1; kl<xhits[ij]; kl++) {
      	    hist_correl[ij][set]->Fill(xpts[ij][jk], xpts[ij][kl]);
      	  }
      	  total_vs_indmul[ij][set]->Fill(nmul, jk);
    	  }
        for (int jk=ij+1; jk<nlayer; jk++) {
      	  mult_correlall[jk][set]->Fill(xhits[ij], xhits[jk]);
      	}
      }

      //Deleting pointers
#if !defined(DATA_ANALYSIS)
        for(int ij=0; ij<nlayer; ij++) {
          if(event->xLayer[ij]){delete event->xLayer[ij];}
          if(event->yLayer[ij]){delete event->yLayer[ij];}
          for(int jk=0; jk<nchannel; jk++) {
            if(event->vxtdc_l[ij][jk]){delete event->vxtdc_l[ij][jk];}
            if(event->vytdc_l[ij][jk]){delete event->vytdc_l[ij][jk];}
          }
        }
#endif



    } //NENTRY LOOP

    // Normalise layer multiplity for each individual value of total multiplicity
    Int_t ntotal = nTrigentry;
    for (int ij=0; ij<nlayer; ij++) {
      //cout<<"ij "<< ij<<endl;
      for (int jk=0; jk<total_vs_indmulall[ij][set]->GetNbinsX(); jk++) {
        double ntt = total_vs_indmulall[ij][set]->GetBinContent(jk+1,1);
        if (ntt<1.0) ntt=1.0;
        for (int kl=0; kl<total_vs_indmulall[ij][set]->GetNbinsY(); kl++) {
          total_vs_indmulall[ij][set]->SetBinContent(jk+1, kl+1, total_vs_indmulall[ij][set]->GetBinContent(jk+1,kl+1)/ntt);
        }
      }

      for (int jk=0; jk<total_vs_indmul[ij][set]->GetNbinsX(); jk++) {
        double ntt = total_vs_indmul[ij][set]->GetBinContent(jk+1,1);
        if (ntt<1.0) ntt=1.0;
        for (int kl=0; kl<total_vs_indmul[ij][set]->GetNbinsY(); kl++) {
          total_vs_indmul[ij][set]->SetBinContent(jk+1, kl+1, total_vs_indmul[ij][set]->GetBinContent(jk+1,kl+1)/ntt);
        }
      }
    }

    double scal = 1./(max(1, ntotal));
    //cout <<"scal "<<scal<<endl;
    n_multall[set]->Scale(scal);
    n_mult[set]->Scale(scal);


    for (int ij=0; ij<nlayer; ij++) {
      layer_hitsall[ij][set]->Scale(scal);
      layer_multall[ij][set]->Scale(scal);
      layer_timeall[ij][set]->Scale(scal);

      mult_correlall[ij][set]->Scale(scal);
      hist_correlall[ij][set]->Scale(scal);

      layer_hits[ij][set]->Scale(scal);
      layer_mult[ij][set]->Scale(scal);
      layer_time[ij][set]->Scale(scal);

      mult_correl[ij][set]->Scale(scal);
      hist_correl[ij][set]->Scale(scal);
    }




  }   //NSET LOOP

  fileOut->cd();

  T1->Write();

  cout<<"Writing Start"<<endl;

  for(unsigned int set=0; set<nset; set++) {
    n_multall[set]->Write(0, TObject::kOverwrite);
    n_mult[set]->Write(0, TObject::kOverwrite);
    x_total_hits[set]->Write(0, TObject::kOverwrite);
    y_total_hits[set]->Write(0, TObject::kOverwrite);
    xtotalhits_vs_layer_corr[set]->Write(0, TObject::kOverwrite);
    ytotalhits_vs_layer_corr[set]->Write(0, TObject::kOverwrite);
    for(unsigned int ij=0; ij<nlayer; ij++) {
      layer_hitsall[ij][set]->Write(0, TObject::kOverwrite);
      layer_hits[ij][set]->Write(0, TObject::kOverwrite);
    }
    for(unsigned int ij=0; ij<nlayer; ij++) {
      layer_multall[ij][set]->Write(0, TObject::kOverwrite);
      layer_mult[ij][set]->Write(0, TObject::kOverwrite);
    }
    for(unsigned int ij=0; ij<nlayer; ij++) {
      layer_timeall[ij][set]->Write(0, TObject::kOverwrite);
      layer_time[ij][set]->Write(0, TObject::kOverwrite);
    }
    for(unsigned int ij=0; ij<nlayer; ij++) {
      total_vs_indmulall[ij][set]->Write(0, TObject::kOverwrite);
      total_vs_indmul[ij][set]->Write(0, TObject::kOverwrite);
    }
    for(unsigned int ij=0; ij<nlayer; ij++) {
      mult_correlall[ij][set]->Write(0, TObject::kOverwrite);
      mult_correl[ij][set]->Write(0, TObject::kOverwrite);
    }
    for(unsigned int ij=0; ij<nlayer; ij++) {
      hist_correlall[ij][set]->Write(0, TObject::kOverwrite);
      hist_correl[ij][set]->Write(0, TObject::kOverwrite);
    }
    for(unsigned int ij=0; ij<nlayer; ij++){
      xlayer_occu[ij][set]->Scale(1./xlayer_occu[ij][set]->GetEntries());
      xlayer_occu[ij][set]->Write(0, TObject::kOverwrite);
      ylayer_occu[ij][set]->Scale(1./ylayer_occu[ij][set]->GetEntries());
      ylayer_occu[ij][set]->Write(0, TObject::kOverwrite);
    }
    for(unsigned int ij=0;ij<nlayer;ij++){
      xlayer_mult[ij][set]->Scale(1./xlayer_mult[ij][set]->GetEntries());
      xlayer_mult[ij][set]->Write(0, TObject::kOverwrite);
      ylayer_mult[ij][set]->Scale(1./ylayer_mult[ij][set]->GetEntries());
      ylayer_mult[ij][set]->Write(0, TObject::kOverwrite);
    }
    for(unsigned int ij=0;ij<nlayer;ij++){
      xlayer_mult_noise[ij][set]->Scale(1./xlayer_mult_noise[ij][set]->GetEntries());
      xlayer_mult_noise[ij][set]->Write(0, TObject::kOverwrite);
      ylayer_mult_noise[ij][set]->Scale(1./ylayer_mult_noise[ij][set]->GetEntries());
      ylayer_mult_noise[ij][set]->Write(0, TObject::kOverwrite);
    }
    for(unsigned int ij=0;ij<nlayer;ij++){
      rawhits_corr_xymul[ij][set]->Write(0, TObject::kOverwrite);
    }
    for(unsigned int ij=0;ij<nlayer;ij++){
      raw_occu[ij][set]->Write(0, TObject::kOverwrite);
    }

    for(unsigned int ij=0;ij<nlayer;ij++){
      rawhits_x_vs_y_corr_mul[ij][set]->Write(0, TObject::kOverwrite);
    }

    for(unsigned int ij=0;ij<nlayer;ij++){
      raw_xmul_vs_time[ij][set]->Write(0, TObject::kOverwrite);
      raw_ymul_vs_time[ij][set]->Write(0, TObject::kOverwrite);
    }

    for(unsigned int ij=0;ij<nlayer;ij++){
      for(unsigned int jk=ij+1;jk<nlayer;jk++){
        rawhits_xlay_corr_mul[ij][jk][set]->Write(0, TObject::kOverwrite);
        rawhits_ylay_corr_mul[ij][jk][set]->Write(0, TObject::kOverwrite);
      }
    }
  }


//  fileIn_Mag->Close();
  fileIn_NonMag->Close();
  fileOut->Close();

//return 0;

}
