#include <string>
#include <vector>
#include <exception>
#include <vector>
#include <utility>  

#include "TROOT.h"
#include "TFile.h"
#include "TDirectory.h"
#include "TChain.h"
#include "TObject.h"
#include "TCanvas.h"
#include "TMath.h"
#include "TLegend.h"
#include "TGraph.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TTree.h"
#include "TF1.h"
#include "TGraphAsymmErrors.h"
#include "TPaveText.h"
#include "TRandom3.h"
#include "TProfile.h"
#include "TDirectory.h"
#include "TCanvas.h"
#include "TProfile.h"
#include "TPaveText.h"
#include "TRandom.h"
#include "TLorentzVector.h"
#include "TH1.h"
#include "TFile.h"

void generate_sample(string OutFilePath="out.txt", string unusedArgument="")
{
     if(OutFilePath == "")return;
     
     FILE* OutputFile = fopen(OutFilePath.c_str(), "w");

     float mu=0.5;
     float sigma=0.1;
     float nb_events=1000000;

     TRandom3 rand(12345);

     TFile* OutputHisto = new TFile((OutFilePath+".root").c_str(),"RECREATE");
     OutputHisto->cd();
     TH2F* HGen         = new TH2F("pdf", "pdf;x;y", 100,0,1, 100, 0, 1 );



     int i=0;
     while(i<nb_events){
        double x = rand.Uniform();
        double y = rand.Uniform();

//       double r=0.5+rand.Uniform()* 0.2; //rand.Gaus(0.2, 0.1);
//       double phi=rand.Uniform() * 2.0 * 3.1415;
//       double x=r*cos(phi);
//       double y=r*sin(phi);

        double eps=0.1;
        double fct=1.0+eps*(2.0*y);
        double max_weight=1.0+2*eps;

//        if(fct>max_weight*rand.Uniform()){
        if(true){        
           fprintf(OutputFile, "%8f %8f  \n", x, y );
           HGen->Fill(x,y);
           i++;
        }
     }


     fclose(OutputFile);

     HGen->Write();
}
