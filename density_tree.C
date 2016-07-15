#include <string>
#include <vector>
#include <exception>
#include <vector>
#include <utility>  
#include <fstream>

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

class box{
//    A box is characterized by 
//      - a list of events in the box (= list of integers)
//      - a set of boudaries (up and down values for each feature)
//          [[down_0,up_0], [down_1, up_1], ...]
//      - perhaphs two daughter boxes
  public:
     std::vector<unsigned int> events;
     std::vector<double> edgesDown;
     std::vector<double> edgesUp;
     std::vector<box> daughters;
     unsigned int nEvents;
     unsigned int nObs;
     int cutVariable;
     double cutValue;
     double deviation;
  public:
    box(){};

    box( const std::vector<double*>& allEvents, std::vector<unsigned int>& events_, unsigned int Nvar=0){
       events = events_;
       for(unsigned int v=0;v<Nvar;v++){ edgesDown.push_back(-1);   edgesUp.push_back(-1);   }
       nEvents = events.size();
    }

    box( const std::vector<double*>& allEvents, unsigned int Nvar=0){
       for(unsigned int e=0;e<allEvents.size();e++){events.push_back(e);}
       for(unsigned int v=0;v<Nvar;v++){ edgesDown.push_back(-1);   edgesUp.push_back(-1);   }
       nEvents = events.size();
    }

    box( const std::vector<double*>& allEvents, std::vector<unsigned int>& events_,  std::vector<double>& edgesDown_, std::vector<double>& edgesUp_){
       events = events_;
       edgesDown = edgesDown_;
       edgesUp = edgesUp_;        
       nEvents = events.size();
       cutVariable = -1;
       //printf("create a box with edge %f %f %f %f\n", edgesDown[0],edgesDown[1],edgesUp[0],edgesUp[1]);
    }

    box(ifstream& file, unsigned int nVars=0, bool first=true){
       if(first){
          for(unsigned int v=0;v<nVars;v++){
             edgesDown.push_back(-1E99);
             edgesUp.push_back(1E99);
          }         
       }

       string line;
       char tag; unsigned int I;  double D;
        streampos pos = file.tellg();
        getline(file,line);
        switch(line[0]){
             case 'C':{ //get a branch
                unsigned int ncharRead;  char* linePtr=(char*)line.c_str();
                sscanf(linePtr, "%c %u %lf %n", &tag, &I, &D, &ncharRead); linePtr+=ncharRead;
                for(unsigned int v=0;v<nVars;v++){sscanf(linePtr, "%lf %lf %n", &edgesDown[v], &edgesUp[v], &ncharRead);  linePtr+=ncharRead;}

                cutVariable = I;
                cutValue = D;
                std::vector<double> edgesDown_a, edgesDown_b;
                std::vector<double> edgesUp_a, edgesUp_b;
                for(unsigned int v=0;v<edgesDown.size();v++){
                   edgesDown_a.push_back(edgesDown[v]);
                   edgesUp_a.push_back(edgesUp[v]);
                   edgesDown_b.push_back(edgesDown[v]);
                   edgesUp_b.push_back(edgesUp[v]);
                }
                edgesUp_a[cutVariable] = cutValue;
                edgesDown_b[cutVariable] = cutValue;
                daughters.push_back(box(file, nVars));  daughters[0].edgesDown = edgesDown_a;   daughters[0].edgesUp = edgesUp_a;
                daughters.push_back(box(file, nVars));  daughters[1].edgesDown = edgesDown_b;   daughters[1].edgesUp = edgesUp_b;
                nEvents =0; for(unsigned int d=0;d<daughters.size();d++){nEvents+=daughters[d].nEvents;}
                return;
                }break;

             case 'L':{ //get a leaf. exit after we filled the information
                sscanf(line.c_str(), "%c %u %lf", &tag, &I, &D);
                nEvents = I;
                return; //don't do anything else
                }break;

             default:{
                //go back to the previous vector position and return;
                file.seekg(pos);
                return;
                }break;
       }
    }



    double vol(){
        double vol=1.0;
        for(unsigned int e=0;e<edgesDown.size();e++){vol*=(edgesUp[e]-edgesDown[e]);}
        return vol;
    }

    double getWeight(double* oneEvent){
       if(daughters.size()==0){
          return events.size() / vol();
       }else{
          if(oneEvent[cutVariable]<cutValue){
             return daughters[0].getWeight(oneEvent);
          }else{
             return daughters[1].getWeight(oneEvent);
          }
       }
    }

    void reccord(std::string path){
        FILE* pFile = fopen(path.c_str(), "w");
        if(pFile){
           reccord(pFile);
           fclose(pFile);
        }
    }

    void reccord(FILE* pFile){
        if(daughters.size()==0){
           fprintf(pFile, "L %u %g \n", nEvents, vol());
        }else{
           fprintf(pFile, "C %u %g ", cutVariable, cutValue);
           for(unsigned int v=0;v<edgesUp.size();v++){fprintf(pFile, "%f %f ", edgesDown[v], edgesUp  [v]);}fprintf(pFile, "\n");
           for(unsigned int d=0;d<daughters.size();d++){  
              daughters[d].reccord(pFile);
           }
        }
    }



   void getBoundaries(const std::vector<double*>& allEvents){
     for(unsigned int v=0;v<edgesUp.size();v++){
         edgesDown[v]=allEvents[events[0]][v];
         edgesUp  [v]=allEvents[events[0]][v];
      }

      for(unsigned int e=0;e<events.size();e++){
         for(unsigned int v=0;v<edgesUp.size();v++){
            if(allEvents[events[e]][v]>edgesUp  [v])edgesUp  [v]=allEvents[events[e]][v];
            if(allEvents[events[e]][v]<edgesDown[v])edgesDown[v]=allEvents[events[e]][v];
         }
      }

      for(unsigned int v=0;v<edgesUp.size();v++){
         edgesUp[v]+=(edgesUp[v]-edgesDown[v])*1e-6;
         edgesDown[v]-=(edgesUp[v]-edgesDown[v])*1e-6;
      }
   }



   unsigned int getHighestGradientVariable(const std::vector<double*>& allEvents){
      //Identify the feature for which the gradient is the largest
     double grad_max=0;
     unsigned int good_var=0;
     for(unsigned int v=0;v<edgesUp.size();v++){
         double down=edgesDown[v];
         double up  =edgesUp[v];
         double step=(up-down)/5.0;
         double a[]={0.0,0.0,0.0,0.0,0.0};
         
         for(unsigned int e=0;e<events.size();e++){
            if(allEvents[events[e]][v]<down+step){
               a[0]+=1;
            }else if(allEvents[events[e]][v]<down+step*2.0){
               a[1]+=1;
            }else if(allEvents[events[e]][v]<down+step*3.0){
               a[2]+=1;
            }else if(allEvents[events[e]][v]<down+step*4.0){
               a[3]+=1;
            }else{
               a[4]+=1;
            }
         }

         double grad=1e-16;
         if(a[0]>0 && a[1]>0)  grad+=fabs((a[1]-a[0])/(a[0]+a[1]));
         if(a[1]>0 && a[2]>0)  grad+=abs((a[2]-a[1])/(a[1]+a[2]));
         if(a[2]>0 && a[3]>0)  grad+=abs((a[3]-a[2])/(a[2]+a[3]));
         if(a[3]>0 && a[4]>0)  grad+=abs((a[4]-a[3])/(a[3]+a[4]));

         //printf("var=%i %f %f %f %f --> grad = %f\n", v, a[0], a[1], a[2], a[3], grad);
 
         if(grad_max<grad){
            grad_max=grad;
            good_var=v;
         }
     }
     return good_var;
   }

   void generateDaughters(const std::vector<double*>& allEvents, int nLayers){
     if(nLayers<1)return;

     unsigned int var=getHighestGradientVariable(allEvents);

      //define two equal size boxes by dividing the current box along the variable "var"
     std::vector<double> valvect;
     for(unsigned int e=0;e<events.size();e++){
        valvect.push_back(allEvents[events[e]][var]);
     }
     std::sort(valvect.begin(), valvect.end());

//     std::map<double, unsigned int>::iterator val2idIterator = val2id.begin();
//     for(unsigned int i=0;i<val2id.size()/2;i++){val2idIterator++;}
//     cutValue    = val2idIterator->first; //simplified with respect to Pierre implementation, don't average between two events.
     cutValue = valvect[(valvect.size()+1)/2];
     cutVariable = var;

     //create the list of events associated with daughters a and b
     std::vector<unsigned int> events_a;
     std::vector<unsigned int> events_b;
     for(unsigned int e=0;e<events.size();e++){
       if(allEvents[events[e]][var]<cutValue){
          events_a.push_back(events[e]);
       }else{
          events_b.push_back(events[e]);
       }
     }

//    create the boundaries associated with daughters a and b
//    PAY ATTENTION THAT THE BOUNDARIES IN self.boundaries, boundaries_a and 
//    boundaries_b DO NOT HAVE THE SAME ADDRESS
     std::vector<double> edgesDown_a, edgesDown_b;
     std::vector<double> edgesUp_a, edgesUp_b;    
     for(unsigned int v=0;v<edgesDown.size();v++){
         edgesDown_a.push_back(edgesDown[v]);
         edgesUp_a.push_back(edgesUp[v]);
         edgesDown_b.push_back(edgesDown[v]);
         edgesUp_b.push_back(edgesUp[v]);
     }
     edgesUp_a[var] = cutValue;
     edgesDown_b[var] = cutValue;

     //printf("mother box with edge %f %f %f %f split on var%i at %f\n", edgesDown[0],edgesDown[1],edgesUp[0],edgesUp[1], var, cutValue);


     daughters.clear(); 
     daughters.push_back(box(allEvents, events_a, edgesDown_a, edgesUp_a));
     daughters.push_back(box(allEvents, events_b, edgesDown_b, edgesUp_b));
     daughters[0].generateDaughters(allEvents, nLayers-1);
     daughters[1].generateDaughters(allEvents, nLayers-1);
   }


   void drawBoxes(const std::vector<double*>& allEvents, int& index, string outputPath="", bool firstBox=true){
      TCanvas* c1 = NULL;

      if(outputPath!=""){
         c1 = new TCanvas("c1", "Tree",500,500);
         c1->SetLogz(true);
      }

      if(firstBox){  
         TH2D* histo = new TH2D("histo", "histo;x;y", 100, edgesDown[0], edgesUp[0], 100, edgesDown[1], edgesUp[1]);        
         for(unsigned int e=0;e<events.size();e++){
            histo->Fill(allEvents[events[e]][0], allEvents[events[e]][1]);
         }
         histo->SetTitle("");
         histo->SetStats(kFALSE);
         histo->Draw("COLZ");
      }

      if(daughters.size()==0){
//         TBox* box1 = new TBox(edgesDown[0],edgesDown[1],edgesUp[0],edgesUp[1]);
         TPaveText* box1 = new TPaveText(edgesDown[0],edgesDown[1],edgesUp[0],edgesUp[1]);         
         char text[256]; sprintf(text, "%i - %i", index, (unsigned int)events.size());
         box1->AddText(text);
         box1->SetTextAlign(22);
         box1->SetTextSize(0.02);
         printf("draw boxes %f %f %f %f c=%i %i\n", edgesDown[0],edgesDown[1],edgesUp[0],edgesUp[1], index, (unsigned int)events.size());
         box1->SetLineColor(1);
         box1->SetLineWidth(2);
         box1->SetFillStyle(0);
         box1->Draw("same l");
         index++;
      }else{
         for(unsigned int d=0;d<daughters.size();d++){
            daughters[d].drawBoxes(allEvents, index, "", false);
         }
      }

      if(outputPath!="")c1->SaveAs(outputPath.c_str());      
   }


    void clearContent(){
       events.clear();
       nObs = 0;
       for(unsigned int d=0;d<daughters.size();d++){
          daughters[d].clearContent();
       }
    }

    void incObsEventBox(double* oneEvent){
      nObs++;
      if(daughters.size()>0){
         int d=(oneEvent[cutVariable]<cutValue)?0:1;
         daughters[d].incObsEventBox(oneEvent);
      }
    }

    void setDeviation(const std::vector<double*>& obsEvents, unsigned int nEventsInTree=0){
       if(nEventsInTree<=0){
          clearContent();
          for(unsigned int e=0;e<obsEvents.size();e++){
             incObsEventBox(obsEvents[e]);
          }
          nEventsInTree = nEvents;
       }

       if(nEventsInTree>0){
          deviation =  1 - ( (double(nEvents)/double(nEventsInTree)) / (double(nObs)/double(obsEvents.size())) );
          for(unsigned int d=0;d<daughters.size();d++){
             daughters[d].setDeviation(obsEvents, nEventsInTree);
          }          
       }
    }

    double getDeviation(double* oneEvent){
      if(daughters.size()>0){
         int d=(oneEvent[cutVariable]<cutValue)?0:1;
         return daughters[d].getDeviation(oneEvent);
      }
      return deviation;            
    }

    unsigned int getDistance(double* eventA, double* eventB){
      if(daughters.size()>0){
         //are they both in the same branch?
         if(eventA[cutVariable]<cutValue && eventB[cutVariable]<cutValue){
            return 1 + daughters[0].getDistance(eventA, eventB);
         }else if(eventA[cutVariable]>cutValue && eventB[cutVariable]>cutValue){
            return 1 + daughters[1].getDistance(eventA, eventB);
         }
         return 1;
      }
      return 0;
    }
};

////////////////////////////////////////////////////////////////////////////////////////////////

class superbox {
  public:
     unsigned int nTrees, nLayers, nVars;
     std::vector<box> boxes;
  public:
    superbox(){};
    superbox(const std::vector<double*>& allEvents, unsigned int nTrees_, unsigned int nLayers_, unsigned int nVars_){
        nTrees  = nTrees_;
        nLayers = nLayers_;
        nVars = nVars_;
        for(unsigned int t=0;t<nTrees;t++){
           std::vector<unsigned int> events;
//         for(unsigned int e=t;e<allEvents.size();e+=nTrees){events.push_back(e);}         
           unsigned int NEvent = allEvents.size()/nTrees;           
           for(unsigned int e=t*NEvent;e<(t+1)*NEvent;e++){events.push_back(e);}           
           box newBox(allEvents, events, nVars);
           newBox.getBoundaries(allEvents); 
           newBox.generateDaughters(allEvents, nLayers);
           boxes.push_back(newBox);           
        }
    }

    superbox(ifstream& file, unsigned int nVars_=0){
       nVars = nVars_;
       nTrees=0;       
       string line;

       getline(file,line);//skip the first line that contains the number of trees
       while(getline(file,line) ){
        switch(line[0]){
             case 'T':{ //get a tree
                nTrees++;
                boxes.push_back(box(file, nVars, true));
                }break;
             default:{
                printf("# # # # # # # # #\nsomething is not ok in the superTree input file\n# # # # # # # # #\n");
                return;
                }break;
        }
      }
    }

    void drawBoxes(const std::vector<double*>& allEvents, string outputPath=""){
        TCanvas* c1 = NULL;
        if(outputPath!=""){
           c1 = new TCanvas("c1", "Tree",800, 800 );
           c1->DivideSquare(std::min(25,(int)boxes.size()));
           c1->SetLogz(true);
        }
        for(unsigned int t=0;t<std::min(25,(int)boxes.size());t++){
           (c1->cd(t+1))->SetLogz(true);
           int index=0;
           boxes[t].drawBoxes(allEvents, index);
        }
        if(outputPath!="")c1->SaveAs(outputPath.c_str());      
    }


     void reccord(std::string path){
        FILE* pFile = fopen(path.c_str(), "w");
        if(pFile){
           reccord(pFile);
           fclose(pFile);
        }
     }

     void reccord(FILE* pFile){
        fprintf(pFile, "Number of trees: %i\n", (unsigned int) boxes.size());
        for(unsigned int t=0;t<boxes.size();t++){
           fprintf(pFile, "T\n");
           boxes[t].reccord(pFile);
        }
     }
    
    double getWeight(double* oneEvent){
       double weight=0.0;
        for(unsigned int t=0;t<boxes.size();t++){
           weight += boxes[t].getWeight(oneEvent);
        }return weight/boxes.size();
    }

    void setDeviation(const std::vector<double*>& obsEvents){
        for(unsigned int t=0;t<boxes.size();t++){
           boxes[t].setDeviation(obsEvents);
        }
    }

    double getDeviation(double* oneEvent){
        double deviation=0;
        for(unsigned int t=0;t<boxes.size();t++){
           deviation+=boxes[t].getDeviation(oneEvent);
        }return deviation/boxes.size();       
    }

    //compute the correlation function for the tree
    double computeTwoPointFct(const std::vector<double*>& obsEvents, TH1* histo){
       double corr=0;

       printf("Progressing Bar                   :0%%       20%%       40%%       60%%       80%%       100%%\n");
       printf("Progressing Bar                   :");
       int TreeStep = obsEvents.size()/50;if(TreeStep==0)TreeStep=1;
       for(unsigned int i=0;i<obsEvents.size();i++){
          if(i%TreeStep==0){printf(".");fflush(stdout);}
          for(unsigned int j=i+1;j<obsEvents.size();j++){
               double distance=0;   double FF=0;
               for(unsigned int t=0;t<boxes.size();t++){
                  int d = nLayers - boxes[t].getDistance(obsEvents[i], obsEvents[j]);
                  distance += d;
                  double ff = boxes[t].getDeviation(obsEvents[i])*boxes[t].getDeviation(obsEvents[j]);
                  FF       += ff;
//                  if(ff<0)printf("distance = %i  corr=%f\n", d, ff);
                  //printf("distance = %i  corr=%f\n", d, ff);
               }
               distance/=(boxes.size()*nLayers);
               FF/=boxes.size();
               if(histo)histo->Fill(distance,FF);
               corr+=distance*FF;
           }
        }printf("\n");
        double factor=2.0/pow(obsEvents.size(), 2);
//        if(histo)histo->Scale(factor);
        return corr*factor;
    }
};


void readEventsFromFile(string inputFile, std::vector<double*>& events){
  events.clear();

  ifstream file (inputFile.c_str());
  if (file.is_open()){
    string line;
    while ( getline (file,line) ){
       double* array = new double[2];
       sscanf(line.c_str(), "%lf %lf", &array[0], &array[1]);
       events.push_back(array);
    }
    file.close();   
  }
}






void density_tree(string inputFilePath="in.txt", string unusedArgument="4")
{
     if(inputFilePath == "")return;

     int nLayers = 4;
     sscanf(unusedArgument.c_str(),"%i", &nLayers);
     printf("will consider %i levels per tree\n", nLayers);
     
     std::vector<double*> allEvents;
     readEventsFromFile(inputFilePath, allEvents);
     printf("EventSize = %i\n",(unsigned int) allEvents.size() );

/*
     box mainBox(allEvents, 2);
     mainBox.getBoundaries(allEvents); 
     mainBox.generateDaughters(allEvents, nLayers);
     mainBox.reccord("savedTree");
//     int index=1;
//     mainBox.drawBoxes(allEvents, index, "test.png");

     ifstream file ("savedTree");
     if (file.is_open()){
        box newBox(file,2);
        newBox.reccord("savedTree2");
     }
*/



       superbox sboxes(allEvents, 100, 6, 2);
       sboxes.drawBoxes(allEvents, "test.png");
       sboxes.reccord("savedTree");
//     ifstream file ("savedTree");
//     if (file.is_open()){
//        superbox newBox(file,2);
//        newBox.reccord("savedTree2");
//     }


       allEvents.clear();
       readEventsFromFile("circle", allEvents);

       std::vector<double*> obsEventsA;  
       std::vector<double*> obsEventsB;  
       std::vector<double*> obsEventsC;  
       for(int i=0;i<min((int)allEvents.size(), 5000);i++){
          if(i<5000)obsEventsA.push_back(allEvents[i]);
          if(i<1000)obsEventsB.push_back(allEvents[i]);
          if(i< 500)obsEventsC.push_back(allEvents[i]);
       }
 
       TCanvas*  c1 = new TCanvas("c1", "Tree",500,500);
//       c1->SetLogy(true);
       TH1D* histoA = new TH1D("histoA", "histo;d;corr", 50, 0, 1);        histoA->SetLineColor(1);
       TH1D* histoB = new TH1D("histoB", "histo;d;corr", 50, 0, 1);        histoB->SetLineColor(4);
       TH1D* histoC = new TH1D("histoC", "histo;d;corr", 50, 0, 1);        histoC->SetLineColor(2);
       histoA->SetTitle("");
       histoA->SetStats(kFALSE);
       sboxes.setDeviation(obsEventsA);
       printf("corr = %f\n", sboxes.computeTwoPointFct(obsEventsA, histoA)  );
       sboxes.setDeviation(obsEventsB);
       printf("corr = %f\n", sboxes.computeTwoPointFct(obsEventsB, histoB)  );
       sboxes.setDeviation(obsEventsC);
       printf("corr = %f\n", sboxes.computeTwoPointFct(obsEventsC, histoC)  );
       histoA->Draw("E1 hist");
       histoB->Draw("E1 hist same");
       histoC->Draw("E1 hist same");
       c1->SaveAs("correlation.png");      
       c1->SaveAs("correlation.C");      
 
   

/*     string OutFilePath = "results.txt";
     TFile* OutputHisto = new TFile((OutFilePath+".root").c_str(),"RECREATE");
     OutputHisto->cd();
     TH2F* HGen         = new TH2F("pdf", "pdf;x;y", 100,0,1, 100, 0, 1 );


     HGen->Write();
*/
}
