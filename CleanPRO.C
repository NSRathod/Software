#include <iostream>
#include <fstream>
#include <TGraph.h>
#include <TCanvas.h>
#include "TLegend.h"
#include <TF1.h>
#include <TMath.h>
#include <vector>
using namespace std;

//#define STEP 0.2  //Smallest step while scanning.....
#define NSL 10.1    //Next Signal..... 
#define step 0.5    //Smallest step while scanning.....

#define max1 32     //Straws in 1 module.....
#define max2 15 

#define ExpXdiff 20.2

Double_t POS0deg1[max1];
Double_t POS0deg2[max1];
Double_t POS22deg[max1];

Double_t Expected_0deg1[max1];
Double_t Expected_0deg2[max1];
Double_t Expected_22deg[max1];

Double_t ODDPlot1[max1];
Double_t EVENPlot1[max1];

Double_t ODDPlot2[max1];
Double_t EVENPlot2[max1];

Double_t Deviation_0deg1[max1];
Double_t Deviation_0deg2[max1];
Double_t Deviation_22deg[max1];

Double_t Difference[max1];
Double_t Yaxis_0deg[max1];
Double_t Yaxis_22deg[max1];
Double_t Yaxis_diff[max1];

int ReadDATA()
{
// open file with data and read data
    TString talyName = "0.5mm_T5s_SW0.11mm_14.12.18-SKAN1.txt";
    ifstream talyInput;
    talyInput.open(talyName);

    if(!talyInput.is_open()){
        cout << "Could not open"<< talyName << endl;
        return kFALSE;
    }

// read in data
TString header;
for(Int_t i=0; i<1; i++){
    header.ReadLine(talyInput);
    }

// 0.5 mm as 1 step while scanning.......

for(Int_t i=0; i<max1; i++){
    
    talyInput >> POS0deg1[i] >> POS22deg[i];
    //cout << POS0deg1[i] <<"  "<< POS22deg[i] <<"\n";

    Expected_0deg1[i] = 49.6331 + (NSL*i);                  //Starting point of every signal has to be manually input i.e.49.4097 
    //Expected_0deg2[i] = 49.4097 + (10.1*i);
    if(i%2==0){Expected_22deg[i] =57.9496+ (NSL*i);}
    else{Expected_22deg[i] = 75.3066 - NSL + (NSL*i);}
//     cout << POS0deg1[i] <<"  "<< Expected_0deg1[i] <<"  "<< POS22deg[i] <<"  "<< Expected_22deg[i] <<"\n";
    
    Deviation_0deg1[i] = (Expected_0deg1[i] - POS0deg1[i])*step;
    //Deviation_0deg2[i] = (Expected_0deg2[i] - POS0deg2[i])*step;
    Deviation_22deg[i] = (Expected_22deg[i] - POS22deg[i])*step;
    //cout << Deviation_0deg1[i] <<"  "<< Deviation_22deg[i] <<"\n";
    
    if(i%2==0){Difference[i] = (POS22deg[i] - POS0deg1[i])*step;}
        else{Difference[i] = (POS22deg[i] - POS0deg1[i])*step;}

    Yaxis_0deg[i] = Difference[i]/TMath::Tan((22.5*3.14159)/180);
    Yaxis_diff[i] = (Yaxis_0deg[i%2==1] + Yaxis_0deg[i%2==0])/2;

    }
talyInput.close();
return 1;
}

Bool_t CleanPRO()
{
ReadDATA();

TFile* hfile = new TFile("TEST.root","RECREATE");

std::vector<double> ODD1;
std::vector<double> EVEN1;
std::vector<double> ODD2;
std::vector<double> EVEN2;

Double_t dummy[max1];
Double_t ST[] = {0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15};
Double_t n[] = {0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31};

TCanvas *c = new TCanvas("c","c");
c->SetGrid();

TCanvas *Diff = new TCanvas("Difference","Difference");
TCanvas *Tom = new TCanvas("Tom","Tom");

// Plotting graph for every point of analysis.....
TGraphErrors *GE0deg1= new TGraphErrors(max1,n,Deviation_0deg1,0,0);
// TGraphErrors *GE0deg2= new TGraphErrors(max1,n,Deviation_0deg2,0,0);
TGraphErrors *GE22deg= new TGraphErrors(max1,n,Deviation_22deg,0,0);

TMultiGraph *mg = new TMultiGraph();
cout << "----------------------------------------------------------" << endl;

GE0deg1->SetMarkerStyle(8);
GE0deg1->SetMarkerSize(1);
GE0deg1->SetMarkerColor(kRed);
GE0deg1->SetLineColor(kRed);
GE0deg1->SetLineWidth(2);

// GE0deg2->SetMarkerStyle(21);
// GE0deg2->SetMarkerSize(1);
// GE0deg2->SetMarkerColor(kGreen);
// GE0deg2->SetLineColor(kBlack);
// GE0deg2->SetLineWidth(2);

GE22deg->SetMarkerStyle(21);
GE22deg->SetMarkerSize(1);
GE22deg->SetMarkerColor(kBlack);
GE22deg->SetLineColor(kBlack);
GE22deg->SetLineWidth(2);

     mg->Add(GE0deg1,"PL");
     //mg->Add(GE0deg2,"PL");
     mg->Add(GE22deg,"PL");
     mg->GetXaxis()->SetRangeUser(0,33);
     mg->GetYaxis()->SetRangeUser(-0.5,1);
     c->cd();
     mg->Draw("AP");
     mg->Write();

     
   TLegend *leg = new TLegend(0.7,0.7,0.9,0.9);
   leg->SetFillColor(1);
   leg->AddEntry(GE0deg1,"0deg-SKAN1","lep");
   //leg->AddEntry(GE0deg2,"0deg-SKAN2","lep");
   leg->AddEntry(GE22deg,"22deg","lep");
   leg->SetFillStyle(0);
   leg->Draw("Deviation measurement");
   leg->SetName("Identity");
   leg->Write();
//c->SetLogx();

mg->GetXaxis()->SetTitle("NUMBER (n)");
mg->GetYaxis()->SetTitle("[ Expected - Meaured ] distance (DELTA_d)");
mg->SetTitle("Deviation in Expected and in measured position of STRAWS");
gPad->Modified();

Diff->cd();
TGraph *gr = new TGraph(max1,n,Yaxis_0deg);

    gr->SetLineColor(2);
    gr->SetLineWidth(4);
    gr->SetMarkerColor(4);
    gr->SetMarkerStyle(21);

    gr->SetTitle("Signal from mean Straw Position");
    gr->GetXaxis()->SetTitle("Straws Number");
    gr->GetYaxis()->SetTitle("Position between Slit and Straws [mm]");
    gr->Draw("AP");

//Difference with average
TGraphErrors *GOKU= new TGraphErrors(max1,n,Yaxis_0deg,0,0);

TGraphErrors *GOHAN= new TGraphErrors(max1,n,Yaxis_diff,0,0);

TMultiGraph *DBZ = new TMultiGraph();
cout << "----------------------------------------------------------" << endl;

GOKU->SetMarkerStyle(8);
GOKU->SetMarkerSize(1);
GOKU->SetMarkerColor(kRed);
GOKU->SetLineColor(kRed);
GOKU->SetLineWidth(2);

GOHAN->SetMarkerStyle(21);
GOHAN->SetMarkerSize(1);
GOHAN->SetMarkerColor(kGreen);
GOHAN->SetLineColor(kBlack);
GOHAN->SetLineWidth(2);

     DBZ->Add(GOKU,"P");
     DBZ->Add(GOHAN,"P");
     DBZ->GetXaxis()->SetRangeUser(0,33);
     DBZ->GetYaxis()->SetRangeUser(5,20);
     Tom->cd();     
     DBZ->Draw("AP");


    TLegend *Vegita = new TLegend(0.7,0.7,0.9,0.9);
    Vegita->SetFillColor(1);
    Vegita->AddEntry(GOKU,"Distance of straws in Yaxis from Slit","lep");
    Vegita->AddEntry(GOHAN,"Diff between Odd-Even straws","lep");
    Vegita->SetFillStyle(0);
    Vegita->Draw();

//c->SetLogx();

    DBZ->GetXaxis()->SetTitle("NUMBER (n)");
    DBZ->GetYaxis()->SetTitle("[ Expected - Meaured ] distance (DELTA_d)");
    DBZ->SetTitle("Deviation in Expected and in measured position of STRAWS");
    gPad->Modified();
    DBZ->Write();
    
//----------------------------------------------------------------------------------------------------------------------------------------//
//Odd-Even X axis fluctuation.......
TH1F *EVENCOUNTS0deg = new TH1F("EVENCOUNTS0deg","Difference in EVEN Straw Position at 0 degree",100,-1,1);
TH1F *ODDCOUNTS0deg = new TH1F("ODDCOUNTS0deg","Difference in ODD Straw Position at 0 degree",100,-1,1);
TH1F *EVENCOUNTS22deg = new TH1F("EVENCOUNTS22deg","Difference in EVEN Straw Position at 22.5 degree",100,-1,1);
TH1F *ODDCOUNTS22deg = new TH1F("ODDCOUNTS22deg","Difference in ODD Straw Position at 22.5 degree",100,-1,1);

ofstream myfile("TEST.dat");
for(int i = 0; i < max1; i+=2){
    cout << POS0deg1[i] <<"  "<< POS22deg[i] <<"\n";
    EVEN1.push_back(POS0deg1[i+2]-POS0deg1[i]);
    ODD1.push_back(POS0deg1[i+3]-POS0deg1[i+1]);
    cout << "EVEN1 :: "<<EVEN1[i/2] << "ODD1 :: "<<ODD1[i/2] << endl; 
    EVEN2.push_back(POS22deg[i+2]-POS22deg[i]);
    ODD2.push_back(POS22deg[i+3]-POS22deg[i+1]);
    //cout << "EVEN2 :: "<<EVEN2[i/2] << "ODD2 :: "<<ODD2[i/2] << endl; 
    }
    
for(Int_t i=0; i<15; i++){
    cout << "EVEN1 :: "<<EVEN1[i] << "ODD1 :: "<<ODD1[i] << endl; 
    ODDPlot1[i] = (ExpXdiff - ODD1[i]);
    EVENPlot1[i] = (ExpXdiff - EVEN1[i]);
    ODDPlot2[i] = (ExpXdiff - ODD2[i]);
    EVENPlot2[i] = (ExpXdiff - EVEN2[i]);
cout << "EVENPlot1 :: "<<EVENPlot1[i] << "ODDPlot1 :: "<<ODDPlot1[i] << "EVENPlot2 :: "<<EVENPlot2[i] << "ODDPlot2 :: "<<ODDPlot2[i] << endl; 
    ODDCOUNTS0deg->Fill(ODDPlot1[i]);
    EVENCOUNTS0deg->Fill(EVENPlot1[i]);
    ODDCOUNTS22deg->Fill(ODDPlot2[i]);
    EVENCOUNTS22deg->Fill(EVENPlot2[i]);
    }
ODDCOUNTS0deg->Draw();
EVENCOUNTS0deg->Write();
ODDCOUNTS22deg->Draw();
EVENCOUNTS22deg->Write();
    
TCanvas *JERRY = new TCanvas("JERRY","JERRY");
JERRY->SetFillColor(42);
JERRY->SetGrid();

TCanvas *MAGNETO = new TCanvas("MAGNETO","MAGNETO");
MAGNETO->SetFillColor(42);
MAGNETO->SetGrid();
    
TGraphErrors *gTalys1 = new TGraphErrors(max2,ST,ODDPlot1,0,0);
TGraphErrors *gTalys2 = new TGraphErrors(max2,ST,EVENPlot1,0,0);
TMultiGraph *mg1 = new TMultiGraph();
cout << "----------------------------------------------------------" << endl;

gTalys1->SetMarkerStyle(8);
gTalys1->SetMarkerSize(1);
gTalys1->SetMarkerColor(kSpring+2);
gTalys1->SetLineColor(kSpring-10);
gTalys1->SetLineWidth(2);

gTalys2->SetMarkerStyle(21);
gTalys2->SetMarkerSize(1);
gTalys2->SetMarkerColor(kBlue+2);
gTalys2->SetLineColor(kBlue-10);
gTalys2->SetLineWidth(2);

     mg1->Add(gTalys1,"PL");
     mg1->Add(gTalys2,"PL");
     mg1->GetXaxis()->SetRangeUser(0,15);
     mg1->GetYaxis()->SetRangeUser(-0.2,1);
     JERRY->cd();
     mg1->Draw("AP");
     mg1->Write();

     
 TLegend *leg1 = new TLegend(0.7,0.7,0.9,0.9);
   leg1->SetFillColor(1);
   leg1->AddEntry(gTalys1,"ODD_NO_STRAW","lep1");
   leg1->AddEntry(gTalys2,"EVEN_NO_STRAW","lep1");
   leg1->SetFillStyle(0);
   leg1->Draw();
   
//c->SetLogx();
mg1->GetXaxis()->SetTitle("NUMBER (n)");
mg1->GetYaxis()->SetTitle("DELTA_d");
//mg->GetXaxis()->SetRange(0,16);
//mg->SetTitle("Total analysis of Counts per Yaxis on Straw Tube [N/L] VS Amplitude [mV]");
     gPad->Modified();

TGraphErrors *gTalys3 = new TGraphErrors(max2,ST,ODDPlot2,0,0);
TGraphErrors *gTalys4 = new TGraphErrors(max2,ST,EVENPlot2,0,0);
TMultiGraph *mg2 = new TMultiGraph();
cout << "----------------------------------------------------------" << endl;

gTalys3->SetMarkerStyle(8);
gTalys3->SetMarkerSize(1);
gTalys3->SetMarkerColor(kSpring+2);
gTalys3->SetLineColor(kSpring-10);
gTalys3->SetLineWidth(2);

gTalys4->SetMarkerStyle(21);
gTalys4->SetMarkerSize(1);
gTalys4->SetMarkerColor(kBlue+2);
gTalys4->SetLineColor(kBlue-10);
gTalys4->SetLineWidth(2);

     mg2->Add(gTalys3,"PL");
     mg2->Add(gTalys4,"PL");
     mg2->GetXaxis()->SetRangeUser(0,15);
     mg2->GetYaxis()->SetRangeUser(-0.2,1);
     MAGNETO->cd();
     mg2->Draw("AP");
     mg2->Write();

 TLegend *leg2 = new TLegend(0.7,0.7,0.9,0.9);
   leg2->SetFillColor(1);
   leg2->AddEntry(gTalys3,"ODD_NO_STRAW","lep1");
   leg2->AddEntry(gTalys4,"EVEN_NO_STRAW","lep1");
   leg2->SetFillStyle(0);
   leg2->Draw();

//c->SetLogx();
mg2->GetXaxis()->SetTitle("NUMBER (n)");
mg2->GetYaxis()->SetTitle("DELTA_d");
//mg->GetXaxis()->SetRange(0,16);
//mg->SetTitle("Total analysis of Counts per Yaxis on Straw Tube [N/L] VS Amplitude [mV]");
     gPad->Modified();

return kTRUE;
}
