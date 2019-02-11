#include <iostream>
#include <fstream>
#include <TGraph.h>
#include <TCanvas.h>
#include "TLegend.h"
#include <TF1.h>
#include <TH1.h>
#include <TMath.h>
#include <TObject.h>
#include "TROOT.h"
#include "TEnv.h"
#include "TBrowser.h"
#include "TMultiGraph.h"
#include "TFile.h"
using namespace std;

#define NS 400      //Number of steps or points in Scanning
#define CH 49       //Number of Channels
#define SLSM 16     //Number of Single Layer Straws Module
#define DLSM 32     //Number of Double Layer Straws Module
#define MD 2        //Number of Module
#define SS 0.5      //Single step in [mm]
#define Time 2      //Sleep time for scanner while measuement
#define ARYABHATA 0  //Owner

Double_t Scalars[CH][NS];


int ReadDATA()
{
// open file with data and read data
    TString talyName = "m20_0deg_y200_0.5mm_T5s_14.12.18-SKAN1.dat";
    ifstream talyInput;
    talyInput.open(talyName);

    if(!talyInput.is_open()){
        cout << "Could not open"<< talyName << endl;
        return 0;
    }

// read in data
    for(Int_t i=0; i<NS; i++){
        for(Int_t j=0; j<CH; j++){
        talyInput >> Scalars[j][i];
        //cout << Scalars[j][i] <<" ";
        }
    cout <<endl;   
    }
return 1;
}

//Fitting function.....
double Dexter(double* x, double* par){
    double X = x[0];
    double var = ((par[0]*par[0])-((X - par[1])*(X - par[1]))/((1 - par[2])*(1-par[2])));
    if(var >= 0) {var = sqrt(var);}
    else {var = 0;}
    return var;
    }


Bool_t Scan_Tech(){
//------------Variable for Analysis----------
    const Int_t max1 = DLSM;
    const Int_t max2 = NS;

    Double_t Signals[NS]; 
    Double_t XSignals[NS]; 
    Double_t YSignals[NS];
    Double_t MeanValues[DLSM]={0};
    Double_t dummy[NS];

// Creating a txt file to store different Mean-Values.....    
    ofstream myfile("MEAN_VALUES_m20_0deg_y200_0.5mm_T5s_14.12.18-SKAN1.dat");

    // Creating a root file to store different plots.....    
    TFile* hfile = new TFile("O-manytest-its-one-of-them.root","RECREATE");

//#######################-----CANVAS-----#########################    
    TCanvas *c1 = new TCanvas("c1","Individual Straw Signal",1500,1200);
    c1->SetFillColor(42);
    c1->SetGrid();
    c1->Divide(8,4);
    
    TCanvas *c2 = new TCanvas("c2","ALL Straws for fitting",1500,1200);
    c2->SetFillColor(42);
    c2->SetGrid();
//###############################################################    
//###############-----GRAPH|DECLARATION-----#####################    
    
    TGraph *g[max1];

    TMultiGraph *mg = new TMultiGraph();
    TMultiGraph *mgeven = new TMultiGraph();
    TMultiGraph *mgodd = new TMultiGraph();
//############################################################### 
       
//Loop creating values for X-axis.....    
for(Int_t k=0; k<NS; k++){
    XSignals[k]= k;
    }

ReadDATA();
Int_t h = 0; 
Int_t j = 0;
for (int i= 0; i<MD; i++){   
    for (int q= (1+(SLSM*i)); q<(SLSM+1+(SLSM*i)); q++){
        h++;
        for (int w= 0; w< NS; w++){
            YSignals[w] = Scalars[(SLSM+1+(SLSM*i))-q+(SLSM*i)][w];
            cout<<YSignals[w]<<"--------"<<i<<"  "<<"Straw#:::"<<h<<"  "<<"j:::"<<j<<endl;
            j++;
        }
    c1->cd(q);
    TGraph *gr = new TGraph(NS,XSignals,YSignals);
    
    gr->SetLineColor(2);
    gr->SetLineWidth(3);
    gr->SetMarkerColor(4);
    gr->SetMarkerStyle(21);
    gr->SetTitle(Form("GRAPH%i",q));
    gr->SetName(Form("GRAPH%i",q));
    gr->SetTitle("Individual Straw Signal");
    gr->GetXaxis()->SetTitle("Straw number");
    gr->GetYaxis()->SetTitle("Scalars");
    gr->Draw("AP");

    c2->cd();
    g[q] = new TGraph(NS,XSignals,YSignals);
    //g[q]->SetLineColor(2);
    //g[q]->SetLineWidth(4);
    g[q]->SetMarkerColor(4);
    g[q]->SetMarkerStyle(21);
    g[q]->SetTitle(Form("GRAPH%i",q));
    g[q]->SetName(Form("GRAPH%i",q));
    
//Fitting all the GRAPHS..... 
    TF1 *Ellipse = new TF1("Ellipse",Dexter,0,400,3);
    Ellipse->SetParameters(735,56+(8*q),0.9);
//    Ellipse->SetParameters(735,140+(24*q),0.9);
//     if(q%2==0){Ellipse->SetParameters(735,140+(30*q),0.09);}
// 	else{Ellipse->SetParameters(735,140+(45*q),0.09);}
    Ellipse->SetParLimits(0,1,73600);
    Ellipse->SetParLimits(1,1,400);
    Ellipse->SetParLimits(2,0.8,1);
    g[q]->Fit(Ellipse,"R");
    MeanValues[q] = Ellipse->GetParameter(1);.
    g[q]->Write();
    if(q%2==0){mgodd->Add(g[q]);}
    else{mgeven->Add(g[q]);}
    mg->Add(g[q]);
    
//    cout << "Mean value for GRAPH%q :::: " << MeanValues[q] << endl;
    myfile << MeanValues[q] << endl;
     }
    mg->Write("All Signals at same time");
}
    mg->Draw("AP");
    mgodd->Draw("AP");
    mgodd->Write("Signals-@-ODD-Straws");
    mgeven->Draw("AP");
    mgeven->Write("Signals-@-EVEN-Straws");

    return kTRUE;  
}
