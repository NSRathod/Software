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

#define NX 1000       //Number of steps in x direction
#define NY 1         //Number of steps in y direction
#define NCH 49       //Number of channels
#define NSTR 32      //Number of straws
#define SX 0.2       //Single step in x direction [mm]
#define SY 200       //Single step in y direction [mm]
#define XSTART 0.0   //Starting position of scanner in x direction [mm]
#define YSTART 0.0   //Starting position of scanner in y direction [mm]
#define RSTR 5.05    //Straw radius in [mm]
#define Time 2       //Sleep time for scanner while measuement
#define ARYABHATA 0  //Owner

Double_t Scalars[NCH][NX];
Double_t Counts[NSTR][NX];

Double_t XSignals[NX]; 
Double_t YSignals[NX];


// data file
TString talyName = "m20_0deg_y200_0.2mm_SW0.11mm_14.12.18-SKAN1.dat";
ifstream talyInput;

int OpenFile()
{
// open file with data 
talyInput.open(talyName);

if(!talyInput.is_open()){
    cout << "Could not open"<< talyName << endl;
    return 0;
}
return 1;
}


int Readscan()
// read in data from one x scan
{
for(Int_t i=0; i<NX; i++){
    for(Int_t j=0; j<NCH; j++){
        talyInput >> Scalars[j][i];
        //cout << Scalars[j][i] <<" ";
    }
    cout <<endl;   
}
return 1;
}

int Translate()
// translate channel to straw number
{
Int_t istr;
for(Int_t ich=0; ich < NCH; ich++){    
    for(Int_t ix=0; ix<NX; ix++){
        // calculate the straw number istr corresponding to channel ich
        istr=-1;
        if((ich>=1) && (ich <=16)) istr=16-ich;
        if((ich>=17) && (ich <=32)) istr=48-ich;
        if((istr>=0) && (istr<NSTR)) {Counts[istr][ix]=Scalars[ich][ix];
    //cout << "counts " << Counts[0][ix] <<"--"<<istr<<"--"<<ix<<endl;
        }
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

//////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////-----Main code-----/////////////////////////////////////////////

bool SCANUPDATE(){
//------------Variable for Analysis----------
Double_t xpos[NX];
Double_t ypos;
Double_t MeanValues[NSTR];
Double_t AMP[NSTR];


//###############-----GRAPH|DECLARATION-----#####################    
    
    TGraph *g[NSTR];
    TMultiGraph *mg = new TMultiGraph();
    TMultiGraph *mgeven = new TMultiGraph();
    TMultiGraph *mgodd = new TMultiGraph();
//############################################################### 
    
// Creating a txt file to store different Mean-Values.....    
ofstream myfile("MEAN_VALUES_TEST.dat");

// Creating a root file to store different plots.....    
TFile* hfile = new TFile("ROOT-TEST.root","RECREATE");
// open data file    
OpenFile();
 
// main loop over scans in the y direction
for(Int_t iy=0; iy<NY; iy++)
{
// read in data for one x-scan
    Readscan();   
// translate channel to straw number
    Translate();
// Calculate x and y positions of the scanner
    for(Int_t ix=0; ix < NX; ix++)
    { 
    xpos[ix]= XSTART+(ix*SX);
    //cout<<"xpos = "<<xpos[ix]<<endl;
    }
    ypos=YSTART+(iy*SY);

// Make graphs of scans and fit them for each straw

TF1 *Ellipse[NSTR];
Int_t h=0;
 for(Int_t ich=0; ich < NSTR; ich++){    

//Fitting all the GRAPHS..... 
    Ellipse[ich] = new TF1(Form("Ellipse%i",ich),Dexter,0,200,3);

    Float_t midval = 29. + (5. * ich );
    if (ich > 1)
    {
        Ellipse[ich]->SetParameter(0, Ellipse[ich-2]->GetParameter(0));
        Ellipse[ich]->SetParLimits(0,1,73600);
        Ellipse[ich]->SetParameter(1, midval);
        Ellipse[ich]->SetParLimits(1,1,200);
        Ellipse[ich]->SetParameter(2, Ellipse[ich-2]->GetParameter(2));
        Ellipse[ich]->SetParLimits(2,0.8,1.2);
    }
    else
    {
        Ellipse[ich]->SetParameter(0, 5000);
        Ellipse[ich]->SetParLimits(0,1,73600);
        Ellipse[ich]->SetParameter(1, midval);
        Ellipse[ich]->SetParLimits(1,1,200);
        Ellipse[ich]->SetParameter(2, 0.89);
        Ellipse[ich]->SetParLimits(2,0.8,1.2);
    }
    //Ellipse->SetPoints(100);

//Plotting TGraph parameterters    
    g[ich] = new TGraph(NX,xpos,Counts[ich]);
    //g[ich]->SetLineColor(2);
    g[ich]->SetLineWidth(4);
    g[ich]->SetMarkerColor(4);
    g[ich]->SetMarkerStyle(21);
    g[ich]->SetTitle(Form("GRAPH%i",ich));
    g[ich]->SetName(Form("GRAPH%i",ich));

    Int_t status = 0;
    Int_t cnt = 0;
    printf("------------------GRAPH--------------- = %d\n", ich);
/*    do {
        TFitResultPtr rptr = g[ich]->Fit(Ellipse[ich],"", "", midval-10, midval+10);
        printf("Fit status = %d\n", int(rptr));
        status = int(rptr);
        --cnt;
        if(status == 0){
        AMP[ich] = Ellipse[ich]->GetParameter(0);
        MeanValues[ich] = Ellipse[ich]->GetParameter(1);
        myfile << AMP[ich] << endl;} 
    } while (cnt > 0 or status != 0);*/
    g[ich]->Write();
    mg->Add(g[ich]);

//     if (ich == 2) break;
   }
}
mg->Draw("APC");
mg->Write();

return true;
}
