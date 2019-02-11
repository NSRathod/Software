#include <iostream> 
#include <fstream>
#include <TMath.h>
#include <TGraph.h>
#include <sstream>
#include <string>
#include "TCanvas.h"
#include "TF1.h"
#include "TH1.h"
#include <TFile.h>

double myfun1(double* x, double* par){
 
double X =x[0];

return par[0]*TMath::Exp((-pow((X-par[1]),3)/pow(par[2],2))-(par[3]*(X-par[1])/par[2]));//TMath::Sqrt(2*TMath::Pi()*par[2]);
//return TMath::Exp(-par[0])*pow(par[0],X)/TMath::Factorial(10);
}


/*
double myfun1(double* x, double* par){
 
double X =x[0];
//return par[0]*TMath::Exp(-pow((X-par[1]),2)/2*pow(par[2],2));//TMath::Sqrt(2*TMath::Pi()*par[2]);
//return TMath::Exp(-par[0])*pow(par[0],X)/TMath::Factorial(10);
return (par[0] * TMath::Sqrt(36-((X-64.6-10.1)*(X-64.6-10.1)/par[1]/par[1])))+200;
}*/
double Straw(double* x, double* par){
    double X =x[0];
    return 7000*TMath::Sqrt((par[0]*par[0] - (par[1]-X)*(par[1]-X)));
//     return (par[1]*X)+par[0];
    }
    
    
double firstexpo(double* x, double* par){
    double X =x[0];
   return TMath::Exp(par[0]+(par[1]*X));
//      return (par[1]*X)+par[0];
    }

double lastexpo(double* x, double* par){
    double X =x[0];
   return TMath::Exp(par[0]+(par[1]*X));
//      return (par[1]*X)+par[0];
        }
        
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
double myfun(double* x, double* par){
    double X =x[0];
    double var = ((par[0]*par[0])-((X - par[1])*(X - par[1]))/((1-par[2])*(1-par[2])));
    if(var >= 0) {var=sqrt(var);}
    else {var = 0;}
    
    return var;
        
}
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
        
        
        
// double myfun(double* x, double* par){
//     double X =x[0];
//     //return (par[0]*TMath::Sqrt(1-par[1]*par[1]) * TMath::Sqrt(1-((X-par[2])*(X-par[2])/par[0]/par[0])))+par[3];
//     //return ((1 / TMath::Sqrt(1-(par[0]*par[0]))) * TMath::Sqrt(par[1]*par[1] - ((X - par[2])*(X - par[2]))))+par[3]; 
//     //return ((1 / TMath::Sqrt(1-(par[0]*par[0]))) * TMath::Sqrt(par[1]*par[1] - ((X - par[2])*(X - par[2]))))+par[3]; 
//     //return (par[0] * TMath::Sqrt(1 -(((X - par[1]) * (X - par[1]))/(par[2] * par[2])))) + par[3];
//     return (TMath::Sqrt((par[0]*par[0])-((X - par[1])*(X - par[1]))/((1-par[2])*(1-par[2]))));
//         
// }

double Ellip(double* x, double* par){
    double X =x[0];

    return (par[0]*TMath::Sqrt(1-((X - par[1])*(X - par[1]))/(par[2]*par[2])));
        
}

double totalfun(double* x, double* par){
    double X =x[0];
    //return (par[0]*TMath::Sqrt(1-par[1]*par[1]) * TMath::Sqrt(1-((X-par[2])*(X-par[2])/par[0]/par[0])))+par[3];
    //return ((1 / TMath::Sqrt(1-(par[0]*par[0]))) * TMath::Sqrt(par[1]*par[1] - ((X - par[2])*(X - par[2]))))+par[3]; 
    //return ((1 / TMath::Sqrt(1-(par[0]*par[0]))) * TMath::Sqrt(par[1]*par[1] - ((X - par[2])*(X - par[2]))))+par[3]; 
    //return (par[0] * TMath::Sqrt(1 -(((X - par[1]) * (X - par[1]))/(par[2] * par[2])))) + par[3];
    return firstexpo(x,par)*myfun(x,&par[2])*lastexpo(x,&par[5]);
}
        

// double Total(double* x){
//     double X3 =x[0];
//     return 
//         }
/*
 double conf_hyperg(double* num, double* mu) {
 double N = num[0];
     
return (TMath::Gamma(N+1)*TMath::Gamma(mu[0]+1)/TMath::Gamma(N-mu[0]+1))*(TMath::Gamma(mu[0]+mu[1]/mu[2])*TMath::Gamma(1/mu[2])*TMath::Gamma(N+(1-mu[1])/mu[2]-mu[0]))/(TMath::Gamma(N+(1/mu[2])) * TMath::Gamma(mu[1]/mu[2]) * TMath::Gamma((1-mu[1])/mu[2]))*75000;
 
 }*/
void Reallycoolfit(){
    for(int i = 1; i<16; i+=2)
    { 
        cout<< "value if the loop :::::::  "<<i <<endl;
    }
    
    TCanvas * c1 = new TCanvas("c1","Canvas");
    //TCanvas * c2 = new TCanvas("c2","Canvas2");


    //canvas = (TCanvas*)fin->Ge//     func->SetParLimits(0,1,73600);
//     func->SetParLimits(1,1,201);
//     func->SetParLimits(2,0.8,1);t("Graph0");  
    //TF1 *fHisto = (TCanvas*)fin->Get("Graph0");



    //TF1 *func =new TF1("func",myfun, 79.6-10.1, 89.7-10.1,3);//[0]*TMath::Exp(-((x[0]-[1])*(x[0]-[1]))/([2]*[2]))
    //----------------------------------------------------------------------------------------------------------------------------------
    //TF1 *func =new TF1("func",myfun, 80.4-10.1-10.1-10.1, 89.68-10.1-10.1-10.1,2);
    /*
    //func->SetParameters(1,12,84.8,200);
    func->SetParLimits(0,1,100);
    func->SetParLimits(1,1,1000);
    func->SetParLimits(2, 79.8, 89.9);
    func->SetParLimits(3,1,26);
    func->SetParNames("Eccentricity","Minor-axis","H-center","V-center");
    //-----------------------------------------------------------------------------------------------------------------------------------
    func->SetParameters(72000,84.8,22,250);
    func->SetParLimits(0,1000,68000);
    //func->SetParLimits(1, 78.8, 89.9);
    func->SetParLimits(1, 1, 100);
    func->SetParLimits(2,1,22);
    func->SetParLimits(3,1,290);
    func->SetParNames("V - Major-axis","H-center","H-centerH - Minor-axis","V-center");*/
    //--------------------------------------------------------------------------------------------------------------------------------
    //TF1 *func =new TF1("func",myfun, 79.8-10.1-10.1, 89.7-10.1-10.1,3);
    //TF1 *func =new TF1("func","pol0[0]+gaus[1]+pol0[3]", 79.546-10.1-10.1, 89.738-10.1-10.1,3);
    //TF1 *func =new TF1("func",myfun, 79.6-10.1-10.1, 89.7-10.1-10.1,3);
    //TF1 *func =new TF1("func",myfun, 50.19-10.1, 58.8-10.1,3);
    
    TF1 *func1 =new TF1("func1",firstexpo, 60, 80, 2);
    func1->SetParameters(-11,2);
//     func1->SetParLimits(0,-265,-200);
//     func1->SetParLimits(1,-2.6,-2.4);
    
    TF1 *func2 =new TF1("func2",lastexpo, 89, 100, 2);
    func2->SetParameters(150,-2.5);
//     func2->SetParLimits(0,1000,600);
//     func2->SetParLimits(1,-2.6,-2.4);
   
TFile* file = new TFile("ROOT-m20_pos_22.5deg_y200_0.2mm_SW0.11mm_14.12.18-SKAN1.root"); 
     
TGraph *fHis = (TGraph*)file->FindObjectAny("GRAPH1");

    TF1 *func =new TF1("func",myfun, 30,50,3);
    func->SetParameters(750,35,0.001);
    func->SetParLimits(0,1,73600);
    func->SetParLimits(1,1,1001);
    func->SetParLimits(2,0.8,1);
    
/*    TF1 *curve =new TF1("curve",Ellip, 81,88,3);
    curve->SetParameters(73500,75,5);
    curve->SetParLimits(0,1,73600);
    curve->SetParLimits(1,50,100);
    curve->SetParLimits(2,4,10);*/    
    
//         TF1 *func =new TF1("func",myfun, 81, 88,3);
//     func->SetParameters(73500,85,0.6);
//     func->SetParLimits(0,1,73600);
//     func->SetParLimits(1,81,90);
//     func->SetParLimits(2,0.1,1);
    
    //func->5(3,1,20);
    //TF1 *total = new TF1("total","func1+func+func2",60,120,7);
//    total->SetParameters(10,2,73600,85,0.6,93,-2.5);
/*    total->SetParLimits(0,72.1, 81.3);
    total->SetParLimits(1,1.9,2.1);
    //total->SetParameters(73600,85,0.6);
    total->SetParLimits(2,1,73600);
    total->SetParLimits(3,81,90);
    total->SetParLimits(4,0.1,1);
    total->SetParLimits(5,88.2, 98);
    total->SetParLimits(6,-2.6,-2.4);*/ 
//     TF1 *func =new TF1("func",myfun, 170, 172.8,3);
//     func->SetParameters(61000,172,0.6);
//     func->SetParLimits(0,1,63000);
//     func->SetParLimits(1,163,175);
//     func->SetParLimits(2,0.1,1);
//     TF1 *bw = new TF1("bw","[0]/((x-[2])*(x-[3])+[1]/4)*800",0,200);
// //bw->SetParameters(700000, 0.5, 10000);
//     bw->SetParLimits(0,1,200);
//     bw->SetParLimits(1,1,200);
//     bw->SetParLimits(2,1,200);
//     bw->SetParLimits(3,1,200);
//     
//     //TF1 *func =new TF1("func","conf_hyperg", 1, 200,3);
//     TF1 *bw = new TF1("bw","[0]/((x-[2])*(x-[3])+[1]/4)*800",0,75000);
// //     bw->SetParameters(750000, 0.5, 10000);
//     bw->SetParLimits(0,1,200);
//     bw->SetParLimits(1,1,200);
//     bw->SetParLimits(2,1,200);
//     bw->SetParLimits(3,1,200);
//     
    //func->SetParNames("Minor-axis","H-center","Eccentricity","V-center");
    //(TMath::Sqrt((par[0]*par[0])-((X - par[1])*(X - par[1]))/((1-par[2])*(1-par[2]))));
//     TF1 *bw = new TF1("bw","TMath::Sqrt(([0]*[1])-((x-[2])*(x-[3])))/((1-[4])*(1-[5]))",10,200);
//     //bw->SetParameters(73500,90,0.6);
//     bw->SetParLimits(0,1,73600);
//     bw->SetParLimits(1,1,73600);
//     bw->SetParLimits(2,1,200);
//     bw->SetParLimits(3,1,200);
//     bw->SetParLimits(4,0.1,1);
//     bw->SetParLimits(5,0.2,2);
    c1->cd();
//       curve->Draw();
    fHis->Draw("AP");
//     fHis->Fit(bw,"PR");
//     fHis->Fit(curve,"R");
// TF1 *f1 = new TF1("f1", "sin(x)", 1, 3.14*20);
// fHis->Fit("f1");
    
//     TF1 *f1 = new TF1("f1","[0]*TMath::Erf((x-[1])/[2])*75000",-0.00005,0.00005);
//     fHis->Fit("f1");
    
       // fHis->Fit(func1,"R+");
       // fHis->Fit(func2,"R+");
        fHis->Fit(func,"R+");
        
        TF1 *fitresult = fHis->GetFunction("func");
        cout << "max value is: " << fitresult->GetMaximumX() << endl;
    
        Double_t p0=func->GetParameter(0);
        Double_t p1=func->GetParameter(1);
        Double_t p2=func->GetParameter(2);
        cout <<"######----Value of p0: " << p0 << "\n" <<"######----Value of p1: " << p1 << "\n" <<"######----Value of p2: " << p2 << endl;  
        
         TF1 *total = new TF1("total",totalfun,0,200,7);
        total->SetParameters(-110,2.1,7350,70,-100000,140, - 2.41);
//         total->SetParameters(-1100,50,1500,-3,735000,2,1);
//          Double_t par[7];
//          func1->GetParameters(&par[0]);
//          func->GetParameters(&par[2]);
//          func2->Get    gr->Fit(bw,"PR");
// Parameters(&par[5]);
//          total->SetParameters(par);


        //2fHis->Fit(total,"R");
//            fHis->Fit("total");

        
//     TF1 *total1 = new TF1("total1",myfun1,1,200,3);
//       total1->SetParameters(7, 800, 10,20,20);
//       total1->SetParLimits(0,1,200);
//       total1->SetParLimits(1,1,200);
//       total1->SetParLimits(2,1,200);
//       total1->SetParLimits(3,1,200);
//       total1->SetParLimits(4,1,200);
//       total1->SetParLimits(5,1,200);
// 
//       fHis->Fit(total1,"PR");

    
    
//     TF1 *Jerry = new TF1("Jerry",Straw,78,91.9,2);
//     Jerry->SetParameters(5.05,85);
    //Jerry->SetParLimits(0,1,10);
//     Jerry->SetParLimits(1,79.8,89.5);
//     Jerry->SetParLimits(0,1,10);
//     bw->SetParLimits(3,1,200);
//     fHis->Fit(Jerry,"PR");

        
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// TF1Convolution *f_conv = new TF1Convolution("firstexpo","myfun","lastexpo",60,100,true);
//    f_conv->SetRange(60.,100.);
//    f_conv->SetNofPointsFFT(201);
//    TF1   *f = new TF1("f",*f_conv, 61., 99., f_conv->GetNpar());
//    //f->SetParameters(1.,-0.3,0.,1.);        
//    fHis->Fit("f");
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    /*
    int num;
    double numx,numy;
    Double_t x, y;
    num = fHis->GetN();
    numx = fHis->GetX();
    numy = fHis->GetY();
    cout<<"numberx "<<numx<<endl;
    for(int i=0; i<num; i++)
    cout << GetPoint(i,x,y) << endl;
    
    Double_t ax[201],ay[201], xpar[201], ypar[201];
    Int_t n;
    double center;
    n=fHis->GetN(); //get ploted array dimention

    cout<<n<<endl;

    for(Int_t i=0; i<n; i++) {
        fHis->GetPoint(i,ax[i],ay[i]); 
            //cout<<i<<"th element of X array: "<<ax[i]<<endl;


        if(ay[i] > 1000 )
        {
            fHis->GetPoint(i,xpar[i],ypar[i]); 
            ypar[i] = ay[i];
            cout<<ax[i]<<"  ------------------------   "<<ay[i]<<" ------ypar------- "<< ypar[i] << " ------xpar------- "<< xpar[i] <<endl;
        }

    }*/


//     TGraph *gr = new TGraph(201,ax,ay);
//     c2->cd();
//     gr->Draw();
//     gr->Fit(Jerry,"R");

    TF1 *bw = new TF1("bw","[0]/((x-[2])*(x-[3])+[1]/4)*800",0,75000);
    // bw->SetParameters(750000, 0.5, 10000);
    bw->SetParLimits(0,1,200);
    bw->SetParLimits(1,1,200);
    bw->SetParLimits(2,1,200);
    bw->SetParLimits(3,1,200);
//     gr->Fit(bw,"PR");
/*
    gr->Fit(func1,"PR+");
    gr->Fit(func2,"PR+");*/

    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    //fHis->Fit(func,"PR");
//     TF1 *fitresult = fHis->GetFunction("func");
//     cout << "max value is: " << fitresult->GetMaximumX() << endl;

    
}
// 
