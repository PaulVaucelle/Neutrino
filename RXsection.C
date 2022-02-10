
{
    gROOT->Reset();/*called for unnamed script to clean the variables*/
    #include <TMath.h>
    #include <TCanvas.h>
    #include <TStyle.h>
    #include "stdio.h"
    #include <TGraph.h>
    
    #define n 10

    /*Needed parameters*/
    float Emin=1.8; /*MeV*/
    float Emax=10.0; /*MeV*/
    /*int n = 10000; Define energy resolution*/
    float pas = Emax/n ;/*pas*/ /* EMax-Emin*/
    float E[n],X[n]; /*Energy, Xsection*/
    float m_e=0.511;/*positron's mass MeV*/
    float D=0.8;/*Delta Mn-Mp MeV*/
    /*X-section computation*/
    E[0]=Emin;
    for (int i=0;i<n;i++)
        {
            X[i]= 0.0952e-42*sqrt((E[i]-D)*(E[i]-D)-m_e*m_e)*(E[i]-D);
            printf("%f \n",X[i]);
            E[i+1]=E[i]+pas;
            printf("%f \n",E[i]);
        }

    /*TGraph part*/
    TGraph *gr1 = new TGraph(n,E,X);
    TCanvas *c1 = new TCanvas("c1","Graph Draw options",200,10,600,400);
    gr1->Draw("ALP");
    gr1->GetXaxis()->SetTitle("Electron Neutrino energy in MeV");
    gr1->GetYaxis()->SetTitle("Cross-section in cm2");
    /*gr1->GetXaxis()->CenterTitle();
    gr1->GetYaxis()->CenterTitle();*/
    /* Library issues*/
    gr1->Draw("ALP");

    auto legend = new TLegend(0.1,0.75,0.35,0.9);/*(gap from left,Gap from the bottom,percentage width filled,gap fro mtop?)
    legend->SetHeader("Legend","C"); /*C to center the title*/
    legend->AddEntry("gr1","Cross section as a function of the neutrino Energy");
    legend->Draw();   
}