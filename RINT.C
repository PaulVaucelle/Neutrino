{
    /**Units**/

    #include <TMath.h>
    #include <TCanvas.h>
    #include <TStyle.h>
    #include "stdio.h"
    #include <TGraph.h>
    #include <TF1.h>
    
    #define m 1000
    #define Nbcore 10

    double Emin=1.81;/*MeV a little bit less in reality, this is an approximation*/
    double Emax=10; /*MeV as the energy of neutrinos coming from NPP reaches only around 10MeV*/
    double baseline[Nbcore]={52750,52840,52420,52510,52120,52210,52760,52630,52320,52200}; /*m could be one of the future entry parameter of our code*/
    double Gth[Nbcore]={2.9e9,2.9e9,2.9e9,2.9e9,2.9e9,2.9e9,4.6e9,4.6e9,4.6e9,4.6e9};

    /*Cross-sections parameters*/

    double m_e=0.511;/*positron's mass MeV*/
    double D=1.29;/*Delta Mn-Mp 1.29 MeV*/

    /*Variables*/
    
    double pas = (Emax-Emin)/m;/*pas*/ 
    double E[m];/*Neutrino energy*/
    double X[m];/* Xsection at  0th order*/
    double Phi[m];

    double Itot=0.;
    double Sumfiei = 0.564*202.36+205.99*0.076+211.12*0.304+214.26*0.056;/*Needed for the cimputation of 13.2*/
    double fi235U=0.564;
    double fi238U=0.076;
    double fi239Pu=0.304;
    double fi241Pu=0.056;
    double epsilon = 1;
    double SiX[m];
    double Const[Nbcore];

    int TimeConverter = 60*60*24;
    /*COmputation of the number of protons*/
    double Tmass= 20000000000;/*g*/
    double wc = 0.8792;
    double wh = 0.1201;
    double Mc = 12;/*g/mol*/
    double Mh = 1;/*g/mol*/
    double mc = wc*Tmass;
    double mh = wh*Tmass;
    double nc = (mc/Mc);
    double nh = (mh/Mh);
    double Na = 6.022e23;
    double nbrC = nc*Na;
    double nbrH = nh*Na;
    double Np = nbrH+6*nbrC;
    /*-------------------------------------*/
    TF1 *fa3 = new TF1("fa3","[0]*exp(0.870-0.160*x-0.091*x*x)*0.0952e-42*sqrt((x-[5])*(x-[5])-[6]*[6])*(x-[5])+[1]*exp(0.896-0.239*x-0.0981*x*x)*0.0952e-42*sqrt((x-[5])*(x-[5])-[6]*[6])*(x-[5])+[3]*exp(0.976-0.162*x-0.0790*x*x)*0.0952e-42*sqrt((x-[5])*(x-[5])-[6]*[6])*(x-[5])+[4]*exp(0.793-0.080*x-0.1085*x*x)*0.0952e-42*sqrt((x-[5])*(x-[5])-[6]*[6])*(x-[5])",Emin,Emax);
    /*TF1 *fa3 = new TF1("fa3","fi235U[0]*exp(0.870-0.160*E-0.091*E*E)*0.0952e-42*sqrt((E-D)*(E-D)-m_e*m_e)*(E-D)+fi239[1]Pu*exp(0.896-0.239*E-0.0981*E*E)*0.0952e-42*sqrt((E-D)*(E-D)-m_e*m_e)*(E-D)+fi238U[3]*exp(0.976-0.162*E-0.0790*E*E)*0.0952e-42*sqrt((E-D)*(E-D)-m_e*m_e)*(E-D)+fi241Pu[4]*exp(0.793-0.080*E-0.1085*E*E)*0.0952e-42*sqrt((E-D)*(E-D)-m_e*m_e)*(E-D)",Emin,Emax);*/
    fa3->SetParameter(0,fi235U);
    fa3->SetParameter(1,fi239Pu);
    fa3->SetParameter(3,fi238U);
    fa3->SetParameter(4,fi241Pu);
    fa3->SetParameter(5,D);
    fa3->SetParameter(6,m_e);
    fa3->SetTitle("function to integrate wrt the neutrino energy");
    fa3->GetXaxis()->SetTitle("Energy of the Neutrino (MeV)");
    fa3->GetYaxis()->SetTitle("Integrated function");
    /*fa3->Draw();*/
    E[0]=Emin;

    
    for (int i=0;i<m;i++)/*sums over bins*/         
        {
        /**Const term computation*/
        for (int r=0;r<Nbcore;r++)/*sum over the cores*/
            {
            Const[r]= Gth[r]*1/(Sumfiei)*TimeConverter*6.42*pow(10,+12)/(4*4*atan(1.0)*baseline[r]*baseline[r]*(10^4))*epsilon*Np;/*Units correction L in cm*/

            /*For the Phi, the fission fraction are respectively: 235U, 239Pu, 238U and 241Pu*/
            E[i+1]=E[i]+pas;/*Neutrino energy incresed by a bin*/
            /*Simpson rule to compute the integral 1/3*/

            
            SiX[i]= fa3->Integral(E[i],E[i+1],1e-12);

            Itot=Itot+Const[r]*SiX[i];/* sum all over the bins*//*13.5*//*Doesn't take into accout the oscillations*/
            }
        }

    double Ip= fa3->Integral(E[0],E[1],1e-12);
    double Ipc= Ip*Const[1];
    printf("Ipc = %f \n", Ipc);
    /*TCanvas *c3 = new TCanvas("c3","Integrated function",200,10,600,400);*/
    TGraph *gr6 = new TGraph(m,E,SiX);/*Si on trace I, on retrouve la forme caluler pour F1 dans le code RFluxocs.C*/
    gr6->SetLineColor(6);
    gr6->SetTitle("Integrated function");
    gr6->GetXaxis()->SetTitle("Energy of the Neutrino (MeV)");
    gr6->GetYaxis()->SetTitle("Integrated function");
    /*gr6->Draw("AL");*/
    printf("Itot = %f \n",Itot);
    

}
