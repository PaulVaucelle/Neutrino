{
    /**Units**/
    gROOT->Reset();/*called for unnamed script to clean the variables*/
    #include <TMath.h>
    #include <TCanvas.h>
    #include <TStyle.h>
    #include "stdio.h"
    #include <TGraph.h>
    
    
    #define m 50
    #define Nbcore 10
    /*Needed parameters*/

   double Emin=2.0;/*MeV a little bit less in reality, this is an approximation*/
    double Emax=8.21; /*MeV as the energy of neutrinos coming from NPP reaches only around 10MeV*/
    double baseline[Nbcore]={52750,52840,52420,52510,52120,52210,52760,52630,52320,52200}; /*m could be one of the future entry parameter of our code*/
    double Gth[Nbcore]={2.9e9,2.9e9,2.9e9,2.9e9,2.9e9,2.9e9,4.6e9,4.6e9,4.6e9,4.6e9};
    double Mbaseline=(baseline[0]+baseline[1]+baseline[2]+baseline[3]+baseline[4]+baseline[5]+baseline[6]+baseline[7]+baseline[8]+baseline[9])/10;
    /*printf("%f \n",Mbaseline);*/
    /*Oscillation parameters Normal hierarchy Nu-fit 2021 */

    double dm_12 = 7.42e-5; /*eV2*/ 
    double dm_13 = 2.515e-3; /*eV2*//*check Nufit*/
    double dm_23 = dm_13-dm_12;/*eV2*/
    double s_12 = 0.304;
    double s_23 = 0.573; 
    double s_13 = 0.02220;
    double T12= 33.44*4*atan(1.0)/180;
    double T23= 49.2*4*atan(1.0)/180;
    double T13= 8.57*4*atan(1.0)/180;
    double s2_13 = sin(2*T13)*sin(2*T13) /*0.087*/;/*remove one of the sx_13 to avoid any issue*/
    double s2_12 = sin(2*T12)*sin(2*T12);

    /* could use a struct for the parametres*/

    /*Oscillation parameters Inverted Hierarchy Nu-fit 2021*/

    double Idm_12 = 7.42e-5; /*eV2*/
    double Idm_13 = -2.498e-3; /*eV2*//*issue here*/
    double Idm_23 = abs(Idm_13)+abs(Idm_12);/*eV2*/
    double Is_12 = 0.304;
    double Is_23 = 0.578;
    double Is_13 = 0.02238;
    double IT13 = 8.60*4*atan(1.0)/180;
    double IT12 = 33.45*4*atan(1.0)/180;
    double Is2_13 = sin(2*IT13)*sin(2*IT13) /*0.087*/;
    double Is2_12 = sin(2*IT12)*sin(2*IT12);
   
    /*Delta functions (contain the L/E dependance for Normal Hierarchy)*/
    double d31[m];
    double d21[m];
    double d32[m];

    /*Delta functions (contain the L/E dependance for Inverted Hierarchy)*/
    double Id31[m];
    double Id21[m];
    double Id32[m];

    /*Cross-sections parameters*/

    double m_e=0.511;/*positron's mass MeV*/
    double D=1.29;/*Delta Mn-Mp 1.27 MeV*/

    /*Variables*/

    double pas = (Emax-Emin)/m;/*pas*/ 
    double E[m];/*Neutrino energy*/
    double X[m];/* Xsection at  0th order*/
    double Phi[m];/*flux of neutrino as [65]*/
    double Eres[m];/*Energy resolution*/
    double ey[m];

        /*Normal Hierarchy*/

    double F1[m];/*F [65] without taking into account the probability*/
    double P21[m];
    double F12[m];/* 1-p21 oscillations*/
    double P31[m];
    double P32[m];
    double P321[m]; /*1-p21-p32-p31,*/
    double F321[m]; /*F [65] just using P321 [65] instead of the formula given in the paper of JP*/
    double F321res[m];

        /*Inverted Hierarchy*/

    double IP21[m];
    double IP31[m];
    double IP32[m];
    double IP321[m]; /*1-p21-p32-p31,*/
    double IF321[m]; 

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
    fa3->Draw();

    /*Proba computation*/
    E[0]=Emin;
    TH1F *hist = new TH1F("hist", "Histogram",m,2.00,8.21);
    TCanvas *c5 = new TCanvas();
    for (int i=0;i<m;i++)
        {
            
             /*delta functions L/E dependance for Normal Hierarchy*/

            d31 [i]= 1.27*dm_13*Mbaseline/(E[i]); /*where does the 1.27 come from?*/
            d21 [i]= 1.27*dm_12*Mbaseline/(E[i]);/* Keep dmxx in eV2 as requested per [65];*/
            d32 [i]= 1.27*dm_23*Mbaseline/(E[i]);

            /*delta functions L/E dependance for Inverted Hierarchy*/

            Id31 [i]= 1.27*Idm_13*Mbaseline/(E[i]);
            Id21 [i]= 1.27*Idm_12*Mbaseline/(E[i]);
            Id32 [i]= 1.27*Idm_23*Mbaseline/(E[i]);

            /*Normal Hierarchy*/

            P21[i]=(1-s_13)*(1-s_13)*s2_12*sin(d21[i])*sin(d21[i]);
            P31[i]=(1-s_12)*s2_13*sin(d31[i])*sin(d31[i]);
            P32[i]=s_12*s2_13*sin(d32[i])*sin(d32[i]);
            P321[i]=1-P21[i]-P31[i]-P32[i];

            /*Inverted Hierarchy*/

            IP21[i]=(1-Is_13)*(1-Is_13)*Is2_12*sin(Id21[i])*sin(Id21[i]);
            IP31[i]=(1-Is_12)*Is2_13*sin(Id31[i])*sin(Id31[i]);
            IP32[i]=Is_12*Is2_13*sin(Id32[i])*sin(Id32[i]);
            IP321[i]=1-IP21[i]-IP31[i]-IP32[i];

            /* Main Parameters */

            X[i]= 0.0952e-42*sqrt((E[i]-D)*(E[i]-D)-m_e*m_e)*(E[i]-D);/*Xsection*/
            Eres[i]=0.03/(sqrt(E[i]));
            E[i+1]=E[i]+pas;/*Neutrino energy*/

            /*flux of neutrinos*/

            Phi[i]=0.564*exp(0.870-0.160*E[i]-0.091*E[i]*E[i])+0.304*exp(0.896-0.239*E[i]-0.0981*E[i]*E[i])+0.076*exp(0.976-0.162*E[i]-0.0790*E[i]*E[i])+0.056*exp(0.793-0.080*E[i]-0.1085*E[i]*E[i]);
            /*For the Phi, the fission fraction are respectively: 235U, 239Pu, 238U and 241Pu*/

            F1[i]=X[i]*Phi[i];/*without oscillations P=1*/
            F321[i]=X[i]*Phi[i]*P321[i]; /* Taking into account small oscillation => Mass Ordering differentiation*/
            F12[i]=X[i]*Phi[i]*(1-P21[i]);/*No Mass Ordering part*/
            IF321[i]=X[i]*Phi[i]*IP321[i];/* Inverted Hierarchy*/
            hist->Fill(E[i],F321[i]);
            F321res[i]=F321[i]*Eres[i];
            ey[i]=0;
        }
    /*hist->Draw("HIST SAME C");*/
    // double scale= 1/hist->Integral();
    hist->GetXaxis()->SetTitle("Neutrino Energy MeV");
    hist->GetYaxis()->SetTitle("Neutrino Spectra");
    hist->SetFillColor(kBlue-7);
    hist->Draw("H");
 
    /*TGraph part*/

    /*Building the Baseline over E array => X axis*/

    double LER[m];
    double Evis[m];
    for (int j=0;j<m;j++)
        {   
            Evis[j]=E[j]-1;
            LER[j]=Mbaseline/E[j];
            
        }
        
        for (int b=0;b<m;b++)
        /*sums over bins*/         
            {for (int r=0;r<Nbcore;r++)
        /**Const term computation*/
            /*sum over the cores*/
                {
                Const[r]= Gth[r]*1/(Sumfiei)*TimeConverter*6.42*pow(10,+12)/(4*4*atan(1.0)*baseline[r]*baseline[r]*(1e4))*epsilon*Np;/*Units correction L in cm*/

                /*printf("Constr = %f \n", Const[r]);*/
                SiX[b]= fa3->Integral(E[b],E[b+1],1e-12);

                /*printf("E[b] = %f  E[b+1] = %f \n",E[b],E[b+1]);*/
                Itot=Itot+Const[r]*SiX[b]*P321[b];/* sum all over the bins*//*13.5*//*Doesn't take into accout the oscillations*/
                }
            }
 
    TCanvas *c3 = new TCanvas("c3","Integrated function",200,10,600,400);
    TGraph *gr6 = new TGraph(m,E,SiX);/*Si on trace I, on retrouve la forme caluler pour F1 dans le code RFluxocs.C*/
    gr6->SetLineColor(6);
    gr6->SetTitle("Integrated function");
    gr6->GetXaxis()->SetTitle("Energy of the Neutrino (MeV)");
    gr6->GetYaxis()->SetTitle("Integrated function");
    gr6->Draw("AL");
    printf("Itot = %f \n",Itot);
    
    /* Neutrino spectrum as a function of E*/
    
    TCanvas *c1 = new TCanvas("c1","First Graph",200,10,600,400);
    TGraph *gr1 = new TGraph(m,Evis,F321);/* don't use small bin number*/
    TGraph *gr8 = new TGraph(m,Evis,IF321);/* don't use small bin number*/
    gr1->SetLineColor(4);
    gr1->SetTitle("Observed neutrino spectrum with oscillations");
    gr1->GetXaxis()->SetTitle("Visible Energy (MeV)");
    gr1->GetYaxis()->SetTitle("Neutrino flux * Xsection * Probability (cm2/MeV/fisson)");
    
    gr1->Draw("AL");
    gr8->Draw("L"); /* ACTIVATE/DESACVTIVATE*/

    /*Neutrino spectrum as a function of L/E*/

    /*Building a graph for each F function*/
    TCanvas *c2 = new TCanvas("c2","Second Graph",200,10,600,400);

    TGraph *gr2 = new TGraph(m,LER,F1);/* don't use small bin number*/
    TGraph *gr3 = new TGraph(m,LER,F321);/* don't use small bin number*/
    TGraph *gr4 = new TGraph(m,LER,F12);/* don't use small bin number*/
    TGraph *gr5 = new TGraph(m,LER,IF321);/* don't use small bin number*/

    gr5->SetLineColor(6);/* 5 is for yellow*/
    gr4->SetLineColor(4);/* 4 is for blue*/
    gr2->SetLineColor(2);/*2 is for red*/
    gr3->SetLineColor(3);/* 3 is for green*/

    /***MultiGraph***/

    TMultiGraph *mg = new TMultiGraph();
    mg->SetTitle("Observed neutrino spectrum with oscillations");
    mg->Add(gr2);
    mg->Add(gr3);
    mg->Add(gr4);
    mg->Add(gr5);

    mg->Draw("AL");/* ACTIVATE/DESACVTIVATE*/

    mg->GetXaxis()->SetTitle("L/E (Meters/MeV)");
    mg->GetYaxis()->SetTitle("Neutrino flux * Xsection * Probability (cm2/MeV/fisson)");

    /*Legend*/

    auto legend = new TLegend(0.9,0.75,0.7,0.9);/*(gap from left,Gap from the bottom,percentage width filled,gap fro mtop?)*/
    legend->SetHeader("Legend","C"); /*C to center the title*/
    legend->AddEntry(gr2,"No oscillation Pee=1");
    legend->AddEntry(gr4,"1-P12");
    legend->AddEntry(gr3,"NH Pee");    
    legend->AddEntry(gr5,"IH Pee");
    legend->Draw();/* ACTIVATE/DESACVTIVATE*/

    TCanvas *c4 = new TCanvas("c4","Fourth Graph",200,10,600,400);
    auto gr = new TGraph(m,Evis,F321);/*use small bin number defined by the constaint of nergy resolution*/
    TGraph *gr9 = new TGraph(m,Evis,IF321);
    gr->SetTitle("Observed Neutrino Spectrum with oscillations");
    gr->GetXaxis()->SetTitle("Visible Energy(MeV)");
    gr->GetYaxis()->SetTitle("Observed Neutrino Spectrum with oscillations");
    gr->SetLineColor(2);
    gr->Draw("AL");
    gr9->Draw("L");

    auto legend2 = new TLegend(0.9,0.75,0.7,0.9);/*(gap from left,Gap from the bottom,percentage width filled,gap fro mtop?)*/
    legend2->SetHeader("Legend","C"); /*C to center the title*/
    legend2->AddEntry(gr,"NH with Energy resolution");
    legend2->AddEntry(gr9,"IH with energy resolution");
    legend2->Draw();
}