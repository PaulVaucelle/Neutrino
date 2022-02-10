
#include "func.h"

 double funX(double E) 
    {
        double D=1.29;
        double m_e=0.511;
        return 0.0952e-42*sqrt((E-D)*(E-D)-m_e*m_e)*(E-D); 
    }

    double funcPhi(double E) 
    {
        double fi235U=0.564;
        double fi238U=0.076;
        double fi239Pu=0.304;
        double fi241Pu=0.056;
        return fi235U*exp(0.870-0.160*E-0.091*E*E)*funX(E)+fi239Pu*exp(0.896-0.239*E-0.0981*E*E)*funX(E)+fi238U*exp(0.976-0.162*E-0.0790*E*E)*funX(E)+fi241Pu*exp(0.793-0.080*E-0.1085*E*E)*funX(E);
    }

int main()
{
    return 0;
}