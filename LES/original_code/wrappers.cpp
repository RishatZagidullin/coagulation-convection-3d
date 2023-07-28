#include "wrappers.h"
#include <cmath>


namespace wrappers
{
    double K(const int & u, const int &v, const double h)
    {
        //simplified kernel for res_full.txt
        //double u1=pow( (u + 1.0) , 2.0/3.0);
        //double v1=pow( (v + 1.0) , 2.0/3.0);
        //double result = u1*v1;

        //ballistic kernel for res_ballistic.txt
        double u1=pow( (u + 1.0) , 1.0/3.0);
        double v1=pow( (v + 1.0) , 1.0/3.0);
        //double result = (u1+v1)*(1./u1+1./v1);
        
        double result = 0.1;
        return result;
    }
}

