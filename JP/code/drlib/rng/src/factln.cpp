#include "edginc/coreConfig.hpp"

#define _NRC2_SOURCE
#include <cstdio>
using namespace std;

CORE_BEGIN_NAMESPACE

RNG_DLL double SC_gammln(double xx);

RNG_DLL double SC_factln(int n)
{
    /* void nrerror(char error_text[]); */
    static double a[101];

    if (n < 0)
        fprintf(stderr,"Negative factorial in routine SC_factln");
    if (n <= 1)
        return 0.0;
    if (n <= 100)
        return a[n] ? a[n] : (a[n]=SC_gammln(n+1.0));
    else
        return SC_gammln(n+1.0);
}
/* (C) Copr. 1986-92 Numerical Recipes Software *0-12'=. */
CORE_END_NAMESPACE
