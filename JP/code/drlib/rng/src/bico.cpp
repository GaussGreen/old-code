
#include "edginc/coreConfig.hpp"

#define _NRC2_SOURCE

#include <cmath>

CORE_BEGIN_NAMESPACE

using namespace std;

RNG_DLL double SC_factln(int n);

RNG_DLL double  SC_bico(int n, int k)
{

    return floor(0.5+exp(SC_factln(n)-SC_factln(k)-SC_factln(n-k)));
}
/* (C) Copr. 1986-92 Numerical Recipes Software *0-12'=. */
CORE_END_NAMESPACE
