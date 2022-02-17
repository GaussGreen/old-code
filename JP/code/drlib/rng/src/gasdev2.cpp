#include "edginc/gasdev2.h"
#include <cmath>

CORE_BEGIN_NAMESPACE

using namespace std;

GasDev2::GasDev2(IUniformRNGGenSP rng, bool _iset, double _gset) :
        uniform(rng), iset(_iset), gset(_gset)
{
    assert(rng.get());
}

GasDev2SP GasDev2::create (IUniformRNGGenSP rng, bool _iset, double _gset)
{
    return GasDev2SP(new GasDev2(rng, _iset, _gset));
}

/* (C) Copr. 1986-92 Numerical Recipes Software *0-12'=. */

double GasDev2::fetch()
{
    double fac,rsq,v1,v2;

    if  (!iset) {
        do {
            v1=2.0*uniform->fetch()-1.0;
            v2=2.0*uniform->fetch()-1.0;
            rsq=v1*v1+v2*v2;
        } while (rsq >= 1.0 || rsq == 0.0);
        fac=sqrt(-2.0*log(rsq)/rsq);
        gset=v1*fac;
        iset=true;
        return v2*fac;
    } else {
        iset=false;
        return gset;
    }
}


SC_gasdev2::SC_gasdev2(ISuperCubeRNGGenSP rng, bool _iset, double _gset) :
        GasDev2(rng, _iset, _gset),
        hooks(rng.get())
{}

SC_gasdev2SP    SC_gasdev2::create(ISuperCubeRNGGenSP rng, bool _iset, double _gset)
{
    return SC_gasdev2SP(new SC_gasdev2(rng, _iset, _gset));
}

/* iset is set to zero by default. If SC_init_gasdev2 is not
called then SC_gasdev2() operates as normal. Note that to
reset SC_ran2 a negative seed MUST be passed in. */

void SC_gasdev2::init()
{
    iset = false;
}

IRNGSP SC_gasdev2::clone() const
{
    return SC_gasdev2SP(new SC_gasdev2(DYNAMIC_POINTER_CAST<ISuperCubeRNGGen>(getGenerator()->clone()), iset, gset));
}

CORE_END_NAMESPACE
