//----------------------------------------------------------------------------
//
//   Group       : Derivatives Research
//
//   Filename    : LossDistribution.hpp
//
//   Description : Convolution Algorithm
//
//   Date        : Fed 2006
//
//
//----------------------------------------------------------------------------
#ifndef EDR_LOSSDISTRIBUTION_HPP
#define EDR_LOSSDISTRIBUTION_HPP

#include "edginc/AtomicArray.hpp"

DRLIB_BEGIN_NAMESPACE

class CONVOLUTION_DLL LossDistribution : public CObject
{
public:
    enum DistributionType { TOPLOSS, BOTTOMLOSS };

    //constructors
    LossDistribution();
    LossDistribution(const int nbLoss, const int maxShortIdx);

    //add a name to the loss distribution
    //- if for the first time, adjustBounds should be true
    //- if not (ie the name has been previously de-convoluted), adjust bounds should be false
    void convolute( const int    nmIdx,
                    const double s,
                    const double weight1,
                    const int    l1,
                    const int    l2);

    //add a distribution (not just a name) to the loss distribution
    void convoluteDistribution(int numPoints,
                               const double * s,
                               const double * l);

    //remove a name from the loss distribution
    //- used for efficient tweaking of CDO portfolios
    //- the untweaked name is first removed, and then the tweaked name is added
    void deconvolute( const int    nmIdx,
                      const double s,
                      const double weight1,
                      const int    l1,
                      const int    l2);

    int getLastRunName();

    DoubleArray& getLossDistribution();

    //----------------
    // CObject methods
    //----------------
    static CClassConstSP const TYPE;
    static void load(CClassSP& clazz);
    static IObject* defaultLossDistribution();

private:

    double maxAtTiny(double); //used to rationalise precision of the loss distributions

    DoubleArray     lossDist;   //the distribution
    int             left;       //how far we have progressed along the distribution
    int             right;      //how far we have progressed along the distribution
    int             zero;       //the position of the "0" loss element cf maxShoftIdx
    int             lastRunName; //which names contribute to this distribution
};

typedef smartPtr<LossDistribution>                  LossDistributionSP;
typedef array<LossDistributionSP, LossDistribution> LossDistributionArray; //note array of smart pointers

DRLIB_END_NAMESPACE
#endif
