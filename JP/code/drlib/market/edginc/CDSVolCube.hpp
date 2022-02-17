//----------------------------------------------------------------------------
//
//   Group       : Credit QR
//
//   Filename    : CDSVolCube.hpp
//
//   Description : Trivial vol cube interface for CDS European option type 
//                 instruments.
//                 The concrete vol cube/smile model classes are 
//                 CDSVolCubeBSImplied and CDSVolCubeMultiQ
//
//   Author      : Charles Morcom
//
//   Date        : 5 January 2006
//
//----------------------------------------------------------------------------
#ifndef QR_ICDSVOLCUBE_HPP
#define QR_ICDSVOLCUBE_HPP

#include "edginc/FORWARD_DECLARE_WRAPPER.hpp"
#include "edginc/ICDSVol.hpp"

DRLIB_BEGIN_NAMESPACE

FORWARD_DECLARE_WRAPPER(ICDSVolCube)

/**Abstract class defining CDS option volatility cube. Note that this doesn't
 * extend CDSVolATMMatrix - that would usually be included as a member in the
 * concrete classes implementing this, but I didn't want to require this in
 * case one wanted to define a vol cube which is based on slices not including
 * ATM, or on some completely different parametric representation. */
class MARKET_DLL ICDSVolCube : public virtual ICDSVol {
public:
    static CClassConstSP const TYPE; // in CDSParSpreads.cpp

    virtual ~ICDSVolCube();
    /**Strikes are absolute numbers*/
    static const string STRIKE_TYPE_ABSOLUTE;
    /**Strikes are defined as a simple ratio to ATM*/
    static const string STRIKE_TYPE_RATIO;
    /**Strikes are defined by their call deltas for ATM vol*/
    static const string STRIKE_TYPE_DELTA;

    /**Smile vols are absolute numbers.*/
    static const string VOL_TYPE_ABSOLUTE;
    /**Smile vols should be added to ATM vols.*/
    static const string VOL_TYPE_ADDITIVE_OFFSET;
    /**Smile vols are ATM vols times the vol number.*/
    static const string VOL_TYPE_MULTIPLICATIVE_OFFSET;
private:
    static void load(CClassSP& clazz);
};



DRLIB_END_NAMESPACE
#endif

