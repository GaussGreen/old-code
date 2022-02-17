//----------------------------------------------------------------------------
//
//   Group       : Credit QR
//
//   Filename    : VNFM.hpp
//
//   Description : Base class for VNFM approximation
//
//   Author      : Keith Jia
//
//
//----------------------------------------------------------------------------

#ifndef _VNFM_HPP
#define _VNFM_HPP

#include "edginc/Object.hpp"
#include "edginc/DECLARE.hpp"

DRLIB_BEGIN_NAMESPACE


class UTIL_DLL VNFM : public virtual IObject 
{
public:
    static CClassConstSP const TYPE;
    virtual ~VNFM()= 0;

    /** prepare for measure change data */
    virtual void       Q3VNFMZero2Swap(double xZeroMat,    /* (I) */
                                       double xZeroRate,   /* (I) */
                                       double& alpha,      /* (O) */
                                       double& power,      /* (O) */
                                       double& zMat) = 0;  /* (O) */

protected:
    VNFM(CClassConstSP clazz = TYPE);

private:
    static void load(CClassSP& clazz);
};

DECLARE(VNFM);

DRLIB_END_NAMESPACE

#endif
