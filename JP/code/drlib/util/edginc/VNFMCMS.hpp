//----------------------------------------------------------------------------
//
//   Group       : Credit QR
//
//   Filename    : VNFMCMS.hpp
//
//   Description : Class for VNFMCMS approximation
//
//   Author      : Keith Jia
//
//
//----------------------------------------------------------------------------

#ifndef _VNFMCMS_HPP
#define _VNFMCMS_HPP

#include "edginc/VNFM.hpp"

DRLIB_BEGIN_NAMESPACE

class UTIL_DLL VNFMCMS: public virtual VNFM,
                        public CObject
{
public:
    static CClassConstSP const TYPE;
    virtual ~VNFMCMS();

    VNFMCMS(double dVolStart,
            double dexpiry,
            double dSwapStart,
            double dSwapMat,
            int    swapFreq,
            double beta,
            bool   isCashSettle,
            double ParFwd,
            double sigATM,
            double fwdAnnuity);

    /** prepare for measure change data */
    virtual void Q3VNFMZero2Swap(double mat, 
                                 double zeroRate,
                                 double& alpha,
                                 double& mpower,
                                 double& zMat);

protected:
    VNFMCMS(CClassConstSP clazz = TYPE);

protected:
    double dVolStart;
    double dexpiry;
    double dSwapStart;
    double dSwapMat;
    int    swapFreq;
    double beta;
    bool   isCashSettle;
    double parFwd;
    double sigATM;
    double fwdAnnuity;

private:
    static void load(CClassSP& clazz);

};

DECLARE(VNFMCMS);

DRLIB_END_NAMESPACE

#endif
