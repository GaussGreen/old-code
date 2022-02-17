//----------------------------------------------------------------------------
//
//   Group       : Credit QR
//
//   Filename    : VNFMCMCDS.hpp
//
//   Description : Class for VNFMCMCDS approximation
//
//   Author      : Mehdi Chaabouni
//
//
//----------------------------------------------------------------------------

#ifndef _VNFMCMCDS_HPP
#define _VNFMCMCDS_HPP

#include "edginc/VNFM.hpp"

DRLIB_BEGIN_NAMESPACE

class UTIL_DLL VNFMCMCDS: public virtual VNFM,
                        public CObject
{
public:
    static CClassConstSP const TYPE;
    virtual ~VNFMCMCDS();

    VNFMCMCDS(double dVolStart,
            double dexpiry,
            double dSwapStart,
            double dSwapMat,
            int    swapFreq,
            double beta,
            bool   isCashSettle,
            double ParFwd,
            double sigATM,
            double fwdAnnuity,
			double fwdSpread,
			double recovery);

    /** prepare for measure change data */
	// Function calculating the zero and the swap vol and correlation
  // Then it calculates the power and alpha coefficients
  // If physically settled uses VNFM to to find them
  // if cash settled alpha = power = 1
  // then  E( Zero | Spread) = alpha * Spread^power (used to adjust density)
    virtual void Q3VNFMZero2Swap(double mat, 
                                 double zeroRate,
                                 double& alpha,
                                 double& mpower,
                                 double& zMat);

protected:
    VNFMCMCDS(CClassConstSP clazz = TYPE);

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
	double fwdSpread;
	double recovery;

private:
    static void load(CClassSP& clazz);

};

DECLARE(VNFMCMCDS);

DRLIB_END_NAMESPACE

#endif
