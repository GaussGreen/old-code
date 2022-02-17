//----------------------------------------------------------------------------
//
//   Group       : Credit QR
//
//   Filename    : VNFMCMS.cpp
//
//   Description : class for CMS VNFM approximation
//
//   Author      : Keith Jia
//
//
//----------------------------------------------------------------------------


#include "edginc/config.hpp"
#include "edginc/Class.hpp"
#include "edginc/VNFMCMS.hpp"

DRLIB_BEGIN_NAMESPACE

#define TINY 1.0e-10
#define BOND_PRICE(c,f,m,y) \
         (fabs(y)<TINY) ? (1 + (m) * (c)) : ((c) / (y) + \
         (1- (c) / (y)) * pow(1 + (y)/(f), -(m) * (f))) 

VNFMCMS::VNFMCMS(CClassConstSP clazz) : CObject(clazz)
{}

VNFMCMS::~VNFMCMS() 
{}

VNFMCMS::VNFMCMS(double dVolStart,
                 double dexpiry,
                 double dSwapStart,
                 double dSwapMat,
                 int    swapFreq,
                 double beta,
                 bool   isCashSettle,
                 double parFwd,
                 double sigATM,
                 double fwdAnnuity) :
    CObject(TYPE), dVolStart(dVolStart), dexpiry(dexpiry), 
    dSwapStart(dSwapStart), dSwapMat(dSwapMat), swapFreq(swapFreq),
    beta(beta), isCashSettle(isCashSettle), parFwd(parFwd),
    sigATM(sigATM), fwdAnnuity(fwdAnnuity)
{}


/** prepare for measure change data */
void VNFMCMS::Q3VNFMZero2Swap(double xZeroMat, 
                              double xZeroRate,
                              double& alpha,
                              double& mpower,
                              double& zMat)
{
    static char routine[] = "VNFMCMS::Q3VNFMZero2Swap";

    zMat = xZeroMat;

    if(xZeroMat < TINY || isCashSettle)
    {
        mpower = 1;
        alpha = 1;
        return;
    }

    /* Adjust swapMat */
    if(swapFreq == 12 || swapFreq == 4 || swapFreq == 2 || swapFreq == 1)
    {
        dSwapMat= floor(swapFreq*dSwapMat + 0.5) / (double) swapFreq;
    }
    else if (swapFreq == 365 || swapFreq == 52 || swapFreq == 26)
    {
        dSwapMat = 1.0/(double) swapFreq;
    }
    else
    {
        throw ModelException(routine, "Invalid frequence" + swapFreq);
    }
    
    /* check swpMat, zeroMat */
    if (dSwapMat < TINY || xZeroMat < TINY) {
        throw ModelException(routine, "Invalid swap maturity or zero maturity");
    }


    /* setup B coefficients */
    double swapPlusBeta, bondPrice1, bondPrice2;
    
    /* swap rate + mean reversion, on the swap rate compounding basis 
     * extra hassle to ensure cvx adj = delay adj when appropriate */
    swapPlusBeta = swapFreq * ((1 + parFwd / swapFreq) * 
                               exp (beta / swapFreq) - 1.);
    
    /* intermediate step: bond calculations */
    bondPrice1 = BOND_PRICE(parFwd, swapFreq, dSwapMat, parFwd);
    bondPrice2 = BOND_PRICE(parFwd, swapFreq, dSwapMat, swapPlusBeta);
    
    /* B coefficient for swap rate */
    double swapBCoeff = (bondPrice1 - bondPrice2) / (fwdAnnuity * beta); 
    
    /* B coefficient for zero rate */
    double zeroBCoeff = (1 + xZeroRate / swapFreq) *
        (1 - exp( -beta * xZeroMat)) / (beta * xZeroMat);


    /* total covariance */
    double covarGauss, betaSum, betaDur;
        
    betaSum = 2 * beta;
    betaDur = (betaSum*(dexpiry - dVolStart)>TINY)?
                (1-exp(-betaSum*(dexpiry - dVolStart))) / 
                (betaSum*(dexpiry - dVolStart)): 1.;

    covarGauss = exp(-betaSum*(dSwapStart - dexpiry)) * betaDur;

    double sVol  = swapBCoeff * swapBCoeff * covarGauss;
    double zVol  = zeroBCoeff * zeroBCoeff * covarGauss;
    double coVar = swapBCoeff * zeroBCoeff * covarGauss;

    double swapVol = sqrt(sVol);
    double zeroVol = sqrt(zVol);
    double zeroSwapCorr = coVar / (swapVol*zeroVol);

    double volRatio = zeroVol / swapVol;
    mpower = zeroSwapCorr * volRatio;
    alpha = xZeroRate / pow(parFwd, mpower) * 
        exp(0.5 * sigATM * sigATM * dexpiry * mpower * (1-mpower));
}


// static 
CClassConstSP const VNFMCMS::TYPE = CClass::registerClassLoadMethod(
    "VNFMCMS", typeid(VNFMCMS), load);



void VNFMCMS::load(CClassSP& clazz)
{
    REGISTER(VNFMCMS, clazz);
    SUPERCLASS(CObject);
    IMPLEMENTS(VNFM);
}

bool VNFMCMSLinkIn() {
    return VNFMCMS::TYPE!=0;
}



DRLIB_END_NAMESPACE

