//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : VolProcessedDVFParam.cpp
//
//   Description : Processed DVF parameterised vols.
//                 Previously, was part of VolParam.cpp
//
//   Author      : Regis Guichard
//
//   Date        : 02 Mai 2001
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/VolProcessedDVFParam.hpp"
#include "edginc/Maths.hpp"

DRLIB_BEGIN_NAMESPACE

const double CVolProcessedDVFParam::SMALL_YEARS_TO_MATURITY = 0.1 / 365.0;
const double CVolProcessedDVFParam::MIN_STRIKE_TWEAK = 0.0001;
const double CVolProcessedDVFParam::MIN_TIME_TWEAK = 0.0001;
const double CVolProcessedDVFParam::MIN_NEXT_DATE = 0.;
const double FP_MIN = 1.0e-10;

// CVolProcessedDVFParam stuff -----------------------------------------------

CVolProcessedDVFParam::~CVolProcessedDVFParam(){}


CClassConstSP const CVolProcessedDVFParam::TYPE =
CClass::registerClassLoadMethod(
    "VolProcessedDVFParam", typeid(CVolProcessedDVFParam), load);

void CVolProcessedDVFParam::load(CClassSP& clazz){
    REGISTER(CVolProcessedDVFParam, clazz);
    SUPERCLASS(CVolProcessedDVF);
}

/////////////////////////
void CVolProcessedDVFParam::Init(   const CVolBase*         vol,
                                    const CVolParamConstSP& volParam,
                                    const CVolRequestDVF* volRequest,
                                    const CAsset* asset,
                                    const VolSurface* VolSurf,
                                    const FXAsset*   fxAsset,
                                    const  Correlation* eqFXCorr)
{
    myVol = CVolBaseConstSP::attachToRef(vol);
    myParamVol = volParam;
    this->asset = AssetConstSP::attachToRef(asset);
    myVolSurf = VolSurfaceConstSP::attachToRef(VolSurf);
    myFXAsset = FXAssetConstSP::attachToRef(fxAsset);
    myEqFXCorr = CorrelationConstSP::attachToRef(eqFXCorr);
    fwdStart = CVolParam::FwdStart(VolSurf->getBaseDate(), volRequest, VolSurf->getTimeMetric(), asset);

    myVolRequest = CVolRequestDVFConstSP(copy(volRequest));

    if (volRequest->isFwdStarting())
    {
        startDate = volRequest->getStartDate();
    }
    else
        startDate = VolSurf->getBaseDate();  //if fwd start must be volRequest->getStartDate() //VolSurf->getBaseDate()
}

CVolProcessedDVFParam::CVolProcessedDVFParam(
    const CClassConstSP& clazz,
    const CVolBase*         vol,
    const CVolParamConstSP& volParam,
    const CVolRequestDVF* volRequest,
    const CAsset* asset,
    const VolSurface* VolSurf,
    const FXAsset*   fxAsset,
    const  Correlation* eqFXCorr):
    CVolProcessedDVF(clazz),
    forwardmid(0),
    growthRate(0),
    maturity(0)
{
    Init(vol,
         volParam,
         volRequest,
         asset,
         VolSurf,
         fxAsset,
         eqFXCorr);
}

CVolProcessedDVFParam::CVolProcessedDVFParam(
    const CVolBase*         vol,
    const CVolParamConstSP& volParam,
    const CVolRequestDVF*    volRequest,
    const CAsset*           asset,
    const VolSurface* VolSurf) :
    CVolProcessedDVF(TYPE),
    forwardmid(0),
    growthRate(0),
    maturity(0)
{
    Init(vol,
         volParam,
         volRequest,
         asset,
         VolSurf);
}

CVolProcessedDVFParam::CVolProcessedDVFParam(
    const CVolBase*         vol,
    const CVolParamConstSP& volParam,
    const CVolRequestDVF*   volRequest,
    const CAsset*           asset,
    const VolSurface*      VolSurf,
    const FXAsset*          fxAsset,
    const Correlation*      eqFXCorr) :
    CVolProcessedDVF(TYPE),
    forwardmid(0),
    growthRate(0),
    maturity(0)
{
    Init(vol,
         volParam,
         volRequest,
         asset,
         VolSurf,
         fxAsset,
         eqFXCorr);
}

/** identifies the market data name of the volatility */
string CVolProcessedDVFParam::getName() const
{
    return myVolSurf->getName();
}

/** calculates the trading time between two dates */
double CVolProcessedDVFParam::calcTradingTime(const DateTime &date1,
                                              const DateTime &date2) const
{
    return GetTimeMetric()->yearFrac(date1, date2);
}


double CVolProcessedDVFParam::computeImpVol(const DateTime& maturity,
                                            double          strike) const {
    if (fwdStart.isFwdStarting()){
        CSliceDouble singlePoint(1);
        DateTimeArray maturities(1, maturity);
        CSliceDouble strikes(1);
        strikes[0] = strike;
        myParamVol->computeFwdStartImpVol(myVol.get(), fwdStart,
                                          strikes, false, // abs strike
                                          maturities, singlePoint);
        return singlePoint[0];
    }
    return myParamVol->ComputeImpVol(myVol.get(), strike, maturity);
}

/** Simple wrapper around myParamVol->ComputeImpVol which selects right
    method to invoke for forward starting */
void CVolProcessedDVFParam::computeImpVol(const CLatticeDouble&    strikes,
                                          const DateTimeArray&     maturities,
                                          CLatticeDouble&          impV) const{
    if (fwdStart.isFwdStarting()){
        myParamVol->computeFwdStartImpVol(myVol.get(), fwdStart,
                                          strikes, false, // abs strike
                                          maturities, impV);
    } else {
        myParamVol->ComputeImpVol(myVol.get(), strikes, maturities, impV);
    }
}

/** retieve time measure for the vol */
TimeMetricConstSP CVolProcessedDVFParam::GetTimeMetric() const
{
    return myVolSurf->getTimeMetric();
}

inline DateTime CVolProcessedDVFParam::tweakDate(
    const DateTime& date, double& timeTweak) const
{
    double yearFrac;
    if (myVolRequest->getUseUnscaledTimeTweak())
        yearFrac = 1.0;
    else
        yearFrac= calcTradingTime(startDate, date);
    timeTweak = Maths::max(MIN_TIME_TWEAK, myVolRequest->getTimeTweakUnscaled() * yearFrac);
    return myVolSurf->getTimeMetric()->impliedTime(date, timeTweak, timeTweak);
}

/** volatility calculator object using a cache of forward values */
class IVolCalculatorDVF : public CVolProcessedDVF::IVolCalculator {
public:
    /** initialise the caching of forward value storing the maturities,
        the forwardmid and growthrate */
    IVolCalculatorDVF(
        const CVolProcessedDVFParam* Ivol,
        const DateTimeArray&   Maturities, // the time axis of the lattice above
        bool                   isIntraDayInterp):
        forwardmid(Maturities.size()-1,0.0),
        growthRate(Maturities.size()-1,0.0),
        maturity(Maturities.size())
        {
            static const string method("CVolProcessedDVFParam::VolCalculatorDVF");
            try {
                vol = CVolProcessedDVFParamSP::attachToRef(Ivol);

                int iDate;

                /* create a copy of the maturity dates */
                for(iDate=0; iDate<maturity.size(); iDate++) {
                    maturity[iDate] = Maturities[iDate];
                }

                /* forward value fpr maturity dates */
                DoubleArray forward(maturity.size());
                vol->asset->fwdValue(maturity,forward);

                DateTimeArray maturityNext(maturity.size()-1);
                double notUse;

                /* creates maturity dates + epsilon */
                for(iDate=0; iDate<maturityNext.size(); iDate++) {
                    if (vol->GetTimeMetric()->impliedTime(
                            maturity[iDate],
                            CVolProcessedDVFParam::MIN_NEXT_DATE,
                            notUse) > maturity[iDate+1])
                    {
                        maturityNext[iDate] = vol->GetTimeMetric()->impliedTime(
                            maturity[iDate],
                            CVolProcessedDVFParam::MIN_NEXT_DATE,
                            notUse);
                    }
                    else {
                        maturityNext[iDate] =  maturity[iDate+1];
                    }
                }

                /* forward values for the maturity dates + epsilon */
                DoubleArray forwardNext(maturityNext.size());
                vol->asset->fwdValue(maturityNext,forwardNext);

                /* fraction of the year between the maturity dates and maturity dates + epsilon */
                DoubleArray yearFrac(maturity.size()-1);
                for(iDate=0; iDate<yearFrac.size(); iDate++) {
                    yearFrac[iDate] = vol->calcTradingTime(maturity[iDate], maturityNext[iDate]);
                }

                /* maturity dates middle */
                DateTimeArray maturitymid(maturity.size()-1);
                double notUsed;
                if (vol->myVolRequest->getUseMidPoint()) {
                    for(iDate=0; iDate<maturitymid.size(); iDate++) {
                        maturitymid[iDate] = vol->myVolSurf->getTimeMetric()->impliedTime(maturity[iDate], yearFrac[iDate] *.5, notUsed);
                    }
                }
                else {
                    for(iDate=0; iDate<maturitymid.size(); iDate++) {
                        maturitymid[iDate] = maturity[iDate];
                    }
                }

                /* forward values for maturity dates middle */
                if (vol->myVolRequest->getUseMidPoint()) {
                    vol->asset->fwdValue(maturitymid,forwardmid);
                }
                else {
                    for(iDate=0; iDate<forwardmid.size(); iDate++) {
                        forwardmid[iDate] = forward[iDate];
                    }
                }

                if( isIntraDayInterp) {
                    DateTimeArray maturityTp1(maturity.size());
                    for(iDate=0; iDate<maturityTp1.size(); iDate++) {
                        maturityTp1[iDate] = maturity[iDate].rollDate(1);
                    }

                    DateTimeArray maturityNextTp1(maturityNext.size());
                    for(iDate=0; iDate<maturityNextTp1.size(); iDate++) {
                        maturityNextTp1[iDate] =  maturityNext[iDate].rollDate(1);
                    }

                    DoubleArray forwardTp1(maturityTp1.size());
                    vol->asset->fwdValue(maturityTp1,forwardTp1);

                    DoubleArray forwardNextTp1(maturityNextTp1.size());
                    vol->asset->fwdValue(maturityNextTp1,forwardNextTp1);

                    for(iDate=0; iDate<forward.size(); iDate++) {
                        forward[iDate] += (forwardTp1[iDate] - forward[iDate]) * maturity[iDate].getTime() /  DateTime::TIME_IN_DAY;
                    }

                    for(iDate=0; iDate<forwardNext.size(); iDate++) {
                        forwardNext[iDate] += (forwardNextTp1[iDate] - forwardNext[iDate]) * maturityNext[iDate].getTime()/ DateTime::TIME_IN_DAY;
                    }
                }

                /* cache growth rates */
                for(iDate=0; iDate<growthRate.size(); iDate++) {
                    if(!Maths::isZero(yearFrac[iDate])) {
                        growthRate[iDate] = (forwardNext[iDate] / forward[iDate] - 1.0) / yearFrac[iDate];
                    }
                }
            }

            catch (exception& e) {
                throw ModelException(e, method);
            }
        }

    /* use the caching of forward values */
    void CalcLocVol(
            const CSliceDouble& strikes,        // Slice of strikes
            int                 iDate,          // index in the array of maturities
            CSliceDouble&       locVol)         // output
        {
            static const string method("CVolProcessedDVFParam::CalcLocVol");
            try {
                if (forwardmid.size()==0) {
                    throw ModelException(method, "the cache for forward values has not been initialized");
                }
                else {

                    DateTime maturityNext;
                    double notUse;

                    if (vol->GetTimeMetric()->impliedTime(maturity[iDate],1/3000,notUse) > maturity[iDate+1]) {
                        maturityNext =  vol->GetTimeMetric()->impliedTime(maturity[iDate],1/3000,notUse);
                    }
                    else {
                        maturityNext= maturity[iDate+1];
                    }

                    /* create yearFrac */
                    double yearFrac = vol->calcTradingTime(maturity[iDate], maturityNext);

                    /* create maturitymid */
                    DateTime maturitymid;
                    double notUsed;
                    if (vol->myVolRequest->getUseMidPoint()) {
                        maturitymid = vol->myVolSurf->getTimeMetric()->impliedTime(maturity[iDate], yearFrac *.5, notUsed);
                    }
                    else {
                        maturitymid = maturity[iDate];
                    }

                    /* create yrsToMat */
                    double yrsToMat = vol->calcTradingTime(vol->startDate, maturitymid);


                    const double* strike = &strikes[0];
                    int sizeStrikeArray = strikes.size();

                    vector<SImpV>  impV(sizeStrikeArray);

                    if (vol->myVolRequest->getUseNextStepDerivs())
                        vol->computeImpVolAndDerivsNextStep(strike,
                                                            sizeStrikeArray,
                                                            maturity[iDate],
                                                            maturityNext,
                                                            yearFrac,
                                                            impV);
                    else {
                        double timeTweakyrs;
                        DateTime endMat = vol->tweakDate(maturitymid, timeTweakyrs);
                        vol->computeImpVolAndDerivsNextStep(strike,
                                                            sizeStrikeArray,
                                                            maturitymid,
                                                            endMat,
                                                            timeTweakyrs,
                                                            impV);
                    }

                    locVol[0] = 0.0;
                    for (int iStrike = 0; iStrike < sizeStrikeArray; ++iStrike) {

                        if (iStrike > 0)
                            locVol[iStrike] = locVol[iStrike - 1];

                        if (vol->computeLocV2(yrsToMat,
                                              strike[iStrike],
                                              forwardmid[iDate],
                                              growthRate[iDate],
                                              impV[iStrike],
                                              locVol[iStrike])){
                            locVol[iStrike] = sqrt(locVol[iStrike]);
                        }
                    }
                }
            }
            catch (exception& e) {
                throw ModelException(e, method);
            }
        }

private:
    CVolProcessedDVFParamSP vol;
    DoubleArray                  forwardmid;
    DoubleArray                  growthRate;
    DateTimeArray                maturity;
};

/** create a volatility calculator with a caching of forward values  */
CVolProcessedDVF::IVolCalculator* CVolProcessedDVFParam::CreateVolCalculator(
    const DateTimeArray&   maturity,
    bool                   isIntraDayInterp) const
{
    return new IVolCalculatorDVF(this, maturity, isIntraDayInterp);
}

void CVolProcessedDVFParam::computeLocVNextStep(
                            const DateTimeArray&  maturity,
                            const double*         strike,
                            const int             sizeStrikeArray,
                            vector<SImpV>&        impV,
                            vector<double>&       locV,
                            bool                  reqVar,
                            bool                  isIntraDayInterp) const
{
    static const string method("CVolProcessedDVFParam::computeLocVNextStep");
    try {
        if (sizeStrikeArray < 1 ) {
            throw ModelException(method, "Must supply at least one strike");
        }
        double notUsed;
        int iStrike;
        double forwardmid;
        DateTime maturitymid;
        DateTime maturityNext;

        DoubleArray forward(2);
        double yearFrac = calcTradingTime(maturity[0], maturity[1]);

        // get local copy of maturity. floor yearFrac if maturity0/1 are ~EOD and ~SOD on neighbor dates
        DateTimeArray matLoc = maturity;
        if( yearFrac < MIN_TIME_TWEAK && maturity[0].getDate() != maturity[1].getDate())
        {
            yearFrac = MIN_TIME_TWEAK;
            matLoc[1] = myVolSurf->getTimeMetric()->impliedTime(matLoc[0], yearFrac, notUsed);
        }

        if (myVolRequest->getUseMidPoint())
            maturitymid = myVolSurf->getTimeMetric()->impliedTime(matLoc[0], yearFrac *.5, notUsed);
        else
            maturitymid = matLoc[0];

        double yrsToMat = calcTradingTime(startDate, maturitymid);

        asset->fwdValue(matLoc,forward);

        if (myVolRequest->getUseMidPoint())
           forwardmid = asset->fwdValue(maturitymid);
        else
           forwardmid = forward[0];

        if( isIntraDayInterp )
        {
           DoubleArray forwardTp1(2);
           DateTimeArray maturityTp1(2);

           maturityTp1[0] =  matLoc[0].rollDate(1);
           maturityTp1[1] =  matLoc[1].rollDate(1);
           asset->fwdValue(maturityTp1,forwardTp1);
           forward[0] += (forwardTp1[0] - forward[0]) * matLoc[0].getTime() / DateTime::TIME_IN_DAY;
           forward[1] += (forwardTp1[1] - forward[1]) * matLoc[1].getTime() / DateTime::TIME_IN_DAY;
        }

        if (myVolRequest->getUseNextStepDerivs())
           computeImpVolAndDerivsNextStep(strike,sizeStrikeArray,
                                          matLoc[0],matLoc[1],yearFrac,impV);
        else
        {
           double timeTweakyrs;
           DateTime endMat = tweakDate(maturitymid, timeTweakyrs);
           computeImpVolAndDerivsNextStep(strike,sizeStrikeArray,
                                          maturitymid, endMat, timeTweakyrs, impV);
        }

        double growthRate = 0.0;

        if(!Maths::isZero(yearFrac))
        {
            growthRate = (forward[1] / forward[0] - 1.0) / yearFrac;
        }

        locV[0] = 0.0;
        for (iStrike = 0; iStrike < sizeStrikeArray; ++iStrike){
            if (iStrike > 0)
                locV[iStrike] = locV[iStrike - 1];

            if (computeLocV2(yrsToMat,
                             strike[iStrike],
                             forwardmid,
                             growthRate,
                             impV[iStrike],
                             locV[iStrike])){
                if(reqVar)
                    locV[iStrike] *= yearFrac;
                else
                    locV[iStrike] = sqrt(locV[iStrike]);
            }
        }
    }
    catch (exception& e) {
        throw ModelException(e, method);
    }
}

double CVolProcessedDVFParam::CalcLocVol(const DateTime&  maturity,
                                         double    strike,
										 bool		isIntraDayInterp) const
{
    vector<double> returnvol(1, 0.0);
    vector<SImpV>  impV(1);

    computeLocVol(maturity,&strike,1,impV,returnvol, isIntraDayInterp);
    return returnvol[0];
}


void CVolProcessedDVFParam::CalcLocVol(
    CLatticeDouble*  strikes,
    DateTimeArray&   maturities,
    CLatticeDouble*  locVol,
	bool			isIntraDayInterp) const {

    int NbMat = maturities.size();
    if (NbMat < 2 ) {
        throw ModelException("CVolProcessedDVFParam::CalcLocVol",
                             "Must supply at least two dates");
    }

    // wrap around single slice version
    DateTimeArray t(2);
    vector < double> vol;
    int n, i, j;

    for (i =0; i<NbMat-1; i++)
    {
        t[0] = maturities[i];
        t[1] = maturities[i+1];
        n = (*strikes)[i].size();
        vol.resize(n);
        // call single slice version
        CalcLocVol(  t,
                     &((*strikes)[i][0]),
                     n,
                     vol,
					 isIntraDayInterp);
        for(j = 0; j<n; j++)
            (*locVol)[i][j] = vol[j];
    }
};

void CVolProcessedDVFParam::CalcLocVar(
    CLatticeDouble*  strikes,
    DateTimeArray&   maturities,
    CLatticeDouble*  locVar,
	bool			isIntraDayInterp) const {

    int NbMat = maturities.size();
    if (NbMat < 2 ) {
        throw ModelException("CVolProcessedDVFParam::CalcLocVol",
                             "Must supply at least two dates");
    }

    // wrap around single slice version
    DateTimeArray t(2);
    vector < double> var;
    int n, i, j;

    for (i =0; i<NbMat-1; i++)
    {
        t[0] = maturities[i];
        t[1] = maturities[i+1];
        n = (*strikes)[i].size();
        var.resize(n);

        // call single slice version
        CalcLocVar(  t,
                     &((*strikes)[i][0]),
                     n,
                     var,
					 isIntraDayInterp);
        for(j = 0; j<n; j++)
            (*locVar)[i][j] = var[j];
    }
}

void CVolProcessedDVFParam::CalcLocVol(const DateTimeArray&  maturity,
                                       const double*         strike,
                                       const int             sizeStrikeArray,
                                       vector<double>&       locVol,
                                       bool					 isIntraDayInterp) const {

    int NbMat = maturity.size();
    if (NbMat != 2 ) {
        throw ModelException("CVolProcessedDVFParam::CalcLocVol",
                             "Must supply at least two dates");
    }

    vector<SImpV>  impV(sizeStrikeArray);

    computeLocVNextStep(maturity,strike,sizeStrikeArray,impV,locVol, false, isIntraDayInterp);
}

void CVolProcessedDVFParam::CalcLocVar(const DateTimeArray&  maturity,
                                       const double*         strike,
                                       const int             sizeStrikeArray,
                                       vector<double>&       locVar,
									   bool					 isIntraDayInterp) const {

    int NbMat = maturity.size();
    if (NbMat != 2 ) {
        throw ModelException("CVolProcessedDVFParam::CalcLocVol",
                             "Must supply at least two dates");
    }

    vector<SImpV>  impV(sizeStrikeArray);

    computeLocVNextStep(maturity,strike,sizeStrikeArray,impV,locVar, true, isIntraDayInterp);
}

void CVolProcessedDVFParam::computeLocVol(const DateTime&  maturity,
                                          const double*         strike,
                                          const int             sizeStrikeArray,
                                          vector<SImpV>&        impV,
                                          vector<double>&       locV,
										  bool					isIntraDayInterp) const
{
    static const string routine = "CVolProcessedDVFParam::computeLocVol";


    /* yrs to mat, date + tweak, tweak  */
    double timeTweak;
    double yrsToMat = calcTradingTime(startDate, maturity);
    DateTime nextMat = tweakDate(maturity, timeTweak);

    double forward = asset->fwdValue(maturity);
    double nextForward = asset->fwdValue(nextMat);

	if( isIntraDayInterp )
	{
		double forwardTp1 = asset->fwdValue(maturity.rollDate(1));
		double nextForwardTp1 = asset->fwdValue(nextMat.rollDate(1));

		/* Interpolate forwards linearly intraday */
		forward += (forwardTp1 - forward) * maturity.getTime() /
			DateTime::TIME_IN_DAY;
		nextForward += (nextForwardTp1 - nextForward) * nextMat.getTime() /
			DateTime::TIME_IN_DAY;
	}

    /* Growth rate */
    double growthRate = (nextForward / forward - 1.0) / timeTweak;

    /* Compute implieds and derivs using next time step */
    computeImpVolAndDerivsNextStep(strike,sizeStrikeArray,
                            maturity, nextMat, timeTweak, impV);

    /* Loop over maturity dates */
    for (int iStrike = 0; iStrike < sizeStrikeArray; ++iStrike){
        if (computeLocV2(yrsToMat,
                          strike[iStrike],
                          forward,
                          growthRate,
                          impV[iStrike],
                          locV[iStrike])){
            locV[iStrike] = sqrt(locV[iStrike]);
        }
    }
}

//impV must be initialized before calling this function
// result vol^2 is in v2 if computed. returns true if v2 update, false if v2 is not computed
bool CVolProcessedDVFParam::computeLocV2(double  yrsToMat,
                                          double  strike,
                                          double  forward,
                                          double  growthRate,
                                          SImpV&  impV,
                                          double  &v2) const
{
    static const string routine = "computeLocV2";

    if (yrsToMat < SMALL_YEARS_TO_MATURITY){
        /* cum prob = Step function.
           Densities should be Diracs...however there are several examples where this
           approximation yields unstable results.  Moreover, if the computation yielded a
           negative vol, previous node's vol was used, which was 0 for T = 0.

           Instead, we now use implied vol.  */

           v2  = impV.vol;
           v2 *= v2;
            return true;

    }

    /* If vol at given strike and yrsToMat is zero, set loc vol to zero */
    if (Maths::isZero(impV.vol)){
        v2 = 0.0;
        return true;
    }

    /* Start the calculation proper. */
    double rootYTM  = sqrt(yrsToMat);
    double stdDev   = impV.vol * rootYTM;
    if (!Maths::isPositive(stdDev)){
        throw ModelException(routine, "stdDev (" + Format::toString(stdDev) +
                             ") is not positive (sqrt(t) = " +
                             Format::toString(rootYTM) + ", vol = " +
                             Format::toString(impV.vol) + ")");
    }

    double x = log(forward / strike) / stdDev + 0.5 * stdDev;

    /* Probability density. */
    double probDensRatio = impV.vol * impV.d2_dStrike2;
    probDensRatio += x * (x - stdDev) * impV.d_dStrike * impV.d_dStrike;
    probDensRatio *= strike * yrsToMat;
    probDensRatio += 2.0 * x * rootYTM * impV.d_dStrike;
    probDensRatio *= strike;
    probDensRatio += 1.0;

    if (Maths::isPositive(probDensRatio) &&
        (probDensRatio < myVolRequest->getProbDensRatioMin())){
        /* Floor the probDensRatio. */
        probDensRatio = myVolRequest->getProbDensRatioMin();
    }

	if (!Maths::isPositive(probDensRatio)){
		return false;
	}

    /* Local volatility */
    double numerator = impV.d_dT + growthRate * strike *  impV.d_dStrike;
    numerator *= 2.0 * yrsToMat;
    numerator += impV.vol;
    numerator *= impV.vol;

    double v = numerator / probDensRatio;

    if (!Maths::isPositive(v)){
        return false;
    }

    v2 = v;
    return true;
}

void CVolProcessedDVFParam::computeImpVolAndDerivsNextStep(
    const double*         strike,
    const int             sizeStrikeArray,
    const DateTime&       maturity,
    const DateTime&       maturityup,
    const double          timeTweakyrs,
    vector<SImpV>&        impV) const
{
    int iStrike;
    vector<double> impVup(sizeStrikeArray);
    vector<double> impVdown(sizeStrikeArray);
    vector<double> impVTup(sizeStrikeArray);
    vector<double> impVTdown(sizeStrikeArray);
    vector<double> impVmid(sizeStrikeArray);
    vector<double> strikeUp(sizeStrikeArray);
    vector<double> strikeDn(sizeStrikeArray);
    vector<double> strikeTweak(sizeStrikeArray);

    for (iStrike = 0; iStrike < sizeStrikeArray; ++iStrike)
    {
        //  Pre-compute strikeTweakUp
        strikeTweak[iStrike] = myVolRequest->getStrikeTweakUnscaled() * strike[iStrike];
        if (Maths::isZero(strikeTweak[iStrike]))
            strikeTweak[iStrike] = MIN_STRIKE_TWEAK;

        // Implied Vol at up strike
        strikeUp[iStrike] = strike[iStrike] + strikeTweak[iStrike];
        strikeDn[iStrike] = strike[iStrike] - strikeTweak[iStrike];
    }

    DateTimeArray maturities(1, maturity);

    CSliceDouble strike_vector(const_cast<double*>(strike), sizeStrikeArray);
    CSliceDouble vol_vector(&impVmid[0], sizeStrikeArray);
    computeImpVol(strike_vector, maturities,vol_vector);

    CSliceDouble StrikeUp_vector(&strikeUp[0], sizeStrikeArray);
    CSliceDouble vol_up(&impVup[0], sizeStrikeArray);
    computeImpVol(StrikeUp_vector, maturities,vol_up);

    CSliceDouble StrikeDn_vector(&strikeDn[0], sizeStrikeArray);
    CSliceDouble vol_dn(&impVdown[0], sizeStrikeArray);
    computeImpVol(StrikeDn_vector, maturities,vol_dn);

    maturities[0] = maturityup;
    CSliceDouble vol_t1(&impVTup[0], sizeStrikeArray);
    computeImpVol(strike_vector, maturities,vol_t1);

    for (iStrike = 0; iStrike < sizeStrikeArray; ++iStrike)
    {
        impV[iStrike].vol = impVmid[iStrike];

        impV[iStrike].d_dStrike = 0.5 * (impVup[iStrike] - impVdown[iStrike])/ strikeTweak[iStrike];

        impV[iStrike].d2_dStrike2 = (impVup[iStrike] - 2.0 * impV[iStrike].vol + impVdown[iStrike])
                                    / Maths::square(strikeTweak[iStrike]);

		// adding robustness condition
		if (Maths::isZero(timeTweakyrs))
		{
			impV[iStrike].d_dT = 0.0;
		}
		else
		{
			impV[iStrike].d_dT = (impVTup[iStrike] - impV[iStrike].vol) / timeTweakyrs;
		}
    }
}

DRLIB_END_NAMESPACE
