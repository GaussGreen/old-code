//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : VolExpSmileSkew.cpp
//
//   Description : Parametrised Vol Surface:
//                 Generated from 4 observed vols and 5 smile slopes
//                 See VolExpSmileSkew spec for full formula
//                 
//
//   Author      : Stephen Hope
//
//   Date        : 3rd July 2002
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/VolBaseParamSurface.hpp"
#include "edginc/Maths.hpp"
#include "edginc/VolatilityDVF.hpp"
#include "edginc/Addin.hpp"
#include "edginc/Calibrator.hpp"

DRLIB_BEGIN_NAMESPACE

class VolExpSmileSkew;
typedef smartPtr<VolExpSmileSkew> VESSSP;

class VolExpSmileSkew: public CVolBaseParamSurface,
                       virtual public IVolatilityBS,
                       virtual public IVolatilityDVF,
                       virtual public DeltaSurface::IShift,
                       virtual public Calibrator::IAdjustable{
private:
    /** Our VolParam class. This class exists primarily to keep the
        ComputeImpVol and spotVolSurfaceFromStrikes methods out of the
        CVolBaseParamSurface interface */
    class VQVolParam: public CVolParam{
    public:
        static CClassConstSP const TYPE;
        // constructor
        VQVolParam(): CVolParam(TYPE){}
        
        virtual void ComputeImpVol(const CVolBase*          vol,
                                   const CLatticeDouble&    strikes,
                                   const DateTimeArray&     maturities,
                                   CLatticeDouble&          impV)const{
            // turn the vol into what we must have
            const VolExpSmileSkew* myVol = 
                static_cast<const VolExpSmileSkew*>(vol);
            // then just pass through the parameterised vol
            myVol->computeImpVol(strikes, maturities, impV);
        }
        
        
        /*** Creates a "Matrix Vol Surface" on provided Strikes  ***/
        virtual VolSurface*       spotVolSurfaceFromStrikes(
            const CVolBase*       vol,
            const CDoubleArray&   strikes) const{
            // turn the vol into what we must have
            const VolExpSmileSkew* myVol = 
                static_cast<const VolExpSmileSkew*>(vol);
            // then just pass through the parameterised vol
            return myVol->spotVolSurfaceFromStrikes(strikes);
        }
        
    private:
        static void load(CClassSP& clazz){
            REGISTER(VQVolParam, clazz);
            SUPERCLASS(CVolParam);
            EMPTY_SHELL_METHOD(defaultCtor);
        }
        static IObject* defaultCtor(){
            return new VQVolParam();
        }       
    };
public:
    friend class VQVolParam;
    static CClassConstSP const TYPE;
    
    void validatePop2Object(){
        static const string method("VolExpSmileSkew::validatePop2Object");
        try {
            if (Maths::isNegative(strikeRef) || Maths::equals(strikeRef, 0.0)) {
                throw ModelException(method,
                                     "strikeRef (" + Format::toString(strikeRef) +
                                     ") must be > 0");
            }
            if (Maths::isNegative(obsVol1) || Maths::isNegative(obsVol2) ||
                Maths::isNegative(obsVol3) 
                || Maths::isNegative(obsVol4)) {
                throw ModelException(method,
                                     "observed vols must be positive ");
            }
            if (Maths::isNegative(matPower)) {
                throw ModelException(method,
                                     "matPower (" + Format::toString(matPower) + 
                                     ") must be >= 0");
            }
            if (Maths::isNegative(shortTermConvFactor) || 
                Maths::equals(shortTermConvFactor, 0.0)){
                throw ModelException(method,
                                     "shortTermConvFactor (" +
                                     Format::toString(shortTermConvFactor) + 
                                     ") must be > 0");
            }
        
            if (Maths::isNegative(mediumTermConvFactor) || 
                Maths::equals(mediumTermConvFactor, 0.0)){
                throw ModelException(method,
                                     "mediumTermConvFactor (" + 
                                     Format::toString(mediumTermConvFactor) + 
                                     ") must be > 0");
            }
        
            if (Maths::isNegative(longTermConvFactor) || 
                Maths::equals(longTermConvFactor, 0.0)){
                throw ModelException(method,
                                     "longTermConvFactor (" + 
                                     Format::toString(longTermConvFactor) + 
                                     ") must be > 0");
            }
            // This bit stops division by zero later
            if (Maths::equals(shortTermConvFactor, mediumTermConvFactor) || 
                Maths::equals(shortTermConvFactor, longTermConvFactor) ||
                Maths::equals(mediumTermConvFactor, longTermConvFactor)){
                throw ModelException(
                    method,
                    "shortTermConvFactor, mediumTermConvFactor "
                    "and longTermConvFactor must all be different");
            }
            if (Maths::isNegative(volFloor)) {
                throw ModelException(method,
                                     "volFloor (" +
                                     Format::toString(volFloor) + 
                                     ") must be >= 0");
            }
            if (Maths::isNegative(slope1Limit) || Maths::isNegative(slope2Limit) ||
                Maths::isNegative(slope3Limit) ||
                Maths::isNegative(slope4Limit) || Maths::isNegative(slope5Limit)) {
                throw ModelException(method,
                                     "slope limits must be postive");
            }
        
            if (Maths::isZero(convexityPower))
            {
                throw ModelException(method,
                                     "convexity power cannot be zero !");
            }
        }
        catch (exception& e) {
            throw ModelException(e,method,"Failed for vol ("+getName()+")");
        }   
    }
    
    /** method that builds a CVolParam. */
    virtual CVolParam* createVolParam() const{
       return new VQVolParam();
    }
      
      void addinCalcFourFactor(double& Zsq, double& Ysq, double& Xsq, double& Wsq, const DateTime& date) {
         baseDate = date;
         HolidaySP simpleHoliday(Holiday::weekendsOnly());
         TimeMetricSP simpleTimeMetric(new TimeMetric(0.05,simpleHoliday.get()));
         this->metric = simpleTimeMetric;
         
         calcFourFactor(Zsq, Ysq, Xsq, Wsq);
      }
      
      string sensName(DeltaSurface* shift) const{
         return getName();
      }
      
      bool sensShift(DeltaSurface* shift){
         double shiftSize = shift->getShiftSize();
         // only bother if non zero 
        if (!Maths::isZero(shiftSize)){
            // store our pre-shift self
            baseVESS = VESSSP(copy(this));

            // then update ourselves AFTER the copy
            useSurfaceDelta = true;
            surfaceDelta = DeltaSurfaceSP(copy(shift));           
            strikeRef *= (1.0 + shiftSize);            
        }
        // No point tweaking backbone as not used
        return false; // no more tweaking required here
    }


private:
    
    void computeImpVol(const CLatticeDouble&      strikes,
                       const DateTimeArray&       maturities,
                       CLatticeDouble&            impV) const {
        static const string routine("VolExpSmileSkew::ComputeImpVol");
        
        try
        {            
            if ((maturities.size() != strikes.size()) ||
                (maturities.size() != impV.size())) 
            {
                throw ModelException(routine, 
                                     "Size mismatch between strikes ("+ 
                                     Format::toString(strikes.size()) +
                                     "), maturities ("+ 
                                     Format::toString(maturities.size())+
                                     ") and impV ("+ 
                                     Format::toString(impV.size())+ ")");
            } 
            // calculate the fwd vols from the observed vols - only need to do this once
            calcFourFactor(Zsq, Ysq, Xsq, Wsq);
            // set up arrays of slopes and slope slimits for interpolation later
            DoubleArray slopeLimits(5);
            DoubleArray smileFactors(5);
            DoubleArray slopes(5);
            slopes[0] = slope1; slopes[1] = slope2; slopes[2] = slope3;
            slopes[3] = slope4; slopes[4] = slope5;
            slopeLimits[0] = slope1Limit;
            slopeLimits[1] = slope2Limit;
            slopeLimits[2] = slope3Limit;
            slopeLimits[3] = slope4Limit;
            slopeLimits[4] = slope5Limit;

            MaturityPeriodSP skewMatP(new MaturityPeriod(skewMaturity));
            MaturityPeriodSP smileMatP(new MaturityPeriod(smileMaturity));

            double smileMat = metric->yearFrac(baseDate, smileMatP->toDate(baseDate));
            double skewMat  = metric->yearFrac(baseDate, skewMatP->toDate(baseDate));
            
            for (int iMat = 0; iMat < maturities.size(); iMat++) 
            {
                if (strikes[iMat].size() != impV[iMat].size())
                {
                    throw ModelException(routine, 
                                         "Size mismatch between strikes"
                                         " & maturities for Mat " +
                                         maturities[iMat].toString() +
                                         " (n "+ Format::toString(iMat) + ")");
                }
                double atmVol;
                // floor yearToDate with 1D to avoid blowing up as t --> 0
               double yearToDate = 
                    Maths::max(1.0/DateTime::DAYS_PER_YEAR,
                               metric->yearFrac(baseDate,maturities[iMat]));

                // calculate atmVol from the observed vols
                atmVol = calcATMvol(yearToDate,
                                    smileMat,
                                    skewMat,
                                    slopeLimits,
                                    slopes);

                for (int iStrike = 0; iStrike < strikes[iMat].size(); iStrike ++)
                {                   
                    impV[iMat][iStrike] = computeVol(atmVol,
                                                     strikes[iMat][iStrike],
                                                     yearToDate,
                                                     smileMat,
                                                     skewMat,
                                                     slopeLimits,
                                                     slopes);
                }
            }                            
        }
        catch (exception& e)
        {
            throw ModelException(e, routine);
        }
    }
    
    void calcFourFactor(double& Zsq, double& Ysq, double& Xsq, double& Wsq)const{
        static const string routine("VolExpSmileSkew::calcFourFactor");
        
        try
        {
            double x[] = {obsVol1*obsVol1,
                          obsVol2*obsVol2,
                          obsVol3*obsVol3,
                          obsVol4*obsVol4};
            

            // convert the maturities to year fractions
            MaturityPeriodSP matP1(new MaturityPeriod(maturity1));
            MaturityPeriodSP matP2(new MaturityPeriod(maturity2));
            MaturityPeriodSP matP3(new MaturityPeriod(maturity3));
            MaturityPeriodSP matP4(new MaturityPeriod(maturity4));

            double mat1 = metric->yearFrac(baseDate, matP1->toDate(baseDate));
            double mat2 = metric->yearFrac(baseDate, matP2->toDate(baseDate));
            double mat3 = metric->yearFrac(baseDate, matP3->toDate(baseDate));
            double mat4 = metric->yearFrac(baseDate, matP4->toDate(baseDate));
 
            // get the matrix values (the coefficients of the fwd vols
            double Gz1, Gz2, Gz3, Gz4, Gy1, Gy2, Gy3, Gy4, Gx1, Gx2, Gx3, Gx4, Gw1, Gw2, Gw3, Gw4;
            calcFourFactorCoeff(mat1, Gz1, Gy1, Gx1, Gw1);
            calcFourFactorCoeff(mat2, Gz2, Gy2, Gx2, Gw2);
            calcFourFactorCoeff(mat3, Gz3, Gy3, Gx3, Gw3);
            calcFourFactorCoeff(mat4, Gz4, Gy4, Gx4, Gw4);

            const int matrixSize = 16;
            const int Mdimn = 4;
            double M[] = {Gz1, Gy1, Gx1, Gw1,
                          Gz2, Gy2, Gx2, Gw2,
                          Gz3, Gy3, Gx3, Gw3,
                          Gz4, Gy4, Gx4, Gw4};
            double A[matrixSize]; // the inverse matrix

            // invert the M matrix
            imsl_d_lin_sol_gen(Mdimn, M, NULL, IMSL_INVERSE_USER, A, IMSL_INVERSE_ONLY, 0);


            // imsl seems to insist on a Matrix name of A and a vector name of x 
            double* result = imsl_d_mat_mul_rect("A*x",
                                                 IMSL_A_MATRIX, Mdimn, Mdimn, A,
                                                 IMSL_X_VECTOR, Mdimn, x,
                                                 0);

            if (!result)
            {
                // just to be paranoid
                throw ModelException(routine,
                                     "imsl failed to multiply matrices");
            }

            Zsq = result[0];
            Ysq = result[1];
            Xsq = result[2];
            Wsq = result[3];

            free(result);

        }
        catch (exception& e)
        {
            throw ModelException(e, routine);
        }
    }
    
    // avoid having this code in 2 places
    double computeVol(double             atmVol,
                      double             strike, 
                      double             yearToDate,
                      double             smileMat,
                      double             skewMat,
                      const DoubleArray& slopeLimits,
                      const DoubleArray& slopes) const {
        static const string routine("VolExpSmileSkew::computeVol");
        try {
            // calculate the smile and skew as a function of this strike
            double skewOfK  = calcSkewOfK(strike);
            double smileOfK = interpSmileFactor(strike, slopeLimits, slopes, smileMat);

            // Now calculate the paramVol
            double paramVol = atmVol + exp(-smileDecay*(yearToDate - smileMat))*smileOfK*(pow(smileMat/yearToDate, matPower))+
                exp(-smileDecay*yearToDate)*(1-exp(smileDecay*smileMat))*skewOfK*(pow(skewMat/yearToDate,matPower))+
                (1-exp(-smileDecay*yearToDate))*skewOfK*(pow(skewMat/yearToDate,matPower));

            if (paramVol < volFloor) {
                paramVol = volFloor;
            }

            return paramVol;
        }
        catch (exception& e)   {
            throw ModelException(e, routine);
        }        
    }

    /*** Creates a "Matrix Vol Surface" on provided Strikes  ***/
    VolSurface* spotVolSurfaceFromStrikes(
        const CDoubleArray&   strikes) const{
        static const string routine("VolExpSmileSkew::spotVolSurfaceFromStrikes");
        try
        {
            const VolSurface* backbone = getBackboneSurface();
            const DateTimeArray& dates = backbone->getDates();
            CDoubleMatrix matrix(strikes.size(), numBMs);
            
            // calculate the fwd vols from the observed vols - only need to do this once
            calcFourFactor(Zsq, Ysq, Xsq, Wsq);
            // set up arrays of slopes and slope slimits for interpolation later
            DoubleArray slopeLimits(5);
            DoubleArray smileFactors(5);
            DoubleArray slopes(5);
            slopes[0] = slope1; slopes[1] = slope2; slopes[2] = slope3;
            slopes[3] = slope4; slopes[4] = slope5;
            slopeLimits[0] = slope1Limit;
            slopeLimits[1] = slope2Limit;
            slopeLimits[2] = slope3Limit;
            slopeLimits[3] = slope4Limit;
            slopeLimits[4] = slope5Limit;

            MaturityPeriodSP skewMatP(new MaturityPeriod(skewMaturity));
            MaturityPeriodSP smileMatP(new MaturityPeriod(smileMaturity));

            double smileMat = metric->yearFrac(baseDate, smileMatP->toDate(baseDate));
            double skewMat =  metric->yearFrac(baseDate, skewMatP->toDate(baseDate));

            for (int iMat = 0; iMat < numBMs; iMat++) {
                double atmVol;
                double yearToDate = metric->yearFrac(baseDate, dates[iMat]);

                // calculate atmVol from the observed vols
                atmVol = calcATMvol(yearToDate,
                                    smileMat,
                                    skewMat,
                                    slopeLimits,
                                    slopes);
                
                for (int iStrike = 0; iStrike < strikes.size(); iStrike ++)
                {
                    matrix[iStrike][iMat] = computeVol(atmVol,
                                                       strikes[iStrike],
                                                       yearToDate,
                                                       smileMat,
                                                       skewMat,
                                                       slopeLimits,
                                                       slopes);
                }                
            }
            
            /** for performance need constructor that takes in
                cached values (to do) */
            VolSurface* volSurf = new VolSurface(backbone, strikes, matrix);

            
            return volSurf;
        }
        catch (exception& e)
        {
            throw ModelException(e, routine);
        }
    }

    double interpSmileFactor(double strike, const DoubleArray& slopeLimits, 
                             const DoubleArray& slopes, double smileMat)const{
        static const string routine("VolExpSmileSkew::interpSmileFactor");

        try
        {
            double smileOfK;
            double percStrike = strike/strikeRef;
            int loIndex=0;
            int hiIndex=0;

            // Find the boundary points to calculate the slope
            for (int i=0; i<slopeLimits.size()-1; i++)
            {
                if (Maths::equals(percStrike,slopeLimits[i]))
                {
                    smileOfK = slopes[i]*(log(percStrike)/log(1.1))*(pow((1/smileMat), matPower));
                    //smileOfK = smileFactors[i];
                    break;
                }
                if (percStrike > slopeLimits[i] && percStrike < slopeLimits[i+1])
                {
                    loIndex = i;
                    hiIndex = i+1;
                }
                if (percStrike < slopeLimits[i] && i==0)
                {
                    loIndex = i;
                    hiIndex = i+1;
                }
                if (percStrike > slopeLimits[slopeLimits.size()-1])
                {
                    loIndex = slopeLimits.size()-2;
                    hiIndex = slopeLimits.size()-1;
                }
            } 

            if (loIndex !=0 || hiIndex!=0)
            {
                double gradGamma;
                double gradSmile;
                double diff;
                if (percStrike < slope1Limit)
                {
                    // Differentiate the smile func wrt strike - see doc for details. 
                    // This has the effect of taking the tangent to the curve at the lowest strike limit and
                    // extrapolating the function linearly to zero strike.
                    gradGamma = (slopes[hiIndex] - slopes[loIndex])/(slope2Limit - slope1Limit);
                    gradSmile = (1./log(1.1))*(pow((1./smileMat), matPower))*
                        ((gradGamma * log(slope1Limit)) + (slopes[loIndex]/slope1Limit));
                    double loLimitSmile = slopes[loIndex]*(log(slope1Limit)/log(1.1))*(pow((1./smileMat), matPower));
                    diff = percStrike - slopeLimits[loIndex];
                    smileOfK = (gradSmile*diff)+ loLimitSmile;

                }
                else if (percStrike > slope5Limit)
                {
                    gradGamma = (slopes[hiIndex] - slopes[loIndex])/(slope5Limit - slope4Limit);
                    gradSmile = (1./log(1.1))*(pow((1./smileMat), matPower))*
                        ((gradGamma * log(slope5Limit)) + (slopes[hiIndex]/slope5Limit));
                    double hiLimitSmile = slopes[hiIndex]*(log(slope5Limit)/log(1.1))*(pow((1./smileMat), matPower));
                    diff = percStrike - slopeLimits[hiIndex];
                    smileOfK = (gradSmile*diff)+ hiLimitSmile;
                }
                else
                {
                    // Interpolate the slope level at this strike
                    double slope = (percStrike - slopeLimits[loIndex])*
                        ((slopes[hiIndex] - slopes[loIndex])/(slopeLimits[hiIndex] - slopeLimits[loIndex])) +
                        slopes[loIndex];
                    smileOfK = slope*(log(percStrike)/log(1.1))*(pow((1./smileMat), matPower));
                }
            }
            
            return smileOfK;
        }
        catch (exception& e)
        {
            throw ModelException(e, routine);
        }
    }

    // calculates (k^p - 1.0) / p
    static double calcStrike2Power(double relStrike, double power){
        if (Maths::isZero(relStrike) && 
            Maths::isZero(power)){
            throw ModelException("VolExpSmileSkew::calcStrike2Power",
                                 "0 to the power 0 is undefined.\n"
                                 "Either the implied volatility should not be evaluated at zero strike;\n"
                                 "or convexityPower should be non-zero");
        }
        // we now know relStrike > 0
        // limiting case: power == 0.0
        if (Maths::isZero(power)){
            return log(relStrike);
        }
        // general case: power != 0.0
        return ((pow(relStrike,power) - 1.0) / power);
    }
    
    double calcSkewOfK(double strike)const{
        static const string routine("VolExpSmileSkew::calcSkewOfK");

        try
        {
            double skewOfK;
            MaturityPeriodSP skewMatP(new MaturityPeriod(skewMaturity));
            double skewMat = metric->yearFrac(baseDate, skewMatP->toDate(baseDate));
            skewOfK = skew
                      * calcStrike2Power(strike/strikeRef, convexityPower)
                      * (1/log(1.1))
                      * (pow((1/skewMat),matPower));
            
            return skewOfK;
        }
        catch (exception& e)
        {
            throw ModelException(e, routine);
        }
    }


    double computeVol(double thisStrike, double paSkew, double paCurvature,
                      double paStrikePctCutDown, double paStrikePctCutUp,
                      double paSmoothWidthDown, double paSmoothWidthUp,
                      double paCubic, double taStrikeCubic, int iMat)const{
        
        static const string routine("VolExpSmileSkew::computeVol");
        try
        {
            return 0;
        }
        catch (exception& e)
        {
            throw ModelException(e, routine);
        }
    }    
    
    double calcATMvol(double             yearToDate,
                      double             smileMat,
                      double             skewMat,
                      const DoubleArray& slopeLimits,
                      const DoubleArray& slopes) const{
        static const string routine("VolExpSmileSkew::calcATMvol");
        
        try
        {
            double atmVol;

            if (!useSurfaceDelta) {
                double Gz, Gy, Gx, Gw;
                
                calcFourFactorCoeff(yearToDate, Gz, Gy, Gx, Gw);
                
                atmVol = sqrt(Gz*Zsq + Gy*Ysq + Gx*Xsq + Gw*Wsq);
            }
            else {
                // if we're doing surface delta we want ATM vol at the new spot
                double spot  = surfaceDelta->getSpot();
                double shift = surfaceDelta->getShiftSize();

                atmVol = baseVESS->calcATMvol(yearToDate,
                                              smileMat,
                                              skewMat,
                                              slopeLimits,
                                              slopes);

                atmVol = baseVESS->computeVol(atmVol,
                                              spot * (1.0 + shift),
                                              yearToDate,
                                              smileMat,
                                              skewMat,
                                              slopeLimits,
                                              slopes);
            }
           
            return atmVol;
        }
        catch (exception& e) {
            throw ModelException(e, routine);
        }
    }
    
    void calcFourFactorCoeff(double T, double& Gz, double& Gy, 
                             double& Gx, double& Gw)const{

        static const string routine("VolExpSmileSkew::calcFourFactorCoeff");

        try
        {
            // For ease of reading
            double st = shortTermConvFactor;
            double mt = mediumTermConvFactor;
            double lt = longTermConvFactor;

            Gz = (1.0 - exp(-st*T))/(st*T);

            Gy = ( (st/(st-mt))*((1-exp(-mt*T))/mt)-(st/(st-mt))*((1-exp(-st*T))/st) )*(1/T);

            Gx = ( (st*mt/((mt - lt)*(st-lt)))*( ((1-exp(-lt*T))/lt)-((1-exp(-st*T))/st) )-
                   ((st*mt)/((mt-lt)*(st-mt)))*( ((1-exp(-mt*T))/mt)-((1-exp(-st*T))/st) ) )*(1/T);

            Gw = (  T -((1-exp(-st*T))/st)-(st/(st-mt))*((1-exp(-mt*T))/mt)+(st/(st-mt))*((1-exp(-st*T))/st)-
                    (st*mt/((mt-lt)*(st-lt)))*((1-exp(-lt*T))/lt)+(st*mt/((mt-lt)*(st-lt)))*((1-exp(-st*T))/st)+
                    (st*mt/((mt-lt)*(st-mt)))*((1-exp(-mt*T))/mt)-(st*mt/((mt-lt)*(st-mt)))*((1-exp(-st*T))/st)  )*(1/T);
        }
        catch (exception& e)
        {
            throw ModelException(e, routine);
        }
    }
    
    
    // registered fields
    double strikeRef;
    double obsVol1;
    double obsVol2;
    double obsVol3;
    double obsVol4;
    // following fields will usually use their default values
    double slope1;
    double slope2;
    double slope3;
    double slope4;
    double slope5;
    double skew;
    mutable double convexityPower;  // if this is passed in as 0.0 need to reset it to 0.001
    double smileDecay;
    double matPower;
    double shortTermConvFactor;
    double mediumTermConvFactor;
    double longTermConvFactor;
    double volFloor;
    string maturity1;
    string maturity2;
    string maturity3;
    string maturity4;
    double slope1Limit;
    double slope2Limit;
    double slope3Limit;
    double slope4Limit;
    double slope5Limit;
    string smileMaturity;
    string skewMaturity;
    
    // transient fields (won't appear in dd interface)
    mutable DateTime    baseDate;
    TimeMetricConstSP   metric;
    int                 numBMs;      // number of benchmarks
    mutable double      Zsq;         // fwd vols calculated from observed vols
    mutable double      Ysq;
    mutable double      Xsq;
    mutable double      Wsq;

    bool                useSurfaceDelta;
    DeltaSurfaceSP      surfaceDelta;
    VESSSP              baseVESS;

    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        clazz->setPublic(); // make visible to EAS/spreadsheet
        REGISTER(VolExpSmileSkew, clazz);
        SUPERCLASS(CVolBaseParamSurface);
        IMPLEMENTS(IVolatilityBS);
        IMPLEMENTS(IVolatilityDVF);
        IMPLEMENTS(DeltaSurface::IShift);
        IMPLEMENTS(Calibrator::IAdjustable);
        EMPTY_SHELL_METHOD(defaultCtor);
        FIELD(strikeRef, "strikeRef");
        FIELD(obsVol1, "observed vol 1");
        FIELD(obsVol2, "observed vol 2");
        FIELD(obsVol3, "observed vol 3");
        FIELD(obsVol4, "observed vol 4");
        FIELD(slope1, "slope 1");
        FIELD(slope2, "slope 2");
        FIELD(slope3, "slope 3");
        FIELD(slope4, "slope 4");
        FIELD(slope5, "slope 5");
        FIELD(skew, "skew");
        FIELD(convexityPower, "skewConvexity");
        FIELD(smileDecay, "smile decay");
        FIELD(matPower, "matPower");
        FIELD(shortTermConvFactor, "shortTermConvFactor");
        FIELD(mediumTermConvFactor, "mediumTermConvFactor");
        FIELD(longTermConvFactor, "longTermConvFactor");
        FIELD(volFloor, "floor when calculating vols");
        // following fields will usually use their default values
        FIELD(maturity1, "maturity1");
        FIELD(maturity2, "maturity2");
        FIELD(maturity3, "maturity3");
        FIELD(maturity4, "maturity4");
        FIELD(slope1Limit, "slope 1 limit");
        FIELD(slope2Limit, "slope 2 limit");
        FIELD(slope3Limit, "slope 3 limit");
        FIELD(slope4Limit, "slope 4 limit");
        FIELD(slope5Limit, "slope 5 limit");
        FIELD(smileMaturity, "smile maturity");
        FIELD(skewMaturity, "skew maturity");
        FIELD(baseDate, "Base Date");
        FIELD(metric, "used to throw at the VolSurface constructor");
        FIELD(numBMs, "Number of Benchmarks");
        FIELD(Zsq, "Z fwd vol");
        FIELD(Ysq, "Y fwd vol");
        FIELD(Xsq, "X fwd vol");
        FIELD(Wsq, "W fwd vol");
        // OPTIONAL fields
        FIELD_MAKE_OPTIONAL(slope1);
        FIELD_MAKE_OPTIONAL(slope2);
        FIELD_MAKE_OPTIONAL(slope3);
        FIELD_MAKE_OPTIONAL(slope4);
        FIELD_MAKE_OPTIONAL(slope5);
        FIELD_MAKE_OPTIONAL(skew);
        FIELD_MAKE_OPTIONAL(convexityPower);
        FIELD_MAKE_OPTIONAL(smileDecay);
        FIELD_MAKE_OPTIONAL(matPower);
        FIELD_MAKE_OPTIONAL(shortTermConvFactor);
        FIELD_MAKE_OPTIONAL(mediumTermConvFactor);
        FIELD_MAKE_OPTIONAL(longTermConvFactor);
        FIELD_MAKE_OPTIONAL(volFloor);
        FIELD_MAKE_OPTIONAL(maturity1);
        FIELD_MAKE_OPTIONAL(maturity2);
        FIELD_MAKE_OPTIONAL(maturity3);
        FIELD_MAKE_OPTIONAL(maturity4);
        FIELD_MAKE_OPTIONAL(slope1Limit);
        FIELD_MAKE_OPTIONAL(slope2Limit);
        FIELD_MAKE_OPTIONAL(slope3Limit);
        FIELD_MAKE_OPTIONAL(slope4Limit);
        FIELD_MAKE_OPTIONAL(slope5Limit);
        FIELD_MAKE_OPTIONAL(smileMaturity);
        FIELD_MAKE_OPTIONAL(skewMaturity);
        FIELD_MAKE_TRANSIENT(baseDate);
        FIELD_MAKE_TRANSIENT(metric);
        FIELD_MAKE_TRANSIENT(numBMs);
        FIELD_MAKE_TRANSIENT(Zsq);
        FIELD_MAKE_TRANSIENT(Ysq);
        FIELD_MAKE_TRANSIENT(Xsq);
        FIELD_MAKE_TRANSIENT(Wsq);

        FIELD(useSurfaceDelta, "");
        FIELD_MAKE_TRANSIENT(useSurfaceDelta);
        FIELD(surfaceDelta, "");
        FIELD_MAKE_TRANSIENT(surfaceDelta);
        FIELD(baseVESS, "");
        FIELD_MAKE_TRANSIENT(baseVESS);

        Calibrator::IAdjustable::registerField(
            clazz, "obsVol1",
            new Range(ClosedBoundary(0.0), Infinity(Infinity::Plus)));
        Calibrator::IAdjustable::registerField(
            clazz, "obsVol2",
            new Range(ClosedBoundary(0.0), Infinity(Infinity::Plus)));
        Calibrator::IAdjustable::registerField(
            clazz, "obsVol3",
            new Range(ClosedBoundary(0.0), Infinity(Infinity::Plus)));
        Calibrator::IAdjustable::registerField(
            clazz, "obsVol4",
            new Range(ClosedBoundary(0.0), Infinity(Infinity::Plus)));
        Calibrator::IAdjustable::registerField(
            clazz, "slope1",
            new InfiniteRange());
        Calibrator::IAdjustable::registerField(
            clazz, "slope2",
            new InfiniteRange());
        Calibrator::IAdjustable::registerField(
            clazz, "slope3",
            new InfiniteRange());
        Calibrator::IAdjustable::registerField(
            clazz, "slope4",
            new InfiniteRange());
        Calibrator::IAdjustable::registerField(
            clazz, "slope5",
            new InfiniteRange());
        Calibrator::IAdjustable::registerField(
            clazz, "skew",
            new InfiniteRange());
        Calibrator::IAdjustable::registerField(
            clazz, "convexityPower",
            new InfiniteRange());   // special care will be taken in the evaluation of the 
                                    // function when convexityPower == 0.0
        Calibrator::IAdjustable::registerField(
            clazz, "smileDecay",
            new InfiniteRange());
        Calibrator::IAdjustable::registerField(
            clazz, "matPower",
            new Range(ClosedBoundary(0.0), Infinity(Infinity::Plus)));
#if 0   // dont make these fields ajustable for now
        // cos they are constrained to one another
        Calibrator::IAdjustable::registerField(
            clazz, "shortTermConvFactor",
            new Range(OpenBoundary(0.0), Infinity(Infinity::Plus)));
        Calibrator::IAdjustable::registerField(
            clazz, "mediumTermConvFactor",
            new Range(OpenBoundary(0.0), Infinity(Infinity::Plus)));
        Calibrator::IAdjustable::registerField(
            clazz, "longTermConvFactor",
            new Range(OpenBoundary(0.0), Infinity(Infinity::Plus)));
#endif
    }
    
    VolExpSmileSkew(): CVolBaseParamSurface(TYPE),
                       strikeRef(0.0), obsVol1(0.0), obsVol2(0.0), 
                       obsVol3(0.0), obsVol4(0.0),
                       slope1(-0.12), slope2(-0.1),  slope3(-0.09), 
                       slope4(-0.07), slope5(-0.05),
                       skew(-1.2), convexityPower(0.001), smileDecay(1.4),
                       matPower(0.5), shortTermConvFactor(11.0), 
                       mediumTermConvFactor(3.00), longTermConvFactor(1.00),
                       volFloor(0.06),
                       maturity1("2W"), maturity2("3M"), maturity3("1Y"), 
                       maturity4("5Y"),
                       slope1Limit(0.7), slope2Limit(0.9), slope3Limit(1.0),
                       slope4Limit(1.1), slope5Limit(1.3),
                       smileMaturity("1M"), skewMaturity("5Y"),
                       numBMs(0), Zsq(0.0), Ysq(0.0), Xsq(0.0), Wsq(0.0),
                       useSurfaceDelta(false)
        {
            // empty
        }
    
    static IObject* defaultCtor(){
        return new VolExpSmileSkew();
    }
    
protected:
    /*** Build the parameterised vol and cache any values **/
    void buildCache() {
        const VolSurface* backbone = getBackboneSurface();
        baseDate   = backbone->getBaseDate();
        metric     = backbone->getTimeMetric();
        numBMs     = backbone->getDates().size();
#if 0
        const DateTimeArray&    dates = backbone->getDates();
        // NB This code was the old caching code - hashed out when the
        // term structure to the power was added. So its reintroduction will
        // need careful analysis.
        scaleCutOff = log(90.0 / 100.0) / cutOffSpeedStart;
        scaleCutOff = N1(scaleCutOff) - 0.5;
        scaleCutOff = 1.0 / scaleCutOff;
        
        atmSlopes = CDoubleArray(dates.size());
        
        for (int i = 0; i <atmSlopes.size(); i++) {
            double yearToDate = baseDate.yearFrac(dates[i]); // calendar time
            double powerToDate;
            if (yearToDate <= 1.0) {
                powerToDate  = powerStart * (1.0 - yearToDate);
                powerToDate += power1Y * yearToDate;
            } else {
                powerToDate  = power1Y  / yearToDate;
                powerToDate += powerEnd * (1.0 - 1.0/yearToDate);
            }
            double tToP = pow(yearToDate, powerToDate);
            double thisScaleCutOff = scaleCutOff;
            atmSlopes[i]  = scaleCutOff * skew1Y90100;
            atmSlopes[i] /= tToP;
        }
#endif
    }

    /** Needed for IAdjustable interface. Returns market data name for vol */
    virtual string getName() const{
        return CVolBaseParamSurface::getName();
    }

};

CClassConstSP const VolExpSmileSkew::TYPE =
CClass::registerClassLoadMethod("VolExpSmileSkew", typeid(VolExpSmileSkew), load);

CClassConstSP const VolExpSmileSkew::VQVolParam::TYPE =
CClass::registerClassLoadMethod("VolExpSmileSkew::VQVolParam", 
                                typeid(VQVolParam), load);


/** Addin to calculate four factors */
class VESScalcFourFactorAddin: public CObject{
    static CClassConstSP const TYPE;

    /**  parameters that our addin takes */
    IObjectSP  vol;
    DateTime   baseDate;
  
    /** the 'addin function' */
    static IObjectSP vessCalcFourFactor(VESScalcFourFactorAddin* params){
        static const string routine = "VESScalcFourFactorAddin:vessCalcFourFactor:";
        try
        {
            CDoubleArraySP fourFactors(new CDoubleArray(4));

         
            if (VolExpSmileSkew::TYPE->isInstance(params->vol.get()))
            {
                double Wsq, Xsq, Ysq, Zsq;
                VolExpSmileSkew* vess = dynamic_cast<VolExpSmileSkew*>(params->vol.get());
     
                vess->addinCalcFourFactor(Zsq, Ysq, Xsq, Wsq, params->baseDate);
                (*fourFactors)[0] = Wsq; (*fourFactors)[1] = Xsq; 
                (*fourFactors)[2] = Ysq; (*fourFactors)[3] = Zsq;
            }
            else
            {
                throw ModelException(routine, "vol object must be of type 'VolExpSmileSkew'");
            }

            return fourFactors;
        }
        catch (exception& e) 
        {
            throw ModelException(e, routine);
        }
    }

    /** for reflection */
    VESScalcFourFactorAddin():  CObject(TYPE){}

    
 /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        REGISTER(VESScalcFourFactorAddin, clazz);
        SUPERCLASS(CObject);
        EMPTY_SHELL_METHOD(defaultVESScalcFourFactorAddin);
        FIELD(vol,"VolExpSmileSkew vol");
        FIELD(baseDate, "base date");
        Addin::registerClassObjectMethod(
            "VESS_CALC_FOUR_FACTOR",
            Addin::UTILITIES,
            "returns the VESS four factors",
            TYPE,
            true,
            Addin::returnHandle,
            (Addin::ObjMethod*)vessCalcFourFactor);
    }

    static IObject* defaultVESScalcFourFactorAddin() {
        return new VESScalcFourFactorAddin();
    }
};

CClassConstSP const VESScalcFourFactorAddin::TYPE = CClass::registerClassLoadMethod(
    "VESScalcFourFactorAddin", typeid(VESScalcFourFactorAddin), load);



// external symbol to allow class to be forced to be linked in
bool VolExpSmileSkewLinkIn(){
    return (VolExpSmileSkew::TYPE != 0);
}

DRLIB_END_NAMESPACE

    


