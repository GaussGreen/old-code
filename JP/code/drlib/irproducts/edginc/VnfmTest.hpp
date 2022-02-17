#ifndef _VNFMTEST_HPP
#define _VNFMTEST_HPP

#include "edginc/Object.hpp"
//#include "edginc/DateSched.hpp"
#include "edginc/DoubleMatrix.hpp"
#include "edginc/IRVol.hpp"
#include "edginc/YieldCurve.hpp"
#include "esl_types.h"
#include "edginc/Optimizer.hpp"

DRLIB_BEGIN_NAMESPACE

class Vnfm : public CObject{
public:
    static CClassConstSP const TYPE;

protected:

    IRVolWrapper        irvol;
    YieldCurveWrapper   yc;
    DayCountConventionSP  daycount;                 // Day count convention
    int                 frequency;                  // Swap frequency
    CIntSP              reference_tenor;            // For spot vol bootsrtapping
    DoubleArray         pC3param;                   // PC3 parameters

    // transient field
    SWAPVOL_DATA        sv_data;                    // Market data container: swap vols 
    CDoubleMatrixSP     yc_data;                    // Market data container: (time, DF, short rate)
    DoubleArraySP       spotvol;                    // Spotvol

public:
    Vnfm(DoubleArray &param, DayCountConventionSP dcc, int freq, CIntSP ref_tenor, 
            IRVolSP ir, YieldCurveSP yc, CClassConstSP const &type=TYPE);

    CDoubleMatrixSP  getSwapVol();                  // calculate model swap vols
    CDoubleArraySP  getSpotVol();                   // do spot vol bootstrapping
    virtual void validatePop2Object();

    static DoubleArraySP    pC3toFix3(DoubleArray   &pC3param);  // calculate fix3 param from pc3 param
    static DoubleArraySP    fix3toPC3(DoubleArray   &fix3param); // calculate pc3 param from fix3 param     
    static CDoubleMatrixSP  pComponents(DoubleArray &pC3param);  // calculate princ. comp.  


protected:
    void set_YIELDCURVE_DATA(YieldCurveSP yc);
    double b_factor(int expiry, int tenor, int index);
    double swap_vol(int expiry, int tenor);
    void set_spotvol(int ref_tenor);

    static CDoubleMatrixSP setrKmatrix( double a1,      // (I) mean reversion for first factor 
                                        double a2,      // (I) mean reversion for second factor
                                        double angle);  // (I) angle in rotation matrix

    static double findAngle(double angle, void* data);

    Vnfm(CClassConstSP const &type=TYPE) : CObject(type) {}

private:
    static void load(CClassSP& clazz);
    static IObject* defaultConstructor(void) { return new Vnfm(); }
};

class VnfmPlus : public Vnfm{
    friend class ObjFunc;
public:
    static CClassConstSP const TYPE;

protected:
    CDoubleMatrix       weightMatrix;               // Swaption volatilities (zero or one elements)     
    CDoubleMatrixSP     rKmatrix;                   // R*K matrix, see document
    CDoubleMatrix       tenorPairs;                 // Tenor pair for PC3 solver
    CDoubleArray        marketCorrelation;          // market swap rate correlations
    CBoolArray          calibrateCorrelation;       // vector of booleans
    string              calibrator;                 // defines calibrator to use
    int                 nbFunc;                     // calibrator specific variable
    int                 nbVar;                      // calibrator specific variable

public:
    VnfmPlus(DoubleArray &param, DayCountConventionSP dcc, int freq, CIntSP ref_tenor, IRVolSP ir, 
        YieldCurveSP yc, CDoubleMatrix &weights, CDoubleMatrix &tp, CDoubleArray &mc, CBoolArray &cc, 
        CClassConstSP const &type=TYPE);
    CDoubleArraySP getCorrelation();                // calculates model swap rate correlations
    CDoubleArraySP getPC3Calibration();             // calibrates pc3 params
    CDoubleArraySP getFix3Calibration();            // calibrates fix3 params

protected:

    void update_rKmatrix();
    double g_function(int expiry, int tenor1, int tenor2, int index);
    double j_function(int expiry, int n, int m, int index);
    double correlation(int expiry, int tenor1, int tenor2);
    void set_factorWeights();
    static double find_q(double q, void* data); 
    DoubleArray rmse_swap_vol(const CDoubleArray&  x);
    CDoubleArraySP getCalibration(string model);

    VnfmPlus(CClassConstSP const &type=TYPE) : Vnfm(type) {}

private:

    static void load(CClassSP& clazz);
    static IObject* defaultConstructor(void) { return new VnfmPlus(); }

};

class VnfmPCMove : public VnfmPlus{
public:
    static CClassConstSP const TYPE;

protected:
    CDoubleArraySP  pC_move;                // move in basis vectors    
    IntArray        pC_order;               // order of the ortogonal vectors
    CDoubleArraySP  pC_bump;                // bump size to produce sensitivities


public:
    VnfmPCMove(DoubleArray &param, DayCountConventionSP dcc, int freq, CIntSP ref_tenor, IRVolSP ir, 
        YieldCurveSP yc, CDoubleMatrix &weights, CDoubleMatrix &tp, CDoubleArray &mc, CBoolArray &cc,
        CDoubleArraySP &pert, IntArray &ortho, CDoubleArraySP &bs);

    CDoubleArraySP  getPCMove();            // given basis move calculates corr. pc3 params
    CDoubleMatrixSP getPCSensitivity();     // claculates principal component sensitivities

private:
    VnfmPCMove() : VnfmPlus(TYPE) {}
    static void load(CClassSP& clazz);
    static IObject* defaultConstructor(void) { return new VnfmPCMove(); }

};


DRLIB_END_NAMESPACE

#endif


