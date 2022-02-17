#ifndef _VNFMTESTTIME_HPP
#define _VNFMTESTTIME_HPP

#include "edginc/Object.hpp"
//#include "edginc/DateSched.hpp"
#include "edginc/DoubleMatrix.hpp"
#include "edginc/IRVol.hpp"
#include "edginc/YieldCurve.hpp"
#include "esl_types.h"
#include "edginc/Optimizer.hpp"
#include "esl.h"

DRLIB_BEGIN_NAMESPACE

class VnfmT : public CObject{
public:
    static CClassConstSP const TYPE;

protected:

    IRVolWrapper         irvol;
    YieldCurveWrapper    yc;
    DayCountConventionSP daycount;                 // Day count convention
    int                  frequency;                  // Swap frequency
    CIntSP               reference_tenor;            // For spot vol bootsrtapping
    CDoubleMatrix        pC3param;                   // PC3 parameters

    // transient field
    SWAPVOL_DATA         sv_data;                    // Market data container: swap vols 
    CDoubleMatrixSP      yc_data;                    // Market data container: (time, DF, short rate)
    DoubleArraySP        spotvol;                    // Spotvol

public:
    VnfmT(CDoubleMatrix &param, DayCountConventionSP dcc, int freq, CIntSP ref_tenor, 
            IRVolSP ir, YieldCurveSP yc, CClassConstSP const &type=TYPE);

    CDoubleMatrixSP  getSwapVol();                  // calculate model swap vols
    CDoubleArraySP  getSpotVol();                   // do spot vol bootstrapping
    virtual void validatePop2Object();

    static CDoubleMatrixSP	pC3toFix3( CDoubleMatrix &pC3param2);  // calculate fix3 param from pc3 param
    static CDoubleMatrixSP  fix3toPC3( CDoubleMatrix &fix3param2); // calculate pc3 param from fix3 param     
    static CDoubleMatrixSP  pComponents(CDoubleMatrix &pC3param);  // calculate princ. comp.  
    static CDoubleArraySP	hump(CDoubleMatrix &pC3param);		   // calculate princ. comp.  

protected:
    void set_YIELDCURVE_DATA(YieldCurveSP yc);
    double b_factor(int expiry, int tenor, int index);
    double swap_vol(int expiry, int tenor);
    void set_spotvol(int ref_tenor);

    static CDoubleMatrixSP setrKmatrix( double a1,      // (I) mean reversion for first factor 
                                        double a2,      // (I) mean reversion for second factor
                                        double angle);  // (I) angle in rotation matrix

    static double findAngle(double angle, void* data);

    VnfmT(CClassConstSP const &type=TYPE) : CObject(type) {}

private:
    static void load(CClassSP& clazz);
    static IObject* defaultConstructor(void) { return new VnfmT(); }
};

class VnfmTPlus : public VnfmT{
    friend class ObjFuncT;
public:
    static CClassConstSP const TYPE;

protected:
    CDoubleMatrix       weightMatrix;               // Swaption volatilities (zero or one elements)     
    CDoubleMatrixSP     rKmatrix;                   // R*K matrix, see document
    CDoubleMatrix       tenorPairs;                 // Tenor pair for PC3 solver
    CDoubleArray        marketCorrelation;          // market swap rate correlations
	CDoubleSP			humptime;						// Location of hump
    CBoolArray          calibrateCorrelation;       // vector of booleans
    string              calibrator;                 // defines calibrator to use
    int                 nbFunc;                     // calibrator specific variable
    int                 nbVar;                      // calibrator specific variable
	int					spotStart;					// reference point for spot vol calibration
	int					spotEnd;					// reference point for spot vol calibration

public:
    VnfmTPlus(CDoubleMatrix &param, DayCountConventionSP dcc, int freq, CIntSP ref_tenor, IRVolSP ir, 
        YieldCurveSP yc, CDoubleMatrix &weights, CDoubleMatrix &tp, CDoubleArray &mc, CBoolArray &cc, 
        CDoubleSP hp, CClassConstSP const &type=TYPE);
    CDoubleArraySP getCorrelation();                // calculates model swap rate correlations
    CDoubleMatrixSP getPC3Calibration();            // calibrates pc3 params

protected:

    void update_rKmatrix();
	static double getAnglefromHumpTime(double angle, void* data);
    double g_function(int start, int end, int expiry, int tenor1, int tenor2, int index);
    double j_function(int start, int end, int expiry, int n, int m, int index);
    double correlation(int expiry, int tenor1, int tenor2);
    void set_factorWeights();
    static double find_q(double q, void* data); 
    DoubleArray rmse_swap_vol(const CDoubleArray&  x);
    CDoubleMatrixSP getCalibration(string model);

    VnfmTPlus(CClassConstSP const &type=TYPE) : VnfmT(type) {}

private:

    static void load(CClassSP& clazz);
    static IObject* defaultConstructor(void) { return new VnfmTPlus(); }

};

class VnfmTPCMove : public VnfmTPlus{
public:
    static CClassConstSP const TYPE;

protected:
    CDoubleArraySP  pC_move;                // move in basis vectors    
    IntArray        pC_order;               // order of the ortogonal vectors
    CDoubleArraySP  pC_bump;                // bump size to produce sensitivities


public:
    VnfmTPCMove(CDoubleMatrix &param, DayCountConventionSP dcc, int freq, CIntSP ref_tenor, IRVolSP ir, 
        YieldCurveSP yc, CDoubleMatrix &weights, CDoubleMatrix &tp, CDoubleArray &mc, CBoolArray &cc,
        CDoubleSP hp, CDoubleArraySP &pert, IntArray &ortho, CDoubleArraySP &bs);

    CDoubleArraySP  getPCMove();            // given basis move calculates corr. pc3 params
    CDoubleMatrixSP getPCSensitivity();     // claculates principal component sensitivities

private:
    VnfmTPCMove() : VnfmTPlus(TYPE) {}
    static void load(CClassSP& clazz);
    static IObject* defaultConstructor(void) { return new VnfmTPCMove(); }

};


DRLIB_END_NAMESPACE

#endif


