//----------------------------------------------------------------------------
//
//   Group       : Convertibles DR
//
//   Filename    : FDUtils.cpp
//
//   Description : Some utility functions useful for FD. Mainly interpolation.
//
//   Author      : Tycho von Rosenvinge
//
//   Date        : June 24, 2002
//
//
//----------------------------------------------------------------------------

#ifndef FD_UTILS_HPP
#define FD_UTILS_HPP

#include "edginc/config.hpp" // needed for precompiled headers
#include "edginc/Schedule.hpp"
#include "edginc/TimeMetric.hpp"
#include "edginc/TimeLine.hpp"
#include "edginc/SparseMatrix.hpp"

DRLIB_BEGIN_NAMESPACE

#define CHEM_EPSILON9 1.e-9
#define SUCCESS 0
#define FAILURE 1
#define DBL_THRESHHOLD (double)1.0e-11    /* good to 10 decimal places */
#define DBL_EQUAL(a,b) (fabs((a) - (b)) < DBL_THRESHHOLD)
#ifndef TRUE
#define TRUE 1
#endif
#ifndef FALSE
#define FALSE 0
#endif
#ifndef MAX
#define MAX(a,b) ((a) > (b) ? (a) : (b))
#endif

#ifndef MIN
#define MIN(a,b) ((a) < (b) ? (a) : (b))
#endif

int FDInterpolation(
    int        nOld,        /* (I) size of the table to interpolate from*/
    double*    xOld,        /* (I) x values of the table to interpolate from */
    double*    yOld,        /* (I) y values of the table to interpolate from */
    int        n,            /* (I) size of the output array */
    double*    x,        /* (I) array of x values at which to interpolate */
    double*    y);        /* (O) array of interpolated values for xvec*/

int FDInterpolationD(
    int        nOld,    /* (I) size of the table to interpolate from*/
    double*    xOld,    /* (I) x values of the table to interpolate from */
    double*    yOld,    /* (I) y values of the table to interpolate from */
    int        n,        /* (I) size of the output array */
    double*    x,        /* (I) array of x values at which to interpolate */
    double*    y,        /* (O) array of interpolated values for xvec*/
    double*    y1,        /* (O) array of first derivatives */
    double*    y2);    /* (O) array of second derivatives */

int FDCubicSplineInterpD(
    int tabsize,            /* (I) size of the table to interpolate from */
    double *xtab,            /* (I) x values of the table to interpolate from */
    double *ytab,            /* (I) y values of the table to interpolate from */
    int vecsize,            /* (I) size of the output array */
    double *xvec,            /* (I) array of x values at which to interpolate */
    double *yvec,            /* (O) array of interpolated values */
    double *deriv1,            /* (O) first derivative */
    double *deriv2);        /* (0) second derivative */

int FDCubicSplineInterp(
    int        tabsize,    /* (I) size of the table to interpolate from*/
    double*    xtab,        /* (I) x values of the table to interpolate from */
    double*    ytab,        /* (I) y values of the table to interpolate from */
    int        vecsize,    /* (I) size of the output array */
    double*    xvec,        /* (I) array of x values at which to interpolate */
    double*    yvec);      /* (O) array of interpolated values for xvec*/

int FDInterpolation1F(
    int        nOld,        // (I) size of the table to interpolate from
    double*    xOld,        // (I) x values of the table to interpolate from
    double*    yOld,        // (I) y values of the table to interpolate from
    int        n,            // (I) size of the output array
    double*    x,            // (I) array of x values at which to interpolate
    double*    y);            // (O) array of interpolated values for xvec


int FDInterpolationD1F(
    int        nOld,        // (I) size of the table to interpolate from
    double*    xOld,        // (I) x values of the table to interpolate from
    double*    yOld,        // (I) y values of the table to interpolate from
    int        n,            // (I) size of the output array
    double*    x,            // (I) array of x values at which to interpolate
    double*    y,            // (O) array of interpolated values for xvec
    double*    y1,            // (O) array of first derivatives
    double*    y2);        // (O) array of second derivatives

int    GridInterpolate1F(
    int        mOld,
    double*    Sold,
    double*    Vold,
    int        m,
    double*    S,
    double*    V,                // (O) vector of option values
    double    upBarrier,
    double    upPayout,
    double    upPayoutDelta,
    double    downBarrier,
    double    downPayout,
    double    downPayoutDelta);


//fts for FD1F and FD2F solver
/** FDUtils algorithm class that supports FDUtils */
class TREE_DLL FDUtils{

public:
    
    static void euroOneStep(
        int bot, int top,
        double * alpha, 
        double * beta, 
        double * gamma,
        double * A, 
        double * B, 
        double * C,
        double* sol_n,     
        double* sol_nplus1,
        bool solveByLine ,
        int which_i,
        int dim2);


    static void euroOneStepWithSource(
        int bot, int top,
        double * alpha,
        double * beta,
        double * gamma,
        double * A,
        double * B,
        double * C,
        double* sol_n,
        double* sol_nplus1,
        double* source,
        bool solveByLine,
        int which_i,
        int dim2);

	static void euroOneStepWithSource(
		int bot, 
		int top, 
		SparseCollection aCurr, 
		SparseCollection aPrev,                       
		double * solCurr, 
		double * solPrev, 
		double * source,
		double interval,
		double theta);

private:

    static void computeRhs(
        int bot, int top,
        double * alpha, 
        double * beta, 
        double * gamma,
        double * A, 
        double * B, 
        double * C,
        double * rhs,                                 
        double* sol_n,                    
        double* sol_nplus1,
        bool solveByLine,
        int which_i,
        int dim2);

    static void computeRhsWithSource(
        int bot, int top,
        double * alpha, 
        double * beta, 
        double * gamma,
        double * A, 
        double * B, 
        double * C,
        double * rhs,     
        double* sol_n,
        double* sol_nplus1, 
        double* source,
        bool solveByLine ,
        int which_i,
        int dim2);

    static void triDiag2D(
        int bot, int top,    
        double * alpha,
        double * beta, 
        double * gamma,
        double * rhs, 
        double* sol,
        bool solveByLine ,
        int which_i,
        int dim2);

	static void sparseLinearSystem(
		int bot, 
		int top, 
		SparseCollection & coeff, 
		double * source, 
		double * solCurr);

//-------------------------------------------------------------------------------
//set more pts at time step around barriers
private:

    //isAddedSeg
    //-1 default value, 0 added +-1d segment, 1 just normal segment
    static void calcSegBarDates(const DateTime& valueDate,
                               const DateTime& matDate,
                               TimeMetricConstSP metric,
                               const DateTimeArray& barDates,
                               DateTimeArray& segDates, 
                               IntArray* isAddedSeg);

    static void SetSegBarDates1Bar(const DateTime& valueDate,
                                       const DateTime& matDate,
                                       TimeMetricConstSP metric,
                                       const ScheduleSP bar, 
                                       DateTimeArray& segDates, IntArray* isAddedSeg);

    static void SetSegBarDates2Bar(const DateTime& valueDate,
                                       const DateTime& matDate,
                                       TimeMetricConstSP metric,
                                       const ScheduleSP ubar, const ScheduleSP lbar, 
                                       DateTimeArray& segDates, IntArray* isAddedSeg);


public:
    static  void SetSegDates(const DateTime& valueDate,
                                        const DateTime& matDate,
                                        TimeMetricConstSP metric,
                                        const ScheduleSP ubar, const ScheduleSP lbar, 
                                        const string& uBarType, const string& lBarType,
                                        DateTimeArray& segDates, IntArray* isAddedSeg);

    static  void SetSegDatesGen(const DateTime& valueDate,
                                    const DateTime& matDate,
                                    TimeMetricConstSP metric,
                                    const vector<ScheduleSP> bar, 
                                    const vector<string>& barType, 
                                    DateTimeArray& segDates, IntArray* isAddedSeg);

//-------------------------------------------------------------------------------

private:

    static double setCritSpacePts( const DateTime& valueDate,
                                const  DateTime& matDate, const  ScheduleSP in);

public: 

    /** dblBarrier special case: may have no more than 2 barriers */
    static void setCritSpacePtsAll(const  DateTime& valueDate,
                                const  DateTime& matDate,
                            const ScheduleSP ubar, const ScheduleSP lbar, 
                            const  string& uBarType, const string& lBarType,
                            DoubleArray& outCritPts);

    /** general case: may have more than 2 barriers */
    static void setCritSpacePtsAllGen(const  DateTime& valueDate,
                                    const  DateTime& matDate,
                                    const  vector<ScheduleSP> bar,
                                    const  vector<string>& barType,
                                    DoubleArray& outCritPts);

//-------------------------------------------------------------------------------

/**add more space pts at critical pt.
actually, max pt is 2. To be generalize */
public:
    static void addMoreSpacePts(vector<double>& v_dxM, 
                            DoubleArray& critSpacePts, int dim1, int minNe,
                            double w5, double w3, double w2,
                            double& lBound5, double& uBound5, 
                            double& lBound3, double& uBound3, 
                            double& lBound2, double& uBound2, 
                            double oneDVolT);

private:
    static void addMoreSpacePtsAt1Pt(vector<double>& v_dxM, 
                        int dim1, int impNe, double w5, double w3, double w2,
                        double& lBound5, double& uBound5, 
                          double& lBound3, double& uBound3, 
                          double& lBound2, double& uBound2, 
                          double barL, double barR);

    static void addMoreSpacePtsAt2Pts(vector<double>& v_dxM, 
                        int dim1, int impNe, double w5, double w3, double w2,
                        double& lBound5, double& uBound5, 
                          double& lBound3, double& uBound3, 
                          double& lBound2, double& uBound2, 
                          double barLowL, double barLowR,
                          double barUpL, double barUpR);

    static void calcLocaldx(vector<double>& v_dxM, double lowB, double upB, int lDim, int uDim);


//public:
    //given interval [lB, uB], with lB <= [barLowL, barUpR] < [barUpL, barUpR] <= uB,
    //total pts: uN - lN + 2 *impNe
    static void addMoreSpacePtsConst(vector<double>& v_dxM, int impNe, const double lB, const double uB, const int lN, const int uN, 
                        double barLowL, double barLowR, double barUpL, double barUpR);

    

};


//for NewFD
//can be used for other products having Barrier feature
//for KO
class TREE_DLL FDUtilsBarrier: public CObject {


public:

    FDUtilsBarrier();
    virtual ~FDUtilsBarrier();

    //filled by product
    double  upBarrier;  
    double  valueAtUpBar;  
    bool needSpecialFDUpBar;  //only true if product is KO, should be same value for all steps

    double  downBarrier;  
    double  valueAtDownBar;  //value at the critical point, if the prod knows the value in advance, set the value at precalcNewFD such as traditional KO
                            //otherwise, the solver will call payoffBCDs to get an approximated value
    bool needSpecialFDDownBar;  //true if product needs to be solved in a smaller grid, value at grids outside solved FD grid can be given by payoffBCDs
                        //otherwise, false


    //true if the input barrier is in the FD grid
    //in the case when product needs barrier treatment, 
    //these two variable tell us if it's active at each step
    bool hasUpBarAtStep;  
    bool hasDownBarAtStep;

    //cloest index smaller or equal to bar
    int botDimBar;
    int topDimBar;

    double  upBarrierGrid; // bar level in the FD grid, ex: if grid is log(S), then upBarrierGrid = log(upBarrier)
                            //this is done by model
                            //used for interp. 

    double  downBarrierGrid; // bar level in the FD grid, ex: if grid is log(S), then upBarrierGrid = log(upBarrier)
                            //this is done by model

    //s is the state variable where we check barrier

    void calcClosestIndex(int& iDown, int& iUp, 
                    double* s, int sSize);


    void adjustCBDs(const vector<double>& gridLevel1, 
                int iDown, int iUp,  //these two index are the cloest index to the critical points
                double* price, double* lastP) ;
    
};

typedef smartPtr<FDUtilsBarrier> FDUtilsBarrierSP;

typedef array<FDUtilsBarrierSP, FDUtilsBarrier >  FDUtilsBarrierArray;
typedef smartPtr<FDUtilsBarrierArray>             FDUtilsBarrierArraySP;

DRLIB_END_NAMESPACE

#endif
