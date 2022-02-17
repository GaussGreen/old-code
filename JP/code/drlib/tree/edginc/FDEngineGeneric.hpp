//----------------------------------------------------------------------------
//
//   Group       : Convertibles DR
//
//   Filename    : FDEngineGeneric.hpp
//
//   Description : One factor generic FD Engine
//
//   Author      : André Segger
//
//   Date        : 14 October 2003
//
//----------------------------------------------------------------------------

#ifndef FDENGINE_GENERIC_HPP
#define FDENGINE_GENERIC_HPP

#include "edginc/config.hpp"
#include "edginc/Payoff1F.hpp"
#include "edginc/FDEngine.hpp"
#include "edginc/Model.hpp"
#include "edginc/FDSolver1FGeneric.hpp"
#include "edginc/FDTermStructure.hpp"

#include <vector>

using namespace std;

DRLIB_BEGIN_NAMESPACE


class TREE_DLL FDPayoff1FGeneric: public FDPayoff1F {
public:

	virtual FDTermStructureSP getDriftTerm(int step, const double* s)
    {
        throw ModelException("FDPayoff1FGeneric::getDriftTerm",
            "This method must be overrriden");
    }

	virtual FDTermStructureSP getDiffusionTerm(int step, const double* s)
    {
        throw ModelException("FDPayoff1FGeneric::getDiffusionTerm",
            "This method must be overrriden");
    }

	virtual FDTermStructureSP getCouponTerm(int step, const double* s)
    {
        throw ModelException("FDPayoff1FGeneric::getCouponTerm",
            "This method must be overrriden");
    }

	virtual FDTermStructureSP getDiscountTerm(int step, const double* s)
    {
        throw ModelException("FDPayoff1FGeneric::getDiscountTerm",
            "This method must be overrriden");
    }

    virtual void preCalcFDGeneric(int step, int idx, int pStart, int pEnd,
                                  double const * const * optionArray,
                                  double const * const * optionArrayOld)
    {  
        preCalcFD(step, idx, pStart, pEnd);
    }
    
    virtual bool doLambdaAdjust() 
    { 
        return false; 
    }
};


class TREE_DLL FDEngine1FGeneric {
public:
    
    /** Simple constructor */
    FDEngine1FGeneric();
    virtual ~FDEngine1FGeneric();

    
    double     *upBarrierNew;
	double     *upPayoutNew;
	double     *upPayoutDeltaNew;
	double     *upBarrierOld;
	double     *upPayoutOld;
	double     *upPayoutDeltaOld;
	double     *downBarrierNew;
	double     *downPayoutNew;
	double     *downPayoutDeltaNew;
	double     *downBarrierOld;
	double     *downPayoutOld;
	double     *downPayoutDeltaOld;

    void loop(FDPayoff1F *engineCBs, double stockNow, double *price, double *divPert, double *irPert);

    void init(int timeSteps, int stockSteps, double *times, double *inForwards, double *inPVs, 
              int inNSegments, int *inSegEnd, double *inSegMax, double *inSegMin, double*inSegStrike,
              bool inUseFwdGrid, double *inIr, double *inDivy, int inNumPriceArrays,
              int inGridType);

    int returnNumPriceArrays() { return numPriceArrays;};

    void initModel(IModelSP modelToUse) {
        model = modelToUse.get();
    }

    void FDLambdaAdjustment (FDSolver1FGeneric& fdSolver,
	                         int			    m,
	                         double		        lambda,
	                         double		        downBarrier,
	                         double		        downPayout,
	                         int			    n,
	                         double*		    V);


private:
    int         nTimeSteps;
    int         nStockSteps;
    double     *dt;
    double     *forwards;
    double     *pvs;
    double     *ir;
    double     *divy;
    
    int         nSegments;
    int        *segEnd;
    double     *segMax; /* maximum stock price needed for the segment */
    double     *segMin; /* minimum stock price needed for the segment */
    double     *segStrike; /* used for sinh grid */

    double     *stockArray;
    double     *stockArrayNew;

    double     *assetInterpLevel;

    TreeSliceGeneral::RangeSP range;
    TreeSliceGeneralContSP    optionArray;
    TreeSliceGeneralContSP    optionOldArray;
    TreeSliceGeneralContSP    optionInterpArray;

    mutable vector<double>  sigma;
    double                 *fwdGrid;
    int                     numPriceArrays;

    IModel *model;

    /** clean up */
    virtual void Clear();

    virtual void AllocateArrays();

    int  varMethod;
    bool useFwdGrid;
    int  gridType;
};

DRLIB_END_NAMESPACE

#endif

