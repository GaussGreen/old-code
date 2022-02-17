//----------------------------------------------------------------------------
//
//   Group       : Convertibles DR
//
//   Filename    : Payoff1F.hpp
//
//   Description : One factor payoff class
//
//   Author      : Tycho von Rosenvinge
//
//   Date        : February 8, 2002
//
//----------------------------------------------------------------------------

#ifndef FDENGINE_HPP
#define FDENGINE_HPP

#include "edginc/config.hpp"
#include "edginc/Payoff1F.hpp"
#include "edginc/TreeSlice.hpp"

#include <vector>

using namespace std;

DRLIB_BEGIN_NAMESPACE


class TREE_DLL FDPayoff1F: public Payoff1F{
public:

    virtual void preCalcFD(int step, int idx, int pStart, int pEnd) {  
        preCalc(step, idx);   
    };
    
    virtual void PayoffAtMatFD(const double * s, int step, int bot, int top, 
        int pStart, int pEnd, double * const * price) {
        PayoffAtMat(s, step, bot, top, pStart, pEnd, price);
    };
    
    virtual void PayoffBeforeMatFD(const double * s, int step, int bot, int top, 
                                   int pStart, int pEnd, double * const * price) {
        PayoffBeforeMat(s, step, bot, top, pStart, pEnd, price);
    };
    
 //   virtual double getVolFD(double stockPrice, int step) = 0;
	virtual void getVolFD(int step, vector<double>& vol, const double * s, int start, int end) = 0;

};


class TREE_DLL FDEngine1F{
public:
    
    /** Simple constructor */
    FDEngine1F();
    virtual ~FDEngine1F();

    
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
    
private:
    int         nTimeSteps;     // No. time steps in entire grid (all segments)
    int         nStockSteps;    // No. stock steps 
    double     *dt;             // Time interval sizes
    double     *forwards;       // Forward prices of underlying asset
    double     *pvs;            // Discount factors to next time step
    double     *ir;             // Forward rates between consecutive time points
    double     *divy;
     
    int     nSegments;          // No. segments
    int    *segEnd;             // Index of the last time point in the segment
    double *segMax;             // maximum stock price needed for the segment 
    double *segMin;             // minimum stock price needed for the segment 
    double *segStrike;          // Used for sinh grid 

    double *stockArray;         // Stock price arrays
    double *stockArrayNew;

    TreeSliceGeneral::RangeSP range;
    TreeSliceGeneralContSP    optionArray;       // Option price arrays
    TreeSliceGeneralContSP    optionOldArray;
    TreeSliceGeneralContSP    optionInterpArray;

    mutable vector<double> sigma;
    double *fwdGrid;
    
    int numPriceArrays;         // No. price arrays / FD layers

    /** clean up */
    virtual void Clear();

    virtual void AllocateArrays();

    int  varMethod;
    bool useFwdGrid;
    int  gridType;              // 0 - sinh, 1 - linear, 2 - exponential
};

DRLIB_END_NAMESPACE

#endif

