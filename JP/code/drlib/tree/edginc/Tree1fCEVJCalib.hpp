//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : Tree1fCEVJ.hpp
//
//   Description : one factor trinomial tree for CEV;bool Jump process.
//
//----------------------------------------------------------------------------

#ifndef EDG_Tree1fCEVJ_CALIB_HPP
#define EDG_Tree1fCEVJ_CALIB_HPP

#include "edginc/Tree1fCEVJ.hpp"
#include "edginc/CEVJProcessed.hpp"
#include "edginc/FXAsset.hpp"
#include "edginc/VolProcessedBS.hpp"
#include "edginc/Correlation.hpp"
#include "edginc/FXVolBase.hpp"
#include "edginc/DoubleMatrix.hpp"

DRLIB_BEGIN_NAMESPACE

// log-normal tree
class TREE_DLL CTree1fCEVJCalib : public CTree1fCEVJ
{
public:
    static CClassConstSP const TYPE;
    static void load(CClassSP&);
    
    CTree1fCEVJCalib();
    virtual ~CTree1fCEVJCalib();

    /** main calibration entry point */
    virtual void Price(CInstrument*  instrument, CControl* control, CResults* results);
    virtual void CallTree1fPrice(CInstrument*  instrument, CControl* control, CResults* results, 
                                double* pInput, double* sum_diff);

    virtual void InitVol();

    /** tree initialisation */
#ifdef RG
    virtual void Setup(const DateTime&      valDate, 
                       const DateTimeArray& segDates, 
                       const vector<int>&   density, 
                       const DateTimeArray* critDates,
                       double               minGap, 
                       bool                 equalTime, 
                       int                  numOfPrice, 
                       int                  numInsertNode,
                       const DateTimeArray* divCritDates = 0,
                       bool                 isCall = true,
                       int                  noExerciseWindow = 0);
#else
    virtual void Setup(const DateTime& valDate, const DateTimeArray& segDates, 
                       const vector<int>& density, 
                       const DateTimeArray* critDates,
                       double minGap, bool equalTime, int numOfPrice, int numInsertNode);
#endif


    CInstrument* this_instrument; // $unregistered
    CControl*    this_control;  // $unregistered
    CResults*    this_results; // $unregistered

protected:
    DoubleMatrix    TargetPrices;   //Target Price Matrix
    DoubleMatrix    TargetVegas;    //Target Vegas Matrix
    DoubleMatrix    WeightMatrix;   //Weight Matrix for Calibration

    //switch param to calibrate
    bool DOATMVol;
    bool DOCEVPower;
    bool DOJumpRate;
    bool DOJumpMean;
    bool DOJumpWidth;

private:
    friend class Tree1fCEVJCalibHelper;

    double  totalWeight;    //Total weight of weight Matrix.  excluded zero-vega points. $unregistered
    int     iMat;           // index for TS params. $unregistered
    bool    isUpdateVol;    //Call InitVol or not. $unregistered
};

DRLIB_END_NAMESPACE
#endif
