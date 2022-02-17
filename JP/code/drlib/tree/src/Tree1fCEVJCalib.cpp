//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : Tree1fCEVJCalib.cpp
//
//   Description : calibrator for CEV+Jump Tree.
//
//----------------------------------------------------------------------------
//#define DEBUG_CEVJ_TREE

#include "edginc/config.hpp"
#include "edginc/Results.hpp"
#include "edginc/Tree1fCEVJCalib.hpp"
#include "edginc/CEVJ.hpp"

#include "edginc/Nrfns.hpp"

DRLIB_BEGIN_NAMESPACE

CTree1fCEVJCalib::CTree1fCEVJCalib():CTree1fCEVJ(TYPE)
{
}

CTree1fCEVJCalib::~CTree1fCEVJCalib()
{
    //Clear();
}


static CTree1fCEVJCalib *modelref; // for accessing by funcCalib
double funcCalib(double *p)
{
    double result =0.0;

    modelref->CallTree1fPrice(modelref->this_instrument,
                            modelref->this_control,
                            modelref->this_results,
                            p,
                            &result);

    return result;   
};

////////////////////////// CTree1f /////////////////////////
/** main model entry point */
void CTree1fCEVJCalib::Price(CInstrument* instrument, 
                    CControl*    control, 
                    CResults*    results)
{
    static const string method = "CTree1fCEVJCalib::Price";
    try {

        const double ftol = 0.03; // 3% tolerance
        //const double min_vega = 0.0001; // minimum vega used for fitting error weighting
        double fitting_error = 0.0;
        double *p, **xi;
        int i, j, ndim;

        ndim=0;
        if (DOATMVol) ndim++;
        if (DOCEVPower) ndim++;
        if (DOJumpRate) ndim++;
        if (DOJumpMean) ndim++;
        if (DOJumpWidth) ndim++;

        p = new double[ndim+1];
        xi = new double*[ndim+1];
        for (i=0; i<ndim+1; i++)
        {
            p[i] = 0.0;
            xi[i] = new double[ndim+1];
            for (j=0; j<ndim+1; j++)
                xi[i][j]= 0.0;
        }

        modelref = this;
        this_instrument = instrument;
        this_control = control; 
        this_results = results;

        isUpdateVol = true;
        CTree1f::Price(instrument,control,results); //Only to get VolCEVJ.  No need run price.

        isUpdateVol = false;
/*
  DoubleArray& taVol = VolCEVJ->GetParamArr("ATM_VOL");
  DoubleArray& taCEVPower = VolCEVJ->GetParamArr("CEV_POWER");
  DoubleArray& taJR = VolCEVJ->GetParamArr("JUMP_RATE");
  DoubleArray& taJM = VolCEVJ->GetParamArr("JUMP_MEAN");
  DoubleArray& taJW = VolCEVJ->GetParamArr("JUMP_WIDTH");
*/
        i=0;
        if (DOATMVol) p[++i] = VolCEVJ->GetParamArr("ATM_VOL")[0];
        if (DOCEVPower) p[++i] = VolCEVJ->GetParamArr("CEV_POWER")[0];
        if (DOJumpRate) p[++i] = VolCEVJ->GetParamArr("JUMP_RATE")[0];
        if (DOJumpMean) p[++i] = VolCEVJ->GetParamArr("JUMP_MEAN")[0];
        if (DOJumpWidth) p[++i] = VolCEVJ->GetParamArr("JUMP_WIDTH")[0];

        bool isCalib = (i==0 ? false : true);
        
        int numTS = VolCEVJ->GetParamArr("ATM_VOL").size();
        int numK = WeightMatrix.numCols();
        
        CDoubleMatrixSP priceMatrix(new DoubleMatrix(numK, numTS));
        
        for (iMat=0; iMat<numTS; iMat++)
        {
            totalWeight = 0.0;
            for (i=0;i<NumOfPrice;i++)
            {//Only same maturity row.
                if (TargetVegas[i][iMat]>FP_MIN)
                    totalWeight += WeightMatrix[i][iMat];
            }

            for (j=0; j<=ndim;j++)
                xi[j][j] = 0.2*p[j];

            if (totalWeight > 0.0 && isCalib)      // Not Calibrate if no weight
            {
                powell(p, xi, ndim, ftol, &i, &fitting_error, funcCalib);
            
                // Store the optimized results to all values for longer than current Maturity.
                for (i=iMat; i<numTS; i++)
                {
                    j = 0;
                    if (DOATMVol) VolCEVJ->VolCEVJ->ATMVolArr[i] = p[++j];
                    if (DOCEVPower) VolCEVJ->VolCEVJ->CEVPowerArr[i] = p[++j];
                    if (DOJumpRate) VolCEVJ->VolCEVJ->JumpRateArr[i] = p[++j];
                    if (DOJumpMean) VolCEVJ->VolCEVJ->JumpMeanArr[i] = p[++j];
                    if (DOJumpWidth) VolCEVJ->VolCEVJ->JumpWidthArr[i] = p[++j];
                }
            }
          
            CTree1f::Price(instrument,control,results); //To get Final Results.
            
            for (int i=0; i<NumOfPrice; i++)
                (*priceMatrix)[i][iMat] = NodePrice[CurrIdx][i][0];     // Store the price

        }

        results->storeGreek(priceMatrix, Results::DEBUG_PACKET, 
                            OutputNameSP(new OutputName("PriceMatrix")));
            
        // Refresh to store the calibration results
        DoubleArraySP outVol = DoubleArraySP(new DoubleArray(VolCEVJ->GetParamArr("ATM_VOL")));
        DoubleArraySP outCEVPower = DoubleArraySP(new DoubleArray(VolCEVJ->GetParamArr("CEV_POWER")));
        DoubleArraySP outJR = DoubleArraySP(new DoubleArray(VolCEVJ->GetParamArr("JUMP_RATE")));
        DoubleArraySP outJM = DoubleArraySP(new DoubleArray(VolCEVJ->GetParamArr("JUMP_MEAN")));
        DoubleArraySP outJW = DoubleArraySP(new DoubleArray(VolCEVJ->GetParamArr("JUMP_WIDTH")));

        results->storeGreek(outVol, Results::DEBUG_PACKET, 
                            OutputNameSP(new OutputName("ATM_VOL")));
        results->storeGreek(outCEVPower, Results::DEBUG_PACKET, 
                            OutputNameSP(new OutputName("CEV_POWER")));
        results->storeGreek(outJR, Results::DEBUG_PACKET, 
                            OutputNameSP(new OutputName("JUMP_RATE")));
        results->storeGreek(outJM, Results::DEBUG_PACKET, 
                            OutputNameSP(new OutputName("JUMP_MEAN")));
        results->storeGreek(outJW, Results::DEBUG_PACKET, 
                            OutputNameSP(new OutputName("JUMP_WIDTH")));

    }
    catch (exception& e) {
        throw ModelException(e, method);
    }
}

////////////////////////// CTree1f /////////////////////////
/** call Tree1f Price */
void CTree1fCEVJCalib::CallTree1fPrice(CInstrument* instrument, 
                    CControl*    control, 
                    CResults*    results,
                    double*      pInput,
                    double*       sum_diff)
{
    static const string method = "CTree1fCEVJCalib::CallTree1fPrice";
    try{//if it failed, return error value
        int i=0;

        if (DOATMVol) VolCEVJ->VolCEVJ->ATMVolArr[iMat] = pInput[++i];
        if (DOCEVPower) VolCEVJ->VolCEVJ->CEVPowerArr[iMat] = pInput[++i];
        if (DOJumpRate) VolCEVJ->VolCEVJ->JumpRateArr[iMat] = pInput[++i];
        if (DOJumpMean) VolCEVJ->VolCEVJ->JumpMeanArr[iMat] = pInput[++i];
        if (DOJumpWidth) VolCEVJ->VolCEVJ->JumpWidthArr[iMat] = pInput[++i];
        
        /*
          static const double MAX_CEVPOWER = 0.999, MIN_CEVPOWER = 0.001;
          static const double MAX_DIFF_VOL = 3.0, MIN_DIFF_VOL = 0.01;
          static const double MAX_JUMP_WIDTH = 0.5, MIN_JUMP_WIDTH = 0.0;
          static const double MAX_JUMP_MEAN = 0.3, MIN_JUMP_MEAN = -0.3;
          static const double MAX_JUMP_RATE = 20.0, MIN_JUMP_RATE = 0.0;
        */
        if (VolCEVJ->VolCEVJ->ATMVolArr[iMat] < FP_MIN  || 
            VolCEVJ->VolCEVJ->ATMVolArr[iMat] > 2.0     ||
            VolCEVJ->VolCEVJ->CEVPowerArr[iMat] > 0.999 ||
            VolCEVJ->VolCEVJ->CEVPowerArr[iMat] < 0.001 ||
            VolCEVJ->VolCEVJ->JumpRateArr[iMat] > 20.0  ||
            VolCEVJ->VolCEVJ->JumpRateArr[iMat] < 0.0   ||
            VolCEVJ->VolCEVJ->JumpMeanArr[iMat] > 0.3   ||
            VolCEVJ->VolCEVJ->JumpMeanArr[iMat] < -0.3  ||
            VolCEVJ->VolCEVJ->JumpWidthArr[iMat] > 20.0 ||
            VolCEVJ->VolCEVJ->JumpWidthArr[iMat] < 0.0  )
        {
            *sum_diff = 9.0E+20;
        }
        else
        {
            CTree1f::Price(instrument,control,results);

            *sum_diff = 0.0;
            for (i=0; i<NumOfPrice; i++)
            {
                if (TargetVegas[i][iMat]>FP_MIN)
                {//only add non zero Vega points.  zero Vega points are ignored
                    *sum_diff += (NodePrice[CurrIdx][i][0] - TargetPrices[i][iMat])
                        *(NodePrice[CurrIdx][i][0] - TargetPrices[i][iMat])
                        *WeightMatrix[i][iMat]/totalWeight/TargetVegas[i][iMat];
                }
            }

            // Penalty when if Jump rate move much from previous time point
            if (DOJumpRate && iMat > 0) 
            {// if diff is less than 1.0, no penalty.  If JR increase more than 1.0, the diff is multiplied as error.
                *sum_diff *= max( VolCEVJ->VolCEVJ->JumpRateArr[iMat] - VolCEVJ->VolCEVJ->JumpRateArr[iMat-1] , 1.0);
            }
        }
    }
    catch (exception) {
        *sum_diff = 9.0E+20;    // return error value
    }
}

/** get processed vol*/
// need to examine if this is general
void CTree1fCEVJCalib::InitVol()
{
    if (isUpdateVol){
        CTree1fCEVJ::InitVol();
    }
    //else
        // not InitVol();
}

/****** overwrite maturity! */
#ifdef RG
void CTree1fCEVJCalib::Setup(const DateTime&      valDate, 
                    const DateTimeArray& segDates, 
                    const vector<int>&   density, 
                    const DateTimeArray* critDatesIn,
                    double               minGap, 
                    bool                 equalTime, 
                    int                  numOfPrice, 
                    int                  numInsertNode,
                    const DateTimeArray* divCritDates,
                    bool                 isCall,
                    int                  noExerciseWindow)
#else
void CTree1fCEVJCalib::Setup(const DateTime& valDate, const DateTimeArray& segDates, 
                    const vector<int>& density, 
                    const DateTimeArray* critDatesIn,
                    double minGap, bool equalTime, int numOfPrice, 
                    int numInsertNode)
#endif
{
    if (isUpdateVol)
    {
        CTree1f::Setup(valDate, segDates, density, critDatesIn, 
                       minGap, equalTime, numOfPrice, numInsertNode);
    }
    else
    {
        DateTimeArray new_segDates;
        new_segDates.resize(2); // to try using more segments

        const DateTimeArray& tmpDates = VolCEVJ->GetBMDates();
        new_segDates[0] = segDates[0];
        new_segDates[1] = tmpDates[iMat];

        CTree1f::Setup(valDate, new_segDates, density, critDatesIn, 
                       minGap, equalTime, numOfPrice, numInsertNode);
    }

}
class Tree1fCEVJCalibHelper{
public:
    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        clazz->setPublic(); // make visible to EAS/spreadsheet
        REGISTER(CTree1fCEVJCalib, clazz);
        SUPERCLASS(CTree1fCEVJ);

        EMPTY_SHELL_METHOD(defaultTree1fCEVJCalib);
        FIELD(TargetPrices,     "Target Price Matrix");
        FIELD(TargetVegas,      "Target Vegas Matrix");
        FIELD(WeightMatrix,      "Weight Matrix for Calibration.");
        FIELD(DOATMVol,      "Need to Calibrate CEVPoer?");
        FIELD(DOCEVPower,    "Need to Calibrate CEVPoer?");
        FIELD(DOJumpRate,    "Need to Calibrate CEVPoer?");
        FIELD(DOJumpMean,    "Need to Calibrate CEVPoer?");
        FIELD(DOJumpWidth,   "Need to Calibrate CEVPoer?");

    }

    static IObject* defaultTree1fCEVJCalib(){
        return new CTree1fCEVJCalib();
    }
};

CClassConstSP const CTree1fCEVJCalib::TYPE = CClass::registerClassLoadMethod(
    "Tree1fCEVJCalib", typeid(CTree1fCEVJCalib), Tree1fCEVJCalibHelper::load);

DRLIB_END_NAMESPACE
