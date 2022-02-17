// -----------------------------------------------------------------------
// This proprietary software has been developed strictly for J.P. Morgan's
// internal use.  Any use or misuse, intentional or otherwise, which
// contradicts or places this policy in jeopardy is strictly forbidden.
//
// Copyright 2000 J.P. Morgan & Co. Incorporated.  All rights reserved.
// -----------------------------------------------------------------------
//
// 1.0 10/23/2003 Afshin Bayrooti
//
// $Header$
//

//#pragma once

#if ! defined(_IR_FIX123_API_)
#define _IR_FIX123_API_

// For now, we assume it is not an issue to use C++ and in particular stl at the API level

#ifdef WIN32
#include <windows.h>  // required for the meantime as windef complains about a typedef of int to INT
#endif
#undef ERROR  // complaints from fix3.h as windows defines an ERROR macro also
#include "IRTraits.h"
#include "fix123.h"
#include "dri_libraryConverters.h"


#include <vector>


namespace IR {

// We use a namespace IR to we don't have to append a prefix to all types and functions.
// This way, we don't polute the global namespace

// Using enums is preferable to using char or int to represent an enum.  For example, reflection information
// can be used to document enum arguments via a hyperlink to the enum definition.
// Validation of such parameters is also automated by enhanced magnet.

// For example, annual and semi-annual curve frequency could be encapsulated in an enum as follows:


//
//   T _ C U R V E _ A P I
//
    
struct T_CURVE_API {

        LG_CLASS_MEMBERS7(
            "Today date",
            Date, TodayDate,
            "Value date",
            Date, ValueDate,
            "Money market basis (360 or 365)",
            int, MMB,
            "Benchmark swap frequency",
            Frequency, SwapFreq,
            "Benchmark swap day count convention",
            YearBasis, SwapDCC,
        
            // Zero curve
            "Zero maturity dates",
            vector<Date>, ZeroDate,
            "Zero rates (Annual ACT/365 basis)",
            vector<double>, Zero )

};

void ToPublicObject( const _T_CURVE* instancePtr, T_CURVE_API* publicInstancePtr );
void FromPublicObject( _T_CURVE* instancePtr, const T_CURVE_API* publicInstancePtr );

LG_DECLARED_CLASS_BEGIN(_T_CURVE)
LG_SET_PUBLIC_CLASS(T_CURVE_API)
LG_DESCRIPTION("This structure stores the zero curve information")
LG_END()


//
//   B A S E V O L _ D A T A
//

// this structure mirrors the relevant parameters in the basevol.dat env input file
struct BaseVolData {

    LG_CLASS_MEMBERS3(
        "Frequency of underlying rate",
        char, Freq,
        "Volatility dates",
        vector<Date>, VolDate,
        "Vol curve",
        vector<double>, Vol)

};

LG_DECLARED_CLASS_BEGIN(BaseVolData)
LG_DESCRIPTION("This structure stores the market volatility data")
LG_END()


//
//   S W A P V O L _ D A T A
//

// this structure mirrors the relevant parameters in the swapvol.dat env input file
struct SwapVolData {

    LG_CLASS_MEMBERS3(
        "Swaption Expiry (in number of months)",
        vector<long>, SwaptionExpiry,
        "Forward starting swap length/Maturity (in years)",
        vector<long>, FwdSwapMaturity,
        "Matrix of swaption volatilities (in percent)",
        Matrix<double>, VolMatrix)

};

LG_DECLARED_CLASS_BEGIN(SwapVolData)
LG_DESCRIPTION("This structure stores the swap volatility data")
LG_END()


// 
//   M K T V O L _ D A T A _ A P I
//

// structure to support the direct specification of relavant market volatilities for the particular 
// instrument

struct MKTVOL_DATA_API {

    LG_CLASS_MEMBERS11(
        // volatility curve
        "Volatility curve base date",
        Date, BaseDate,
        "Number of points in the vol curve",
        int, NbVol,
        "List of volatility dates",
        vector<Date>, VolDate,
        "Vol curve values (in decimal format)",
        vector<double>, Vol,
        "TRUE (=1) if vol used in calibration - note integer type for the meantime but would like to change to vector of bools",
        vector<int>, VolUsed,

        // forward swap information
        "Coupon frequency of underlying swap",
        char, Freq,
        "Day count convention of underlying swap",
        char, DCC,
        "Underlying swap start",
        vector<Date>, SwapSt,
        "Underlying swap maturity",
        vector<Date>, SwapMat,

        // calibration flags
        "Index calibration flag",
        bool, CalibFlag,
        "Skip calibration failure points",
        bool, SkipFlag)

};

void ToPublicObject( const _MKTVOL_DATA* instancePtr, MKTVOL_DATA_API* publicInstancePtr );
void FromPublicObject( _MKTVOL_DATA* instancePtr, const MKTVOL_DATA_API* publicInstancePtr );

LG_DECLARED_CLASS_BEGIN(_MKTVOL_DATA)
LG_SET_PUBLIC_CLASS(MKTVOL_DATA_API)
LG_DESCRIPTION("Structure for the market volatility data relevant to the instrument")
LG_END()


//
//    M O D E L P A R A M E T E R S _ D A T A
//

// this structure mirrors the relevant parameters in the modelParameters.dat env input file
// to stay under 21 parameters, I've implemented some of the inputs as vectors
struct ModelParametersData {

    LG_CLASS_MEMBERS18(

    // one factor tree settings
    "One Factor Tree Mean Reversion",
    double, oneFactorMR,
    "One Factor Tree Volatility",
    double, oneFactorVol,
    "One Factor Tree Periods Per Year",
    double, oneFactorPPY,

    // two factor tree settings
    "Two Factor Tree Mean Reversions",
    vector<double>, twoFactorMR,
    "Two Factor Tree Volatilities",
    vector<double>, twoFactorVol,
    "Two Factor Tree Correlation of Factors One and Two",
    double, twoFactorCorrelation,
    "Two Factor Tree Periods Per Year",
    double, twoFactorPPY,

    // three factor tree settings
    "Three Factor Tree Mean Reversions",
    vector<double>, threeFactorMR,
    "Three Factor Tree Volatilities",
    vector<double>, threeFactorVol,
    "Three Factor Tree Correlation of Factors One and Two",
    double, threeFactorCorrelation12,
    "Three Factor Tree Correlation of Factors One and Three",
    double, threeFactorCorrelation13,
    "Three Factor Tree Correlation of Factors Two and Three",
    double, threeFactorCorrelation23,
    "Three Factor Tree Periods Per Year",
    double, threeFactorPPY,

    // 2 Q parameters
    "2Q Distribution qleft Parameter",
    double, qLeft,
    "2Q Distribution qright Parameter",
    double, qRight,

    "FwdShift",
    double, FwdShift,
    "Number of iterations in Calibration CET",
    int, CetNbIter,
    "Backbone",
    double, backBone)  // make an optional parameter
/*
    LG_X_CLASS_MEMBERS1(
    "Backbone",
    double, backBone, =999.0);  */

};

LG_DECLARED_CLASS_BEGIN(ModelParametersData)
LG_DESCRIPTION("This structure stores the default model parameters for the environment")
LG_END()


//
//    M A R K E T _ D A T A
//

// this is a container for all of the market data structures
// replicating all information in the wrapper env directories
// ?? should this have a public/private representation in order to put the code
// ?? required to convert basevols and swaption vols from percentage to decimal (and
// ?? back again) in the fromPublic function in the same way as the others have
// ?? even though the public and private would be contain exactly the same elements

struct MARKET_DATA {

    LG_X_CLASS_MEMBERS4(
        "Swaption volatility matrix data",
        SwapVolData, swapVol_data, LG_MANDATORY,
        "Base volatility data",
        BaseVolData, baseVol_data, LG_MANDATORY,
        "Zero/discount curve data",
        vector<T_CURVE>, t_curve, LG_MANDATORY,
        "Model Parameters Data",
        //shared_ptr<ModelParametersData>, modelParameters_data, =shared_ptr<ModelParametersData>())
        shared_ptr<ModelParametersData>, modelParameters_data, LG_MANDATORY)
        /*"Flag to use supplied modelParameters Data",
        bool, useModelParameters,)*/

};

LG_DECLARED_CLASS_BEGIN(MARKET_DATA)
LG_DESCRIPTION("This structure stores the market data")
LG_END()


//
//   T R E E _ D A T A _ A P I
//

struct TREE_DATA_API {

    LG_CLASS_MEMBERS4( 
        "Curve to Diffuse",
        int, CvDiff,
        "Discount curve",
        int, CvDisc,

        "Nb of std dev to cut the tree ",
        int, NbSigmaMax,
        "Number of factors",
        int, NbFactor)

};

void ToPublicObject( const _TREE_DATA* instancePtr, TREE_DATA_API* publicInstancePtr );
void FromPublicObject( _TREE_DATA* instancePtr, const TREE_DATA_API* publicInstancePtr );

LG_DECLARED_CLASS_BEGIN(_TREE_DATA)
LG_SET_PUBLIC_CLASS(TREE_DATA_API)
LG_DESCRIPTION("This structure stores the tree data")
LG_END()


//
//   M O D E L _ D A T A
//

// used as overwrites for the default model parameters (set strings to nil)
// or if default model parameters are not available
class MODEL_DATA {
public:
    
    LG_X_CLASS_MEMBERS7(
        "Number of periods per year in tree",
        string, treePPY, LG_MANDATORY,
        "Volatility Calibration Index",
        string, calibrationIndex, LG_MANDATORY,
        "Model type and distribution Parameters",
        string, modelType, LG_MANDATORY, 
        "Factor Weight Overrides",
        string, factorWeights, LG_MANDATORY,
        "Mean reversion coefficients",
        string, meanReversion, LG_MANDATORY,
        "Correlations",
        string, correlations, LG_MANDATORY,
        "Backbone value",
        string, backbone, LG_MANDATORY)

};

LG_CLASS_BEGIN(MODEL_DATA)
LG_DESCRIPTION("This structure defines the product specific model values and overrides")
LG_END()


//
//   O P T _ O U T _ D A T A _ A P I
//

struct OPT_OUT_DATA_API {

    LG_CLASS_MEMBERS2(
        "Output option price",
        double, Option,
        "Additional outputs ([0] = Fwd Swap Price, [1] = Fix leg price, [2] = Forward Leg Price [3-6] = statistics if turned on",
        vector<double>, Price)

};

void ToPublicObject( const _OPT_OUT_DATA* instancePtr, OPT_OUT_DATA_API* publicInstancePtr );
void FromPublicObject( _OPT_OUT_DATA* instancePtr, const OPT_OUT_DATA_API* publicInstancePtr );

LG_DECLARED_CLASS_BEGIN(_OPT_OUT_DATA)
LG_SET_PUBLIC_CLASS(OPT_OUT_DATA_API)
LG_DESCRIPTION("This structure stores the output data")
LG_END()


//
//   P O P U L A T E _ M K T V O L _ D A T A
//

void PopulateMktVolData(
    const MARKET_DATA &market_data,
    const string& calibrationIndex,
    MKTVOL_DATA &mktvol_data);

LG_CLASS_BEGIN(POPULATE_MKTVOL_DATA)
LG_EXECUTABLE(mktvol_data)
LG_DESCRIPTION("Function to populate the mktvol_data structure from the supplied environment and calibration index")
LG_END()

class POPULATE_MKTVOL_DATA {
public:
    void execute() {
        PopulateMktVolData(market_data, calibrationIndex, mktvol_data);
    }

private:

    LG_CLASS_MEMBERS2(
        "Market environment data",
        MARKET_DATA, market_data,
        "Calibration index",
        string, calibrationIndex)

        MKTVOL_DATA mktvol_data;
};

//
//    F I X 3 E r r o r
//

// static/global class to store fix3 error messages
class FIX3Error
{
    private:

    static string _errorMessage;
    static bool _errorDetected;
    FIX3Error();

    public:

    operator bool()
    {
        return _errorDetected;
    }

    static void Initialize()
    {
        _errorMessage = "";
        _errorDetected = false;
    }
    static void SetError(string errMsg)
    {
        _errorMessage = errMsg;
        _errorDetected = true;
    }
    static void AppendError(string errMsg)
    {
        _errorMessage += errMsg;
        _errorDetected = true;
    }
    static string GetError()
    {
        return _errorMessage;
    }
    static bool IsError() 
    {
        return _errorDetected;
    }
};


void CheckArraySize( string name, int size, string vectorName, int vectorSize );


// These two processing functions should be part of the fix3 library
// Had to be pulled out and modified as they mixed code for reading from environment files
// and calculating parameters and populating structures.  Now the environment input data is
// all in the vol and market_data sructures.
// No problem throwing errors from these two functions directly as there is no heap
// memory allocation and they should always be called before any library functions that allocate memory


// Process and populate the market volatility structure from the "raw" market data sent to
// the pricer
// The calibration index specific to the instrument is required for this process also
void CalculateMktVol(
    const string &calibrationIndex,
    const T_CURVE &t_curve,
    const BaseVolData &baseVol_data,
    const SwapVolData &swapVol_data,
    MKTVOL_DATA &mktvol_data);

class MODEL_DATA;

// set the model parameters from either the input modelParameters structure or take
// the defaults supplied with the instrument model_data sructure.  If modelParameters supplied
// apply any overwrites supplied in the instrument's model_data structure (ie. where != "nil")
void Fix3_Param_Input (
    MKTVOL_DATA   *mktvol_data,   
    FIX3_TREE_DATA     *tree_data,     
    const MODEL_DATA& model_data,  // model data overrides supplied by the instrument
    const MARKET_DATA &market_data);

} // IR



#endif
