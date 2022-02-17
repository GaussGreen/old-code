//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : Tree1fLVGD.cpp
//
//   Description : local vol tree with gamma scaling
//
//   Author      : Ning Shen
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/Maths.hpp"
#include "edginc/Results.hpp"
#include "edginc/Tree1fLVGD.hpp"

DRLIB_BEGIN_NAMESPACE

/** class methods **/

CTree1fLVGD::CTree1fLVGD(): CTree1fLV(TYPE),
                        gammaStart(0.0),
                        gammaWidth(1.0),
                        scalingRange(1.0),
                        isRepricing(false){}

CTree1fLVGD::CTree1fLVGD(CClassConstSP clazz): CTree1fLV(clazz),
                        gammaStart(0.0),
                        gammaWidth(1.0),
                        scalingRange(1.0),
                        isRepricing(false){}

void CTree1fLVGD::load(CClassSP& clazz){
        clazz->setPublic(); // make visible to EAS/spreadsheet
        REGISTER(CTree1fLVGD, clazz);
        SUPERCLASS(CTree1fLV);
        EMPTY_SHELL_METHOD(defaultTree1fLVGD);
        FIELD(gammaStart, "default=0, the starting gamma value for scaling to kick in");
        FIELD(gammaWidth, "gamma width for vol decay/increase");
        FIELD(scalingRange, "default=1, the range of vol scaling factor, btween 1.0 to 10.0");

        // transient
        FIELD(scalingFactors, "");
        FIELD_MAKE_TRANSIENT(scalingFactors);
        FIELD(isRepricing, "");
        FIELD_MAKE_TRANSIENT(isRepricing);
}

IObject* CTree1fLVGD::defaultTree1fLVGD(){
    return new CTree1fLVGD();
}

CClassConstSP const CTree1fLVGD::TYPE = CClass::registerClassLoadMethod(
    "Tree1fLVGD", typeid(CTree1fLVGD), load);


/** input data validation */
void CTree1fLVGD::validatePop2Object()
{
    const string method = "CTree1fLVGD::validatePop2Object";

    //need to switch off ctrl_var
    DEBUG_UseCtrlVar = false;

    if (scalingRange< 1.0 || scalingRange >100.0)
        throw ModelException(method, "scalingRange must be >=1.0 and <=100.0");

    if (gammaWidth <=0.0)
        throw ModelException(method, "gammaWidth must be >0.0");
}

/***override the entry point 
    this is to allow possible use of an un-damped price array to estimate local gamma values
    main model entry point *********/
void CTree1fLVGD::Price(CInstrument* instrument, 
                    CControl*    control, 
                    CResults*    results)
{
    // do scaling range first
    if (isTreeRebuilt(control)){
        isRepricing = false;
        CResults dummy; // not used
        // call base method
        bool useGeom = UseSameGridGeometry;
        UseSameGridGeometry = false;
        CTree1fLV::Price(instrument, control, &dummy);
        UseSameGridGeometry = useGeom;
    }
    // pricing
    isRepricing = true;
    CTree1fLV::Price(instrument, control, results);
}

// override to scale vol according to gamma
int CTree1fLVGD::GetStepVol(int step, vector<double>& vol, const double* s_inp, int start, int end)
{
    int j;
    double u, d, gamma, factor;
    // values at next slice
    const double* nextPrice = NodePrice[1- CurrIdx][0];
    const double* nextStock = Stock[1- CurrIdx];

    // call base method first
    CTree1fLV::GetStepVol(step, vol, s_inp, start, end);
    // compute local gamma and scale loc vol, use index[0] array
    // no re-branching considered
    // no ccy struck adjustment !!

    start = Maths::max(start, -BotClip[1- CurrIdx]); // needed as penultimate smoothing go to real end
    end = Maths::min(end, TopClip[1- CurrIdx]);

    if (isRepricing)
    {
        for (j=start; j<=end; j++)
        {
            vol[j-start] *= scalingFactors[step][j-start];
        }
    }
    else // compute scaling factors
    {
        if (scalingFactors.size() == 0) {
            scalingFactors.resize(timeLine->NumOfStep + 1, DoubleArray());
        }
        if (scalingFactors[step].size() == 0) {
            scalingFactors[step].resize(end - start + 1);
            for (j=start; j<=end; j++)
            {
                d = (nextPrice[j] - nextPrice[j-1])/(nextStock[j] - nextStock[j-1]); // lower delta
                u = (nextPrice[j+1] - nextPrice[j])/(nextStock[j+1] - nextStock[j]); // upper delta
                gamma = (u-d)/(nextStock[j+1]-nextStock[j-1])/0.5;
            
                if (gamma > 0.0)
                    gamma = Maths::max(0.0, gamma - gammaStart);
                else
                    gamma = Maths::min(0.0, gamma + gammaStart);

                // compute adjustment factor
                factor = fabs(gamma/gammaWidth);
                factor = 1.0/scalingRange + (1.0 - 1.0/scalingRange)*exp(-factor);
                if (gamma <0.0)
                    factor = 1.0/factor;

                scalingFactors[step][j-start] = factor;
            }
        }
    }
    return (end-start+1);
}

extern bool Tree1fLVGDLoad()
{
    return (true && CTree1fLVGD::TYPE);
}

DRLIB_END_NAMESPACE
