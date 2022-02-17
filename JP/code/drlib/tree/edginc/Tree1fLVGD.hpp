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
#ifndef TREE1F_LVGD
#define TREE1F_LVGD

#include "edginc/Tree1fLV.hpp"

DRLIB_BEGIN_NAMESPACE

class TREE_DLL CTree1fLVGD : public CTree1fLV //, public IAggregateModel
{
public:
    friend class Tree1fLVGDHelper;
    static CClassConstSP const TYPE;
    static void load(CClassSP&);
    static IObject* defaultTree1fLVGD();

    /***override the entry point 
        this is to allow possible use of an un-damped price array to estimate local gamma values
        main model entry point *********/
    virtual void Price(CInstrument* instrument, 
                        CControl*    control, 
                        CResults*    results);

    /* returns number of vol calculated - one (flat for all node) or num 
    overriden here to scale local vol first */
    //virtual bool CalcStepDriftAndVar(const double* s, int start, int end, vector<double>& dvar_dt,    vector<double>* drift);
    
    /* returns number of vol calculated - one (flat for all node) or num 
    overriden here to scale local vol first */
    virtual int GetStepVol(int step, vector<double>& vol, const double* s_inp, int start, int end);
    
    /** input data validation */
    virtual void validatePop2Object();

    virtual ~CTree1fLVGD(){};
    CTree1fLVGD(CClassConstSP clazz);

private:
    CTree1fLVGD();

    double      gammaStart;
    double      gammaWidth;
    double      scalingRange;

    // trnasient
    // local scaling factor storage
    DoubleArrayArray    scalingFactors;
    bool                isRepricing;
};

DRLIB_END_NAMESPACE

#endif
