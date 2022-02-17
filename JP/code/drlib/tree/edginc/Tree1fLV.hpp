//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : Tree1fLV.cpp
//
//   Description : one factor trinomial tree for local vol process.
//
//   Author      : Ning Shen
//
//----------------------------------------------------------------------------

#ifndef EDG_TREE1F_LV_HPP
#define EDG_TREE1F_LV_HPP

#include "edginc/Tree1f.hpp"
#include "edginc/MarketDataFetcherLN.hpp"
#include "edginc/VolProcessedDVF.hpp"
#include "edginc/VolProcessedBS.hpp"
#include "edginc/FXAsset.hpp"
#include "edginc/Correlation.hpp"
#include "edginc/FXVolBase.hpp"

DRLIB_BEGIN_NAMESPACE

// log-normal tree
class TREE_DLL CTree1fLV : public CTree1f
{
public:
    friend class Tree1fLVHelper;
    static CClassConstSP const TYPE;

    static const string SPEEDUP_METHOD_NONE;
    static const string SPEEDUP_METHOD_GRID;
    
    CTree1fLV(CClassConstSP type = TYPE);
    virtual ~CTree1fLV(){};

    /** Override default createMDF in order to set the right MDF */
    virtual MarketDataFetcherSP createMDF() const;

    // from new FDModel interface
    virtual void finaliseModel(CControl*    control);

    /**
     * Whether to enable RiskMapping when computing sensitivities for
     * instruments priced using this model
     *
     * Returns riskMappingIrrelevant: this is a plain old log-normal model.
     * See IModel::wantsRiskMapping().
     */

    virtual IModel::WantsRiskMapping wantsRiskMapping() const;

protected:
    /** set up variance array */
    virtual    void PostSetup();

    /** calculate drift and variance in unit of fwd (ie. exp(v_dt)-1 for LN) for current step
        they must be either one per step or one per node (from -BotClip to TopClip)
        returns true if constant probabilities can be used for all nodes (fast roll) */
    virtual bool CalcStepDriftAndVar(const double* s, int start, int end, vector<double>& dvar_dt,    vector<double>* drift);
    
    /** calculate vol for current step for a set of spot levels
        returns number of vol calculated - one (flat for all node) or num */
    virtual int GetStepVol(int step, vector<double>& vol, const double* s, int start, int end);
    
    /** calculate node spacing for each time segment, usually just one segment */
    virtual    void CalcNodeSpacing();

    /** set up tree nodes at the given step */
    virtual void NodeSetup(int step, int idx);
    
    /**   tree node space set up */
    virtual void NodeSetupSpace(double* s, int bot, int top, 
                  double* sLast, int botLast, int topLast, 
                  double f0, double spacing, int step, double driftNext);
    
    /** determines if cached rebranching info should be used for same grid tweaks */
    virtual bool UseRebranchingCache(){
        return true;
    }

    // total variance from step 0
    CDoubleArray Variance; // $unregistered
    
    // some helpers
    /* get processed vol for tree1f */
    virtual void InitVol();
    
    /** calculate a term structure vol^2 */
    virtual void CalcV2Term(const DateTime& valDate, const DateTime& startDate,
                            const DateTime& matDate, CTermStructure& v2);

//private: // put this back once CTree1fLVGD::GetStepVol() goes to base
    // log-normal vol
    CVolProcessedDVFSP VolLV; // $unregistered
    string volType;
    
    bool            useTweakingForTimeDerivs;
    double          tweakStrikeUnscaled;
    double          tweakTimeUnscaled;
    double          probDensRatioMin;
    bool            useMidPoint;

    string          speedup;

    FXAssetConstSP fxAsset; // $unregistered
 
    const Correlation* eqFXCorr; // $unregistered

    CVolProcessedBSSP volFXBS; // $unregistered

    refCountPtr<CVolProcessedDVF::IVolCalculator> volCalculator; // $unregistered
};

DRLIB_END_NAMESPACE
#endif
