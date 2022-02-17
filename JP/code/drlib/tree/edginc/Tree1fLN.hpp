//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : Tree1fLN.hpp
//
//   Description : one factor trinomial tree for log-normal process.
//
//   Author      : Ning Shen
//
//----------------------------------------------------------------------------

#ifndef EDG_TREE1FLN_HPP
#define EDG_TREE1FLN_HPP

#include "edginc/Tree1f.hpp"
#include "edginc/MarketDataFetcherLN.hpp"
#include "edginc/VolProcessedBS.hpp"

#include "edginc/FXAsset.hpp"
#include "edginc/Correlation.hpp"
#include "edginc/FXVolBase.hpp"


DRLIB_BEGIN_NAMESPACE

// log-normal tree
class TREE_DLL CTree1fLN : public CTree1f
{
public:
    friend class Tree1fLNHelper;
    static CClassConstSP const TYPE;
    static void load(CClassSP&);
    
    CTree1fLN();
    virtual ~CTree1fLN(){};

    /** Override default createMDF in order to set the right MDF */
    virtual MarketDataFetcherSP createMDF() const;

    /** make a Tree1fLN with most things safely defaulted */
    static CTree1fLN* make(const string& volType, bool useDivTreatment = true);

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

    /**   tree node space set up */
    virtual void NodeSetupSpace(double* s, int bot, int top, 
                  double* sLast, int botLast, int topLast, 
                  double f0, double spacing, int step, double driftNext);

    // total variance from step 0
    CDoubleArray Variance; // $unregistered

//********* some helpers
    /** get processed vol for tree1f */
    virtual void InitVol();
    /** calculate a term structure vol^2 */
    virtual void CalcV2Term(const DateTime& valDate, const DateTime& startDate,
                        const DateTime& matDate, CTermStructure& v2);
    /** overwrite control of sameGrid from instrument, mainly for fwdStarting */
    virtual void controlSameGridFwdStart(const string& ccyTreatment);

private:
    // log-normal vol
    CVolProcessedBSSP VolLN; // $unregistered

    string volType;
};

typedef CTree1fLN Tree1fLN;

typedef smartConstPtr<Tree1fLN> Tree1fLNConstSP;
typedef smartPtr<Tree1fLN> Tree1fLNSP;


DRLIB_END_NAMESPACE
#endif
