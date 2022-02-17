#ifndef QLIB_FD1D_LN_HPP
#define QLIB_FD1D_LN_HPP

#include "edginc/FD1D.hpp"
#include "edginc/MarketDataFetcherLN.hpp"
#include "edginc/MultiFactors.hpp"
#include "edginc/VolProcessedBS.hpp"

DRLIB_BEGIN_NAMESPACE

class TREE_DLL FD1DLN : public FD1D
{
public:
    static const CClassConstSP TYPE;
    static void load( CClassSP & type );

    FD1DLN( const CClassConstSP & type = TYPE );
    ~FD1DLN();

    static IObject * defaultFD1DLN()
    {
        return new FD1DLN();
    }

    /** validate some inputs */
    virtual void validatePop2Object();

    /** Override default createMDF in order to set the right MDF */
    virtual MarketDataFetcherSP createMDF() const;

    /**
     * Whether to enable RiskMapping when computing sensitivities for
     * instruments priced using this model
     *
     * Returns riskMappingIrrelevant: this is a plain old log-normal model.
     * See IModel::wantsRiskMapping().
     */

    virtual IModel::WantsRiskMapping wantsRiskMapping() const;

    CVolProcessedBSConstSP getProcessedVol();
    
    /** get the initial conditions and corresponding indexes in the slices
        for forward induction */
    virtual void getInitialConditions(IntArray& initialIndex, DoubleArray& initialValue) const{};
    
    
protected:
    /** retrieving market data from MDF */    
    virtual void retrieveFactor();
    
    virtual void initModel();

    /**-----------------------------------------------------
       model initialisation  (second part)                 
       most memory init should be here
       finalise model initialisation once after initModel() and product init() 
       called by validate 
       -----------------------------------------------------*/
    
    virtual void finaliseModel(CControl*    control);

    /**-----------------------------------------------------
        each derived model can set diff. boundaries 
        based on the dynamics of underlying,
        alpha is the input truncation1D (nb of std) 
        outLowB, outUpB are low and up boundaries for FD 
    -----------------------------------------------------*/
    virtual void setFdBounds(
        double alpha,
        double & outLowB,
        double & outUpB ) const;

    /**----------------------------------------------------------------
        PDE 1 factor is
        U_t + a*U_x + c*U_xx + f*U + Jump(x) = g;
    ----------------------------------------------------------------*/
    virtual void pdeCoeff(
        int step,
        TreeSliceEQ & a,
        TreeSliceEQ & c,
        TreeSliceEQ & f,
        TreeSliceEQ & g ) const;

    /**----------------------------------------------------------------
        calculate vol for current step for a set of spot levels
        returns number of vol calculated - one (flat for all node) or num 
    ----------------------------------------------------------------*/

    int GetStepVol(int step, vector<double>& vol, const double* s, int start, int end) const;

    // ------------------ data ------------------------------------------------------

    // log-normal vol
    string volType; // $unregistered

    CVolProcessedBSSP volLN1; // $unregistered

private:
    void generateDriftAndVols(const DateTimeArray& futurePathDates);
    
    /** keep the same structure as FD2DLN */  
    vector<DoubleArray>       vols; // $unregistered
    vector<CVolRequestLNSP> volRequests;    // [numAssets]  $unregistered
    vector<DoubleArray>       drifts;     /* x[numAssets] each with // $unregistered
                                              NbVolReq(iAsset) cols x
                                              NbSteps rows */
};

DRLIB_END_NAMESPACE

#endif
