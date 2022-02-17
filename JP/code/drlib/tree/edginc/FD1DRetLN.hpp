//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : FD1DRetLN.hpp
//
//   Description : one factor finite difference engine for 1 LN
//                    asset factors   dX_i = drift_i*dt + sigma_i * dW_i with i=1,
//
//   Author      : Xiaolan Zhang
//
//   Date        : Feb 14, 2006
//
//----------------------------------------------------------------------------

#ifndef FD1D_RET_LN_HPP
#define FD1D_RET_LN_HPP
#include "edginc/FD1DRet.hpp"
#include "edginc/MarketDataFetcherLN.hpp"
#include "edginc/MultiFactors.hpp"

#include "edginc/VolProcessedBS.hpp"

DRLIB_BEGIN_NAMESPACE

class TREE_DLL FD1DRetLN : public FD1DRet{
public:
    static CClassConstSP const TYPE;
    static void load(CClassSP& );
    static IObject* defaultFD1DRetLN(){
        return new FD1DRetLN();
    }
    
    friend class FD1DRetSolver;

    FD1DRetLN();
    virtual ~FD1DRetLN();

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
    
    virtual  void setFdBounds(double& volForBound1, double alpha1, double& outLowB1, double& outUpB1);

    /**-----------------------------*/
    /** update payoffIndex          */
    /** ----------------------------*/
   
    virtual void payoffIndexUpdate (int& step, FDProduct::UpdateType);
    
    /**----------------------------------------------------------------
       setup the coefficients of the PIDE needed to be solved
       Assume the most general PDE is 
       U_t + a*U_x + b*U_y + c*U_xx + d*U_yy + e*U_xy + f*U + Jump(x,y) = 0;
       this ft is called at each time step, so they're time depd.
       ----------------------------------------------------------------*/
    
    /** assuming that DimX (vertical) is stock and DimY is Variance*/
    
    void pdeCoeff(int step, double** coeff, int bot1, int top1);

    /**----------------------------------------------------------------
        calculate vol for current step for a set of spot levels
        returns number of vol calculated - one (flat for all node) or num 
    ----------------------------------------------------------------*/

    int GetStepVol(int step, vector<double>& vol, const double* s, int start, int end);

    /** override FDModel::makeProduct to retrieve spot index product */
    virtual FDProductSP makeProduct( const IProdCreatorSP & creator );

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

typedef smartPtr<FD1DRetLN> FD1DRetLNSP;

DRLIB_END_NAMESPACE
#endif
