//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : FD2DSVCJ.hpp
//
//   Description : two factor finite difference engine for SVCJ vol
//                    asset factor   dX = drift*dt + sqrt(V) dW + Merton jump in X
//                    vol factor  dV = k(V_0 - V)dt + vol*sqrt(V) dZ + Merton jump in V
//                    V_0 is allowed to be a gamma distribution of mean v0.
//
//   Author        Xiaolan Zhang
//
//   Date        : November 29, 2006
//
//----------------------------------------------------------------------------

#ifndef FD2DSVCJ_HPP
#define FD2DSVCJ_HPP
#include "edginc/FD2DSV.hpp"
#include "edginc/VolSVCJ.hpp"
#include "edginc/MarketDataFetcherLN.hpp"

DRLIB_BEGIN_NAMESPACE

class TREE_DLL FD2DSVCJ : public FD2DSV{
public:
    static CClassConstSP const TYPE;
    static void load(CClassSP& );
    static IObject* defaultFD2DSVCJ(){
        return new FD2DSVCJ();
    }
    
    friend class FD2DSolver;
    friend class FD2DSolverJumps;

    FD2DSVCJ( const CClassConstSP & type = TYPE );


//----------------------------------------------------------------

    virtual ~FD2DSVCJ();

    /** Create a MarketDataFetcher which will be used for retrieving market data etc */
    MarketDataFetcherSP createMDF() const;
        
    /**
     * Whether to enable RiskMapping when computing sensitivities for
     * instruments priced using this model
     *
     * Returns allowRiskMapping ? riskMappingAllowed : riskMappingDisallowed.
     * See IModel::wantsRiskMapping().
     */

    virtual IModel::WantsRiskMapping wantsRiskMapping() const;


    /** get the initial conditions and corresponding indexes in the slices
        for forward induction */
    //virtual void getInitialConditions(IntArray& initialIndex, DoubleArray& initialValue) const;

protected:

    /** retrieving market data from MDF */
    virtual void retrieveFactor();

    
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

//    virtual  void setFdBounds(double& volForBound1, double alpha1, double& outLowB1, double& outUpB1,
//                              double& volForBound2, double alpha2, double& outLowB2, double& outUpB2);
    
    //virtual   void SetSpaceGridData();
    
    /**----------------------------------------------------------------
       setup the coefficients of the PIDE needed to be solved
       Assume the most general PDE is 
       U_t + a*U_x + b*U_y + c*U_xx + d*U_yy + e*U_xy + f*U + Jump(x,y) = 0;
       this ft is called at each time step, so they're time depd.
       ----------------------------------------------------------------*/
    
    /** assuming that DimX (vertical) is stock and DimY is Variance*/
    
    void pdeCoeff(int step, double*** coeff,  
                  int bot1, int top1, int bot2, int top2 );
    
    /**-------------------------------------
       for SVCJ, not in the base FD2D model 
       -------------------------------------*/
    //jumps parts
    double jumpProDen(double x, double y);
    //temporary, to remove
    double getJumpSizeMean(double y);
    double getJumeSizeVar();

    
    VolSVCJConstSP          volSVCJ; // $unregistered

    // transient fields
    double                   cCrashRatedt; //commonCrashRate*dt

    //temporary, only for test purpose.

    int whichM; 

};

typedef smartPtr<FD2DSVCJ> FD2DSVCJSP;

DRLIB_END_NAMESPACE
#endif
