//----------------------------------------------------------------------------
//
//   Group       : QR Equities London
//
//   Filename    : VanillaGridMulti.cpp
//
//   Description : Multi Factor Version of VanillaGrid
//
//----------------------------------------------------------------------------
#include "edginc/config.hpp"
#include "edginc/GenericNFBase.hpp"
#include "edginc/LinearStrikeTSVolRequest.hpp"
#include "edginc/MonteCarlo.hpp"
#include "edginc/IAggregate.hpp"

DRLIB_BEGIN_NAMESPACE

typedef array<IAggregateMakerSP, IAggregateMaker> IAggregateMakerArray;
typedef smartPtr<IAggregateMakerArray> IAggregateMakerArraySP;

typedef array<IDoubleArrayModifierMakerSP, IDoubleArrayModifierMaker> IDoubleArrayModifierMakerArray;
typedef smartPtr<IDoubleArrayModifierMakerArray> IDoubleArrayModifierMakerArraySP;

/** forward declaration */
class VanillaGridMultiInstrumentCollection;

/** VanillaGridMulti product */
class PRODUCTS_DLL VanillaGridMulti: public GenericNFBase, 
                                     virtual public IMCIntoProduct {

    friend class VanillaGridMultiInstrumentCollection;

protected:
    /** registered fields */ 
    DateTimeArray                           avgOutDates;
    IAggregateMakerArraySP                  baskets;         // specification of weighting for all dates
    IDoubleArrayModifierMakerArraySP        perfTypes;    // specification of performance for each date

    /** transient field */
    int nbPerfPerMat;

public:
    static CClassConstSP const TYPE;
    friend class VanillaGridMultiMC;

    /** validation */
    virtual void validatePop2Object(); 
   
    /** mandatory, since derivation from GenericNFBase */ 
    virtual const DateTimeArray samplingDates() const; 

    /** implementation of MonteCarlo::IntoProduct interface */
    virtual IMCProduct* createProduct(const MonteCarlo* model) const; 

    /** additional constructor */
    VanillaGridMulti(VanillaGridMulti*                  grid,
                     const DateTimeArray&               avgOutDate,
                     IAggregateMakerArraySP             basket,
                     IDoubleArrayModifierMakerArraySP   perfType);    

private:
    VanillaGridMulti(); 
    
    static IObject* defaultVanillaGridMulti(); 

    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz); 
};

typedef smartPtr<VanillaGridMulti> VanillaGridMultiSP;
typedef smartConstPtr<VanillaGridMulti> VanillaGridMultiConstSP;


/** MC product class for VanillaGridMulti */
class PRODUCTS_DLL VanillaGridMultiMC : public IMCProduct,
                                        virtual public IMCProductLN {
private:
    const VanillaGridMulti*         inst;           // reference to original instrument
    
    int                             nbAssets;       // nb of assets    
    CDoubleMatrixSP                 payoffMatrix;   // nb rows = nb maturities, nb = nb 
    CDoubleMatrixSP                 payoffMatrixSq; // for calc of std error
    CDoubleMatrixSP                 payoffMatrixSqHelper;   // for calc of std error

    DoubleArray				        refLevel;	    // refLevel for each asset
    DoubleArray                     refLevelSoFar;
    
    vector<IAggregateSP>            basket;
    SimpleDoubleArray		        basketHelper;
    vector<IDoubleArrayModifierSP>  performance;
    TrivialDoubleArray              performanceHelper;  

    double                          countPaths;     // manual counting of nb of paths
    int                             futureBeginIdx; //

public:
    
    /** Use this opportunity to do any LogNormal driven initialisation
        of the instrument before the main MC loop. e.g closed form barrier adjustment */
    void initialiseLN(const  MCPathGenerator*  pathGen) const; 

    /** equivalent to InstIntoMCProduct */
    VanillaGridMultiMC(const VanillaGridMulti*  inst,
                       const SimSeriesSP&       simSeries); 

    /** Called within the simulation loop */
    void payoff(const MCPathGenerator*  pathGen,
                IMCPrices&                prices);

    /** invoked after final simulated path is run. */
    virtual void recordExtraOutput(CControl*     control,
                                   Results*      results,
                                   const IMCPrices& prices) const;

    /** need interpolation level for the LogNormal path generator */ 
    CVolRequestLNArray getVolInterp(const MCPathGenerator* pathGen,
                                    int                     iAsset) const; 
};

DRLIB_END_NAMESPACE