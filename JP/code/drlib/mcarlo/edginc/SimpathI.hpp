//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : SimpathI.hpp
//
//   Description : Sampras Monte Carlo model
//
//
//
//----------------------------------------------------------------------------

#ifndef EDR_SIMPATHI_HPP
#define EDR_SIMPATHI_HPP
#include "edginc/FORWARD_DECLARE.hpp"
#include "edginc/Class.hpp"
//#include "edginc/Instrument.hpp"
#include "edginc/Model.hpp"
#include "edginc/RegressionTest.hpp" // to enable testing via models.exe
#include "edginc/ClientRunnable.hpp"

//#include "edginc/MarketDataFetcher.hpp" // MarketDataFetcherSP
//#include "edginc/MCPathGenerator.hpp" // MCPathGenerator.hpp
#include "edginc/Control_forward.hpp"
#include "edginc/Results_forward.hpp"
#include "edginc/MarketData_forward.hpp"

DRLIB_BEGIN_NAMESPACE

class IMCPathConfig;
class CInstrument;

FORWARD_DECLARE(IMCPathConfig);
FORWARD_DECLARE(MarketDataFetcher);
FORWARD_DECLARE_REF_COUNT(MCPathGenerator);
FORWARD_DECLARE_REF_COUNT(IMCProduct);
FORWARD_DECLARE_REF_COUNT(IMCPricesGeneric);

class MCARLO_DLL SimpathI : public CModel, 
                            public virtual IRegressionTest,
                            public virtual ClientRunnable
{
public:
    static CClassConstSP const TYPE;
    virtual ~SimpathI();

    // override the default deep copy semantics
    virtual IObject* clone() const { return const_cast<SimpathI*>(this); }

    // return all state var interfaces
    bool Initialize();
    
    IObjectSP GeneratePath(int idx); // return prices from 1 path
    IObjectSP GeneratePath(int idx, IMCPrices& prices); // return and update prices
    
    
    virtual void Price(CInstrument*  instrument,
                       CControl*     control,
                       CResults*     results)
    {
        throw ModelException("Price called on SimpathI model!");
    }

    /**
     * Whether to enable RiskMapping when computing sensitivities for
     * instruments priced using this model
     *
     * See IModel::wantsRiskMapping().
     */

    virtual IModel::WantsRiskMapping wantsRiskMapping() const {
        return riskMappingIrrelevant;
    }

    int getNbIter(); // access nbIter
    virtual IObjectSP runTest() const; // Implements IRegressionTest interface
    virtual IObjectSP run() {
        return runTest();
    }
protected:
    SimpathI(CClassConstSP clazz);
    MarketDataFetcherSP createMDF() const;

    // market
    CMarketDataSP market;

    // model
    IMCPathConfigSP pathConfig;
    MCPathGeneratorSP pathGen; // $unregistered

    // instrument
    CInstrumentSP instrument;
    IMCProduct* instrumentMC; // $unregistered

    int nbIter;

private:
    SimpathI();
    SimpathI(const SimpathI& rhs);
    SimpathI& operator=(const SimpathI& rhs);
    static IObject* defaultSimpathI();
    static void load(CClassSP& clazz);
public:
	void serialize(string fn);
};

DECLARE(SimpathI);

//////////////////////////////////////////////////////
class MCARLO_DLL SimpathIController : public CObject {
public:
    static CClassConstSP const TYPE;
    virtual ~SimpathIController() {};

    bool Initialize();

    IObjectSP GeneratePath();

	IObjectSP runTest();

protected:
    SimpathIController(CClassConstSP clazz);

    SimpathI* model;
    int idx;

private:
    SimpathIController();
    static IObject* defaultSimpathIController();
    static void load(CClassSP& clazz);
};


//////////////////////////////////////////////////////
class MCARLO_DLL SimpathIShell : public CObject,
                                 public virtual ClientRunnable
{
public:
    static CClassConstSP const TYPE;
    virtual ~SimpathIShell() {};

    IObjectSP run();

protected:
    SimpathIShell(CClassConstSP clazz);

    SimpathI* simpathi;

private:
    SimpathIShell();
    static IObject* defaultSimpathIShell();
    static void load(CClassSP& clazz);
};


DRLIB_END_NAMESPACE


#endif



