//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : PastPathGenerator.hpp
//
//   Description : Monte Carlo path generator for historic 'pass'
//
//   Date        : Nov 6 2001
//
//
//----------------------------------------------------------------------------

#ifndef EDR_PASTPATHGENERATOR_HPP
#define EDR_PASTPATHGENERATOR_HPP
#include <map>
#include "edginc/MCPathConfig.hpp"
#include "MCProductClient.hpp"

#define STATE_VARIABLES

#ifdef STATE_VARIABLES
#include "edginc/SVGenSpot.hpp"
#endif


DRLIB_BEGIN_NAMESPACE
/** Dymmy past path generator, used for Sampras purpose */
class MCARLO_DLL DummyPastPathGenerator: virtual public MCPathGenerator{
public:
    virtual void            generatePath(int /*pathIdx*/) {};

    /** Returns true if this path generator is being used to 'simulate'
        the past */
    virtual bool            doingPast() const { return true; };
    
    /** Returns  1.0 */
    virtual double          maxDriftProduct(int /*iAsset*/) const { return 1.0; };

    /** to do: check that this is the number of simulated assets
        rather than the number of assets that the payoff is expecting */
    virtual int             NbSimAssets() const { return 0; };

    virtual const double*   Path(int /*iAsset*/, int /*iPath*/) const { return 0; };

    /** Returns the reference level for given asset and path */
    virtual double refLevel(int /*iAsset*/, int /*iPath*/) const { return 1.0; };

    /** Indicates which values in array returned by Path() is being
        simulated. Clients should use a loop of the form:
        "for (iStep = begin(); iStep < end(); iStep++)" */
    virtual int begin(int /*iAsset*/) const { return 0; };

    //// see  begin(int iAsset)
    virtual int end(int /*iAsset*/) const { return 0; };

    /** Returns true if there are historic simulation dates */
    virtual bool hasPast() const { return false; };

    virtual int getPathIndex() const { return 0; };

    //// Constructor
    DummyPastPathGenerator() {};

    /** Creates a state variable corresponding to a generator */
    virtual IStateVariableSP create(const IStateVariableGen* /*svGen*/) 
        { return IStateVariableSP(   ); };

private:
};


/** Class for dealing with past - hopefully we can reuse this for different
    generators. Supports different dates per asset */
class MCARLO_DLL PastPathGenerator: virtual public MCPathGenerator{
public:
    virtual void            generatePath(int pathIdx);

    /** Returns true if this path generator is being used to 'simulate'
        the past */
    virtual bool            doingPast() const;
    
    /** Returns  1.0 */
    virtual double          maxDriftProduct(int iAsset) const;

    /** to do: check that this is the number of simulated assets
        rather than the number of assets that the payoff is expecting */
    virtual int             NbSimAssets() const;

    virtual const double*   Path(int iAsset, int iPath) const;

    /** Returns the reference level for given asset and path */
    virtual double refLevel(int iAsset, int iPath) const;

    /** Indicates which values in array returned by Path() is being
        simulated. Clients should use a loop of the form:
        "for (iStep = begin(); iStep < end(); iStep++)" */
    virtual int begin(int iAsset) const;

    //// see  begin(int iAsset)
    virtual int end(int iAsset) const;

    /** Returns true if there are historic simulation dates */
    virtual bool hasPast() const;

    virtual int getPathIndex() const;

    //// Constructor
    PastPathGenerator(const IMCProduct*  prod);

    /** Creates a state variable corresponding to a generator */
    virtual IStateVariableSP create(const IStateVariableGen* svGen);

private:
    int                      numAssets;  
    bool                     simHasPast;
    DoubleArray              refLevels;  /* ref level per asset */
    vector<DoubleMatrix>     paths;  /* indexed by [numAssets][0][NbPastSteps]
                                        - this is to get a double*    */
};


#ifdef STATE_VARIABLES

class PastPathGenSpot;
typedef refCountPtr<PastPathGenSpot> PastPathGenSpotSP;

class PastPathGenQuadVar;
typedef refCountPtr<PastPathGenQuadVar> PastPathGenQuadVarSP;

class PastPathGenSqrtAnnualQuadVar;
typedef refCountPtr<PastPathGenSqrtAnnualQuadVar> PastPathGenSqrtAnnualQuadVarSP;

/** Class for dealing with past - hopefully we can reuse this for different
    generators. Supports different dates per asset */
class MCARLO_DLL PastPathGen: virtual public MCPathGenerator,
                   virtual public IMCStatelessGen,
                   virtual public IStateVariableGen::IStateGen {
public:
    virtual void            generatePath(int pathIdx);

    virtual void            advance();

    /** Returns true if this path generator is being used to 'simulate'
        the past */
    virtual bool            doingPast() const;
    
    /** Returns  1.0 */
    virtual double          maxDriftProduct(int iAsset) const;

    /** to do: check that this is the number of simulated assets
        rather than the number of assets that the payoff is expecting */
    virtual int             NbSimAssets() const;

    virtual const double*   Path(int iAsset, int iPath) const;

    /** Returns the reference level for given asset and path */
    virtual double refLevel(int iAsset, int iPath) const;

    /** Indicates which values in array returned by Path() is being
        simulated. Clients should use a loop of the form:
        "for (iStep = begin(); iStep < end(); iStep++)" */
    virtual int begin(int iAsset) const;

    //// see  begin(int iAsset)
    virtual int end(int iAsset) const;

    /** Returns true if there are historic simulation dates */
    virtual bool hasPast() const;

    virtual int getPathIndex() const;

    /** Constructor */
    PastPathGen(const MCProductClient* prodClient);

    /** Returns past path generator for spot */
    PastPathGenSpotSP getPathGenSpot() const;

    /** Returns past path generator for spot */
    PastPathGenQuadVarSP getPathGenQuadVar() const;

    /** Returns past path generator for spot */
    PastPathGenSqrtAnnualQuadVarSP getPathGenSqrtAnnualQuadVar() const;

    /** Returns the state variable corresponding to generator.
        Part of the IStateVariableGen::IStateGen IFace */
    virtual IStateVariableSP create(const IStateVariableGen* svGen);

private:
    DateTime                        today;
    StateVarDBase                   svDBase;                    //!< Collection of Generators + statevars
    PastPathGenSpotSP               pathGenSpot;                //!< Spot generator
    PastPathGenQuadVarSP            pathGenQuadVar;             //!< QuadVar generator
    PastPathGenSqrtAnnualQuadVarSP  pathGenSqrtAnnualQuadVar;   //!< SqrtAnnualQuadVar generator
};

typedef refCountPtr<PastPathGen> PastPathGenSP;


/** Class for dealing with past - hopefully we can reuse this for different
    generators. Supports different dates per asset */
class MCARLO_DLL PastPathGenSpot: virtual public MCPathGen {
public:
    /** Constructor */
    PastPathGenSpot(
        const MCProductClient* prodClient, 
        const SVGenSpotArray&     spotGenArray,
        StateVarDBase&         svDBase);

    /** Part of the MCPathGen interface */
    virtual void generatePath(int pathIdx);

    virtual void advance();

    /** Indicates whether there is past */
    bool hasPast() const;

    bool doingPast() const;

private:
    int                         numAssets;  
    bool                        simHasPast;
    vector<DoubleArray>         paths;          // [numAssets][NbPastSteps]
    // XXXX, giant hack here
    std::vector<MCPath::IStateVar*> spotSvSet;
};

/** Class for dealing with past - hopefully we can reuse this for different
    generators. Supports different dates per asset */
class MCARLO_DLL PastPathGenQuadVar: virtual public MCPathGen {
public:
    /** Constructor */
    PastPathGenQuadVar(
        const MCProductClient* prodClient, 
        const MCQuadVarArray& quadVarGenArray,
        StateVarDBase& svDBase);

    /** Part of the MCPathGen interface */
    virtual void generatePath(int pathIdx);

    /** Indicates whether there is past */
    bool hasPast() const;

    bool doingPast() const;

private:
    int                         numAssets;  
    bool                        simHasPast;
    vector<DoubleArray>         paths;          // [numAssets][NbPastSteps]
};

/** Class for dealing with past - hopefully we can reuse this for different
    generators. Supports different dates per asset */
class MCARLO_DLL PastPathGenSqrtAnnualQuadVar: virtual public MCPathGen {
public:
    /** Constructor */
    PastPathGenSqrtAnnualQuadVar(
        const MCProductClient*          prodClient, 
        const MCSqrtAnnualQuadVarArray& sqrtAnnualQuadVarGenArray,
        StateVarDBase&                  svDBase);

    /** Part of the MCPathGen interface */
    virtual void generatePath(int pathIdx);

    /** Indicates whether there is past */
    bool hasPast() const;

    bool doingPast() const;

private:
    int                         numAssets;  
    bool                        simHasPast;
    vector<DoubleArray>         paths;          // [numAssets][NbPastSteps]
};

#endif

DRLIB_END_NAMESPACE
#endif
