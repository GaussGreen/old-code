//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : MCPathGenerator.hpp
//
//   Description : Monte Carlo interface
//
//   Author      : Oliver Brockhaus
//
//   Date        : April 2003
//
//
//----------------------------------------------------------------------------

#ifndef EDR_MCPATHGENERATOR_HPP
#define EDR_MCPATHGENERATOR_HPP

#include "edginc/StateVariableGen.hpp"

DRLIB_BEGIN_NAMESPACE

//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
class IMCStatelessGen {
public:
    virtual void advance() = 0;
};

//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
class MCSliceGen {
public:
    virtual void generateSlice(int sliceIdx) = 0;
    virtual bool doingPast() const = 0;
    virtual ~MCSliceGen() {};
};

class MCStatelessSliceGen : public MCSliceGen, 
                            virtual public IMCStatelessGen {
};

/////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////


/** New interface for state variables approach. The path generator, 
    from the Monte Carlo point of view, can generate paths. */
class MCPathGen {
public:
    /** Requests path generator to create a new path (using new random
        variables etc). The monte carlo engine invokes this before each
        call the product's payoff function. The pathIdx parameter 
        indicates which path should be generated - this is only relevant
        if storePaths() has been invoked */
    virtual void generatePath(int pathIdx) = 0;

    /** Returns true if this path generator is being used to 'simulate'
        the past */
    virtual bool doingPast() const = 0;

    /** Destructor */
    virtual ~MCPathGen() {};
};

class  MCStatelessPathGen : virtual public MCPathGen,
                            virtual public IMCStatelessGen {
};

////////////////////////////////////////////////////////////////////

/** Pre state variables approach: Path generator can generate
    paths AND provide access to state variables i.e. spot paths
    and reflevels. */
class MCPathGenerator: virtual public MCPathGen,
                       virtual public IStateVariableGen::IStateGen {
public:
    /** Returns the number of simulated assets */
    virtual int             NbSimAssets() const = 0;

    /** returns the path for the given asset and path index. Use
        begin() and end() to see which values are being simulated. Values
        before begin() are valid and will be historic values */
    virtual const double*   Path(int iAsset, int iPath) const = 0; 

    /** Returns the reference level for given asset and path. For
        the PastPathGenerator then any future average values will
        be 'estimated' (current methodology is current spot). For
        a future path generator, then this must return 'estimated'
        values too (eg same values as the Past generator) until
        generatePath is invoked. At which point it will contain
        the correct value for the current path */
    virtual double          refLevel(int iAsset, int iPath) const = 0;

    /** Returns the MAX(1,d_0)*MAX(1,d_1)*...*MAX(1,d_n) where d_i is
        the drift between simulation date i and simulation date i+1.
        (With d_0 = simulation start date). Range of dates is restricted to
        those in the future so the past path generator will
        return 1.0 */
    virtual double          maxDriftProduct(int iAsset) const = 0;

    /** Indicates which values in array returned by Path() is being
        simulated. Clients should use a loop of the form:
        "for (iStep = begin(); iStep < end(); iStep++)" */
    virtual int             begin(int iAsset) const = 0;

    //// see  begin(int iAsset)
    virtual int             end(int iAsset) const = 0;

    /** Returns true if there are historic simulation dates */
    virtual bool            hasPast() const = 0;

    /** Returns true if this path generator is being used to 'simulate'
        the past */
    virtual bool            doingPast() const = 0;

    /** Requests path generator to create a new path (using new random
        variables etc). The monte carlo engine invokes this before each
        call the product's payoff function. The pathIdx parameter 
        indicates which path should be generated - this is only relevant
        if storePaths() has been invoked */
    virtual void            generatePath(int pathIdx) = 0;

    /** Provide access to the current 'pathIdx' value as per generatePath() */
    virtual int             getPathIndex() const = 0;

    virtual ~MCPathGenerator() {};
};

typedef MCPathGenerator IMCPathGenerator;
typedef MCPathGenerator IPathGenerator;
typedef refCountPtr<MCPathGenerator> MCPathGeneratorSP;

DRLIB_END_NAMESPACE

#endif



