#ifndef SVGenDemo_HPP
#define SVGenDemo_HPP

#include "edginc/ElemStateVariableGen.hpp"


DRLIB_BEGIN_NAMESPACE

class MCARLO_DLL SVGenDemo: virtual public IElemStateVariableGen,
                                 public virtual VirtualDestructorBase
{
public:

    /** Interface for the state variable that SVGenExpectedSpot produces. This is
        the type that products deal with in the payoff. The payoff obtains
        it by calling the getExpSpotSV() method below. Note support here
        for only one asset per state variable (unlike SVGenSpot). May
        need to reconsider at some point */
    class MCARLO_DLL IStateVar: public virtual IStateVariable{
    public:
        virtual ~IStateVar() {};
        virtual double* getSpot() = 0;
    };
    typedef smartPtr<IStateVar> IStateVarSP;

    class MCARLO_DLL PastSV : public IStateVar {
    public:
        double* getSpot() { return 0; }
        bool doingPast() const { return true; }
    };

    // Is this obsolete?  It is not being called anywhere...
    virtual IStateVariableSP create(
        IStateVariableSP oldStateVar,
        IStateVariableGen::IStateGen* pathGen) const;

    IStateVarSP getDemoStateVar(IStateVariableGen::IStateGen* pathGen) const;

    // used with past path generator
    IStateVarSP pastSV() const;


    ////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////
    class MCARLO_DLL INewStateVar: public virtual IStateVariable {
    public:
        virtual ~INewStateVar() {};
        virtual double getSpot() = 0;
    };
    typedef smartPtr<INewStateVar> INewStateVarSP;

    class MCARLO_DLL PastNewSV : public INewStateVar {
    public:
        double getSpot() { return 0; }
        bool doingPast() const { return true; }
    };

    INewStateVarSP getNewDemoStateVar(IStateVariableGen::IStateGen* pathGen) const;

    INewStateVarSP pastNewSV() const;

    /** Constructor */
    SVGenDemo() {};
    void attachSVGen(IElemStateVariableGenVisitor*) const;


private:
};

typedef smartPtr<SVGenDemo> SVGenDemoSP;

DRLIB_END_NAMESPACE

#endif
