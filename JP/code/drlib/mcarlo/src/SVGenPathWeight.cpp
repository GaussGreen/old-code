//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : SVGenPathWeight.cpp
//
//   Description : A Generator of MC Path Weight State Variables
//                 (SV as they can develop to become instrument specific)
//
//   Date        : 22 Nov 2006
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/SVGenPathWeight.hpp"
#include "edginc/IElemStateVariableGenVisitor.hpp"
#include "edginc/IQMCStateVariableBase.hpp"

DRLIB_BEGIN_NAMESPACE

SVGenPathWeight::SVGenPathWeight(const DateTimeArray& _dates) : dates(_dates)
{}


/** Create the corresponding State Variable for this State
    Variable Generator (from IStateVariableGen interface). The
    previous IStateVariableSP (may be null) should be passed in.  The
    return object may or may not be the same as oldStateVar. */
IStateVariableSP SVGenPathWeight::create(IStateVariableSP             oldStateVar,
                                   IStateVariableGen::IStateGen* pathGen) const{
    return getSVPathWeight(pathGen);
}


/** Returns a MC Path Weight state variable which then
    provides access to the path etc. */
SVPathWeightSP SVGenPathWeight::getSVPathWeight(IStateVariableGen::IStateGen* pathGen) const
{
    SVPathWeightSP pathWeightSV(&dynamic_cast<SVPathWeight&>(*pathGen->create(this)));
    return pathWeightSV;
}

class SVGenPathWeight::DeterminsticSV: public SVPathWeight, public IQMCStateVariableSpot {
public:

    /** Returns true if this state variable is being used to 'simulate'
        the past. This is a useful method for users of state variables - 
        need to see how hard it is to implement */
    virtual bool doingPast() const{
        return doingThePast;
    }
    /** All elements are identically 1.0 for the past */
    virtual double element(int idx) const {
        return 1.0;
    }

    virtual double getWeight(int i)const { return element(i);}
    virtual void prepare(bool) {}

    DeterminsticSV(bool                 doingThePast):
    IQMCStateVariableSpot(DateTimeArray()),
    doingThePast(doingThePast)
    {}

private:
    bool           doingThePast;
};

/** For use by Path Generators that want to use deterministic rates. */
SVPathWeight* SVGenPathWeight::determinsticSV(bool doingPast) const{
    return new DeterminsticSV(doingPast);
}

/** implementing 'visitor' model */
void SVGenPathWeight::attachSVGen(IElemStateVariableGenVisitor* sv) const
{
    sv->processSVGen(this);
}

DRLIB_END_NAMESPACE
