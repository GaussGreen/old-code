//----------------------------------------------------------------------------
//
//   Description : demo state variable for monte carlo 
//
//   Author      : Jay Z Wang
//
//-----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/SVGenDemo.hpp"
#include "edginc/MCPathConfigSRMGenSV.hpp"

DRLIB_BEGIN_NAMESPACE

IStateVariableSP SVGenDemo::create(
    IStateVariableSP oldStateVar,
    IStateVariableGen::IStateGen* pathGen) const{
    return getDemoStateVar(pathGen);
}

SVGenDemo::IStateVarSP SVGenDemo::getDemoStateVar(
    IStateVariableGen::IStateGen* pathGen) const
{
    IStateVarSP DemoStateVar(
        &dynamic_cast<SVGenDemo::IStateVar&>(*pathGen->create(this)));
    return DemoStateVar;
}


SVGenDemo::IStateVarSP SVGenDemo::pastSV() const
{
    return IStateVarSP(new SVGenDemo::PastSV());
}

SVGenDemo::INewStateVarSP SVGenDemo::getNewDemoStateVar(
    IStateVariableGen::IStateGen* pathGen) const
{
    INewStateVarSP DemoStateVar(
        &dynamic_cast<SVGenDemo::INewStateVar&>(*pathGen->create(this)));
    return DemoStateVar;
}


SVGenDemo::INewStateVarSP SVGenDemo::pastNewSV() const
{
    return INewStateVarSP(new SVGenDemo::PastNewSV());
}

/** implementing 'visitor' model */
void SVGenDemo::attachSVGen(IElemStateVariableGenVisitor* sv) const
{
    sv->processSVGen(this);
}


DRLIB_END_NAMESPACE


