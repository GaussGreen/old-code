//----------------------------------------------------------------------------
//
//   Description : Energy Future state variable for monte carlo 
//
//   Filename    : SVGenEnergyFuturePrice.cpp
//
//   Author      : Spyridon Schismenos
//
//   Date        : June 22, 2005
//
//-----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/SVGenEnergyFuturePrice.hpp"
//#include "edginc/MCPathConfigSRMGenSV.hpp"
#include "edginc/IElemStateVariableGenVisitor.hpp"

DRLIB_BEGIN_NAMESPACE

IStateVariableSP SVGenEnergyFuturePrice::create(IStateVariableSP oldStateVar,
                                             IStateVariableGen::IStateGen* pathGen) const
{
    return getEnergyFuturePriceSV(pathGen);
}

SVGenEnergyFuturePrice::IStateVarSP SVGenEnergyFuturePrice::getEnergyFuturePriceSV(IStateVariableGen::IStateGen* pathGen) const
{
    IStateVarSP EnergyFuturePriceSV(&dynamic_cast<SVGenEnergyFuturePrice::IStateVar&>(*pathGen->create(this)));
    return EnergyFuturePriceSV;
}


SVGenEnergyFuturePrice::IStateVarSP SVGenEnergyFuturePrice::pastSV() const
{
    return IStateVarSP(new SVGenEnergyFuturePrice::PastSV());
}

SVGenEnergyFuturePrice::SVGenEnergyFuturePrice(const DateTimeArray& dates,
										 const EnergyContractLabelArray& labels): 
							   priceDates(dates),contractLabels(labels){}


DateTimeArray SVGenEnergyFuturePrice::getPriceDates() const {
	return priceDates;
}
/** implementing 'visitor' model */
void SVGenEnergyFuturePrice::attachSVGen(IElemStateVariableGenVisitor* sv) const
{
    sv->processSVGen(this);
}

DRLIB_END_NAMESPACE


