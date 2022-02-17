//----------------------------------------------------------------------------
//
//   Group       : GCCG Derivatives Research
//
//   Filename    : EnergyFixing.hpp
//
//   Description : Energy Fixing info
//
//   Author      : Sean Chen
//
//   Date        : Aug. 29, 2005
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/EnergyFixing.hpp"

DRLIB_BEGIN_NAMESPACE

EnergyFixing::EnergyFixing() : CObject(TYPE), rate(0.0), nearbyRel(0)
{
}

EnergyFixing::~EnergyFixing()
{
}

EnergyFixing::EnergyFixing(const EnergyFixing& v) : CObject(TYPE)
{
	*this = v;
}

EnergyFixing& EnergyFixing::operator=(const EnergyFixing& v)
{
	if (&v!=this)
	{
		energyUnderlyer = v.energyUnderlyer;
		rate = v.rate;
		nearbyRel = v.nearbyRel;
		fixingDate = v.fixingDate;
		contractLabelForFixing = v.contractLabelForFixing;
	}

	return *this;
}


void EnergyFixing::setEnergyUnderlyer(const EnergyUnderlyerSP& underlyer)
{
	energyUnderlyer = underlyer;
}

void EnergyFixing::setFixingRate(double fixingRate)
{
	rate = fixingRate;
}

void EnergyFixing::setNearbyRel(int nearby)
{
	nearbyRel = nearby;
}

// date of contract is labeled by its expiry date
void EnergyFixing::setFixingDate(const DateTime& theFixingDate)
{
	fixingDate = theFixingDate;
}

void EnergyFixing::setContractLabelForFixing(const EnergyContractLabel& label)
{
		contractLabelForFixing = label;
}

EnergyUnderlyerSP EnergyFixing::getEnergyUnderlyer() const
{
	return energyUnderlyer;
}

double EnergyFixing::getFixingRate() const
{
	return rate;
}

int EnergyFixing::getNearbyRel() const
{
	return nearbyRel;
}

DateTime EnergyFixing::getFixingDate() const
{
	return fixingDate;
}

EnergyContractLabel EnergyFixing::getContractLabelForFixing() const
{
	return contractLabelForFixing;
}

/*******
int EnergyFixing::GetFixingContractNumber(DRDate today) const
{
	return energyUnderlyer->calculateContractNumber(today, fFixingLabel, DRCommodityIndex::kIndex) + nearbyRel;
}
***/

// resolve label into the expiry date contract - but can't apply nearby rules
DateTime EnergyFixing::getContractDateForFixing() const
{
	return energyUnderlyer->expiryDate(contractLabelForFixing, EnergyUnderlyer::INDEX);
}

CClassConstSP const EnergyFixing::TYPE = CClass::registerClassLoadMethod(
    "EnergyFixing", typeid(EnergyFixing),load);

bool  EnergyFixingLoad() { return (EnergyFixing::TYPE != 0);   }

DRLIB_END_NAMESPACE

