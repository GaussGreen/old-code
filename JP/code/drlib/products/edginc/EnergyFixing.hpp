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

#ifndef _ENERGYFIXING_H_
#define _ENERGYFIXING_H_

#include "edginc/config.hpp"
#include "edginc/Class.hpp"
#include "edginc/Object.hpp"
#include "edginc/EnergyUnderlyer.hpp"
#include "edginc/EnergyContractLabel.hpp"

#include <string>
using namespace std;

DRLIB_BEGIN_NAMESPACE

class PRODUCTS_DLL EnergyFixing: public CObject
{
    public:

		static CClassConstSP const TYPE ;

        EnergyFixing();
        ~EnergyFixing();
        EnergyFixing(const EnergyFixing& v);
        EnergyFixing& operator=(const EnergyFixing& v);

        // sets
        void setEnergyUnderlyer(const EnergyUnderlyerSP& underlyer);
        void setNearbyRel(int nearbyRel);
        void setFixingDate(const DateTime& contractExpiryDate);
		void setFixingRate(double fixingRate);
        void setContractLabelForFixing(const EnergyContractLabel& label);

        // gets
        EnergyUnderlyerSP getEnergyUnderlyer() const;
        int getNearbyRel() const;
        DateTime getFixingDate() const;
		double getFixingRate() const;
        EnergyContractLabel getContractLabelForFixing() const;
        bool isFixingLabelSet() const;

		DateTime getContractDateForFixing() const;

		static void load(CClassSP& classToLoad)
		{
           // empty for now
		}
        
    private:

        EnergyUnderlyerSP            energyUnderlyer; // $unregistered
		double                     rate; // $unregistered
        int                        nearbyRel;  // relative to current contract. $unregistered
        
        DateTime                   fixingDate; // $unregistered
        EnergyContractLabel        contractLabelForFixing; // $unregistered
        
};

    
DRLIB_END_NAMESPACE

#endif
