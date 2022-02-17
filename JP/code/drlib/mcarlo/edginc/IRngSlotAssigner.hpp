//----------------------------------------------------------------------------
//
//
//   Filename    : IRngSlotAssigner.hpp
//
//   Description : Interface to allow non-trivial assignment of random number slots for different assets.
//                 Mainly targets QSRM, where we want most of the assets to use the same random numbers even if new assets are aded.
//
//
//----------------------------------------------------------------------------

#ifndef IRNGSLOTASSIGNER_HPP
#define IRNGSLOTASSIGNER_HPP

#include "edginc/DECLARE.hpp"
#include <string>


DRLIB_BEGIN_NAMESPACE

using std::string;

/** Products may implement this interface (see SimpathIInstrumentMC) in order to guarantee random number slots for each asset.
 */

class MCARLO_DLL IRngSlotAssigner
{
public:
	virtual int getRngIdx(int globalIdx, const string& assetName) const = 0;
	virtual ~IRngSlotAssigner() {}
};

DECLARE_REF_COUNT(IRngSlotAssigner);

/** RngSlotAssignerLinear assigns factors linearly without gaps */
class RngSlotAssignerLinear : public IRngSlotAssigner
{
public:
	virtual int getRngIdx(int globalIdx, const string& /*assetName*/ ) const {
		return globalIdx;
	}
};

DRLIB_END_NAMESPACE
#endif


