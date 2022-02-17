#ifndef __drdatein__h
#define __drdatein__h

#include "drplatdep.h"
#include IOSTREAM_H

extern "C" {
#include "bastypes.h"
#include "busday.h"
}

//	This class wraps a TDateInterval.
//	My motivations are to get constructors and be able to name them in the 
//  global namespace as constants.  Thus, I can write stuff like
//	ONE_MONTH, etc anywhere.

class DRDateIn : public TDateAdjIntvl{
public:
	DRDateIn(int numPeriods, char prdType, TBoolean isBusDaysP, char* holidayFileP = "NONE");

	friend DRDateIn operator*(int, const DRDateIn&);

	friend ostream& operator<<(ostream&, const DRDateIn&);
};

// Some useful constants

const DRDateIn ONE_MONTH = DRDateIn(1, 'M', FALSE);
const DRDateIn ONE_YEAR = DRDateIn(1, 'A', FALSE);
const DRDateIn ONE_DAY = DRDateIn(1, 'D', FALSE);
const DRDateIn ONE_SEMI = DRDateIn (1, 'S', FALSE);
const DRDateIn ONE_QUARTER = DRDateIn(1, 'Q', FALSE);
const DRDateIn ONE_BUS_DAY = DRDateIn(1, 'D', TRUE);

#endif
