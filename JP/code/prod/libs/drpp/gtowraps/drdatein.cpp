#include "drdatein.h"

ostream& operator<<(ostream& s, const DRDateIn& a)
{
	s << "(" << a.interval.prd << "," ;
	s << a.interval.prd_typ << "," ;
	s << a.isBusDays << "," ;
	s << a.holidayFile << ")";
	return s;
}

DRDateIn::DRDateIn (int numPeriods, char prdType, TBoolean isBusDaysP, char* holidayFileP)
{
        SET_TDATE_INTERVAL (interval, numPeriods, prdType);
        isBusDays = isBusDaysP;
        holidayFile = holidayFileP;
    badDayConv = GTO_BAD_DAY_NONE;
}

DRDateIn operator*(int num, const DRDateIn& a)
{
        int newPrds = num * a.interval.prd;

        return DRDateIn(newPrds, a.interval.prd_typ, a.isBusDays, a.holidayFile);
}
