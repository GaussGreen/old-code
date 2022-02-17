/****************************************************************************/
/*      Declaration file for structure of data.                             */
/****************************************************************************/

#ifndef	_DDMAP_H_
#define	_DDMAP_H_

#include "kplatform.h"
#include "General/General.h"
#include "bastypes.h"        /* TCouponDates */
#include "cdate.h"           /* TDate TDateInterval */
#include "fltrate.h"  /* Float rate */
#include "kdate.h"   






class DDMap:public std::map<KDate,double>, public CM::Object
{

public:
	DDMap(){}
	DDMap(std::string spec, KDate endDate);
	DDMap(int numPoints, KDate *dates, double *rates);
	// get amount at date
	double get_amount(KDate date) const;
	//get amount between date1 and date2
	double get_amount(KDate date1, KDate date2) const;
	double get_average_amount(KDate date1, KDate date2) const;
	DDMap  get_subMap(KDate date1, KDate date2) const;
	KVector(KDate) get_datelist(KDate date1=0, KDate date2=0) const;
};

#endif
