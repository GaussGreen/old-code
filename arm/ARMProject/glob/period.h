/*!
 * Copyright (c) NATIXIS  April 2007 Paris
 *
 *	\file period.h
 *	\author  Francois Poitou
 */

#ifndef _PERIOD_H
#define _PERIOD_H

/// gpbase

/// kernel
/*
	reminder : itsUnit can be 
	- K_DAILY
	- K_ANNUAL
	- K_SEMIANNUAL
	- K_QUARTERLY
	- K_BIMONTHLY
	- K_MONTHLY
	- K_WEEKLY
	- K_DAILY


*/

//CC_BEGIN_NAMESPACE( ARM )

/*	inline int IS_ZERO(double x) 
	{
		if ( fabs(x) < 1e-8 )
		{
		   return(1);
		}
		else
		{
		   return(0);
		}
	}
*/
	class Period /* : public ?? */{
		public:
			Period(): itsLength(0), itsUnit(K_DAILY)
			{}
			Period(int length, int unit): itsLength(length), itsUnit(unit)
			{}
			Period(const Period& per): itsLength(per.itsLength), itsUnit(per.itsUnit)
			{}
			virtual ~Period()
			{}

			int GetLength() const {return itsLength ;}
			int GetUnit() const {return itsUnit ;}
		private :
			int itsLength ;
			int itsUnit ;

	};

//CC_END_NAMESPACE()

#endif