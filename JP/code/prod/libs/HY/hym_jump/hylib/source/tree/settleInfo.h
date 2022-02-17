// settleInfo.h: interface for the SettlementInfo class.
//
//////////////////////////////////////////////////////////////////////

#ifndef	_SETTLEINFO_H_
#define	_SETTLEINFO_H_

#include "kdate.h"  

class SettleInfo  
{
	KDate m_settlementDate;

public:
	SettleInfo();
	SettleInfo(KDate settlementDate);

	void set_settlementDate(KDate settleDate);
	KDate get_settlementDate() const;

};

#endif
