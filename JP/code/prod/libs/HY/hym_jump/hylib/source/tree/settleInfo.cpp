// SettlementInfo.cpp: implementation of the SettlementInfo class.
//
//////////////////////////////////////////////////////////////////////

#include "settleInfo.h"
#include "kdate.h"  


SettleInfo::SettleInfo(){}

SettleInfo::SettleInfo(KDate settlementDate):m_settlementDate(settlementDate){}


void SettleInfo::set_settlementDate(KDate settlementDate)
{
	m_settlementDate = settlementDate;
}


KDate SettleInfo::get_settlementDate() const{return m_settlementDate;}
