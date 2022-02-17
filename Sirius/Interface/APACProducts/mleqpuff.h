// MlEqPuff.h: interface for the MlEqPuff class.
//
//////////////////////////////////////////////////////////////////////

#ifndef _MLEQPUFF_H
#define _MLEQPUFF_H

#pragma once


#include "MlEqObjects.h"
#include "smart.h"

class MlEqPuff : public RCObject
{
public:
	MlEqPuff();

	double								GetStrike() const;
	void								PutStrike(double f);
	long								GetMaturity() const;
	void								PutMaturity(long n);
	double								GetMultiplier() const;
	void								PutMultiplier(double f);
	MlEqAssetHandle						GetUnderlying() const;
	void								PutUnderlying(MlEqAssetHandle h);

	double								GetPrice(void) const;

protected:
	MlEqAssetHandle						m_hUnderlying;
	double								m_fMultiplier;
	double								m_fStrike;
	long								m_nDate;
};

#endif
