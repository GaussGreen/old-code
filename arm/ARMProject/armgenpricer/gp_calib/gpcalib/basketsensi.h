/*!
 *
 * Copyright (c) CDC IXIS CM November 2003 Paris
 *
 *	\file basket.h
 *
 *  \brief base class for basket calibration
 *
 *	\author  J. Messines
 *	\version 1.0
 *	\date August 2006
 */


#ifndef _INGPCALIB_BASKETSENSI_H
#define _INGPCALIB_BASKETSENSI_H

/// this header has to come first
/// as env.h defines pre-processor constants for checking and validation
/// namely for std::vector<double> and ARM_Matrix

/// gpbase
#include "gpbase/env.h"
#include "gpbase/port.h"
#include "gpbase/typedef.h"
#include "gpbase/countedptr.h"

#include "gpbase/gpvector.h"
#include "gpbase/gplinalgtypedef.h"

/// kernel
#include <glob/dates.h>

class ARM_ZeroCurve;
class ARM_SpreadOption;
class ARM_SwapLeg;


CC_BEGIN_NAMESPACE( ARM )

class ARM_STRIPPER;

struct ARM_VanillaSpreadOptionArg;

class ARM_BasketSensi
{
private:
	double				itsMktPrice;
	double				itsWeight;
	std::deque<bool> /* FIXMEFRED: mig.vc8 (23/05/2007 14:27:19): dont use std::vector<bool>. see item 18 in "Effective STL" from Scott Meyers*/		itsRompu;
protected:
	ARM_IntVector		itsCpnIndex;
	std::vector<double>		itsPayDates;
	std::vector<double>		itsNotional;
	std::vector<double>		itsResetDates;
	std::vector<double>		itsStartDates;
	std::vector<double>		itsDF;
	std::vector<double>		itsFees;

	ARM_Date			itsEndDate;
	ARM_Date			itsRealEndDate;

public:
	ARM_BasketSensi(double w, double p)
		: itsMktPrice(p),itsWeight(w){};
	~ARM_BasketSensi(){};

	virtual double GetCpn(int i) = 0;
	virtual double GetDCpn(int i, const ARM_STRIPPER& stripper) = 0;
	virtual void Init(int i, const ARM_STRIPPER& stripper)=0;

	int GetCpnSize()				{return itsResetDates.size();}
	double Notional(int i)			{return itsNotional[i];}
	double GetPay(int i)			{return itsPayDates[i];}
	double GetDF(int i)				{return itsDF[i];}
	double GetW()					{return itsWeight;}
	
	const ARM_Date& GetEndDate()		{return itsEndDate;};
	const ARM_Date& GetRealEndDate()	{return itsRealEndDate;};

	virtual void FillCpnIndex(const std::vector<double>& exerDates,const std::vector<double>& startDates);
	int GetCpnIndex(int exerIndex)	{return itsCpnIndex[exerIndex];}
	
	void SetRompu(int exerIndex, bool status)	{itsRompu[exerIndex]=status;};
	bool IsRompu(int exerIndex)					{return itsRompu[exerIndex];};
	
	void SetFee(int exerIndex, double fee)		{itsFees[exerIndex]=fee;};
	double GetFee(int exerIndex)				{return itsFees[exerIndex];}
	virtual double ComputeFee(int exerIndex,double startDate);
	
};

class ARM_BasketSensiSO : public ARM_BasketSensi
{
private:
	ARM_VanillaSpreadOptionArg*	itsArg;

	std::vector<double> itsCpnValue;
	std::vector<double> itsCpnLongSensi;
	std::vector<double> itsCpnShortSensi;

	std::vector<double> itsCpnCoverages;

	std::vector<double> itsLong0;
	std::vector<double> itsShort0;

public:
	ARM_BasketSensiSO(ARM_ZeroCurve* curve,ARM_SpreadOption* so,double w, double p);
	~ARM_BasketSensiSO(){};

	virtual double	GetCpn(int i);
	virtual double	GetDCpn(int i, const ARM_STRIPPER& stripper);
	virtual void	Init(int i, const ARM_STRIPPER& stripper);
};

class ARM_BasketSensiCSO : public ARM_BasketSensi
{
private:
	ARM_VanillaSpreadOptionArg*	itsArg;

	std::vector<double> itsFixValues;			//taux fixe
	std::vector<double> itsFixedCpnValue;		//proba fixe
	std::vector<double> itsFixedCpnLongSensi;
	std::vector<double> itsFixedCpnShortSensi;

	std::vector<double> itsPayIndexFwd;		//levier index
	std::vector<double> itsPayIndexMult;		//fwd index
	std::vector<double> itsVarCpnValue;		//proba variable	
	std::vector<double> itsVarCpnLongSensi;
	std::vector<double> itsVarCpnShortSensi;

	std::vector<double> itsCpnCoverages;

	std::vector<double> itsLong0;
	std::vector<double> itsShort0;
	std::vector<double> itsPay0;
	bool		  itsIsVariable;

public:
	ARM_BasketSensiCSO(ARM_ZeroCurve* curve,ARM_SpreadOption* so,double w, double p);
	~ARM_BasketSensiCSO(){};

	virtual double	GetCpn(int i);
	virtual double	GetDCpn(int i, const ARM_STRIPPER& stripper);
	virtual void	Init(int i, const ARM_STRIPPER& stripper);

	bool			IsVariableCpn(size_t cpnIndex) const {return itsIsVariable;};
};

class ARM_BasketSensiFunding : public ARM_BasketSensi
{
private:
	std::vector<double> itsFwdRates;
	std::vector<double> itsFundingSpreads;
	std::vector<double> itsCpnCoverages;
public:
	ARM_BasketSensiFunding(ARM_ZeroCurve* curve,ARM_SwapLeg* leg,double w, double p);
	~ARM_BasketSensiFunding(){};

	virtual double	GetCpn(int i);
	virtual double	GetDCpn(int i, const ARM_STRIPPER& stripper);
	virtual void	Init(int i, const ARM_STRIPPER& stripper);
	virtual double ComputeFee(int exerIndex,double startDate);
};

class ARM_BasketSensiFixLeg : public ARM_BasketSensi
{
private:
	std::vector<double> itsCpnCoverages;
	std::vector<double> itsFixRates;

public:
	ARM_BasketSensiFixLeg(ARM_ZeroCurve* curve,ARM_SwapLeg* leg,double w, double p);
	~ARM_BasketSensiFixLeg(){};

	virtual double	GetCpn(int i);
	virtual double	GetDCpn(int i, const ARM_STRIPPER& stripper){return 0.;};
	virtual void	Init(int i, const ARM_STRIPPER& stripper);
	virtual double ComputeFee(int exerIndex,double startDate);

};


CC_END_NAMESPACE()

#endif
//-----------------------------------------------------------------------------
/*---- End of file ----*/
