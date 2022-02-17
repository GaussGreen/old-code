//	MlEqStrike.h : strike object class
//
/////////////////////////////////////////////////////////////////////////////

#ifndef _MLEQSTRIKE_H_
#define _MLEQSTRIKE_H_

#include "smart.h"


/////////////////////////////////////////////////////////////////////////////
// We expose a vector of strikes to the end-user.
// Hence the following structure.
//
class MlEqStrikes : public RCObject
{
public:
	std::vector<MlEqStrikeHandle>		m_hStrikes;
};


/////////////////////////////////////////////////////////////////////////////
// Generic strike class.
//
class MlEqStrike : public RCObject
{	

protected:

	MlEqConstDateHandle					m_startDate;
	long								m_endDate;						// may not be set
	bool								m_isFloating;

public:

	bool GetIsFloating(void) const		{return m_isFloating;}
	void SetEndDate(long endDate){m_endDate = endDate;};
	long GetEndDate()const {return m_endDate;};
	void SetStartDateHandle(MlEqConstDateHandle startDate){m_startDate = startDate;};
	MlEqConstDateHandle GetStartDateHandle()const {return m_startDate;};

	static void convertStrikes(MlEqStrike& OutputStrike,const MlEqStrike& Inputstrike);

public:
	
	MlEqStrike(MlEqConstDateHandle startDate=NULL);
	MlEqStrike(double strike,MlEqConstDateHandle startDate=NULL);
	MlEqStrike(const MlEqAsset& asset,double strike,MlEqConstDateHandle startDate);

	void initialize(double strike);
	void initialize(double strike,MlEqConstDateHandle startDate);

	virtual MlEqStrikeHandle copy() const{ return new MlEqStrike( *this ); }
	virtual void reinitializeStrike(const MlEqAsset& asset);
	
	virtual void deriv(CVector& res)const;
	
	double					m_strike;
	StrikesTypeEnum			m_strikeType;	
};		
	
	
	
	
/******************************************************************
**
**	Normalized Strike Class
**
******************************************************************/
	
class  MlEqNormalizedStrike : public MlEqStrike
{


public:
	
	
	MlEqNormalizedStrike();


	MlEqNormalizedStrike(const MlEqAsset& asset,double strike, long nMaturity,
						 bool isFloating = true,
						 MlEqConstDateHandle startDate= NULL,double volscale=1.0,double elasticity=1.0,
						 double maturityShift=0.0);


	void initialize(const MlEqAsset& asset,double strike, long nMaturity,
					bool isFloating ,MlEqConstDateHandle startDate,
					double volscale=1.0,double elasticity=1.0,double maturityShift=0.0);

	MlEqNormalizedStrike(double strike, long nMaturity,double vol,double fwd,
						 bool isFloating = true,
						 MlEqConstDateHandle startDate= NULL,double volscale=1.0,double elasticity=1.0,
						 double maturityShift=0.0);


	void initialize(double strike, long nMaturity,double vol,double fwd,
					bool isFloating = true,MlEqConstDateHandle startDate= NULL,
					double volscale=1.0,double elasticity=1.0,double maturityShift=0.0);

	
	virtual MlEqStrikeHandle copy() const{ return new MlEqNormalizedStrike( *this ); }
	virtual void reinitializeStrike(const MlEqAsset& asset);


	virtual void deriv(CVector& res)const;
	
	// MEMBER DATA	
	double m_volscale;
	double m_maturityshift;
	double m_elasticity;
	double m_fwdparameter;
	double m_vol;
	double m_maturity;	
	
};		

		
/******************************************************************
**
**	ForwardBased Strike Class
**
******************************************************************/
	
			
class  MlEqForwardBasedStrike : public MlEqStrike
{	

	
	public:

	MlEqForwardBasedStrike();
	MlEqForwardBasedStrike(const MlEqAsset& asset,double strike,long enddate,MlEqConstDateHandle startDate=NULL,bool isFloating=true);
	MlEqForwardBasedStrike(double strike,double fwdparam,long enddate,MlEqConstDateHandle startDate=NULL,	bool isFloating=true);

	
	virtual MlEqStrikeHandle copy() const{ return new MlEqForwardBasedStrike( *this ); }

	double m_fwdparameter;  // ie, forward price, but do I need to know about the term it is based?
	virtual void reinitializeStrike(const MlEqAsset& asset);
	
	virtual void deriv(CVector& res)const;

};
	

/******************************************************************
**
**	LogBased Strike Class
**
******************************************************************/
	
			
class  MlEqLogStrike : public MlEqStrike
{	
	
	public:

	MlEqLogStrike();
	MlEqLogStrike(double strike,double refVal,long enddate,MlEqConstDateHandle startDate=NULL,bool isFloating=true);
	
	virtual MlEqStrikeHandle copy() const{ return new MlEqLogStrike( *this ); }
	MlEqLogStrike(const MlEqAsset& asset,double strike,long enddate,MlEqConstDateHandle startDate,bool isFloating=true);


	virtual void deriv(CVector& res)const;

	virtual void reinitializeStrike(const MlEqAsset& asset);

	double m_fwdparameter;// must be forward value  

};

/******************************************************************
**
**	SpotBased Strike Class
**
******************************************************************/

class  MlEqSpotBasedStrike : public MlEqStrike
{	
	
	public:

	MlEqSpotBasedStrike();
	MlEqSpotBasedStrike(double strike,double spotParam,MlEqConstDateHandle startDate=NULL,bool isFloating=true); // spotParam is Spot Price

	MlEqSpotBasedStrike(const MlEqAsset& asset,double strike,MlEqConstDateHandle startDate,bool isFloating=true);

	virtual MlEqStrikeHandle copy() const{ return new MlEqSpotBasedStrike( *this ); }
	void reinitializeStrike(const MlEqAsset& asset);


	virtual void deriv(CVector& res)const;

	double m_spotparameter;  
};	

	
	



#endif
