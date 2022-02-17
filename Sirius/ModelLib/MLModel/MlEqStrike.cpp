//	MlEqStrike.cpp : Implementation of the Dividend Schedule
//
/////////////////////////////////////////////////////////////////////////////

#include "stdafx.h"
#include "MlEqStrike.h"
#include "MlEqObjects.h"


/****************************************************************
**	Class  : MlEqStrike 
**	Routine: MlEqStrike
**	Returns: nothing
**	Action : creates  Strike class
**           
****************************************************************/
MlEqStrike::MlEqStrike(MlEqConstDateHandle startDate)
{
	m_strikeType	= Fixed;
	m_strike		= 0.0;

	m_startDate		= startDate;	
	m_endDate		= -123;
	m_isFloating	= false;

}

/****************************************************************
**	Class  : MlEqStrike 
**	Routine: MlEqStrike
**	Returns: nothing
**	Action : creates  Strike class
**           
****************************************************************/
MlEqStrike::MlEqStrike(double strike,MlEqConstDateHandle startDate)
{
	m_strikeType	= Fixed;
	m_strike		= strike;
	
	m_startDate		= startDate;	
	m_endDate		= -123;
	m_isFloating	= false;

}

/****************************************************************
**	Class  : MlEqStrike 
**	Routine: MlEqStrike
**	Returns: nothing
**	Action : creates  Strike class
**           
****************************************************************/
MlEqStrike::MlEqStrike(const MlEqAsset& asset, double strike,MlEqConstDateHandle startDate)
{
	m_strikeType	= Fixed;
	m_strike		= strike;

	if ( !!startDate ){
		m_startDate		= startDate;
	}
	else{

		MlEqConstDateHandle mdate = asset.GetDateHandle();
		m_startDate    = new MlEqDate(*mdate);
	}
}



/****************************************************************
**	Class  : MlEqStrike 
**	Routine: MlEqStrike
**	Returns: nothing
**	Action : creates  Strike class
**           
****************************************************************/
void MlEqStrike::initialize(double strike)
{
	m_strikeType	= Fixed;
	m_strike		= strike;

	m_startDate		= NULL;
	m_endDate		= -123;
	m_isFloating	= false;
}

/****************************************************************
**	Class  : MlEqStrike 
**	Routine: MlEqStrike
**	Returns: nothing
**	Action : creates  Strike class
**           
****************************************************************/

void MlEqStrike::initialize(double strike,	MlEqConstDateHandle startDate)
{
	m_strikeType	= Fixed;
	m_strike		= strike;
	m_startDate		= startDate;
	m_isFloating	= false;
}

/****************************************************************
**	Class  : MlEqStrike 
**	Routine: MlEqStrike
**	Returns: nothing
**	Action : creates  Strike class
**           
****************************************************************/

void MlEqStrike::deriv(CVector& res)const
{
	res[0] = 1.0;
	res[1] = 0.0;
}

/****************************************************************
**	Class  : MlEqStrike 
**	Routine: MlEqStrike
**	Returns: nothing
**	Action : creates  Strikes
**           
****************************************************************/

void MlEqStrike::reinitializeStrike(const MlEqAsset& asset)
{

	if ( !m_isFloating ){
		return;
	}


	long nStart = m_startDate->GetDate();

	long nToday = asset.GetCurrentDate();
	if ( nStart < nToday ){
		nStart = nToday;
	}
}


/****************************************************************
**	Class  : MlEqForwardBasedStrike 
**	Routine: MlEqForwardBasedStrike
**	Returns: nothing
**	Action : creates forwardbased Strikes
**           
****************************************************************/
MlEqForwardBasedStrike::MlEqForwardBasedStrike()
: MlEqStrike( 1.0 )
{
	m_strikeType	= ForwardBased;
	m_fwdparameter	= 1.0;
}


/****************************************************************
**	Class  : MlEqForwardBasedStrike 
**	Routine: MlEqForwardBasedStrike
**	Returns: nothing
**	Action : creates forwardbased Strikes
**           
****************************************************************/

MlEqForwardBasedStrike::MlEqForwardBasedStrike(double strike,double fwdparam,long enddate,MlEqConstDateHandle startDate,	bool isFloating)
: MlEqStrike( strike )
{
	m_strikeType	= ForwardBased;	
	m_fwdparameter	= fwdparam;
	m_startDate		= startDate;
	m_endDate		= enddate;
	m_isFloating	= isFloating;

}

/****************************************************************
**	Class  : MlEqForwardBasedStrike 
**	Routine: MlEqForwardBasedStrike
**	Returns: nothing
**	Action : creates forwardbased Strikes
**           
****************************************************************/

MlEqForwardBasedStrike::MlEqForwardBasedStrike(const MlEqAsset& asset,double strike,long enddate,MlEqConstDateHandle startDate,bool isFloating)
:MlEqStrike(asset,strike,startDate)
{

	m_strikeType	= ForwardBased;	
	m_fwdparameter	= asset.GetForward(startDate->GetDate(), enddate, true);

	m_endDate		= enddate;
	m_isFloating	= isFloating;

}


/****************************************************************
**	Class  : MlEqNormalizedStrike 
**	Routine: deriv
**	Returns: nothing
**	Action : creates  normalized Strike class
**           
****************************************************************/

void MlEqForwardBasedStrike::deriv(CVector& res)const
{
	double deriv;
	deriv = 1.0/m_fwdparameter;

	res[0] = deriv;
	res[1] = 0.0;

}

/****************************************************************
**	Class  : MlEqForwardBasedStrike 
**	Routine: MlEqForwardBasedStrike
**	Returns: nothing
**	Action : creates forwardbased Strikes
**           
****************************************************************/

void MlEqForwardBasedStrike::reinitializeStrike(const MlEqAsset& asset)
{
	if ( ! m_isFloating ){
		return;
	}

	MlEqStrike::reinitializeStrike(asset);

	long currentDate = asset.GetDateHandle()->GetDate();

	if ( m_endDate < currentDate ){
		throw("this ain't make sense: endDate of normalised strike lower than currentDate of asset");
	}

	long nStart = m_startDate->GetDate();
	long nToday = asset.GetCurrentDate();
	if ( nStart < nToday ){
		nStart = nToday;
	}

	m_fwdparameter = asset.GetForward(nStart,m_endDate,false); 
	MlEqStrike strike(m_fwdparameter);	
}


/****************************************************************
**	Class  : LogBased 
**	Routine: LogBased
**	Returns: nothing
**	Action : creates LogBased Strikes
**           
****************************************************************/

MlEqLogStrike::MlEqLogStrike()
: MlEqStrike( 1.0 )
{
	m_strikeType	= LogBased;	
	m_fwdparameter	= 1.0;
}

/****************************************************************
**	Class  : MlEqForwardBasedStrike 
**	Routine: MlEqForwardBasedStrike
**	Returns: nothing
**	Action : creates forwardbased Strikes
**           
****************************************************************/

MlEqLogStrike::MlEqLogStrike(double strike,double refValue,long enddate,MlEqConstDateHandle startDate,bool isFloating)
: MlEqStrike( strike )
{
	m_strikeType	= LogBased;
	m_fwdparameter	= refValue;
	m_startDate		= startDate;
	m_endDate		= enddate;
	m_isFloating	= isFloating;
}

/****************************************************************
**	Class  : MlEqForwardBasedStrike 
**	Routine: MlEqForwardBasedStrike
**	Returns: nothing
**	Action : creates forwardbased Strikes
**           
****************************************************************/


MlEqLogStrike::MlEqLogStrike(const MlEqAsset& asset,double strike,long enddate,MlEqConstDateHandle startDate,bool isFloating)
:MlEqStrike(asset,strike,startDate)
{

	m_strikeType	= LogBased;	
	m_fwdparameter	= asset.GetForward(startDate->GetDate(), enddate, true);
	m_endDate		= enddate;
	m_isFloating	= isFloating;

}


/****************************************************************
**	Class  : MlEqForwardBasedStrike 
**	Routine: MlEqForwardBasedStrike
**	Returns: nothing
**	Action : creates forwardbased Strikes
**           
****************************************************************/

void MlEqLogStrike::reinitializeStrike(const MlEqAsset& asset)
{
	if ( ! m_isFloating ){
		return;
	}

	MlEqStrike::reinitializeStrike(asset);
	long currentDate = asset.GetDateHandle()->GetDate();
	if ( m_endDate < currentDate ){
		throw("this ain't make sense: endDate of normalised strike lower than currentDate of asset");
	}

	long nStart = m_startDate->GetDate();
	long nToday = asset.GetCurrentDate();
	if ( nStart < nToday ){
		nStart = nToday;
	}

	m_fwdparameter = asset.GetForward(nStart,m_endDate,false); 
}

/****************************************************************
**	Class  : MlEqNormalizedStrike 
**	Routine: deriv
**	Returns: nothing
**	Action : creates  normalized Strike class
**           
****************************************************************/

void MlEqLogStrike::deriv(CVector& res)const
{

	double deriv;
	MlEqStrike fixed;
	convertStrikes(fixed,*this);
	deriv = 1.0/fixed.m_strike;

	res[0] = deriv;
	res[1] = -deriv*deriv;

}

/****************************************************************
**	Class  : MlEqSpotBasedStrike 
**	Routine: MlEqSpotBasedStrike
**	Returns: nothing
**	Action : creates spotbased Strikes
**           
****************************************************************/

MlEqSpotBasedStrike::MlEqSpotBasedStrike()
: MlEqStrike( 1.0 )
{
	m_strikeType	= SpotBased;
	m_spotparameter	= 1.0;
}

/****************************************************************
**	Class  : MlEqSpotBasedStrike
**	Routine: MlEqSpotBasedStrike
**	Returns: nothing
**	Action : creates forwardbased Strikes
**           
****************************************************************/

MlEqSpotBasedStrike::MlEqSpotBasedStrike(double strike,double spotparam,MlEqConstDateHandle startDate,bool isFloating)
: MlEqStrike( strike, startDate)
{
	m_strikeType	= SpotBased;
	m_spotparameter	= spotparam;	
	m_isFloating	= isFloating;

}

/****************************************************************
**	Class  : MlEqSpotBasedStrike
**	Routine: MlEqSpotBasedStrike
**	Returns: nothing
**	Action : creates forwardbased Strikes
**           
****************************************************************/

MlEqSpotBasedStrike::MlEqSpotBasedStrike(const MlEqAsset& asset,double strike,MlEqConstDateHandle startDate,bool isFloating)
: MlEqStrike(asset, strike ,startDate)
{
	m_strikeType	= SpotBased;
	m_spotparameter	= asset.GetSpot(asset.GetDateHandle()->GetDate());
	m_isFloating	= isFloating;
}

/****************************************************************
**	Class  : MlEqForwardBasedStrike 
**	Routine: MlEqForwardBasedStrike
**	Returns: nothing
**	Action : creates forwardbased Strikes
**           
****************************************************************/

void MlEqSpotBasedStrike::reinitializeStrike(const MlEqAsset& asset)
{

	if ( ! m_isFloating ){
		return;
	}

	MlEqStrike::reinitializeStrike(asset);

	long nStart = m_startDate->GetDate();

	long nToday = asset.GetCurrentDate();
	if ( nStart < nToday ){
		nStart = nToday;
	}

	if ( nStart > asset.GetDateHandle()->GetDate()){
		m_spotparameter = asset.GetForward(nStart,false); 
	}
	else{
		m_spotparameter = asset.GetSpot(nStart);
	}

	
}

/****************************************************************
**	Class  : MlEqNormalizedStrike 
**	Routine: deriv
**	Returns: nothing
**	Action : creates  normalized Strike class
**           
****************************************************************/

void MlEqSpotBasedStrike::deriv(CVector& res)const
{
	double deriv = 1.0/m_spotparameter;

	res[0] = deriv;
	res[1] = 0.0;

}


/****************************************************************
**	Class  : MlEqNormalizedStrike 
**	Routine: MlEqNormalizedStrike
**	Returns: nothing
**	Action : creates  normalized Strike class
**           
****************************************************************/

MlEqNormalizedStrike::MlEqNormalizedStrike()
:	MlEqStrike(0.0)
{
	m_strikeType	= Normalised;

	m_vol			= 0;
	m_maturity		= 0;
	m_volscale		= 1.0;
	m_elasticity	= 1.0;
	m_maturityshift = 0.0;
	m_fwdparameter	= 1.0;
}

/****************************************************************
**	Class  : MlEqNormalizedStrike 
**	Routine: deriv
**	Returns: nothing
**	Action : creates  normalized Strike class
**           
****************************************************************/

void MlEqNormalizedStrike::deriv(CVector& res)const
{

	double deriv;
	MlEqStrike fixed;
	convertStrikes(fixed,*this);

	if ( fabs(fixed.m_strike) < 1e-6 ){
		throw("infinite strike derivative encountered");
	}

	deriv = 1.0/fixed.m_strike*m_strike*log(m_fwdparameter/fixed.m_strike);

	res[0] = deriv;
	res[1] = 1.0/fixed.m_strike*deriv;


}

/****************************************************************
**	Class  : MlEqNormalizedStrike 
**	Routine: MlEqNormalizedStrike
**	Returns: nothing
**	Action : reinitializes market data in strikes
**           
****************************************************************/


void MlEqNormalizedStrike::reinitializeStrike(const MlEqAsset& asset)
{
	MlEqStrike::reinitializeStrike(asset);

	long nToday = asset.GetCurrentDate();
	long nStart = m_startDate->GetDate();

	if ( nStart < nToday ){
		nStart = nToday;
	}

	m_fwdparameter = asset.GetForward(nStart,m_endDate,false); 
	MlEqStrike strike(m_fwdparameter);
	m_vol = asset.GetVolatility(strike,nStart,m_endDate); 
	m_maturity	= m_startDate->GetYearFraction(m_endDate);
	
}

/****************************************************************
**	Class  : MlEqNormalizedStrike 
**	Routine: MlEqNormalizedStrike
**	Returns: nothing
**	Action : creates  normalized Strike class
**           
****************************************************************/

void MlEqNormalizedStrike::initialize(const MlEqAsset& asset,double strike, long nMaturity,
									  bool isFloating ,MlEqConstDateHandle startDate,
									  double volscale,double elasticity,double maturityShift)

{

	MlEqStrike::initialize( strike,startDate);

	m_strikeType	= Normalised;
	m_volscale		= volscale;
	m_elasticity	= elasticity;
	m_maturityshift = maturityShift;
	m_isFloating	= isFloating;
	
	m_startDate		= startDate;
	m_fwdparameter	= asset.GetForward(startDate->GetDate(),nMaturity,false); 
	MlEqStrike xstrike(m_fwdparameter);
	m_vol			= asset.GetVolatility(xstrike,startDate->GetDate(),nMaturity); 
	m_maturity		= startDate->GetYearFraction(nMaturity);
	m_endDate		= nMaturity;

}


/****************************************************************
**	Class  : MlEqNormalizedStrike 
**	Routine: MlEqNormalizedStrike
**	Returns: nothing
**	Action : creates  normalized Strike class
**           
****************************************************************/


MlEqNormalizedStrike::MlEqNormalizedStrike(const MlEqAsset& asset,double strike, long nMaturity,
											 bool isFloating,
											 MlEqConstDateHandle startDate,double volscale,double elasticity,
											 double maturityShift)
{

	initialize(asset,strike, nMaturity,
			   isFloating,startDate,
			   volscale,elasticity,maturityShift);

}

/****************************************************************
**	Class  : MlEqNormalizedStrike 
**	Routine: MlEqNormalizedStrike
**	Returns: nothing
**	Action : creates  normalized Strike class
**           
****************************************************************/

MlEqNormalizedStrike::MlEqNormalizedStrike(double strike, long nMaturity,double vol,double fwd,
						 bool isFloating ,
						 MlEqConstDateHandle startDate,double volscale,double elasticity,
						 double maturityShift)
{


	initialize(strike,nMaturity,vol,fwd,
			   isFloating,startDate,
			   volscale,elasticity,maturityShift);

}

/****************************************************************
**	Class  : MlEqNormalizedStrike 
**	Routine: MlEqNormalizedStrike
**	Returns: nothing
**	Action : creates  normalized Strike class
**           
****************************************************************/

void MlEqNormalizedStrike::initialize(double strike, long nMaturity,double vol,double fwd,
						bool isFloating ,MlEqConstDateHandle startDate,
						double volscale,double elasticity,double maturityShift)
{



	MlEqStrike::initialize( strike,startDate);

	m_strikeType	= Normalised;
	m_volscale		= volscale;
	m_elasticity	= elasticity;
	m_maturityshift = maturityShift;
	m_isFloating	= isFloating;
	
	m_startDate	    = startDate;
	m_vol			= vol;
	m_maturity		= startDate->GetYearFraction(nMaturity);
	m_endDate		= nMaturity;
	m_fwdparameter  = fwd;
}


/****************************************************************
**	Class  :  
**	Routine: convertStrikes
**	Returns: nothing
**	Action : converts one strike class into another
**           
****************************************************************/
void MlEqStrike::convertStrikes(MlEqStrike& OutputStrike,const MlEqStrike& InputStrike)
{	

	if(&OutputStrike == &InputStrike)
		return;


	// Given any strike type, first convert to fixed strike, and then convert to output strike type.
	double FixedStrike = 0.0;	

	switch( InputStrike.m_strikeType )
	{
	case Fixed:
		{
			FixedStrike		= InputStrike.m_strike;
			break;
		}
	case ForwardBased:
		{
			MlEqForwardBasedStrike* IPtr = (MlEqForwardBasedStrike*) (&InputStrike);
			FixedStrike		= IPtr->m_fwdparameter * InputStrike.m_strike;
			break;
		}
	case SpotBased:
		{
			MlEqSpotBasedStrike* IPtr = (MlEqSpotBasedStrike*) (&InputStrike);
			FixedStrike		= IPtr->m_spotparameter * InputStrike.m_strike;
			break;
		}
	case Normalised:
		{
			MlEqNormalizedStrike* IPtr = (MlEqNormalizedStrike*) (&InputStrike);

			FixedStrike		= IPtr->m_fwdparameter*exp(InputStrike.m_strike*
								(IPtr->m_volscale*pow(IPtr->m_vol/IPtr->m_volscale,IPtr->m_elasticity)*
								pow(IPtr->m_maturity*IPtr->m_maturity+IPtr->m_maturityshift,0.25)));
			break;
		}
	case LogBased:
		{
			MlEqLogStrike* IPtr = (MlEqLogStrike*) (&InputStrike);
			FixedStrike		= IPtr->m_fwdparameter * exp(InputStrike.m_strike);
			break;
		}
	default :
		throw(  "unknown input strikeType encountered" );

	}

	// Now, convert to new strike
	// reuse FixedStrike variable
	switch ( OutputStrike.m_strikeType )
	{
	case Fixed:
		break;

	case ForwardBased:
		{
			MlEqForwardBasedStrike* OutStrikePtr = (MlEqForwardBasedStrike*) (&OutputStrike);

			FixedStrike /= OutStrikePtr->m_fwdparameter;
			break;
		}

	case SpotBased:
		{
			MlEqSpotBasedStrike* OutStrikePtr = (MlEqSpotBasedStrike*) (&OutputStrike);

			FixedStrike /= OutStrikePtr->m_spotparameter;
			break;
		}
	case Normalised:
		{
			MlEqNormalizedStrike* OPtr = (MlEqNormalizedStrike*) (&OutputStrike);

			FixedStrike = log(FixedStrike/OPtr->m_fwdparameter)/
				(OPtr->m_volscale*pow(OPtr->m_vol/OPtr->m_volscale,OPtr->m_elasticity)*
				 pow(OPtr->m_maturity*OPtr->m_maturity+OPtr->m_maturityshift,0.25));

			break;
		}
	case LogBased:
		{

			MlEqLogStrike* OPtr = (MlEqLogStrike*) (&OutputStrike);			
			
			FixedStrike	= log(FixedStrike/OPtr->m_fwdparameter);
			break;
		}
	default:
			throw(  "unknown output strikeType encountered" );
	}

	OutputStrike.m_strike = FixedStrike;

}

