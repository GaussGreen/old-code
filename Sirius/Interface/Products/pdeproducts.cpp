
#include "stdafx.h"
#include "pdeproducts.h"
#include "MlEqObjects.h"




/****************************************************************
**	Class  : knockOutBarrierHelper 
**	Routine: strategy_function
**	Returns: 
**	Action : 
**           
****************************************************************/

plainVanillaOptionHelper::plainVanillaOptionHelper( CVector& CallPut,GVector<MlEqStrikeHandle>& strike,double spot)
{
	initialize( CallPut,strike,spot);
} 


void plainVanillaOptionHelper::initialize( CVector& CallPut,GVector<MlEqStrikeHandle>& strike,double spot)
{

		m_CallPut		= CallPut;
		m_spot			= spot;

		m_strikes.resize(strike.getsize());
		for ( int i = 0;  i < strike.getsize(); i++ )
		{
			MlEqStrike stk;
			MlEqStrike::convertStrikes(stk,*strike[i]);

			m_strikes[i] = stk.m_strike;
		}
} 


/****************************************************************
**	Class  : plainVanillaOptionHelper 
**	Routine: initial_condition
**	Returns: 
**	Action : 
**           
****************************************************************/

void plainVanillaOptionHelper::initial_condition(void* ptr,int nd,int nx,CMatrix& u,MlEqPdeDriver* pde )
{

	if ( nd != m_strikes.getsize()){
		throw("inconsistent pde set up");
	}

	int ix;
	for (ix = 0; ix < nx;ix++ )
	{
		double x = pde->pde_x(ix);
		for ( int id = 0 ; id < nd; id++ ){
			u[id][ix] = MlEqMaths::Max(m_CallPut[id]*(m_spot*exp(x)-m_strikes[id]),0.0);
		}
	}

}	

/****************************************************************
**	Class  : plainVanillaOptionHelper 
**	Routine: boundary_condition
**	Returns: 
**	Action : 
**           
****************************************************************/

void plainVanillaOptionHelper::boundary_condition(int id,int it,double t, void* ptr,pde_boundary_condition& lower,pde_boundary_condition& upper,MlEqPdeDriver* pde)
{

	lower.coord = PdeEnums::PDE_COORDINATES_EXP;
	lower.scale = 1.0;
	lower.shift = 0.0;

	lower.a = 1.0;
	lower.b = 0.0;
	lower.c = 0.0;
    lower.d = 0.0;


	upper.coord = PdeEnums::PDE_COORDINATES_EXP;
	upper.scale = 1.0;
	upper.shift = 0.0;

	upper.a = 1.0;
	upper.b = 0.0;
	upper.c = 0.0;
	upper.d = 0.0;

	
}



/****************************************************************
**	Class  : knockOutBarrierHelper 
**	Routine: strategy_function
**	Returns: 
**	Action : 
**           
****************************************************************/
americanOptionHelper::americanOptionHelper(double spot,double CallPut,double strike)
{
	initialize(spot,CallPut,strike);
}

void americanOptionHelper::initialize(double spot,double CallPut,double strike)

{		m_spot			= spot;
		m_CallPut		= CallPut;
		m_strike		= strike;
} 

/****************************************************************
**	Class  : americanOptionHelper 
**	Routine: initial_condition
**	Returns: 
**	Action : 
**           
****************************************************************/

void americanOptionHelper::initial_condition(void* ptr,int nd,int nx,CMatrix& u,MlEqPdeDriver* pde )
{

	int ix;
	for (ix = 0; ix < nx;ix++ )
	{
		double x = pde->pde_x(ix);
		u[0][ix] = MlEqMaths::Max(m_CallPut*(m_spot*exp(x)-m_strike),0.0);	
	}
}	


/****************************************************************
**	Class  : knockOutBarrierHelper 
**	Routine: boundary_condition
**	Returns: 
**	Action : 
**           
****************************************************************/

void americanOptionHelper::boundary_condition(int id,int it,double t, void* ptr,pde_boundary_condition& lower,pde_boundary_condition& upper,MlEqPdeDriver* pde)
{

	lower.coord = PdeEnums::PDE_COORDINATES_EXP;
	lower.scale = 1.0;
	lower.shift = 0.0;

	lower.a = 1.0;
	lower.b = 0.0;
	lower.c = 0.0;
    lower.d = 0.0;


	upper.coord = PdeEnums::PDE_COORDINATES_EXP;
	upper.scale = 1.0;
	upper.shift = 0.0;

	upper.a = 1.0;
	upper.b = 0.0;
	upper.c = 0.0;
	upper.d = 0.0;

	lower.b = -1.;	// in this coordinates gamma( S ) = 0 <=> gamma( x ) = delta( x )
	upper.b = -1.;

	
}

/****************************************************************
**	Class  : knockOutBarrierHelper 
**	Routine: strategy_function
**	Returns: 
**	Action : 
**           
****************************************************************/

PdeEnums::PdeMode americanOptionHelper::strategy_function(int it, double t,  void* ptr, 

                   int nd, int nx, CMatrix& u_it,MlEqPdeDriver* pde)

{

	int ix;
	for (ix = 0; ix < nx; ix++ )
	{
		double x = pde->pde_x(ix);
		double intrinsic = MlEqMaths::Max(m_CallPut*(m_spot*exp(x)-m_strike),0.0);
		u_it[0][ix] = MlEqMaths::Max( intrinsic, u_it[0][ix] );	
	}

	return PdeEnums::PDE_NORMAL_MODE;
} 




/****************************************************************
**	Class  : knockOutBarrierHelper 
**	Routine: strategy_function
**	Returns: 
**	Action : 
**           
****************************************************************/
knockOutBarrierHelper::knockOutBarrierHelper(double UpperBarrier,double LowerBarrier,
											double UpperRebate,double LowerRebate,
											double spot,double CallPut,double strike,
											bool payEnd, const CVector& df)
{
	initialize(UpperBarrier,LowerBarrier,UpperRebate,LowerRebate,spot,CallPut,strike,payEnd, df);
}


void knockOutBarrierHelper::initialize(double UpperBarrier,double LowerBarrier,
									   double UpperRebate,double LowerRebate,
									   double spot,double CallPut,double strike,
									   bool payEnd, const CVector& df)

{
		m_UpperBarrier	= UpperBarrier;
		m_LowerBarrier	= LowerBarrier;
		m_UpperRebate	= UpperRebate;
		m_LowerRebate	= LowerRebate;
		m_spot			= spot;
		m_CallPut		= CallPut;
		m_strike		= strike;
		m_logBarrUp		= log(m_UpperBarrier/m_spot);
		m_logBarrDo		= log(m_LowerBarrier/m_spot);

		m_payEnd		= payEnd;
		m_df			= df;
} 


/****************************************************************
**	Class  : knockOutBarrierHelper 
**	Routine: initial_condition
**	Returns: 
**	Action : 
**           
****************************************************************/

void knockOutBarrierHelper::initial_condition(void* ptr,int nd,int nx,CMatrix& u,MlEqPdeDriver* pde )
{

	int ix;
	for (ix = 0; ix < nx;ix++ )
	{
		double x = pde->pde_x(ix);

		if ( x >= m_logBarrUp - 1e-9 )
			u[0][ix] = m_UpperRebate;
		else if ( x <= m_logBarrDo + 1e-9 )
			u[0][ix] = m_LowerRebate;
		else
			u[0][ix] = MlEqMaths::Max(m_CallPut*(m_spot*exp(x)-m_strike),0.0);	
	}

	m_lowSpot = u[0][0];
	m_highSpot = u[0][ix-1];

}	


/****************************************************************
**	Class  : knockOutBarrierHelper 
**	Routine: boundary_condition
**	Returns: 
**	Action : 
**           
****************************************************************/

void knockOutBarrierHelper::boundary_condition(int id,int it,double t, void* ptr,pde_boundary_condition& lower,pde_boundary_condition& upper,MlEqPdeDriver* pde)
{

	lower.coord = PdeEnums::PDE_COORDINATES_EXP;
	lower.scale = 1.0;
	lower.shift = 0.0;

	lower.a = 0.0;
	lower.b = 0.0;
	lower.c = 1.0;
    

	upper.coord = PdeEnums::PDE_COORDINATES_EXP;
	upper.scale = 1.0;
	upper.shift = 0.0;

	upper.a = 0.0;
	upper.b = 0.0;
	upper.c = 1.0;

	if( m_payEnd )
	{
		double df = m_df[it];
		lower.d = m_LowerRebate * df;
		upper.d = m_UpperRebate * df;
	}
	else
	{
		lower.d = m_LowerRebate;
		upper.d = m_UpperRebate;
	}

	
}

/****************************************************************
**	Class  : callableNote00Helper 
**	Routine: initialize
**	Returns: 
**	Action : 
**           
****************************************************************/


void callableNote00Helper::initialize(	double spot, 
										double refSpot,
										MlEqConstDateHandle hDate,	// this is today...
										const CVector& coupons, 
										const GVector<long>& couponDates,
										const GVector<long>& pdeDates,
										double CallPut,double strike, double gearing, double fixed_coupon)
{
		m_spot			= spot;
		m_factor		= refSpot / spot;
		m_CallPut		= CallPut;
		m_strike		= strike;
		m_gearing		= gearing;
		m_fixed_coupon	= fixed_coupon;

// now build coupon schedule

		int nt = pdeDates.getsize();
		int nct = couponDates.getsize();
		if( nct > nt )	throw" too few steps in the pde";
		m_coupons.resize(nt);
		m_isCouponDate.resize(nt);

		int it = 0;
		long nextCDate, lastCDate = couponDates[0];

		for(int cdate=1; cdate<nct; ++cdate)
		{
			nextCDate = couponDates[cdate];			
			double tprec = hDate->GetYearFraction(lastCDate);
			double yFrac = hDate->GetYearFraction(nextCDate) - tprec;

			while( it < nt && pdeDates[it] <= nextCDate )
			{
				double x = (hDate->GetYearFraction(pdeDates[it]) - tprec) / yFrac;
				m_coupons[it] = coupons[cdate-1] * x;
				m_isCouponDate[it] = false;
				it++;
			}
			if( it > 0 ){ // otherwise problem for past fixings...
				m_isCouponDate[it-1] = true;
			}
			lastCDate = nextCDate;
		}

		it = nt-1;
		lastCDate = couponDates[nct-1];

		while( it >= 0 && pdeDates[it] > lastCDate ) // not necessary, but safer
		{
			m_isCouponDate[it] = false;
			m_coupons[it] = 0.0;
			it--;
		}

}

/****************************************************************
**	Class  : callableNote00Helper 
**	Routine: initial_condition
**	Returns: 
**	Action : 
**           
****************************************************************/


void callableNote00Helper::initial_condition(void* ptr,int nd,int nx,CMatrix& u,MlEqPdeDriver* pde )
{
	int nt = m_isCouponDate.getsize() - 1;
	double dcoupon	= (m_isCouponDate[nt])?m_coupons[nt]:0.0;

	int ix;
	for (ix = 0; ix < nx;ix++ )
	{
		double x = pde->pde_x(ix);
		u[0][ix] = dcoupon + m_fixed_coupon + m_gearing * MlEqMaths::Max(m_CallPut*(m_factor*exp(x)-m_strike),0.0);	
//		u[1][ix] = 1.;
	}
}	


/****************************************************************
**	Class  : callableNote00Helper 
**	Routine: boundary_condition
**	Returns: 
**	Action : 
**           
****************************************************************/


void callableNote00Helper::boundary_condition(int id,int it,double t, void* ptr,pde_boundary_condition& lower,pde_boundary_condition& upper,MlEqPdeDriver* pde)
{

	lower.coord = PdeEnums::PDE_COORDINATES_EXP;
	lower.scale = 1.0;
	lower.shift = 0.0;

	lower.a = 1.0;
	lower.c = 0.0;
    lower.d = 0.0;


	upper.coord = PdeEnums::PDE_COORDINATES_EXP;
	upper.scale = 1.0;
	upper.shift = 0.0;

	upper.a = 1.0;
	upper.c = 0.0;
	upper.d = 0.0;

	lower.b = -1.;	// in this coordinates gamma( S ) = 0 <=> gamma( x ) = delta( x )
	upper.b = -1.;

	
}


/****************************************************************
**	Class  : callableNote00Helper 
**	Routine: strategy_function
**	Returns: 
**	Action : 
**           
****************************************************************/


PdeEnums::PdeMode callableNote00Helper::strategy_function(int it, double t,  void* ptr, 

                   int nd, int nx, CMatrix& u_it,MlEqPdeDriver* pde)

{
	double coupon_t	= m_coupons[it];
	double dcoupon	= (m_isCouponDate[it])?m_coupons[it]:0.0;

	int ix;
	for (ix = 0; ix < nx; ix++ )
	{
		double x = pde->pde_x(ix);

		double intrinsic = MlEqMaths::Max(m_CallPut*(m_factor*exp(x)-m_strike),0.0);
		intrinsic = m_fixed_coupon + m_gearing * intrinsic;
		intrinsic += coupon_t ;

		u_it[0][ix] = MlEqMaths::Max( intrinsic, u_it[0][ix] + dcoupon ) ;	

	}

	return PdeEnums::PDE_NORMAL_MODE;
} 





/****************************************************************
**	Class  : callableNote01Helper 
**	Routine: initialize
**	Returns: 
**	Action : 
**           
****************************************************************/


void callableNote01Helper::initialize(	double spot, 
										double refSpot,
										MlEqConstDateHandle hDate,	// this is today...
										const CVector& coupons, 
										const GVector<long>& couponDates,
										const CVector& rebates, 
										const GVector<long>& exerciseDates,
										const GVector<long>& pdeDates,
										double CallPut,double strike, double gearing, double fixed_coupon)
{
		m_spot			= spot;
		m_factor		= refSpot / spot;
		m_CallPut		= CallPut;
		m_strike		= strike;
		m_gearing		= gearing;
		m_fixed_coupon	= fixed_coupon;

// now build coupon schedule

		int nt = pdeDates.getsize();
		int nct = couponDates.getsize();
		if( nct > nt )	throw" too few steps in the pde";
		m_coupons.resize(nt);
		m_isCouponDate.resize(nt);
	

		int it = 0;
		long nextCDate, lastCDate = couponDates[0];

		for(int cdate=1; cdate<nct; ++cdate)
		{
			nextCDate = couponDates[cdate];			
			double tprec = hDate->GetYearFraction(lastCDate);
			double yFrac = hDate->GetYearFraction(nextCDate) - tprec;

			while( it < nt && pdeDates[it] <= nextCDate )
			{
				double x = (hDate->GetYearFraction(pdeDates[it]) - tprec) / yFrac;
				m_coupons[it] = coupons[cdate-1] * x;
				m_isCouponDate[it] = false;
				it++;
			}
			if( it > 0 ){ // otherwise problem for past fixings...
				m_isCouponDate[it-1] = true;
			}
			lastCDate = nextCDate;
		}

		it = nt-1;
		lastCDate = couponDates[nct-1];

		while( it >= 0 && pdeDates[it] > lastCDate ) // not necessary, but safer
		{
			m_isCouponDate[it] = false;
			m_coupons[it] = 0.0;
			it--;
		}

		nct = exerciseDates.getsize();
		if( nct > nt )	throw" too few steps in the pde";
		m_isCallableDate.resize(nt);
		m_rebates.resize(nt);

		it = 0;
		lastCDate = pdeDates[0] ;
		
		for(int cdate=0; cdate<nct; ++cdate)
		{
			nextCDate = exerciseDates[cdate];			
	
			while( it < nt && pdeDates[it] <= nextCDate ){
				m_isCallableDate[it++] = false;
			}
			if( it > 0 ){ // otherwise problem for past fixings...
				m_isCallableDate[it-1] = true;
				m_rebates[it-1] = rebates[cdate];
			}
			lastCDate = nextCDate;
		}

		it = nt-1;
		lastCDate = exerciseDates[nct-1];

		while( it >= 0 && pdeDates[it] > lastCDate ) {	// not necessary, but safer
			m_isCallableDate[it--] = false;
		}
}



/**************************************************************
**
**	Class  : callableNote00Helper 
**	Routine: initial_condition
**	Returns: 
**	Action : 
**           
***************************************************************/


void callableNote01Helper::initial_condition(void* ptr,int nd,int nx,CMatrix& u,MlEqPdeDriver* pde )
{
	int nt = m_isCouponDate.getsize() - 1;
	double dcoupon	= (m_isCouponDate[nt])?m_coupons[nt]:0.0;

	bool exit = m_isCallableDate[nt];
	double rebate = m_rebates[nt];

	int ix;
	for (ix = 0; ix < nx;++ix )
	{
		double x = pde->pde_x(ix);
		double payout = dcoupon + m_fixed_coupon + m_gearing * std::max(m_CallPut*(m_factor*exp(x)-m_strike),0.0);	
		if( exit ){
			payout = std::max( payout, rebate );
		}		
		u[0][ix] = payout ;
	}
}	


/****************************************************************
**	Class  : callableNote00Helper 
**	Routine: boundary_condition
**	Returns: 
**	Action : 
**           
****************************************************************/


void callableNote01Helper::boundary_condition(int id,int it,double t, void* ptr,pde_boundary_condition& lower,pde_boundary_condition& upper,MlEqPdeDriver* pde)
{

	lower.coord = PdeEnums::PDE_COORDINATES_EXP;
	lower.scale = 1.0;
	lower.shift = 0.0;

	lower.a = 1.0;
	lower.c = 0.0;
    lower.d = 0.0;


	upper.coord = PdeEnums::PDE_COORDINATES_EXP;
	upper.scale = 1.0;
	upper.shift = 0.0;

	upper.a = 1.0;
	upper.c = 0.0;
	upper.d = 0.0;

	lower.b = -1.;	// in this coordinates gamma( S ) = 0 <=> gamma( x ) = delta( x )
	upper.b = -1.;

	
}


/****************************************************************
**	Class  : callableNote00Helper 
**	Routine: strategy_function
**	Returns: 
**	Action : 
**           
****************************************************************/


PdeEnums::PdeMode callableNote01Helper::strategy_function(int it, double t,  void* ptr, 

                   int nd, int nx, CMatrix& u_it,MlEqPdeDriver* pde)

{
	if( !m_isCallableDate[it] && !m_isCouponDate[it] ){
		return PdeEnums::PDE_NORMAL_MODE;
	}



	double dcoupon	= (m_isCouponDate[it])?m_coupons[it]:0.0;
	double rebate = m_rebates[it] + m_coupons[it];

	for (int ix = 0; ix < nx; ix++ ){
		u_it[0][ix] = std::max( rebate, u_it[0][ix] + dcoupon ) ;	
	}

	return PdeEnums::PDE_NORMAL_MODE;
} 






/****************************************************************
**	Class  : knockOutBarrierHelper01 
**	Routine: knockOutBarrierHelper01
**	Returns: 
**	Action : 
**           
****************************************************************/



knockOutBarrierHelper01::knockOutBarrierHelper01(double UpperBarrier,double LowerBarrier,
											double UpperRebate,double LowerRebate, bool isUp, bool isDown,
											const CVector& dfDn, const CVector& dfUp,
											MlEqInterpolatorHandle hFinalValue)
{
	initialize(UpperBarrier,LowerBarrier,UpperRebate,LowerRebate,isUp, isDown, dfDn, dfUp, hFinalValue);
}


/****************************************************************
**	Class  : knockOutBarrierHelper01 
**	Routine: knockOutBarrierHelper01
**	Returns: 
**	Action : 
**           
****************************************************************/


knockOutBarrierHelper01::knockOutBarrierHelper01(BarrierFeatures* bfeat, 
												 const CVector& dfDn, const CVector& dfUp,
												MlEqInterpolatorHandle hFinalValue)
{
		m_UpperBarrier	= bfeat->barrUp;
		m_LowerBarrier	= bfeat->barrDown;
		m_UpperRebate	= bfeat->rebateUp;
		m_LowerRebate	= bfeat->rebateDown;
		m_isUpBarrier	= bfeat->isBarrierUp;
		m_isDownBarrier	= bfeat->isBarrierDown;	
	
		m_dfDn			= dfDn;
		m_dfUp			= dfUp;
		m_hFinalValue	= hFinalValue->copy();
}


/****************************************************************
**	Class  : knockOutBarrierHelper01 
**	Routine: knockOutBarrierHelper01
**	Returns: 
**	Action : 
**           
****************************************************************/


void knockOutBarrierHelper01::initialize(double UpperBarrier,double LowerBarrier,
									   double UpperRebate,double LowerRebate, bool isUp, bool isDown,
									   const CVector& dfDn, const CVector& dfUp,
									   MlEqInterpolatorHandle hFinalValue)

{
		m_UpperBarrier	= UpperBarrier;
		m_LowerBarrier	= LowerBarrier;
		m_UpperRebate	= UpperRebate;
		m_LowerRebate	= LowerRebate;
		m_isUpBarrier	= isUp;
		m_isDownBarrier	= isDown;

	
		m_dfDn			= dfDn;
		m_dfUp			= dfUp;

		m_hFinalValue	= hFinalValue;
} 



/****************************************************************
**	Class  : knockOutBarrierHelper 
**	Routine: initial_condition
**	Returns: 
**	Action : 
**           
****************************************************************/


void knockOutBarrierHelper01::initial_condition(void* ptr,int nd,int nx,CMatrix& u,MlEqPdeDriver* pde )
{
	for(int id=0; id<nd; ++id)
	{		
		for (int ix = 0; ix < nx; ++ix)
		{
			double x = pde->pde_x(ix);
			u[id][ix] = m_hFinalValue->getValue( x, id );	
		}
	}

}

/****************************************************************
**	Class  : knockOutBarrierHelper 
**	Routine: boundary_condition
**	Returns: 
**	Action : 
**           
****************************************************************/

void knockOutBarrierHelper01::boundary_condition(int id,int it,double t, void* ptr,pde_boundary_condition& lower,pde_boundary_condition& upper,MlEqPdeDriver* pde)
{

	lower.coord = PdeEnums::PDE_COORDINATES_EXP;
	lower.scale = 1.0;
	lower.shift = 0.0;   

	upper.coord = PdeEnums::PDE_COORDINATES_EXP;
	upper.scale = 1.0;
	upper.shift = 0.0;

	lower.a = 1.0;
	lower.c = 0.0;
	lower.d = 0.0;

	upper.a = 1.0;
	upper.c = 0.0;
	upper.d = 0.0;

	lower.b = -1.;		// default boundary gamma = 0
	upper.b = -1.;


	if( m_isUpBarrier ) 
	{
		upper.a = 0.0;
		upper.b = 0.0;
		upper.c = 1.0;
		upper.d = m_UpperRebate * m_dfUp[it];;
	}
	if( m_isDownBarrier ) 
	{			
		lower.a = 0.0;
		lower.b = 0.0;
		lower.c = 1.0;
		lower.d = m_LowerRebate * m_dfDn[it];;

	}	
	
}

/****************************************************************
**	Class  : knockOutBarrierHelper01 
**	Routine: finalize
**	Returns: 
**	Action : 
**           
****************************************************************/

void knockOutBarrierHelper01::finalize( int id0, MlEqPdeDriver* pde )
{
	CMatrix u = (pde->pde_full_solution())[0]; // final solution
	int nx = u.cols();
	int nd = u.rows();


	for(int id=0; id<nd; ++id)
	{
		std::vector<double> xValues, yValues;
		double x = pde->pde_x(0);


		if( m_isDownBarrier ) 
		{
			x *= (x>0.)? 0.5:1.5;
			if( fabs(x) < 1e-6 ){
				x = -1e-2;
			}

			xValues.push_back( x );
			yValues.push_back( m_LowerRebate * m_dfDn[0] );
		}

		for (int ix = 0; ix < nx;ix++ )
		{
			x = pde->pde_x(ix);
			xValues.push_back( x );
			yValues.push_back( u[id][ix] );	
		}
		
		if( m_isUpBarrier ) 
		{
			x *= (x<0.)? 0.5:1.5;
			if( fabs(x) < 1e-6 ){
				x = 1e-2;
			}
			xValues.push_back( x );
			yValues.push_back( m_UpperRebate * m_dfUp[0] );
		}


		int new_nx = xValues.size();
		CVector xv(new_nx), yv(new_nx);

		for(int i=0 ;i<new_nx; ++i)
		{
			xv[i] = xValues[i];
			yv[i] = yValues[i];
		}

		m_hFinalValue->reinitialize( xv, yv , id);
	}

	
}



/****************************************************************
**	Class  : BarrierPricer 
**	Routine: BarrierPricer
**	Returns: 
**	Action : 
**           
****************************************************************/

BarrierPricer::BarrierPricer(MlEqAssetHandle			hUnderlying,
							 GVector<BarrierFeatures>	barrierFeatures,
							 long	maturityDate)
{
	initialize(hUnderlying,	barrierFeatures, maturityDate);
}



/****************************************************************
**	Class  : BarrierPricer 
**	Routine: initialize
**	Returns: 
**	Action : 
**           
****************************************************************/


void BarrierPricer::initialize(	 MlEqAssetHandle			hUnderlying,
								 GVector<BarrierFeatures>	barrierFeatures,
								 long						maturityDate)
{
	m_hUnderlying		=	hUnderlying;
	m_maturityDate		=	maturityDate;
	m_storeResult		=	false;	// for debug purposes

	std::vector<long> dates;
	std::set<long>	sdates;
	int nfeat = barrierFeatures.getsize();

	long nToday = m_hUnderlying->GetDateHandle()->GetDate();
	sdates.insert(nToday);

	for(int i=0; i<nfeat; ++i)
	{
		sdates.insert(barrierFeatures[i].startDate);
		sdates.insert(barrierFeatures[i].endDate);
	}
	sdates.insert(m_maturityDate);

	std::set<long> ::iterator  it = sdates.begin();
	int nslice = 0;
	for(it = sdates.begin(); it!= sdates.end(); it++)
	{
		dates.push_back(*it);
		nslice++;
	}

	m_barrierFeatures.resize(nslice-1);

	m_priceKnockIn = false;

	for(int slice=0; slice<nslice-1; ++slice)
	{
		BarrierFeatures tmp;
		tmp.startDate = dates[slice];
		tmp.endDate = dates[slice+1];	

		int fcount = 0;
	
		for(int i=0; i<nfeat; i++)
		{
			if(barrierFeatures[i].startDate <= dates[slice] && dates[slice+1] <= barrierFeatures[i].endDate )
			{
				fcount++;
				bool isup = barrierFeatures[i].isBarrierUp;
				bool isdn = barrierFeatures[i].isBarrierDown;

				if( isdn )
				{
					tmp.isBarrierDown	 = true;
					tmp.barrDown		 = barrierFeatures[i].barrDown ;
					tmp.rebateDown		 = barrierFeatures[i].rebateDown ;
					tmp.isKnockOutDown	 = barrierFeatures[i].isKnockOutDown ;
					tmp.payDateDn		 = barrierFeatures[i].payDateDn;
				}
				if( isup )
				{
					tmp.isBarrierUp	 = true;
					tmp.barrUp		 = barrierFeatures[i].barrUp ;
					tmp.rebateUp	 = barrierFeatures[i].rebateUp ;
					tmp.isKnockOutUp = barrierFeatures[i].isKnockOutUp ;
					tmp.payDateUp	 = barrierFeatures[i].payDateUp;
				}
			}
		}
		if( fcount > 2 )
			throw"too many features overlapping";

		m_barrierFeatures[slice] = tmp;

		if( !tmp.isKnockOutDown || !tmp.isKnockOutUp ){
			m_priceKnockIn = true;
		}
	}

	m_knockInRebate = 0.;

	if( m_priceKnockIn )
	{			
		for(int slice=nslice-2; slice>=0; --slice)
		{
			BarrierFeatures tmp = m_barrierFeatures[slice];
			MlEqZeroCurveHandle hZC = m_hUnderlying->GetPayZeroCurve(true);

			if( tmp.isBarrierDown && !tmp.isKnockOutDown )
			{
				long payDn = tmp.payDateDn;
				if( payDn == 0){
					payDn = m_maturityDate;
				}
				double rebDn = tmp.rebateDown ;

				if( payDn < m_maturityDate ){
					m_knockInRebate = rebDn / hZC->GetDiscountFactor(payDn, m_maturityDate);
				}
				else if( payDn > m_maturityDate ){
					m_knockInRebate = rebDn * hZC->GetDiscountFactor(m_maturityDate, payDn);
				}
				else{
					m_knockInRebate = rebDn ;
				}				
				break; // only one rebate...
			}
			if( tmp.isBarrierUp && !tmp.isKnockOutUp )
			{
				long payUp = tmp.payDateUp;
				if( payUp == 0){
					payUp = m_maturityDate;
				}
				double rebUp = tmp.rebateUp ;

				if( payUp < m_maturityDate ){
					m_knockInRebate = rebUp / hZC->GetDiscountFactor(payUp, m_maturityDate);
				}
				else if( payUp < m_maturityDate ){
					m_knockInRebate = rebUp * hZC->GetDiscountFactor(m_maturityDate, payUp);
				}
				else{
					m_knockInRebate = rebUp;
				}
				break; // only one rebate...
			}		
		}
	}
	if( fabs(m_knockInRebate) > 1e-6 )
		throw "Knock In rebates not supported";


	MlEqConstDateHandle	dateh = m_hUnderlying->GetDateHandle();
	double fwd = m_hUnderlying->GetForward(m_maturityDate, false);
	double mat = dateh->GetYearFraction(m_maturityDate);	
	double sdev = 3.5;
	double vol = m_hUnderlying->GetVolatility(MlEqStrike(fwd), m_maturityDate);		

	m_upper = fwd * exp( sdev*vol*sqrt(mat) );
	m_lower = fwd * exp( -1.5*sdev*vol*sqrt(mat) );

	m_nt = 100;	// number of steps a year..
	m_nx = int(m_nt * mat);

	m_loanspread_bump = 0.0001;	// 1 bp
}


/****************************************************************
**	Class  : BarrierPricer 
**	Routine: solve
**	Returns: 
**	Action : 
**           
****************************************************************/


void BarrierPricer::solve( CVector& result )
{
	initialCondition();	

	m_KnockOutInterpolator = m_hValueAtSlice ;
	if( m_priceKnockIn )
	{
		m_KnockInInterpolator = m_hValueAtSlice ;
//		addKnockInRebate( m_KnockOutInterpolator );
	}

	int nblock = m_barrierFeatures.getsize();

	for(int iblock=nblock-1; iblock>=0; --iblock)
	{
		buildTimeGrid( iblock );
	
		if( m_priceKnockIn )
		{
			m_hValueAtSlice = m_KnockInInterpolator;
			solveKnockInBlock( iblock, m_KnockInInterpolator );			
		}

		BarrierFeatures tmpFeatures = m_barrierFeatures[iblock];
		m_hValueAtSlice = m_KnockOutInterpolator;
		solveBlock( tmpFeatures, m_KnockOutInterpolator );
	}

	computeValueAtSlice( m_KnockOutInterpolator, m_KnockInInterpolator );	
	
	getGreeks( result );
}


/****************************************************************
**	Class  : BarrierPricer 
**	Routine: getGreeks
**	Returns: 
**	Action : 
**           
****************************************************************/


void BarrierPricer::getGreeks( CVector& result )
{
	result.resize(5);

	double eps = 0.001;
	double price	= m_hValueAtSlice->getValue(0.);
	double price_	= m_hValueAtSlice->getValue(eps);
	double _price	= m_hValueAtSlice->getValue(-eps);

	long nToday = m_hUnderlying->GetDateHandle()->GetDate();
	double spot = m_hUnderlying->GetSpot(nToday);

	double deltaX = (price_-_price) / (2.*eps);
	double gammaX = (price_+_price - 2.*price) / (eps*eps);

	result[0] = price;
	result[1] = deltaX / spot;
//	result[2] = ( gammaX - deltaX ) / (spot*spot) ;
	result[2] = ( gammaX - deltaX ) / spot * 0.01 ;	// delta move for 1% spot move

	double dvega = m_hValueAtSlice->getValue(0., 1);
	result[3] = dvega - price;

	double dlsbump = m_hValueAtSlice->getValue(0., 2);
	result[4] = dlsbump - price;
}



/****************************************************************
**	Class  : BarrierPricer 
**	Routine: solveBlock
**	Returns: 
**	Action : 
**           
****************************************************************/

void BarrierPricer::solveBlock(BarrierFeatures bfeatures, MlEqInterpolatorHandle& retInterpolator)
{
	double upper = m_upper;
	double lower = m_lower;	
	if( bfeatures.isBarrierUp )		upper = bfeatures.barrUp;	
	if( bfeatures.isBarrierDown )	lower = bfeatures.barrDown;

	int nx = int( m_nx * log(upper/lower)/log(m_upper/m_lower) ) ;

	DupireLocalVolHandle lv = new DupireLocalVol; 
	lv->initialize(m_hUnderlying, lower, upper, m_timeGrid,  nx, 0.01, 0);

	CVector tmp(1);

	RCPtr<knockOutBarrierHelper01> barrier 
		= new knockOutBarrierHelper01( &bfeatures, m_dfDn, m_dfUp, m_hValueAtSlice);

	pdeLocalVolBump pde;
	pde.init( lv,&(*barrier), m_loanspread_bump );
	pde.pde_integrate(tmp);

	barrier->finalize(0, &pde);
	retInterpolator = barrier->getFinalValue();
}

/****************************************************************
**	Class  : BarrierPricer 
**	Routine: solveKnockInBlock
**	Returns: 
**	Action : 
**           
****************************************************************/



void BarrierPricer::solveKnockInBlock(int iblock, MlEqInterpolatorHandle& KnockInInterpolator)
{
	/*	Here remove all the "Knock In" barriers and price only the remaining (if any)
		Knock Out barriers...next step for this block is to solve all the barriers as
		Out and get the difference of the 2 options as the actual price...the rebate 
		paid is the first to occur (you can have only one "In" rebate)
	*/

	BarrierFeatures tmpFeatures = m_barrierFeatures[iblock];

	if( m_barrierFeatures[iblock].isBarrierUp && m_barrierFeatures[iblock].isBarrierDown )
	{
		tmpFeatures.isBarrierDown = m_barrierFeatures[iblock].isKnockOutDown;
		tmpFeatures.isBarrierUp   = m_barrierFeatures[iblock].isKnockOutUp;
	}
	else if( m_barrierFeatures[iblock].isBarrierDown && !m_barrierFeatures[iblock].isKnockOutDown ) // no barrier up but KI
	{
		tmpFeatures.isBarrierDown = false;
	}
	else if( m_barrierFeatures[iblock].isBarrierUp && !m_barrierFeatures[iblock].isKnockOutUp ) // no barrier dn but KI
	{
		tmpFeatures.isBarrierUp = false;
	}
	
	solveBlock( tmpFeatures, KnockInInterpolator );	

}


/****************************************************************
**	Class  : BarrierPricer 
**	Routine: computeValueAtSlice
**	Returns: 
**	Action : 
**           
****************************************************************/


void BarrierPricer::computeValueAtSlice( MlEqInterpolatorHandle	KnockOutInterpolator,
										 MlEqInterpolatorHandle	KnockInInterpolator )
{
	if( !m_priceKnockIn )
	{
		m_hValueAtSlice = KnockOutInterpolator ;
	}
	else
	{
		CVector xo = KnockOutInterpolator->getXData();
		CVector xi = KnockInInterpolator->getXData();

		std::set<double> xval;
		int no = xo.getsize();
		int ni = xi.getsize();

		for(int i=0; i<no; ++i){
			xval.insert(xo[i]);
		}
		for(int i=0; i<ni; ++i){
			xval.insert(xi[i]);
		}

		int ntot = xval.size();
		CVector X(ntot), Y(ntot);

		int nd = KnockOutInterpolator->getNumberOfSlices();	// check...

		for(int id=0; id<nd; ++id)
		{

			std::set<double>::iterator it = xval.begin();

			for(int i=0; i<ntot; ++i)
			{
				double x = *it++;
				double y = KnockInInterpolator->getValue( x, id );
				y		-= KnockOutInterpolator->getValue( x, id );

				X[i] = x;
				Y[i] = y;	// std::max( y , 0.) ; // to avoid noise in degenerate cases...
			}

			m_hValueAtSlice->reinitialize( X, Y, id );

		}

	}

}


void BarrierPricer::addKnockInRebate( MlEqInterpolatorHandle&	KnockOutInterpolator )
{
	if( !m_priceKnockIn || fabs(m_knockInRebate) < 1e-6 )
	{
		m_knockInRebate = 0.0;
		return;		
	}
	else
	{
		CVector x = KnockOutInterpolator->getXData();
		int nx = x.getsize();
		CVector y(nx);

		int nd = KnockOutInterpolator->getNumberOfSlices();

		for(int id=0; id<nd; ++id)
		{
			for(int i=0; i<nx; ++i){
				y[i] = KnockOutInterpolator->getValue( x[i], id ) - m_knockInRebate;
			}

			KnockOutInterpolator->reinitialize( x, y, id );
		}		
	}

}


/****************************************************************
**	Class  : BarrierPricer 
**	Routine: buildTimeGrid
**	Returns: 
**	Action : 
**           
****************************************************************/

void BarrierPricer::buildTimeGrid(int iblock)
{
	long startDate	= m_barrierFeatures[iblock].startDate;
	long endDate	= m_barrierFeatures[iblock].endDate;
	long payDateDn	= m_barrierFeatures[iblock].payDateDn;
	long payDateUp	= m_barrierFeatures[iblock].payDateUp;

	bool dpu = (payDateUp > 0);
	bool dpd = (payDateDn > 0);

	MlEqConstDateHandle	dateh = m_hUnderlying->GetDateHandle();

	double T = dateh->GetYearFraction(endDate) - dateh->GetYearFraction(startDate) ;
	int nt = int(T*m_nt) + 1;
	nt = std::max( nt, 2 );

	double dt = double(endDate-startDate) / (nt-1.);
	long  date = startDate;

	m_timeGrid.resize(nt);
	m_dfUp.resize(nt, 1.);
	m_dfDn.resize(nt, 1.);
	MlEqZeroCurveHandle hZeroCurve = m_hUnderlying->GetPayZeroCurve(true);

	for(int i=0; i<nt; ++i)
	{
		date = long( startDate + i*dt );
		m_timeGrid[i] = date;

		if(dpu){
			m_dfUp[i] = hZeroCurve->GetDiscountFactor(date, payDateUp);	
		}
		if(dpd){
			m_dfDn[i] = hZeroCurve->GetDiscountFactor(date, payDateDn);	
		}
	}
}



/****************************************************************
**	Class  : BarrierPricerTest 
**	Routine: BarrierPricerTest
**	Returns: 
**	Action : 
**           
****************************************************************/


BarrierPricerTest::BarrierPricerTest(MlEqAssetHandle			hUnderlying,
									 GVector<BarrierFeatures>	barrierFeatures,
									 long	maturityDate, double strike, double cp)
:BarrierPricer(hUnderlying,	barrierFeatures, maturityDate),
m_strike(strike),
m_cp(cp)
{
}



/****************************************************************
**	Class  : BarrierPricerTest 
**	Routine: initialCondition
**	Returns: 
**	Action : 
**           
****************************************************************/


void BarrierPricerTest::initialCondition()
{
	// To do: use the barrier boundary conditions in the KO case...
	
	long nToday = m_hUnderlying->GetDateHandle()->GetDate();
	double spot = m_hUnderlying->GetSpot(nToday);

	double upper = m_upper, lower = m_lower;

	int nblock = m_barrierFeatures.getsize();
	BarrierFeatures feat = m_barrierFeatures[nblock-1] ;
	bool addRebate = false;

	int nx = m_nx/2;
	if( addRebate ) nx++;

	CVector u(nx), x(nx);

	if( m_cp > 0. )	// call 
	{
		double divide = addRebate? (nx-3.):(nx-2.);
		double dx = log(upper/m_strike) / divide;
		x[0] = log(lower/spot);			u[0] = 0.;
		x[1] = log(m_strike/spot);		u[1] = 0.;

		for(int i=2; i<nx; ++i)
		{
			x[i] = x[1] + (i-1) * dx ;
			u[i] = spot * exp(x[i]) - m_strike ;
		}

		if( addRebate ) {
			u[nx-1] = u[nx-2] = feat.rebateUp;
		}
	}
	else
	{
		double divide = addRebate? (nx-3.):(nx-2.);
		double dx = -log(lower/m_strike) / divide;
		x[nx-1] = log(upper/spot);			u[nx-1] = 0.;
		x[nx-2] = log(m_strike/spot);		u[nx-2] = 0.;

		for(int i=3; i<=nx; ++i){
			x[nx-i] = x[nx-2] - (i-2) * dx ;
			u[nx-i] = m_strike - spot * exp(x[nx-i]);
		}

		if( addRebate ) {
			u[0] = u[1] = feat.rebateDown;
		}
	}


	int nd = 3;
	GVector<CVector> X(1), U(nd);
	X[0] = x;
	for(int i=0; i<nd; ++i)	{
		U[i] = u;
	}

	int			addTanhWings = 1;	// as in mc model info
	double		cL = 0.01;	
	double		cR = 3.;			
	int			yPower = 2;	
	
	m_hValueAtSlice = new MlEqMonotonicSplineInterpolator( 	X, U, addTanhWings, cL, cR, yPower );
}




