
#ifndef __PDEPRODUCTS_
#define __PDEPRODUCTS_

#include "mleqpde.h"





class plainVanillaOptionHelper : public MlEqPdeHelper
{
public:
		plainVanillaOptionHelper(){}

		plainVanillaOptionHelper( CVector& CallPut,GVector<MlEqStrikeHandle>& strike,double spot);
		~plainVanillaOptionHelper(){}


//		double m_CallPut;
//		double m_strike;

		CVector m_CallPut;
		CVector m_strikes;

		virtual void initial_condition(void* ptr,int nd,int nx,CMatrix& u ,MlEqPdeDriver* pde);
		virtual void boundary_condition(int id,int it,double t, void* ptr,pde_boundary_condition& lower,pde_boundary_condition& upper,MlEqPdeDriver* pde);

		void initialize(CVector& CallPut,GVector<MlEqStrikeHandle>& strike,double spot);

};



class americanOptionHelper : public MlEqPdeHelper
{
public:
		americanOptionHelper(){}
		americanOptionHelper(double spot,double CallPut,double strike);
		~americanOptionHelper(){}


		double m_CallPut;
		double m_strike;

		virtual PdeEnums::PdeMode strategy_function(int it, double t,  void* ptr, int nd, int nx, CMatrix& u_it,MlEqPdeDriver* pde);
		virtual void initial_condition(void* ptr,int nd,int nx,CMatrix& u ,MlEqPdeDriver* pde);
		virtual void boundary_condition(int id,int it,double t, void* ptr,pde_boundary_condition& lower,pde_boundary_condition& upper,MlEqPdeDriver* pde);
		void initialize(double spot,double CallPut,double strike);
};



class knockOutBarrierHelper : public MlEqPdeHelper
{
protected:
		double m_logBarrUp;
		double m_logBarrDo;
		CVector m_df;

		bool	m_isUpBarrier;
		bool	m_isLowBarrier;

public:
		double m_lowSpot;
		double m_highSpot;

		double m_UpperBarrier;
		double m_LowerBarrier;
		double m_UpperRebate;
		double m_LowerRebate;
		bool   m_payEnd;
		double m_CallPut;
		double m_strike;


//		virtual PdeEnums::PdeMode strategy_function(int it, double t,  void* ptr, int nd, int nx, CMatrix& u_it,MlEqPdeDriver* pde);
		virtual void initial_condition(void* ptr,int nd,int nx,CMatrix& u ,MlEqPdeDriver* pde);
		virtual void boundary_condition(int id,int it,double t, void* ptr,pde_boundary_condition& lower,pde_boundary_condition& upper,MlEqPdeDriver* pde);
		void initialize(double UpperBarrier,double LowerBarrier,
						double UpperRebate,double LowerRebate,
						double spot,double CallPut,double strike,
						bool payEnd, const CVector& df);

public:
	knockOutBarrierHelper(){}
	knockOutBarrierHelper(double UpperBarrier,double LowerBarrier,
						double UpperRebate,double LowerRebate,
						double spot,double CallPut,double strike,
						bool payEnd, const CVector& df);
	~knockOutBarrierHelper(){}

};





class callableNote00Helper : public MlEqPdeHelper
{
protected:
		CVector m_coupons;
		GVector<bool> m_isCouponDate;

public:
		double m_factor;
		double m_lowSpot;
		double m_highSpot;

		double m_fixed_coupon;
		double m_gearing;
		double m_CallPut;
		double m_strike;


		virtual PdeEnums::PdeMode strategy_function(int it, double t,  void* ptr, int nd, int nx, CMatrix& u_it,MlEqPdeDriver* pde);
		virtual void initial_condition(void* ptr,int nd,int nx,CMatrix& u ,MlEqPdeDriver* pde);
		virtual void boundary_condition(int id,int it,double t, void* ptr,pde_boundary_condition& lower,pde_boundary_condition& upper,MlEqPdeDriver* pde);

		void initialize(double spot, 
						double refSpot,
						MlEqConstDateHandle hDate, 
						const CVector& coupons, 
						const GVector<long>& couponDates, 
						const GVector<long>& pdeDates,
						double CallPut,double strike, double gearing, double fixed_coupon);
};



class callableNote01Helper : public MlEqPdeHelper
{
protected:
		CVector m_coupons;
		GVector<bool> m_isCouponDate;

		CVector m_rebates;
		GVector<bool> m_isCallableDate;

public:
		double m_factor;
		double m_lowSpot;
		double m_highSpot;

		double m_fixed_coupon;
		double m_gearing;
		double m_CallPut;
		double m_strike;


		virtual PdeEnums::PdeMode strategy_function(int it, double t,  void* ptr, int nd, int nx, CMatrix& u_it,MlEqPdeDriver* pde);
		virtual void initial_condition(void* ptr,int nd,int nx,CMatrix& u ,MlEqPdeDriver* pde);
		virtual void boundary_condition(int id,int it,double t, void* ptr,pde_boundary_condition& lower,pde_boundary_condition& upper,MlEqPdeDriver* pde);

		void initialize(double spot, 
						double refSpot,
						MlEqConstDateHandle hDate, 
						const CVector& coupons, 
						const GVector<long>& couponDates, 
						const CVector& rebates, 
						const GVector<long>& exerciseDates, 
						const GVector<long>& pdeDates,
						double CallPut,double strike, double gearing, double fixed_coupon);
};





struct BarrierFeatures
{
	double	barrUp;
	double	barrDown;
	double	rebateUp;
	double	rebateDown;
	bool	isBarrierUp;
	bool	isBarrierDown;
	bool	isKnockOutUp;
	bool	isKnockOutDown;

	long	payDateUp;
	long	payDateDn;

	long	startDate;
	long	endDate;

	BarrierFeatures()
	{
		barrUp = barrDown = rebateUp = rebateDown = 0.;
		isBarrierUp = isBarrierDown = false;
		isKnockOutUp = isKnockOutDown = true;
		payDateUp = payDateDn = startDate = endDate = 0L;
	}
};


class knockOutBarrierHelper01 : public MlEqPdeHelper	// this will be a pde slice...
{
protected:
		double m_logBarrUp;
		double m_logBarrDo;
		CVector m_dfDn;
		CVector m_dfUp;

		bool	m_isUpBarrier;
		bool	m_isDownBarrier;

		MlEqInterpolatorHandle m_hFinalValue;

public:
		double m_lowSpot;
		double m_highSpot;

		double m_UpperBarrier;
		double m_LowerBarrier;
		double m_UpperRebate;
		double m_LowerRebate;
		double m_CallPut;
		double m_strike;


		virtual void initial_condition(void* ptr,int nd,int nx,CMatrix& u ,MlEqPdeDriver* pde);
		virtual void boundary_condition(int id,int it,double t, void* ptr,pde_boundary_condition& lower,pde_boundary_condition& upper,MlEqPdeDriver* pde);
		void initialize(double UpperBarrier,double LowerBarrier,
						double UpperRebate,double LowerRebate,
						bool isUp, bool isDown,
						const CVector& dfDn,const CVector& dfUp,
						MlEqInterpolatorHandle hFinalValue);

		void finalize( int id, MlEqPdeDriver* pde );

		MlEqInterpolatorHandle getFinalValue() const { return m_hFinalValue;}

public:
	knockOutBarrierHelper01(){}
	knockOutBarrierHelper01(double UpperBarrier,double LowerBarrier,
							double UpperRebate,double LowerRebate, bool isUp, bool isDown,
							const CVector& dfDn, const CVector& dfUp,
							MlEqInterpolatorHandle hFinalValue);
	knockOutBarrierHelper01(BarrierFeatures* bfeat, const CVector& dfDn, const CVector& dfUp, MlEqInterpolatorHandle hFinalValue);

	~knockOutBarrierHelper01(){}

};

struct PdeFullSolution
{
	GVector<CVector>	u;
	GVector<CVector>	x;
	GVector<long>		t;
};

class BarrierPricer : public RCObject	
{
protected:
		RCPtr< knockOutBarrierHelper01 >	m_barrierSlice;		
		MlEqInterpolatorHandle				m_hValueAtSlice;

		MlEqInterpolatorHandle				m_KnockInInterpolator ;
		MlEqInterpolatorHandle				m_KnockOutInterpolator ;

		long								m_maturityDate;
		GVector<BarrierFeatures>			m_barrierFeatures;
		GVector<long>						m_timeGrid;
		CVector								m_dfUp;
		CVector								m_dfDn;
		double								m_knockInRebate;
		double								m_lower;
		double								m_upper;
		int									m_nx;
		int									m_nt;	// number of steps a year
		MlEqAssetHandle						m_hUnderlying;
		PdeFullSolution						m_u_x_t ;	
		bool								m_storeResult;
		bool								m_priceKnockIn;

		void			buildTimeGrid(int slice);
public:
		void			getGreeks( CVector& result );
		double			m_loanspread_bump;	// add whatever you like...

protected:
		void			solveBlock(BarrierFeatures bfeatures, MlEqInterpolatorHandle& retInterpolator);
		void			solveKnockInBlock(int iblock, MlEqInterpolatorHandle& KnockInInterpolator);
		void			computeValueAtSlice(MlEqInterpolatorHandle	ko, MlEqInterpolatorHandle	ki);
		void			addKnockInRebate( MlEqInterpolatorHandle&	KnockOutInterpolator );

		virtual void	initialCondition() = 0;

public:
		BarrierPricer(){}
		~BarrierPricer(){}
		BarrierPricer(	 MlEqAssetHandle			hUnderlying,
						 GVector<BarrierFeatures>	barrierFeatures,
						 long	maturityDate);

		void initialize( MlEqAssetHandle			hUnderlying,
						 GVector<BarrierFeatures>	barrierFeatures,
						 long	maturityDate);

		void solve( CVector& result );

		const PdeFullSolution& getFullResult(){ return m_u_x_t ;}
};


class BarrierPricerTest	:	public	BarrierPricer
{
public:
	BarrierPricerTest(){}
	~BarrierPricerTest(){}
	BarrierPricerTest(	 MlEqAssetHandle			hUnderlying,
						 GVector<BarrierFeatures>	barrierFeatures,
						 long	maturityDate, double strike, double cp);
protected:
	double m_strike;
	double m_cp;

	virtual void initialCondition();
};





#endif