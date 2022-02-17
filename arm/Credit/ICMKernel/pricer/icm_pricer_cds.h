#ifndef _ICM_Pricer_CDS_H
#define _ICM_Pricer_CDS_H

#include "ICMKernel/pricer/icm_pricer_security.h"
#include "ICMKernel/crv/icm_defaultcurve.h"

/*********************************************************************************/
/*! \class  ICM_Pricer_Cds icm_pricer_cds.h "icm_pricer_cds.h"
 *  \author Damien Pouponneau 
 *	\version 1.0
 *	\date   June 2004 
 *	\file   icm_pricer_cds.h
 *	\brief  Pricer for ICM_Cds security */
/***********************************************************************************/

class ICM_Poly; 

class ICM_Pricer_Cds : public ICM_Pricer_Security
{
	friend void objfun_cds(void*,double,double&) ; 
	friend void objfun_cds_acc(void*,double,double&) ; 
	class disc_index 
	{
	public:
		double lag; 
		qCredit_Leg_Style legStyle ;
	public:
		bool operator<(const disc_index&ref) const
		{
			if (lag<ref.lag) return true; 
			else if (legStyle<ref.legStyle) return true; 
			return false; 
		}
	}; 
	typedef std::map<disc_index,ARM_Vector> discs_t ;
private :
	
	double		itsRefNotional;				//we pay (itsRefNotional-Recovery)*Notional in default case
											//default is itsRefNotional=1
	ARM_ZeroCurve*	itsDiscountCurve;		//  association - Discount Curve
	ARM_ZeroCurve*	itsCouponCurve;			//	association - Curve for rates computation
	const ICM_DefaultCurve* itsDefaultCurve;		//	association - 
	discs_t  itsIntegrationSchedules; 
public :
	ICM_Pricer_Cds(void) { Init();}
	virtual ~ICM_Pricer_Cds() {}
	virtual void PropagateModel(ARM_Model *model) 
	{
		LoadMarketData(model,model->GetAsOfDate()) ;
	}

	void Init() ;

	void Set(ARM_Security*sec,ARM_Object*mod,const ICM_Parameters&params,const ARM_Date&asof) ;

protected:
	virtual double FeeLegPV ();
	virtual double DefLegPV ();
	virtual double Accrued();
	virtual double ComputeDuration(void);
private:
	// former production method, summit like
	double DefLegPV_Summit();	
	double FeeLegPV_Summit();
private:
	// new production method, non summit like
	double DefLegPV_NonSummit();
	double FeeLegPV_NonSummit();
private:
	// former "non summit" pricing. 
	double DefLegPV_Old();
	double FeeLegPV_Old();
public:

	virtual double ExpectedLossTranche(const double& yearterm,vector<double>& losses)
	{ return 1. - itsDefaultCurve->SurvivalProba(yearterm);}

	virtual void Reset(void)
	{
		ICM_Pricer_Security::Reset();
	}

	// API Implementation 
	virtual double ComputeSpread(const double& MtM = 0.);
	virtual double Compute_Fwd_Spread(const ARM_Date & Mty, const ARM_Date & CDS_ExpiryDate, double& FlatRbp_Fwd);
	

	protected:
	virtual double ComputeSensitivity(const qSENSITIVITY_TYPE& typesensi, 
		const std::string& plot,
		const std::string& label,
		double epsilon, double epsilonGamma =0 ); 
	public:


	private:
		void getSchedule_Summit(const qCredit_Leg_Style& Leg_Type, 
									const double& lag , ARM_Vector& res);
		void GenerateMergeSchedule_Summit(const qCredit_Leg_Style& Leg_Type, 
									const double& lag , ARM_Vector& res);
	public:

	void CptPolyIrCurve(ICM_Poly& poly,double& alpha, double yf1,double yf2,double Lambda = 0.0, double lag=0.);

	// To deal with ZC
	// void	DummyIRInterpolation(ARM_ZeroCurve*& TheIRCurve, double x, double& TheZCValue);


	// ---------------------------------------------------------------
	// Full Taylor Expansion for expression Exp(a0+a1*X+a2*X²) on [t1,t2] returns alpha & Px 
	// for Exp(alpha*X)*Px
	// ---------------------------------------------------------------
	// inline void TaylorExpansion(const double& a0,
	// 							const double& a1,
	// 							const double& a2,
	// 							const double& t1,
	// 							const double& t2,
	// 							double& alpha,
	// 							ICM_Poly& Px);
	
	virtual ICM_Pricer* CloneOnlyParams(ARM_Security* sec,ARM_Model* mod) ;

	inline double GetRefNotional() {return itsRefNotional;}
	inline void SetRefNotional(const double& value) {itsRefNotional=value;}

	void LoadMarketData(ARM_Object* object,const ARM_Date&AsOf) ;

};




#endif

