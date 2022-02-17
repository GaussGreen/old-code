#ifndef _ARMBSCONVADJUSTREP_H
#define _ARMBSCONVADJUSTREP_H

#include "bsconvadjust.h"
#include "volcurv.h"
#include "volint.h"
#include "gpclosedforms/inverse.h"
#include "gpclosedforms/gaussian_integrals.h"

using ARM::DoubleToDoubleFunc;
using ARM::LegendreInt;
using ARM::brentSolve;

inline double KSI(double x) {return fabs(x) < K_DOUBLE_TOL ? 1. : (1. - exp(-x))/x;};

class ARM_CMSVolInterSpline
{
private:
	ARM_Vector	itsX;
	ARM_Vector	itsY;

	ARM_Vector	itsMK;
	ARM_Vector	itsAK;
	ARM_Vector	itsBK;

	bool		itsIsInit;

public:
	ARM_CMSVolInterSpline(void) : itsX(0), itsY(0)
	{
		itsIsInit = false;
	}

	~ARM_CMSVolInterSpline(void)
	{}

	ARM_CMSVolInterSpline(const ARM_Vector& x, const ARM_Vector& y);

	ARM_CMSVolInterSpline(const ARM_CMSVolInterSpline& rhs);

public:

	double operator()(double x) const;

private:

	void	buildSpline();
	void	findMK(double * fk, double * hk);
	double	LeftExtrapol(double x) const;
	double	RightExtrapol(double x) const;
};

inline double ARM_CMSVolInterSpline::operator ()(double x) const
{
	if(itsIsInit == false)
		throw Exception(__LINE__, __FILE__, ERR_INVALID_OPERATION,
                        " ARM_CMSVolInterSpline : spline not initialised");


	if(x < itsX[0]) return LeftExtrapol(x);
	if(x > itsX[itsX.size()-1]) return RightExtrapol(x);

	int idx = itsX.LookupPrevIndex(x);
	
	double xd = x - itsX[idx], xu = itsX[idx+1] - x;

	return itsMK[idx+1] * xd*xd*xd / (6. * (itsX[idx+1] - itsX[idx]))
		+ itsMK[idx] * xu*xu*xu / (6. * (itsX[idx+1] - itsX[idx]))
		+ itsAK[idx] * xd + itsBK[idx];
}

inline double ARM_CMSVolInterSpline::LeftExtrapol(double x) const
{
	return itsAK[0] * (x - itsX[0]) + itsBK[0];
}

inline double ARM_CMSVolInterSpline::RightExtrapol(double x) const
{
	return itsY[itsX.size()-1] + itsAK[itsX.size()-2] * (x - itsX[itsX.size()-1]);
}

class ARM_BSConvAdjustRep : public ARM_BSConvAdjust
{
private:
	
	ARM_VolLInterpol 			itsCvxCMSVol;	// les vols pour l'ajustement de convexité naturel
	vector<ARM_VolLInterpol>	itsCMSVol;
	ARM_Vector					itsStdDev;
	int							itsIdxATM;
	int							itsNbPointsForRepliq;
	ARM_Model *					itsUsedModel;
	ARM_VolLInterpol			itsMR;
	bool						itsFullRepliq;
	double						itsUpperProba;
	double						itsLowerProba;

public:
	ARM_BSConvAdjustRep(ARM_Model* UsedModel, const ARM_VolLInterpol& swoptVolCurve, const ARM_Vector& stddev, const ARM_VolLInterpol& MR, int NbPointsForRepliq = 8, bool fullRepliq = false,
						double upperProba = 0.999, double lowerProba = 0.001);// Mode: K_CONV_ADJ_EXP (ARM formulae)
    ARM_BSConvAdjustRep(const ARM_BSConvAdjustRep& rhs);
    ARM_BSConvAdjustRep& operator= (const ARM_BSConvAdjustRep& rhs);


	virtual ~ARM_BSConvAdjustRep(void);
        
 
	// return the natural convexity adjustment of CMS 
	virtual double NaturalAdjstCMS(ARM_Model* Model, 
                               const NaturalAdjustData& Input, 
                               StoreFwdRateInfo* StoreInfo = NULL);

	// compute the cms volatility used in the cms adjustment
	virtual double ComputeCMSVolatility(ARM_Model* Model, const NaturalAdjustData& Input, 
                                    StoreFwdRateInfo* StoreInfo = NULL);

private:

	void doColReplication(int col);

	double	ImpCMSVol(double fwd, double mat, double strike, int CallPut, double price);

	class CMSReplicator : public DoubleToDoubleFunc
	{
	private:
		ARM_Model *		itsUsedModel;
		double			itsResetDate;
		double			itsStartDate;
		double			itsEndDate;
		double			itsFwdRate;
		double			itsAnnuity;
		int				itsCallPut;
		ARM_Vector		itsZCFwd;
		ARM_Vector		itsvarZC;
		ARM_Vector		itsracVarZC;
		ARM_Vector		itsDayCount;
		ARM_Vector		itsTheta;

		double			itsAlpha;
		double			itsBeta;
		double			itsMR;
		
		double			itsUBound;
		double			itsLBound;

	public:
		CMSReplicator(ARM_Model * model, ARM_Date& StartDate, ARM_Date& EndDate, int DayCount, int Freq, 
			int CallPut, double mr, double upperProba, double lowerProba, double strike = 0.);

		double	operator()(double x) const;

		double	limitFromStrike(double strike);

		double	extraTerms(double xlim, double strike) const;

		void	setType(int CallPut) {itsCallPut = CallPut;};

		void	setCoeffs(double alpha, double beta) {itsAlpha = alpha; itsBeta = beta;};

		double	GetHWVol(double target, double strike);

		double	GetFwdRate() const {return itsFwdRate;};

		double	GetUBound() const {return itsUBound;};
		double	GetLBound() const {return itsLBound;};

		double	BoundFromProba(double proba);

		class SwapFromFactor : public DoubleToDoubleFunc
		{
		private:
			CMSReplicator * _this;
			const ARM_Vector *	itsVars;
			const ARM_Vector *	itsracVars;

		public:
			SwapFromFactor(CMSReplicator * replicator, const ARM_Vector * vars, const ARM_Vector * racVars) : _this(replicator),
				itsVars(vars), itsracVars(racVars)
			{}

			double operator()(double x) const;
		};

		class ReplicationCoeffs
		{
		public:
			double	itsSwapRate;
			double	itsDerivSwapRate;
			double	itsBeta;
			double	itsRatio;
			double	itsDerivRatio;

		public:
			ReplicationCoeffs(double x, const ARM_Vector& zcFwd, const ARM_Vector& varZC, const ARM_Vector& racVarZC, const ARM_Vector& dcf, double a, double b);
		};

		class HWSwopt : public DoubleToDoubleFunc
		{
		private:
			CMSReplicator * _this;
			double			itsMR;
			double			itsMat;
			double			itsStrike;
			ARM_Vector		itsFlows;
			ARM_Vector		itsVars;
			ARM_Vector		itsracVars;
			int				itsCallPut;

		public:
			HWSwopt(CMSReplicator * replicator, double MR, double Mat, double K, int CallPut);

			double operator()(double x) const;
		};

		class EquivStrikeForProba : public DoubleToDoubleFunc
		{
		private:
			double	itsReset;
			double	itsStart;
			double	itsEnd;
			double	itsFwd;
			double	itsAnnuity;

			ARM_Model * itsUsedModel;

		public:
			EquivStrikeForProba(double reset, double start, double end, double fwd, double annuity, ARM_Model * model) :
			  itsReset(reset), itsStart(start), itsEnd(end), itsFwd(fwd), itsAnnuity(annuity), itsUsedModel(model)
			  {}

			double operator()(double x) const
			{
				double bin = itsUsedModel->EuroSwaption(itsReset, itsStart, itsEnd, -1, itsFwd, x + 0.0001, itsAnnuity)
						   - itsUsedModel->EuroSwaption(itsReset, itsStart, itsEnd, -1, itsFwd, x - 0.0001, itsAnnuity);

				return fabs(bin) / 2. / 0.0001 / itsAnnuity;
			}
		};

		friend class ReplicationCoeffs;
		friend class SwapFromFactor;
		friend class HWSwopt;
	};

	class EquivConvAdjVol : public DoubleToDoubleFunc
	{
	private:
		double	itsFwd;
		double	itsTheta;
		double	itsN;
		double	itsResetLag;

	public:
		EquivConvAdjVol(double fwd, double theta, double N, double resetLag) : 
		  itsFwd(fwd), itsTheta(theta), itsN(N), itsResetLag(resetLag)
		{}

		double operator()(double x) const
		{
			double lambda=1.0-(itsTheta*itsN*itsFwd)/((1.0+itsTheta*itsFwd)*(pow((1.0+itsTheta*itsFwd),itsN)-1.0));
			double exp1=exp(x*x*itsResetLag)-1.0;
			double ki = lambda*exp1;
			return itsFwd*(ki+1.);
		}
	};

	class EquivCMSVol : public DoubleToDoubleFunc
	{
	private:
		double	itsFwd;
		double	itsResetLag;
		double	itsStrike;
		int		itsCallPut;

	public:
		EquivCMSVol(double fwd, double reset, double strike, int callput) : 
		  itsFwd(fwd), itsResetLag(reset), itsStrike(strike), itsCallPut(callput)
		{}

		double operator()(double x) const
		{
			return bsOption(itsFwd, itsStrike, x, 0., 0., itsResetLag, itsCallPut);
		}
	};
};


#endif