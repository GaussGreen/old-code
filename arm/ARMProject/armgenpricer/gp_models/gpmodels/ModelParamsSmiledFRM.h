/*
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 */

#ifndef _INGPMODELS_MODELPARAMSSMILED_H
#define _INGPMODELS_MODELPARAMSSMILED_H

#include "gpinfra/modelparams.h"
#include "gpinfra/curvemodelparam.h"
#include "gpinfra/modelparamsvec.h"
#include "gpbase/curve.h"
#include "typedef.h"

CC_BEGIN_NAMESPACE( ARM )


class ARM_ModelParamsSmiled : public ARM_ModelParams 
{
public:
	enum CalibProxy
	{
		ATM = 0,
		AtmBlack,
		MomentMatching,
		LocalVolatility,
		LocalVolatilityBlack,
		LocalVolatilityWithRescaling,
		LocalVolatilityBlackWithRescaling,
		EffectiveSkew,
		EffectiveSkewWithRescaling,
		GaussBasketAtm,
		GaussBasketMoneyness
	};
	enum CorrelType
	{
		Beta = 0,
		Theta,
		BetaWithRecorrel,
		CorrelMatrix,
		Fwd
	};
private :	
	void	BuildVolatilityStructure();
	void	UpdateVolatilityFromHump();
	void	UpdateVolatilityFromVolatilites();
	void	UpdateIntegratedCovarianceModel();
	void	UpdateIntegratedCovarianceCalib();

	double	InstantaneousCorrelation(double T1, double T2) const;

	double	IntegratedCovarianceFromZero(	double t, 
											size_t i, 
											size_t j ) const;


public:
	ARM_ModelParamsSmiled(	const ARM_ModelParamVector& params, 
								size_t dim = 1, 
								CorrelType correlTYpe = Beta, 
								bool allowInterpol=false,
								double recorrel=0.);

	ARM_ModelParamsSmiled( const ARM_ModelParamsSmiled& rhs);
	virtual ~ARM_ModelParamsSmiled();

	inline void SetNbFactors(size_t nbFactors) { itsNbFactors = nbFactors; }
	inline	CorrelType GetCorrelType() const { return itsCorrelType; }
	inline void SetCorrelType(CorrelType correlType) { itsCorrelType = correlType; }
	inline void SetAllowInterpol(bool allowInterpol) { itsAllowInterpol = allowInterpol; }

	inline virtual size_t		FactorCount() const				{ return itsNbFactors;}
	inline virtual void			SetFactorCount(size_t fc)		{ itsNbFactors  = fc;}

	inline double				meanReversion(size_t i) const   { return 0.;	}

	inline ARM_Curve*			GetVolCurve(size_t i) const		 { return itsVol[i];	}
	inline ARM_Curve*			GetVolCurveCalib(size_t i) const { return itsVolCalib[i];	}

	inline ARM_Curve*			GetCorrelCurve() const			{ return ((ARM_CurveModelParam&)GetModelParam(ARM_ModelParamType::BetaCorrelation)).GetCurve();	}
	inline bool					UsingRecorrel() const			{ return DoesModelParamExist(ARM_ModelParamType::ReCorrelation); }
	inline ARM_Curve*			GetReCorrelCurve() const		{ return ((ARM_CurveModelParam&)GetModelParam(ARM_ModelParamType::ReCorrelation)).GetCurve();	}

	virtual void		PreProcessing(ARM_ModelFitter& modelFitter,int factorNb=0) {};
    virtual void		PostProcessing(const ARM_ModelFitter& modelFitter,ARM_PricingModel* model,int factorNb=0 ) {};
	void				setVolatilitiesAndCorrelMatrix(ARM_GP_Matrix& volatilities, ARM_MatrixVector& matrixVector );
	void				setVolatilities(ARM_GP_Matrix& volatilities);
	void				setResetDates( const ARM_GP_VectorPtr& resetDates);
	void				AdviseCurrentCalib(double from, double to);
	void				ForOptimize();
    void				AdviseSwaptionApprox(CalibProxy proxy)	;
    
	virtual iterator	SetModelParamValue( int paramType, size_t i, double value, double time, double tenor = 0.0, size_t factorNb=0 );

	double	GetEquivalentVolCalib(double evalTime, double swapResetTime, size_t i, size_t j, const std::vector<double>& weight, const std::vector<double>& coeff) const;
	double	GetEquivalentVolModel(double evalTime, double swapResetTime, size_t i, size_t j, const std::vector<double>& weight, const std::vector<double>& coeff) const;
	double	GetEquivalentVol(double evalTime, double swapResetTime, size_t i, size_t j, const std::vector<double>& weight, const std::vector<double>& coeff, const ARM_MatrixVector& cov) const;
	double	GetEquivalentVol(double evalTime, double swapResetTime, size_t i, size_t j, const std::vector<double>& weight, const std::vector<double>& coeff, CalibProxy proxy) const;

	double	GetEquivalentCov(double evalTime, double swapResetTime, const std::vector<double>& weight1, const std::vector<double>& weight2, const std::vector<double>& coeff, const ARM_MatrixVector& cov) const;
	double	GetEquivalentCov(double evalTime, double swapResetTime, const std::vector<double>& weight1, const std::vector<double>& weight2, const std::vector<double>& coeff, CalibProxy proxy) const;


	double	IntegratedCovariance(	double s, 
									double t, 
									size_t i, 
									size_t j ) const;

	double	IntegratedVariance( double s, double t, size_t i) const {return IntegratedCovariance(s, t, i, i);};

	double	StationaryConstraint( double s, double t);

	/// Standard ARM object support
	virtual ARM_Object* Clone() const	{ return new ARM_ModelParamsSmiled(*this);};
	virtual string toString(const string& indent="",const string& nextIndent="") const;

	void	PreComputeIntegratedCovariance();

	const std::vector<double>&	GetGlobCovTimes() const {return itsGlobalCovTimes;};

	double	InstantCorrel(	double t,
							size_t i,
							size_t j,
							size_t index) const;
	
	double	InstantCorrel(double t, int i, int j) const;
	
private:
	typedef vector<ARM_Curve*>	SBGM_VolCurveVector;

	SBGM_VolCurveVector			itsVol;
	SBGM_VolCurveVector			itsVolCalib;
	ARM_VectorPtr				itsResetDates;
	ARM_MatrixVector			itsGlobalCovariance;
	ARM_MatrixVector			itsGlobalCovarianceCalib;
	std::vector<double>&				itsGlobalCovTimes;
	size_t						itsNbFwds;
	size_t						itsNbFactors;
	size_t						itsIndexFrom;
	size_t						itsIndexTo;

	ARM_ModelParamsSmiled::CorrelType	itsCorrelType;
	ARM_GP_MatrixPtr		itsVolatilities;
	ARM_MatrixVector			itsCorrelMatrix;
	double						itsRecorrel;

	bool						itsAllowInterpol;
	bool						itsUseModifiedCovForCalib;

};


inline	size_t IndexOfLastLowerEqInVector(double x, const std::vector<double>& v)
{
	size_t size = v.size();

	if (v.size()>0)
	{
		size_t index=0;
		if (x >= v[size-1])
			index = size-1;
		else
		{
			while (v[index]<=x)
				index++;
			index-=1;
		}
		if (index==-1)
			ARM_THROW( ERR_INVALID_ARGUMENT, "ARM_ModelParamsSmiledFRM::IndexOfLastLowerEqInVector : strictly lower than first");
		return index;
	}
	else
		ARM_THROW( ERR_INVALID_ARGUMENT, "ARM_ModelParamsSmiledFRM::IndexOfLastLowerEqInVector : vector of size 0");
};

inline	size_t IndexOfFirstHigherEqInVector_DefaultLast(double x, const std::vector<double>& v)
{
	size_t size = v.size();

	if (v.size()>0)
	{
		size_t index=0;
		if (x > v[size-1])
			index = size-1;
		else
		{
			while (v[index]<x)
				index++;
		}
		return index;
	}
	else
		ARM_THROW( ERR_INVALID_ARGUMENT, "ARM_ModelParamsSmiledFRM::IndexOfFirstHigherEqInVector : vector of size 0");
};

inline	size_t IndexOfFirstHigherInVector_DefaultLast(double x, const std::vector<double>& v)
{
	size_t size = v.size();

	if (v.size()>0)
	{
		size_t index=0;
		if (x >= v[size-1])
			index = size-1;
		else
		{
			while (v[index]<=x)
				index++;
		}
		return index;
	}
	else
		ARM_THROW( ERR_INVALID_ARGUMENT, "ARM_ModelParamsSmiledFRM::IndexOfFirstHigherEqInVector : vector of size 0");
};

inline	size_t IndexOfLastLowerEqInVector_DefaultFirst(double x, const std::vector<double>& v)
{
	size_t size = v.size();

	if (v.size()>0)
	{
		size_t index=0;
		if (x >= v[size-1])
			index = size-1;
		else
		{
			while (v[index]<=x)
				index++;
			index-=1;
		}
		if (index==-1)
			return 0;
		return index;
	}
	else
		ARM_THROW( ERR_INVALID_ARGUMENT, "ARM_ModelParamsSmiledFRM::IndexOfLastLowerEqInVector : vector of size 0");
};

CC_END_NAMESPACE()

#endif

/*----------------------------------------------------------------------------*/
/*---- End of file ----*/
