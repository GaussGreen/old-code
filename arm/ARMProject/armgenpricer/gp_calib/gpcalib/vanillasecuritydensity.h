/*!
 *
 * Copyright (c) CDC IXIS CM November 2003 Paris
 *
 *	\file vanillasecuritydensity.h
 *
 *  \brief vanilla security density
 *
 *	\author  A. Schauly
 *	\version 1.0
 *	\date August 2005
 */

#ifndef _INGPCALIB_VANILLADENSITY_H
#define _INGPCALIB_VANILLADENSITY_H

/// this header has to come first
/// as env.h defines pre-processor constants for checking and validation
/// namely for std::vector<double> and ARM_Matrix
#include "gpbase/env.h"
#include "gpbase/rootobject.h"
#include "gpbase/assignop.h"
#include "gpinfra/pricingfunctionequity.h"
#include "gpinfra/typedef.h"
#include "gpcalib/typedef.h"
#include "densityfunctors.h"

CC_BEGIN_NAMESPACE( ARM )

class ARM_DateStrip;

class ARM_VanillaSecurityDensity : public ARM_RootObject
{
public:
	/// Constructors and Destructors
	ARM_VanillaSecurityDensity( double resetDate = 0,
								double startDate = 0,
								double endDate	 = 0,
								ARM_DensityFunctorPtr densityFunctor = ARM_DensityFunctorPtr(NULL), 
								int frequency	 = GETDEFAULTVALUE,
								int daycount	 = GETDEFAULTVALUE,
								int stubRule	 = GETDEFAULTVALUE,
								double weight	 = 1.,
								double adjfwdadd = 0.,
								double adjfwdmult = 1.,
								ARM_ZeroCurvePtr ZeroCurve = ARM_ZeroCurvePtr(NULL)
								) ;

	ARM_VanillaSecurityDensity(const ARM_VanillaSecurityDensity& rhs);

	ASSIGN_OPERATOR(ARM_VanillaSecurityDensity);
	virtual ~ARM_VanillaSecurityDensity();

	inline double getResetDate() const { return itsResetDate; }
	inline void setResetDate( const double& resetDate ) { itsResetDate = resetDate; }
	
	inline double getStartDate() const { return itsStartDate; }
	inline void setStartDate( const double& startDate ) { itsStartDate = startDate; }

	inline double getEndDate() const { return itsEndDate; }
	inline void setEndDate( const double& endDate ) { itsEndDate = endDate; }
	
	inline ARM_ZeroCurvePtr getZeroCurve() const { return itsZeroCurve; }
	inline void setZeroCurve(const ARM_ZeroCurvePtr& zeroCurve) { itsZeroCurve = zeroCurve; }
	
	inline ARM_DensityFunctorPtr getDensityFunctor() const { return itsDensityFunctor; }
	inline void setDensityFunctor( const ARM_DensityFunctorPtr& densityFunctor ) { itsDensityFunctor=densityFunctor; }

	inline ARM_DateStripPtr getCalibSchedule() const { return itsCalibSchedule; }
	virtual void AssociateCalibSchedule(ARM_DateStripPtr calibSchedule );

	inline double getWeight() const {return itsWeight;};
	inline void setWeight(double weight) {itsWeight = weight;};

	/// some other accessors
	inline const int& getFrequency() const {return itsFrequency;}
	inline const std::vector<double>& getInterestTerms() const {return itsInterestTerms;}
	inline const ARM_IntVector& getPayDatesRelIndexes() const {return itsPayDatesRelIndexes;}
	inline const ARM_IntVector& getPayDatesAbsIndexes() const {return itsPayDatesAbsIndexes;}

	/// Standard ARM Object Support
	virtual ARM_Object* Clone() const {return new ARM_VanillaSecurityDensity(*this);};
	virtual string ExportShortName() const { return "LVSEC";}
	virtual string toString(const string& indent="",const string& nextIndent="") const ;

	/// What it is for
	virtual void UpdateStates( const ARM_PricingStatesPtr& states, const ARM_GP_VectorPtr& probabilities, size_t lineIdx ) const ;
		
	/// forward and level (libor or swaprate)
	virtual double ComputeForwardRate() const {return ComputeForwardRate(itsStartDate, itsEndDate, itsPayDatesAbsIndexes, itsInterestTerms);};
	virtual double ComputeLevel() const {return ComputeLevel(itsStartDate, itsEndDate, itsPayDatesAbsIndexes, itsInterestTerms);};


private:

	
protected: 
	ARM_DensityFunctorPtr itsDensityFunctor; /// Density of itsUnderlying
	ARM_ZeroCurvePtr itsZeroCurve; /// ZeroCurve to build its ForwardValues

	ARM_DateStripPtr itsCalibSchedule; /// ASSOC (owned by numerical model fitter)
	
	/// underlying descr (libor or swaprate)
	/// required to be included within itsCalibSchedule
	/// !!! these are dates, not times.... !!!
	double itsResetDate;
	double itsStartDate;
	double itsEndDate;

	int itsFrequency;
	int itsDayCount;
	int itsStubRule;

	double itsWeight;
	double itsAdjFwdAdd;
	double itsAdjFwdMult;

	std::vector<double> itsInterestTerms;
	ARM_IntVector itsPayDatesRelIndexes; /// w.r.t start date
	ARM_IntVector itsPayDatesAbsIndexes; /// w.r.t schedule

protected:

	double	ComputeForwardRate(double startdate, double enddate, const ARM_IntVector& PayDatesAbsIndexes, const std::vector<double>& InterestTerms) const;
	double	ComputeLevel(double startdate, double enddate, const ARM_IntVector& PayDatesAbsIndexes, const std::vector<double>& InterestTerms) const;

	/// compute indexes
	virtual void IndexesInCalibSchedule(double StartDate, double EndDate, size_t& resetAnStartIndex, size_t& endIndex) const;
	virtual void IndexesInCalibSchedule(size_t& resetAnStartIndex, size_t& endIndex) const {IndexesInCalibSchedule(itsStartDate, itsEndDate, resetAnStartIndex, endIndex);};


	virtual void AssociateCalibSchedule(ARM_DateStripPtr calibSchedule, double StartDate, double EndDate, int Frequency, int Daycount, std::vector<double>& InterestTerms, 
					ARM_IntVector& PayDatesRelIndexes, ARM_IntVector& PayDatesAbsIndexes);
};

class ARM_VanillaSecurityDensitySpread : public ARM_VanillaSecurityDensity
{
private:
	double			itsStartDate2;
	double			itsEndDate2;
	int				itsFrequency2;
	int				itsDayCount2;

	std::vector<double>	itsInterestTerms2;
	ARM_IntVector	itsPayDatesRelIndexes2; /// w.r.t start date
	ARM_IntVector	itsPayDatesAbsIndexes2; /// w.r.t schedule

public:

	ARM_VanillaSecurityDensitySpread( double resetDate = 0,
								double startDate = 0,
								double endDate	 = 0,
								double startDate2 = 0,
								double endDate2 = 0,
								ARM_DensityFunctorPtr densityFunctor = ARM_DensityFunctorPtr(NULL), 
								int frequency	 = GETDEFAULTVALUE,
								int daycount	 = GETDEFAULTVALUE,
								int frequency2 = GETDEFAULTVALUE,
								int daycount2 = GETDEFAULTVALUE,
								int stubRule	 = GETDEFAULTVALUE,
								ARM_ZeroCurvePtr ZeroCurve = ARM_ZeroCurvePtr(NULL) ) ;

	ARM_VanillaSecurityDensitySpread(const ARM_VanillaSecurityDensitySpread& rhs);

	ASSIGN_OPERATOR(ARM_VanillaSecurityDensitySpread);
	virtual ~ARM_VanillaSecurityDensitySpread();
	
	inline double getStartDate1() const { return itsStartDate; }
	inline void setStartDate1( const double& startDate ) { itsStartDate = startDate; }

	inline double getEndDate1() const { return itsEndDate; }
	inline void setEndDate1( const double& endDate ) { itsEndDate = endDate; }
	
	inline double getStartDate2() const { return itsStartDate2; }
	inline void setStartDate2( const double& startDate ) { itsStartDate2 = startDate; }

	inline double getEndDate2() const { return itsEndDate2; }
	inline void setEndDate2( const double& endDate ) { itsEndDate2 = endDate; }

	/// some other accessors
	inline const int& getFrequency1() const {return itsFrequency;}
	inline const std::vector<double>& getInterestTerms1() const {return itsInterestTerms;}
	inline const ARM_IntVector& getPayDatesRelIndexes1() const {return itsPayDatesRelIndexes;}
	inline const ARM_IntVector& getPayDatesAbsIndexes1() const {return itsPayDatesAbsIndexes;}

	inline const int& getFrequency2() const {return itsFrequency2;}
	inline const std::vector<double>& getInterestTerms2() const {return itsInterestTerms2;}
	inline const ARM_IntVector& getPayDatesRelIndexes2() const {return itsPayDatesRelIndexes2;}
	inline const ARM_IntVector& getPayDatesAbsIndexes2() const {return itsPayDatesAbsIndexes2;}

	virtual void AssociateCalibSchedule(ARM_DateStripPtr calibSchedule );

	/// Standard ARM Object Support
	virtual ARM_Object* Clone() const {return new ARM_VanillaSecurityDensitySpread(*this);};
	virtual string toString(const string& indent="",const string& nextIndent="") const ;

	/// What it is for	
	/// forward and level (libor or swaprate)
	virtual double ComputeForwardRate1() const {return ComputeForwardRate(itsStartDate, itsEndDate, itsPayDatesAbsIndexes, itsInterestTerms);};
	virtual double ComputeLevel1() const {return ComputeLevel(itsStartDate, itsEndDate, itsPayDatesAbsIndexes, itsInterestTerms);};

	virtual double ComputeForwardRate2() const {return ComputeForwardRate(itsStartDate2, itsEndDate2, itsPayDatesAbsIndexes2, itsInterestTerms2);};
	virtual double ComputeLevel2() const {return ComputeLevel(itsStartDate2, itsEndDate2, itsPayDatesAbsIndexes2, itsInterestTerms2);};
};

class ARM_VanillaSecurityDensityFX : public ARM_VanillaSecurityDensity
{
	public:
	/// Constructors and Destructors
	ARM_VanillaSecurityDensityFX( double resetDate = 0,
								ARM_DensityFunctorPtr densityFunctor = ARM_DensityFunctorPtr(NULL),
								ARM_ZeroCurvePtr DomZeroCurve = ARM_ZeroCurvePtr(NULL),
								ARM_ZeroCurvePtr FgnZeroCurve = ARM_ZeroCurvePtr(NULL),
								double FXSpot = 1.0);
	ARM_VanillaSecurityDensityFX(const ARM_VanillaSecurityDensityFX& rhs);

	ASSIGN_OPERATOR(ARM_VanillaSecurityDensityFX);
	virtual ~ARM_VanillaSecurityDensityFX();

	inline ARM_ZeroCurvePtr getFgnZeroCurve() const { return itsFgnZeroCurve; }
	inline void setFgnZeroCurve(const ARM_ZeroCurvePtr& fgnZeroCurve) { itsFgnZeroCurve = fgnZeroCurve;}
	inline double getFXSpot() const { return itsFXSpot; }
	inline void setFXSpot(double FXSpot) { itsFXSpot = FXSpot; }
	inline double getSettlementDate() const { return itsSettlementDate; }
	inline void setSettlementDate(double settlementDate) { itsSettlementDate = settlementDate; }

	// We don't need calib schedule for FX
	virtual void AssociateCalibSchedule(ARM_DateStripPtr calibSchedule ){};

	/// Standard ARM Object Support
	virtual ARM_Object* Clone() const {return new ARM_VanillaSecurityDensityFX(*this);};
	virtual string ExportShortName() const { return "LFSEC";}
	virtual string toString(const string& indent="",const string& nextIndent="") const;

	/// forward FX Rate
	virtual double ComputeForwardRate() const;

private:
	double itsFXSpot;
	double itsSettlementDate;
	ARM_ZeroCurvePtr itsFgnZeroCurve;
};


CC_END_NAMESPACE()
#endif