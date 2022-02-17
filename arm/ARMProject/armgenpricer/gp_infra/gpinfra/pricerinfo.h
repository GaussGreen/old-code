/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 *	\file pricerinfo.h
 *
 *  \brief 
 *
 *	\author  E. Benhamou
 *	\version 1.0
 *	\date February 2004
 */


#ifndef _INGPINFRA_PRICERINFO_H
#define _INGPINFRA_PRICERINFO_H

/// this header comes firts as it includes some preprocessor constants!
#include "gpbase/removeidentifiedwarning.h"

/// env.h defines pre-processor constants for checking and validation
#include "gpbase/env.h"
#include "gpbase/port.h"
#include "gpbase/rootobject.h"
#include "gpbase/gpmatrix.h"

#include "gramfunctorargdict.h"
#include "pricingadviser.h"
#include <map>
CC_USING_NS( std, map )
CC_USING_NS( std, less )
#include <vector>
CC_USING_NS(std,vector)
#include <gpbase/gpmatrix.h>
#include <gpbase/timer.h>

#include "gpnumlib/stddevfunc.h"


CC_BEGIN_NAMESPACE( ARM )

/// forward declaration
class ARM_PricingModel;
class ARM_PricingStates;
class ARM_MCMethod;

/// pricer info is a pure interface class to store
/// pricing information ...
/// this allows the GenPricer object to be used as a single
//
// however, it needs to use clone. We took this opportunity to define the clone
/// a const method!

struct ARM_PricerInfo : public ARM_RootObject, public ARM_Timer
{
	ARM_PricerInfo();
	ARM_PricerInfo( const ARM_PricerInfo& rhs );
	ARM_PricerInfo& operator=(const ARM_PricerInfo& rhs );
	virtual ~ARM_PricerInfo();

    /// ================== Standard ARM object support ========================
    virtual string toString(const string& indent="", const string& nextIndent="") const;

	virtual void StoreCurrentLoopPrice( 
		const ARM_PricingStatesPtr& states, 
		size_t bucketSize, size_t currentBucketIndex, 
		ARM_PricingAdviser::PricingDir pricingDir, int colIdx=0  )			= 0;
	virtual void PostProcess(const ARM_TimeInfoPtrVector& timeInfo,const ARM_NumMethodPtr& numMethod) = 0;
	virtual void reset()                                        = 0;
	virtual string GetName() const								= 0;
	virtual double GetPrice() const								= 0;
	virtual ARM_GramFctorArgDict GetContents(const string& colName="") const			= 0;
};

/// class for single loop pricing .... like trees!
struct ARM_PInfo_SingleLoop : public ARM_PricerInfo
{
	ARM_PInfo_SingleLoop();
	ARM_PInfo_SingleLoop( const ARM_PInfo_SingleLoop& rhs );
	ARM_PInfo_SingleLoop& operator=( const ARM_PInfo_SingleLoop& rhs );
	virtual ~ARM_PInfo_SingleLoop();

    /// ================== Standard ARM object support ========================
	virtual ARM_Object* Clone() const;

	virtual void StoreCurrentLoopPrice( 
		const ARM_PricingStatesPtr& states, 
		size_t bucketSize, size_t currentBucketIndex, 
		ARM_PricingAdviser::PricingDir pricingDir, int colIdx=0 );
	virtual void PostProcess(const ARM_TimeInfoPtrVector& timeInfo,const ARM_NumMethodPtr& numMethod);
    virtual void reset() {};
	virtual string GetName() const { return "Single Loop"; }
	virtual double GetPrice() const;
	virtual ARM_GramFctorArgDict GetContents(const string& colName="") const;

private:
	double itsPrice;
	ARM_VectorPtr       itsIntermediatePrices;
    ARM_GP_MatrixPtr    itsPriceSnapshots;
    ARM_GP_MatrixPtr    itsSpotProbabilities;
};

/// class for single loop pricing .... like PDE!
struct ARM_PInfo_SingleLoopPDE : public ARM_PInfo_SingleLoop
{
	ARM_PInfo_SingleLoopPDE();
	ARM_PInfo_SingleLoopPDE( size_t priceIndex ) : itsPriceIndex( priceIndex ),ARM_PInfo_SingleLoop() {}
	ARM_PInfo_SingleLoopPDE( const ARM_PInfo_SingleLoopPDE& rhs );
	ARM_PInfo_SingleLoopPDE& operator=( const ARM_PInfo_SingleLoopPDE& rhs );
	virtual ~ARM_PInfo_SingleLoopPDE();

    /// ================== Standard ARM object support ========================
	virtual ARM_Object* Clone() const;

	virtual void StoreCurrentLoopPrice( 
		const ARM_PricingStatesPtr& states, 
		size_t bucketSize, size_t currentBucketIndex, 
		ARM_PricingAdviser::PricingDir pricingDir, int colIdx=0 );
	virtual void PostProcess(const ARM_TimeInfoPtrVector& timeInfo,const ARM_NumMethodPtr& numMethod);
    virtual void reset() {};
	virtual string GetName() const { return "Single Loop For PDE"; }
	virtual double GetPrice() const;
	virtual ARM_GramFctorArgDict GetContents(const string& colName="") const;

private:
	double itsPrice;
	size_t itsPriceIndex;
	ARM_VectorPtr       itsIntermediatePrices;
    ARM_GP_MatrixPtr    itsPriceSnapshots;
    ARM_GP_MatrixPtr    itsSpotProbabilities;
};

/// multiple loop pricer info This allows to use pricing that requires multiple loop
/// like in Monte Carlo for instance!
struct ARM_PInfo_MultipleLoop : public ARM_PricerInfo
{
	virtual ~ARM_PInfo_MultipleLoop();

    /// ================== Standard ARM object support ========================
	virtual ARM_Object* Clone() const;

	virtual void StoreCurrentLoopPrice( 
		const ARM_PricingStatesPtr& states, 
		size_t bucketSize, size_t currentBucketIndex, 
		ARM_PricingAdviser::PricingDir pricingDir, int colIdx=0);
	virtual void PostProcess(const ARM_TimeInfoPtrVector& timeInfo,const ARM_NumMethodPtr& numMethod);
	virtual void reset();
	virtual string GetName() const { return "Multiple Loop"; }
	virtual double GetPrice() const;
	virtual ARM_GramFctorArgDict GetContents(const string& colName="") const;
	inline ARM_MomentFuncPtr GetMomentFunc() const { return itsMomentFunc; };
	inline int GetRunningAvgStep() const { return itsRunningAvgStep; }

private:
	double itsPrice;
	ARM_VectorPtr itsIntermediatePrices;
	double itsStdDev;
	ARM_VectorPtr itsIntermediateStdDevs;
	ARM_VectorPtr itsValues;
	ARM_GP_MatrixPtr itsIntermediateValues;
	ARM_VectorPtr itsRunningAverage;
	int itsRunningAvgStep;
	ARM_GP_MatrixPtr itsIntermediateRunningAverages;
	double itsEmpiricalStdDev;
	ARM_VectorPtr itsIntermediateEmpiricalStdDevs;
	ARM_VectorPtr itsBuckets;
	ARM_VectorPtrVectorPtr itsValuesPerBucket;
	ARM_VectorPtr itsPricesPerBucket;
	ARM_VectorPtr itsStdDevPerBucket;
	size_t itsValuesIndex;
	ARM_MomentFuncPtr itsMomentFunc;
	ARM_MomentFuncPtrVector itsInterMomentFunc;
	bool itsIntermediateArrayInitialized;

	/// to prevent creation outside the factory method!
	ARM_PInfo_MultipleLoop(const ARM_PricingModel& mod, const ARM_MomentFuncPtr& momentFunc = ARM_MomentFuncPtr( new ARM_StdMomentFunc ) );
	ARM_PInfo_MultipleLoop( const ARM_PInfo_MultipleLoop& rhs );
	ARM_PInfo_MultipleLoop& operator=( const ARM_PInfo_MultipleLoop& rhs );
	friend class ARM_MCMethod; /// make MC method friend to allow only this method to create it
};


/// forward declaration
class ARM_PInfo_Map;


/// pricer info that stores the control variate 
struct ARM_PInfo_Cov : public ARM_PricerInfo
{
private:
	double itsPrice;
	double itsStdDev;
	double itsTheoreticalStdDev;
	ARM_VectorPtr itsValues;
	ARM_VectorPtr itsRunningAverage;
	int	itsRunningAvgStep;

	ARM_VectorPtr itsBetas;
	ARM_GP_MatrixPtr itsVarCovCV;
	ARM_GP_VectorPtr itsCovRefCV;
	ARM_GP_VectorPtr itsCorrRefCV;

	string itsRefColumn;
	ARM_StringVectorPtr itsCVColumns;
	ARM_GP_VectorPtr itsCVPrices;

	/// various function for the constructor
	void ComputeBetas( ARM_PInfo_Map* pinfoMap, ARM_PInfo_MultipleLoop* pinfoMultiLoop, const std::vector<double>& userBetas );
	void PostProcessCV( ARM_PInfo_Map* pinfoMap, ARM_PInfo_MultipleLoop* pinfoMultiLoop, const std::vector<double>& userBetas );

public:
	ARM_PInfo_Cov( ARM_PInfo_Map* pinfoMap, const string& refPriceColumn, const ARM_StringVector& CVColumns, const std::vector<double>& CVPrices, const std::vector<double>& Betas );
	ARM_PInfo_Cov( const ARM_PInfo_Cov & rhs );
	ARM_PInfo_Cov & operator=(const ARM_PInfo_Cov & rhs );
	virtual ~ARM_PInfo_Cov(){};

	virtual void StoreCurrentLoopPrice( const ARM_PricingStatesPtr& states, 
		size_t bucketSize, size_t currentBucketIndex, ARM_PricingAdviser::PricingDir pricingDir, int colIdx=0 );
	virtual void PostProcess(const ARM_TimeInfoPtrVector& timeInfo,const ARM_NumMethodPtr& numMethod);
	virtual ARM_GramFctorArgDict GetContents(const string& colName="") const;
	virtual double GetPrice() const { return itsPrice; }
    virtual void reset() {};
	virtual ARM_Object* Clone() const { return new ARM_PInfo_Cov(*this); }
	virtual string GetName() const { return "Covariance pricer info"; }
	inline ARM_VectorPtr GetBetas() {return itsBetas;}
};


/// The pricer info will be used to values few columns of the same secuirty
class ARM_PInfo_Map : public ARM_PricerInfo
{
public:
	ARM_PInfo_Map(const vector<string>& columnNames, const vector<ARM_PricerInfo*>& pricerInfo);
	virtual ~ARM_PInfo_Map() {};

    /// ================== Standard ARM object support ========================
	virtual ARM_Object* Clone() const;
	virtual string toString(const string& indent="", const string& nextIndent="") const;

	virtual void StoreCurrentLoopPrice(
	const ARM_PricingStatesPtr& states, 
	size_t bucketSize, size_t currentBucketIndex, 
	ARM_PricingAdviser::PricingDir pricingDir, int colIdx=0);
	virtual void PostProcess(const ARM_TimeInfoPtrVector& timeInfo,const ARM_NumMethodPtr& numMethod);
	virtual void reset();
	virtual string GetName() const { return "Map"; }
	virtual double GetPrice() const;
	virtual ARM_GramFctorArgDict GetContents(const string& colName="") const;
	inline pair<string, ARM_PricerInfoPtr> operator[]( size_t i ) const { return itsPricerInfoMap[i]; }
	inline pair<string, ARM_PricerInfoPtr>& operator[]( size_t i ) { return itsPricerInfoMap[i]; }
	inline void push_back( const string& name, const ARM_PricerInfoPtr& pinfo ) { itsPricerInfoMap.push_back( pair<string, ARM_PricerInfoPtr>( name, pinfo ) ); }

private:
	typedef vector< pair<string, ARM_PricerInfoPtr> > PricerInfoMap;
	PricerInfoMap itsPricerInfoMap;
	ARM_PInfo_Map( const ARM_PInfo_Map& rhs );
	ARM_PInfo_Map& operator=( const ARM_PInfo_Map& rhs );
};

CC_END_NAMESPACE()

#endif

/*----------------------------------------------------------------------------*/
/*---- End of file ----*/
