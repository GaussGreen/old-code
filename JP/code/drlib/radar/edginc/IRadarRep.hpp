#ifndef IRADAR_REP_H
#define IRADAR_REP_H

#include "edginc/DECLARE.hpp"
#include "edginc/Format.hpp"
#include "edginc/AtomicArray.hpp"
#include "edginc/Object.hpp"
#include "edginc/RadarRepUtil.hpp"
#include "edginc/IDataPartition.hpp"
#include "edginc/IFunctionBasis.hpp"
#include "edginc/IFittingVarTransform.hpp"

#include <vector>
#include <map>
#include <numeric>

DRLIB_BEGIN_NAMESPACE

using namespace std;

// a class that calculates the value of a radar deal on a specific date for a specific data partition element
class RADAR_DLL IRadarRep {
protected:
	virtual double getValue(const FittingArray& x ) const = 0;
public:
	double operator()(const FittingArray& x ) const { 
        return getValue(x);
    }
	virtual ~IRadarRep() {}
};

DECLARE_REF_COUNT(IRadarRep);

class RADAR_DLL RadarRep : public virtual IRadarRep 
{
public:
	RadarRep(const RegressionCoeff& coef,
		  IFunctionBasisSP f,
		  IFittingVarTransformSP transform ) :
			m_coef( coef ),
			m_f( f),
            m_transform( transform )
	{
		if ( coef.size() != f->getNumBasisFunctions())
		  throw string("Number of basis functions '") + Format::toString(f->getNumBasisFunctions()) +
		    "' does not equal number of coefficients '" + Format::toString(coef.size()) ;
	}

	// returns sum( coef[i]*f_i(x) ); i.e., the value of the radar deal for this set of fitting variables.
	double getValue(const FittingArray& x ) const;

    RegressionCoeff getCoef(void) const;
    
private:
	RegressionCoeff m_coef;
	IFunctionBasisSP m_f;
	IFittingVarTransformSP m_transform;    
};

DECLARE_REF_COUNT(RadarRep);

// a class that calculates the value of a radar deal on a given date.
// essentially this is a partition for the date and a vector of Proxies, 
// one for each element of the partition
class RADAR_DLL RadarRepAtTimeT : public virtual IRadarRep {
public:
    RadarRepAtTimeT(
	       IDataPartitionSP partition,
	       const vector<IRadarRepSP>& radar ) : m_partition( partition ), m_radar(radar) {}
  
    // input:   a vector of fitting variables
    // output:  the value of the radar deal  
    double getValue( const FittingArray& x ) const {
        double result = (*(m_radar[ getBucket(x) ]))(x);
        return result;
    }
private:

	// input:	a vector of fitting variables
	// output:	the index (=bucket) of the corresponding partition element
	size_t getBucket( const FittingArray& x ) const {
		size_t bucket = m_partition->classify(x);
		// ASSERT(part < m_partition->getSize());
		return bucket;
	}

	IDataPartitionSP m_partition;
	vector<IRadarRepSP> m_radar;
};
DECLARE_REF_COUNT(RadarRepAtTimeT);

// a class which transforms market observables into the variables of basis functions. 
// one example is the transformation which maps prices into [0,1] by, say, x --> cumnorm(log(x))
// (formerly known as BasisTransformation)
class RADAR_DLL IFittingVariable {
public:
	// input:	a vector a market observables (eventually, may generalize this to some kind of
	//          market environment structure
	// output:	the value of this particular fitting variable
	virtual double getValue( TDate t, const vector<double>& mktObs ) = 0;
	virtual ~IFittingVariable() {}
};

DECLARE_REF_COUNT(IFittingVariable);

// Interface that bridges "market observables" and "fitting variables". Not needed inside Qlib, but can be needed if we go from SimDesk (via files, for example).
class RADAR_DLL IMarketFitter {
public:
	virtual FittingArray operator() (TDate t, const MarketArray& mkt) const = 0;
	virtual ~IMarketFitter() {}
};
DECLARE_REF_COUNT(IMarketFitter);

// Simple implementation of IMarketFitter
class RADAR_DLL MarketFitter : public IMarketFitter {
public:
	MarketFitter(
		const vector<IFittingVariableSP>& fittingVar) : m_fittingVar(fittingVar)
	{}
	virtual FittingArray operator() (TDate t, const MarketArray& mktObs) const
	{
		FittingArray result( m_fittingVar.size() );
		for ( size_t i = 0; i < result.size(); ++i ) {
		  result[i] = m_fittingVar[i] -> getValue( t, mktObs );
		}
		return result;
	}
private:
	vector<IFittingVariableSP> m_fittingVar;
};

DECLARE_REF_COUNT(IMarketFitter);

DRLIB_END_NAMESPACE

#endif  
