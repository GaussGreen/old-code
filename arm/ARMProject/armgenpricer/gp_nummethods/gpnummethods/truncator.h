/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 *	\file truncator.h
 *
 *  \brief
 *
 *	\author  JM Prie, E. Benhamou
 *	\version 1.0
 *	\date November 2004
 */


#ifndef _INGPNUMMETHODS_TRUNCATOR_H
#define _INGPNUMMETHODS_TRUNCATOR_H

#include "gpbase/port.h"
#include "gpbase/rootobject.h"
#include "gpbase/typedef.h"
#include "gpbase/gpvector.h"
#include "gpbase/gpmatrix.h"

#include "typedef.h"

CC_BEGIN_NAMESPACE( ARM )

/// forward declaration
struct ARM_SamplerBase;
struct ARM_ReconnectorBase;
struct ARM_NodeBase;
struct ARM_SliceBase;
class ARM_Truncator1D;
class ARM_TruncatorND;
class ARM_Truncator1DArrowDebreu;
class ARM_TruncatorNDArrowDebreu;
class ARM_TreeIndex;
class ARM_PricingModel;

struct ARM_TruncatorBase : public ARM_RootObject
{
    static const int NO_MAXMAX_INDEX;
    static const size_t MAX_NBSTEP_EXPANSION [];

    enum TruncatorType
	{
		StandardDeviation = 0,
		ArrowDebreu
    };

    ARM_TruncatorBase() : itsReconnector(NULL) {}

    ARM_TruncatorBase( const ARM_TruncatorBase& rhs );
    ARM_TruncatorBase& operator=(const ARM_TruncatorBase& rhs );
    virtual ~ARM_TruncatorBase();

    /// downcast to the appriate type
    virtual ARM_Truncator1D* ToTruncator1D();
    virtual const ARM_Truncator1D* ToTruncator1D() const;

    virtual ARM_TruncatorND* ToTruncatorND();
    virtual const ARM_TruncatorND* ToTruncatorND() const;

    virtual ARM_Truncator1DArrowDebreu* ToTruncator1DArrowDebreu();
    virtual const ARM_Truncator1DArrowDebreu* ToTruncator1DArrowDebreu() const;

    virtual ARM_TruncatorNDArrowDebreu* ToTruncatorNDArrowDebreu();
    virtual const ARM_TruncatorNDArrowDebreu* ToTruncatorNDArrowDebreu() const;

	/// for hybrid creation of the corresponding Truncator1D
    virtual ARM_Truncator1D* CorrespondingTruncator1D(size_t dim) const;

    /// accessor
    const ARM_ReconnectorBase* GetReconnector() const { return itsReconnector; }
    void SetReconnector(const ARM_ReconnectorBase* reconnector);

    /// Init
    virtual void Init(const ARM_SliceVectorPtr& slices, const ARM_SamplerBase* sampler, const ARM_IntMatrixPtr& center=ARM_IntMatrixPtr(NULL)) = 0;

private:
    ARM_ReconnectorBase* itsReconnector;
};



////////////////////////////////////////////
/// ARM_Truncator1D
/// 1D version of a standard deviation based truncator
/////////////////////////////////////////////
class ARM_Truncator1D : public ARM_TruncatorBase
{
private:
    double itsStdDevRatio;

    /// Node index to detect truncation on down hedge
    ARM_IntVector itsMinIndex;

    /// Node index to detect truncation on up hedge
    ARM_IntVector itsMaxIndex;

protected:
    ARM_IntMatrixPtr ComputeCenter(const ARM_SamplerBase* sampler, const ARM_SliceVectorPtr& slices) const;

public:
    ARM_Truncator1D( double stdDevRatio=5.0 )
    : ARM_TruncatorBase(), itsStdDevRatio( stdDevRatio ), itsMinIndex(), itsMaxIndex() {}

    ARM_Truncator1D( const ARM_Truncator1D& rhs );
    ARM_Truncator1D& operator=(const ARM_Truncator1D& rhs );
    virtual ~ARM_Truncator1D();

    virtual ARM_Truncator1D* ToTruncator1D() { return this; }
    virtual const ARM_Truncator1D* ToTruncator1D() const { return this; }

    /// Accessors
    double GetStdDevRatio() const { return itsStdDevRatio; }
    double GetMinIndex(size_t i) const { return itsMinIndex[i]; }
    double GetMaxIndex(size_t i) const { return itsMaxIndex[i]; }

    virtual int GetMaxMaxIndex() const { return ARM_TruncatorBase::NO_MAXMAX_INDEX; }

    /// Init
    virtual void Init(const ARM_SliceVectorPtr& slices, const ARM_SamplerBase* sampler, const ARM_IntMatrixPtr& center=ARM_IntMatrixPtr(NULL));

    /// Truncation detection
    virtual bool NeedTruncation(double sliceTime, ARM_SliceBase* tmpSlice, size_t stateIdx, ARM_SliceBase* tmpNextSlice, int nextNodeIdx) const;

    /// Reconnection is truncation occurs
    virtual void UpdateNode( double sliceTime, ARM_SliceBase* tmpSlice, size_t stateIdx, ARM_SliceBase* tmpNextSlice, double meanRel, double varRel ) const;

    /// ARM_Object
    virtual string toString(const string& indent="", const string& nextIndent="") const;
    virtual ARM_Object* Clone() const { return new ARM_Truncator1D(*this); }
};


////////////////////////////////////////////
/// ARM_Truncator1DArrowDebreu
/// 1D version of a truncation based on standard
/// deviation but also on Arrow-Debreu prices 
/// beyond a time horizon: 
/////////////////////////////////////////////
class ARM_Truncator1DArrowDebreu : public ARM_Truncator1D
{
private:
    size_t itsMaxMaxIndex;
    double itsArrowDebreuThreshold;
    double itsTimeThreshold;

    ARM_IntMatrixPtr itsCenter;

public:
    ARM_Truncator1DArrowDebreu( double stdDevRatio=5.0, size_t maxMaxIndex=25, 
        double arrowDebreuThreshold=1.0e-7, double timeThreshold=10.0 )
    :   ARM_Truncator1D( stdDevRatio ), itsMaxMaxIndex( maxMaxIndex ),
        itsArrowDebreuThreshold( arrowDebreuThreshold ), itsTimeThreshold( timeThreshold ), itsCenter(NULL) {}

    ARM_Truncator1DArrowDebreu( const ARM_Truncator1DArrowDebreu& rhs );
    ARM_Truncator1DArrowDebreu& operator=(const ARM_Truncator1DArrowDebreu& rhs );
    virtual ~ARM_Truncator1DArrowDebreu();

    virtual ARM_Truncator1DArrowDebreu* ToTruncator1DArrowDebreu() { return this; }
    virtual const ARM_Truncator1DArrowDebreu* ToTruncator1DArrowDebreu() const { return this; }

    virtual int GetMaxMaxIndex() const { return itsMaxMaxIndex; }

    /// Init
    virtual void Init(const ARM_SliceVectorPtr& slices, const ARM_SamplerBase* sampler, const ARM_IntMatrixPtr& center=ARM_IntMatrixPtr(NULL));

    /// Truncation detection
    virtual bool NeedTruncation(double sliceTime, ARM_SliceBase* tmpSlice, size_t stateIdx, ARM_SliceBase* tmpNextSlice, int nextNodeIdx) const;

    /// Reconnection is truncation occurs
    virtual void UpdateNode( double sliceTime, ARM_SliceBase* tmpSlice, size_t stateIdx, ARM_SliceBase* tmpNextSlice, double meanRel, double varRel ) const;

    /// ARM_Object
    virtual string toString(const string& indent="", const string& nextIndent="") const;
    virtual ARM_Object* Clone() const  { return new ARM_Truncator1DArrowDebreu(*this); }
};



////////////////////////////////////////////
/// ARM_TruncatorND
/// ND version of a standard deviation based truncator
/////////////////////////////////////////////
class ARM_TruncatorND : public ARM_TruncatorBase
{
private:
    std::vector<double> itsStdDevRatio;
    ARM_IntMatrix itsMinIndex;
    ARM_IntMatrix itsMaxIndex;

protected:
    ARM_IntMatrixPtr ComputeCenter(const ARM_SamplerBase* sampler, const ARM_SliceVectorPtr& slices) const;

public:
    ARM_TruncatorND(size_t dim, double stdDev = 5.0)
    :	ARM_TruncatorBase(), itsStdDevRatio(dim,stdDev), itsMinIndex(), itsMaxIndex() {}
    ARM_TruncatorND(const std::vector<double>& stdDevs)
    :	ARM_TruncatorBase(), itsStdDevRatio(stdDevs), itsMinIndex(), itsMaxIndex() {}
    
	ARM_TruncatorND( const ARM_TruncatorND& rhs )
	:	ARM_TruncatorBase(rhs), itsStdDevRatio( rhs.itsStdDevRatio ), 
	itsMinIndex( rhs.itsMinIndex ),itsMaxIndex( rhs.itsMaxIndex )
	{}

    ARM_TruncatorND& operator=(const ARM_TruncatorND& rhs )
	{
		if( this != &rhs )
		{
			ARM_TruncatorBase::operator=(rhs);
			itsStdDevRatio	= rhs.itsStdDevRatio;
			itsMinIndex		= rhs.itsMinIndex;
			itsMaxIndex		= rhs.itsMaxIndex;
		}
		return *this;
	}
    virtual ~ARM_TruncatorND() {};

    virtual ARM_TruncatorND* ToTruncatorND() { return this; }
    virtual const ARM_TruncatorND* ToTruncatorND() const { return this; }

    virtual void GetMaxMaxIndex(ARM_IntVector& maxMaxIndex) const { maxMaxIndex.resize(dim(),ARM_TruncatorBase::NO_MAXMAX_INDEX); }

    /// Accessors
    const std::vector<double>& GetStdDevRatio() const { return itsStdDevRatio; }
    double GetStdDevRatio(size_t i) const { return itsStdDevRatio[i]; }
    double GetMinIndex(size_t i, size_t j) const { return itsMinIndex(i,j); }
    double GetMaxIndex(size_t i, size_t j) const { return itsMaxIndex(i,j); }

    size_t dim() const { return itsStdDevRatio.size(); }

    /// Init
    virtual void Init(const ARM_SliceVectorPtr& slices, const ARM_SamplerBase* sampler, const ARM_IntMatrixPtr& center=ARM_IntMatrixPtr(NULL));

    /// Truncation detection
    virtual ARM_BoolVector NeedTruncation(double sliceTime, ARM_SliceBase* tmpSlice, size_t stateIdx, ARM_SliceBase* tmpNextSlice, const ARM_TreeIndex& nextNodeIdx) const;

    /// Reconnection is truncation occurs
    virtual void UpdateNode( double sliceTime, ARM_SliceBase* tmpSlice, size_t stateIdx, ARM_SliceBase* tmpNextSlice, const std::vector<double>& meanRel, const std::vector<double>& varRel ) const;

    virtual string toString(const string& indent="", const string& nextIndent="") const;
    virtual ARM_Object* Clone() const { return new ARM_TruncatorND(*this); }

	/// hybrid support
    virtual ARM_Truncator1D* CorrespondingTruncator1D(size_t dim) const;
};



////////////////////////////////////////////
/// ARM_TruncatorNDArrowDebreu
/// ND version of a truncation based on standard
/// deviation but also on Arrow-Debreu prices 
/// beyond a time horizon: 
/////////////////////////////////////////////
class ARM_TruncatorNDArrowDebreu : public ARM_TruncatorND
{
private:
    ARM_IntVector itsMaxMaxIndex;
    double itsArrowDebreuThreshold;
    double itsTimeThreshold;

    ARM_IntMatrixPtr itsCenter;

public:
	/// standard constructor
    ARM_TruncatorNDArrowDebreu(const std::vector<double>& stdDevs, const ARM_IntVector& maxMaxIndex,
        double arrowDebreuThreshold=1.0e-7, double timeThreshold=10.0 )
    :   ARM_TruncatorND(stdDevs), itsMaxMaxIndex(maxMaxIndex),
        itsArrowDebreuThreshold( arrowDebreuThreshold ), itsTimeThreshold( timeThreshold ), itsCenter(NULL) {}

	/// constructor with dimension
    ARM_TruncatorNDArrowDebreu(size_t dim, double stdDev = 5.0 , size_t maxMaxIndex = 25,
        double arrowDebreuThreshold=1.0e-7, double timeThreshold=10.0 )
    :   ARM_TruncatorND(dim,stdDev), itsMaxMaxIndex(dim,maxMaxIndex),
        itsArrowDebreuThreshold( arrowDebreuThreshold ), itsTimeThreshold( timeThreshold ), itsCenter(NULL) {}


    ARM_TruncatorNDArrowDebreu( const ARM_TruncatorNDArrowDebreu& rhs );
    ARM_TruncatorNDArrowDebreu& operator=(const ARM_TruncatorNDArrowDebreu& rhs );
    virtual ~ARM_TruncatorNDArrowDebreu();

    virtual ARM_TruncatorNDArrowDebreu* ToTruncatorNDArrowDebreu() { return this; }
    virtual const ARM_TruncatorNDArrowDebreu* ToTruncatorNDArrowDebreu() const { return this; }

    virtual void GetMaxMaxIndex(ARM_IntVector& maxMaxIndex) const { maxMaxIndex=itsMaxMaxIndex; }

    /// Init
    virtual void Init(const ARM_SliceVectorPtr& slices, const ARM_SamplerBase* sampler, const ARM_IntMatrixPtr& center=ARM_IntMatrixPtr(NULL));

    /// Truncation detection
    virtual ARM_BoolVector NeedTruncation(double sliceTime, ARM_SliceBase* tmpSlice, size_t stateIdx, ARM_SliceBase* tmpNextSlice, const ARM_TreeIndex& nextNodeIdx) const;

    /// Reconnection is truncation occurs
    virtual void UpdateNode( double sliceTime, ARM_SliceBase* tmpSlice, size_t stateIdx, ARM_SliceBase* tmpNextSlice, const std::vector<double>& meanRel, const std::vector<double>& varRel ) const;

    /// ARM_Object
    virtual string toString(const string& indent="", const string& nextIndent="") const;
    virtual ARM_Object* Clone() const { return new ARM_TruncatorNDArrowDebreu(*this); }

	/// hybrid support
    virtual ARM_Truncator1D* CorrespondingTruncator1D(size_t dim) const;
};



CC_END_NAMESPACE()

#endif

/*----------------------------------------------------------------------------*/
/*---- End of file ----*/

