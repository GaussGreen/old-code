/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 *	\file correlmatparam.h
 *
 *  \brief 
 *	\author  E Benhamou
 *	\version 1.0
 *	\date January 2004
 */

#ifndef _INGPINFRA_CORRELMATPARAM_H
#define _INGPINFRA_CORRELMATPARAM_H

/// this headers has to come first
/// as it contains pragma for stupid warning
#include "gpbase/removeidentifiedwarning.h"

#include "gpbase/port.h"
#include "gpbase/gpmatrix.h"  /// to access to the destructor
#include "curvemodelparam.h"
#include "typedef.h"

CC_BEGIN_NAMESPACE( ARM )

/// forward declaration
class ARM_DateStrip;


///////////////////////////////////////////////////
/// \class ARM_CorrelMatParam
/// \brief class for a multi correlations model param
/////////////////////////////////////////////////////
class ARM_CorrelMatParam : public ARM_CurveModelParam
{
public:
	/// constructor, copy assignment and destructor
    ARM_CorrelMatParam( 
		ParamType type,
		std::vector<double>* breakPointTimes,
		std::vector<double>* values,
		const string& interpolatorName,
		std::vector<double>* lowerBound   = NULL,
		std::vector<double>* upperBound   = NULL,
		bool adviseBreakPointTimes	= false );
	ARM_CorrelMatParam( const ARM_CorrelMatParam& rhs );
	ARM_CorrelMatParam& operator=( const ARM_CorrelMatParam& rhs );
    virtual ~ARM_CorrelMatParam();
	
	inline const ARM_MultiCurve* GetMultiCurve() const { return itsCurves; }
	inline ARM_MultiCurve* GetMultiCurve() { return itsCurves; }
	void SetMultiCurve(ARM_MultiCurve* mCurve);

	virtual void SetFactorCount( size_t factorCount ){};
	void UpdateRealizedCorrel();

	std::vector<double>* GetData( ARM_DataType type, long& rows, long& cols ) const;

	virtual ARM_Object* Clone() const { return new ARM_CorrelMatParam(*this); }
    virtual string toString(const string& indent="", const string& nextIndent="") const;

private:
    ARM_MultiCurve* itsCurves;
	ARM_GP_Matrix* itsRealizedCorrelMatrix;
	void CreateMultiCurves( const string& interpolatorName );
};


///////////////////////////////////////////////////
/// \class ARM_CorrelMatParam
/// \brief correlation matrix that does the ACP
/////////////////////////////////////////////////////

class ARM_CorrelACPMatParam : public ARM_CorrelMatParam
{
public:
	ARM_CorrelACPMatParam( 
		std::vector<double>* breakPointTimes,
		std::vector<double>* values,
		const string& interpolatorName,
		std::vector<double>* lowerBound   = NULL,
		std::vector<double>* upperBound   = NULL,
		bool adviseBreakPointTimes	= false );

	ARM_CorrelACPMatParam( const ARM_CorrelACPMatParam& rhs );
	ARM_CorrelACPMatParam& operator=(const ARM_CorrelACPMatParam& rhs );
	virtual ~ARM_CorrelACPMatParam();
	virtual ARM_Object* Clone() const { return new ARM_CorrelACPMatParam(*this); }
    virtual string toString(const string& indent="", const string& nextIndent="") const;
	
	void SetInputCorrelMatrix( ARM_GP_Matrix* inputCorrelMatrix );
	void UpdateMultiCurves( std::vector<double>* breakPointTimes );
	virtual void SetFactorCount( size_t factorCount );

private:
    ARM_GP_Matrix* itsInputCorrelMatrix;
    ARM_GP_Matrix* itsACPMatrix;
	int itsFactorCount;
	void PrintMatrix( CC_Ostringstream& os, const string& matrixName, ARM_GP_Matrix* matrix, const string& indent, const string& nextIndent ) const;
};




///////////////////////////////////////////////////
/// \class ARM_CorrelTrigoMatParam
/// \brief
/// ARM_CorrelTrigoMatParam is a trigonometric matrix
///		parametrized only in one coefficient theta
/////////////////////////////////////////////////////
class ARM_CorrelTrigoMatParam : public ARM_CorrelACPMatParam
{
public:
	ARM_CorrelTrigoMatParam( 
		double theta,
		const ARM_DateStrip& datestrip,
		double asOf,
		const string& interpolatorName,
		std::vector<double>* lowerBound   = NULL,
		std::vector<double>* upperBound   = NULL,
		bool adviseBreakPointTimes	= false );

	ARM_CorrelTrigoMatParam(
		double theta,
		const std::vector<double>* resetDates,
		double asOf,
		const string& interpolatorName,
		std::vector<double>* lowerBound,
		std::vector<double>* upperBound,
		bool adviseBreakPointTimes );

	ARM_CorrelTrigoMatParam( const ARM_CorrelTrigoMatParam& rhs );
	ARM_CorrelTrigoMatParam& operator=(const ARM_CorrelTrigoMatParam& rhs );

    virtual void MergeModelParam( ARM_ModelParam* NewValue );
	virtual ~ARM_CorrelTrigoMatParam();
	virtual ARM_Object* Clone() const { return new ARM_CorrelTrigoMatParam(*this); }
    virtual string toString(const string& indent="", const string& nextIndent="") const;
	virtual void SetFactorCount( size_t factorCount );

private:
	double itsTheta;
	std::vector<double>* itsCorrelBreakPointTimes;
};



///////////////////////////////////////////////////
/// \class ARM_CorrelMatNoCalibParam
/// \brief
/// correlation matrix without any calibration
/////////////////////////////////////////////////////
class ARM_CorrelMatNoCalibParam : public ARM_CorrelACPMatParam
{
public:
	ARM_CorrelMatNoCalibParam( 
		std::vector<double>* breakPointTimes,
		std::vector<double>* values,
		const string& interpolatorName,
		std::vector<double>* lowerBound   = NULL,
		std::vector<double>* upperBound   = NULL,
		bool adviseBreakPointTimes	= false );

	ARM_CorrelMatNoCalibParam( const ARM_CorrelMatNoCalibParam& rhs );
	ARM_CorrelMatNoCalibParam& operator=(const ARM_CorrelMatNoCalibParam& rhs );
	~ARM_CorrelMatNoCalibParam() {};
	virtual ARM_Object* Clone() const { return new ARM_CorrelMatNoCalibParam(*this); }
	virtual void SetFactorCount( size_t factorCount );
    virtual void MergeModelParam( ARM_ModelParam* NewValue );
};



///////////////////////////////////////////////////
/// \class ARM_CorrelMatParam
/// \brief
/// ARM_BrownianCorrelationParam is a 2D vector of correlation
///		parametrized with theta(i) (leads to cos(Theta(i)) sin(Theta(i))
/////////////////////////////////////////////////////
class ARM_BrownianCorrelationParam : public ARM_CorrelMatParam
{
public:
	ARM_BrownianCorrelationParam( 
		std::vector<double>* breakPointTimes,
		std::vector<double>* values,
		const string& interpolatorName,
		std::vector<double>* lowerBound   = NULL,
		std::vector<double>* upperBound   = NULL,
		bool adviseBreakPointTimes	= false );
	ARM_BrownianCorrelationParam( const ARM_BrownianCorrelationParam& rhs);
	ARM_BrownianCorrelationParam& operator=( const ARM_BrownianCorrelationParam& rhs );
	virtual ~ARM_BrownianCorrelationParam();
	virtual ARM_Object* Clone() const { return new ARM_BrownianCorrelationParam(*this); }
	virtual void SetFactorCount( size_t factorCount );
    virtual void MergeModelParam( ARM_ModelParam* NewValue );
private:
	void ComputeMultiCurves( std::vector<double>* values );
};



CC_END_NAMESPACE()

#endif

/*----------------------------------------------------------------------------*/
/*---- End of file ----*/

