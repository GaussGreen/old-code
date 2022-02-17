/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 * \file gramfunctorarg.h
 *
 *  \brief 
 *
 *	\author  E. Benhamou
 *	\version 1.0
 *	\date October 2003
 */



#ifndef _INGPINFRA_GRAMFUNCTORARG_H
#define _INGPINFRA_GRAMFUNCTORARG_H

/// this header comes first as it includes some preprocessor constants!
#include "gpbase/removeidentifiedwarning.h"
#include "gpbase/port.h"
#include "typedef.h"
#include "gramargs.h"
#include <glob/dates.h>
#include "gpbase/curve.h"
#include "gpbase/datestrip.h"
#include <string>
CC_USING_NS( std, string )

CC_BEGIN_NAMESPACE( ARM )


/////////////////////////////////////////////////////////////////
/// \class	ARM_GramFctorArg
/// \brief	class to describe general argument of the grammar
/// can be double, date, string or vector of double
/////////////////////////////////////////////////////////////////

struct ARM_GramFctorArg
{
	/// for shorter name!
	typedef GramFuncArgType Type; 

	/// helper function for types
	static string TypeToString( Type t) { return GramFuncArgTypeToString( t ); }
	static Type GetResultType( Type );

	/// Accessors Get and Set
	void SetDouble( double d )											{ itsDouble=d;          itsType=GFAT_DOUBLE_TYPE;			}
	void SetString( const string& s )									{ itsString=s;          itsType=GFAT_STRING_TYPE;			}
	void SetDate( const ARM_Date& d )									{ itsDate=d;            itsType=GFAT_DATE_TYPE;				}
	void SetVector( const ARM_VectorPtr& v )							{ itsVector=v;          itsType=GFAT_VECTOR_TYPE;			}
    void SetMatrix( const ARM_GP_MatrixPtr& m )							{ itsMatrix=m;          itsType=GFAT_MATRIX_TYPE;			}
    void SetCurve( const ARM_GP_CurvePtr& c )							{ itsCurve=c;           itsType=GFAT_CURVE_TYPE;			}
    void SetStringVector( const ARM_StringVectorPtr& sv)				{ itsStringVector=sv;   itsType=GFAT_STRINGVEC_TYPE;		}
    void SetStringVectorTrans( const ARM_StringVectorPtr& sv)			{ itsStringVector=sv;   itsType=GFAT_STRINGVECTRANS_TYPE;	}
  	void SetDateStrip( const ARM_DateStripPtr& sched)					{ itsDateStrip=sched;   itsType=GFAT_DATESTRIP_TYPE;		}

    inline Type GetType() const											{ return itsType; }
	inline double GetDouble() const										{ CheckType(GFAT_DOUBLE_TYPE);			return itsDouble;		}
	inline string GetString() const										{ CheckType(GFAT_STRING_TYPE);			return itsString;		}
	inline ARM_Date GetDate() const										{ CheckType(GFAT_DATE_TYPE);			return itsDate;			}
	inline ARM_VectorPtr GetVector() const								{ CheckType(GFAT_VECTOR_TYPE);			return itsVector;		}
	inline ARM_VectorPtr GetVector()									{ CheckType(GFAT_VECTOR_TYPE);			return itsVector;		}
	inline ARM_GP_MatrixPtr GetMatrix() const							{ CheckType(GFAT_MATRIX_TYPE);			return itsMatrix;		}
	inline ARM_GP_MatrixPtr GetMatrix()									{ CheckType(GFAT_MATRIX_TYPE);			return itsMatrix;		}
    inline ARM_GP_CurvePtr GetCurve()									{ CheckType(GFAT_CURVE_TYPE);			return itsCurve;		}
    inline ARM_GP_CurvePtr GetCurve() const								{ CheckType(GFAT_CURVE_TYPE);			return itsCurve;		}
    inline ARM_StringVectorPtr GetStringVector()						{ CheckType(GFAT_STRINGVEC_TYPE);		return itsStringVector; }
    inline ARM_StringVectorPtr GetStringVector() const					{ CheckType(GFAT_STRINGVEC_TYPE);		return itsStringVector; }
    inline ARM_StringVectorPtr GetStringVectorTrans()					{ CheckType(GFAT_STRINGVECTRANS_TYPE);  return itsStringVector; }
    inline ARM_StringVectorPtr GetStringVectorTrans() const				{ CheckType(GFAT_STRINGVECTRANS_TYPE);  return itsStringVector; }
	inline ARM_DateStripPtr GetDateStrip()								{ CheckType(GFAT_DATESTRIP_TYPE);		return itsDateStrip;	}
    inline ARM_DateStripPtr GetDateStrip() const						{ CheckType(GFAT_DATESTRIP_TYPE);		return itsDateStrip;	}

	/// constructor with arguments
	/// because of the default d=0.0 a default ARM_GramFctorArg
	/// is a double with value 0
	ARM_GramFctorArg( double d = 0.0);
	ARM_GramFctorArg( const string& s);
	ARM_GramFctorArg( const ARM_Date& d );
	ARM_GramFctorArg( const ARM_VectorPtr& v );
	ARM_GramFctorArg( const ARM_GP_MatrixPtr& m );
    ARM_GramFctorArg( const ARM_GP_CurvePtr& c );
    ARM_GramFctorArg( const ARM_StringVectorPtr& sv );
    ARM_GramFctorArg( const ARM_DateStripPtr& sched );


	/// for debugging and exception and so on
	string toString() const;

	/// memory management
	void Reset();

	/// copy constructor, assignment and destructor
	ARM_GramFctorArg( const ARM_GramFctorArg& rhs );
	ARM_GramFctorArg& operator=( const ARM_GramFctorArg& rhs );
	~ARM_GramFctorArg();
	ARM_GramFctorArg Duplicate() const;
	
private:

	/// instead of using union that would need to 
	/// give standart type (basically with no copy
	/// constructor), we use for each data type
	/// a member to be able to use smartPointors!
	double              itsDouble;
	ARM_VectorPtr       itsVector;
	ARM_GP_MatrixPtr    itsMatrix;
    ARM_GP_CurvePtr     itsCurve;
	ARM_StringVectorPtr itsStringVector;
	string              itsString;
	ARM_Date            itsDate;
    ARM_DateStripPtr    itsDateStrip;

	Type                itsType;

	void CopyNoCleanUp( const ARM_GramFctorArg& rhs );
	void CheckType( Type t ) const;
};


CC_END_NAMESPACE()

#endif

/*----------------------------------------------------------------------------*/
/*---- End of file ----*/

