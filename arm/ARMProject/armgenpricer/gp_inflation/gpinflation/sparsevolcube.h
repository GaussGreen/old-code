/*
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 * $Log: sparsevolcube.h,v $
 * Revision 1.8  2003/11/20 15:24:05  ebenhamou
 * added last known date
 *
 * Revision 1.7  2003/11/19 20:00:19  ebenhamou
 * remove inappropriate struct
 *
 * Revision 1.6  2003/10/30 20:31:33  ebenhamou
 * change for DotNet compilation
 *
 * Revision 1.5  2003/09/23 17:32:26  ebenhamou
 * change spotTime to public
 *
 * Revision 1.4  2003/09/22 13:51:25  ebenhamou
 * stricter using directive
 *
 * Revision 1.3  2003/09/11 10:43:18  ebenhamou
 * meaningful static const name
 *
 * Revision 1.2  2003/09/09 12:56:40  ebenhamou
 * move constant to static variable
 *
 * Revision 1.1  2003/09/09 11:42:05  ebenhamou
 * Initial revision
 *
 *
 */

    
/*----------------------------------------------------------------------------*/

/*! \file sparsevolcube.h
 *
 *  \brief sparse vol cube
 *
 *	\author  Eric Benhamou
 *	\version 1.0
 *	\date September 2003
 */

 
/*----------------------------------------------------------------------------*/


/*! \class   ARM_SparseVolCube
 *	\brief  defines an object that fills only partially a vol cube
 *			a sparse vol cube has the ability to be filled according
 *			to various dimension of the cube with consistency check!
 *	\author  Eric Benhamou
 *	\version 1.0
 *	\date August 2003
 */

#ifndef _INGPINFLATION_SPARSEVOLCUBE_H
#define _INGPINFLATION_SPARSEVOLCUBE_H

#include <vector>
#include "volint.h"

/// gpbase
#include "gpbase/port.h"
#include "gpbase/rootobject.h"
#include "gpbase/gplinalgtypedef.h"

CC_USING_NS( std, vector )
CC_USING_NS( std, pair )
//CC_USING_NS_BI_T( std, pair, int, int )

typedef pair< int, int > pairIntInt;
/// forward declaration
class ARM_VolCube;

CC_BEGIN_NAMESPACE( ARM )

//// a sparse vol cube is just a little object to help to 
//// build a vol cube
//// it has almost no functionality
class ARM_SparseVolCube : public ARM_RootObject
{
public:
	static const double spotTime;
private :
	vector< ARM_VolCurve* >* itsVols;
	ARM_GP_Vector* itsUnderlyings;
	ARM_Date itsLastKnownDate;
	ARM_Date itsAsOf;
	
	void CleanUp();
	void CopyNoCleanUp( const ARM_SparseVolCube& rhs );
	void Init();
	void AddOneVolCurve( ARM_VolCurve* volCurv, double underlying );

	/// some helper function that can be used very generally!
	static vector< pairIntInt > intersectPos( const ARM_GP_Vector& v1, const ARM_GP_Vector& v2 );
	static void compare( ARM_VolCurve* nVolCrv, ARM_VolCurve* oVolCrv, double tenor );
	static vector< pairIntInt > mergePos( const ARM_GP_Vector& v1, const ARM_GP_Vector& v2, ARM_GP_Vector& writeResult );
	static ARM_VolCurve* merge( ARM_VolCurve* nVolCrv, ARM_VolCurve* oVolCrv, double tenor );
	void Validation( ARM_GP_Vector* dim1, ARM_GP_Vector* strikes, int dim1Type, int otherDimType );

public:
	ARM_SparseVolCube( 
		const ARM_Date& asOf,
		const ARM_Date& lastKnownDate,
		const string& indexName,
		int dim1Type,
		ARM_GP_Vector* dim1,
		int otherDimType,
		double otherDimValue,
		ARM_GP_Vector* strikes,
		ARM_Matrix* volatilities,
		int strikeType,
		int volType );

	ARM_SparseVolCube( const ARM_SparseVolCube& rhs );
	
	ARM_SparseVolCube& operator=( const ARM_SparseVolCube& rhs );

	virtual ~ARM_SparseVolCube();

	/// ARM_Obnject stuff
	virtual ARM_Object* Clone() const;
	virtual void View(char* id = NULL, FILE* ficOut = NULL) const;
    virtual string toString(const string& indent="", const string& nextIndent="") const { return "ARM_SparseVolCube"; }


	/// fill methods
	void AddVolCurves(
		int dim1Type,
		ARM_GP_Vector* dim1,
		int otherDimType,
		double otherDimValue,
		ARM_GP_Vector* strikes,
		ARM_Matrix* volatilities );

	/// accessors
	vector < ARM_VolCurve* > * GetVols() const;
	ARM_GP_Vector*	GetUnderlyings() const;

	inline ARM_Date GetAsOf() const { return itsAsOf; }
	inline ARM_Date GetLastKnownDate() const { return itsLastKnownDate; }

	/// conversion to a volcube
	ARM_VolCube* ConvertToVolCube() const;
};

CC_END_NAMESPACE()

#endif
/*----------------------------------------------------------------------------*/
/*---- End Of File ----*/
