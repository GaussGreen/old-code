/*
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 * $Log: world.h,v $
 * Revision 1.1  2004/03/01 07:52:17  ebenhamou
 * Initial revision
 *
 */


/*! \file .h
 *
 *  \brief
 *	\author  E. Benhamou
 *	\version 1.0
 *	\date March 2004
 */


#ifndef _INGPBASE_DATESTRIPCOMBINER_H
#define _INGPBASE_DATESTRIPCOMBINER_H

#include "port.h"
#include "typedef.h"
#include "datestrip.h"

#include <string>
CC_USING_NS(std,string)

CC_BEGIN_NAMESPACE( ARM )

class ARM_DateStripCombiner : public ARM_RootObject
{
public:
	static const double DateStripCombiner_BlankData;

    ARM_DateStripCombiner( const ARM_DateStripVector& dateStrips = ARM_DateStripVector(), const string& MergeDateFuncString = "ResetDate" );
	ARM_DateStripCombiner( const ARM_DateStripCombiner& rhs );
	ARM_DateStripCombiner& operator=( const ARM_DateStripCombiner& rhs );
	virtual ~ARM_DateStripCombiner();

	/// Accessor
	size_t size() const { return itsContents.size(); }
	ARM_DateStripPtr GetDateStrip( size_t i ) const { return itsContents[i]; }
	ARM_VectorPtr GetMergeData() const { return itsMergeData; }

	/// standard ARM Object support
	/// because of the clone implementation in terms of the copy constructor
	/// there is no need to define BitwiseCopy and Copy!
    virtual string toString(const string& indent="",const string& nextIndent="") const;
	virtual ARM_Object* Clone() const;

private:
	typedef std::vector<double>* (ARM_DateStrip::*MergeFunc)() const;
	MergeFunc GetMergeFuncFromString( const string& MergeDateFuncString ) const;

	vector<ARM_DateStripPtr> itsContents;
	ARM_VectorPtr itsMergeData;

	/// function to merge all the data
	void MergeData( const MergeFunc& func );

	/// function for formatting
	static string PrintOneDateStripDateElem( double elem );
	static string PrintOneDateStripDoubleElem( double elem );


};

CC_END_NAMESPACE()

#endif

/*----------------------------------------------------------------------------*/
/*---- End of file ----*/

