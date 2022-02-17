/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 *	\file dealdescription.h
 *
 *  \brief 
 *
 *	\author  E. Benhamou
 *	\version 1.0
 *	\date October 2003
 */

#ifndef _INGPINFRA_DEALDESCRIPTION_H
#define _INGPINFRA_DEALDESCRIPTION_H

/// this headers has to come first
/// as it contains pragma for stupid warning
#include "gpbase/removeidentifiedwarning.h"

#include "gpbase/port.h"
#include "gpbase/rootobject.h"
#include "gpbase/valuetype.h"

#include <vector>
CC_USING( std::vector )
#include <string>
CC_USING( std::string )

#include "typedef.h"


CC_BEGIN_NAMESPACE( ARM )

/// \class ARM_DealDescription
/// \brief little class to help to get the deal description
class ARM_DealDescription : public ARM_RootObject
{
private:
	ARM_StringVector itsText;
	vector< ARM_GP_VALUE_TYPE > itsFormat;
	size_t itsRowsNb;
	size_t itsColsNb;
	void FilterBlank();
	string toStringCommon(const string& indent, const string& nextIndent, bool detailMode) const;
	vector< int > itsPricedColumns;
	vector< string > itsPricedColumnNames;
	string PrintOneColumn( size_t jMin, size_t jMax, bool detailMode ) const;

public:
	ARM_DealDescription( const ARM_StringVector& txt, const vector< ARM_GP_VALUE_TYPE >& format, size_t RowsNb, size_t ColsNb,  const ARM_StringVector& PricedColumns = ARM_StringVector(0));
	virtual ~ARM_DealDescription();
	ARM_DealDescription( const ARM_DealDescription& );
	ARM_DealDescription& operator=( const ARM_DealDescription& );

	/// vector of string are written left to right up to bottom
	ARM_DealDescriptionPtr GetRow( int RowNb ) const;

    ARM_DealDescriptionPtr GetSubDescription(size_t rowBegin, size_t rowEnd, size_t nbCols, const ARM_StringVector& PricedColumns = ARM_StringVector(0) ) const;

	/// various accessors
	inline size_t size() const { return itsRowsNb * itsColsNb; }
	inline size_t GetRowsNb() const { return itsRowsNb; }
	inline size_t GetColsNb() const { return itsColsNb; }
	inline string GetElem( size_t row, size_t col = 0 ) const { return itsText[ row * itsColsNb + col ]; }
	inline ARM_GP_VALUE_TYPE GetElemFormat( size_t row, size_t col = 0 ) const { return itsFormat[ row * itsColsNb + col ]; }
	inline ARM_StringVector GetText() const { return itsText; }
	inline vector< ARM_GP_VALUE_TYPE > GetFormat() const { return itsFormat; }
    void SetElem(size_t row, size_t col,const string& elem,ARM_GP_VALUE_TYPE elemFormat) {itsText[ row * itsColsNb + col ]=elem;itsFormat[ row * itsColsNb + col ]=elemFormat;}


	size_t GetColIndex(const string& colName) const;

	// Priced Columns accessors
	inline size_t GetNbPricedColumns() const { return itsPricedColumns.size(); }
	inline int GetPricedColumn(size_t idx) const { return itsPricedColumns[idx]; }
	inline const vector<string>& GetPricedColumnNames() const { return itsPricedColumnNames; }

	/// standard ARM Object support
	/// because of the clone implementation in terms of the copy constructor
	/// there is no need to define BitwiseCopy and Copy!
	virtual ARM_Object* Clone() const;
    virtual string toString(const string& indent="", const string& nextIndent="") const;
	string detailledString(const string& indent, const string& nextIndent ) const;
};


CC_END_NAMESPACE()

#endif

/*----------------------------------------------------------------------------*/
/*---- End of file ----*/

