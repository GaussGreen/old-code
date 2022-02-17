/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 *	\file dealdescription.cpp
 *
 *  \brief
 *
 *	\author  E. Benhamou
 *	\version 1.0
 *	\date October 2003
 */

/// gpinfra
#include "gpinfra/dealdescription.h"
#include "gpinfra/retcppcode.h"

#include <glob/expt.h>
#include <glob/dates.h>

/// gpbase
#include "gpbase/ostringstream.h"
#include "gpbase/gpvector.h"
#include "gpbase/stringmanip.h"

#ifndef unix
	#include <iomanip>
#endif


CC_BEGIN_NAMESPACE( ARM )


////////////////////////////////////////////////////
///	Class  : ARM_DealDescription
///	Routine: Constructor
///	Returns: nothing
///	Action : Constructor
////////////////////////////////////////////////////
ARM_DealDescription::ARM_DealDescription( const ARM_StringVector& txt, const vector< ARM_GP_VALUE_TYPE >& format,
	size_t RowsNb, size_t ColsNb,  const ARM_StringVector& pricedColumns)
: ARM_RootObject(), itsText( txt ), itsRowsNb( RowsNb ), itsColsNb( ColsNb ), itsFormat( format )
{
	/// some validation
	if( itsText.size() != itsRowsNb * itsColsNb )
		throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, 
		"ARM_DealDescription(): invalid vector of string compared to the nb of cols and rows" );
	
	if( itsFormat.size() != itsRowsNb * itsColsNb )
		throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, 
		"ARM_DealDescription(): invalid vector of format compared to the nb of cols and rows" );

	FilterBlank();
	SetName( ARM_DEALDES );

	if (pricedColumns.size() != 0)
	{
		itsPricedColumnNames.resize(pricedColumns.size());
		itsPricedColumns.resize(pricedColumns.size());
		for (size_t i = 0; i < pricedColumns.size(); ++i)
		{
			itsPricedColumnNames[i] = pricedColumns[i];
			itsPricedColumns[i] = GetColIndex(pricedColumns[i]);
		}
	}
	else
	{
		itsPricedColumns.resize(1);
		itsPricedColumns[0] = itsColsNb-1;
	}
};


////////////////////////////////////////////////////
///	Class  : ARM_DealDescription
///	Routine: FilterBlank
///	Returns: vector of string with no blank
///	Action : 
////////////////////////////////////////////////////
void ARM_DealDescription::FilterBlank()
{
	/// we filter blank!
	/// do the checking on the first column to get nbofRows
	/// and on the first row to get nbofCols
	size_t i = 0, j =0 ;
	while( i<itsRowsNb  && itsFormat[i*itsColsNb] != ARM_MISSING_TYPE )
		++i;
	int tmpRowsNb = i;

	while( j<itsColsNb  && itsFormat[j] != ARM_MISSING_TYPE )
		++j;
	int tmpColsNb = j;


	/// test that these newtmpCols and tmpRows are
	/// making sense
	for( i=0; i<itsRowsNb; ++i )
	{
		int jStart = 0;
		if(i< tmpRowsNb)
			jStart = tmpColsNb;
		
		for(j=jStart; j<itsColsNb; ++j )
		{
			if( itsFormat[ i*itsColsNb+j ] != ARM_MISSING_TYPE )
			{
				CC_Ostringstream os;
				os << "ARM_DealDescription(): tried to filter blank "
					<< "but it seems that there are some problems with cell "  
					<< "row: " << i << " col: " << j << " ";
				throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, os.str() );
			}
		}
	}

	/// get the correct Text and Format
	ARM_StringVector tmpText( tmpRowsNb * tmpColsNb );
	vector< ARM_GP_VALUE_TYPE > tmpFormat( tmpRowsNb * tmpColsNb );

	for( i=0; i<tmpRowsNb; ++i )
		for( j=0; j<tmpColsNb; ++j )
		{
			tmpText[i*tmpColsNb+j]	 = itsText[ i*itsColsNb+j ];
			tmpFormat[i*tmpColsNb+j] = itsFormat[ i*itsColsNb+j ];
		}
	
	itsFormat.swap( tmpFormat );
	itsText.swap( tmpText );
	itsRowsNb = tmpRowsNb;
	itsColsNb = tmpColsNb;
}


////////////////////////////////////////////////////
///	Class  : ARM_DealDescription
///	Routine: GetRow
///	Returns: vector of string are written left to right up to bottom
///	Action : Accessor to a given row
////////////////////////////////////////////////////
ARM_CountedPtr< ARM_DealDescription > ARM_DealDescription::GetRow( int RowNb ) const 
{
	/// a little bit clumsy but returns the correct smart Pointor!
	ARM_DealDescription* rawPtr = new ARM_DealDescription( ARM_StringVector( itsText.begin() + itsColsNb* RowNb,
		itsText.begin() + itsColsNb * (RowNb+1) ),
		vector< ARM_GP_VALUE_TYPE >( itsFormat.begin() + itsColsNb* RowNb,
		itsFormat.begin() + itsColsNb * (RowNb+1) ),
		1, itsColsNb );

	return ARM_CountedPtr< ARM_DealDescription >( rawPtr );
}

////////////////////////////////////////////////////
///	Class  : ARM_DealDescription
///	Routine: GetSubDescription
///	Returns: pointor to a deal description
///	Action : select a deal description made of lines
///          from rowBegin to rowEnd and with the nbCols
///          first columns
////////////////////////////////////////////////////
ARM_DealDescriptionPtr ARM_DealDescription::GetSubDescription(size_t rowBegin, size_t rowEnd, size_t nbCols, const ARM_StringVector& PricedColumns ) const
{
    if(nbCols>itsColsNb || rowBegin>rowEnd || rowEnd>=itsRowsNb)
        throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, "Bad sub selection in the deal description");

    size_t descSize=nbCols*(rowEnd-rowBegin+2);
    ARM_StringVector subText(descSize);
    vector< ARM_GP_VALUE_TYPE > subFormat(descSize);

    /// Names
    for(int j=0;j<nbCols;++j)
    {
        subText[j]=itsText[j];
        subFormat[j]=itsFormat[j];
    }

    /// Values
    size_t offset=rowBegin*itsColsNb;
    size_t subOffset=nbCols;
    for(size_t i=rowBegin;i<=rowEnd;++i)
    {
        for(size_t j=0;j<nbCols;++j)
        {
            subText[subOffset+j]=itsText[offset+j];
            subFormat[subOffset+j]=itsFormat[offset+j];
        }
        subOffset += nbCols;
        offset += itsColsNb;
    }

    return ARM_DealDescriptionPtr(new ARM_DealDescription(subText,subFormat,rowEnd-rowBegin+2,nbCols, PricedColumns ) );
}

////////////////////////////////////////////////////
///	Class  : ARM_DealDescription
///	Routine: toString,detailledString,toStringCommon
///	Returns: stringify the object
///	Action : 
////////////////////////////////////////////////////
string ARM_DealDescription::PrintOneColumn( size_t jMin, size_t jMax, bool detailMode ) const
{
	size_t i,j;
	CC_Ostringstream os;
	const int spaceMin			= 9;
	const int nbBlankToSepare	= 3;
	ARM_IntVector txtWidth( jMax-jMin, spaceMin );

	for(i=0; i<itsRowsNb; ++i)
	{
		for(j=jMin; j<jMax; ++j)
		{
			/// for non date type get the max of text
			if( itsFormat[i*itsColsNb+j] != ARM_DATE_TYPE )
			{
				if( itsText[i*itsColsNb+j].size() > txtWidth[j-jMin] )
					txtWidth[j-jMin] = itsText[i*itsColsNb+j].size();
			}
			/// for date generates the date and takes the max of its length
			else
			{
				string dateStr = ARM_Date( atof( itsText[i*itsColsNb+j].c_str() ) ).toString();
				if( dateStr.size() > txtWidth[j-jMin] )
					txtWidth[j-jMin] = dateStr.size();
			}

			if( ARM_TYPE::Name( itsFormat[i*itsColsNb+j] ).size() > txtWidth[j-jMin] )
				txtWidth[j-jMin] = ARM_TYPE::Name( itsFormat[i*itsColsNb+j] ).size();
		}
	}

	os << CC_NS( std, endl );
	for(i=0; i<itsRowsNb; ++i)
	{
		for(j=jMin; j<jMax; ++j)
		{
			os  << CC_NS( std, left ) << CC_NS( std, setfill )(' ') 
				<< CC_NS( std, setw)( txtWidth[j-jMin] );

			if( itsFormat[i*itsColsNb+j] == ARM_DATE_TYPE )
			{
				string dateStr = ARM_Date( atof( itsText[i*itsColsNb+j].c_str() ) ).toString();
				os << dateStr << "\t";
			}
			else
				os << itsText[i*itsColsNb+j]  << "\t";
		}

		os << CC_NS( std, endl );

		if( detailMode )
		{
			for(j=jMin; j<jMax; ++j)
			{
				os  << CC_NS( std, left ) << CC_NS( std, setfill )(' ') 
					<< CC_NS( std, setw)( txtWidth[j-jMin] );
				os << ARM_TYPE::Name( itsFormat[i*itsColsNb+j] ) << "\t";
			}

			os << CC_NS( std, endl );
		}
	}
	return os.str();
}



////////////////////////////////////////////////////
///	Class  : ARM_DealDescription
///	Routine: toString,detailledString,toStringCommon
///	Returns: stringify the object
///	Action : 
////////////////////////////////////////////////////
string ARM_DealDescription::toString(const string& indent, const string& nextIndent ) const
{ 
#ifdef _DEBUG
	return toStringCommon(indent,nextIndent,true); 
#else
	return toStringCommon(indent,nextIndent,false); 
#endif
}

string ARM_DealDescription::detailledString(const string& indent, const string& nextIndent ) const
{ return toStringCommon(indent,nextIndent,true); }

string ARM_DealDescription::toStringCommon(const string& indent, const string& nextIndent, bool detailMode) const
{

	CC_Ostringstream os;
	os << "Deal Description Object\n" ;

	if( !itsPricedColumnNames.empty() )
	{
		os << indent << " -Columns " << itsPricedColumnNames.size() << "\n";
		for( size_t i=0; i<itsPricedColumnNames.size(); ++i )
			os << indent << " " << itsPricedColumnNames[i] << "\n";
		os << "\n";
	}
	
	/// get first the textWidth!
	int columnToPrint = 10, i;
	for( i=0; i<itsColsNb/columnToPrint; ++i )
	{
		os << PrintOneColumn(i*columnToPrint, (i+1)*columnToPrint, detailMode );
	}
	os << CC_NS( std, endl );

	if( itsColsNb> i*columnToPrint )
		os << PrintOneColumn( i*columnToPrint, itsColsNb, detailMode );

	return os.str();
}


////////////////////////////////////////////////////
///	Class  : ARM_DealDescription
///	Routine: Clone method
///	Returns: a new copy of this
///	Action : Copy this
////////////////////////////////////////////////////
ARM_Object* ARM_DealDescription::Clone() const
{
	return new ARM_DealDescription( *this ); 
}


////////////////////////////////////////////////////
///	Class  : ARM_DealDescription
///	Routine: Destructor
///	Returns: 
///	Action : 
////////////////////////////////////////////////////
ARM_DealDescription::~ARM_DealDescription()
{}


////////////////////////////////////////////////////
///	Class  : ARM_DealDescription
///	Routine: Copy Constructor
///	Returns: 
///	Action : 
////////////////////////////////////////////////////
ARM_DealDescription::ARM_DealDescription( const ARM_DealDescription& rhs)
:    ARM_RootObject( rhs ),
	 itsText( rhs.itsText ), 
	 itsRowsNb( rhs.itsRowsNb ), 
	 itsColsNb( rhs.itsColsNb ),
	 itsFormat( rhs.itsFormat ), 
	 itsPricedColumns(rhs.itsPricedColumns), 
	 itsPricedColumnNames(rhs.itsPricedColumnNames)
{}

////////////////////////////////////////////////////
///	Class  : ARM_DealDescription
///	Routine: Assignment operator
///	Returns: 
///	Action : 
////////////////////////////////////////////////////
ARM_DealDescription& ARM_DealDescription::operator=( const ARM_DealDescription& rhs )
{
	if( this !=	 &rhs )
	{
	    ARM_RootObject::operator = ( rhs );
		itsText				=  rhs.itsText;
		itsFormat			=  rhs.itsFormat;
		itsRowsNb			=  rhs.itsRowsNb;
		itsColsNb			=  rhs.itsColsNb;
		itsPricedColumns	= rhs.itsPricedColumns;
		itsPricedColumnNames = rhs.itsPricedColumnNames;
	}
	return *this;
}

size_t ARM_DealDescription::GetColIndex(const string& colName) const
{
	for (int i = 0; i < itsColsNb; ++i)
		if ( stringGetUpper(itsText[i]) == stringGetUpper(colName) )
			return i;

	CC_Ostringstream os;
	os << "ARM_DealDescription(): the colum name "
		<< colName << " doesn't exist in the deal description";
	throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, os.str() );
}

CC_END_NAMESPACE()


/*---------------------------------------------------------------------------*/
/*---- End of file ----*/

