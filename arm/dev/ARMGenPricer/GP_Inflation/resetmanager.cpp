/*
 * $Log: resetmanager.cpp,v $
 * Revision 1.7  2003/11/06 13:14:32  ebenhamou
 * change to have day display in date in view
 *
 * Revision 1.6  2003/09/29 09:50:26  ebenhamou
 * filtering of zero and correction in sprintf
 *
 * Revision 1.5  2003/08/21 08:42:17  ebenhamou
 * added method GetReset
 *
 * Revision 1.4  2003/08/20 08:42:29  ebenhamou
 * added constructor from excel
 *
 * Revision 1.3  2003/07/29 08:52:21  ebenhamou
 * change following code review of MAB
 *
 * Revision 1.1  2003/07/16 07:01:37  ebenhamou
 * Initial revision
 *
 *
 *
 */


#include <glob/firsttoinc.h>
#include "gpinflation/resetmanager.h"

/// standard libs
#include <utility>
#include <stdarg.h>
#include <algorithm>

///	flag for interp
#include <glob/armdef.h>

CC_BEGIN_NAMESPACE( ARM )

////////////////////////////////////////////////////
/// Note that when using iterator
/// cannot use iter->method but rather
/// (*iter).method for Solaris compiler to work
////////////////////////////////////////


////////////////////////////////////////////////////
////////////////////////////////////////////////////
///  single index reset
////////////////////////////////////////////////////
////////////////////////////////////////////////////


////////////////////////////////////////////////////
///	Class  : OneResetHistory
///	Routine: Constructor
///	Returns: Built object
///	Action : a one reset History is basically a simple class 
///				that contains a vector of dates transformed into
///				double and a vector of values
////////////////////////////////////////////////////

OneResetHistory::OneResetHistory( const string& indexName,
	const vector<double>& dates, 
	const vector<double>& values )
	: itsIndexName( indexName )
{
	for(int i=0; i<values.size(); i++ )
		itsResets.insert( pair<const double,double>( dates[i], values[i] ) );
	
	doubleDoubleMap::const_iterator firstpos = itsResets.begin();
	doubleDoubleMap::const_iterator lastpos  = itsResets.end();

	/// get the first
	/// and last reset 
	/// only if the history 
	/// is not empty
	if( firstpos != lastpos )
	{
		lastpos--;
		/// see the explanation above about the fact that Solaris compiler requires
		/// to use (*iter).first rather than iter->first
		itsfirstDate	= (*firstpos).first;		
		itslastDate		= (*lastpos).first;
	}
}


////////////////////////////////////////////////////
///	Class  : OneResetHistory
///	Routine: OneResetHistory
///	Returns: Built object
///	Action : copy constructor
////////////////////////////////////////////////////

OneResetHistory::OneResetHistory( const OneResetHistory& rhs )
: itsIndexName( rhs.itsIndexName), itsResets( rhs.itsResets ),
	itsfirstDate( rhs.itsfirstDate ), itslastDate( rhs.itslastDate )
{}


////////////////////////////////////////////////////
///	Class  : OneResetHistory
///	Routine: operator =
///	Returns: 
///	Action : operator =
////////////////////////////////////////////////////

OneResetHistory& OneResetHistory::operator = ( const OneResetHistory& rhs )
{
	if (this == &rhs )
		return *this;
	else
	{
		itsIndexName	= rhs.itsIndexName;
		itsResets		= rhs.itsResets;
		itsfirstDate	= rhs.itsfirstDate;
		itslastDate		= rhs.itslastDate;
	}
	return *this;
}


////////////////////////////////////////////////////
///	Class  : OneResetHistory
///	Routine: GetIndexName
///	Returns: string
///	Action : Get Index Name (accessor)
////////////////////////////////////////////////////

string OneResetHistory::GetIndexName() const
{
	return itsIndexName;
}



////////////////////////////////////////////////////
///	Class  : OneResetHistory
///	Routine: GetReset
///	Returns: double
///	Action : Get Reset (accessor)
////////////////////////////////////////////////////

double OneResetHistory::GetReset( double julianDate ) const
{
	doubleDoubleMap::const_iterator iter = itsResets.find( julianDate );
	if(iter == itsResets.end())
	{
		ARM_Date tmpDate( julianDate );
		char dateChar[20];
		tmpDate.JulianToStrDateDay( dateChar );
		char msg[255];
		sprintf( msg, "date : in its julian format  %f  corresponding to %s not found in the data with index %s",
			julianDate,dateChar, itsIndexName.c_str() );
		throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, msg );
	}

	return 	(*iter).second;
}


////////////////////////////////////////////////////
////////////////////////////////////////////////////
///  multiple indexes
////////////////////////////////////////////////////
////////////////////////////////////////////////////



////////////////////////////////////////////////////
///	Class  : ARM_ResetManager
///	Routine: ARM_ResetManager
///	Returns: 
///	Action : Constructor
////////////////////////////////////////////////////

ARM_ResetManager::ARM_ResetManager ( const rawinputMap& data )
{
	rawinputMap::const_iterator pos;

	for( pos = data.begin(); pos != data.end(); pos++)
	{
		oneMarketRawInputMap::const_iterator permktpos;
		oneMarketDataMap	oneMarketData;
		
		for( permktpos= ((*pos).second).begin(); permktpos != ((*pos).second).end(); permktpos++)
		{
			OneResetHistory singleReset( 
				(*permktpos).first, 
				((*permktpos).second).first, 
				((*permktpos).second).second );
			oneMarketData.insert( pair<const string, OneResetHistory >( 
				(*permktpos).first, 
				singleReset ) );
		}

		itsMktDataResets.insert( pair<const string, oneMarketDataMap >( 
			(*pos).first, 
			oneMarketData ) );
	}

	/// type is ARM_RESETMANAGER
	SetName( ARM_RESETMANAGER );
}


////////////////////////////////////////////////////
///	Class  : ARM_ResetManager
///	Routine: ARM_ResetManager
///	Returns: 
///	Action : Constructor using the data from excel
////////////////////////////////////////////////////

ARM_ResetManager::ARM_ResetManager( const vector<string>& data, int nbRows, int nbColumns )
{
	/// data are rows by rows
	int beginPos	= -1;
	int endPos		= 0;
	bool hasFound = false;

	for(;;)
	{
		/// beginPos is the column corresponding to "MKT"
		++beginPos;
		while( beginPos<nbColumns && strcmp(ARM_UPPERCASE( const_cast<char*>(data[beginPos].c_str())),"MKT")!= 0 )
			++beginPos;

		/// are we done if os exit the infinite loop!
		if( beginPos == nbColumns )
			break;

		/// find the next MKT column
		endPos = beginPos+2;
		while( endPos<nbColumns && strcmp(ARM_UPPERCASE(const_cast<char*>(data[endPos].c_str())),"MKT")!= 0 )
			++endPos;

		hasFound = true;
		string marketTag = data[beginPos+1];

		/// first get the dates
		int row;

		vector<double> dates;
		dates.reserve(nbRows-2);
		for( row = 2; row<nbRows; row++ )
			dates.push_back( atof( data[row*nbColumns+beginPos].c_str()) );

		oneMarketDataMap	oneMarketData;
		///then fills the marketDataMap
		int col;
		for( col= beginPos+1; col< endPos; col++ )
		{
			string tag = data[nbColumns+col];
			vector<double> values;
			values.reserve(nbRows-2);

			/// filter zero only for Inflation Index!
			double singleValue;
			for( row = 2; row<nbRows &&  ( ( singleValue = atof(data[row*nbColumns+col].c_str() ) ) != 0 || (marketTag == "IR") ) ; row++)
				values.push_back( singleValue );
			
			/// take the corresponding dates
			vector< double > singleDates( dates.begin(), dates.begin() +row-2 );
			OneResetHistory singleReset( tag, singleDates, values );

			oneMarketData.insert( pair<const string, OneResetHistory >( 
				tag, singleReset ) );
		}
		
		// insert the Mkt data
		itsMktDataResets.insert( pair<const string, oneMarketDataMap >( 
			marketTag, 
			oneMarketData ) );
	}
	
	if( !hasFound )
		throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, "data do not contain any market!" );
	
	/// type is ARM_RESETMANAGER
	SetName( ARM_RESETMANAGER );
}



////////////////////////////////////////////////////
///	Class  : ARM_ResetManager
///	Routine: ARM_ResetManager
///	Returns: 
///	Action : Copy constructor
////////////////////////////////////////////////////

ARM_ResetManager::ARM_ResetManager(const ARM_ResetManager& src)
: ARM_Object(src), itsMktDataResets( src.itsMktDataResets )
{}


////////////////////////////////////////////////////
///	Class  : ARM_ResetManager
///	Routine: Operator = 
///	Returns: 
///	Action : Operator = 
////////////////////////////////////////////////////

ARM_ResetManager& ARM_ResetManager::operator = ( const ARM_ResetManager& src )
{
	if( this != &src )
	{
		ARM_Object::operator=( src );
		itsMktDataResets = src.itsMktDataResets;
	}
	return *this;
}


////////////////////////////////////////////////////
///	Class  : ARM_ResetManager
///	Routine: Clone
///	Returns: 
///	Action : clone the object
////////////////////////////////////////////////////

ARM_Object* ARM_ResetManager::Clone(void)
{
	return new ARM_ResetManager( *this );
}

////////////////////////////////////////////////////
///	Class  : ARM_ResetManager
///	Routine: ~ARM_ResetManager
///	Returns: 
///	Action : Destructor
////////////////////////////////////////////////////

ARM_ResetManager::~ARM_ResetManager()
{}


////////////////////////////////////////////////////
///	Class  : ARM_ResetManager
///	Routine: GetReset
///	Returns: double
///	Action : accessor
////////////////////////////////////////////////////

double ARM_ResetManager::GetReset( double julianDate, const string& indexName, const string& market ) const
{
	MktMap::const_iterator mktIter = itsMktDataResets.find( market );
	if( mktIter == itsMktDataResets.end() )
	{
		char msg[255];
		sprintf( msg, "Mkt %s  not found in the reset Manager", market.c_str() );
		throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, msg );
	}

	oneMarketDataMap::const_iterator oneResetIter = ((*mktIter).second).find( indexName );
	if( oneResetIter == ((*mktIter).second).end() )
	{
		char msg[255];
		sprintf( msg, "indexName %s  not found in the reset Manager", indexName.c_str() );
		throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, msg );
	}

	return ((*oneResetIter).second).GetReset( julianDate );
}


////////////////////////////////////////////////////
///	Class  : ARM_ResetManager
///	Routine: GetReset
///	Returns: double
///	Action : accessor for just one index
////////////////////////////////////////////////////

double ARM_ResetManager::GetReset( double julianDate) const
{
	return ((((itsMktDataResets.begin())->second).begin())->second).GetReset( julianDate );
}



////////////////////////////////////////////////////
///	Class  : ARM_ResetManager
///	Routine: View
///	Returns: 
///	Action : view the details of the curve
///				a Mkt is a map with a market tag
///					INF and OneMktData
///					OneMktData is a map with
///					Indexes and doubleDoubleMap
///					a doubleDoubleMap is indeed a collection of dates 
///					and values corresponding
////////////////////////////////////////////////////

void ARM_ResetManager::View(char* id, FILE* ficOut)
{

	/// to be consistent with other interface
	/// need to use fprintf
    FILE* fOut;
    char fOutName[40];
	
	/// do we have already a file opened?
    if ( ficOut == NULL )
    {
		ARM_GetViewFile(ARM_VIEW_FILE, id, fOutName);
		fOut = fopen(fOutName, "w");
    }
    else
		fOut = ficOut;

    fprintf(fOut, "\n\n\t =====> Reset Manager \n\n");

	MktMap::const_iterator iter		= itsMktDataResets.begin();
	MktMap::const_iterator iterEnd	= itsMktDataResets.end();

	/// loop over markets
	for( ; iter != iterEnd; ++iter )
	{
		fprintf(fOut, "Markets:\t%s\n", ((*iter).first).c_str() );

		// is it empty?
		if( ((*iter).second).empty() )
			break;

		oneMarketDataMap::const_iterator oneMktDataIter		= ((*iter).second).begin();
		oneMarketDataMap::const_iterator oneMktDataIterEnd	= ((*iter).second).end();

		/// loop over indexes for a given market
		for( ; oneMktDataIter != oneMktDataIterEnd; ++oneMktDataIter )
		{
			OneResetHistory resetHistory = (*oneMktDataIter).second;
			// is it empty?
			if( resetHistory.empty() )
				break;

			doubleDoubleMap	resets = resetHistory.GetResets();
			fprintf(fOut, "\nIndex:\t%s\n", ((*oneMktDataIter).first).c_str() );
			fprintf(fOut, "Dates\t Value\n");
			doubleDoubleMap::const_iterator resetHistoryIter	= resets.begin();
			doubleDoubleMap::const_iterator resetHistoryIterEnd = resets.end();

			char date[20];
			
			// loop over dates and value for one index
			for( ; resetHistoryIter != resetHistoryIterEnd; ++resetHistoryIter )
			{
				((ARM_Date) ((*resetHistoryIter).first)).JulianToStrDateDay(date);
				fprintf(fOut, "%s\t%f\n", date, (*resetHistoryIter).second );
			}

			fprintf(fOut, "\n" );
		}
	}

	if ( ficOut == NULL )
		fclose(fOut);
}


CC_END_NAMESPACE()

/*---------------------------------------------------------------*/
/*---- End Of File ----*/

