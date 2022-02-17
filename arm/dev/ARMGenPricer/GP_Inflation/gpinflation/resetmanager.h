/*
 * $Log: resetmanager.h,v $
 * Revision 1.8  2003/09/29 08:17:15  ebenhamou
 * use the macro CC_USING_NS
 *
 * Revision 1.7  2003/09/22 13:51:58  ebenhamou
 * more strict using directive
 *
 * Revision 1.6  2003/08/21 08:42:53  ebenhamou
 * more explicit message
 *
 * Revision 1.4  2003/08/05 08:30:41  ebenhamou
 * change Getits into Get
 *
 * Revision 1.3  2003/07/29 08:52:28  ebenhamou
 * change following code review of MAB
 *
 * Revision 1.1  2003/07/16 07:01:28  ebenhamou
 * Initial revision
 *
 * Revision 1.4  2003/06/30 17:14:03  ebenhamou
 * dos2unix
 *
 * Revision 1.3  2003/06/30 16:11:13  ebenhamou
 * discounting version
 *
 */

/*!
	
    Declaration of the reset manager which is defined
	as a map per market of a map of OneResetHistory
*----------------------------------------------------------------------------*/ 

#ifndef _INGPINFLATION_INFRESETMANAGER_H
#define _INGPINFLATION_INFRESETMANAGER_H

#include <glob/firsttoinc.h>

#include <string>
#include <map>
#include <vector>

#include <glob/dates.h>		// because of the use of date

#include "gpbase/port.h"
CC_USING_NS( std, map )
CC_USING_NS( std, pair )
CC_USING_NS( std, string )
CC_USING_NS( std, make_pair )
CC_USING_NS( std, less )

typedef map< double, double, less< double> > doubleDoubleMap;

CC_BEGIN_NAMESPACE( ARM )


/*
 * for one single reset
 */
class OneResetHistory
{
private:
	string			itsIndexName;
	doubleDoubleMap	itsResets;
	double			itsfirstDate;
	double			itslastDate;
public:
	OneResetHistory( const string& indexName, 
		const vector<double>& dates, 
		const vector<double>& values );

	double GetFirstDate() const { return itsfirstDate;}
	double GetLastDate() const { return itslastDate;}
	OneResetHistory( const OneResetHistory& src );
	OneResetHistory& operator = ( const OneResetHistory& src );
	string GetIndexName() const;
	doubleDoubleMap GetResets() const { return itsResets; }
	bool empty() const { return itsResets.empty(); }
	double GetReset( double julianDate ) const;
};


typedef pair< vector<double>, vector<double> >					pairvdouble;
typedef map< string, pairvdouble, less<string> >				oneMarketRawInputMap;
typedef map< string, OneResetHistory, less<string>  >			oneMarketDataMap;
typedef map< string, oneMarketRawInputMap, less<string> >		rawinputMap;


/*
 * class to hold the various resets
 * history
 */

class ARM_ResetManager: public ARM_Object 
{
private:
	// only used judiciously hence private
	ARM_ResetManager();
	typedef map< string, oneMarketDataMap, less< string> > MktMap;
	MktMap itsMktDataResets;

public:
	ARM_ResetManager( const vector<string>& data, int columns, int rows  );
	ARM_ResetManager( const rawinputMap& data  );
	ARM_ResetManager( const ARM_ResetManager& src );
	ARM_ResetManager& operator = ( const ARM_ResetManager& src );
	virtual ~ARM_ResetManager();

	/// function to compute a reset
	double GetReset( double julianDate, const string& indexName, const string& market ) const;
	double GetReset( double julianDate) const; // for only one index in the resetManager

	OneResetHistory GetResets() const { return (((itsMktDataResets.begin())->second).begin())->second; }
    
	/// standard ARM_Object support
	virtual void View(char* id = NULL, FILE* ficOut = NULL);
	virtual ARM_Object* Clone();

};

CC_END_NAMESPACE()

#endif
