/*
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 * $Log: correlmanager.cpp,v $
 * Revision 1.18  2004/06/14 15:08:30  emezzine
 * check intramarket's tenor and intramarket'currencies separatly.
 *
 * Revision 1.17  2004/02/17 14:11:40  mcampet
 * MC add validation at CorrelMatrix level
 *
 * Revision 1.16  2004/02/13 14:38:46  mcampet
 *  MC add ARMCorrelMatrix builders
 *
 * Revision 1.15  2003/11/07 11:42:50  ebenhamou
 * added default correl
 *
 * Revision 1.14  2003/11/06 13:14:25  ebenhamou
 * change to have day display in date in view
 *
 * Revision 1.13  2003/10/27 07:17:06  ebenhamou
 *  change for .Net compatibility
 *
 * Revision 1.11  2003/10/23 09:36:00  ebenhamou
 * more checking to avoid non working release xll
 *
 * Revision 1.10  2003/10/07 15:08:43  ebenhamou
 * remove ugly const_cast because added the correct const accessor in appropriate class
 *
 * Revision 1.8  2003/09/29 08:19:17  ebenhamou
 * using macro CC_USING_NS
 *
 * Revision 1.7  2003/09/26 17:05:12  ebenhamou
 * version with namespace handled by macro
 *
 * Revision 1.6  2003/09/25 07:36:49  ebenhamou
 * inequality has to be non strict
 *
 * Revision 1.5  2003/09/22 18:08:14  ebenhamou
 * added NOCHECK functionality
 *
 * Revision 1.3  2003/09/19 18:18:35  ebenhamou
 * more validation
 *
 * Revision 1.2  2003/09/18 17:10:10  ebenhamou
 * more validation
 *
 * Revision 1.1  2003/09/17 18:23:37  ebenhamou
 * Initial revision
 *
 *
 */

 
/*----------------------------------------------------------------------------*
    correlmanager.cpp
	Copyright (c) CDC IXIS CM July 2003 Paris
*----------------------------------
------------------------------------------*/

#include "correlmanager.h"
#include "dates.h"
#include "linalg.h"
#include "fromto.h"
#include <algorithm>
#include <set>

/// for unlink support
#ifdef unix
#include <sys/types.h>
#include <unistd.h>
#endif

using std::pair;
using std::string;
using std::vector;
using std::set;


/////////////////////////////////////////////////////////////////
/// function to returns the upper case version of a string
/////////////////////////////////////////////////////////////////
string stringToUpper( const string& s )
{
    // because of UNIX compatibility

    const char* strChar = s.c_str();
    
    char buf[100];

    strcpy(buf, strChar);

    (void) ARM_UPPERCASE(buf);

	// std::transform( a.begin(), a.end(), a.begin(), toupper );
    
    string res(buf);

    return(res);
}



/////////////////////////////////////////////////////////////////
//// function to check that a string is really in uppercase
/////////////////////////////////////////////////////////////////
void StringUpperCaseValidateWith( const string& s, const string& stringName )
{
	if( stringToUpper( s ) != s ) 
	{
		char msg[255];
		sprintf( msg, "the %s tag %s, is not in upper case! Please advise", stringName.c_str(), s.c_str() );
		throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, msg );
	}
}

bool splitAndFindString( const string& s, const string& delimiter, const string& sToFind)
{
	set<string> subTags;
	string::size_type pos = 0, prev_pos = 0;
	while( (pos = s.find_first_of(delimiter, pos ) ) != string::npos )
	{
		subTags.insert(s.substr( prev_pos, pos-prev_pos ) );
		prev_pos = ++pos;
	}
	subTags.insert(s.substr( prev_pos, pos-prev_pos ) );
	
	set<string>::iterator iter = subTags.find(sToFind);

	return ( iter != subTags.end() );
}



/*!
 * Default value centralised CENTRALISED
 */
double CorrelDataDefault::defaultCorrel = 0.7;




/*!
 * Local function to convert Correlation tenor to double
 */
double ConvertCOTenorToDouble(const char* matu)
{
    double resMatu = 0.0;
    int inMatu;

    char buf[20];
    char buf0[20];

    buf[0] = ' ';
    buf0[0] = ' ';
    buf0[1] = ' ';

    sscanf(matu, "%[A-Z,a-z]%d%s",buf0, &inMatu, buf);

    if (( toupper(buf[0]) == 'J' ) || ( toupper(buf[0]) == 'D' ))
    {
       resMatu = inMatu/365.0;
    }
    else if ( toupper(buf[0]) == 'M' )
    {
       resMatu = inMatu/12.0;
    }
    else if (( toupper(buf[0]) == 'A' ) || ( toupper(buf[0]) == 'Y' ))
    {
       resMatu = inMatu;
    }
    else
    {
       throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
         "===> Correlation tenor must look like : CO3D, CO2M, CO5Y,...EURIBOR3M,...");
    }

    return(resMatu);
}

/*!
 * Init method
 */
void ARM_CorrelMatrix::Init(void)
{
	SetName( ARM_CORRELMATRIX );
}


/*!
 * Constructor
 */
ARM_CorrelMatrix::ARM_CorrelMatrix( const ARM_Date& asOf, 
	ARM_Vector* X, ARM_Vector* Y, ARM_Matrix* Z )
:	ARM_VolLInterpol( asOf, X, Y, Z, K_STK_TYPE_PRICE, K_ATMF_VOL )
{
	Init( );
}

ARM_CorrelMatrix::ARM_CorrelMatrix( const ARM_VolLInterpol* correl )
:	ARM_VolLInterpol(*correl)
{
	Init( );
}

/*!
 * Constructor with a volint 
 * pseudo conversion operator
 * warning this is defined as explicit to avoid silent conversion
 */
ARM_CorrelMatrix::ARM_CorrelMatrix( const ARM_VolLInterpol& rhs )
:	ARM_VolLInterpol( rhs )
{
	Init();
}

ARM_CorrelMatrix::ARM_CorrelMatrix( const ARM_VolFlat* correl )
:   ARM_VolLInterpol(correl->GetAsOfDate(),
                     (ARM_Vector*) (new ARM_Vector(1,1))->Clone(),
                     (ARM_Vector*) (new ARM_Vector(1,((ARM_VolFlat*) correl)->GetVolatility()))->Clone() 
                     )
{
	 
    Init();
}


ARM_CorrelMatrix::ARM_CorrelMatrix( const ARM_VolCube* correl )
:   ARM_VolLInterpol(*correl)
{	 
    Init();
}


/*!
 * Copy Constructor
 */
ARM_CorrelMatrix::ARM_CorrelMatrix( const ARM_CorrelMatrix& rhs )
:	ARM_VolLInterpol( rhs )
{
	Init();
}


/*!
 * Operator = 
 */
ARM_CorrelMatrix& ARM_CorrelMatrix::operator = (const ARM_CorrelMatrix &rhs )
{
	if( this !=	 &rhs )
		ARM_VolLInterpol::operator = ( rhs );

	return *this;
}

/*!
 * Destructor: itsInfFwdCurve is not deleted because it is shared
 */
ARM_CorrelMatrix::~ARM_CorrelMatrix()
{
}

/*!
 * BitwiseCopy
 */
void ARM_CorrelMatrix::BitwiseCopy(const ARM_Object* src)
{
}

	
/*!
 * Copy
 */
void ARM_CorrelMatrix::Copy(const ARM_Object* src)
{
	ARM_VolLInterpol::Copy(src);
	BitwiseCopy(src);
}



/*!
 * Clone
 */
ARM_Object* ARM_CorrelMatrix::Clone()
{
	return new ARM_CorrelMatrix( *this );
}


/*!
 * View method ... very similar to an ARM_VolInt
 * with the difference of concept of strike and atm vol!
 */
void ARM_CorrelMatrix::View(char* id, FILE* ficOut)
{
    FILE* fOut;

    char fOutName[200];
    char strDate[50];
 
 
    if ( ficOut == NULL )
    {
       ARM_GetViewFile(ARM_VIEW_FILE, id, fOutName);
       (void) unlink(fOutName);
       fOut = fopen(fOutName, "w");
    }
    else
       fOut = ficOut;

    fprintf(fOut, "\n Correlation matrix\n");
    fprintf(fOut, "\n\t Currency : %s \n", GetCurrency()->GetCcyName());
    if ( GetInterpType() == K_DIAG_INTERPOL )
       fprintf(fOut, "\n\t Interpolation Mode : DIAG \n\n");
    else
       fprintf(fOut, "\n\t Interpolation Mode : LINEAR \n\n");

    GetAsOfDate().JulianToStrDateDay(strDate);
    fprintf(fOut, "\n AsOfDate  : %s \n\n", strDate);

    ARM_Vector* terms = GetExpiryTerms();
    ARM_Matrix* vols = GetVolatilities();
    int nbLg  = terms->GetSize();
    int nbCol = GetStrikes()->GetSize();

    int i, j;

    fprintf(fOut,"             ");
    for (j = 0; j < nbCol; j++)
        fprintf(fOut," %10.5lf", GetStrikes()->Elt(j));

    fprintf(fOut,"\n\n");

    for (i = 0; i < nbLg; i++)
    {
        fprintf(fOut," %10.5lf  ", terms->Elt(i));
        for (j = 0; j < nbCol; j++)
            fprintf(fOut," %10.5lf", vols->Elt(i, j));

        fprintf(fOut,"\n");
    }

	/// to allow to have nested view
    if ( ficOut == NULL )
       fclose(fOut);
}


/* !
 * Accessor on the individual correlData with some computation
 */
double ARM_CorrelMatrix::ComputeCorrelData( double x, double y ) const
{
	//// we do simple interpolation without any tenor
	/// reference hence the third term equal to zero!

	///////// AAAAAAAAAAAAAAAAAAAAAAARRRRRRGGGGGGGG
	///////// the poor const correctness of ARM is so painful!
	return const_cast<ARM_CorrelMatrix*>( this )->ComputeVolatility( x, y, 0.0 );
}

/* !
 * Accessor on the correlData matrix
 */
ARM_Matrix* ARM_CorrelMatrix::GetCorrelData( ) const
{
	return GetVolatilities();
}

/* !
 * Accessor on the dimension X
 */
ARM_Vector* ARM_CorrelMatrix::GetX() const
{
	return GetStrikes();
}


/* !
 * Accessor on the dimension Y
 */
ARM_Vector* ARM_CorrelMatrix::GetY() const
{
	return GetExpiryTerms();
}


///////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////
///////////// ARM_CorrelManager PART	///////////////////////////
///////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////

/*!
 * Constructor
 */
ARM_CorrelManager::ARM_CorrelManager( const string& mktTag, const string& intraMktTag, const ARM_Date& asOf, 
									  ARM_Vector* X, ARM_Vector* Y, ARM_Matrix* Z )
:	itsMktData()
{
	SetName( ARM_CORRELMANAGER );
	ARM_CorrelMatrix correlMatrix( asOf, X, Y, Z );
	Fill( mktTag, intraMktTag, correlMatrix );
}


/*!
 * Constructor with an ARM_VolLInterpol
 */
ARM_CorrelManager::ARM_CorrelManager( const string& mktTag, const string& intraMktTag, ARM_CorrelMatrix* correlMatrix )
{
	SetName( ARM_CORRELMANAGER );
	Fill( mktTag, intraMktTag, *correlMatrix );
}


ARM_CorrelManager::ARM_CorrelManager( const string& mktTag, const string& intraMktTag, ARM_VolFlat* volFlat )
{
	SetName( ARM_CORRELMANAGER );
	ARM_CorrelMatrix	correlMatrix(volFlat); 
	Fill( mktTag, intraMktTag, correlMatrix );
}

ARM_CorrelManager::ARM_CorrelManager( const string& mktTag, const string& intraMktTag, ARM_VolLInterpol* volInterp )
{
	SetName( ARM_CORRELMANAGER );
	ARM_CorrelMatrix	correlMatrix(volInterp); 
	Fill( mktTag, intraMktTag, correlMatrix );
}


ARM_CorrelManager::ARM_CorrelManager( const string& mktTag, const string& intraMktTag, ARM_VolCube* volCube )
{
	SetName( ARM_CORRELMANAGER );
	ARM_CorrelMatrix	correlMatrix(volCube); 
	Fill( mktTag, intraMktTag, correlMatrix );
}


 


/*!
 * Copy Constructor
 */
ARM_CorrelManager::ARM_CorrelManager( const ARM_CorrelManager& rhs )
:	ARM_Object( rhs )
{
	CopyNoCleanUp( rhs);
}


/*!
 * Operator = 
 */
ARM_CorrelManager& ARM_CorrelManager::operator = (const ARM_CorrelManager &rhs )
{
	if( this !=	 &rhs )
	{
		ARM_Object::operator = ( rhs );
		CleanUp();
		CopyNoCleanUp( rhs );
	}

	return *this;
}

/*!
 * Destructor: itsInfFwdCurve is not deleted because it is shared
 */
ARM_CorrelManager::~ARM_CorrelManager()
{
	CleanUp();
}



void ARM_CorrelManager::Init(void)
{
    SetName(ARM_CORRELMANAGER);
}


/*!
 * BitwiseCopy
 */
void ARM_CorrelManager::BitwiseCopy(const ARM_Object* src)
{
	ARM_CorrelManager* rhs = dynamic_cast<ARM_CorrelManager*>( const_cast<ARM_Object*>( src ) );
	if( rhs )
	{
		CleanUp();
		CopyNoCleanUp( *rhs );
	}
}

	
/*!
 * CleanUp
 */
void ARM_CorrelManager::CleanUp()
{
	/// nothing
}


/*!
 * CopyNoCleanUp
 */
void ARM_CorrelManager::CopyNoCleanUp( const ARM_CorrelManager& rhs )
{
	itsMktData= rhs.itsMktData;
}


/*!
 * Copy
 */
void ARM_CorrelManager::Copy(const ARM_Object* src)
{
	ARM_Object::Copy(src);
	BitwiseCopy(src);
}


/*!
 * Clone
 */
ARM_Object* ARM_CorrelManager::Clone()
{
	return new ARM_CorrelManager( *this );
}



/*
 * Fill method
 * after validation, it inserts the data
 * by design decision, we do not allow to reenter the same data!
 */
void ARM_CorrelManager::Fill( const string& tmpMktTag, const string& tmpIntraMktTag, const ARM_CorrelMatrix& correlMatrix)
{
	/// Does the uppercase conversion!
	string mktTag = stringToUpper(tmpMktTag );
	string intraMktTag	= stringToUpper(tmpIntraMktTag );

	/// validation and insertion: Validate only if there is already something to validate!
	ValidateData( mktTag, intraMktTag, correlMatrix );
	intraMarketCorrelsMap intraMarketCorrels = GetIntraMarketCorrels( mktTag );
	
	/// model comp
	if( ValidateCorrelMatrix(mktTag,intraMktTag,correlMatrix) == COMP )
		const_cast<ARM_CorrelMatrix&>(correlMatrix).SetInterpType( K_DIAG_INTERPOL );

	/// this is new tag, so use insert!
	intraMarketCorrels.insert( pair< const string, ARM_CorrelMatrix >( intraMktTag, correlMatrix ) );
	/// this can be a new or old tag, so use []!
	itsMktData[ mktTag ] = intraMarketCorrels;

	
}


/*!
 * Accesor to the AsOfDate
 * Action: get the asOfDate of the first Correlation matrix
 */
ARM_Date ARM_CorrelManager::GetAsOfDate() const
{
	if( itsMktData.empty() )
		throw Exception(__LINE__, __FILE__, ERR_INVALID_DATA, "you should never use this function with empty mktData, please advise!" );

	/// gets the asOfDate of the first index
	AllMarketCorrelsMap::const_iterator MktIter			= itsMktData.begin();
	intraMarketCorrelsMap::const_iterator intraMktIter	= ((*MktIter).second).begin();
	return ((*intraMktIter).second).GetAsOfDate();
};


/*!
 * Fill method with an ARM_VolLInterpol
 *
 */
void ARM_CorrelManager::Fill( const string& mktTag, const string& intraMktTag, ARM_CorrelMatrix* correlMatrix )
{
	/// check the asOfDate
	ARM_Date asOfDate = GetAsOfDate();
	ARM_Date correlAsOfDate = correlMatrix->GetAsOfDate();

	if( correlAsOfDate != asOfDate )
	{
		/// meaningful error message
		char msg[255];
		char asOfDateChar[20];
		char correlAsOfDateChar[20];
		asOfDate.JulianToStrDateDay( asOfDateChar );
		correlAsOfDate.JulianToStrDateDay( correlAsOfDateChar );
		sprintf( msg, "the volcurve taken has a different asOfDate %s than previous data %s, please advise", correlAsOfDateChar, asOfDateChar );
		throw Exception(__LINE__, __FILE__, ERR_INVALID_DATA, msg );
	}

	Fill(  mktTag, intraMktTag, *correlMatrix );
}



void ARM_CorrelManager::Fill( const string& mktTag, const string& intraMktTag, ARM_Vector* X, ARM_Vector* Y, ARM_Matrix* Z )
{
	ARM_Date asOf = GetAsOfDate();
	ARM_CorrelMatrix correlMatrix( asOf, X, Y, Z );
	Fill(  mktTag, intraMktTag, correlMatrix );
}


/*!
 * accessor for intraMarketCorrel for a given mktTag
 *
 */
ARM_CorrelManager::intraMarketCorrelsMap ARM_CorrelManager::GetIntraMarketCorrels( const string& tmpMktTag ) const
{
	/// empty or do not exist!
	if( itsMktData.empty() )
		return intraMarketCorrelsMap();

	/// else do a find with first the uppercase conversion
	string mktTag = stringToUpper(tmpMktTag );

	/// find algo of the STL
	AllMarketCorrelsMap::const_iterator mktIter = itsMktData.find( mktTag );
	if( mktIter == itsMktData.end() )
		return intraMarketCorrelsMap();
	else
		return (*mktIter).second;
};



/*!
 * Accessor that returns NULL if not find
 */
ARM_CorrelMatrix* ARM_CorrelManager::GetIntraMarketCorrel( const string& tmpMktTag, const string& tmpIntraMktTag ) const
{
	/// upper case conversion
	string mktTag = stringToUpper(tmpMktTag );
	string intraMktTag	= stringToUpper(tmpIntraMktTag );

	//// do the search
	AllMarketCorrelsMap::const_iterator mktIter = itsMktData.find( mktTag );

	/// do we have the correct data per mkt!
	if( mktIter == itsMktData.end() )
		return NULL;
	
	/// if so goahead with the intraMktTag
	intraMarketCorrelsMap::const_iterator oneMktIter = ( (*mktIter).second).find( intraMktTag );

	/// test that we could find something
	if( oneMktIter == ( (*mktIter).second).end() )
		return NULL;

	/// ok there is really something
	/// a little bit ugly but required to const_cast!
	return &( const_cast<ARM_CorrelMatrix&>( (*oneMktIter).second ) );
}

/*!
 * ComputeCorrelData does a lookup through the market data
 * to get the various data
 *
 */
double ARM_CorrelManager::ComputeCorrelData( const string& mktTag, const string& intraMktTag, double x, double y ) const
{
    //// in case somehow uses something like itsCorrelManager->ComputeCorrelData( blahblah ... ) with itsCorrelManager = NULL
    if( !this )
		throw Exception(__LINE__, __FILE__, ERR_INVALID_DATA, "you are trying to compute correl with an unitialised correl manager, please advise!" );

    
    if( itsMktData.size() == 1 )
    {
        AllMarketCorrelsMap::const_iterator mktIter = itsMktData.begin();
        
        if ( ((*mktIter).second).size() == 1)
        {
            intraMarketCorrelsMap::const_iterator oneMktIter = ( (*mktIter).second).begin();
            return  (const_cast<ARM_CorrelMatrix*> ( &(*oneMktIter).second))->ComputeCorrelData( x, y);
            
        }
    }
	return ComputeCorrelData( mktTag, intraMktTag )->ComputeCorrelData( x, y);
}



double ARM_CorrelManager::ComputeIndexIndexCorrelData(char* strccy1, const string& index1Name, 
                                                      char* strccy2, const string& index2Name,
                                                      double x1, double x2, bool aAcceptLiborMode)
{
    double vCorrel = 0.0;

	try
	{

		string intraMktTag;
    
		string ccy1(strccy1);
		string ccy2(strccy2);


		if ( ccy1 == ccy2 )
		{
		   if ( ConvertCOTenorToDouble(index1Name.c_str()) 
				< ConvertCOTenorToDouble(index2Name.c_str()) 
			  ) 
		   {
			  intraMktTag = ccy1+"_"+index1Name+"_"+ccy1+"_"+index2Name +"_"+"NOCHECK";

			  // because of Summit the correl matrix is (LIBOR, CMS) and
			  // because of the correl manager ordering :

			  if ((( index1Name[0] == 'C' ) && ( index1Name[1] == 'M' ) && ( index1Name[2] == 'S' ))
				  &&
				  (( index2Name[0] == 'E' ) || ( index2Name[0] == 'L' )) // EURIBOR or LIBOR
				 )
			  {
				 // don't forget here to swap x and y in that case

				 double tmpVal = x1;

				 x1 = x2;

				 x2 = tmpVal;
			  }
		   }
		   else
		   {
			  intraMktTag = ccy1+"_"+index2Name+"_"+ccy1+"_"+index1Name+"_"+"NOCHECK";
    
			  // because of Summit, the correl matrix is (LIBOR, CMS) and
			  // because of the correl manager ordering :

			  if (!((( index1Name[0] == 'C' ) && ( index1Name[1] == 'M' ) && ( index1Name[2] == 'S' ))
				   &&
				   (( index2Name[0] == 'E' ) || ( index2Name[0] == 'L' )) // EURIBOR or LIBOR
				  )
				 )
			  {
				 // don't forget here to swap x and y in that case

				 double tmpVal = x1;

				 x1 = x2;

				 x2 = tmpVal;
			  }
		   }
		}
		else if ( string(ccy1) < string(ccy2) )
		{
		   intraMktTag = ccy1+"_"+index1Name+"_"+ccy2+"_"+index2Name+"_"+"NOCHECK";
		}
		else
		{
		   intraMktTag = ccy2+"_"+index2Name+"_"+ccy1+"_"+index1Name+"_"+"NOCHECK";

		   // because of Summit the correl matrix is (LIBOR, CMS) and
		   // because of the correl manager ordering :

		   // don't forget here to swap x and y in that case

		   double tmpVal = x1;

		   x1 = x2;

		   x2 = tmpVal;
		}
		
		vCorrel = ComputeCorrelData("IR/IR", intraMktTag, x1, x2);

		return(vCorrel);
	}
	catch(...)
	{

		string	intraMktTag;

		string	ccy1(strccy1);
		string	ccy2(strccy2);

		string	vIndex1Name = index1Name;
		string	vIndex2Name = index2Name;

		if( vIndex1Name.find("CMS") != string::npos )
		{
			vIndex1Name.replace(0, 3, "CO");

			if( vIndex1Name.at(vIndex1Name.size() - 1) != 'Y' )
				vIndex1Name += "Y";
		}
		else if( vIndex1Name.find("LIBOR") != string::npos )
			vIndex1Name.replace(0, 5, "CO");
		else if( vIndex1Name.find("EURIBOR") != string::npos )
			vIndex1Name.replace(0, 7, "CO");

		if( vIndex2Name.find("CMS") != string::npos )
		{
			vIndex2Name.replace(0, 3, "CO");

			if( vIndex2Name.at(vIndex2Name.size() - 1) != 'Y' )
				vIndex2Name += "Y";
		}
		else if( vIndex2Name.find("LIBOR") != string::npos )
			vIndex2Name.replace(0, 5, "CO");
		else if( vIndex2Name.find("EURIBOR") != string::npos )
			vIndex2Name.replace(0, 7, "CO");

		if ( ccy1 == ccy2 )
		{
		   if ( ConvertCOTenorToDouble(vIndex1Name.c_str()) < ConvertCOTenorToDouble(vIndex2Name.c_str()) ) 
		   {
			  intraMktTag = ccy1 + "_" + vIndex1Name + "_" + ccy1 + "_" + vIndex2Name +"_"+"NOCHECK";

			  // because of Summit the correl matrix is (LIBOR, CMS) and
			  // because of the correl manager ordering :
			  if( (index1Name.substr(0, 3) == "CMS") &&
				  (index2Name.substr(0, 3) != "CMS") ) // EURIBOR or LIBOR

			  {
				 double tmpVal = x1;
				 x1 = x2;
				 x2 = tmpVal;
			  }
		   }
		   else
		   {
			  intraMktTag = ccy1 + "_" + vIndex2Name + "_" + ccy1 + "_" + vIndex1Name +"_"+"NOCHECK";

			  // because of Summit, the correl matrix is (LIBOR, CMS) and
			  // because of the correl manager ordering :
			  if( (index1Name.substr(0, 3) == "CMS") &&
				  (index2Name.substr(0, 3) != "CMS") ) // EURIBOR or LIBOR
			  {	
				 double tmpVal = x1;
				 x1 = x2;
				 x2 = tmpVal;
			  }
		   }
		}
		else if ( string(ccy1) < string(ccy2) )
		{
		   intraMktTag = ccy1 + "_" + vIndex1Name + "_" + ccy2 + "_" + vIndex2Name + "_" + "NOCHECK";
		}
		else
		{
		   intraMktTag = ccy2 + "_" + vIndex2Name + "_" + ccy1 + "_" + vIndex1Name + "_" + "NOCHECK";

		   // because of Summit the correl matrix is (LIBOR, CMS) and
		   // because of the correl manager ordering :

		   double tmpVal = x1;
		   x1 = x2;
		   x2 = tmpVal;
		}
		try
		{
			vCorrel = ComputeCorrelData("IR/IR", intraMktTag, x1, x2);
		}
		catch(...)
		{
			int pos = intraMktTag.find("NOCHECK");
			intraMktTag.replace(pos,pos+7,"COMP");
			vCorrel = ComputeCorrelData("IR/IR", intraMktTag, x1, x2);
		}

		return	vCorrel;
	}

}



/*!
 * ComputeCorrelData does a lookup through the market data
 * to get the appropriate Correlation matrix
 *
 */
ARM_CorrelMatrix* ARM_CorrelManager::ComputeCorrelData( const string& tmpMktTag, const string& tmpIntraMktTag ) const
{
	/// upper case conversion
	string mktTag = stringToUpper(tmpMktTag );
	string intraMktTag	= stringToUpper(tmpIntraMktTag );

	//// do the search
	AllMarketCorrelsMap::const_iterator mktIter = itsMktData.find( mktTag );

	/// do we have the correct data per mkt!
	if( mktIter == itsMktData.end() )
	{
		char msg[255];
		sprintf( msg, "Could not find any mkt data for the mkt tag %s, Please advise", mktTag.c_str() );
		throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, msg );
	}
	
	intraMarketCorrelsMap::const_iterator oneMktIter = ( (*mktIter).second).find( intraMktTag );
	if( oneMktIter == ( (*mktIter).second).end() )
	{
		char msg[255];
		sprintf( msg, "Could not find any mkt data for the mkt tag %s, and intraMktTag %s Please advise", mktTag.c_str(), intraMktTag.c_str() );
		throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, msg );
	}

	//// forced to const_cast because we are in a const function
	//// fine here!
	return &( const_cast<ARM_CorrelMatrix&>( ( (*oneMktIter).second) ) );
}



/*!
 * View method for inf cap and floor
 */
void ARM_CorrelManager::View(char* id, FILE* ficOut)
{
    FILE* fOut;
    char fOutName[200];

    if ( ficOut == NULL )
    {
       ARM_GetViewFile(ARM_VIEW_FILE, id, fOutName);
       (void) unlink(fOutName);
       fOut = fopen(fOutName, "w");
    }
    else
       fOut = ficOut;

    fprintf(fOut, "\n       >>>>>>>>>>> Correlation Manager <<<<<<<<<<<\n");

	AllMarketCorrelsMap::const_iterator iter	= itsMktData.begin();
	AllMarketCorrelsMap::const_iterator iterEnd	= itsMktData.end();

	/// loop over markets
	for( ; iter != iterEnd; ++iter )
	{
		fprintf(fOut, "Markets:\t%s\n", ((*iter).first).c_str() );

		intraMarketCorrelsMap::const_iterator oneMktDataIter		= ((*iter).second).begin();
		intraMarketCorrelsMap::const_iterator oneMktDataIterEnd		= ((*iter).second).end();

		/// loop over indexes for a given market
		for( ; oneMktDataIter != oneMktDataIterEnd; ++oneMktDataIter )
		{
			fprintf(fOut, "\nIndex tag:\t%s\n", ((*oneMktDataIter).first).c_str() );

			///view of the correlation matrix
			( const_cast<ARM_CorrelMatrix&>( ( (*oneMktDataIter).second) ) ).View( id, fOut );
			fprintf(fOut, "\n" );
		}
	}

    if ( ficOut == NULL )
       fclose(fOut);
}

/*!
 * validation for the data at its own level and in the CorrelManager's level
 *
 * 1) regex on the mktTag and intraMktTag. 
 *  mktTag has to be of the form mkt1/mkt2
 *  with mkt1 < mkt2 (sorted in alphabetical order)
 *	intraMktTag should be of the form
 *		- ccy1_Tenor1_Ccy2_DIAG
 *		- ccy1_Tenor1_Ccy2_Tenor2_COMP
 *		- ccy1_Tenor1_Ccy2_Tenor2_NOCHECK
 *		- ccy1/ind1_ccy2/ind2
 *
 * 2) each correl data has to be between -1 and 1 (up to the resize factor which is defined in ARM_Constants::correlBase
 * 3) matrix in complet mode has to be 
 *			- general
 *				- squared
 *				- same X and Y
 *			- if using the same subtags
 *				- symetric 
 *				- with the diagonal with 1.0 (up to the resize factor which is defined in ARM_Constants::correlBase )
 *
 */

 ARM_CorrelManager::MatType ARM_CorrelManager::ValidateCorrelMatrix( const string& mktTag, const string& intraMktTag, const ARM_CorrelMatrix& correlMatrix ) const
{
	/// 1) regex on the mktTag and intraMktTag. 
	///  mktTag has to be of the form mkt1/mkt2
	///  with mkt1 < mkt2 (sorted in alphabetical order)
	///	intraMktTag should be of the form
	///		- ccy1_Tenor1_Ccy2_DIAG
	///		- ccy1_Tenor1_Ccy2_Tenor2_COMP
	///		- ccy1/ind1_ccy2/ind2

	/// case  of the mktTag
	vector< string > subTags = splitStringIntoPieces( mktTag, "/" );
	if( subTags.size() != 2 )
		throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, "mktTag should be of the form mkt1/mkt2" );
	if( subTags[ 0 ] > subTags[ 1 ] )
		throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, 
			"mktTag should be of the form mkt1/mkt2 with mkt1 < mkt2 according to alphabetical order" );

	/// case of the intraMktTag
	subTags = splitStringIntoPieces( intraMktTag, "_" );
	bool fine = false;

	MatType validationType = NOCHECK;
	
	///	intraMktTag should be of the form
	///		- ccy1_Tenor1_Ccy2_DIAG
	///		- ccy1_Tenor1_Ccy2_Tenor2_COMP
	///		- ccy1/ind1_ccy2/ind2 sorted in alphabetical order
	switch( subTags.size() )
	{
	case 2:	
		/// checks that it is correctly sorted
		fine = subTags[ 0 ] <= subTags[ 1 ] ;
		break;
	case 4:
		{
			fine = ( subTags[3] == "DIAG" );
			validationType = DIAG;
			break;
		}
	case 5:
		{
            bool comp_nocheck_fine = ( subTags[4] == "COMP" ||  subTags[4] == "NOCHECK");
            if(!comp_nocheck_fine)
                throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"intraMktTag should be ended by  sorted COMP or NOCHECK" );
            
            bool tenor_fine     = ConvertCOTenorToDouble(subTags[1].c_str()) <= ConvertCOTenorToDouble(subTags[3].c_str()) ? true : false;
            bool equal_ccy_fine = (subTags[0] == subTags[2]);
            if(comp_nocheck_fine && !tenor_fine)
                    throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"intraMktTag's tenors should be sorted" );
            else if (subTags[0] > subTags[2])
                throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"intraMktTag's currencies should be sorted" );
            else
                fine = true;

			validationType = subTags[4] == "COMP" ? COMP : NOCHECK;
			break;
		}
	default:
		fine = false;
	}
	if( !fine )
		throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, 
			"intraMktTag should be ccy1_Tenor1_Ccy2_DIAG, or ccy1_Tenor1_Ccy2_Tenor2_COMP ccy1_Tenor1_Ccy2_Tenor2_NOCHECK or ccy1|ind1_ccy2|ind2 sorted in alphabetical order" );

    /// 2) each correl data has to be between -1 and 1 (up to the resize factor which is defined in ARM_Constants::correlBase )
	ARM_Matrix* data= correlMatrix.GetCorrelData();

	unsigned int i, j;
	unsigned int imax = data->GetNumLines();
	unsigned int jmax = data->GetNumCols();

	/// some validation on size to avoid infinite loop
	unsigned int MAXSIZE = 10000;
	if( imax > MAXSIZE )
	{
		// char msg[255];
		// sprintf( msg, "the row size found for the correlation matrix is %d which is bigger than the max size authorised %d! Please advise", imax, MAXSIZE);
		// throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, msg );
		ARMTHROW(ERR_INVALID_ARGUMENT,"the row size found for the correlation matrix is "<<imax
			<<" which is bigger than the max size authorised "<<MAXSIZE
			<<" ! Please advise") ; 
	}

	if( jmax > MAXSIZE )
	{
		// char msg[255];
		// sprintf( msg, "the column size found for the correlation matrix is %d which is bigger than the max size authorised %d! Please advise", jmax, MAXSIZE);
		// throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, msg );
		ARMTHROW(ERR_INVALID_ARGUMENT,"the column size found for the correlation matrix is "
			<<jmax<<" which is bigger than the max size authorised "
			<<MAXSIZE<<"! Please advise"); 
	}

	for( i=0; i<imax; ++i )
		for( j=0; j<jmax; ++j )
		{
			if( data->Elt(i,j) > ARM_Constants::correlBase || data->Elt(i,j) < - ARM_Constants::correlBase )
			{
				// char msg[255];
				// sprintf( msg, "the data %s for row %d and column %d is not consistent, Correlation has to be between -%d and %d up to the rescaling factor! Please advise", 
				// 	ARM_Constants::correlBase, ARM_Constants::correlBase, i, j );
				// throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, msg );
				ARMTHROW(ERR_INVALID_ARGUMENT,"the data "<<ARM_Constants::correlBase
					<<" for row " <<i <<"  and column "<< j<<" is not consistent, Correlation has to be between "
					<<"-"<<ARM_Constants::correlBase <<"and "<<ARM_Constants::correlBase <<"up to the rescaling factor! Please advise"); 
			}
		}
	
    /// 3) matrix in complet mode has to be 
	///			- general
	///				- squared
	///				- same X and Y
	///			- if using the same subtags
	///					- symetric 
	///					- with the diagonal with 1.0 (up to the resize factor which is defined in ARM_Constants::correlBase )
	if( validationType == COMP )
	{
		if( imax != jmax )
			throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, "matrix has to be squared!, please advise" );

		ARM_Vector* X = correlMatrix.GetX();
		ARM_Vector* Y = correlMatrix.GetY();

		/// squared with the same tenor!
		for( i=0; i<imax; ++i )
			if( fabs( X->Elt(i) - Y->Elt(i) ) > PRECISION )
			{
				char msg[255];
				sprintf( msg, "the tenor for X and Y should be the same.. apparently not for i = j = %d value for i %f for j %f",
					i, X->Elt(i), Y->Elt(i) );
				throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, msg );
			}

			
		/// tests that the subtype are the same
		if( ( subTags[ 0 ] + "_" + subTags[ 1 ] ) == ( subTags[ 2 ] + "_" + subTags[ 3 ] ) )
		{
			for( i=0; i<imax; ++i )
			{
				if( fabs( data->Elt(i,i ) - ARM_Constants::correlBase ) > PRECISION )
				{
					char msg[255];
					sprintf( msg, "diagonal terms should be %.0f, this is not the case for index %d with value %.3f", 
						ARM_Constants::correlBase, i, data->Elt(i,i ) );
					throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, msg);
				}

				for( j=0; j<i; ++j )
				{
					if( fabs( data->Elt(i,j ) - data->Elt(j ,i ) ) > PRECISION )
					{
						char msg[255];
						sprintf( msg, "the matrix has to be symmetric, this is not the case for index i= %d and j= %d, value %.3f compared to its symmetric one %.3f", 
							i, j, data->Elt(i,j ), data->Elt(j,i ) );
						throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, msg);
					}
				}
			}
		}
	}


    return validationType;
	
}

/*!
 * validation for the data at its own level and in the CorrelManager's level
 *
 * preliminary checks that keys are upper case
 * 1) checks CorrelMatrix data
 * 2) if the tag already exists, throw an exception...does not allow overwritting
 * 3) mod DIAG has to be consistent with the corresponding COMP mod!
 * 4) mod COMP has to be consistent with the corresponding DIAG mod!
 *
 *
 */
void ARM_CorrelManager::ValidateData( const string& mktTag, const string& intraMktTag, const ARM_CorrelMatrix& correlMatrix ) const
{
	/// preliminary checks that keys are upper case
	StringUpperCaseValidateWith( mktTag, "mkTag" );
	StringUpperCaseValidateWith( intraMktTag, "intraMkTag" );

    vector< string > subTags = splitStringIntoPieces( intraMktTag, "_" );
    /// 1) checks CorrelMatrix data
    MatType validationType = ValidateCorrelMatrix(mktTag,intraMktTag,correlMatrix);

	/// 2) checks that the tag does not already exist!
	intraMarketCorrelsMap intraMarketCorrels = GetIntraMarketCorrels( mktTag );
	if( (intraMarketCorrels.find( intraMktTag ) != intraMarketCorrels.end() ) && (intraMarketCorrels.size()>1))
	{
		char msg[255];
		sprintf( msg, "the mkt data %s, with intraMktTag %s already exists! Please advise", mktTag.c_str(), intraMktTag.c_str() );
		throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, msg );
	}
	
   	
	ARM_Matrix* data= correlMatrix.GetCorrelData();
	unsigned int i, j;
	unsigned int imax = data->GetNumLines();
	unsigned int jmax = data->GetNumCols();

    /// 3) mod DIAG has to be consistent with the corresponding COMP mod!
	if( validationType == DIAG )
	{
		ARM_Vector* X = correlMatrix.GetX();
		ARM_Vector* Y = correlMatrix.GetY();
		imax = X->GetSize();
		
		for( i=0; i<imax; ++i )
		{
			string subTag1	= subTags[0] + "_" + subTags[1];
			string tenor(YearTermToStringMatu( X->Elt(i)));
			string subTag2	= subTags[2] + "_" + string( tenor );
			string tag = ( subTag1 < subTag2 ? subTag1 + "_" + subTag2 : subTag2 + "_" + subTag1 ) 
				+ "_" + "COMP";
			
			ARM_CorrelMatrix* compCorrelMatrix = GetIntraMarketCorrel( mktTag, tag );
			if( compCorrelMatrix )
			{
				for( j = 0; j<jmax; ++j )
				{
					double ourValue  = correlMatrix.ComputeCorrelData( X->Elt(i), Y->Elt(j) );
					double compValue = compCorrelMatrix->ComputeCorrelData( Y->Elt(j),Y->Elt(j) );
					if( fabs( compValue - ourValue ) > PRECISION )
					{
						// char msg[255];
						// sprintf( msg, "unconsistency between COMP Matrix and DIAG one, found that for intraMktTag %s for tenor %s for fixing %f found value %f as opposed to the one in comp %f", tag, tenor, Y->Elt(j), ourValue, compValue );
						// throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, msg);
						ARMTHROW(ERR_INVALID_ARGUMENT,"unconsistency between COMP Matrix and DIAG one, found that for intraMktTag "
							<< tag <<" for tenor "<< tenor <<" for fixing "<< Y->Elt(j) 
							<<" found value "<<ourValue <<" as opposed to the one in comp "<<compValue); 
					}
				}
			}
		}
	}
	
	/// 4) mod COMP has to be consistent with the corresponding DIAG mod!
	if( validationType == COMP )
	{

		ARM_Vector* X = correlMatrix.GetX();
		ARM_Vector* Y = correlMatrix.GetY();
		imax = X->GetSize();

		/// search for the two corresponding DIAG matrix
		vector< string > TagsToLook;
		vector< string > TenorsToLook;

		TagsToLook.push_back( subTags[0] + "_" + subTags[1] + "_" + subTags[2] + "_" + string("DIAG") );
		TenorsToLook.push_back( subTags[3] );

		string subTag = subTags[2] + "_" + subTags[3] + "_" + subTags[0] + "_" + string("DIAG");
		if( subTag != TagsToLook[0] )
		{
			TagsToLook.push_back( subTag );
			TenorsToLook.push_back( subTags[1] );
		}

		std::vector< string >::const_iterator tagIter	= TagsToLook.begin();
		std::vector< string >::const_iterator tagEnd	= TagsToLook.end();
		std::vector< string >::const_iterator tenorIter	= TenorsToLook.begin();


		for( ; tagIter!= tagEnd; ++tagIter, ++tenorIter )
		{
			ARM_CorrelMatrix* compCorrelMatrix = GetIntraMarketCorrel( mktTag, *tagIter );
			if( compCorrelMatrix )
			{
				char* tenor = (char*) (*tenorIter).c_str();
				double tenorDble = StringMatuToYearTerm( tenor );
				delete tenor;
				for( i=0; i <imax; ++i )
				{
					double ourValue  = correlMatrix.ComputeCorrelData( X->Elt(i), Y->Elt(i) );
					double compValue = compCorrelMatrix->ComputeCorrelData( X->Elt(i), tenorDble );
					if( fabs( compValue - ourValue ) > PRECISION )
					{
						// char msg[255];
						// sprintf( msg, "unconsistency between COMP Matrix and DIAG one, found that for intraMktTag %s for tenor %s for fixing %f found value %f as opposed to the one in comp %f", *tagIter, *tenorIter, Y->Elt(j), ourValue, compValue );
						// throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, msg);
						ARMTHROW(ERR_INVALID_ARGUMENT,
							"unconsistency between COMP Matrix and DIAG one, found that for intraMktTag "<<*tagIter
							<<"  for tenor "<< *tenorIter<<" for fixing "<< Y->Elt(j)
							<<" found value "<< ourValue<<" as opposed to the one in comp"
							<< compValue);  
					}
				}
			}
		}
	}
};


void ARM_CorrelManager::BumpVolatility( const vector<string> mktTag, 
										const vector<string> intraMktTag, 
										double value, 
										int nthLine, 
										int nthCol,
										int isCumul, 
										int isAbsolute )
{
	for (int k=0; k<mktTag.size(); k++)
	{
		ARM_CorrelMatrix* correlMatrix = GetIntraMarketCorrel(mktTag[k], intraMktTag[k]);

		if (correlMatrix == NULL)
		{
			char msg[255];
			sprintf(msg, "tags do not correspond to a curve of Correl Manager\n");
			throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, msg);
		}

		correlMatrix->BumpVolatility(value, nthLine, nthCol, isCumul, isAbsolute);

		//validate values
		for (int i=0; i<correlMatrix->GetExpiryTerms()->GetSize(); i++)
		{
			for (int j=0; j<correlMatrix->GetStrikes()->GetSize(); j++)
			{
				if (correlMatrix->GetVolatilities()->Elt(i, j) > 100.0)
					correlMatrix->GetVolatilities()->Elt(i, j) = 100.0;

				if (correlMatrix->GetVolatilities()->Elt(i, j) < -100.0)
					correlMatrix->GetVolatilities()->Elt(i, j) = -100.0;
			}
		}
	}
}

ARM_CorrelMatrix* ARM_CorrelManager::GetCorrelData( const string& tmpMktTag, const string& tmpIntraMktTag ) const
{
	/// upper case conversion
	string mktTag = stringToUpper(tmpMktTag );
	string intraMktTag	= stringToUpper(tmpIntraMktTag );

	//// do the search
	AllMarketCorrelsMap::const_iterator mktIter = itsMktData.find( mktTag );

	/// do we have the correct data per mkt!
	if( mktIter == itsMktData.end() )
		return NULL;
	
	/// if so goahead with the intraMktTag
	intraMarketCorrelsMap::const_iterator oneMktIter = ( (*mktIter).second).find( intraMktTag );

	/// test that we could find something
	if( oneMktIter == ( (*mktIter).second).end() )
		return NULL;

	/// ok there is really something
	/// a little bit ugly but required to const_cast!
	return &( const_cast<ARM_CorrelMatrix&>( (*oneMktIter).second ) );
}

//add public 


/*******************************************************************
///////////////////////////////////////////////////////////////////*
///////////////////////////////////////////////////////////////////*
///////////////////////////////////////////////////////////////////*
///////////////////////////////////////////////////////////////////*
///////////// ARM_CorrelatorManager PART	///////////////////////*
///////////////////////////////////////////////////////////////////*
///////////////////////////////////////////////////////////////////*
///////////////////////////////////////////////////////////////////*
/*******************************************************************
*/
ARM_CorrelatorManager::ARM_CorrelatorManager()
{
	Init();
};

ARM_CorrelatorManager::~ARM_CorrelatorManager()
{
	CleanUp();
};


ARM_CorrelatorManager& ARM_CorrelatorManager::operator= ( const ARM_CorrelatorManager& src )
{
	if( this != &src )
	{
		CleanUp();	
		ARM_Object::operator = ( src );
		BitwiseCopy( &src );
	}

	return *this;
};

ARM_CorrelatorManager::ARM_CorrelatorManager( const ARM_CorrelatorManager& src )
//:ARM_Object( src )
{
	Copy( &src );
};

void ARM_CorrelatorManager::BitwiseCopy(const ARM_Object* src)
{

	ARM_CorrelatorManager* CorrSrc = dynamic_cast< ARM_CorrelatorManager* >( const_cast< ARM_Object* >( src ) );

	if( CorrSrc )
	{

		map< string , ARM_IndexIndexCorrelCube* >::iterator vIndexIterator = CorrSrc->itsIndexData.begin();
		for( ; vIndexIterator != CorrSrc->itsIndexData.end() ; vIndexIterator++ )
			itsIndexData[vIndexIterator->first] = ( ARM_IndexIndexCorrelCube* )(vIndexIterator->second)->Clone();

		map< string , ARM_HyperCube* >::iterator vDiagIterator = CorrSrc->itsDiagData.begin();
		for( ; vDiagIterator != CorrSrc->itsDiagData.end() ; vDiagIterator++ )
			itsDiagData[vDiagIterator->first] = ( ARM_HyperCube* )(vDiagIterator->second)->Clone();

		map< string , ARM_VolCurve*>::iterator vCorrIterator = CorrSrc->itsCorrData.begin();
		for( ; vCorrIterator != CorrSrc->itsCorrData.end() ; vCorrIterator++ )
			itsCorrData[vCorrIterator->first] = ( ARM_VolCurve* )(vCorrIterator->second)->Clone();

		map< string , ARM_VolCurve*>::iterator vMarketIterator = CorrSrc->itsSimpleModeData.begin();
		for( ; vMarketIterator != CorrSrc->itsSimpleModeData.end() ; vMarketIterator++ )
			itsSimpleModeData[vMarketIterator->first] = ( ARM_VolCurve* )(vMarketIterator->second)->Clone();

		vDiagIterator = CorrSrc->itsIRVolData.begin();
		for( ; vDiagIterator != CorrSrc->itsIRVolData.end() ; vDiagIterator++ )
			itsIRVolData[vDiagIterator->first] = ( ARM_HyperCube* )(vDiagIterator->second)->Clone();

		vDiagIterator = CorrSrc->itsVolVolData.begin();
		for( ; vDiagIterator != CorrSrc->itsVolVolData.end() ; vDiagIterator++ )
			itsVolVolData[vDiagIterator->first] = ( ARM_HyperCube* )(vDiagIterator->second)->Clone();

		map< string , ARM_VolCube*>::iterator vFXIterator = CorrSrc->itsVolFXData.begin();
		for( ; vFXIterator != CorrSrc->itsVolFXData.end() ; vFXIterator++ )
			itsVolFXData[vFXIterator->first] = ( ARM_VolCube* )(vFXIterator->second)->Clone();
	}

};


ARM_Object* ARM_CorrelatorManager::Clone()
{
	ARM_CorrelatorManager* theClone = new ARM_CorrelatorManager();

	theClone->Copy(this);

	return theClone;
};

void ARM_CorrelatorManager::Copy(const ARM_Object* src)
{
	ARM_Object::Copy(src);

	BitwiseCopy(src);
	
};

void ARM_CorrelatorManager::CleanUp()
{
	try
	{
		if( !itsIndexData.empty() )
		{
			map< string , ARM_IndexIndexCorrelCube* >::iterator vIndexIterator = itsIndexData.begin();
			for( ; vIndexIterator != itsIndexData.end() ; vIndexIterator++ )
			{	
				if( vIndexIterator->second )
				{
					delete vIndexIterator->second;
					vIndexIterator->second = NULL;
				}
			}
		}

		if( !itsDiagData.empty() )
		{
			map< string , ARM_HyperCube* >::iterator vDiagIterator = itsDiagData.begin();
			for( ; vDiagIterator != itsDiagData.end() ; vDiagIterator++ )
			{	
				if( vDiagIterator->second )
				{
					delete vDiagIterator->second;
					vDiagIterator->second = NULL;
				}
			}
		}

		if( !itsCorrData.empty() )
		{
			map< string , ARM_VolCurve*>::iterator vMarketIterator = itsCorrData.begin();
			for( ; vMarketIterator != itsCorrData.end() ; vMarketIterator++ )
			{	
				if( vMarketIterator->second )
				{
					delete vMarketIterator->second;
					vMarketIterator->second = NULL;
				}
			}
		}

		if( !itsSimpleModeData.empty() )
		{
			map< string , ARM_VolCurve*>::iterator vMarketIterator = itsSimpleModeData.begin();
			for( ; vMarketIterator != itsSimpleModeData.end() ; vMarketIterator++ )
			{	
				if( vMarketIterator->second )
				{
					delete vMarketIterator->second;
					vMarketIterator->second = NULL;
				}
			}
		}

		if( !itsIRVolData.empty() )
		{
			map< string , ARM_HyperCube* >::iterator vDiagIterator = itsIRVolData.begin();
			for( ; vDiagIterator != itsIRVolData.end() ; vDiagIterator++ )
			{	
				if( vDiagIterator->second )
				{
					delete vDiagIterator->second;
					vDiagIterator->second = NULL;
				}
			}
		}

		if( !itsVolVolData.empty() )
		{
			map< string , ARM_HyperCube* >::iterator vDiagIterator = itsVolVolData.begin();
			for( ; vDiagIterator != itsVolVolData.end() ; vDiagIterator++ )
			{	
				if( vDiagIterator->second )
				{
					delete vDiagIterator->second;
					vDiagIterator->second = NULL;
				}
			}
		}

		if( !itsVolFXData.empty() )
		{
			map< string , ARM_VolCube* >::iterator vDiagIterator = itsVolFXData.begin();
			for( ; vDiagIterator != itsVolFXData.end() ; vDiagIterator++ )
			{	
				if( vDiagIterator->second )
				{
					delete vDiagIterator->second;
					vDiagIterator->second = NULL;
				}
			}
		}
	}
	catch(...)
	{
		throw Exception(__LINE__, __FILE__, ERR_INVALID_DATA,  
						"error in ARM_CorrelatorManager::CleanUp()");
	}
};

void ARM_CorrelatorManager::Init()
{
    SetName(ARM_CORRELMANAGER);
};

void ARM_CorrelatorManager::Add( const string& name, const ARM_Object* src, MapType mtype)
{
	string theName = stringToUpper( name );
	//ARM_Object* theObject = ( ARM_Object* )src;

	switch( mtype )
	{
		case IndexIndexCorr :
		{
			ARM_IndexIndexCorrelCube* dest = ( ARM_IndexIndexCorrelCube* )src ;
			
			map< string , ARM_IndexIndexCorrelCube* >::iterator iter = itsIndexData.find(theName);
			if( iter != itsIndexData.end() )
				delete iter->second;
			itsIndexData[theName] = ( ARM_IndexIndexCorrelCube* )dest->Clone();

			break;
		}

		case DiagCorr :
		{
			ARM_HyperCube* dest = ( ARM_HyperCube* )src ;
			
			map< string , ARM_HyperCube* >::iterator iter = itsDiagData.find(theName);
			if( iter != itsDiagData.end() )
				delete iter->second;
			itsDiagData[theName] = ( ARM_HyperCube* )dest->Clone();

			break;
		}

		case Corr:
		{
			ARM_VolCurve* dest = ( ARM_VolCurve* )src;
			
			map< string , ARM_VolCurve*>::iterator iter = itsCorrData.find(theName);
			if( iter != itsCorrData.end() )
				delete iter->second;
			itsCorrData[theName] = ( ARM_VolCurve* )dest->Clone();

			break;
		}

		case IRVolCorr :
		{
			ARM_HyperCube* dest = ( ARM_HyperCube* )src ;
			
			map< string , ARM_HyperCube* >::iterator iter = itsIRVolData.find(theName);
			if( iter != itsIRVolData.end() )
				delete iter->second;
			itsIRVolData[theName] = ( ARM_HyperCube* )dest->Clone();

			break;
		}

		case VolVolCorr :
		{
			ARM_HyperCube* dest = ( ARM_HyperCube* )src ;
			
			map< string , ARM_HyperCube* >::iterator iter = itsVolVolData.find(theName);
			if( iter != itsVolVolData.end() )
				delete iter->second;
			itsVolVolData[theName] = ( ARM_HyperCube* )dest->Clone();

			break;
		}

		case FXVolCorr :
		{
			ARM_VolCube* dest = ( ARM_VolCube* )src ;
			
			map< string , ARM_VolCube* >::iterator iter = itsVolFXData.find(theName);
			if( iter != itsVolFXData.end() )
				delete iter->second;
			itsVolFXData[theName] = ( ARM_VolCube* )dest->Clone();

			break;
		}

		default:
		{
			ARM_VolCurve* dest = ( ARM_VolCurve* )src;
			
			map< string , ARM_VolCurve*>::iterator iter = itsSimpleModeData.find(theName);
			if( iter != itsSimpleModeData.end() )
				delete iter->second;
			itsSimpleModeData[theName] = ( ARM_VolCurve* )dest->Clone();

			break;
		}

	}
};

ARM_IndexIndexCorrelCube* ARM_CorrelatorManager::GetCorrelIndexIndexCube( const string ccy ) const
{
	map< string , ARM_IndexIndexCorrelCube* >::const_iterator iter = itsIndexData.find( stringToUpper( ccy ) );
	if( iter == itsIndexData.end() )
	{
		char vErrorMsg[50];

		sprintf(vErrorMsg, "Couldn't find %s IndexIndex Correl cube", ccy.c_str());

		throw Exception(__LINE__, __FILE__, ERR_INVALID_DATA, vErrorMsg, 
						"ARM_CorrelatorManager::GetCorrelIndexIndexCube(...): Correl not found");
	}

	return iter->second;
	
};
	
ARM_HyperCube* ARM_CorrelatorManager::GetCorrelHyperCube( const string ccy ) const
{
	map< string , ARM_HyperCube* >::const_iterator iter =  itsDiagData.find( stringToUpper( ccy ) );
	if( iter == itsDiagData.end() )
	{
		char vErrorMsg[50];

		sprintf(vErrorMsg, "Couldn't find %s HyperCube correl", ccy.c_str());

		throw Exception(__LINE__, __FILE__, ERR_INVALID_DATA, vErrorMsg, 
						"ARM_CorrelatorManager::GetCorrelHyperCube(...): Correl not found");
	}

	return iter->second;
};

ARM_VolCurve* ARM_CorrelatorManager::GetCorrModeCorrelCurve( const string name ) const
{
	map< string , ARM_VolCurve* >::const_iterator iter = itsCorrData.find( stringToUpper( name ) );
	if( iter == itsCorrData.end() )
	{
		char vErrorMsg[50];

		sprintf(vErrorMsg, "Couldn't find %s VolCurve ", name.c_str());

		throw Exception(__LINE__, __FILE__, ERR_INVALID_DATA, vErrorMsg, 
						"ARM_CorrelatorManager::GetCorrModeCorrelCurve(...): Correl not found");
	}

	return iter->second;
};

ARM_VolCurve* ARM_CorrelatorManager::GetSimpleModeCorrelCurve( const string name ) const
{
	map< string , ARM_VolCurve* >::const_iterator iter = itsSimpleModeData.find( stringToUpper( name ) );
	if( iter == itsSimpleModeData.end() )
	{
		char vErrorMsg[50];

		sprintf(vErrorMsg, "Couldn't find %s VolCurve ", name.c_str());

		throw Exception(__LINE__, __FILE__, ERR_INVALID_DATA, vErrorMsg, 
						"ARM_CorrelatorManager::GetSimpleModeCorrelCurve(...): Correl not found");
	}

	return iter->second;
};


ARM_HyperCube* ARM_CorrelatorManager::GetIRVolCorrelHyperCube( const string ccy ) const
{
	map< string , ARM_HyperCube* >::const_iterator iter =  itsIRVolData.find( stringToUpper( ccy ) );
	if( iter == itsIRVolData.end() )
	{
		char vErrorMsg[50];

		sprintf(vErrorMsg, "Couldn't find IR/Vol Correl for %s", ccy.c_str());

		throw Exception(__LINE__, __FILE__, ERR_INVALID_DATA, vErrorMsg, 
						"ARM_CorrelatorManager::GetIRVolCorrelHyperCube(...) : Correl not found");
	}

	return iter->second;
}


ARM_HyperCube* ARM_CorrelatorManager::GetVolVolCorrelHyperCube( const string ccy ) const
{
	map< string , ARM_HyperCube* >::const_iterator iter =  itsVolVolData.find( stringToUpper( ccy ) );
	if( iter == itsVolVolData.end() )
	{
		char vErrorMsg[50];

		sprintf(vErrorMsg, "Couldn't find Vol/Vol Correl for %s", ccy.c_str());

		throw Exception(__LINE__, __FILE__, ERR_INVALID_DATA, vErrorMsg, 
						"ARM_CorrelatorManager::GetVolVolCorrelHyperCube(...) : Correl not found");
	}

	return iter->second;
}

ARM_VolCube* ARM_CorrelatorManager::GetFXCorrelCube( const string ccy ) const
{
	map< string , ARM_VolCube* >::const_iterator iter =  itsVolFXData.find( stringToUpper( ccy ) );
	if( iter == itsVolFXData.end() )
	{
		char vErrorMsg[50];

		sprintf(vErrorMsg, "Couldn't find FX Correl for %s", ccy.c_str());

		throw Exception(__LINE__, __FILE__, ERR_INVALID_DATA, vErrorMsg, 
						"ARM_CorrelatorManager::GetFXCorrelHyperCube(...) : Correl not found");
	}

	return iter->second;
}

void ARM_CorrelatorManager::View(char* id, FILE* ficOut)
{
    FILE* fOut;
    char fOutName[200];

    if ( ficOut == NULL )
    {
       ARM_GetViewFile(ARM_VIEW_FILE, id, fOutName);
       (void) unlink(fOutName);
       fOut = fopen(fOutName, "w");
    }
    else
       fOut = ficOut;

    fprintf(fOut, "\n       >>>>>>>>>>> Correlation Manager <<<<<<<<<<<\n");
	

	if ( ! itsDiagData.empty() )
	{
		map< string , ARM_HyperCube* >::iterator iterDiag = itsDiagData.begin();
		map< string , ARM_HyperCube* >::iterator iterDiagEnd = itsDiagData.end();
		fprintf(fOut, "\n  xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx \n");
		fprintf(fOut, "\n\n >>>>>>>>>>>>>>>>>>>>>> Begin HyperCube Correl Map <<<<<<<<<<<<<<<<<<<<<<<\n\n");
	
		for( ; iterDiag != iterDiagEnd; ++iterDiag )
		{
			fprintf(fOut, "Currency:\t%s\n", ((*iterDiag).first).c_str() );
			(*iterDiag).second->View( id, fOut );

		}

		fprintf(fOut, "\n\n <<<<<<<<<<<<<<<<<<<<<< End HyperCube Correl Map >>>>>>>>>>>>>>>>>>>>>>>\n\n");
		fprintf(fOut, "\n xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx \n");
	}

	if ( ! itsIndexData.empty() )
	{
		map< string , ARM_IndexIndexCorrelCube* >::iterator iterIndex = itsIndexData.begin();
		map< string , ARM_IndexIndexCorrelCube* >::iterator iterIndexEnd = itsIndexData.end();
		fprintf(fOut, "\n  xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx \n");
		fprintf(fOut, "\n\n >>>>>>>>>>>>>>>>>>>>>> Begin Index-Index Correl Map <<<<<<<<<<<<<<<<<<<<<<<\n\n");
	
		for( ; iterIndex != iterIndexEnd; ++iterIndex )
		{
			fprintf(fOut, "Currency:\t%s\n", ((*iterIndex).first).c_str() );
			(*iterIndex).second->View( id, fOut );

		}

		fprintf(fOut, "\n\n <<<<<<<<<<<<<<<<<<<<<< End Index-Index Correl Map >>>>>>>>>>>>>>>>>>>>>>>\n\n");
		fprintf(fOut, "\n xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx \n");
	}

	if ( ! itsCorrData.empty() )
	{
		map< string , ARM_VolCurve* >::iterator iterMarket = itsCorrData.begin();
		map< string , ARM_VolCurve* >::iterator iterMarketEnd = itsCorrData.end();
		fprintf(fOut, "\n  xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx \n");
		fprintf(fOut, "\n\n >>>>>>>>>>>>>>>>>>>>>> Begin VolCurve Correl Map <<<<<<<<<<<<<<<<<<<<<<<\n\n");

		for( ; iterMarket != iterMarketEnd; ++iterMarket )
		{
			fprintf(fOut, "Currency:\t%s\n", ((*iterMarket).first).c_str() );
			(*iterMarket).second->View( id, fOut );

		}

		fprintf(fOut, "\n\n <<<<<<<<<<<<<<<<<<<<<< End VolCurve Correl Map >>>>>>>>>>>>>>>>>>>>>>>\n\n");
		fprintf(fOut, "\n xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx \n");
	}

	if ( ! itsSimpleModeData.empty() )
	{
		map< string , ARM_VolCurve* >::iterator iterMarket = itsSimpleModeData.begin();
		map< string , ARM_VolCurve* >::iterator iterMarketEnd = itsSimpleModeData.end();
		fprintf(fOut, "\n  xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx \n");
		fprintf(fOut, "\n\n >>>>>>>>>>>>>>>>>>>>>> Begin VolCurve Index Map <<<<<<<<<<<<<<<<<<<<<<<\n\n");

		for( ; iterMarket != iterMarketEnd; ++iterMarket )
		{
			fprintf(fOut, "Currency:\t%s\n", ((*iterMarket).first).c_str() );
			(*iterMarket).second->View( id, fOut );

		}

		fprintf(fOut, "\n\n <<<<<<<<<<<<<<<<<<<<<< End VolCurve Index Map >>>>>>>>>>>>>>>>>>>>>>>\n\n");
		fprintf(fOut, "\n xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx \n");
	}

	if ( ! itsIRVolData.empty() )
	{
		map< string , ARM_HyperCube* >::iterator iterDiag = itsIRVolData.begin();
		map< string , ARM_HyperCube* >::iterator iterDiagEnd = itsIRVolData.end();
		fprintf(fOut, "\n  xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx \n");
		fprintf(fOut, "\n\n >>>>>>>>>>>>>>>>>>>>>> Begin IR/VOL Correl Map <<<<<<<<<<<<<<<<<<<<<<<\n\n");
	
		for( ; iterDiag != iterDiagEnd; ++iterDiag )
		{
			fprintf(fOut, "Currency:\t%s\n", ((*iterDiag).first).c_str() );
			(*iterDiag).second->View( id, fOut );

		}

		fprintf(fOut, "\n\n <<<<<<<<<<<<<<<<<<<<<< End IR/VOL Correl Map >>>>>>>>>>>>>>>>>>>>>>>\n\n");
		fprintf(fOut, "\n xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx \n");
	}

	if ( ! itsVolVolData.empty() )
	{
		map< string , ARM_HyperCube* >::iterator iterDiag = itsVolVolData.begin();
		map< string , ARM_HyperCube* >::iterator iterDiagEnd = itsVolVolData.end();
		fprintf(fOut, "\n  xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx \n");
		fprintf(fOut, "\n\n >>>>>>>>>>>>>>>>>>>>>> Begin VOL/VOL Correl Map <<<<<<<<<<<<<<<<<<<<<<<\n\n");
	
		for( ; iterDiag != iterDiagEnd; ++iterDiag )
		{
			fprintf(fOut, "Currency:\t%s\n", ((*iterDiag).first).c_str() );
			(*iterDiag).second->View( id, fOut );

		}

		fprintf(fOut, "\n\n <<<<<<<<<<<<<<<<<<<<<< End VOL/VOL Correl Map >>>>>>>>>>>>>>>>>>>>>>>\n\n");
		fprintf(fOut, "\n xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx \n");
	}

	if ( ! itsVolFXData.empty() )
	{
		map< string , ARM_VolCube* >::iterator iterDiag = itsVolFXData.begin();
		map< string , ARM_VolCube* >::iterator iterDiagEnd = itsVolFXData.end();
		fprintf(fOut, "\n  xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx \n");
		fprintf(fOut, "\n\n >>>>>>>>>>>>>>>>>>>>>> Begin FX Correl Map <<<<<<<<<<<<<<<<<<<<<<<\n\n");
	
		for( ; iterDiag != iterDiagEnd; ++iterDiag )
		{
			fprintf(fOut, "Currency:\t%s\n", ((*iterDiag).first).c_str() );
			(*iterDiag).second->View( id, fOut );

		}

		fprintf(fOut, "\n\n <<<<<<<<<<<<<<<<<<<<<< End FX Correl Map >>>>>>>>>>>>>>>>>>>>>>>\n\n");
		fprintf(fOut, "\n xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx \n");
	}


    if ( ficOut == NULL )
       fclose(fOut);
};

double ARM_CorrelatorManager::ComputeIndexIndexCorrel(const string& ccy, const string& tenor1, const string& tenor2,
														  double expiry1, double expiry2) const
{	
	string ccyStr( ccy );
	string tenor2Str( tenor2 );
	string tenor1Str( tenor1 );

	return GetCorrelIndexIndexCube( ccyStr )->ComputeIndexIndexCorrel(tenor1Str, expiry1,tenor2Str,expiry2); 
};

double ARM_CorrelatorManager::ComputeHyperCorrel( const string& ccy, const string& tenor1, const string& tenor2,
												  double expiry, double strike) const
{
	string ccyStr( ccy );

	double dTenor1 =StringMatuToYearTerm( (char*)tenor1.c_str() );
	double dTenor2 =StringMatuToYearTerm( (char*)tenor2.c_str() );

	return GetCorrelHyperCube( ccyStr )->ComputeHyperCorrel(expiry, dTenor1, dTenor2, strike );
	
};

double ARM_CorrelatorManager::ComputeCorrelByExpiry(const string& ccy, const string& tenor1, const string& tenor2,
													double expiry) const
{
	string ccyStr( ccy );

	double dTenor1 =StringMatuToYearTerm( (char*)tenor1.c_str() );
	double dTenor2 =StringMatuToYearTerm( (char*)tenor2.c_str() );

	return GetCorrelHyperCube( ccyStr )->ComputeCorrelByExpiry( expiry, dTenor1, dTenor2 );
};


double ARM_CorrelatorManager::ComputeCorrModeCorrel( const string& ccy, const string& tenor, double expiry, double strike ) const 
{

	string ccyStr( ccy );
	double underlying =StringMatuToYearTerm( (char*)tenor.c_str() );

	return GetCorrModeCorrelCurve( ccyStr )->ComputeVolatility( expiry, underlying );

	
};

double ARM_CorrelatorManager::ComputeSimpleModeCorrel( const string& ccy, const string& tenor, double expiry, double strike ) const 
{

	string ccyStr( ccy );
	double underlying =StringMatuToYearTerm( (char*)tenor.c_str() );

	return GetSimpleModeCorrelCurve( ccyStr )->ComputeVolatility( expiry, underlying );

	
};

double ARM_CorrelatorManager::ComputeCorrelData(const string& mktTag, const string& intraMktTag, double x, double y ) const
{
	vector<string> subTags = splitStringIntoPieces(intraMktTag,"_");
	string	tenor1,tenor2;
	double	correl, dtenor1, dtenor2;

	switch( subTags.size() )
	{
		case 5:
		{
			if( ( subTags[4] != "COMP" ) &&  ( subTags[4] != "NOCHECK" ) )

				throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, 
								"intraMktTag should be ccy1_Tenor1_Ccy2_DIAG, or ccy1_Tenor1_Ccy2_Tenor2_COMP ccy1_Tenor1_Ccy2_Tenor2_NOCHECK or ccy1|ind1_ccy2|ind2 sorted in alphabetical order" );

			// jusqu'à mnt on suppose que dans l'expression ccy1_Tenor1_Ccy2_Tenor2_COMP ou ccy1_Tenor1_Ccy2_Tenor2_NOCHECK ccy1==ccy2
			// il faut prévoir un traitement pour ccy1 != ccy2
			dtenor1 = ConvertCOTenorToDouble(subTags[1].c_str());
			dtenor2 = ConvertCOTenorToDouble(subTags[3].c_str());
			tenor1 = ConvertYearTermToStringMatu(dtenor1);
			tenor2 = ConvertYearTermToStringMatu(dtenor2);
			correl = ComputeIndexIndexCorrel(subTags[0],tenor1,tenor2,x,y);
			break;
		}
		case 4:
		{
			if( subTags[3] != "DIAG" )

				throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, 
								"intraMktTag should be ccy1_Tenor1_Ccy2_DIAG, or ccy1_Tenor1_Ccy2_Tenor2_COMP ccy1_Tenor1_Ccy2_Tenor2_NOCHECK or ccy1|ind1_ccy2|ind2 sorted in alphabetical order" );

			dtenor1 = ConvertCOTenorToDouble(subTags[1].c_str());
			tenor1 = ConvertYearTermToStringMatu(dtenor1);
			tenor2 = ConvertYearTermToStringMatu(y);

			if(mktTag == "IR/VOL")
				correl = ComputeIRVolCorrel(subTags[0], tenor1, tenor2, x);
			else if(mktTag == "VOL/VOL")
				correl = ComputeVolVolCorrel(subTags[0], tenor1, tenor2, x);
			else
				correl = ComputeHyperCorrel(subTags[0], tenor1, tenor2, x);
			break;
		}
		case 2:
		{
			// dans ce cas on traite les ccy_index ou index_ccy ou ccy/index1_ccy/index2
			// le dernier cas on l'assimile avec les autres et on prend tout le intraMktTag comme clé.
		
			tenor2 = ConvertYearTermToStringMatu(y);

			if( subTags[0] == "CORR" )
			{
				correl = ComputeCorrModeCorrel(subTags[1],tenor2,x);
			}
			else if( subTags[1] == "CORR" )
			{
				correl = ComputeCorrModeCorrel(subTags[0],tenor2,x);
			}
			else
			{
				string tempstr = subTags[0];
				int len = tempstr.size();

				if( ( tempstr.find("CMS") != -1 ) || ( tempstr.find("BOR") != -1 ) )
					correl = ComputeSimpleModeCorrel(subTags[1],tenor2,x);
				else
					correl = ComputeSimpleModeCorrel(subTags[0],tenor2,x);
			}
			break;
		}
		default:
			throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, 
								"intraMktTag should be ccy1_Tenor1_Ccy2_DIAG, or ccy1_Tenor1_Ccy2_Tenor2_COMP ccy1_Tenor1_Ccy2_Tenor2_NOCHECK or ccy1|ind1_ccy2|ind2 sorted in alphabetical order" );
	}

	return correl;
}


double ARM_CorrelatorManager::ComputeIndexIndexCorrelData(char* strccy1, const string& index1Name, 
														  char* strccy2, const string& index2Name,
														  double x1, double x2, bool aAcceptLiborMode)
{
	double	vCorrel = 0.0;

	if(aAcceptLiborMode)
	{
		try
		{
			vCorrel = ARM_CorrelManager::ComputeIndexIndexCorrelData(strccy1, index1Name, 
																	 strccy2, index2Name,
																	 x1, x2);
			return	vCorrel;
		}
		catch(...)
		{
			vCorrel = 0.0;
		}
	}

	string	intraMktTag;

	string	ccy1(strccy1);
	string	ccy2(strccy2);

	string	vIndex1Name = index1Name;
	string	vIndex2Name = index2Name;

	if( vIndex1Name.find("CMS") != string::npos )
	{
		vIndex1Name.replace(0, 3, "CO");

		if( vIndex1Name.at(vIndex1Name.size() - 1) != 'Y' )
			vIndex1Name += "Y";
	}
	else if( vIndex1Name.find("LIBOR") != string::npos )
		vIndex1Name.replace(0, 5, "CO");
	else if( vIndex1Name.find("EURIBOR") != string::npos )
		vIndex1Name.replace(0, 7, "CO");

	if( vIndex2Name.find("CMS") != string::npos )
	{
		vIndex2Name.replace(0, 3, "CO");

		if( vIndex2Name.at(vIndex2Name.size() - 1) != 'Y' )
			vIndex2Name += "Y";
	}
	else if( vIndex2Name.find("LIBOR") != string::npos )
		vIndex2Name.replace(0, 5, "CO");
	else if( vIndex2Name.find("EURIBOR") != string::npos )
		vIndex2Name.replace(0, 7, "CO");

	if ( ccy1 == ccy2 )
	{
	   if ( ConvertCOTenorToDouble(vIndex1Name.c_str()) < ConvertCOTenorToDouble(vIndex2Name.c_str()) ) 
	   {
		  intraMktTag = ccy1 + "_" + vIndex1Name + "_" + ccy1 + "_" + vIndex2Name +"_"+"NOCHECK";

		  // because of Summit the correl matrix is (LIBOR, CMS) and
		  // because of the correl manager ordering :
		  if( (index1Name.substr(0, 3) == "CMS") &&
			  (index2Name.substr(0, 3) != "CMS") ) // EURIBOR or LIBOR

		  {
			 double tmpVal = x1;
			 x1 = x2;
			 x2 = tmpVal;
		  }
	   }
	   else
	   {
		  intraMktTag = ccy1 + "_" + vIndex2Name + "_" + ccy1 + "_" + vIndex1Name +"_"+"NOCHECK";

		  // because of Summit, the correl matrix is (LIBOR, CMS) and
		  // because of the correl manager ordering :
		  if( (index1Name.substr(0, 3) == "CMS") &&
			  (index2Name.substr(0, 3) != "CMS") ) // EURIBOR or LIBOR
		  {	
			 double tmpVal = x1;
			 x1 = x2;
			 x2 = tmpVal;
		  }
	   }
	}
	else if ( string(ccy1) < string(ccy2) )
	{
	   intraMktTag = ccy1 + "_" + vIndex1Name + "_" + ccy2 + "_" + vIndex2Name + "_" + "NOCHECK";
	}
	else
	{
	   intraMktTag = ccy2 + "_" + vIndex2Name + "_" + ccy1 + "_" + vIndex1Name + "_" + "NOCHECK";

	   // because of Summit the correl matrix is (LIBOR, CMS) and
	   // because of the correl manager ordering :

	   double tmpVal = x1;
	   x1 = x2;
	   x2 = tmpVal;
	}
	
	vCorrel = ComputeCorrelData("IR/IR", intraMktTag, x1, x2);

    return	vCorrel;
}


double ARM_CorrelatorManager::ComputeIRVolCorrel(const string& ccy, const string& tenor1, 
												  const string& tenor2, double expiry) const
{
	double	vTenor1 = StringMatuToYearTerm( (char*)tenor1.c_str() );
	double	vTenor2 = StringMatuToYearTerm( (char*)tenor2.c_str() );

	return	GetIRVolCorrelHyperCube(ccy)->ComputeHyperCorrel(expiry, vTenor1, vTenor2);
}


double ARM_CorrelatorManager::ComputeVolVolCorrel(const string& ccy, const string& tenor1, 
												  const string& tenor2, double expiry) const
{
	double	vTenor1 = StringMatuToYearTerm( (char*)tenor1.c_str() );
	double	vTenor2 = StringMatuToYearTerm( (char*)tenor2.c_str() );

	return	GetVolVolCorrelHyperCube(ccy)->ComputeHyperCorrel(expiry, vTenor1, vTenor2);
}

double ARM_CorrelatorManager::ComputeFXCorrel(const string& ccy, const string& ccy2, 
											  const string& tenor, double expiry) const
{
	double	vTenor = StringMatuToYearTerm( (char*)tenor.c_str() );

	ARM_VolCube* cube = GetFXCorrelCube(ccy);

	vector<ARM_VolCurve*>* vVolCurve = cube->GetVols();
	vector<ARM_VolCurve*>::iterator iter = vVolCurve->begin();
	vector<ARM_VolCurve*>::iterator iterEnd = vVolCurve->end();
	
	while( (iter != iterEnd) && (strcmp(ccy2.c_str(),((ARM_VolCurve*)(*iter))->GetCurrency()->GetCcyName())))
		iter++;

	if(iter == iterEnd)
	{
		char vErrorMsg[50];
		sprintf(vErrorMsg, "Couldn't find FX Correl for %s_%s", ccy.c_str(),ccy2.c_str());
		throw Exception(__LINE__, __FILE__, ERR_INVALID_DATA, vErrorMsg, 
						"ARM_CorrelatorManager::ComputeFXCorrel(...) : Correl not found");
	}

	return ((ARM_VolCurve*)(*iter))->ComputeVolatility(expiry, vTenor);

}

void ARM_CorrelatorManager::BumpVolatility(const string& ccy, MapType type, double value, int nthLine, int nthCol,
										   int isCumul, int isAbsolute)
{
	BumpCorrel(ccy, type, value, nthLine, nthCol, isCumul, isAbsolute);
}


void ARM_CorrelatorManager::BumpCorrel(const string& ccy, MapType type, double value, int nthLine, int nthCol,
									   int isCumul, int isAbsolute)
{
	string ccystr( ccy );
	switch( type )
	{
		case IdCorr:
		{
			GetSimpleModeCorrelCurve( ccystr )->BumpVolatility(value, nthLine, nthCol, isCumul, isAbsolute);
			break;
		}
		
		case IndexIndexCorr:
		{
			GetCorrelIndexIndexCube( ccystr )->BumpVolatility(value, nthLine, nthCol, isCumul, isAbsolute);
			break;
		}
		
		case Corr:
		{
			GetCorrModeCorrelCurve( ccystr )->BumpVolatility(value, nthLine,nthCol,isCumul,isAbsolute);
			break;
		}

		case IRVolCorr:
		{
			GetIRVolCorrelHyperCube( ccystr )->BumpVolatility(value, nthLine,nthCol,isCumul,isAbsolute);
			break;
		}

		case VolVolCorr:
		{
			GetVolVolCorrelHyperCube( ccystr )->BumpVolatility(value, nthLine,nthCol,isCumul,isAbsolute);
			break;
		}

		case FXVolCorr:
		{
			vector<string> subTags = splitStringIntoPieces(ccystr,"_");
			ARM_VolCube* cube = GetFXCorrelCube( subTags[0] );
			switch(subTags.size())
			{
			case 1:
				{
					cube->BumpSmile(value,0,nthLine, nthCol, isCumul, isAbsolute);
					break;
				}
			case 2:
				{
					vector<ARM_VolCurve*>* vVolCurve = cube->GetVols();
					vector<ARM_VolCurve*>::iterator iter = vVolCurve->begin();
					vector<ARM_VolCurve*>::iterator iterEnd = vVolCurve->end();
					
					while( (iter != iterEnd) && (strcmp(subTags[1].c_str(),((ARM_VolCurve*)(*iter))->GetCurrency()->GetCcyName())))
						iter++;

					if(iter == iterEnd)
					{
						// char vErrorMsg[50];
						// sprintf(vErrorMsg, "Couldn't find FX Correl for %s_%s", subTags[0].c_str(),subTags[1].c_str());
						// throw Exception(__LINE__, __FILE__, ERR_INVALID_DATA, vErrorMsg, 
						// 				"ARM_CorrelatorManager::BumpCorrel(...) : Correl not found");
						ARMTHROW(ERR_INVALID_DATA,
							"ARM_CorrelatorManager::BumpCorrel(...) : Correl not found: Couldn't find FX Correl for "
							<<subTags[0]<<"_"<<subTags[1]); 
					}

					((ARM_VolCurve*)(*iter))->BumpVolatility(value, nthLine,nthCol,isCumul,isAbsolute);

					break;
				}
			default:
				{
					char vErrorMsg[50];
					sprintf(vErrorMsg, "Couldn't find FX Correl for %s", ccystr);
					throw Exception(__LINE__, __FILE__, ERR_INVALID_DATA, vErrorMsg, 
										"ARM_CorrelatorManager::BumpCorrel(...) : Correl not found");
				}

			}
			break;
		}
		
		default:
		{
			GetCorrelHyperCube( ccystr )->BumpVolatility(value, nthLine, nthCol, isCumul, isAbsolute);
		}
	}
	
};

/*---------------------------------------------------------------------------*/
/*---- End of file ----*/

