/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 *	\file mktdataretriever.cpp
 *
 *  \brief
 *
 *	\author  E. Benhamou
 *	\version 1.0
 *	\date March 2004
 */

#include "gpbase/removeidentifiedwarning.h"
#include "gpinfra/mktdataretriever.h"

/// gpinfra
#include "gpinfra/mktdatamanager.h"
#include "gpinfra/argconvdefault.h"

/// gpbase
#include "gpbase/singleton.h"
#include "gpbase/autocleaner.h"
#include "gpbase/env.h"
#include "gpbase/stringconvert.h"
#include "gpbase/gplinalgconvert.h"
#include "gpbase/gpvector.h"
#include "gpbase/gpmatrix.h"

/// kernel
#include <ccy/currency.h>
#include <crv/volint.h>
#include <crv/volcube.h>
#include <crv/zeroint.h>
#include <mod/bsmodel.h>
/// STL
#include <iomanip>
#include <memory>
#include <cstdio>	/// for fopen... and other file access ... 
					/// we did not use ifstream for performance reason
					/// as the first version shows that this was slower!
CC_BEGIN_NAMESPACE( ARM )


/////////////////////////////////////////////////////////////////
///	Class  : ARM_MarketData_RetrieverImp
///	Routine: GetSummitFile
///	Returns: an ifstream*
///	Action : finds file according to the commonFolderName and the 
///		commonFileName (try the standard one and the archive file)
///		returns a pointor to avoid destroying the object
/////////////////////////////////////////////////////////////////

FILE* ARM_MarketData_RetrieverImp::GetSummitFile(
	const string& commonFolderName,
	const string& commonFileName,
	const ARM_Date& asOfDate ) const
{
	/// design is to try two folders: 
	///		-the standard one 
	///		-the archive folder ( which have on top of the standard folder the extra argument SUMMIT_FILE_LOCATION_HISTO )
	CC_Ostringstream fileName_os;
	fileName_os << commonFileName;
	
	/// US have a different convention for date hence a different storing!
	if (strcmp(ARM_DEFAULT_COUNTRY,"USD") == 0)
	{
		/// convention is yyyy/mm/dd
		/// if the month is less than 10, fills it with a zero and same for day
		fileName_os << const_cast<ARM_Date&>(asOfDate).GetYear() 
			<< CC_NS(std,setw)(2) << CC_NS(std,setfill)('0') << const_cast<ARM_Date&>(asOfDate).GetMonth() 
			<< CC_NS(std,setw)(2) << CC_NS(std,setfill)('0') << const_cast<ARM_Date&>(asOfDate).GetDay();
	}
	else
	{
		/// convention is dd/mm/yyyy
		fileName_os << CC_NS(std,setw)(2) << CC_NS(std,setfill)('0') << const_cast<ARM_Date&>(asOfDate).GetDay() 
			<< CC_NS(std,setw)(2) << CC_NS(std,setfill)('0') << const_cast<ARM_Date&>(asOfDate).GetMonth();

		/// in Europe, we truncate the year to the last two digits!
		CC_Ostringstream fullYearString;
		fullYearString << const_cast<ARM_Date&>(asOfDate).GetYear();
		string yearString( fullYearString.str() );
		yearString = yearString.substr(2);
		fileName_os << yearString;
	}
	
	/// extension of the file!
	fileName_os << ".000";
	
	/// get the two corresponding file name one for the standard folder and the one for the archive folder
	string stdFileName( commonFolderName + fileName_os.str() );
	string archiveFileName( commonFolderName + SUMMIT_FILE_LOCATION_HISTO + fileName_os.str() );

	FILE* dataFile = fopen( stdFileName.c_str(), "r" );
	if (dataFile == NULL )
	{
		dataFile = fopen(archiveFileName.c_str(), "r" );
		if( dataFile == NULL )
		{
			delete dataFile;
			CC_Ostringstream os;
			os << ARM_USERNAME << ": could not open std file " << stdFileName
				<< " and archive " << archiveFileName;
			throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, os.str() );
		}
	}
	return dataFile;
}


/////////////////////////////////////////////////////////////////
///	Class  : ARM_MarketData_RetrieverImp
///	Routine: ConvertDateOrMatuString
///	Returns: ARM_Date
///	Action : convert a dateString or a maturity according to the default setting
/////////////////////////////////////////////////////////////////

ARM_Date ARM_MarketData_RetrieverImp::ConvertDateOrMatuString( const string& matuOrDate,
	const string& ccy,
	const ARM_Date& asOfDate,
	int Rule ) const
{
	/// test if it is a standard period
	long nb,freqId;
	char cMatu;
	
	sscanf(matuOrDate.c_str(), "%ld%c", &nb, &cMatu);
	cMatu = toupper(cMatu);

	/// sorted in their order of frequent use
	if ( cMatu == 'Y' )
		freqId = K_ANNUAL;
	else if ( cMatu == 'M' ) 
		freqId = K_MONTHLY;
	else if ( cMatu == 'D' )
		freqId = K_DAILY;
	else if ( cMatu == 'W' )  
		freqId = K_WEEKLY;
	else
	/// this is not a period and should be converted from string to date
		return ConvertDateStringToDate( matuOrDate );

	ARM_Date tmpDate(asOfDate);

	if( freqId == K_DAILY )
		tmpDate.NextBusinessDay( nb, (char*) ccy.c_str() );
	else
	{
		tmpDate.AddPeriodMult( freqId, nb, (char*) ccy.c_str() );
		if( Rule )
			tmpDate.AdjustToBusDate( (char*) ccy.c_str(), Rule );
	}
	
	return tmpDate;
}


/////////////////////////////////////////////////////////////////
///	Class  : ARM_MarketData_RetrieverImp
///	Routine: ConvertVectorFromMatuOrDateToYearFrac
///	Returns: ARM_GP_Vector* 
///	Action : convert a vector of dateString or a maturity according to a vector of year frac in double
/////////////////////////////////////////////////////////////////
ARM_GP_Vector* ARM_MarketData_RetrieverImp::ConvertVectorFromMatuOrDateToYearFrac( const vector<string>& vecData, 
	const ARM_Date& asOfDate,
	const string& ccy ) const
{
	ARM_GP_Vector* result = new ARM_GP_Vector(vecData.size());
	for( size_t i=0; i<vecData.size(); ++i )
		(*result)[i]= (double) (ConvertDateOrMatuString(vecData[i],ccy,asOfDate)-const_cast<ARM_Date&>(asOfDate) )/ K_YEAR_LEN;
	return result;
}


/////////////////////////////////////////////////////////////////
///	Class  : ARM_MarketData_RetrieverImp
///	Routine: ConvertVectorFromMatuOrDateWRefDateToYearFrac
///	Returns: ARM_GP_Vector*
///	Action : convert a vector of dateString or a maturity according to a vector of year frac in double
/////////////////////////////////////////////////////////////////

ARM_GP_Vector* ARM_MarketData_RetrieverImp::ConvertVectorFromMatuOrDateWRefDateToYearFrac( const vector<string>& vecData, 
	const ARM_Date& asOfDate,
	const ARM_Date& refDate,
	const string& ccy ) const
{
	ARM_GP_Vector* result = new ARM_GP_Vector(vecData.size());
	for( size_t i=0; i<vecData.size(); ++i )
		(*result)[i]= (double) (ConvertDateOrMatuString(vecData[i],ccy,refDate, K_FOLLOWING)-const_cast<ARM_Date&>(asOfDate) )/ K_YEAR_LEN;
	return result;
}


/////////////////////////////////////////////////////////////////
///	Class  : ARM_MarketData_RetrieverImp
///	Routine: GetZCCurveFromSummit
///	Returns: ARM_ZeroCurve*
///	Action : get a curve according to its type!
/////////////////////////////////////////////////////////////////
ARM_ZeroCurve* ARM_MarketData_RetrieverImp::GetZCCurveFromSummit( 
	const string& index,
	const string& ccy,
	const string cvName,
	const ARM_Date& asOfDate ) const
{
	/// check that the asOfDate is meaningfull
	ValidateAsOfDate( asOfDate );

	/// get the files!
	CC_Ostringstream commonFolder_os;
	commonFolder_os << SUMMIT_DATA_FOLDER << SUMMIT_ZC_FILE_LOCATION;
	CC_Ostringstream commonFileName_os;
	commonFileName_os <<  "ZC_" << ccy <<  "_" << index << "_" << cvName << ".";

	FILE* dataFile = GetSummitFile( commonFolder_os.str(), commonFileName_os.str(), asOfDate );
	
	/// read the file! reserve 30 for performance reason!
	vector<double> matu;
	matu.reserve(30);
	vector<double> rate;
	rate.reserve(30);

	char ignoreString[50],dateString[50];
	double value;
	int readFileReturn = 0;

	/// order is date (ignoreString), dateString, value
	while( readFileReturn != EOF )
	{
		fscanf( dataFile, "%s", ignoreString );
		fscanf( dataFile, "%s", dateString );
		ARM_Date tmpDate = ConvertDateStringToDate( dateString );
		matu.push_back( (double) (tmpDate-const_cast<ARM_Date&>(asOfDate) )/K_YEAR_LEN );
		
		readFileReturn = fscanf( dataFile, "%lf", &value );
		rate.push_back( value );
	}

	fclose( dataFile );

	/// uses auto_ptr for exception safety!
	CC_NS(std,auto_ptr)<ARM_Vector> mktData( new ARM_Vector( rate.size(), &rate[0] ) );
	CC_NS(std,auto_ptr)<ARM_Vector> mktMatu( new ARM_Vector( matu.size(), &matu[0] ) );
	CC_NS(std,auto_ptr)<ARM_Currency> ccyObj( new ARM_Currency( (char*) ccy.c_str() ) );
	ARM_ZeroCurve* curve = new ARM_ZeroLInterpol( const_cast<ARM_Date&>(asOfDate), mktMatu.get(), mktData.get(), 0, 0, K_LINEAR );
	curve->SetCurrencyUnit(ccyObj.get());
	return curve;
}

/////////////////////////////////////////////////////////////////
///	Class  : ARM_MarketData_RetrieverImp
///	Routine: GetVolFromSummit_Common
///	Returns: ARM_VolLInterpol*
///	Action : get a vol curve according to its type!
/////////////////////////////////////////////////////////////////

ARM_VolLInterpol* ARM_MarketData_RetrieverImp::GetVolFromSummit_Common(
	const string& commonFolderName,
	const string& commonFileName,
	const string& ccy,
	const ARM_Date& asOfDate,
	int volType ) const
{
	/// check that the asOfDate is meaningfull
	ValidateAsOfDate( asOfDate );
	
	/// flat file parsing
	/// the result!
	ARM_VolLInterpol* volCrv = NULL;

	/// get the files!
	FILE* dataFile = GetSummitFile( commonFolderName, commonFileName, asOfDate );

	/// read the file!
	ARM_GP_Vector* vMaturitiesOrStrikes = new ARM_GP_Vector;
	vMaturitiesOrStrikes->reserve( 30 );
	vector<string> vYearTerms;
	vYearTerms.reserve(30);
	vector<double> vVols;
	vVols.reserve(250);

	char ignoreString[50],dateOrMatuString[50];
	double value;
	
	/// data are read row by row...assumption is that there
	/// is no hole in the matrix!
	/// hence to get column data we just need to read the first row
	bool readingFirstRow= true;
	
	if( volType == K_ATMF_VOL )
	{
		vector<string> vmaturitiesString;
		vmaturitiesString.reserve( 30 );

		/// order is maturities (string), year terms (string), vols (double)
		int readFileReturn = fscanf( dataFile, "%s", ignoreString );
		while( readFileReturn != EOF )
		{
			fscanf( dataFile, "%s", dateOrMatuString );

			if( vmaturitiesString.empty() )
				vmaturitiesString.push_back( dateOrMatuString );
			else 
				/// are we not at the beginning of a new row?
				if( vmaturitiesString[vmaturitiesString.size()-1] != dateOrMatuString )
				{
					vmaturitiesString.push_back( dateOrMatuString );
					readingFirstRow = false;
				}
				
			fscanf( dataFile, "%s", dateOrMatuString );

			/// because the first row is representative of all column
			/// data, we just need to read the first row!
			if( readingFirstRow )
				vYearTerms.push_back( dateOrMatuString );

			/// read all the vols
			fscanf( dataFile, "%lf", &value );
			vVols.push_back( value );

			/// skip one string
			readFileReturn = fscanf( dataFile, "%s", ignoreString );
		}

		/// Changes data from string type to double maturity type
		vMaturitiesOrStrikes = ConvertVectorFromMatuOrDateToYearFrac( vmaturitiesString, asOfDate, ccy );
	}
	else if( volType == K_SMILE_VOL )
	{
		/// order is year terms (string), strikes (double), vols (double)
		int readFileReturn = fscanf( dataFile, "%s", ignoreString );
		while( readFileReturn != EOF )
		{
			fscanf( dataFile, "%s", dateOrMatuString );

			if( vYearTerms.empty() )
				vYearTerms.push_back( dateOrMatuString );
			else 
				/// are we not at the beginning of a new row?
				if( vYearTerms[vYearTerms.size()-1] != dateOrMatuString )
				{
					vYearTerms.push_back( dateOrMatuString );
					readingFirstRow = false;
				}

			fscanf( dataFile, "%lf", &value );

			/// because the first row is representative of all column
			/// data, we just need to read the first row!
			if( readingFirstRow )
				vMaturitiesOrStrikes->push_back( value );

			/// read all the vols
			fscanf( dataFile, "%lf", &value );
			vVols.push_back( value );
	
			/// skip one string
			readFileReturn = fscanf( dataFile, "%s", ignoreString );
		}
	}
	else
		throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, "Unknown vol case!" );


	/// closes the file
	fclose(dataFile);

	/// process year tearms
	CC_NS(std,auto_ptr)<ARM_Currency> ccyObj( new ARM_Currency( (char*) ccy.c_str() ) );
	long spotDays = ccyObj->GetSpotDays();
	ARM_Date settleDate(asOfDate);
	settleDate.NextBusinessDay(spotDays, (char*) ccy.c_str());
	ARM_GP_Vector* tmpYearTerms = ConvertVectorFromMatuOrDateWRefDateToYearFrac( vYearTerms, asOfDate, settleDate, ccy );

	/// creates the vol matrix
	ARM_Vector* yearTerms = To_pARM_Vector(tmpYearTerms);
	delete tmpYearTerms;

	ARM_Vector* maturitiesOrStrikes = To_pARM_Vector( vMaturitiesOrStrikes );
	delete vMaturitiesOrStrikes;

	ARM_GP_Matrix tmpVols( yearTerms->size(), maturitiesOrStrikes->size(), &vVols[0] );

	/// in the ATM case, the vol are stored in the opposite order
	/// hence we need to reassign!
	if( volType == K_ATMF_VOL )
		for (size_t i=0;i<tmpVols.rows(); ++i )
			for (size_t j=0; j<tmpVols.cols(); ++j )
				tmpVols(i,j) = vVols[j*tmpVols.rows()+i];

	ARM_Matrix* vols = To_pARM_Matrix( &tmpVols );
		
	/// creates the vol object!
	volCrv = new ARM_VolLInterpol( asOfDate, yearTerms, maturitiesOrStrikes, vols, 1, volType );
	volCrv->SetCurrencyUnit(ccyObj.get() );

	return volCrv;
}


/////////////////////////////////////////////////////////////////
///	Class  : ARM_MarketData_RetrieverImp
///	Routine: GetVolATMFromSummit
///	Returns: ARM_VolLInterpol*
///	Action : get an ATM vol curve according to its type!
/////////////////////////////////////////////////////////////////

ARM_VolLInterpol* ARM_MarketData_RetrieverImp::GetVolATMFromSummit(
	const string& index,
	const string& ccy,
	const string& cvName,
	const ARM_Date& asOfDate,
	const string& volMktType ) const
{
	string commonFolder		= SUMMIT_DATA_FOLDER + SUMMIT_VOL_FILE_LOCATION;
	string commonFileName	= string("IRFWDVOL_") + ccy + "_" + index + "_" + volMktType + "_" + cvName  +".";
	return GetVolFromSummit_Common(	commonFolder, commonFileName, ccy, asOfDate, K_ATMF_VOL );
}


/////////////////////////////////////////////////////////////////
///	Class  : ARM_MarketData_RetrieverImp
///	Routine: GetVolSmileFromSummit
///	Returns: ARM_VolLInterpol*
///	Action : get a smile vol curve according to its type!
/////////////////////////////////////////////////////////////////

ARM_VolLInterpol* ARM_MarketData_RetrieverImp::GetVolSmileFromSummit(
	const string& index,
	const string& ccy,
	const string& cvName,
	const ARM_Date& asOfDate,
	const string& volMktType,
	const string& matuIndex ) const
{
	string commonFolder		= SUMMIT_DATA_FOLDER + SUMMIT_SMILE_FILE_LOCATION;
	string commonFileName	= string("SMILE_") + ccy + "_" + index + "_" + volMktType + "_" + matuIndex + "_C_" + cvName + ".";
	return GetVolFromSummit_Common(	commonFolder, commonFileName, ccy, asOfDate, K_SMILE_VOL );
}	


/////////////////////////////////////////////////////////////////
///	Class  : ARM_MarketData_RetrieverImp
///	Routine: GetVolCubeFromSummit
///	Returns: ARM_VolCube*
///	Action : get a smile vol curve according to its type!
/////////////////////////////////////////////////////////////////
ARM_VolCube* ARM_MarketData_RetrieverImp::GetVolCubeFromSummit(
	const string& index,
	const string& ccy,
	const string& cvName,
	const ARM_Date& asOfDate,
	const string& volMktType ) const
{
	ARM_VolCurve* ATMVolCrv	= GetVolATMFromSummit( index, ccy, cvName, asOfDate, volMktType );
	ARM_MarketData::VolMktType mktType	= (ARM_MarketData::VolMktType ) ARM_ArgConv_VolMktType.GetNumber(volMktType);

	vector<string> vTenorsString;
	vector<double> vTenors;
	const string* tenorsString;
	size_t i, totalSize ;
	vector<ARM_VolCurve*> inVols;

	/// FIX FIX should be more detailled with specification by index and currency
	switch( mktType )
	{
	case ARM_MarketData::MKT_CAPORCAPLET_VOL:
		{
			tenorsString = CAPLET_TENORS_TABLE;
			totalSize	 = CAPLET_TENORS_TABLE_SIZE;
		}
		break;
	case ARM_MarketData::MKT_SWAPTION_VOL:
		{
			tenorsString = SWO_TENORS_TABLE;
			totalSize	 = SWO_TENORS_TABLE_SIZE;
		}
		break;
	default:
		throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, "Unknown vol mkt type!" );
	}

	vTenorsString.resize(totalSize);
	vTenors.resize(totalSize);
	inVols.resize(totalSize);

	for(i=0; i<totalSize; ++i)
	{
		vTenorsString[i]= tenorsString[i];
		vTenors[i]		= StringMaturityToYearTerm(vTenorsString[i]);			
		inVols[i]		= GetVolSmileFromSummit( index, ccy, cvName, asOfDate, volMktType, vTenorsString[i] ); 
	}

	ARM_Vector tenors( totalSize,&vTenors[0] );
	return new ARM_VolCube(ATMVolCrv, &inVols[0], totalSize, &tenors );
}


/////////////////////////////////////////////////////////////////
///	Class  : ARM_MarketData_RetrieverImp
///	Routine: GetVolATMFromSummit
///	Returns: ARM_VolLInterpol*
///	Action : get an ATM vol curve according to its type!
/////////////////////////////////////////////////////////////////

ARM_VolLInterpol* ARM_MarketData_RetrieverImp::GetVolFromSummit(
	const string& index,
	const string& ccy,
	const string& cvName,
	const ARM_Date& asOfDate,
	ARM_MarketData::VolMktType volMktType,
	ARM_MarketData::VolType volType ) const
{
	switch( volType )
	{
	case ARM_MarketData::MKT_ATM_VOL:
		return GetVolATMFromSummit( index, ccy, cvName, asOfDate, ARM_ArgConvReverse_VolMktType.GetString(volMktType) );
	case ARM_MarketData::MKT_ATM_AND_SMILE_VOL:
		return GetVolCubeFromSummit( index, ccy, cvName, asOfDate, ARM_ArgConvReverse_VolMktType.GetString(volMktType) );
	default:
		throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, ARM_USERNAME + ": type unknown!" );
	}

}

/////////////////////////////////////////////////////////////////
///	Class  : ARM_MarketData_RetrieverImp
///	Routine: GetVolMktModel
///	Returns: ARM_Model*
///	Action : get a model for a given vol!
/////////////////////////////////////////////////////////////////

ARM_Model* ARM_MarketData_RetrieverImp::GetVolMktModel( 
		const string& indexName, 
		const string& ccy, 
		const string& cvName, 
		const ARM_Date& asOf, 
		const string& volMktType,
		const string& source ) const
{
	ARM_VolCurve* volCurv		= ARM_TheMarketData_Manager.Instance()->GetVolCurve(
		indexName, ccy, cvName, asOf, volMktType, "ATMANDSMILE", source );
	ARM_ZeroCurve* zcCurve	= ARM_TheMarketData_Manager.Instance()->GetZCCurve(
		indexName, ccy, cvName, asOf, source );
	return new ARM_BSModel( const_cast<ARM_Date&>(asOf), 0, zcCurve	, zcCurve, volCurv );
}


/// Creation of the market data retriever
ARM_SingletonHolder<ARM_MarketData_RetrieverImp> ARM_TheMarketDataRetriever;


CC_END_NAMESPACE()

/*---------------------------------------------------------------------------*/
/*---- End of file ----*/

