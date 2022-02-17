/*
 *
 * Copyright (c) CDC IXIS CM November 2003 Paris
 *
 * $Log: seasonmanager.cpp,v $
 * Revision 1.3  2003/12/15 15:21:51  arm
 * Correction of Syntax error in:
 * sprintf(msg, "%s is not a valdid month name.", month.c_str());
 *
 * Revision 1.2  2003/11/25 16:06:04  rguillemot
 * Inflation Seasonality Manager
 *
 * Revision 1.1  2003/11/25 10:15:44  rguillemot
 * Initial revision
 *
 *
 */

/*----------------------------------------------------------------------------*
    seasonmanager.cpp
	Copyright (c) CDC IXIS CM November 2003 Paris
*----------------------------------
------------------------------------------*/

#include "gpinflation/seasonmanager.h"

/// gpbase
#include "gpbase/gpvector.h"
#include "gpbase/gplinalgtypedef.h"
#include "gpbase/gplinalginterpol.h"
#include "gpbase\gpmatrix.h"

#include <glob/dates.h>

CC_BEGIN_NAMESPACE( ARM )


////////////////////////////////////////////////////
///	Class  : ARM_SeasonalityManager
///	Routine: ValidateData
///	Returns: void
///	Action : Validation of the input data of the seasonality manager
////////////////////////////////////////////////////

void ARM_SeasonalityManager::ValidateData(const vector<string>& monthList, const ARM_GP_Vector& seasonSpread, const vector<double>& horizonList)
{
	int numberOfHorizons(0);
	numberOfHorizons = (horizonList.size()) ? horizonList.size() : 1;

	// The month list and the season spread list must have the same size
	if (numberOfHorizons*monthList.size() != seasonSpread.size())
		throw Exception(__LINE__, __FILE__, ERR_INVALID_DATA, "The season spread list must be of size : months number * horizons number." );

	// The month list and the season spread list must have the same size
	if (numberOfHorizons*monthList.size() != itsSeasonData.size())
		throw Exception(__LINE__, __FILE__, ERR_INVALID_DATA, "The month list size must contains : 12 months and at least 1 horizon." );
}



////////////////////////////////////////////////////
///	Class  : ARM_SeasonalityManager
///	Routine: CheckDataIntegrity
///	Returns: void
///	Action : Test if all the season data have been inputed
////////////////////////////////////////////////////

void ARM_SeasonalityManager::CheckDataIntegrity()
{
	int i = 0;
	for (i = 0; i < itsHasBeenStored.size(); ++i)
	{
		if (!itsHasBeenStored[i])
		{
			char msg[255];
			sprintf(msg, "Season spread for %s is missing.", ARM_Date::GetLongMonthName(i+1));
			throw Exception(__LINE__, __FILE__, ERR_INVALID_DATA, msg);
		}
	}
}


////////////////////////////////////////////////////
///	Class  : ARM_SeasonalityManager
///	Routine: Init
///	Returns: void
///	Action : Init method
////////////////////////////////////////////////////
void ARM_SeasonalityManager::Init(const vector<double>& horizonList)
{
	SetName(ARM_SEASONMANAGER);

	// The season spread data contain 12 entries one per month per horizon
	int numberOfHorizons(0);
	numberOfHorizons = (horizonList.size()) ? horizonList.size() : 1;

	itsSeasonData.resize(12*numberOfHorizons);
	itsHasBeenStored.resize(itsSeasonData.size());

	// At the beginning no data have been inputed
	int i = 0;
	for (i = 0; i < itsSeasonData.size(); ++i)
		itsHasBeenStored[i] = false;

}


////////////////////////////////////////////////////
///	Class  : ARM_SeasonalityManager
///	Routine: Factorize
///	Returns: void
///	Action : Factorize
////////////////////////////////////////////////////
void ARM_SeasonalityManager::Factorize(
		const std::vector<double>& horizonList, 
		const std::vector<std::string>& monthList,
		const ARM_GP_Vector& seasonSpreadList)
{
	int numberOfHorizons = horizonList.size();
	int numberOfMonths	 = monthList.size();
	itsSeasonTermData.resize( numberOfMonths, numberOfHorizons); 
	

	for (int i = 0; i< numberOfMonths; ++i)
	{
		for (int j = 0; j< numberOfHorizons; ++j)
		{
			AddMonthSeasonSpread(monthList[i], seasonSpreadList.Elt(j+i*numberOfHorizons), j*numberOfMonths);
			itsSeasonTermData.Elt(i, j)	= seasonSpreadList.Elt(j+i*numberOfHorizons);
		}
	}
	
}


////////////////////////////////////////////////////
///	Class  : ARM_SeasonalityManager
///	Routine: Constructor
///	Returns: void
///	Action : Constructor
////////////////////////////////////////////////////

ARM_SeasonalityManager::ARM_SeasonalityManager(const std::vector<std::string>& monthList, 
	const ARM_GP_Vector& seasonSpreadList, 
	CorrectionMode mode )
:
	itsMode(mode), itsHorizonTerms(NULL), itsSeasonTermData(NULL)
{
	Init(itsHorizonTerms);
	// Do some data validation
	ValidateData(monthList,seasonSpreadList, itsHorizonTerms);

	int i = 0;
	int numberOfData = monthList.size();

	for (i = 0; i < numberOfData; ++i)
		AddMonthSeasonSpread(monthList[i], seasonSpreadList.Elt(i), 0);

	CheckDataIntegrity();
}

////////////////////////////////////////////////////
///	Class  : ARM_SeasonalityManager
///	Routine: Overloaded Constructor
///	Returns: void
///	Action : Overloaded Constructor
////////////////////////////////////////////////////

ARM_SeasonalityManager::ARM_SeasonalityManager( const std::vector<std::string>& monthList,
	const ARM_GP_Vector& seasonSpreadList,
	const std::vector<double>& horizonList,
	CorrectionMode mode)
:
	itsMode(mode), itsHorizonTerms (horizonList)
{
	Init(horizonList);
	// Cehck the dimensions of the tables and the references
	ValidateData(monthList,seasonSpreadList, horizonList);

	/// Distribute the horizons to the data vector 
	Factorize(horizonList, monthList, seasonSpreadList);
	
	// Do some data valida
	CheckDataIntegrity();
	
}



////////////////////////////////////////////////////
///	Class  : ARM_SeasonalityManager
///	Routine: Copy Constructor
///	Returns: 
///	Action : Copy Constructor
////////////////////////////////////////////////////

ARM_SeasonalityManager::ARM_SeasonalityManager(const ARM_SeasonalityManager& rhs)
:
	ARM_Object(rhs),
	itsSeasonData( rhs.itsSeasonData ),
	itsSeasonTermData( rhs.itsSeasonTermData ),
	itsHorizonTerms( rhs.itsHorizonTerms ),
	itsMode( rhs.itsMode )
{}



////////////////////////////////////////////////////
///	Class  : ARM_SeasonalityManager
///	Routine: Destructor
///	Returns: 
///	Action : Destructor
////////////////////////////////////////////////////
ARM_SeasonalityManager::~ARM_SeasonalityManager()
{}


////////////////////////////////////////////////////
///	Class  : ARM_SeasonalityManager
///	Routine: Operator = 
///	Returns: 
///	Action : Operator = 
////////////////////////////////////////////////////

ARM_SeasonalityManager& ARM_SeasonalityManager::operator=(const ARM_SeasonalityManager &rhs)
{
	if( this != &rhs )
	{
		ARM_Object::operator = ( rhs );
		itsSeasonData		= rhs.itsSeasonData;
		itsSeasonTermData	= rhs.itsSeasonTermData;
		itsHorizonTerms		= rhs.itsHorizonTerms;
		itsMode				= rhs.itsMode;
	}

	return *this;
}


////////////////////////////////////////////////////
///	Class  : ARM_SeasonalityManager
///	Routine: Clone
///	Returns: 
///	Action : Clone
////////////////////////////////////////////////////

ARM_Object* ARM_SeasonalityManager::Clone()
{
	return new ARM_SeasonalityManager( *this );
}



////////////////////////////////////////////////////
///	Class  : ARM_SeasonalityManager
///	Routine: AddMonthSeasonSpread
///	Returns: 
///	Action : Add a new season spread in the manager for a corresponding month
////////////////////////////////////////////////////

void ARM_SeasonalityManager::AddMonthSeasonSpread(const string& month, double SeasonSpread, int horizon)
{	
	int monthNumber = GetNumMonthFromStr((char*) month.c_str());

	// To prevent wrong month name
	if ((monthNumber <= 0) || (monthNumber > 12))
	{
		char msg[255];
		sprintf(msg, "%s is not a valdid month name.", month.c_str());
		throw Exception(__LINE__, __FILE__, ERR_INVALID_DATA, msg);
	}

	// To prevent to input twice the same month
	if (itsHasBeenStored[monthNumber-1+horizon])
	{
		char msg[255];
		sprintf(msg, "The data for %s has been inputed twice.", ARM_Date::GetLongMonthName(monthNumber));
		throw Exception(__LINE__, __FILE__, ERR_INVALID_DATA, msg);
	}

	itsSeasonData[monthNumber-1+horizon] = SeasonSpread;
	itsHasBeenStored[monthNumber-1+horizon] = true;
}



////////////////////////////////////////////////////
///	Class  : ARM_SeasonalityManager
///	Routine: Destructor
///	Returns: 
///	Action : view the details of the seasonality manager
///				Display the season spread curve
////////////////////////////////////////////////////

void ARM_SeasonalityManager::View(char* id, FILE* ficOut)
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

    fprintf(fOut, "\n\n\t =====> Seasonality Manager \n\n");
    fprintf(fOut, "\n\n\t   Mode = ");
	switch( itsMode )
	{
	case K_ADDITIVE:
	    fprintf(fOut, "Additive correction\n\n");
		break;
	case K_MULTIPLICATIVE:
	    fprintf(fOut, "Multiplicative correction\n\n");
		break;
	default:
		throw Exception(__LINE__, __FILE__, ERR_INVALID_DATA, "Unknown correction mode for the manager, supported is additive and multiplicative!" );
	}

	int i,j;
	if (!itsSeasonTermData.size())
	{
		for(i = 0; i < itsSeasonData.size(); ++i)
			fprintf(fOut, "%12s \t %+6.6f\n", ARM_Date::GetLongMonthName(i+1), itsSeasonData[i]);
	}
	else
	{
		fprintf(fOut, "    Horizons");
		for(j = 0; j < itsHorizonTerms.size(); ++j)
			fprintf(fOut, "\t %1.0f Y \t", itsHorizonTerms[j]);
		fprintf(fOut, "\n");

		for(i = 0; i < itsSeasonTermData.GetRowsNb(); ++i)
		{
			fprintf(fOut, "%12s \t", ARM_Date::GetLongMonthName(i+1));

			for(j = 0; j < itsSeasonTermData.GetColsNb(); ++j)
				fprintf(fOut, "%+6.6f\t", itsSeasonTermData.Elt(i,j));

			fprintf(fOut, "\n");
		}
		
	}



	if ( ficOut == NULL )
		fclose(fOut);
}




////////////////////////////////////////////////////
///	Class  : ARM_SeasonalityManager
///	Routine: SeasonSpreadCorrect
///	Returns: 
///	Action : Correct the interpolated CPI value by the season spread
///				interpolDateLag is the date of interpolation
///				refDate1Lag is the previous reference date
///				refDate2Lag is the next reference date
///					
///				in an additive scheme,		CPI = CPI (not correction ) + SeasonBump - linCorrect; 
///				in an multiplicative scheme, CPI = CPI (not correction ) * SeasonBump / linCorrect; 
////////////////////////////////////////////////////

void ARM_SeasonalityManager::SeasonSpreadCorrect( 
	ARM_Date asOfDate,
	double& CPIResult, 
	double interpolDateLag, 
	double refDate1Lag, 
	double refDate2Lag) const
{
	/// Recovers the necessary datas
	ARM_Date refDate1(refDate1Lag);
	ARM_Date refDate2(refDate2Lag);
	ARM_Date interpolDate(interpolDateLag);

	int refDate1Month		= refDate1.GetMonth();
	int refDate2Month		= refDate2.GetMonth();
	int interpolDateMonth	= interpolDate.GetMonth();

	double seasonBump;
	double seasonBump1;
	double seasonBump2;
	double linCorrect;

	if (!itsHorizonTerms.size())
	{
		/// In the case of one static correction : no horizon specification
		seasonBump		= itsSeasonData[interpolDateMonth-1];
		seasonBump1		= itsSeasonData[refDate1Month-1];
		seasonBump2		= itsSeasonData[refDate2Month-1];
		linCorrect		= LinearInterpolation<double>(interpolDateLag, refDate1Lag, seasonBump1, refDate2Lag, seasonBump2);
	}
	else
	{
		/// In the case of different bumps by horizon 
		vector<double> exactSeasonTerms;
		exactSeasonTerms.resize(12);

		double interpolHorizon = interpolDate.GetYear()-asOfDate.GetYear()+1;
		double ref1Horizon = refDate1.GetYear()-asOfDate.GetYear()+1;
		double ref2Horizon = refDate2.GetYear()-asOfDate.GetYear()+1;
		double HorizonLeft = 0;
		double HorizonRight = 0;

		if (interpolHorizon >= itsHorizonTerms[itsHorizonTerms.size()-1])
		{
			interpolHorizon = itsHorizonTerms[itsHorizonTerms.size()-1];
			seasonBump2		= itsSeasonTermData.Elt(refDate2Month-1,itsHorizonTerms.size()-1);
		}

		for (int i =0; i< itsHorizonTerms.size()-1; i++)
		{
			if(interpolHorizon >= itsHorizonTerms[i] && interpolHorizon <= itsHorizonTerms[i+1])
			{
				/// Bounding the interpolHorizon
				HorizonLeft			= itsHorizonTerms[i];
				HorizonRight		= itsHorizonTerms[i+1];

				for (int j=0; j<exactSeasonTerms.size(); j++)
					exactSeasonTerms[j] = LinearInterpolation<double>(
											interpolHorizon,  
											HorizonLeft,
											itsSeasonTermData.Elt(j,i),
											HorizonRight,
											itsSeasonTermData.Elt(j,i+1));

				seasonBump	= exactSeasonTerms[interpolDateMonth-1];

				if (interpolHorizon == ref1Horizon)
				{
					seasonBump1 = exactSeasonTerms[refDate1Month-1];
				}
				else
				{
					if (ref1Horizon >= itsHorizonTerms[itsHorizonTerms.size()-1])
					{
						seasonBump1	= itsSeasonTermData.Elt(refDate1Month-1,itsHorizonTerms.size()-1);
					}
					else
					{
						for (int j =0; j< itsHorizonTerms.size()-1; j++)
						{
							if(ref1Horizon >= itsHorizonTerms[j] && ref1Horizon <= itsHorizonTerms[j+1])
							{
								/// Bounding the interpolHorizon
								HorizonLeft  = itsHorizonTerms[j];
								HorizonRight = itsHorizonTerms[j+1];

								seasonBump1 = LinearInterpolation<double>(
														ref1Horizon,  
														HorizonLeft,
														itsSeasonTermData.Elt(refDate1Month-1,j),
														HorizonRight,
														itsSeasonTermData.Elt(refDate1Month-1,j+1));

								j = itsHorizonTerms.size()-1;
							}
						}
					}
				}

				if (interpolHorizon == ref2Horizon)
				{
					seasonBump2 = exactSeasonTerms[refDate2Month-1];
				}
				else
				{
					if (ref2Horizon >= itsHorizonTerms[itsHorizonTerms.size()-1])
					{
						seasonBump2	= itsSeasonTermData.Elt(refDate2Month-1,itsHorizonTerms.size()-1);
					}
					else
					{
						for (int j =0; j< itsHorizonTerms.size()-1; j++)
						{
							if(ref2Horizon >= itsHorizonTerms[j] && ref2Horizon <= itsHorizonTerms[j+1])
							{
								/// Bounding the interpolHorizon
								HorizonLeft  = itsHorizonTerms[j];
								HorizonRight = itsHorizonTerms[j+1];

								seasonBump2 = LinearInterpolation<double>(
														ref2Horizon,
														HorizonLeft,
														itsSeasonTermData.Elt(refDate2Month-1,j),
														HorizonRight,
														itsSeasonTermData.Elt(refDate2Month-1,j+1));

								j = itsHorizonTerms.size()-1;
							}
						}
					}
				}
				i = itsHorizonTerms.size()-1;
			}
		}

		//seasonBump2 = exactSeasonTerms[refDate2Month-1];
		linCorrect	= LinearInterpolation<double>(interpolDateLag, refDate1Lag, seasonBump1, refDate2Lag, seasonBump2);
	}

	switch( itsMode )
	{
	case K_ADDITIVE:
		CPIResult += seasonBump - linCorrect;
		break;
	case K_MULTIPLICATIVE:
		CPIResult *= seasonBump / linCorrect;
		break;
	default:
		throw Exception(__LINE__, __FILE__, ERR_INVALID_DATA, "Unknown correction mode for the manager, supported is additive and multiplicative!" );
	}
}

CC_END_NAMESPACE()

/*---------------------------------------------------------------*/
/*---- End Of File ----*/

