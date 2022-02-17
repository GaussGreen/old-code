/*
 *
 * Copyright (c) CDC IXIS CM November 2003 Paris
 *
 * $Log: seasonmanager.h,v $
 * Revision 1.2  2003/11/25 16:06:30  rguillemot
 * Inflation Seasonality Manager
 *
 * Revision 1.1  2003/11/25 10:16:13  rguillemot
 * Initial revision
 *
 *
 */

/*----------------------------------------------------------------------------*
    seasonmanager.h
	Copyright (c) CDC IXIS CM November 2003 Paris
*----------------------------------
------------------------------------------*/

/*----------------------------------------------------------------------------*/

/*! \file seasonmanager.h
 *
 *  \brief Seasonality manager
 *
 *	\author  Richard GUILLEMOT
 *	\version 1.0
 *	\date November 2003
 */

/*----------------------------------------------------------------------------*/

#ifndef _INGPINFLATION_SEASONMANAGER_H
#define _INGPINFLATION_SEASONMANAGER_H

#include "firsttoinc.h"
#include "gpbase/port.h"

/// gpbase
#include "gpbase/gplinalgtypedef.h"
#include "gpbase\gpmatrix.h"

/// kernel
#include <glob/armglob.h>

/// STL
#include <string>
CC_USING_NS( std, string )
#include <vector>
CC_USING_NS( std, vector )

class ARM_Date;

/*!
 * the Solaris compiler seems to be using std namespace
 * so the using directive are not necessary
 * using std::string; 
 * using std::vector;
 */

CC_BEGIN_NAMESPACE( ARM )


/*! \class   ARM_SeasonalityManager
 *	\brief  This object that is responsible for providing various season adjustment to add on the 
 *	interplated CPI values
 *  It is basicaly a vector of season spreads corresponding to each month. It also contains a vector
 *	of boolean to check if no data is missing or inputed twice. The method used by the inflation curve
 *  to compute the season adjustment is ComputeSeasonSpread.
 *	\author  Richard GUILLEMOT
 *	\version 1.0
 */

 class ARM_SeasonalityManager : public ARM_Object
{
public:
	/// enum to specify the mode
	enum CorrectionMode
	{
		K_ADDITIVE		= K_SEASONADJ_ADDITIVE, 
		K_MULTIPLICATIVE= K_SEASONADJ_MULTIPLICATIVE
	};

private:
	vector<double> itsSeasonData;
	ARM_GP_Matrix  itsSeasonTermData;
	vector<double> itsHorizonTerms;
	std::deque<bool> /* FIXMEFRED: mig.vc8 (23/05/2007 14:27:19): dont use std::vector<bool>. see item 18 in "Effective STL" from Scott Meyers*/ itsHasBeenStored;
	CorrectionMode itsMode;

	/// Validation of the input data of the seasonality manager
	void ValidateData(const vector<string>& monthList, const ARM_GP_Vector& seasonSpread, const vector<double>& horizonList);
	// Test if all the season data have been inputed
	void CheckDataIntegrity();
	/// Init Method
	void Init(const vector<double>& horizonList);
	/// Factorize Method
	void Factorize(const vector<double>& horizonList, const vector<string>& monthList, const ARM_GP_Vector& seasonSpreadList);

public:
	/// Main constructor
	ARM_SeasonalityManager(const vector<string>& monthList, const ARM_GP_Vector& seasonSpreadList, CorrectionMode mode=K_ADDITIVE);
	ARM_SeasonalityManager(const vector<string>& monthList, const ARM_GP_Vector& seasonSpreadList, const vector<double>& horizonList, CorrectionMode mode=K_ADDITIVE);
	ARM_SeasonalityManager(const ARM_SeasonalityManager& rhs);
	virtual ~ARM_SeasonalityManager();
	ARM_SeasonalityManager& operator=(const ARM_SeasonalityManager &rhs);

	//// standard ARM Object support
	virtual ARM_Object* Clone();
	virtual void View(char* id = NULL, FILE* ficOut = NULL);

	/// Add a new season spread in the manager for a corresponding month
	void AddMonthSeasonSpread(const string& month, double SeasonSpread, int horizon);
	/// Correct the CPI result by the seasonality
	void SeasonSpreadCorrect( ARM_Date asOfDate, double& CPIResult, double interpolDateLag, double refDate1Lag, double refDate2Lag) const;
};

CC_END_NAMESPACE()

#endif
