/*
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 *
 * $Log: correlmanager.h,v $
 * Revision 1.9  2004/02/17 14:11:16  mcampet
 * MC add validation at CorrelMatrix level
 *
 * Revision 1.8  2004/02/13 14:37:22  mcampet
 * MC add ARMCorrelMatrix builders
 *
 * Revision 1.7  2003/11/07 11:43:01  ebenhamou
 * added default correl
 *
 * Revision 1.6  2003/09/29 08:19:28  ebenhamou
 * using macro CC_USING_NS
 *
 * Revision 1.5  2003/09/26 17:05:31  ebenhamou
 * version with namespace handled by macro
 *
 * Revision 1.4  2003/09/22 18:08:25  ebenhamou
 *  added NOCHECK functionality
 *
 * Revision 1.2  2003/09/19 18:18:58  ebenhamou
 * more validation
 *
 * Revision 1.1  2003/09/17 18:23:47  ebenhamou
 * Initial revision
 *
 *
 *
 */

    
/*----------------------------------------------------------------------------*/

/*! \file correlmanager.h
 *
 *  \brief Correlation manager
 *
 *	\author  Eric Benhamou
 *	\version 1.0
 *	\date September 2003
 */

/*----------------------------------------------------------------------------*/


#ifndef _CORRELMANAGER_H
#define _CORRELMANAGER_H

#include <firsttoinc.h>
#include <volint.h>
#include <volflat.h>
#include <volcube.h>
#include <string>
#include <map>
#include "IndexIndexCorrelCube.h"

class ARM_Date;
class ARM_Vector;
class ARM_Matrix;
class ARM_IndexIndexCorrelCube;
class ARM_HyperCube;

/*!
 * the Solaris compiler seems to be using std namespace
 * so the using directive are not necessary
 * using std::map; 
 * using std::pair; 
 * using std::string; 
 * using std::make_pair;
 */

using std::map;
using std::string;
using std::less;
using std::vector;


/*! \class   ARM_CorrelMatrix
 *	\brief  object that looks like an ARM_VolCurv more precisely an ARM_VolLInterpol
 * but with meaningful name and specific view. 
 * It basically wraps an ARM volcurve to become a correlmatrix
 * all the function names are defiend with more meaningful name to avoid 
 * manipulating correlation matrix still using vol curve convention
 *
 * the object is defined as a struct since everything is public and
 * it has not other data than an ARM_VolInterpol
 *
 *	\author  Eric Benhamou
 *	\version 1.0
 */
struct ARM_CorrelMatrix : public ARM_VolLInterpol
{
	void Init();
	ARM_CorrelMatrix( const ARM_Date& asOf, ARM_Vector* X, ARM_Vector* Y, ARM_Matrix* Z );
	ARM_CorrelMatrix( const ARM_CorrelMatrix& rhs );
    ARM_CorrelMatrix( const ARM_VolLInterpol* correl );
    ARM_CorrelMatrix( const ARM_VolFlat* correl );
	ARM_CorrelMatrix( const ARM_VolCube* correl );
	ARM_CorrelMatrix& operator =( const ARM_CorrelMatrix& rhs );

	/// conversion operator from a volint
	ARM_CorrelMatrix( const ARM_VolLInterpol& rhs );

	//// standard ARM Object support
	virtual ~ARM_CorrelMatrix();
	void BitwiseCopy(const ARM_Object* src);
	virtual void Copy(const ARM_Object* src);
	virtual ARM_Object* Clone();
	virtual void View(char* id = NULL, FILE* ficOut = NULL);
	double ComputeCorrelData( double x, double y ) const;
	ARM_Matrix* GetCorrelData( ) const;
	ARM_Vector* GetX() const;
	ARM_Vector* GetY() const;
};



struct CorrelDataDefault
{
	static double defaultCorrel;
};


/*! \class   ARM_CorrelManager
 *	\brief  object that is responsible for providing various correlation information
 *	It is basically a dictionary of Correlation Data
 *  with two levels of tags: The first one is the market tag while the second one
 *  specifies the intra-market information
 *  the mktTag is basically mkt1/mkt2 in alphabetical order! (checked when validating)
 *	the intraMktTag can only be of the form
 *		- ccy1_Tenor1_Ccy2_DIAG
 *		- ccy1_Tenor1_Ccy2_Tenor2_COMP
 *		- ccy1_Tenor1_Ccy2_Tenor2_NOCHECK
 *		- ccy1/ind1_ccy2/ind2
 *
 * The validation of the data has the following rules
 * 1) if the tag already exists, throw an exception...does not allow overwritting
 *
 * 2) each correl data has to be between -1 and 1 (up to the resize factor which is defined in CC_NS( ARM_Constants, correlBase )
 *
 * 3) matrix in complet mode has to be 
 *			- general
 *				- squared
 *				- same X and Y
 *			- if using the same subtags
 *				- symetric 
 *				- with the diagonal with 1.0 (up to the resize factor which is defined in CC_NS( ARM_Constants, correlBase ) )
 * 4) mod DIAG has to be consistent with the corresponding COMP mod!
 *
 * 5) mktTag that has to be in alphabetical order!
 * DIAG mod is for intra Market correlation and defines with
 * intraMktTag ccy1/ind1_Tenor1_Ccy2_DIAG and is for a given index defined by its ccy1/ind1 and tenor1
 *
 * 6) regex on the intraMktTag. Could be of the form
 *		- ccy1/ind1_Tenor1_Ccy2_DIAG
 *		- ccy1_Tenor1_Ccy2_Tenor2_COMP
 *		- ccy1_Tenor1_Ccy2_Tenor2_NOCHECK
 *		- ccy1/ind1_ccy2/ind2
 *
 *
 *	\author  Eric Benhamou
 *	\version 1.0
 */
class ARM_CorrelManager: public ARM_Object
{
public:
	typedef map< string, ARM_CorrelMatrix, less< string > >			intraMarketCorrelsMap;
	typedef map< string, intraMarketCorrelsMap, less< string > >	AllMarketCorrelsMap;
private:
	AllMarketCorrelsMap itsMktData;
	
	void ValidateData( const string& mktTag, const string& intraMktTag, const ARM_CorrelMatrix& correlMatrix ) const;
//	void Fill( const string& tmpMktTag, const string& tmpIntraMktTag, const ARM_CorrelMatrix& correlMatrix );

	intraMarketCorrelsMap GetIntraMarketCorrels( const string& mktTag ) const;
	ARM_CorrelMatrix* GetIntraMarketCorrel( const string& tmpMktTag, const string& tmpIntraMktTag ) const;
	ARM_Date GetAsOfDate() const;

	/// to avoid code duplication
	void CleanUp( );
	void CopyNoCleanUp( const ARM_CorrelManager& rhs );

	/// this enum defines the validation done by the correlation manager
	/// NOCHECK does no check as it sounds
	/// DIAG/COMP has to be consistent with other data
	/// by default, the matType is initialised to NOCHECK!
	enum MatType
	{
		COMP = 0,
		DIAG,
		NOCHECK
	};

    MatType ValidateCorrelMatrix( const string& mktTag, const string& intraMktTag, const ARM_CorrelMatrix& correlMatrix ) const;
	
public:

	void Fill( const string& tmpMktTag, const string& tmpIntraMktTag, const ARM_CorrelMatrix& correlMatrix );

	ARM_CorrelManager( const string& mktTag, const string& intraMktTag, const ARM_Date& asOf, ARM_Vector* X, ARM_Vector* Y, ARM_Matrix* Z );
	ARM_CorrelManager( const string& mktTag, const string& intraMktTag, ARM_CorrelMatrix* volInt );

	ARM_CorrelManager( const string& mktTag, const string& intraMktTag, ARM_VolFlat* volFlat );
	ARM_CorrelManager( const string& mktTag, const string& intraMktTag, ARM_VolLInterpol* volInterp );
	ARM_CorrelManager( const string& mktTag, const string& intraMktTag, ARM_VolCube* volCube );

	//// the object can be created and filled after creation
	void Fill( const string& mktTag, const string& intraMktTag, ARM_Vector* X, ARM_Vector* Y, ARM_Matrix* Z );
	void Fill( const string& mktTag, const string& intraMktTag, ARM_CorrelMatrix* volInt );

    void Init(void);
	ARM_CorrelManager( const ARM_CorrelManager& rhs );
	ARM_CorrelManager& operator= ( const ARM_CorrelManager& rhs );

	//// standard ARM Object support
	ARM_CorrelManager(){};
	virtual ~ARM_CorrelManager();

    void BitwiseCopy(const ARM_Object* src);
	virtual void Copy(const ARM_Object* src);
	virtual ARM_Object* Clone();
	virtual void View(char* id = NULL, FILE* ficOut = NULL);

	AllMarketCorrelsMap GetMktData() const { return itsMktData;}
	/// functions to get the data
	ARM_CorrelMatrix* ComputeCorrelData( const string& mktTag, const string& intraMktTag ) const;

	ARM_CorrelMatrix* GetCorrelData( const string& mktTag, const string& intraMktTag ) const;
    
    virtual double ComputeCorrelData(const string& mktTag, const string& intraMktTag, double x, double y ) const;

    virtual double ComputeIndexIndexCorrelData(	char* strccy1, const string& index1Name, 
												char* strccy2, const string& index2Name,
												double x, double y, bool aAcceptLiborMode = false);

	void BumpVolatility(const vector<string> mktTag, const vector<string> intraMktTag, 
						double value, int nthLine = 0, int nthCol = 0,
						int isCumul = K_NO, int isAbsolute = K_YES);
};

/*******************************************************************
///////////////////////////////////////////////////////////////////*
///////////////////////////////////////////////////////////////////*
///////////////////////////////////////////////////////////////////*
///////////////////////////////////////////////////////////////////*
///////////// Class ARM_CorrelatorManager 	///////////////////////*
///////////////////////////////////////////////////////////////////*
///////////////////////////////////////////////////////////////////*
///////////////////////////////////////////////////////////////////*
/*******************************************************************
*/

class ARM_CorrelatorManager : public ARM_CorrelManager
{
public:
	enum MapType
	{
		DiagCorr = 0,
		IndexIndexCorr,
		Corr,
		IdCorr,
		IRVolCorr,
		VolVolCorr,
		FXVolCorr
	};


	ARM_CorrelatorManager();


	virtual ~ARM_CorrelatorManager();

	void Init();

	ARM_CorrelatorManager( const ARM_CorrelatorManager& rhs );

	ARM_CorrelatorManager& operator= ( const ARM_CorrelatorManager& src );

// FIXMEFRED: mig.vc8 (21/05/2007 10:58:59): missing return type
	void Add( const string& name, const ARM_Object* src, MapType mtype);

	ARM_IndexIndexCorrelCube* GetCorrelIndexIndexCube( const string ccy ) const;
	
	ARM_HyperCube* GetCorrelHyperCube( const string ccy ) const;

	ARM_VolCurve* GetCorrModeCorrelCurve( const string name ) const;

	ARM_VolCurve* GetSimpleModeCorrelCurve( const string name ) const;

	ARM_HyperCube* GetIRVolCorrelHyperCube( const string ccy ) const;

	ARM_HyperCube* GetVolVolCorrelHyperCube( const string ccy ) const;

	ARM_VolCube* GetFXCorrelCube( const string ccy ) const;


	void CleanUp();

	void BitwiseCopy( const ARM_Object* src);
	
	virtual void Copy(const ARM_Object* src);

	virtual ARM_Object* Clone();

	virtual void View(char* id = NULL, FILE* ficOut = NULL);

	double ComputeIndexIndexCorrel(const string& ccy, const string& tenor1, const string& tenor2,
                                   double expiry1, double expiry2) const;
	
	double ComputeHyperCorrel( const string& ccy, const string& tenor1, const string& tenor2,
							   double expiry, double strike = 0.0) const;

	double ComputeCorrelByExpiry(const string& ccy, const string& tenor1, const string& tenor2, double expiry) const;

	double ComputeCorrModeCorrel( const string& ccy, const string& tenor, double expiry, double strike = 0.0 ) const;

	double ComputeSimpleModeCorrel( const string& ccy, const string& tenor, double expiry, double strike = 0.0 ) const;

	virtual double ComputeCorrelData(const string& mktTag, const string& intraMktTag, double x, double y ) const;

    virtual double ComputeIndexIndexCorrelData(	char* strccy1, const string& index1Name, 
												char* strccy2, const string& index2Name,
												double x, double y, bool aAcceptLiborMode = false);

	double ComputeIRVolCorrel(const string& ccy, const string& tenor1, const string& tenor2, double expiry) const;

	double ComputeVolVolCorrel(const string& ccy, const string& tenor1, const string& tenor2, double expiry) const;

	double ComputeFXCorrel(const string& ccy, const string& ccy2, const string& tenor, double expiry) const;


	void BumpVolatility(const string& ccy, MapType type, double value, int nthLine = 0, int nthCol = 0,
                            int isCumul = K_NO, int isAbsolute = K_YES);

	void BumpCorrel(const string& ccy, MapType type, double value, int nthLine = 0, int nthCol = 0,
                    int isCumul = K_NO, int isAbsolute = K_YES);


private:
	map< string , ARM_IndexIndexCorrelCube* > itsIndexData;
	map< string , ARM_HyperCube* >	itsDiagData;
	map< string , ARM_VolCurve* >	itsCorrData;
	map< string , ARM_VolCurve* >	itsSimpleModeData;
	map< string , ARM_HyperCube* >	itsIRVolData;
	map< string , ARM_HyperCube* >	itsVolVolData;
	map< string , ARM_VolCube* >  itsVolFXData;
};

#endif
/*----------------------------------------------------------------------------*/
/*---- End of file ----*/
