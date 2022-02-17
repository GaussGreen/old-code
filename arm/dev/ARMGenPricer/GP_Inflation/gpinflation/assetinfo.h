/*
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 * $Log: assetinfo.h,v $
 * Revision 1.1  2004/06/08 07:52:17  ebenhamou
 * Initial revision
 *
 */

/*! \file assetinfo.h
 *
 *  \brief
 *	\author  E. Benhamou
 *	\version 1.0
 *	\date June 2004
 */


#ifndef _INGPINFLATION_ASSETINFO_H
#define _INGPINFLATION_ASSETINFO_H

#include "gpbase/port.h"	/// port.h
#include "gpbase/env.h"

/// kernel
#include <glob/armglob.h>

#include <string>
CC_USING_NS(std,string)
#include <vector>
CC_USING_NS(std,vector)



CC_BEGIN_NAMESPACE( ARM )


/// to keep track of asset information
/// single asset
class ARM_SingleAssetInfo
{
private:
	string itsName;
	int itsRcvOrPay;
	int itsType;
	double itsFwd;
	double itsVol;
	double itsTenor;
	double itsExpiry;
	bool itsHasBeenComputed;
	double itsStrikeForVol;

public:
	ARM_SingleAssetInfo( const string& name, int rcvOrPay, int type )
	:	itsName(name), itsRcvOrPay( rcvOrPay ), itsType( type ), itsFwd( 0.0 ), itsVol( 0.0), itsTenor( 0.0), itsExpiry( 0.0 ), itsHasBeenComputed( false ), itsStrikeForVol(0.0)
	{}

	string toString( const string& indent = "" ) const;
	inline void SetTenor( double tenor ){ itsTenor = tenor; }
	inline void SetExpiry( double expiry){ itsExpiry= expiry; }
	inline void SetFwdAndVol( double fwd, double vol ){ itsFwd = fwd; itsVol = vol;}
	inline void SetVol( double vol ){ itsVol = vol; itsHasBeenComputed = true; }
	inline void SetFwd( double fwd ){ itsFwd = fwd;}
	inline void SetStrikeForVol(double strike) { itsStrikeForVol=strike; }
	bool isFixedAsset( ) const;
	
	inline bool HasBeenComputed( ) const { return itsHasBeenComputed; }
	inline int GetType() const { return itsType; }
	inline double GetTenor() const { return itsTenor; }
	inline double GetExpiry() const { return itsExpiry; }
	inline double GetVol() const { return itsVol; }
	inline double GetFwd() const { return itsFwd; }
	inline double GetStrikeForVol() const { return itsStrikeForVol; }
};


/// two assets information
class ARM_TwoAssetsInfo : public ARM_Object
{
private:
	ARM_SingleAssetInfo itsFirstAsset;
	ARM_SingleAssetInfo itsSecondAsset;
	double itsCorrelation;
	double itsAnnuity;
	double itsOptionMaturity;
	double itsStrike;
	double itsPricingStrike;
	int	   itsInfType;
	double itsFirstAssetCoupon; 
	

public:
	ARM_TwoAssetsInfo( const ARM_SingleAssetInfo& infAsset, 
		const ARM_SingleAssetInfo& secondAsset, double strike )
	:	
		itsFirstAsset(infAsset), itsSecondAsset(secondAsset), 
		itsCorrelation(0.0), itsAnnuity(0.0), itsOptionMaturity(0.0),itsStrike( strike ), itsPricingStrike( strike ) 
	{}

	ARM_TwoAssetsInfo( const ARM_TwoAssetsInfo& rhs )
	:
		ARM_Object( rhs ),
		itsFirstAsset( rhs.itsFirstAsset ),
		itsSecondAsset( rhs.itsSecondAsset ),
		itsCorrelation( rhs.itsCorrelation ),
		itsAnnuity( rhs.itsAnnuity ),
		itsOptionMaturity( rhs.itsOptionMaturity ),
		itsStrike( rhs.itsStrike )
	{}

	ARM_TwoAssetsInfo& operator=(const ARM_TwoAssetsInfo& rhs )
	{
		if( this != &rhs )
		{
			ARM_Object::operator =( rhs );
			itsFirstAsset		= rhs.itsFirstAsset;
			itsSecondAsset		= rhs.itsSecondAsset;
			itsCorrelation		= rhs.itsCorrelation;
			itsAnnuity			= rhs.itsAnnuity;
			itsOptionMaturity	= rhs.itsOptionMaturity;
			itsStrike			= rhs.itsStrike;
		}
	};
	virtual ~ARM_TwoAssetsInfo() {};

	string toString( const string& indent = "" ) const;

	/// standard ARM support
	virtual ARM_Object* Clone();
	virtual void View(char* id = NULL, FILE* ficOut = NULL);



	/// the various set function!
	void SetFwdAndVolAsset1( double fwd, double vol ){ itsFirstAsset.SetFwdAndVol( fwd, vol ); }
	void SetFwdAndVolAsset2( double fwd, double vol ){ itsSecondAsset.SetFwdAndVol( fwd, vol );} 
	void SetVolAsset1( double vol ){ itsFirstAsset.SetVol( vol ); }
	void SetVolAsset2( double vol ){ itsSecondAsset.SetVol( vol ); }
	void SetTenorAsset1( double tenor ){ itsFirstAsset.SetTenor( tenor ); }
	void SetTenorAsset2( double tenor ){ itsSecondAsset.SetTenor( tenor ); }
	void SetExpiryAsset1( double expiry ){ itsFirstAsset.SetExpiry( expiry ); }
	void SetExpiryAsset2( double expiry ){ itsSecondAsset.SetExpiry( expiry ); }
	void SetCorrelation( double correlation ) { itsCorrelation=correlation; }
	void SetAnnuity( double annuity ) { itsAnnuity=annuity; }
	void SetOptionMaturity( double OptionMaturity ) { itsOptionMaturity=OptionMaturity; }
	void SetPricingStrike( double pricingStrike ) { itsPricingStrike=pricingStrike;}
	void SetFirstAssetType( int firstAssetType ) { itsInfType = firstAssetType;}
	void SetFirstAssetCoupon( double firstAssetCoupon) { itsFirstAssetCoupon = firstAssetCoupon;}
	void SetFwdAsset1( double fwd){ itsFirstAsset.SetFwd(fwd); }
	void SetFwdAsset2( double fwd){ itsSecondAsset.SetFwd(fwd); }
	void SetStrikeForVol1(double strike) { itsFirstAsset.SetStrikeForVol(strike); }
	void SetStrikeForVol2(double strike) { itsSecondAsset.SetStrikeForVol(strike); }
	
	/// get function
	double GetOptionMaturity() const { return itsOptionMaturity; }
	double GetTenorAsset1() const { return itsFirstAsset.GetTenor(); }
	double GetTenorAsset2() const { return itsSecondAsset.GetTenor(); }
	double GetExpiryAsset1() const { return itsFirstAsset.GetExpiry(); }
	double GetExpiryAsset2() const { return itsSecondAsset.GetExpiry(); }
	double GetVolAsset1() const { return itsFirstAsset.GetVol(); }
	double GetVolAsset2() const { return itsSecondAsset.GetVol(); }
	double GetCorrelation() const { return itsCorrelation; }
	double GetFwdAsset1() const { return itsFirstAsset.GetFwd(); }
	double GetFwdAsset2() const { return itsSecondAsset.GetFwd(); }
	double GetStrikeForVol1() const { return itsFirstAsset.GetStrikeForVol(); }
	double GetStrikeForVol2() const { return itsSecondAsset.GetStrikeForVol(); }
	

	
	
	/// function to know if it has been computed
	bool HasBeenComputedTwice()  const { return itsFirstAsset.HasBeenComputed() && itsSecondAsset.HasBeenComputed(); }
	bool HasBeenComputed1()  const { return itsFirstAsset.HasBeenComputed();}
	bool HasBeenComputed2() const { return itsSecondAsset.HasBeenComputed();}
	int GetOptionType()		const { return itsSecondAsset.GetType(); }
	int GetFirstAssetType() const { return itsInfType; }
	double GetFirstAssetCoupon() const { return itsFirstAssetCoupon; }
	double GetPricingStrike() const { return itsPricingStrike; }
};


CC_END_NAMESPACE()

#endif

/*----------------------------------------------------------------------------*/
/*---- End of file ----*/



