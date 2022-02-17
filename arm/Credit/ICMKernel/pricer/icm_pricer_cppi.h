
#ifndef _ICM_PRICER_CPPI_H
#define _ICM_PRICER_CPPI_H

#include "ICMKernel/glob/icm_enums.h"
#include "ICMKernel\pricer\icm_pricer.h"

/*********************************************************************************/
/*! \class  ICM_Pricer_Cppi icm_pricer_cppi.h "icm_pricer_cppi.h"
 *  \author Damien Pouponneau
 *	\version 1.0
 *	\date   April 2005
 *	\file   icm_pricer_cppi.h
 *	\brief  Pricer for CPPI  */
/***********************************************************************************/

class ICM_Pricer_Cppi : public ICM_Pricer
{
private:
	// Rendement placement
	double	itsFunding;
	double  itsFinancement;
	double	itsRunning;
										
	// Processus de poisson 
	double itsLambda;

	// Calcul de la soult
	double itsBidMid;
	double itsSwitch;
	
	// Modele Vasicek
	double	itsR0;
	double itsIRVol;
	double itsA;

	// Monte Carlo Data
	int itsNbSimul;
	double itsStep;
	int itsIndexIndice;

	// Vol Implied Mez
	double itsVolImplied;
	double itsMeanImplied;
	double itsSpeedReversionImplied;

	vector<double> itsZeroCoupon;


	// Coupon et prime upfront
	double itsSpreadCoupon;
	double itsPrimeUpFront;
	double itsFactorCoupon;
	vector<double> itsNbCouponPayed;
	double itsStopMaturity;
	int itsFixedEurib;

	// BarriereDown
	double itsLockedTime;
	double itsBarriereDown;

	// Widener
	double itsInitWidener;
	double itsMuWidener;
	double itsSigmaWidener;
	double itsPartWidener;

	// Data for Test
	vector<double> itsZImplied;
	vector<double> itsImplied;
	vector<double> itsE;
	vector<double> itsf;
	vector<double> itsN;
	vector<double> itsV;
	vector<double> itsP;
	vector<double> itsAlphaMin;
	vector<double> itsAlphaMax;
	vector<double> itsAlpha;
	vector<int> itsNbRebal;
	vector<double> itsDur;
	vector<double> itsSoult;
	vector<double> itsR;
	vector<double> itsS;
	vector<double> itsTimeToRoll;
	vector<double> itsJump;
	vector<double> itsSaut;
	vector<double> itsSigmaV;	// itsSigmaV[i] = Vol de V apres la ième année.
	vector<double> itsMeanV;	// itsMeanV[i] = Mean de V apres la ième année.
	vector<double> itsWidener;
	vector<double> itsBarriereDownTime;

	ICM_QMatrix<double>* itsDistribTRI;

public :
	ICM_Pricer_Cppi(void) { Init();}

	/** ICM_Pricer_Cppi(ARM_Security* sec, ARM_Object* datamkt,const ICM_Parameters&parameters,const ARM_Date&asof)
	{
		Init();
		Set(sec, datamkt,parameters,asof);
	} **/ 

	~ICM_Pricer_Cppi()
	{
		itsZImplied.clear();
		itsImplied.clear();
		itsE.clear();
		itsf.clear();
		itsN.clear();
		itsV.clear();
		itsP.clear();
		itsAlphaMin.clear();
		itsAlphaMax.clear();
		itsAlpha.clear();
		itsNbRebal.clear();
		itsZeroCoupon.clear();
		itsDur.clear();
		itsSoult.clear();
		itsR.clear();
		itsS.clear();
		itsTimeToRoll.clear();
		itsJump.clear();
		itsSaut.clear();
		itsSigmaV.clear();
		itsMeanV.clear();
		itsWidener.clear();
		itsNbCouponPayed.clear();
		itsBarriereDownTime.clear();

		if (itsDistribTRI)
			delete itsDistribTRI;
		itsDistribTRI = NULL;
	}

	void Init() ;
	vector<double> GetImplied()			{ return itsImplied;}
	vector<double> GetE()				{ return itsE;}
	vector<double> Getf()				{ return itsf;}
	vector<double> GetN()				{ return itsN;}
	vector<double> GetV()				{ return itsV;}
	vector<double> GetP()				{ return itsP;}
	vector<double> GetAlphaMin()		{ return itsAlphaMin;}
	vector<double> GetAlphaMax()		{ return itsAlphaMax;}
	vector<double> GetAlpha()			{ return itsAlpha;}
	vector<double> GetSigmaV()			{ return itsSigmaV;}
	vector<double> GetMeanV()			{ return itsMeanV;}
	vector<double> GetWidener()			{ return itsWidener;}
	vector<double> GetBarriereDownTime(){ return itsBarriereDownTime;}
	vector<double>	   GetNbCouponPayed()	{ return itsNbCouponPayed;}
	
	int GetRowDistrib()		{ return itsDistribTRI->Getnbrows();}
	int GetColDistrib()		{ return itsDistribTRI->Getnbcols();}
	int GetNbRebal(int k)	{ if (k<itsNbRebal.size())	{return itsNbRebal[k];} else { return -999.;}}
	
	double GetZImplied(int k)			{ if (k<itsZImplied.size()) {return itsZImplied[k];} else { return -999.;}}
	double GetImplied(int k)			{ if (k<itsImplied.size())	{return itsImplied[k];} else { return -999.;}}
	double GetE(int k)					{ if (k<itsE.size())		{return itsE[k];} else { return -999.;}}
	double Getf(int k)					{ if (k<itsf.size())		{return itsf[k];} else { return -999.;}}
	double GetN(int k)					{ if (k<itsN.size())		{return itsN[k];} else { return -999.;}}
	double GetV(int k)					{ if (k<itsV.size())		{return itsV[k];} else { return -999.;}}
	double GetP(int k)					{ if (k<itsP.size())		{return itsP[k];} else { return -999.;}}
	double GetAlphaMin(int k)			{ if (k<itsAlphaMin.size()) {return itsAlphaMin[k];} else { return -999.;}}
	double GetAlphaMax(int k)			{ if (k<itsAlphaMax.size()) {return itsAlphaMax[k];} else { return -999.;}}
	double GetAlpha(int k)				{ if (k<itsAlpha.size())	{return itsAlpha[k];} else { return -999.;}}; 
	double GetDur(int k)				{ if (k<itsDur.size())		{return itsDur[k];} else { return -999.;}};
	double GetSoult(int k)				{ if (k<itsSoult.size())	{return itsSoult[k];} else { return -999.;}};
	double GetR(int k)					{ if (k<itsR.size())		{return itsR[k];} else { return -999.;}};
	double GetS(int k)					{ if (k<itsS.size())		{return itsS[k];} else { return -999.;}};
	double GetTimeToRoll(int k)			{ if (k<itsTimeToRoll.size()) {return itsTimeToRoll[k];} else { return -999.;}};
	double GetJump(int k)				{ if (k<itsJump.size())		{return itsJump[k];} else { return -999.;}};
	double GetSaut(int k)				{ if (k<itsSaut.size())		{return itsSaut[k];} else { return -999.;}};
	double GetSigmaV(int k)				{ if (k<itsSigmaV.size())	{return itsSigmaV[k];} else { return -999.;}}
	double GetMeanV(int k)				{ if (k<itsMeanV.size())	{return itsMeanV[k];} else { return -999.;}}
	double GetWidener(int k)			{ if (k<itsWidener.size())	{return itsWidener[k];} else { return -999.;}}
	double GetNbCouponPayed(int k)			{ if (k<itsNbCouponPayed.size())	{return itsNbCouponPayed[k];} else { return -999.;}}
	double GetDistribTRI(int i,int j)		{ return (*itsDistribTRI)(i,j);}
	double GetBarriereDownTime(int k)	{ if (k<itsBarriereDownTime.size()) {return itsBarriereDownTime[k];} else { return -999.;}};

	void Set(ARM_Security *sec, ARM_Object *mod,const ICM_Parameters &parameters,const ARM_Date&asof);
	
	virtual void Reset(void) {}
	// Faire la methode de RESET

	virtual double ComputePrice(qCMPMETH measure) ;

	virtual void MarketData(ARM_Security* sec,vector<string>& DMKT){}

protected:
	virtual double ComputeSensitivity(const qSENSITIVITY_TYPE& typesensi, 
		const std::string& plot,
		const std::string& label,
		double epsilon, double epsilonGamma =0 ) ; 

private:
	virtual double Accrued() ; 
	virtual double FeeLegPV () ;
	virtual double DefLegPV () ;
	virtual double ComputeDuration(void) ;
	virtual double ComputeSpread(const double& MtM = 0.) ;
 	virtual	double ComputeImpliedVol(const double& Price) ;
	virtual double Compute_Fwd_Spread(const ARM_Date &Mty, const ARM_Date &CDS_ExpiryDate,double& dur) ;
};

#endif

