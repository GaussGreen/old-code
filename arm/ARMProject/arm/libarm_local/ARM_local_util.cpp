#pragma warning(disable : 4541)

#include "firstToBeIncluded.h"
#include "arm_local_glob.h"
#include "arm_local_zccurve.h"
#include "arm_local_persistent.h"
#include "arm_local_mod.h"
#include "arm_local_util.h"

#include <ARM\libarm_frometk\ARM_local_parsexml.h>
#include <ARM\libarm_frometk\ARM_local_parsexml_nt_pf.h>
#include <ARM\libarm_frometk\ARM_local_etoolkit.h>
#include <ARM\libarm_local\ARM_local_volcrv.h>
#include <ARM\libarm_local\ARM_local_class.h>

#include <ARM\local_xlarm\ARM_local_interglob.h>

#include <inst\swap.h>
#include <inst\optionportfolio.h>
#include <inst\spreadoption.h>
#include <util\fromto.h>
#include <crv\zerointimp.h>

double ARM_MATU_COLS[MATU_COLS_SIZE] = {0.25, 0.5, 1.0, 2.0, 3.0, 4.0, 5.0,
										7.0, 10.0, 12.0, 15.0, 20.0, 30.0}; 

double ARM_MATU_LINES[MATU_LINES_SIZE] = {0.25, 0.5, 1.0, 2.0, 3.0, 4.0, 5.0,
											7.0, 10.0, 12.0, 15.0, 20.0, 30.0}; 

ARM_Swap* curSwap;
FILE* Fp;


long prepareData(ARM_Date* Date1, ARM_Date* Date2, ARM_Vector* Data1, ARM_Vector* Data2, long nbDate1, long nbDate2)
{
	ARM_Date tmpDate1 (Date1[0]);
	ARM_Date tmpDate2 (Date2[0]);

	ARM_Vector* newData1 = new ARM_Vector(Data1->GetSize());
	ARM_Vector* newData2 = new ARM_Vector(Data2->GetSize());

	long i(0), j(0);
	long nbDates (0);
	while ( (i < nbDate1) && (j < nbDate2) )
	{
		if (tmpDate1 == tmpDate2)
		{
			newData1->Elt(nbDates) = Data1->Elt(i);
			newData2->Elt(nbDates) = Data2->Elt(j);

			i++;
			j++;
			nbDates++;

			tmpDate1 = Date1[i];
			tmpDate2 = Date1[j];
		}
		else if (tmpDate1 < tmpDate2)
		{
			i++;
			tmpDate1 = Date1[i];
		}
		else
		{
			j++;
			tmpDate2 = Date2[j];
		}
	}

	Data1->Resize(nbDates);
	Data2->Resize(nbDates);

	for (int i=0;i<nbDates;i++)
	{
		Data1->Elt(i) = newData1->Elt(i);
		Data2->Elt(i) = newData2->Elt(i);
	}

	if (newData1)
		delete newData1;
	newData1 = NULL;

	if (newData2)
		delete newData2;
	newData2 = NULL;

	return 0;
}


long prepareDataMoy(ARM_Date* Date1, ARM_Date* Date2, ARM_Vector** Data1, ARM_Vector** Data2, long nbDate1, long nbDate2, long lExpiry1, long lExpiry2)
{
	ARM_Date tmpDate1 (Date1[0]);
	ARM_Date tmpDate2 (Date2[0]);

	ARM_Vector** newData1 = new ARM_Vector* [lExpiry1];
	ARM_Vector** newData2 = new ARM_Vector* [lExpiry2];
    int compteur;
	for (compteur = 0; compteur<lExpiry1; compteur++)
	{
		newData1[compteur] = new ARM_Vector(Data1[compteur]->GetSize());
	}

	for (compteur = 0; compteur<lExpiry2; compteur++)
	{
		newData2[compteur] = new ARM_Vector(Data2[compteur]->GetSize());
	}

	long i(0), j(0);
	long nbDates (0);
	while ( (i < nbDate1) && (j < nbDate2) )
	{
		if (tmpDate1 == tmpDate2)
		{
			for (compteur=0;compteur<lExpiry1;compteur++)
			{
				newData1[compteur]->Elt(nbDates) = Data1[compteur]->Elt(i);
			}
			for (compteur=0;compteur<lExpiry2;compteur++)
			{
				newData2[compteur]->Elt(nbDates) = Data2[compteur]->Elt(i);
			}
			i++;
			j++;
			nbDates++;

			tmpDate1 = Date1[i];
			tmpDate2 = Date1[j];
		}
		else if (tmpDate1 < tmpDate2)
		{
			i++;
			tmpDate1 = Date1[i];
		}
		else
		{
			j++;
			tmpDate2 = Date2[j];
		}
	}

	for (compteur=0;compteur<lExpiry1;compteur++)
	{
		Data1[compteur]->Resize(nbDates);
	}
	for (compteur=0;compteur<lExpiry2;compteur++)
	{
		Data2[compteur]->Resize(nbDates);
	}

	for (i=0;i<nbDates;i++)
	{
		for (compteur=0;compteur<lExpiry1;compteur++)
		{
			Data1[compteur]->Elt(i) = newData1[compteur]->Elt(i);
		}
		for (compteur=0;compteur<lExpiry2;compteur++)
		{
			Data2[compteur]->Elt(i) = newData2[compteur]->Elt(i);
		}
	}

	for (compteur=0;compteur<lExpiry1;compteur++)
	{
		delete newData1[compteur];
	}
	for (compteur=0;compteur<lExpiry2;compteur++)
	{
		delete newData2[compteur];
	}

	if (newData1)
		delete [] newData1;
	newData1 = NULL;

	if (newData2)
		delete [] newData2;
	newData2 = NULL;

	return 0;
}




long transformData(ARM_Date* Date1,ARM_Vector* Data1, ARM_Vector* Data2, long Type) 
{
	ARM_Vector* newData1 = new ARM_Vector(Data1->GetSize()-1);
	ARM_Vector* newData2 = new ARM_Vector(Data2->GetSize()-1);

	for (int i=0; i<Data1->GetSize()-1; i++)
	{
		// LogNormal
		if (Type == 1)
		{
	//		newDate1[i] = Date1[i+1];
			newData1->Elt(i) = log(Data1->Elt(i+1)/Data1->Elt(i));
			newData2->Elt(i) = log(Data2->Elt(i+1)/Data2->Elt(i));
		}
		// Normal
		else
		{
	//		newDate1[i] = Date1[i+1];
			newData1->Elt(i) = Data1->Elt(i+1) - Data1->Elt(i);
			newData2->Elt(i) = Data2->Elt(i+1) - Data2->Elt(i);
		}
	}

	Data1->Resize(Data1->GetSize()-1);
	Data2->Resize(Data2->GetSize()-1);

	for (int i=0; i<Data1->GetSize()-1; i++)
	{
		Data1->Elt(i) = newData1->Elt(i);
		Data2->Elt(i) = newData2->Elt(i);
	}

	delete newData1;
	delete newData2;

	return 0;
}


long transformDataMoy(ARM_Date* Date1,ARM_Vector** Data1, ARM_Vector** Data2, long Type, long lExpiry1, long lExpiry2)
{
	ARM_Vector** newData1 = new ARM_Vector* [lExpiry1];
	ARM_Vector** newData2 = new ARM_Vector* [lExpiry2];
	int compteur;
	for (compteur = 0; compteur<lExpiry1; compteur++)
	{
		newData1[compteur] = new ARM_Vector(Data1[compteur]->GetSize()-1);
	}
	for (compteur = 0; compteur<lExpiry2; compteur++)
	{
		newData2[compteur] = new ARM_Vector(Data2[compteur]->GetSize()-1);
	}

	for (compteur = 0; compteur<lExpiry1; compteur++)
	{
		for (int i=0; i<Data1[compteur]->GetSize()-1; i++)
		{
			// LogNormal
			if (Type == 1)
			{
		//		newDate1[i] = Date1[i+1];
				newData1[compteur]->Elt(i) = log(Data1[compteur]->Elt(i+1)/Data1[compteur]->Elt(i));
			}
			// Normal
			else
			{
		//		newDate1[i] = Date1[i+1];
				newData1[compteur]->Elt(i) = Data1[compteur]->Elt(i+1) - Data1[compteur]->Elt(i);
			}
		}
	}
	for (compteur = 0; compteur<lExpiry2; compteur++)
	{
		for (int i=0; i<Data2[compteur]->GetSize()-1; i++)
		{
			// LogNormal
			if (Type == 1)
			{
				newData2[compteur]->Elt(i) = log(Data2[compteur]->Elt(i+1)/Data2[compteur]->Elt(i));
			}
			// Normal
			else
			{
				newData2[compteur]->Elt(i) = Data2[compteur]->Elt(i+1) - Data2[compteur]->Elt(i);
			}
		}
	}


	for (compteur = 0; compteur<lExpiry1; compteur++)
	{
		Data1[compteur]->Resize(Data1[compteur]->GetSize()-1);

		for (int i=0; i<Data1[compteur]->GetSize(); i++)
		{
			Data1[compteur]->Elt(i) = newData1[compteur]->Elt(i);
		}
	}
	for (compteur = 0; compteur<lExpiry2; compteur++)
	{
		Data2[compteur]->Resize(Data2[compteur]->GetSize()-1);

		for (int i=0; i<Data2[compteur]->GetSize(); i++)
		{
			Data2[compteur]->Elt(i) = newData2[compteur]->Elt(i);
		}
	}

	for (compteur=0;compteur<lExpiry1;compteur++)
	{
		delete newData1[compteur];
	}
	for (compteur=0;compteur<lExpiry2;compteur++)
	{
		delete newData2[compteur];
	}

	if (newData1)
		delete [] newData1;
	newData1 = NULL;

	if (newData2)
		delete [] newData2;
	newData2 = NULL;

	return 0;
}

double calcMoy(ARM_Vector* Data1)
{
	double sum = 0;

	for (int i=0; i<Data1->GetSize(); i++)
		sum+= Data1->Elt(i);

	return sum/((double)Data1->GetSize()-1.);
}


long calcMoy2(ARM_Vector* Data1, ARM_Vector* Data2, double* moy1, double* moy2)
{
	if (Data1->GetSize() != Data2->GetSize())
		return -1;

	*moy1 = 0.;
	*moy2 = 0.;

	for (int i=0; i<Data1->GetSize(); i++)
	{
		*moy1+= Data1->Elt(i);
		*moy2+= Data2->Elt(i);
	}
	
	*moy1 /= (Data1->GetSize()-1.);
	*moy2 /= (Data2->GetSize()-1.);

	return 0;
}


long calcCorrelInst(ARM_Vector* Data1, ARM_Vector* Data2, double* correl, double* vol1, double* vol2, long* nbData)
{
	if (Data1->GetSize() != Data2->GetSize())
		return -1;

//	double var1(0.), var2(0.);
//	double moy1(0.), moy2(0.);

	double sumXi(0.), sumYi(0.), sumXi2(0.), sumYi2(0.), sumXiYi(0.);

	*correl = 0.;

//	if (calcMoy2(Data1, Data2, &moy1, &moy2) == -1)
//		return -1;

//	double tmp1, tmp2;
	for (int i=0; i<Data1->GetSize(); i++)
	{
		sumXi += Data1->Elt(i);
		sumXi2 += (Data1->Elt(i)*Data1->Elt(i));
		sumYi += Data2->Elt(i);
		sumYi2 += (Data2->Elt(i)*Data2->Elt(i));
		sumXiYi += (Data1->Elt(i)*Data2->Elt(i));

/*		tmp1 = Data1->Elt(i) - moy1;
		tmp2 = Data2->Elt(i) - moy2;
		*correl += tmp1 * tmp2;
		var1 += tmp1*tmp1;
		var2 += tmp2*tmp2;
*/	}

	*vol1 = sumXi2/(Data1->GetSize()-1.) - (Data1->GetSize()-2.)*(sumXi*sumXi)/((Data1->GetSize()-1.)*(Data1->GetSize()-1.)*(Data1->GetSize()-1.));
//	*vol1 = var1/(Data1->GetSize()-1.);
	*vol2 = sumYi2/(Data1->GetSize()-1.) - (Data1->GetSize()-2.)*(sumYi*sumYi)/((Data1->GetSize()-1.)*(Data1->GetSize()-1.)*(Data1->GetSize()-1.));
//	*vol2 = var2/(Data2->GetSize()-1.);
	*correl = sumXiYi/(Data1->GetSize()-1.) - (Data1->GetSize()-2.)*(sumXi*sumYi)/((Data1->GetSize()-1.)*(Data1->GetSize()-1.)*(Data1->GetSize()-1.));
	*nbData = Data1->GetSize();

	if ( (*vol1 == 0.) && (*vol2== 0) )
		*correl = 1.;
	else if ( (*vol1 == 0.) || (*vol2 == 0) )
		*correl = 0.;
	else
		// Pas besoin de diviser par N : simplification
		*correl /= sqrt((*vol1)*(*vol2));

	return 0.;
}

/*
long calcCorrelMoy(ARM_Vector** Data1, ARM_Vector** Data2, double* correl, double* vol1, double* vol2, long* nbData, long lExpiry1, long lExpiry2)
{
	double tmp1, tmp2;

	long lExpiry = MIN(lExpiry1,lExpiry2);

	for (int compteur=0;compteur<lExpiry;compteur++)
	{
		if (Data1[compteur]->GetSize() != Data2[compteur]->GetSize())
			return -1;
	}

	double* covar = new double[lExpiry];
	double* var1 = new double[lExpiry1];
	double* var2 = new double[lExpiry2];
	double* moy1 = new double[lExpiry1];
	double* moy2 = new double[lExpiry2];

	for (compteur=0;compteur<lExpiry;compteur++)
		covar[compteur]=0.;

	for (compteur=0;compteur<lExpiry1;compteur++)
	{
		var1[compteur]=0.;
		moy1[compteur]=0.;
	}

	for (compteur=0;compteur<lExpiry2;compteur++)
	{
		var2[compteur]=0.;
		moy2[compteur]=0.;
	}

	*correl = 0.;
	*vol1 = 0.;
	*vol2 = 0.;


	// Calcul des variances et covariances jusqu'au min des 2
	for (compteur = 0;compteur<lExpiry;compteur++)
	{
		if (calcMoy2(Data1[compteur], Data2[compteur], &(moy1[compteur]), &(moy2[compteur])) == -1)
			return -1;

		for (int i=0; i<Data1[compteur]->GetSize(); i++)
		{
			tmp1 = Data1[compteur]->Elt(i) - moy1[compteur];
			tmp2 = Data2[compteur]->Elt(i) - moy2[compteur];
			covar[compteur] += tmp1 * tmp2;
			var1[compteur] += tmp1*tmp1;
			var2[compteur] += tmp2*tmp2;
		}
//		covar[compteur] /= (Data1[compteur]->GetSize());

//		var1[compteur] /= (Data1[compteur]->GetSize());
//		var2[compteur] /= (Data2[compteur]->GetSize());
	}

	// Calcul des variances pour la série dont le ténor est le plus grand
	for (compteur = lExpiry;compteur < MAX(lExpiry1,lExpiry2);compteur++)
	{
		if (lExpiry1 > compteur)
		{
			moy1[compteur] = calcMoy(Data1[compteur]);

			for (int i=0; i<Data1[compteur]->GetSize(); i++)
			{
				tmp1 = Data1[compteur]->Elt(i) - moy1[compteur];
				var1[compteur] += tmp1*tmp1;
			}

//			var1[compteur] /= Data1[compteur]->GetSize();
		}
		if (lExpiry2 > compteur)
		{
			moy2[compteur] = calcMoy(Data2[compteur]);

			for (int i=0; i<Data2[compteur]->GetSize(); i++)
			{
				tmp2 = Data2[compteur]->Elt(i) - moy2[compteur];
				var2[compteur] += tmp2*tmp2;
			}
//			var2[compteur] /= Data2[compteur]->GetSize();
		}
	}

	*nbData = Data1[0]->GetSize();

	for (compteur=0;compteur<lExpiry;compteur++)
	{
		*vol1 +=var1[compteur];
		*vol2 +=var2[compteur];
		*correl+=covar[compteur];
	}

	if (lExpiry1 > compteur)
	{
		for (compteur=lExpiry;compteur<lExpiry1;compteur++)
		{
			*vol1 +=var1[compteur];
		}
	}	
		
	if (lExpiry2 > compteur)
	{
		for (compteur=lExpiry;compteur<lExpiry2;compteur++)
		{
			*vol2 +=var2[compteur];
		}
	}	

	*correl /= sqrt((*vol1)* (*vol2));
*/
/*
	if ( (var1 == 0.) && (var2 == 0) )
		*correl = 1.;
	else if ( (var1 == 0.) || (var2 == 0) )
		*correl = 0.;
	else
	{
		//estimateur biaisé
		*correl /= Data1->GetSize()-1.;
		*correl /= (sqrt(1./(Data1->GetSize()-1.)*var1))*(sqrt(1./(Data2->GetSize()-1.)*var2));
	}
*/

//	return 0.;
//}

long calcCorrelMoy(ARM_Vector** Data1, ARM_Vector** Data2, double* correl, double* vol1, double* vol2, long* nbData, long lExpiry1, long lExpiry2)
{
	long lExpiry = MIN(lExpiry1,lExpiry2);
	int compteur;

	for (compteur=0;compteur<lExpiry;compteur++)
	{
		if (Data1[compteur]->GetSize() != Data2[compteur]->GetSize())
			return -1;
	}

	double* sumXiYi = new double[lExpiry];
	double* sumXi = new double[lExpiry1];
	double* sumYi = new double[lExpiry2];
	double* sumXi2 = new double[lExpiry1];
	double* sumYi2 = new double[lExpiry2];

	for (compteur=0;compteur<lExpiry;compteur++)
		sumXiYi[compteur]=0.;

	for (compteur=0;compteur<lExpiry1;compteur++)
	{
		sumXi[compteur]=0.;
		sumXi2[compteur]=0.;
	}

	for (compteur=0;compteur<lExpiry2;compteur++)
	{
		sumYi[compteur]=0.;
		sumYi2[compteur]=0.;
	}

	*correl = 0.;
	*vol1 = 0.;
	*vol2 = 0.;


	// Calcul des variances et covariances jusqu'au min des 2
	for (compteur = 0;compteur<MAX(lExpiry1,lExpiry2);compteur++)
	{
//		if (calcMoy2(Data1[compteur], Data2[compteur], &(moy1[compteur]), &(moy2[compteur])) == -1)
//			return -1;

		for (int i=0; i<Data1[compteur]->GetSize(); i++)
		{
			if (compteur < lExpiry1)
			{
				sumXi[compteur] += Data1[compteur]->Elt(i);
				sumXi2[compteur] += (Data1[compteur]->Elt(i)*Data1[compteur]->Elt(i));
			}
			if (compteur < lExpiry2)
			{
				sumYi[compteur] += Data2[compteur]->Elt(i);
				sumYi2[compteur] += (Data2[compteur]->Elt(i)*Data2[compteur]->Elt(i));
			}
			if (compteur < lExpiry)
			{
				sumXiYi[compteur] += (Data1[compteur]->Elt(i)*Data2[compteur]->Elt(i));
			}

		}
	}

	*nbData = Data1[0]->GetSize();

	for (compteur=0;compteur<lExpiry;compteur++)
	{
		*vol1 += sumXi2[compteur]/(Data1[compteur]->GetSize()-1.) - (Data1[compteur]->GetSize()-2.)*(sumXi[compteur]*sumXi[compteur])/((Data1[compteur]->GetSize()-1.)*(Data1[compteur]->GetSize()-1.)*(Data1[compteur]->GetSize()-1.));
		*vol2 += sumYi2[compteur]/(Data2[compteur]->GetSize()-1.) - (Data2[compteur]->GetSize()-2.)*(sumYi[compteur]*sumYi[compteur])/((Data2[compteur]->GetSize()-1.)*(Data2[compteur]->GetSize()-1.)*(Data2[compteur]->GetSize()-1.));
		*correl+= sumXiYi[compteur]/(Data1[compteur]->GetSize()-1.) - (Data1[compteur]->GetSize()-2.)*(sumXi[compteur]*sumYi[compteur])/((Data1[compteur]->GetSize()-1.)*(Data1[compteur]->GetSize()-1.)*(Data1[compteur]->GetSize()-1.));
	}

	if (lExpiry1 > compteur)
	{
		for (compteur=lExpiry;compteur<lExpiry1;compteur++)
		{
			*vol1 += sumXi2[compteur]/(Data1[compteur]->GetSize()-1.) - (Data1[compteur]->GetSize()-2.)*(sumXi[compteur]*sumXi[compteur])/((Data1[compteur]->GetSize()-1.)*(Data1[compteur]->GetSize()-1.)*(Data1[compteur]->GetSize()-1.));
		}
	}	
		
	if (lExpiry2 > compteur)
	{
		for (compteur=lExpiry;compteur<lExpiry2;compteur++)
		{
			*vol2 += sumYi2[compteur]/(Data2[compteur]->GetSize()-1.) - (Data2[compteur]->GetSize()-2.)*(sumYi[compteur]*sumYi[compteur])/((Data2[compteur]->GetSize()-1.)*(Data2[compteur]->GetSize()-1.)*(Data2[compteur]->GetSize()-1.));
		}
	}	

	*correl /= sqrt((*vol1)* (*vol2));

	return 0.;
}



long MakeGenSwap(ARM_Date AsOf, ARM_Date tmpDate1, long Plot, char* ccyName)
{
	ARM_Date swapStartDate, swapEndDate;

	ARM_Currency ccy(ccyName);
	double spotDays (ccy.GetSpotDays());

	swapStartDate = tmpDate1;
	swapStartDate.NextBusinessDay(spotDays,ccyName);

	//swapStartDate.AdjustToBusDate(ccyName,K_MOD_FOLLOWING);

	swapEndDate = swapStartDate;
	swapEndDate.AddMonths(Plot);

	double fixedRate (10.);
	double spread(0.);

	int resetFreq(-1);
	int payFreq(-1);

	ARM_INDEX_TYPE CCY_INDEX = GetDefaultIndexFromCurrency(ccyName);

	if (curSwap)
		delete curSwap;

	curSwap = new ARM_Swap(swapStartDate, swapEndDate, CCY_INDEX,
								spread, fixedRate, K_RCV, resetFreq, payFreq,
								&ccy);
	return 0;
}


double calcFwdYield(char* ccy, ARM_YCModel* yc, ARM_Date tmpDate1, long Plot)
{
	double fwdYield;

	ARM_Currency myCCY(ccy);

	ARM_Date startDate(tmpDate1);

	// compute the Fwd Rate
	if (Plot < 12)
	{
		startDate.NextBusinessDay(myCCY.GetSpotDays(),ccy);
//		startDate.AdjustToBusDate(ccy,K_MOD_FOLLOWING);

		ARM_Date tmpDate2 (startDate);
		tmpDate2.AddMonths(Plot);
		if (tmpDate2.IsBusinessDay(ccy) == 0)
			tmpDate2.NextBusinessDay(1,ccy);

		if ( (strcmp((const char*)ccy,"GBP") == 0)
			|| (strcmp((const char*)ccy,"JPY") == 0) )
		{
			fwdYield = yc->GetZeroCurve()->ForwardYield(startDate, tmpDate2,-1);
		}
		else
		{
			fwdYield = yc->GetZeroCurve()->ForwardYield(startDate, tmpDate2,-1)*360./365.;
		}

/*		char sSettleDate[12];
		(zc->GetAsOfDate()).JulianToStrDate(sSettleDate);
		char sStartDate[12];
		startDate.JulianToStrDate(sStartDate);
		char sTmpDate2[12];
		tmpDate2.JulianToStrDate(sTmpDate2);
		fprintf(Fp,"%s %s %s %s %lf\n",sSettleDate,ccy,sStartDate,sTmpDate2,fwdYield);
*/	}
	// Compute the Swap Forward Rate
	else
	{
		MakeGenSwap(yc->GetZeroCurve()->GetAsOfDate(), tmpDate1, Plot, ccy);

		double swapPrice(0.);
 
		curSwap->SetModel(yc);

		fwdYield = curSwap->PriceToRate(yc->GetZeroCurve()->GetAsOfDate(),swapPrice);

/*		char sSettleDate[12];
		(zc->GetAsOfDate()).JulianToStrDate(sSettleDate);
		char sStartDate[12];
		(curSwap->GetStartDate()).JulianToStrDate(sStartDate);
		char sEndDate[12];
		(curSwap->GetEndDate()).JulianToStrDate(sEndDate);
		fprintf(Fp,"%s %s %s %s %lf\n",sSettleDate,ccy,sStartDate,sEndDate,fwdYield);
*/	}

	return fwdYield;
}

long transformInMonth(const CCString& pStr)
{
	char* term = (char*)pStr;

    char matu;
	int Nb;

	sscanf(term, "%d%c", &Nb, &matu);
	matu = toupper(matu);

	if (term)
	{
		delete term;
		term = NULL;
	}

	if (matu == 'M')
		return Nb;

	if (matu == 'Y')
		return Nb*12;

	return 0;
}


long ARMLOCAL_GetCorrelInst(double date1,
							double date2,
							const CCString& ccy1,
							const CCString& index1,
							const CCString& expiry1,
							const CCString& tenor1,
							const CCString& ccy2,
							const CCString& index2,
							const CCString& expiry2,
							const CCString& tenor2,
							const CCString& curve1_ccy1,
							const CCString& curve2_ccy1,
							const CCString& Snbmonths_curve1_ccy1,
							const CCString& curve1_ccy2,
							const CCString& curve2_ccy2,
							const CCString& Snbmonths_curve1_ccy2,
							long typeId,
							double lambda,
							double precision,
							const CCString& ccy,
							ARM_result& result)
{
	long lExpiry1, lExpiry2, lTenor1, lTenor2;
	lExpiry1 = transformInMonth(expiry1);
	lExpiry2 = transformInMonth(expiry2);
	lTenor1 = transformInMonth(tenor1);
	lTenor2 = transformInMonth(tenor2);

	long nbmonths_curve1_ccy1=transformInMonth(Snbmonths_curve1_ccy1); 
	long nbmonths_curve1_ccy2=transformInMonth(Snbmonths_curve1_ccy2); 

	char sDate1[11];
	char sDate2[11];

	Local_XLDATE2ARMDATE(date1,sDate1);
	Local_XLDATE2ARMDATE(date2,sDate2);

	ARM_Date tmpdate1(sDate1);
	ARM_Date tmpdate2(sDate2);

	long nbjours = (long)(date2-date1+1);

	ARM_Vector* FwdData1 = new ARM_Vector(nbjours);
	ARM_Vector* FwdData2 = new ARM_Vector(nbjours);

	ARM_Date* DateData1 = new ARM_Date[nbjours];
	ARM_Date* DateData2 = new ARM_Date[nbjours];

	long nbDevises(1);

	if (strcmp((const char*)ccy1,(const char*)ccy2) )
		nbDevises = 2;

	ARM_Date curDate(tmpdate1);
	ARM_Date tmpDate;
	ARM_Date settleDate;
	long j1(0);
	long j2(0);

	double fwd;
	long zcId(-1);
	long zcSmoothId(-1);
	long modId(-1);

	ARM_ZeroLInterpol* zc = NULL;
	ARM_ZeroInterpolation* zcSmooth = NULL;
	ARM_YCModel* yc = NULL;

	long objId(-1); //a supprimer

	char* sCCY = (char*) ccy;

	Fp = fopen("correl.out","w");

	while (curDate <= tmpdate2)
	{
        char buf[30];
		curDate.JulianToStrDate(buf);
		double dDate = Local_ARMDATE2XLDATE(buf);
		
		long retCode;
		if (zcId != -1)
		{
			if (lExpiry1 <= nbmonths_curve1_ccy1)
			retCode = ARMLOCAL_GetZCFromSummit (index1,
												 ccy1,
												 curve1_ccy1,
												 dDate,
												 K_LINEAR,
												 result,
												 zcId);
			else
			retCode = ARMLOCAL_GetZCFromSummit (index1,
												 ccy1,
												 curve2_ccy1,
												 dDate,
												 K_LINEAR,
												 result,
												 zcId);
		}
		else
		{
			if (lExpiry1 <= nbmonths_curve1_ccy1)
			retCode = ARMLOCAL_GetZCFromSummit (index1,
												 ccy1,
												 curve1_ccy1,
												 dDate,
												 K_LINEAR,
												 result);
			else
			retCode = ARMLOCAL_GetZCFromSummit (index1,
												 ccy1,
												 curve2_ccy1,
												 dDate,
												 K_LINEAR,
												 result);

			if (retCode == ARM_OK)
				zcId = result.getLong();
			else
				zcId = -1;
		}


		if (zcId != -1)
		{
			if (zcSmoothId != -1)
				retCode = ARMLOCAL_ZCINTSMOOTH(zcId,
											   lambda,
											   precision,
											   result,
											   zcSmoothId);
			else
			{
				retCode = ARMLOCAL_ZCINTSMOOTH(zcId,
											   lambda,
											   precision,
											   result);

				if (retCode == ARM_OK)
					zcSmoothId = result.getLong();
				else
					zcSmoothId = -1;
			}

			if (zcSmoothId != -1)
			{
				if (modId != -1)
					retCode = ARMLOCAL_ycmod(zcSmoothId,
											 ARM_NULL_OBJECT,
											 result,
											 modId);
				else
				{
					retCode = ARMLOCAL_ycmod(zcSmoothId,
											 ARM_NULL_OBJECT,
											 result);

					if (retCode == ARM_OK)
						modId = result.getLong();
					else
						modId = -1;
				}
			}

			yc = (ARM_YCModel*)LOCAL_PERSISTENT_OBJECTS->GetObject(modId);

			if (yc != NULL)
			{
				settleDate = (yc->GetZeroCurve())->GetAsOfDate();
				DateData1[j1] = settleDate;

				tmpDate = settleDate;
				tmpDate.AddMonths(lExpiry1);
				if (tmpDate.IsBusinessDay(sCCY) == 0)
					tmpDate.NextBusinessDay(1,sCCY);

				fwd = calcFwdYield(ccy1, yc, tmpDate, lTenor1);
				FwdData1->Elt(j1) = fwd;

				if (nbDevises == 1)
				{
					DateData2[j1] = settleDate;
					tmpDate = settleDate;

					tmpDate.AddMonths(lExpiry2);
					if (tmpDate.IsBusinessDay(ccy1) == 0)
						tmpDate.NextBusinessDay(1,ccy1);

					fwd = calcFwdYield(ccy2, yc, tmpDate, lTenor2);
					FwdData2->Elt(j1) = fwd;

					j2++;
				}
				j1++;
			}

			if (nbDevises == 2)
			{
				ARM_ZeroLInterpol* zc2 = NULL;
				ARM_ZeroInterpolation* zcSmooth2 = NULL;
				ARM_YCModel* yc2 = NULL;

				if (objId != -1)
				{
					if (lExpiry2 <= nbmonths_curve1_ccy2)
					retCode = ARMLOCAL_GetZCFromSummit (index2,
														 ccy2,
														 curve1_ccy2,
														 dDate,
														 K_LINEAR,
														 result,
														 objId);
					else
					retCode = ARMLOCAL_GetZCFromSummit (index2,
														 ccy2,
														 curve2_ccy2,
														 dDate,
														 K_LINEAR,
														 result,
														 objId);

				}
				else
				{
					if (lExpiry2 <= nbmonths_curve1_ccy2)
					retCode = ARMLOCAL_GetZCFromSummit (index2,
														 ccy2,
														 curve1_ccy2,
														 dDate,
														 K_LINEAR,
														 result);
					else
					retCode = ARMLOCAL_GetZCFromSummit (index2,
														 ccy2,
														 curve2_ccy2,
														 dDate,
														 K_LINEAR,
														 result);

					objId = result.getLong();
				}

				if (objId != -1)
				{
					zc2 = (ARM_ZeroLInterpol*)(LOCAL_PERSISTENT_OBJECTS->GetObject(objId));
					if (zc2 != NULL)
					{
						zcSmooth2 = new ARM_ZeroInterpolation(zc2, lambda, precision);
						if (zcSmooth2 != NULL)
						{
							yc2 = new ARM_YCModel(zcSmooth2);
						}
					}

					if (yc2 != NULL)
					{
						settleDate = zc2->GetAsOfDate();
						DateData2[j2] = settleDate;

						tmpDate = settleDate;
						tmpDate.AddMonths(lExpiry2);
						if (tmpDate.IsBusinessDay(ccy2) == 0)
							tmpDate.NextBusinessDay(1,ccy2);

						fwd = calcFwdYield(ccy2, yc2, tmpDate, lTenor2);
						FwdData2->Elt(j2) = fwd;

						j2++;
					}
				}

				if (zcSmooth2)
					delete zcSmooth2;
				zcSmooth2 = NULL;

				if (yc2)
					delete yc2;
				yc2 = NULL;
			}
		}
		curDate.NextBusinessDay(sCCY);
	}

	fclose(Fp);

	if (sCCY)
		delete sCCY;
	sCCY = NULL;

	double correl;
	double vol1, vol2;
	long nbData;

	long i = prepareData(DateData1,DateData2,FwdData1,FwdData2,j1,j2);

	i = transformData(DateData1,FwdData1,FwdData2,typeId);

	i = calcCorrelInst(FwdData1, FwdData2, &correl, &vol1, &vol2, &nbData);

	if (FwdData1)
	{
		delete FwdData1;
		FwdData1 = NULL;
	}

	if (FwdData2)
	{
		delete FwdData2;
		FwdData2 = NULL;
	}

	if (DateData1)
	{
		delete [] DateData1;
		DateData1 = NULL;
	}

	if (DateData2)
	{
		delete [] DateData2;
		DateData2 = NULL;
	}

	result.setArray(correl,0);
	result.setArray(vol1,1);
	result.setArray(vol2,2);
	result.setArray(nbData,3);

	return ARM_OK;
}


long ARMLOCAL_GetMoyCorrel(double date1,
						   double date2,
						   const CCString& ccy1,
						   const CCString& index1,
						   const CCString& expiry1,
						   const CCString& tenor1,
						   const CCString& ccy2,
						   const CCString& index2,
						   const CCString& expiry2,
						   const CCString& tenor2,
						   const CCString& curve1_ccy1,
						   const CCString& curve2_ccy1,	
						   const CCString& Snbmonths_curve1_ccy1,
						   const CCString& curve1_ccy2,
						   const CCString& curve2_ccy2,	
						   const CCString& Snbmonths_curve1_ccy2,
						   long typeId,
						   double lambda,
						   double precision,
						   const CCString& ccy,
						   ARM_result& result)
{
	long lExpiry1, lExpiry2, lTenor1, lTenor2;
	lExpiry1 = transformInMonth(expiry1);
	lExpiry2 = transformInMonth(expiry2);
	lTenor1 = transformInMonth(tenor1);
	lTenor2 = transformInMonth(tenor2);

	long nbmonths_curve1_ccy1=transformInMonth(Snbmonths_curve1_ccy1); 
	long nbmonths_curve1_ccy2=transformInMonth(Snbmonths_curve1_ccy2); 

	long lExpiry = MIN(lExpiry1,lExpiry2);

	char sDate1[11];
	char sDate2[11];

	Local_XLDATE2ARMDATE(date1,sDate1);
	Local_XLDATE2ARMDATE(date2,sDate2);

	ARM_Date tmpdate1(sDate1);
	ARM_Date tmpdate2(sDate2);

	long nbjours = (long)(date2-date1+1);
	int compteur;

	ARM_Vector** FwdData1 = new ARM_Vector*[lExpiry1];
	ARM_Vector** FwdData2 = new ARM_Vector*[lExpiry2];

	for (compteur=0;compteur<lExpiry;compteur++)
	{
		FwdData1[compteur] = new ARM_Vector(nbjours);
		FwdData2[compteur] = new ARM_Vector(nbjours);
	}

	for (compteur=lExpiry;compteur<lExpiry1;compteur++)
	{
		FwdData1[compteur] = new ARM_Vector(nbjours);
	}

	for (compteur=lExpiry;compteur<lExpiry2;compteur++)
	{
		FwdData2[compteur] = new ARM_Vector(nbjours);
	}

//	ARM_Vector* FwdData1 = new ARM_Vector(nbjours);
//	ARM_Vector* FwdData2 = new ARM_Vector(nbjours);

	ARM_Date* DateData1 = new ARM_Date[nbjours];
	ARM_Date* DateData2 = new ARM_Date[nbjours];

	long nbDevises(1);

	if (strcmp((const char*)ccy1,(const char*)ccy2) )
		nbDevises = 2;

	ARM_Date curDate(tmpdate1);
	ARM_Date tmpDate;
	ARM_Date settleDate;
	long j1(0);
	long j2(0);

	double fwd;
	long objId(-1);

	long objId1(-1);
	long objId2(-1);

	ARM_ZeroLInterpol* zc = NULL;
	ARM_ZeroInterpolation* zcSmooth = NULL;
	ARM_YCModel* yc = NULL;

	ARM_ZeroLInterpol* zc_1 = NULL;
	ARM_ZeroInterpolation* zcSmooth_1 = NULL;
	ARM_YCModel* yc_1 = NULL;

	ARM_ZeroLInterpol* zc_2 = NULL;
	ARM_ZeroInterpolation* zcSmooth_2 = NULL;
	ARM_YCModel* yc_2 = NULL;

	char* sCCY = (char*) ccy;

	Fp = fopen("correl.out","w");

	while (curDate <= tmpdate2)
	{
        char buf[30];
		curDate.JulianToStrDate(buf);

		double dDate = Local_ARMDATE2XLDATE(buf);

		long retCode;
		if ((objId1 != -1) && (objId2 != -1))
		{
			retCode = ARMLOCAL_GetZCFromSummit (index1,
												 ccy1,
												 curve1_ccy1,
												 dDate,
												 K_LINEAR,
												 result,
												 objId1);
			retCode = ARMLOCAL_GetZCFromSummit (index1,
												 ccy1,
												 curve2_ccy1,
												 dDate,
												 K_LINEAR,
												 result,
												 objId2);
		}
		else
		{
			retCode = ARMLOCAL_GetZCFromSummit (index1,
												 ccy1,
												 curve1_ccy1,
												 dDate,
												 K_LINEAR,
												 result);

			if (retCode == ARM_OK)
				objId1 = result.getLong();
			else
				objId1 = -1;

			retCode = ARMLOCAL_GetZCFromSummit (index1,
												 ccy1,
												 curve2_ccy1,
												 dDate,
												 K_LINEAR,
												 result);

			if (retCode == ARM_OK)
				objId2 = result.getLong();
			else
				objId2 = -1;
		}

        int i;
		if ((objId1 != -1) && (objId2 != -1))
		{
			zc_1 = (ARM_ZeroLInterpol*)(LOCAL_PERSISTENT_OBJECTS->GetObject(objId1));
			zc_2 = (ARM_ZeroLInterpol*)(LOCAL_PERSISTENT_OBJECTS->GetObject(objId2));

			if ((zc_1 != NULL) && (zc_2 != NULL))
			{
				zcSmooth_1 = new ARM_ZeroInterpolation(zc_1, lambda, precision);
				zcSmooth_2 = new ARM_ZeroInterpolation(zc_2, lambda, precision);

				if ((zcSmooth_1 != NULL) && (zcSmooth_2 != NULL))
				{
					yc_1 = new ARM_YCModel(zcSmooth_1);
					yc_2 = new ARM_YCModel(zcSmooth_2);
				}
			}

			if ((yc_1 != NULL) && (yc_2 != NULL))
			{
				settleDate = zc_1->GetAsOfDate();
				DateData1[j1] = settleDate;
				
				if (nbDevises == 1)
					DateData2[j1] = settleDate;

				for (i=0;i<MAX(lExpiry1,lExpiry2);i++)
				{
					tmpDate = settleDate;
					tmpDate.AddMonths(i+1);

					if (tmpDate.IsBusinessDay(sCCY) == 0)
						tmpDate.NextBusinessDay(1,sCCY);

					if (i < lExpiry1)
					{
						if ((i+1) <= nbmonths_curve1_ccy1)
						fwd = calcFwdYield(ccy1, yc_1, tmpDate, lTenor1);
						else
						fwd = calcFwdYield(ccy1, yc_2, tmpDate, lTenor1);

						FwdData1[i]->Elt(j1) = fwd;
					}

					if (nbDevises == 1)
					{		
						if (i < lExpiry2)
						{
							if ((i+1) <= nbmonths_curve1_ccy2)
								fwd = calcFwdYield(ccy2, yc_1, tmpDate, lTenor2);
							else
								fwd = calcFwdYield(ccy2, yc_2, tmpDate, lTenor2);

							// meme tmpdate
							FwdData2[i]->Elt(j1) = fwd;
						}
					}
				}
				i++;
			}
			j1++;

			if (zcSmooth_1)
				delete zcSmooth_1;
			zcSmooth_1 = NULL;

			if (zcSmooth_2)
				delete zcSmooth_2;
			zcSmooth_2 = NULL;

			if (yc_1)
				delete yc_1;
			yc_1 = NULL;

			if (yc_2)
				delete yc_2;
			yc_2 = NULL;

			if (nbDevises == 2)
			{

				ARM_ZeroLInterpol* zc2_1 = NULL;
				ARM_ZeroInterpolation* zcSmooth2_1 = NULL;
				ARM_YCModel* yc2_1 = NULL;

				ARM_ZeroLInterpol* zc2_2 = NULL;
				ARM_ZeroInterpolation* zcSmooth2_2 = NULL;
				ARM_YCModel* yc2_2 = NULL;

				if ((objId1 != -1) && (objId2 != -1))
				{
					retCode = ARMLOCAL_GetZCFromSummit (index2,
														 ccy2,
														 curve1_ccy2,
														 dDate,
														 K_LINEAR,
														 result,
														 objId1);

					retCode = ARMLOCAL_GetZCFromSummit (index2,
														 ccy2,
														 curve2_ccy2,
														 dDate,
														 K_LINEAR,
														 result,
														 objId2);

				}
				else
				{
					retCode = ARMLOCAL_GetZCFromSummit (index2,
														 ccy2,
														 curve1_ccy2,
														 dDate,
														 K_LINEAR,
														 result);
					objId1 = result.getLong();

					retCode = ARMLOCAL_GetZCFromSummit (index2,
														 ccy2,
														 curve2_ccy2,
														 dDate,
														 K_LINEAR,
														 result);
					objId2 = result.getLong();
				}

				if ((objId1 != -1) && (objId2 != -1))
				{
					zc2_1 = (ARM_ZeroLInterpol*)(LOCAL_PERSISTENT_OBJECTS->GetObject(objId1));
					zc2_2 = (ARM_ZeroLInterpol*)(LOCAL_PERSISTENT_OBJECTS->GetObject(objId2));

					if ((zc2_1 != NULL) && (zc2_2 != NULL))
					{
						zcSmooth2_1 = new ARM_ZeroInterpolation(zc2_1, lambda, precision);
						zcSmooth2_2 = new ARM_ZeroInterpolation(zc2_2, lambda, precision);

						if ((zcSmooth2_1 != NULL) && (zcSmooth2_2 != NULL))
						{
							yc2_1 = new ARM_YCModel(zcSmooth2_1);
							yc2_2 = new ARM_YCModel(zcSmooth2_2);
						}
					}

					if ((yc2_1 != NULL) && (yc2_2 != NULL))
					{
						settleDate = zc2_1->GetAsOfDate();
						DateData2[j2] = settleDate;
						for (int i=0;i<lExpiry2;i++)
						{
							tmpDate = settleDate;
							tmpDate.AddMonths(i+1);

							if (tmpDate.IsBusinessDay(sCCY) == 0)
								tmpDate.NextBusinessDay(1,sCCY);

							if ((i+1) <= nbmonths_curve1_ccy2)
								fwd = calcFwdYield(ccy2, yc2_1, tmpDate, lTenor2);
							else
								fwd = calcFwdYield(ccy2, yc2_2, tmpDate, lTenor2);

							FwdData2[i]->Elt(j2) = fwd;
						}
						i++;
					}
					j2++;
				}

				if (zcSmooth2_1)
					delete zcSmooth2_1;
				zcSmooth2_1 = NULL;

				if (yc2_1)
					delete yc2_1;
				yc2_1 = NULL;

				if (zcSmooth2_2)
					delete zcSmooth2_2;
				zcSmooth2_2 = NULL;

				if (yc2_2)
					delete yc2_2;
				yc2_2 = NULL;
			}
		}
		curDate.NextBusinessDay(sCCY);
	}

	fclose(Fp);

	if (sCCY)
		delete sCCY;
	sCCY = NULL;

	double correl;
	double vol1, vol2;
	long nbData;

	if (nbDevises==1)
		j2 = j1;

	long i = prepareDataMoy(DateData1,DateData2,FwdData1,FwdData2,j1,j2,lExpiry1,lExpiry2);

	i = transformDataMoy(DateData1,FwdData1,FwdData2,typeId,lExpiry1,lExpiry2);

	i = calcCorrelMoy(FwdData1, FwdData2, &correl, &vol1, &vol2, &nbData, lExpiry1, lExpiry2);

	for(compteur=0;compteur<lExpiry1;compteur++)
	{
		delete FwdData1[compteur];
	}

	for(compteur=0;compteur<lExpiry2;compteur++)
	{
		delete FwdData2[compteur];
	}

	if (FwdData1)
	{
		delete [] FwdData1;
		FwdData1 = NULL;
	}

	if (FwdData2)
	{
		delete [] FwdData2;
		FwdData2 = NULL;
	}

	if (DateData1)
	{
		delete [] DateData1;
		DateData1 = NULL;
	}

	if (DateData2)
	{
		delete [] DateData2;
		DateData2 = NULL;
	}

	result.setArray(correl,0);
	result.setArray(vol1,1);
	result.setArray(vol2,2);
	result.setArray(nbData,3);

	return ARM_OK;
}



long ARMLOCAL_GetCorrelQuanto(double date1,
							  double date2,
							  const CCString& ccy,
							  const CCString& index,
							  const CCString& expiry,
							  const CCString& tenor,
							  const CCString& cvname1,
							  const CCString& cvname2,
							  const CCString& switchinmonth,
							  const CCString& domccy,
							  const CCString& domindex,
							  const CCString& forccy,
							  const CCString& forindex,
							  long typeId,
							  double lambda,
							  double precision,
							  const CCString& calccy,
							  long fwdOrNot,
							  ARM_result& result)
{
	long lExpiry, lTenor;
	lExpiry = transformInMonth(expiry);
	lTenor = transformInMonth(tenor);

	long lSwitchInMonth=transformInMonth(switchinmonth); 

	char sDate1[11];
	char sDate2[11];

	Local_XLDATE2ARMDATE(date1,sDate1);
	Local_XLDATE2ARMDATE(date2,sDate2);

	ARM_Date tmpdate1(sDate1);
	ARM_Date tmpdate2(sDate2);

	long nbjours = (long)(date2-date1+1);

	if (nbjours < 2)
	{
		return ARM_KO;
	}

	vector<lCurve> lCurves;

	ARM_Vector* FwdData1 = new ARM_Vector(nbjours);
	ARM_Vector* FwdData2 = new ARM_Vector(nbjours);

	ARM_Date* DateData1 = new ARM_Date[nbjours];
	ARM_Date* DateData2 = new ARM_Date[nbjours];

	ARM_Date curDate(tmpdate1);
	ARM_Date tmpDate;
	ARM_Date settleDate;
	long j1(0);
	long j2(0);

	double fwd;

	long zcId(-1);
	long zcSmoothId(-1);
	long modId(-1);
	ARM_ZeroInterpolation* zcsmooth = NULL;
	ARM_YCModel* yc = NULL;

	long zcDomId(-1);
	long zcSmoothDomId(-1);

	long zcForId(-1);
	long zcSmoothForId(-1);

	char* sCCY = (char*) ccy;
	double FX;

	FILE* FpQuanto = fopen("C:\\Program Files\\ARM\\correl.out","a+");

	long domIndexCurve, forIndexCurve;

	// Courbes pour la données de taux
	lCurve lYield;
	lYield.Currency = ccy;
	if (lExpiry <= lSwitchInMonth)
		lYield.CvName = cvname1;
	else
		lYield.CvName = cvname2;
	lYield.Index = index;
	lYield.id = -1;

	lCurves.push_back(lYield);

	// Courbes pour le change
	long trouveDom = 0;
	for (int i=0; i<lCurves.size(); i++)
	{
		lCurve curveI = lCurves[i];
		if (
			(strcmp(domccy,(const char*)curveI.Currency) == 0) &&
			(strcmp(domindex,(const char*)curveI.Index) == 0) &&
			(strcmp("MO",(const char*)curveI.CvName) == 0)
			)
		{
			domIndexCurve = i;
			trouveDom = 1;
		}
	}

	if (trouveDom == 0)
	{
		lCurve lDom;
		lDom.Currency = domccy;
		lDom.CvName = "MO";
		lDom.Index = domindex;
		lDom.id = -1;
		domIndexCurve = lCurves.size();
		lCurves.push_back(lDom);
	}

	long trouveFor = 0;
	for (int i=0; i<lCurves.size(); i++)
	{
		lCurve curveI = lCurves[i];
		if (
			(strcmp(forccy,(const char*)curveI.Currency) == 0) &&
			(strcmp(forindex,(const char*)curveI.Index) == 0) &&
			(strcmp("MO",(const char*)curveI.CvName) == 0)
			)
		{
			forIndexCurve = i;
			trouveFor = 1;
		}
	}

	if (trouveFor == 0)
	{
		lCurve lFor;
		lFor.Currency = forccy;
		lFor.CvName = "MO";
		lFor.Index = forindex;
		lFor.id = -1;
		forIndexCurve = lCurves.size();
		lCurves.push_back(lFor);
	}

	while (curDate <= tmpdate2)
	{
		char sSettleDate[30];
		curDate.JulianToStrDate(sSettleDate);
		double dDate = Local_ARMDATE2XLDATE(sSettleDate);

		long retCode = ARM_OK;

		for (int i=0;i<lCurves.size();i++)
		{
			retCode += ARMLOCAL_GetZCFromSummit (lCurves[i].Index,
												 lCurves[i].Currency,
												 lCurves[i].CvName,
												 dDate,
												 K_LINEAR,
												 result);
			retCode += ARMLOCAL_ZCINTSMOOTH(result.getLong(),
											lambda,
											precision,
											result,
											lCurves[i].id);

			if ((retCode == ARM_OK) && (lCurves[i].id == -1))
				lCurves[i].id = result.getLong();
		}
		
		// Récupération du forward de change
		retCode += ARMLOCAL_FxConvert(domccy,forccy,dDate,1.,"MO",result);

		if (retCode == ARM_OK)
			FX = result.getDouble();
		
		if (retCode == ARM_OK)
		{
			zcsmooth = (ARM_ZeroInterpolation *)LOCAL_PERSISTENT_OBJECTS->GetObject(lCurves[0].id);
			yc = new ARM_YCModel(zcsmooth);

			if (yc != NULL)
			{
				settleDate = (yc->GetZeroCurve())->GetAsOfDate();
				DateData1[j1] = settleDate;

				tmpDate = settleDate;
				tmpDate.AddMonths(lExpiry);
				if (tmpDate.IsBusinessDay(sCCY) == 0)
					tmpDate.NextBusinessDay(1,sCCY);

				fwd = calcFwdYield(ccy, yc, tmpDate, lTenor);
				FwdData1->Elt(j1) = fwd;

				j1++;
			}
			// Fin de la récupération du taux

			delete yc;

			// Calcul de la date de forward sur FX : J+2 de la date de fixing sur le calendrier croisé des 2 devises
			retCode = ARMLOCAL_NextBusinessDay(dDate, domccy+forccy, 2, result);

			double J2 = Local_ARMDATE2XLDATE(result.getString());

			double J2Frac = (J2-dDate)/365.;

			retCode = ARMLOCAL_DiscountPrice(lCurves[domIndexCurve].id,J2Frac,result);
			double tmpVal1 = result.getDouble();

			retCode = ARMLOCAL_DiscountPrice(lCurves[forIndexCurve].id,J2Frac,result);
			double tmpVal2 = result.getDouble();

			// rustine pour contrôler que le smooth ne déconne pas sur la partie courte de la courbe
			if ( (tmpVal2/tmpVal1 < 0.5) || (tmpVal2/tmpVal1 > 2.) )
			{
				FX = -1.;
				j1--;
			}
			else
				FX *= tmpVal2/tmpVal1;

			if (FX != -1)
			{
				if (fwdOrNot == 1)
				{
					char* sTmpDate = new char[11];
					tmpDate.JulianToStrDate(sTmpDate);

					double fixingDate = Local_ARMDATE2XLDATE(sTmpDate);

					// Calcul de la date de forward sur FX : J+2 de la date de fixing sur le calendrier croisé des 2 devises
					retCode = ARMLOCAL_NextBusinessDay(fixingDate, domccy+forccy, 2, result);

					if (sTmpDate)
						delete [] sTmpDate;

					sTmpDate = result.getString();

					fixingDate = Local_ARMDATE2XLDATE(result.getString());

					double yearFrac = (fixingDate-dDate)/365.;

					retCode = ARMLOCAL_DiscountPrice(lCurves[domIndexCurve].id,yearFrac,result);
					double val1 = result.getDouble();

					retCode = ARMLOCAL_DiscountPrice(lCurves[forIndexCurve].id,yearFrac,result);
					double val2 = result.getDouble();

			//		Fp = fopen("C:\\Program Files\\ARM\\correl.out","a+");
					fprintf(FpQuanto,"%s %lf %s %s %s %lf %lf %lf %lf %lf %lf %lf\n",sSettleDate,fwd,(const char*)domccy,(const char*)forccy, sTmpDate, FX *tmpVal1/tmpVal2,tmpVal1,tmpVal2,FX, val1, val2, FX * val1 / val2);
			//		fclose(Fp);

					if (sTmpDate)
						delete [] sTmpDate;
					sTmpDate = NULL;

					FX *= val1/val2;
				}

				DateData2[j2] = settleDate;

				FwdData2->Elt(j2) = FX;

				j2++;
				// Fin de la récupération du forward de change
			}
		}
		curDate.NextBusinessDay(sCCY);
	}

	fclose(FpQuanto);

	if (sCCY)
		delete sCCY;
	sCCY = NULL;

	double correl;
	double vol1, vol2;
	long nbData;
    int i;
	i = prepareData(DateData1,DateData2,FwdData1,FwdData2,j1,j2);

	i = transformData(DateData1,FwdData1,FwdData2,typeId);

	i = calcCorrelInst(FwdData1, FwdData2, &correl, &vol1, &vol2, &nbData);

	if (FwdData1)
	{
		delete FwdData1;
		FwdData1 = NULL;
	}

	if (FwdData2)
	{
		delete FwdData2;
		FwdData2 = NULL;
	}

	if (DateData1)
	{
		delete [] DateData1;
		DateData1 = NULL;
	}

	if (DateData2)
	{
		delete [] DateData2;
		DateData2 = NULL;
	}

	result.setArray(correl,0);
	result.setArray(vol1,1);
	result.setArray(vol2,2);
	result.setArray(nbData,3);

	return ARM_OK;
}


long ARMLOCAL_GetMoyCorrelQuanto(double date1,
								 double date2,
								 const CCString& ccy,
								 const CCString& index,
								 const CCString& expiry,
								 const CCString& tenor,
								 const CCString& cvname1,
								 const CCString& cvname2,
								 const CCString& switchinmonth,
								 const CCString& domccy,
								 const CCString& domindex,
								 const CCString& forccy,
								 const CCString& forindex,
								 long typeId,
								 double lambda,
								 double precision,
								 const CCString& calccy,
								 long fwdOrNot,
								 ARM_result& result)
{
	long lExpiry, lTenor;
	lExpiry = transformInMonth(expiry);
	lTenor = transformInMonth(tenor);

	long lSwitchInMonth=transformInMonth(switchinmonth); 

	char sDate1[11];
	char sDate2[11];

	Local_XLDATE2ARMDATE(date1,sDate1);
	Local_XLDATE2ARMDATE(date2,sDate2);

	ARM_Date tmpdate1(sDate1);
	ARM_Date tmpdate2(sDate2);

	long nbjours = (long)(date2-date1+1);

	if (nbjours < 2)
		return ARM_KO;

	ARM_Vector** FwdData1 = new ARM_Vector*[lExpiry];
	ARM_Vector** FwdData2 = new ARM_Vector*[lExpiry];

	for (int compteur=0;compteur<lExpiry;compteur++)
	{
		FwdData1[compteur] = new ARM_Vector(nbjours);
		FwdData2[compteur] = new ARM_Vector(nbjours);
	}

	ARM_Date* DateData1 = new ARM_Date[nbjours];
	ARM_Date* DateData2 = new ARM_Date[nbjours];

	ARM_Date curDate(tmpdate1);
	ARM_Date tmpDate;
	ARM_Date settleDate;
	long j1(0);
	long j2(0);

	double fwd;

	long zcId(-1);
	long zcSmoothId(-1);
	long modId(-1);
	ARM_ZeroInterpolation* zcsmooth1 = NULL;
	ARM_ZeroInterpolation* zcsmooth2 = NULL;
	ARM_YCModel* yc1 = NULL;
	ARM_YCModel* yc2 = NULL;

	long zcDomId(-1);
	long zcSmoothDomId(-1);

	long zcForId(-1);
	long zcSmoothForId(-1);

	char* sCCY = (char*) calccy;

	double FX = 0.0;

	FILE* FpQuanto = fopen("C:\\Program Files\\ARM\\correl.out","a+");

	vector<lCurve> lCurves;
	long domIndexCurve, forIndexCurve;
	long yIndexCurve1,yIndexCurve2;

	// Courbes pour la données de taux
	lCurve lYield1;
	lYield1.Currency = ccy;
	lYield1.CvName = cvname1;
	lYield1.Index = index;
	lYield1.id = -1;
	lCurves.push_back(lYield1);
	yIndexCurve1 = 0;

	if(strcmp((const char*) cvname1,(const char*) cvname2) != 0)
	{
		lCurve lYield2;
		lYield2.Currency = ccy;
		lYield2.CvName = cvname2;
		lYield2.Index = index;
		lYield2.id = -1;
		lCurves.push_back(lYield2);
		yIndexCurve2 = 1;
	}
	else
		yIndexCurve2 = yIndexCurve1;

	// Courbes pour le change
	long trouveDom = 0;
	for (int i=0; i<lCurves.size(); i++)
	{
		lCurve curveI = lCurves[i];
		if (
			(strcmp(domccy,(const char*)curveI.Currency) == 0) &&
			(strcmp(domindex,(const char*)curveI.Index) == 0) &&
			(strcmp("MO",(const char*)curveI.CvName) == 0)
			)
		{
			domIndexCurve = i;
			trouveDom = 1;
		}
	}

	if (trouveDom == 0)
	{
		lCurve lDom;
		lDom.Currency = domccy;
		lDom.CvName = "MO";
		lDom.Index = domindex;
		lDom.id = -1;
		domIndexCurve = lCurves.size();
		lCurves.push_back(lDom);
	}

	long trouveFor = 0;
	for (int i=0; i<lCurves.size(); i++)
	{
		lCurve curveI = lCurves[i];
		if (
			(strcmp(forccy,(const char*)curveI.Currency) == 0) &&
			(strcmp(forindex,(const char*)curveI.Index) == 0) &&
			(strcmp("MO",(const char*)curveI.CvName) == 0)
			)
		{
			forIndexCurve = i;
			trouveFor = 1;
		}
	}

	if (trouveFor == 0)
	{
		lCurve lFor;
		lFor.Currency = forccy;
		lFor.CvName = "MO";
		lFor.Index = forindex;
		lFor.id = -1;
		forIndexCurve = lCurves.size();
		lCurves.push_back(lFor);
	}

	while (curDate <= tmpdate2)
	{
		char sSettleDate[30];
		curDate.JulianToStrDate(sSettleDate);

		double dDate = Local_ARMDATE2XLDATE(sSettleDate);

		long retCode=ARM_OK;

		// Récupération des courbes de taux
		for (int i=0;i<lCurves.size();i++)
		{
			retCode += ARMLOCAL_GetZCFromSummit (lCurves[i].Index,
												 lCurves[i].Currency,
												 lCurves[i].CvName,
												 dDate,
												 K_LINEAR,
												 result);
			retCode += ARMLOCAL_ZCINTSMOOTH(result.getLong(),
											lambda,
											precision,
											result,
											lCurves[i].id);

			if ((retCode == ARM_OK) && (lCurves[i].id == -1))
				lCurves[i].id = result.getLong();
		}


		// Récupération du forward de change
		retCode += ARMLOCAL_FxConvert(domccy,forccy,dDate,1.,"MO",result);


		if (retCode == ARM_OK)
		{
			FX = result.getDouble();

			// Calcul de la date de forward sur FX : J+2 de la date de fixing sur le calendrier croisé des 2 devises
			retCode = ARMLOCAL_NextBusinessDay(dDate, domccy+forccy, 2, result);

			double J2 = Local_ARMDATE2XLDATE(result.getString());

			double J2Frac = (J2-dDate)/365.;

			retCode = ARMLOCAL_DiscountPrice(lCurves[domIndexCurve].id,J2Frac,result);
			double tmpVal1 = result.getDouble();

			retCode = ARMLOCAL_DiscountPrice(lCurves[forIndexCurve].id,J2Frac,result);
			double tmpVal2 = result.getDouble();

			// rustine pour contrôler que le smooth ne déconne pas sur la partie courte de la courbe
			if ( (tmpVal2/tmpVal1 < 0.5) || (tmpVal2/tmpVal1 > 2.) )
				FX = -1.;
			else
				FX *= tmpVal2/tmpVal1;

			if (FX != -1.)
			{
				zcsmooth1 = (ARM_ZeroInterpolation *)LOCAL_PERSISTENT_OBJECTS->GetObject(lCurves[yIndexCurve1].id);
				yc1 = new ARM_YCModel(zcsmooth1);

				if (yIndexCurve1 != yIndexCurve2)
				{
					zcsmooth2 = (ARM_ZeroInterpolation *)LOCAL_PERSISTENT_OBJECTS->GetObject(lCurves[yIndexCurve2].id);
					yc2 = new ARM_YCModel(zcsmooth2);
				}
				else
				{
					yc2 = yc1;
				}
			}

			if ((yc1 != NULL) && (yc2 != NULL))
			{
				settleDate = (yc1->GetZeroCurve())->GetAsOfDate();
				DateData1[j1] = settleDate;

				for (int i=0;i<lExpiry;i++)
				{
					tmpDate = settleDate;
					tmpDate.AddMonths(i+1);

					if (tmpDate.IsBusinessDay(sCCY) == 0)
						tmpDate.NextBusinessDay(1,sCCY);

					if ((i+1) <= lSwitchInMonth)
						fwd = calcFwdYield(ccy, yc1, tmpDate, lTenor);
					else
						fwd = calcFwdYield(ccy, yc2, tmpDate, lTenor);

					FwdData1[i]->Elt(j1) = fwd;

					if (fwdOrNot == 1)
					{
						char* sTmpDate = new char[11];
						tmpDate.JulianToStrDate(sTmpDate);

						double fixingDate = Local_ARMDATE2XLDATE(sTmpDate);

						// Calcul de la date de forward sur FX : J+2 de la date de fixing sur le calendrier croisé des 2 devises
						retCode = ARMLOCAL_NextBusinessDay(fixingDate, domccy+forccy, 2, result);

						if (sTmpDate)
							delete [] sTmpDate;

						sTmpDate = result.getString();

						fixingDate = Local_ARMDATE2XLDATE(result.getString());

						double yearFrac = (fixingDate-dDate)/365.;

						retCode = ARMLOCAL_DiscountPrice(lCurves[domIndexCurve].id,yearFrac,result);
						double val1 = result.getDouble();

						retCode = ARMLOCAL_DiscountPrice(lCurves[forIndexCurve].id,yearFrac,result);
						double val2 = result.getDouble();

						fprintf(FpQuanto,"%s %lf %s %s %s %lf %lf %lf %lf %lf %lf %lf\n",sSettleDate,fwd,(const char*)domccy,(const char*)forccy, sTmpDate, FX*tmpVal1/tmpVal2,tmpVal1,tmpVal2,FX, val1, val2, FX * val1 / val2);
						
						if (sTmpDate)
							delete [] sTmpDate;
						sTmpDate = NULL;

						fwd = FX * val1/val2;
					}
					else
					{
						fwd = FX;
					}

					DateData2[j1] = settleDate;

					FwdData2[i]->Elt(j1) = fwd;
				}

				if (yc1)
					delete yc1;
				yc1 = NULL;

				if (yIndexCurve1 != yIndexCurve2)
				{
					if (yc2)
						delete yc2;
					yc2 = NULL;
				}

				j1++;
			}
		}
		curDate.NextBusinessDay(sCCY);
	}

	fclose(FpQuanto);

	if (sCCY)
		delete sCCY;
	sCCY = NULL;

	double correl;
	double vol1, vol2;
	long nbData;
    int i;
	i = prepareDataMoy(DateData1,DateData2,FwdData1,FwdData2,j1,j1,lExpiry,lExpiry);

	i = transformDataMoy(DateData1,FwdData1,FwdData2,typeId,lExpiry,lExpiry);

	i = calcCorrelMoy(FwdData1, FwdData2, &correl, &vol1, &vol2, &nbData, lExpiry, lExpiry);

	for(int compteur=0;compteur<lExpiry;compteur++)
	{
		delete FwdData1[compteur];
		delete FwdData2[compteur];
	}

	if (FwdData1)
	{
		delete [] FwdData1;
		FwdData1 = NULL;
	}

	if (FwdData2)
	{
		delete [] FwdData2;
		FwdData2 = NULL;
	}

	if (DateData1)
	{
		delete [] DateData1;
		DateData1 = NULL;
	}

	if (DateData2)
	{
		delete [] DateData2;
		DateData2 = NULL;
	}

	result.setArray(correl,0);
	result.setArray(vol1,1);
	result.setArray(vol2,2);
	result.setArray(nbData,3);

	return ARM_OK;
}


long ARMLOCAL_GetDealsFromSummitFilter (const CCString& filter,
										VECTOR<CCString>& listDeals,
										VECTOR<CCString>& listTypes,
										ARM_result& result)
{

	if (GetDataRetrieverVersion() == FFRETRIEVER)
	{
		result.setMsg ("ARM_ERR: Function not implemented without ETK");
		return ARM_KO;
	}

	CCString msg (" ");

	try
	{
		listDeals = ARMLOCAL_GetDealByFilter(filter,listTypes);
	}
	catch (Exception& x)
	{
		x.DebugPrint();

		ARM_RESULT();
	}

	if (listDeals.size() == 0)
		return ARM_KO;

	return ARM_OK;
}


double ARM_CptFwd(ARM_Date asOf,
				  const CCString& ccy,
				  const CCString& index,
				  const CCString& cvName,
				  double expiry,
				  double matu)
{
	double res;
	ARM_result result;

	ARM_ZeroLInterpol* createdZcLin = NULL;
	
	if (GetDataRetrieverVersion () >= ETKRETRIEVER)
	{
		CCString xmlResponse;

		xmlResponse = etoolkit_getXMLMYAndZCFromSummit(index,ccy,cvName,asOf);

		int interp = ARMLOCAL_ParseXMLForInterpolator(xmlResponse);

		xmlResponse = etoolkit_getXMLZCFromSummit(index,ccy,cvName,asOf);

		createdZcLin = ARMLOCAL_ParseXMLForZC(xmlResponse, asOf,(const char*) ccy, interp);
	}
	else
	{
		createdZcLin = GetZCFromSummitNoETK(index,
											ccy,
											cvName,
											asOf,
											K_LINEAR,
											result);
	}

	if (createdZcLin == NULL)
	{
		throw Exception(__LINE__, __FILE__, ERR_INVALID_INPUT,
						"check your parameters for getting ZcCurve");
	}

	double fwdYield;

	if (matu < 1.0)
	{
		// Compute fwd Rate
		fwdYield = createdZcLin->ForwardYield(expiry, expiry + matu, -1);

		if ( (strcmp((const char*)ccy,"GBP") != 0)
			&& (strcmp((const char*)ccy,"JPY") != 0) )
		{
			fwdYield *= 360.0/365.0;
		}
		res = fwdYield;
	}
	else
	{
		// Compute swap fwd Rate
		char currency[4];
		strcpy(currency,(const char*) ccy);
		ARM_Currency isoccy(currency);
		ARM_Swap* curSwap = MakeGenSwap(asOf, expiry, matu, &isoccy);

		ARM_YCModel mod(createdZcLin);

		curSwap->SetModel(&mod);

		res = curSwap->PriceToRate(asOf,0.0);

		delete curSwap;
	}

	delete createdZcLin;

	return res;
}


long ARMLOCAL_GetAsOfVolOrRate(double date1,
							   const CCString& ccy,
							   const CCString& index,
							   const CCString& cvName,
							   const CCString& expiry,
							   const CCString& matu,
							   long yieldOrValId,
							   long calcModId,
							   const CCString& volType,
							   ARM_result& result)
{

	double curMatu;
	double curExpiry;
	double res;

	curMatu = FromStrMatuToDouble((const char*)matu);
	curExpiry = FromStrMatuToDouble((const char*)expiry);

	if (curMatu <= 0.0)
	{
		result.setMsg ("ARM_ERR: Maturity must be > 0");
		return ARM_KO;
	}

	if (curExpiry < 0.0)
	{
		result.setMsg ("ARM_ERR: Expiry must be >= 0");
		return ARM_KO;
	}

	char sDate[11];
	Local_XLDATE2ARMDATE(date1,sDate);
	ARM_Date asof(sDate);

	CCString msg (" ");

	try
	{
		if (yieldOrValId == ARM_YIELD)
		{
			res = ARM_CptFwd(asof, ccy, index, cvName, curExpiry,curMatu);
		}
		else
		{
			double vol;
			double fwdyield;

			fwdyield = ARM_CptFwd(asof, ccy, index, cvName, curExpiry,curMatu);

			ARM_VolLInterpol* VolATM = NULL;

			if (GetDataRetrieverVersion() >= ETKRETRIEVER)
			{
				VolATM = etoolkit_GetVolATMFromSummit(index,
													  ccy,
													  cvName,
													  asof,
													  volType);
			}
			else
			{
				VolATM = ARMLOCAL_GetVolATMFromSummit(index,
													  ccy,
													  cvName,
													  asof,
													  volType,
													  result);

				if ( (VolATM == NULL) && (GetFallBackDataRetrieverVersion() != 0))
				{
					VolATM = etoolkit_GetVolATMFromSummit(index,
														  ccy,
														  cvName,
														  asof,
														  volType);
				}
			}

			if (VolATM == NULL)
			{
				throw Exception(__LINE__, __FILE__, ERR_INVALID_INPUT,
								"check your parameters for getting VolCurve");
			}

			vol = VolATM->ComputeVolatility(curExpiry,curMatu);

			res = vol;

			switch(calcModId)
			{
			case ARM_NOR_MOD :
				{
					res = vol * fwdyield;
				};
				break;

			case ARM_MOD_SQR :
				{
					res = sqrt(fwdyield) * vol;
				}
				break;
			}

			delete VolATM;

		}

		result.setDouble(res);
	}
	catch (Exception& x)
	{
		x.DebugPrint();

		ARM_RESULT();
	}

	return ARM_OK;

}


long ARMLOCAL_GetFwdRatesMatrix(double asOf,
								const CCString& ccy,
								const CCString& index,
								const CCString& cvName,
								ARM_result& result)
{

	ARM_ZeroLInterpol* zc = obj_getZcFromSummit(index,ccy,cvName,asOf,K_LINEAR,result);

	for (int i = 1; i <= MATU_COLS_SIZE; i++)
		result.setArray(ARM_MATU_COLS[i-1],i);

	for (int i = 1; i <= MATU_LINES_SIZE; i++)
		result.setArray(ARM_MATU_LINES[i-1],i*(MATU_COLS_SIZE+1));

	for (int i = 0; i < MATU_LINES_SIZE; i++)
	{
		for (int j = 0; j < MATU_COLS_SIZE; j++)
		{
			double fwd = ApproximatedForwardRate(ARM_MATU_LINES[i],ARM_MATU_COLS[j],zc,1);

			result.setArray(fwd,(i+1)*(MATU_COLS_SIZE+1)+j+1);
		}
	}

	if (zc)
		delete zc;
	zc = NULL;

	return ARM_OK;

}


long ARMLOCAL_ARM_GetInfo (long secId,
						   const CCString& secClass,
						   const CCString& type,
						   ARM_result& result,
						   long objId)
{
	CCString msg (" ");

	try
	{
		if (secClass == LOCAL_SPREAD_OPTION_CLASS)
		{
			ARM_SpreadOption *spro = (ARM_SpreadOption *) LOCAL_PERSISTENT_OBJECTS->GetObject(secId);

			if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(spro, ARM_SPREADOPTION) == 0)
			{
				result.setMsg ("ARM_ERR: object is not a spreadoption");
				return ARM_KO;
			}

			if (type == "WEIGHT1")
				result.setDouble(spro->GetWeight1());
			else if (type == "WEIGHT2")
				result.setDouble(spro->GetWeight2());
			else if (type == "SLOPEFLAG")
				result.setDouble(spro->GetSlopeFlag());
			else if (type == "TENOR1")
				result.setString(YearTermToStringMatu(spro->GetUnderMatSpread(1)).c_str());
			else if (type == "TENOR2")
				result.setString(YearTermToStringMatu(spro->GetUnderMatSpread(2)).c_str());
			else if (type == "PAYTENOR")
				result.setString(YearTermToStringMatu(spro->GetUnderMatPayIndex()).c_str());
			else if (type == "STRIKE")
			{
				ARM_ReferenceValue* strike = NULL;
				long strikeId;

				strike = (ARM_ReferenceValue*) spro->GetStrikes()->Clone();
				if ( objId == -1 )
				{
					CREATE_GLOBAL_OBJECT();

					strikeId = LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)strike);

					if ( strikeId == RET_KO )
					{
						if (strike)
							delete strike;
						strike = NULL;

						result.setMsg ("ARM_ERR: Pb with inserting object");				
						return ARM_KO;
					}

					result.setLong(strikeId);

					return ARM_OK;
				}
				else
				{
					ARM_ReferenceValue* oldStrike = (ARM_ReferenceValue *) LOCAL_PERSISTENT_OBJECTS->GetObject(objId);

					if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(oldStrike, ARM_REFERENCE_VALUE) == 1)
					{
						if (oldStrike)
						{
							delete oldStrike;
							oldStrike = NULL;
						}

						LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)strike, objId);

						result.setLong(objId);
						return ARM_OK;
					}
					else
					{
						if (strike)
							delete strike;
						strike = NULL;

						result.setMsg ("ARM_ERR: previous object is not of a good type");
						return ARM_KO;
					}
				}
			}
			else if ( (type == "LEG1") || (type == "LEG2") || (type == "PAYLEG") )
			{
				ARM_SwapLeg* leg = NULL;
				long legId;

				if (type == "LEG1")
					leg = (ARM_SwapLeg*) spro->GetSpreadLeg()->GetFirstLeg()->Clone();
				else if (type == "LEG2")
					leg = (ARM_SwapLeg*) spro->GetSpreadLeg()->GetSecondLeg()->Clone();
				else if ( (type == "PAYLEG") && (spro->IsCorridorSpreadOption() || spro->IsDigitalFLT()) )
					leg = (ARM_SwapLeg*) spro->GetPayIndexLeg()->Clone();
				else
				{
					result.setMsg ("ARM_ERR: you can not ask the Payleg for this SO");
					return ARM_KO;
				}

				if ( objId == -1 )
				{
					CREATE_GLOBAL_OBJECT();

					legId = LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)leg);

					if ( legId == RET_KO )
					{
						if (leg)
							delete leg;
						leg = NULL;

						result.setMsg ("ARM_ERR: Pb with inserting object");				
						return ARM_KO;
					}

					result.setLong(legId);

					return ARM_OK;
				}
				else
				{
					ARM_SwapLeg* oldLeg = (ARM_SwapLeg *) LOCAL_PERSISTENT_OBJECTS->GetObject(objId);

					if (LocalPersistent::LOCAL_IS_SWAPLEG_OK(oldLeg) == 1)
					{
						if (oldLeg)
						{
							delete oldLeg;
							oldLeg = NULL;
						}

						LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)leg, objId);

						result.setLong(objId);
						return ARM_OK;
					}
					else
					{
						if (leg)
							delete leg;
						leg = NULL;

						result.setMsg ("ARM_ERR: previous object is not of a good type");
						return ARM_KO;
					}
				}
			}
			else
			{
				result.setMsg ("ARM_ERR: Bad Type, valid are WEIGHT1, WEIGHT2, SLOPEFLAG, STRIKE, LEG1, LEG2 or PAYLEG");
				return ARM_KO;
			}
		}
		else if (secClass == LOCAL_SWAP_CLASS)
		{
			ARM_Swap *swap = (ARM_Swap *) LOCAL_PERSISTENT_OBJECTS->GetObject(secId);

			if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(swap, ARM_SWAP) == 0)
			{
				result.setMsg ("ARM_ERR: object is not a swap");
				return ARM_KO;
			}
			if ( (type == "LEG1") || (type == "LEG2") )
			{
				ARM_SwapLeg* leg = NULL;
				long legId;

				if (type == "LEG1")
					leg = (ARM_SwapLeg*) swap->GetFirstLeg()->Clone();
				else
					leg = (ARM_SwapLeg*) swap->GetSecondLeg()->Clone();

				if ( objId == -1 )
				{
					CREATE_GLOBAL_OBJECT();

					legId = LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)leg);

					if ( legId == RET_KO )
					{
						if (leg)
							delete leg;
						leg = NULL;

						result.setMsg ("ARM_ERR: Pb with inserting object");				
						return ARM_KO;
					}

					result.setLong(legId);

					return ARM_OK;
				}
				else
				{
					ARM_SwapLeg* oldLeg = (ARM_SwapLeg *) LOCAL_PERSISTENT_OBJECTS->GetObject(objId);

					if (LocalPersistent::LOCAL_IS_SWAPLEG_OK(oldLeg) == 1)
					{
						if (oldLeg)
						{
							delete oldLeg;
							oldLeg = NULL;
						}

						LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)leg, objId);

						result.setLong(objId);
						return ARM_OK;
					}
					else
					{
						if (leg)
							delete leg;
						leg = NULL;

						result.setMsg ("ARM_ERR: previous object is not of a good type");
						return ARM_KO;
					}
				}
			}
			else
			{
				result.setMsg ("ARM_ERR: Bad Type, valid are LEG1 or LEG2");
				return ARM_KO;
			}
		}
		else if (secClass == LOCAL_OPTION_PORTFOLIO_CLASS)
		{
			ARM_ReferenceValue* strike = NULL;
			long strikeId;

			ARM_OptionPortfolio *opPf = (ARM_OptionPortfolio *) LOCAL_PERSISTENT_OBJECTS->GetObject(secId);

			if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(opPf, ARM_OPTIONPORTFOLIO) == 0)
			{
				result.setMsg ("ARM_ERR: object is not an Option Portfolio");
				return ARM_KO;
			}
			if ( (type == "BARRIERUP") || (type == "BARRIERDOWN"))
			{
				ARM_CorridorLeg* leg = (ARM_CorridorLeg*) (opPf->GetPtf()->GetAsset(0));

				if (type == "BARRIERUP")
					strike = (ARM_ReferenceValue*) leg->GetUpBarriers()->Clone();
				else
					strike = (ARM_ReferenceValue*) leg->GetDownBarriers()->Clone();

				if ( objId == -1 )
				{
					CREATE_GLOBAL_OBJECT();

					strikeId = LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)strike);

					if ( strikeId == RET_KO )
					{
						if (strike)
							delete strike;
						strike = NULL;

						result.setMsg ("ARM_ERR: Pb with inserting object");				
						return ARM_KO;
					}

					result.setLong(strikeId);

					return ARM_OK;
				}
				else
				{
					ARM_ReferenceValue* oldStrike = (ARM_ReferenceValue *) LOCAL_PERSISTENT_OBJECTS->GetObject(objId);

					if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(oldStrike, ARM_REFERENCE_VALUE) == 1)
					{
						if (oldStrike)
						{
							delete oldStrike;
							oldStrike = NULL;
						}

						LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)strike, objId);

						result.setLong(objId);
						return ARM_OK;
					}
					else
					{
						if (strike)
							delete strike;
						strike = NULL;

						result.setMsg ("ARM_ERR: previous object is not of a good type");
						return ARM_KO;
					}
				}
			}
			else if (type == "MRMIN")
			{
				result.setDouble(opPf->GetCraPricing()->Elt(2));
			}
			else if (type == "MRMAX")
			{
				result.setDouble(opPf->GetCraPricing()->Elt(3));
			}
			else
			{
				result.setMsg ("ARM_ERR: Bad Type, valid are BARRIERUP,BARRIERDOWN, MRMIN or MRMAX");
				return ARM_KO;
			}
		}
		else
		{
			result.setMsg ("ARM_ERR: Bad Instrument, valid is SpreadOption,Swap or OPTIONPORTFOLIO");
			return ARM_KO;
		}
	}
	catch (Exception& x)
	{
		x.DebugPrint();

		ARM_RESULT();
	}

	return ARM_OK;
}

long ARMLOCAL_ARM_GetModelFactor(double AsOfDate,
								 const CCString& model,
								 const CCString& type,
								 const CCString& factorName,
								 const CCString& ccy,
								 const CCString& index,
								 const CCString& cvName,
								 long calcModId,
								 ARM_result& result,
								 long objId)
{
	if (GetDataRetrieverVersion () == FFRETRIEVER)
	{
		result.setMsg ("ARM_ERR: This function is not implemented without ETK");
		return ARM_KO;
	}

	long curveId;

	ARM_ReferenceValue* modelFactor = NULL;
	ARM_ReferenceValue* modelFactorOld = NULL;

	CCString msg (" ");
	
	CCString xmlResponse;

	char sDate[11];
	Local_XLDATE2ARMDATE(AsOfDate,sDate);
	ARM_Date myDate(sDate);

	long retCode;

	try
	{
		modelFactor = etoolkit_getXMLModelParamFromSummit(myDate,
														model,
														type,
														factorName,
														ccy,
														index,
														cvName,
														calcModId);

		if ( modelFactor == NULL )
		{
		   result.setMsg("ARM_ERR: Curve is null");
			
           return(ARM_KO);
		}

		if(objId == -1)
		{
			CREATE_GLOBAL_OBJECT();

			curveId = LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)modelFactor);

			if (curveId == RET_KO)
			{
				if (modelFactor)
					delete modelFactor;
				modelFactor = NULL;

				result.setMsg ("ARM_ERR: Pb with inserting object");				
				return ARM_KO;
			}

			result.setLong(curveId);

			return ARM_OK;
		}
		else
		{
			modelFactorOld = (ARM_ReferenceValue*) LOCAL_PERSISTENT_OBJECTS->GetObject(objId);

			if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(modelFactorOld, ARM_REFERENCE_VALUE) == 1)
			{
				if (modelFactorOld)
				{
					delete modelFactorOld;
					modelFactorOld = NULL;
				}

				LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)modelFactor, objId);

				return ARM_OK;
			}
			else
			{
				result.setMsg ("ARM_ERR: previous object is not of a good type");
				return ARM_KO;
			}
		}
	}
	
    catch (Exception& x)
	{
		if (modelFactor)
			delete modelFactor;
		modelFactor = NULL;

		x.DebugPrint();

		ARM_RESULT();
	}

	return retCode;
}
