#include "ARMKernel/glob/firsttoinc.h"
#include "icm_str_models.h"

//------------------------------
#ifndef unix 
#define ITOA _itoa
#else
#define ITOA intToStr
#endif


// -----------------------------------------------------------------------------------------------------
// Fonction prenant en parametre un tableau de courbe de defaut et renvoyant la courbe de défaut moyenne
// -----------------------------------------------------------------------------------------------------
ICM_DefCurvStr* AverageCurveStr(const std::vector<const ICM_DefaultCurve*>& MultiCurve )
{
	int i,j, k =0;
	int sizeYT = 0;
	int indiceYT = 0;

	// Pour le cas ou certains plots ne seraient renséignés On recupere la courbe possédant le plus de plot
	for (i=0;i<MultiCurve.size();i++)
	{
		if (MultiCurve[i]->GetYearTerms()->GetSize() > sizeYT)
		{
			indiceYT = i;
			sizeYT = MultiCurve[i]->GetYearTerms()->GetSize();
		}
	}

	ARM_Vector* YF_plot = (ARM_Vector*)unconst( * MultiCurve[indiceYT]->GetYearTerms()).Clone();
	
	// Moyenne des recoveries de chaque tranche
	double Recovery = 0.;
	int NbRecov = 0;
	for (i=0;i<MultiCurve.size();i++)
	{
		if (MultiCurve[i]->GetRecovery())
		{
			Recovery += MultiCurve[i]->GetRecovery();
			NbRecov ++;
		}
	}
	Recovery /= NbRecov;

	// On recupere la courbe de taux (on suppose que c'est la meme pour toutes les courbes.
	ARM_ZeroCurve* ZCurv = (ARM_ZeroCurve*) MultiCurve[0]->GetZeroCurve()->Clone();

	// On recupere la moyenne des spreads au niveau des plots 
	ARM_Vector Spread = ARM_Vector(sizeYT - 1,0.); // -1 parce qu'on ne s'interesse pas au plot fictif en 0
	ARM_Vector NbCurvebyPlot = ARM_Vector(sizeYT - 1,0.);
	for (i=0;i<MultiCurve.size();i++)
	{
		// On recupere un vecteur de dimension = 5 + 1 (plot 0)
		ARM_Vector* YF_plotTemp = (ARM_Vector*) unconst(* MultiCurve[i]->GetYearTerms()).Clone();
		for (j=1;j<YF_plotTemp->GetSize();j++)
		{
			// on ne part pas du plot 0 parce qu'il correspond au plot fictif (plot(0) = 0.)
			k = 1;
			while (k<sizeYT)
			{ 
				if (fabs(YF_plotTemp->Elt(j) - YF_plot->Elt(k))<1E-8)
				{
					Spread.Elt(k-1) += MultiCurve[i]->GetRate(j);	// + 1 parce que dans la construction de la defcurve on a rajouté un plot fictif en 0
					NbCurvebyPlot.Elt(k-1) ++;
					k = sizeYT;
				}
				k++;
			}
		}
		if (YF_plotTemp)
		delete YF_plotTemp;
			YF_plotTemp = NULL;
	}

	for (int il=0;il<sizeYT-1;il++)
		Spread.Elt(il) /= NbCurvebyPlot.Elt(il);
		
		
	char* Label = "AVERAGE_CRB";
	
	// On ne recupere que la partie non nulle du vecteur de plot.
	ARM_Vector YearTerm = ARM_Vector(sizeYT-1,1.);
	 memcpy(YearTerm.GetElt(),YF_plot->GetElt()+1 ,sizeof(double)* (sizeYT-1));

	ICM_DefCurvStr* AVGCurve  = new ICM_DefCurvStr(&Spread,
 												   &YearTerm,
												   Recovery,
												   ZCurv,
												   Label);


	// Delete
	MyDelete(YF_plot);
	MyDelete(ZCurv);

	return AVGCurve;
}
// -----------------------------------------------------------------------------------------------------
// Fonction prenant en parametre un tableau de courbe de defaut et renvoyant la courbe de défaut moyenne
// -----------------------------------------------------------------------------------------------------
ICM_DefaultCurve* AverageCurvePWC(ICM_DefaultCurve** MultiCurve, int NbCurve)
{
	ARM_Date AsOf;
	AsOf.Today();
	int i,j, k =0;
	int sizeYT = 0;
	int indiceYT = 0;

	// Pour le cas ou certains plots ne seraient renséignés On recupere la courbe possédant le plus de plot
	for (i=0;i<NbCurve;i++)
	{
		if (MultiCurve[i]->GetYearTerms()->GetSize() > sizeYT)
		{
			indiceYT = i;
			sizeYT = MultiCurve[i]->GetYearTerms()->GetSize();
		}
	}

	ARM_Vector* YF_plot = (ARM_Vector*) unconst(* MultiCurve[indiceYT]->GetYearTerms()).Clone();
	
	// Moyenne des recoveries de chaque tranche
	double Recovery = 0.;
	int NbRecov = 0;
	for (i=0;i<NbCurve;i++)
	{
		if (MultiCurve[i]->GetRecovery())
		{
			Recovery += MultiCurve[i]->GetRecovery();
			NbRecov ++;
		}
	}
	Recovery /= NbRecov;

	// On recupere la courbe de taux (on suppose que c'est la meme pour toutes les courbes.
	ARM_ZeroCurve* ZCurv = (ARM_ZeroCurve*) MultiCurve[0]->GetZeroCurve()->Clone();
	string ccy = string(ZCurv->GetCurrencyUnit()->GetCcyName());
	// On recupere la moyenne des spreads au niveau des plots 
	ARM_Vector Spread = ARM_Vector(sizeYT - 1,0.); // -1 parce qu'on ne s'interesse pas au plot fictif en 0
	ARM_Vector NbCurvebyPlot = ARM_Vector(sizeYT - 1,0.);
	for (i=0;i<NbCurve;i++)
	{
		// On recupere un vecteur de dimension = 5 + 1 (plot 0)
		ARM_Vector* YF_plotTemp = (ARM_Vector*) unconst(* MultiCurve[i]->GetYearTerms()).Clone();
		for (j=1;j<YF_plotTemp->GetSize();j++)
		{
			// on ne part pas du plot 0 parce qu'il correspond au plot fictif (plot(0) = 0.)
			k = 1;
			while (k<sizeYT)
			{ 
				if (fabs(YF_plotTemp->Elt(j) - YF_plot->Elt(k))<1E-8)
				{
					Spread.Elt(k-1) += MultiCurve[i]->GetRate(j);	// + 1 parce que dans la construction de la defcurve on a rajouté un plot fictif en 0
					NbCurvebyPlot.Elt(k-1) ++;
					k = sizeYT;
				}
				k++;
			}
		}
		MyDelete(YF_plotTemp);
	}

	for (int il=0;il<sizeYT-1;il++)
		Spread.Elt(il) /= NbCurvebyPlot.Elt(il);
		
		
	char* Label = "AVERAGE_CRB";
	
	// On ne recupere que la partie non nulle du vecteur de plot.
	vector<string>  Term (sizeYT-1); 
	for (i=0;i<sizeYT-1;i++)
	{
		char buffer[65];
		ITOA(YF_plot->Elt(i+1),buffer,10);
		Term[i] = string(buffer) + string("Y");

	}

	ICM_Constant_Piecewise* AVGCurve = new ICM_Constant_Piecewise(AsOf,
							 Term,
							 (ARM_Vector*)Spread.Clone(),
							 Recovery,
							 (ARM_ZeroCurve*)ZCurv,
							 K_ADJUSTED,
							 K_ADJUSTED,
							 qCredit_Adjust20,
							 ccy, /*ARM_DEFAULT_COUNTRY,*/
							 Label,
							 false,
							 NULL,
							 //2NULL,
							 K_QUARTERLY,
							 MultiCurve[0]->GetCalibrationAlgo(),
							 MultiCurve[0]->GetCalibrationData(),
							 ARM_Currency(ccy.c_str()).GetCreditStartDateLag(),
							 ICM_Parameters()
							 );


	// Delete
	MyDelete(YF_plot);
	MyDelete(ZCurv);

	return AVGCurve;
}