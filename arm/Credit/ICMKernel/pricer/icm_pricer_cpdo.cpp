#include "firsttoinc.h"
#include "icm_pricer_cpdo.h"
#include "ICMKernel/inst/icm_cpdo.h"
#include "ICMKernel/util/icm_utils.h"
#include "ICMKernel/mod/modelmulticurves.h"

void ICM_Pricer_CPDO::Init()
{
	
	ICM_Pricer::Init();

	itsNbSimulations = 1;
	itsNbPath = 0; 

	itsRunning = 0.;
	itsBidMid = 0.000025;
	itsBidMidFloor = 0.000025;
	itsBidMidCap = 0.0001;
	itsRecovery = 0.4; 

	itsInitialLoss	= 0.05;
	itsSpreadSlope = 0.725; 
	itsStressRoll = 0.0005;

	itsAverageCashInDay = 0.; 
	itsAverageCashOutDay = 0.; 
	itsCountRebalance =0.;
	itsAverageFinalVal =0.;
	itsCashInPerc = 0.;
	itsCashOutPerc = 0.;

	itsIsDetailedSimulation = true; 
	itsIsRandomIR = true;
	itsIsRandomSpreads = true;
	itsIsRandomDefaults = true;
	
	itsDFs = NULL;
	itsDates = NULL; 

	SetName(ICM_PRICER_CPDO);
	
}

void ICM_Pricer_CPDO::Set(ARM_Security *sec, ARM_Object *mod, const ICM_Parameters &params, const ARM_Date &asof) 
{
	ICM_Pricer::Set(sec,mod,params,&asof);

	//if (params)
	//{

		ARM_Vector* a = params.GetColVect("a");
		if (a)	{itsa		= a->Elt(0);}
		
		ARM_Vector* IRVol = params.GetColVect("IRVol");
		if (IRVol)	{itsIRVol		= IRVol->Elt(0);}

		ARM_Vector* SpreadVol = params.GetColVect("SpreadVol");
		if (SpreadVol)	{itsSpreadVol		= SpreadVol->Elt(0);}

		ARM_Vector* NbSimulations = params.GetColVect("NbSimulations");
		if (NbSimulations)	{itsNbSimulations		= NbSimulations->Elt(0);}

		ARM_Vector* SimulationPath = params.GetColVect("SimulationPath");
		if (SimulationPath)	{itsSimulationPath		= SimulationPath->Elt(0);}

		ARM_Vector* Running = params.GetColVect("Running");
		if (Running)	{itsRunning		= Running->Elt(0);}

		ARM_Vector* BidMid = params.GetColVect("BidMid");
		if (BidMid)	{itsBidMid		= BidMid->Elt(0);}

		ARM_Vector* BidMidFloor = params.GetColVect("BidMidFloor");
		if (BidMidFloor)	{itsBidMidFloor		= BidMidFloor->Elt(0);}

		ARM_Vector* BidMidCap = params.GetColVect("BidMidCap");
		if (BidMidCap)	{itsBidMidCap		= BidMidCap->Elt(0);}

		ARM_Vector* SpreadSlope = params.GetColVect("SpreadSlope");
		if (SpreadSlope)	{itsSpreadSlope		= SpreadSlope->Elt(0);}
		
		ARM_Vector* InitialLoss = params.GetColVect("InitialLoss");
		if (InitialLoss)	{itsInitialLoss		= InitialLoss->Elt(0);}

		ARM_Vector* vStressRoll = params.GetColVect("StressRoll");
		if(vStressRoll)		{itsStressRoll		= vStressRoll->Elt(0);}

		double dIsDetailedSimulation ; 
		ARM_Vector* IsDetailedSimulation = params.GetColVect("IsDetailedSimulation");
		if (IsDetailedSimulation)	{dIsDetailedSimulation		= IsDetailedSimulation->Elt(0);}
		if (fabs(dIsDetailedSimulation) > DB_TOL)
			itsIsDetailedSimulation = true; 
		else
			itsIsDetailedSimulation = false;
		
		ARM_Vector* NumDetailedSimulation = params.GetColVect("NumDetailedSimulation");
		if (NumDetailedSimulation)	{itsNumDetailedSimulation		= NumDetailedSimulation->Elt(0);}
		
		double dIsRandomIR ; 
		ARM_Vector* IsRandomIR = params.GetColVect("IsRandomIR");
		if (IsRandomIR)	{dIsRandomIR		= IsRandomIR->Elt(0);}
		if (fabs(dIsRandomIR) > DB_TOL)
			itsIsRandomIR = true; 
		else
			itsIsRandomIR = false; 

		double dIsRandomSpreads ; 
		ARM_Vector* IsRandomSpreads = params.GetColVect("IsRandomSpreads");
		if (IsRandomSpreads)	{dIsRandomSpreads		= IsRandomSpreads->Elt(0);}
		if (fabs(dIsRandomSpreads) > DB_TOL)
			itsIsRandomSpreads = true; 
		else
			itsIsRandomSpreads = false; 

		double dIsRandomDefaults ; 
		ARM_Vector* IsRandomDefaults = params.GetColVect("IsRandomDefaults");
		if (IsRandomDefaults)	{dIsRandomDefaults		= IsRandomDefaults->Elt(0);}
		if (fabs(dIsRandomDefaults) > DB_TOL)
			itsIsRandomDefaults = true;
		else
			itsIsRandomDefaults = false;

		double dIsSRSched ; 
		ARM_Vector* IsSRSched = params.GetColVect("ISSRSched");
		if (IsSRSched)	{dIsSRSched		= IsSRSched->Elt(0);}
		if (fabs(dIsSRSched) > DB_TOL)
			itsIsSRSched = true; 
		else
			itsIsSRSched = false; 


		ARM_Vector* Duration = params.GetColVect("Duration");
		if (Duration)	{itsDuration		= Duration->Elt(0);}

		ARM_Vector* Recovery = params.GetColVect("Recovery");
		if (Recovery)	{itsRecovery		= Recovery->Elt(0);}


		// Used to enter deterministic paths 
		ARM_Vector* ImpliedSpread = params.GetColVect("ImpliedSpread");
		if (ImpliedSpread)	{itsDetImpliedSpread		= *ImpliedSpread;}
		

		ARM_Vector* InstRate = params.GetColVect("InstRate");
		if (InstRate)	
		{
			if (itsIsSRSched)
				itsDetInstRate		= *InstRate;
			else
				itsDetXt		= *InstRate;
		}


		ARM_Vector* Defaults = params.GetColVect("Defaults");
		if (Defaults)	{itsDetDefaults		= *Defaults;}

		/*
		ARM_Vector* ZCoupon = parameters->GetColVect("CRB_TAUX");
		if (ZCoupon)
		{
			int size = ZCoupon->GetSize();
			itsZeroCoupon.resize(size);
			for (int i = 0; i< size;i++)
			{itsZeroCoupon[i] = ZCoupon->Elt(i);}
		}*/
	//}
	
	ICM_cpdo* cpdo = dynamic_cast<ICM_cpdo*>(GetSecurity());
	itsMaturityDays = cpdo->GetCPDOMaturity() - asof; 
	itsSimulationSchedule = GenerateIntSch(0., itsMaturityDays, itsSimulationPath);
	
	itsNbPath = itsSimulationSchedule->GetSize(); 
	
}


void ICM_Pricer_CPDO::View(char* id, FILE* ficOut)
{	
	FILE* fOut;
	char  fOutName[200];

	if ( ficOut == NULL )
	{
	ARM_GetViewFile(ARM_VIEW_FILE, id, fOutName);
	
	   (void) unlink(fOutName);

       fOut = fopen(fOutName, "w"); 
    }
	else
	{
	fOut = ficOut;
	} 

	int size =0;

	ICM_Pricer::View(id, fOut);

	fprintf(fOut, "\t\t\t ------------------ Random Parameters ----------------- \n\n\n");
	fprintf(fOut, "\n");
	fprintf(fOut, "\t\t\t ------------------ Interest Rates ----------------- \n\n\n");

	fprintf(fOut, " a : %f \n",itsa);
	fprintf(fOut, " IR Vol : %f \n",itsIRVol);

	fprintf(fOut, "\n");
	fprintf(fOut, "\t\t\t ------------------ Spreads ----------------- \n\n\n");
	fprintf(fOut, " Spread Vol : %f \n",itsSpreadVol);

	fprintf(fOut, "\n\n");
	fprintf(fOut, "\t\t\t ------------------ Simulation Parameters ----------------- \n\n\n");
	fprintf(fOut, " Nb of simulations : %d \n",itsNbSimulations);
	fprintf(fOut, " nb days until maturity : %d \n",itsMaturityDays);
	fprintf(fOut, " Nb of paths in one simulation : %d \n",itsNbPath);
	fprintf(fOut, " Initial path of one simulation : %d \n",itsSimulationPath);


	fprintf(fOut, "\n\n");
	fprintf(fOut, "\t\t\t ------------------ Additional parameters ----------------- \n\n\n");
	fprintf(fOut, " Running Spread : %f \n",itsRunning);
	fprintf(fOut, " Bid Mid : %f \n",itsBidMid);
	fprintf(fOut, " Bid Mid Floor : %f \n",itsBidMidFloor);
	fprintf(fOut, " Bid Mid Cap : %f \n",itsBidMidCap);
	fprintf(fOut, " Recovery : %f \n",itsRecovery);

	fprintf(fOut, "\n\n");
	fprintf(fOut, "\t\t\t ------------------ Global Results ----------------- \n\n\n");
	fprintf(fOut, " Nb of Cash In : %d \n",itsNbCashIn);
	fprintf(fOut, " Nb of Cash Out : %d \n",itsNbCashOut);
	fprintf(fOut, " Perc of Cash In : %f \n",itsCashInPerc);
	fprintf(fOut, " Perc of Cash Out : %f \n",itsCashOutPerc);

	fprintf(fOut, " Nb of no in - no out : %d \n",itsNbNoInNoOut);

	fprintf(fOut, " Average Day of Cash In : %f \n",itsAverageCashInDay);
	fprintf(fOut, " Average Day of Cash Out : %f \n",itsAverageCashOutDay);
	fprintf(fOut, " Average Nb of rebalancing (for all simulations) : %f \n",itsCountRebalance);
	fprintf(fOut, " Average final value if no cash in and no cash out : %f \n",itsAverageFinalVal);

	// Just if we want the details of one simulation 
	fprintf(fOut, "\n\n");
	fprintf(fOut, "\t\t\t ------------------ Outputs for one simulation ----------------- \n\n\n");

    fprintf(fOut, "%14s\t%14s\t%14s\t%14s\t\n",
         "Risky Coupon Date",
		 "Non Risky Coupon Date ", 
		 "Roll Date",
		 "DF") ;
	
	int MaxSize = MAX(MAX(itsRiskyCouponDate.GetSize(),itsNoRiskyCouponDate.GetSize()),itsRollDate.GetSize());
	MaxSize = MAX(MaxSize,itsDF0.GetSize());

	for (int i = 0; i < MaxSize ; i++)
	{
		double dd1,dd2,dd3,dd4;
		if (i<itsRiskyCouponDate.GetSize()) dd1 = itsRiskyCouponDate.Elt(i); else dd1=-999;
		if (i<itsNoRiskyCouponDate.GetSize()) dd2 = itsNoRiskyCouponDate.Elt(i); else dd2=-999;
		if (i<itsRollDate.GetSize()) dd3 = itsRollDate.Elt(i); else dd3=-999;
		if (i<itsDF0.GetSize()) dd4 = itsDF0.Elt(i); else dd4=-999;

		fprintf(fOut, "%10.4lf\t%10.4lf\t%10.4lf\t%10.4lf\t\n",
         dd1,
		 dd2, 
		 dd3,
		 dd4) ;
	
	}

	fprintf(fOut, "\n\n\n\n\n\n");

	if (itsIsDetailedSimulation)
	{	
    fprintf(fOut, "%14s\t%14s\t%14s\t%14s\t	%14s\t%14s\t%14s\t%14s\t%14s\t%14s\t%14s\t%14s\t%14s\t%14s\t%14s\t%14s\t\n",
         "Simulation Dates",
		 "Risky Coupon Date ", 
		 "Non Risky Coupon Date",
		 "Roll Date",
		 "Forward Rate",
		 "Xt",
		 "Implied Spread",
		 "Inst Rate ",
		 "Default",
		 "Target Value",
		 "Risky Value",
		 "MTM",
		 "Cash",
		 "Notional",
		 "PNL",
		 "DF");

	double d1,d2,d3,d4,d5,d6,d7,d8,d9,d10,d11,d12,d13,d14,d15,d16;
		
		for (i = 0; i < itsSimulationSchedule->GetSize(); i++)
		{
			if (itsSimulationSchedule->GetSize()>0) d1=itsSimulationSchedule->Elt(i); else d1 = 0; 
			if (itsIsRiskyCouponDate.GetSize()>0) d2=itsIsRiskyCouponDate.Elt(i); else d2 = 0; 
			if (itsIsNoRiskyCouponDate.GetSize()>0) d3=itsIsNoRiskyCouponDate.Elt(i); else d3 = 0; 
			if (itsIsRollDate.GetSize()>0) d4=itsIsRollDate.Elt(i); else d4 = 0; 
			if (itsFwdRate.GetSize()>0) d5=itsFwdRate.Elt(i); else d5 = 0; 
			if (itsXt.GetSize()>0) d6=itsXt.Elt(i); else d6 = 0; 
			if (itsImpliedSpread.GetSize()>0) d7=itsImpliedSpread.Elt(i); else d7 = 0; 
			if (itsInstRate.GetSize()>0) d8=itsInstRate.Elt(i); else d8 = 0; 
			if (itsDefaults.GetSize()>0) d9=itsDefaults.Elt(i); else d9 = 0; 
			if (itsTargetValue.GetSize()>0) d10=itsTargetValue.Elt(i); else d10 = 0; 
			if (itsRiskyValue.GetSize()>0) d11=itsRiskyValue.Elt(i); else d11 = 0; 
			if (itsMTM.GetSize()>0) d12=itsMTM.Elt(i); else d12 = 0; 
			if (itsCash.GetSize()>0) d13=itsCash.Elt(i); else d13 = 0; 
			if (itsNotional.GetSize()>0) d14=itsNotional.Elt(i); else d14 = 0; 
			if (itsPNL.GetSize()>0) d15=itsPNL.Elt(i); else d15 = 0; 
			if (itsDF.GetSize()>0) d16=itsDF.Elt(i); else d16 = 0; 
			
			fprintf(fOut,"%10.4lf\t%10.4lf\t%10.4lf\t%10.4lf\t%10.4lf\t%10.10lf\t%10.4lf\t%10.4lf\t%10.4lf\t%10.4lf\t%10.4lf\t%10.4lf\t%10.4lf\t%10.4lf\t%10.4lf\t%10.10lf\t\n", 
				d1,	d2, d3,	d4,	d5,
				d6,	d7,	d8,	d9,	d10,
				d11,d12,d13,d14,d15,d16);

			for (int j = 0; j < itsDates[i]->GetSize()  ; j++)
				fprintf(fOut,"%10.4lf\t", itsDates[i]->Elt(j));	

			fprintf(fOut, "\n");

			for (j = 0; j < itsDFs[i]->GetSize()  ; j++)
				fprintf(fOut,"%10.10lf\t", itsDFs[i]->Elt(j));	

			fprintf(fOut, "\n");


		}
	}

	if ( ficOut == NULL )
	{
		fclose(fOut);
	}


}


// Used to merge the schedule present in security "sec" with itsSimulationSchedule 
void ICM_Pricer_CPDO::AddSchedule(const ICM_Security& sec, const ARM_Date& dAsOf) 
{

	ARM_Vector Schedule = sec.GetAccEndDates();
	int size=0;
	size = Schedule.GetSize();
	ARM_Vector dAsOfs(size,dAsOf.GetJulian());

	// Transforming dates in Nb of Days from AsOf  
	Schedule -= dAsOfs; 

	ARM_Vector* Result = NULL;
    MergeDates(&Result, *itsSimulationSchedule, Schedule);

	if (itsSimulationSchedule) 
		delete itsSimulationSchedule;

	itsSimulationSchedule = Result; 

}

void ICM_Pricer_CPDO::GenerateCpnRollDates(const ICM_Security& sec, 
										   const ARM_Date& dAsOf, 
										   const ARM_Date& dPastFixingDate, 
										   ARM_Vector& IsDate, 
										   ARM_Vector& CpnDate) 
{

	ARM_Vector Schedule = sec.GetAccEndDates();
	int size=0;
	size = Schedule.GetSize();
	ARM_Vector dAsOfs(size,dAsOf.GetJulian());

	// Transforming dates in Nb of Days from AsOf  
	Schedule -= dAsOfs; 

	IsDate.Resize(itsNbPath); 
	CpnDate.Resize(size+1); 

	CpnDate.Elt(0) = dAsOf.GetJulian() - dPastFixingDate.GetJulian(); 
	for(int j=0; j<size; j++)
		CpnDate.Elt(j+1) = Schedule.Elt(j);	

	int p=0; 
	for(int i=0 ; i<itsNbPath; i++)
	{
		if (p<size)
		{
			if (fabs(itsSimulationSchedule->Elt(i) - Schedule.Elt(p))<DB_TOL) 
			{
				IsDate[i] = 1;
				p++; 
			}
			else
				IsDate[i] = 0;

		}
		else
			IsDate[i] = 0;
	}

}

// Creates the simulation schedule from 
// MonteCarlo Schedule and
// the schedules of coupons and Roll dates 
// Initializes all the schedules with the good size itsNbPath
void ICM_Pricer_CPDO::CreateSimulationSchedule()
{
	// Getting security
	ICM_cpdo* cpdo = dynamic_cast<ICM_cpdo*>(GetSecurity());

	// Getting maturity of CPDO
	ARM_Date dAsOf = GetAsOfDate() ; 

	
	const ICM_Security& securityRisky = cpdo->GetRiskyLeg().GetCreditInfosRef();
	const ICM_Security& securityRoll = cpdo->GetRollLeg().GetCreditInfosRef();
	const ICM_Security& securityNoRisky = cpdo->GetNoRiskyLeg().GetCreditInfosRef();
	
	// Merging the schedules of coupon and roll dates with the schedule of Monte Carlo 
	AddSchedule(securityRisky, dAsOf);	
	AddSchedule(securityRoll, dAsOf);	
	AddSchedule(securityNoRisky, dAsOf);	

	// Storing Schedule in cache
	adoptVectorValue("SIMULATION_SCHEDULE", (ARM_Vector*)itsSimulationSchedule->Clone());

	// Initialization of itsNbPath
	itsNbPath = itsSimulationSchedule->GetSize(); 

	//Populate Coupon and Roll Schedules 
	GenerateCpnRollDates(securityRisky,dAsOf,dAsOf,itsIsRiskyCouponDate, itsRiskyCouponDate);
	GenerateCpnRollDates(securityRoll,dAsOf,dAsOf,itsIsRollDate, itsRollDate);
	GenerateCpnRollDates(securityNoRisky,dAsOf,dAsOf,itsIsNoRiskyCouponDate, itsNoRiskyCouponDate);

	// Storing Schedules of Risky and non risky security in cache
	adoptVectorValue("SCHED_RISKY", (ARM_Vector*)itsRiskyCouponDate.Clone());
	adoptVectorValue("SCHED_ROLL", (ARM_Vector*)itsRollDate.Clone());
	adoptVectorValue("SCHED_NORISKY", (ARM_Vector*)itsNoRiskyCouponDate.Clone());

	adoptVectorValue("SCHED_ISRISKY", (ARM_Vector*)itsIsRiskyCouponDate.Clone());
	adoptVectorValue("SCHED_ISROLL", (ARM_Vector*)itsIsRollDate.Clone());
	adoptVectorValue("SCHED_ISNORISKY", (ARM_Vector*)itsIsNoRiskyCouponDate.Clone());


	if (itsIsDetailedSimulation) 
	{
		if (itsNumDetailedSimulation > itsNbSimulations)
			itsNumDetailedSimulation = 1;
	}

	itsDFs = new ARM_Vector*[itsNbPath];
	itsDates = new ARM_Vector*[itsNbPath];
	for (int i = 0; i<itsNbPath; i++) 
	{
		itsDFs[i] = NULL;
		itsDates[i] = NULL;
	}

	// Initialization of vectors of size itsNbPath
	itsFwdRate.Resize(itsNbPath); 
	itsXt.Resize(itsNbPath);

	itsImpliedSpread.Resize(itsNbPath);
	itsInstRate.Resize(itsNbPath);
	itsDefaults.Resize(itsNbPath);

	//if (!(itsIsRandom))
	//	PopulateVectors();

	itsRiskyValue.Resize(itsNbPath);		
	itsTargetValue.Resize(itsNbPath);
	itsMTM.Resize(itsNbPath);
	itsCash.Resize(itsNbPath);
	itsPNL.Resize(itsNbPath);
	itsNotional.Resize(itsNbPath);


}

// This is to avoid that a user tries to input a vector which has a size different 
// from the Simulation Schedule size
void ICM_Pricer_CPDO::TestSizeNonRandomVectors()
{
	if (! itsIsRandomIR) 
		if (itsIsSRSched)
		{
			if (itsDetInstRate.GetSize() != itsSimulationSchedule->GetSize())
			{
				unsetVectorFlg("SIMULATION_SCHEDULE"); 			
				unsetVectorFlg("SCHED_RISKY"); 			
				unsetVectorFlg("SCHED_ISRISKY"); 			
				unsetVectorFlg("SCHED_ROLL"); 			
				unsetVectorFlg("SCHED_ISROLL"); 			
				unsetVectorFlg("SCHED_NORISKY"); 			
				unsetVectorFlg("SCHED_ISNORISKY"); 			
				ICMTHROW(ERR_INVALID_ARGUMENT,"ICM_Pricer_CPDO :: pb in size of Interest Rate path");
			}
		}
		else
		{
			if (itsDetXt.GetSize() != itsSimulationSchedule->GetSize())
			{
				unsetVectorFlg("SIMULATION_SCHEDULE"); 			
				unsetVectorFlg("SCHED_RISKY"); 			
				unsetVectorFlg("SCHED_ISRISKY"); 			
				unsetVectorFlg("SCHED_ROLL"); 			
				unsetVectorFlg("SCHED_ISROLL"); 			
				unsetVectorFlg("SCHED_NORISKY"); 			
				unsetVectorFlg("SCHED_ISNORISKY"); 			
				ICMTHROW(ERR_INVALID_ARGUMENT,"ICM_Pricer_CPDO :: pb in size of Interest Rate path");
			}
		}

	if (! itsIsRandomSpreads) 
		if (itsDetImpliedSpread.GetSize() != itsSimulationSchedule->GetSize())
		{
			unsetVectorFlg("SIMULATION_SCHEDULE"); 			
			unsetVectorFlg("SCHED_RISKY"); 			
			unsetVectorFlg("SCHED_ISRISKY"); 			
			unsetVectorFlg("SCHED_ROLL"); 			
			unsetVectorFlg("SCHED_ISROLL"); 			
			unsetVectorFlg("SCHED_NORISKY"); 			
			unsetVectorFlg("SCHED_ISNORISKY"); 			
			ICMTHROW(ERR_INVALID_ARGUMENT,"ICM_Pricer_CPDO :: pb in size of Spreads path");
		}

	if (! itsIsRandomDefaults) 
		if (itsDetDefaults.GetSize() != itsSimulationSchedule->GetSize())
		{
			unsetVectorFlg("SIMULATION_SCHEDULE"); 			
			unsetVectorFlg("SCHED_RISKY"); 			
			unsetVectorFlg("SCHED_ISRISKY"); 			
			unsetVectorFlg("SCHED_ROLL"); 			
			unsetVectorFlg("SCHED_ISROLL"); 			
			unsetVectorFlg("SCHED_NORISKY"); 			
			unsetVectorFlg("SCHED_ISNORISKY"); 			
			ICMTHROW(ERR_INVALID_ARGUMENT,"ICM_Pricer_CPDO :: pb in size of Defaultspath");
		}

}


void ICM_Pricer_CPDO::ComputeInstRate()
{
	// Compute inst rate 
	for (int i=0; i<itsInstRate.GetSize(); i++)
		itsInstRate.Elt(i) = 0.; 
	
	itsInstRate += itsXt  ; 
	itsInstRate += itsFwdRate; 

}

void ICM_Pricer_CPDO::ComputeXt()
{
	// Compute inst rate 
	for (int i=0; i<itsXt.GetSize(); i++)
		itsXt.Elt(i) = 0.; 
	
	itsXt += itsInstRate  ; 
	itsXt -= itsFwdRate; 
 
}

/*
void ICM_Pricer_CPDO::SetRandoms()
{
	// Compute Xt

	// Compute inst rate 
	itsInstRate = itsXt + itsFwdRate ; 
}
*/

void ICM_Pricer_CPDO::InitSimulation()
{
	const ARM_Vector& cSimulationSchedule = getVectorValue("SIMULATION_SCHEDULE");

	ICM_ModelMultiCurves* model = dynamic_cast<ICM_ModelMultiCurves*>(GetModel());

	if (model)
	{	
		ARM_ZeroCurve* DFCurve = model->GetZeroCurve(); 
		
		itsDF.Resize(itsNbPath); 
		itsDF0.Resize(itsNoRiskyCouponDate.GetSize()); 

		itsDF.Elt(0) = 1; 
		itsDF0.Elt(0) = 1; 
		int p = 1; 

		// Gets the vector of DF for the days of Coupons 
		// Computing a vector of DF for all days of simulation itsDF vector 
		
		for (int i = 1; i < itsNbPath ; i++)
		{	
			double NbDaysFromAsOf = cSimulationSchedule.Elt(i);
			itsDF.Elt(i) = DFCurve->DiscountPrice(NbDaysFromAsOf / 365.); 
			if (fabs(NbDaysFromAsOf - itsNoRiskyCouponDate.Elt(p)) < DB_TOL) 
			{
				if (p < itsNoRiskyCouponDate.GetSize()) 
					itsDF0.Elt(p) = DFCurve->DiscountPrice(NbDaysFromAsOf / 365.);
				else 
					ICMTHROW(ERR_INVALID_ARGUMENT,"ICM_Pricer_CPDO :: pb in size of DF curve");

				p++; 
			}
		}

		// Computing a vector of Fwd for all days of simulation
		for (i = 0; i<itsNbPath-1; i++)
		{
			itsFwdRate.Elt(i) = - 365./(cSimulationSchedule.Elt(i+1)-cSimulationSchedule.Elt(i)) * log(itsDF.Elt(i+1)/itsDF.Elt(i)) ; 
		}
		itsFwdRate.Elt(itsNbPath - 1) = itsFwdRate.Elt(itsNbPath - 2); 

		itsDates[0] = (ARM_Vector*) itsNoRiskyCouponDate.Clone(); 
		itsDFs[0] = (ARM_Vector*) itsDF0.Clone(); 


	}
	else
		ICMTHROW(ERR_INVALID_ARGUMENT,"ICM_Pricer_CPDO :: unable to get model");

}

/*
double Curve::ComputeDF(const int& NumSimul)
{

	// used to compute a Discount Factor for each path of simulation 
	int Index = 0 ; 
	if (its NbJours < itsMatu[itsSizeCurve - 1])
	{
		while (itsMatu[Index] <= NbJours) 
			Index ++ ; 
		
		return itsDF[Index-1] + 
			(double)(NbJours - itsMatu[Index - 1]) / 
				(double)(itsMatu[Index] - itsMatu[Index - 1]) * (itsDF[Index] - itsDF[Index - 1]);
	}
	else
	{
		return itsDF[itsSizeCurve - 1]; 
	}

}
*/

void ICM_Pricer_CPDO::Simulate()
{
	itsNbCashIn = 0; 
	itsNbCashOut = 0; 
	itsNbNoInNoOut = 0;

	itsCountRebalance = 0.; 
	itsAverageCashInDay = 0.; 
	itsAverageCashOutDay = 0.; 
	itsAverageFinalVal = 0.;

	TestSizeNonRandomVectors();

	ICM_cpdo* cpdo = dynamic_cast<ICM_cpdo*>(GetSecurity());
	
	itsSimulDayTarget.Resize(itsNbSimulations);	
	itsSimulPNL.Resize(itsNbSimulations);	
	itsSimulIsCashIn.Resize(itsNbSimulations);
	itsSimulIsCashOut.Resize(itsNbSimulations);	
	itsValo.Resize(itsNbSimulations);


	if (cpdo) 
	{
		double Target = cpdo->GetTarget();
		double UFFees = cpdo->GetUFFees();
		double RunningFees = cpdo->GetRunningFees();
		double InitialValo = cpdo->GetInitialValo(); 
		double ExpoOnV = cpdo->GetVExpo();
		double ExpoOnV0 = cpdo->GetV0Expo();
		double Alpha = cpdo->GetAlpha();
		double Beta = cpdo->GetBeta(); 
		double Desactivation = cpdo->GetDesactivation();
		int NbAssets = cpdo->GetNbAssets(); 
		
		InitSimulation();

		//Loop on simulations 
		for (int i = 0 ; i<itsNbSimulations; i++)
		{
			if (itsIsRandomIR)
			{
				// to be changed to generate Random Xt and then Inst Rate
				ComputeInstRate();
			}
			else
			{
				if (itsIsSRSched)
				{
					itsInstRate = itsDetInstRate; 
					ComputeXt();  
				}
				else
				{
					itsXt = itsDetXt; 
					ComputeInstRate();  
				}
			}

			if (itsIsRandomSpreads) 
				itsImpliedSpread = itsDetImpliedSpread; // to be changed to generate Random Implied Spreads
			else
				itsImpliedSpread = itsDetImpliedSpread;

			if (itsIsRandomDefaults) 
				itsDefaults = itsDetDefaults; // to be changed to generate Random defaults 
			else
				itsDefaults = itsDetDefaults; 

			SimulateOnePath(i, Target, UFFees, RunningFees, 
				InitialValo, ExpoOnV, ExpoOnV0, 
				Alpha, Beta, Desactivation,NbAssets); 
		}

		ComputeGlobalOutputs(); 

	}
	else
		ICMTHROW(ERR_INVALID_ARGUMENT,"ICM_Pricer_CPDO :: unable to simulate , No security associated");


}

// Function used to generate one single simulation 
void ICM_Pricer_CPDO::SimulateOnePath(const int& NumSimu, 
									  const double& Target,
									  const double& UFFees,
									  const double& RunningFees,
									  const double& InitialValo,
									  const double& ExpoOnV,
									  const double& ExpoOnV0,
									  const double& Alpha,
									  const double& Beta,
									  const double& Desactivation,
									  const int& NbAssets)

{
	const ARM_Vector& cSimulationSchedule = getVectorValue("SIMULATION_SCHEDULE");

	double TargetPV = 0.; 
	double PNL = 0.; 
	double PeriodLen = 0.; 
	double RefMat = 4.75; 

	int NbDefaultNames=0;

	int NumRiskyCpn = 1; 
	int NumNoRiskyCpn = 1;
	int NumRoll = 0; 

	double NoRiskyPeriod; 
	double RiskyPeriod;
	double RollPeriod; 

	bool IsNoRiskyCpnDate;

	NoRiskyPeriod = itsNoRiskyCouponDate.Elt(NumNoRiskyCpn) - itsNoRiskyCouponDate.Elt(NumNoRiskyCpn-1);
	RiskyPeriod = itsRiskyCouponDate.Elt(NumRiskyCpn) - itsRiskyCouponDate.Elt(NumRiskyCpn-1);
	RollPeriod = 0.5 ; // (itsRollDate.Elt(NumRoll) - itsRollDate.Elt(NumRoll-1)) / 365. ;

	double NonRiskyCpn = 0.03 ; 
	double Running  = itsRunning ;
	double BidMid = itsBidMid ; 
	double BidMidCap = itsBidMidCap ;			// Cap Bid Mid : 1 bp
	double BidMidFloor = itsBidMidFloor;	// Floor Bid Mid : 0.25 bp


	itsSimulDayTarget.Elt(NumSimu) = 0; 

	// Set randoms for one simulation or set paths entered by user 
	// SetRandoms();

	TargetPV = ComputeTarget(NoRiskyPeriod,
								false, 
								InitialValo,
								Target,
								RunningFees,
								itsNoRiskyCouponDate,
								itsDF0); 


	double InstRate = itsInstRate.Elt(0); 

	// Getting UpFront Fees
	double V = InitialValo * (1 - UFFees); 
	PNL += InitialValo * UFFees; 

	// Taking into account an initial loss
	V -= itsInitialLoss * InitialValo ;

	double ExpoCap = MIN(ExpoOnV * V, ExpoOnV0 * InitialValo) ;
	double Ntheo = MIN(ExpoCap, Alpha * (TargetPV - V) + Beta) ; 

	double TimeToRoll = itsRollDate.Elt(NumRoll) / 365. ; //PeriodLen * RollDates[NumRoll] ; 
	double ImpliedSpread = itsImpliedSpread.Elt(0) * (1- (RollPeriod  - TimeToRoll)/2 * (1-itsSpreadSlope));

	double Duration = (1 - exp(-(InstRate + ImpliedSpread) * (TimeToRoll + RefMat))) / (InstRate + ImpliedSpread); 

	// Initialization of variables 
	// We assume we start the simulation on a coupon date so that Accrued = 0 
	double Notional = Ntheo ; 
	double Accrued = 0.; 
	double MTM = (Running - ImpliedSpread) * Duration * Notional + Accrued ; 
		
	// Initial Cash minus UpFrontFees (UFFee) minus execution costs  
	double Cash = V + (ImpliedSpread - BidMid - Running) * Duration * Notional - Accrued;

	// Rebalancing parameters
	bool IsRebalance = false; 
	int CountRebalance = 0 ; 
	double TestRebal = Notional / Ntheo; 

	itsIsCashIn = false; 
	itsIsCashOut = false; 

	int DayCashIn = 0; 
	int DayCashOut = 0; 

	// if there is at least one default on the first period 
	if (itsDefaults.Elt(0)) 
	{
		double Jump = 1; 
		for (int ploop = 0 ; ploop < itsDefaults.Elt(0) ; ploop ++)
			Jump *= (1-1/(double)(NbAssets - NbDefaultNames - ploop)) ; 
		Notional *= Jump; 
		NbDefaultNames += itsDefaults.Elt(0); 
	}


	// Loop on all paths 
	// Main point of simulation 

	// For logs 
	if ((itsIsDetailedSimulation) && itsNumDetailedSimulation == NumSimu+1) 
	{
		itsRiskyValue.Elt(0) = V; 
		itsTargetValue.Elt(0) = TargetPV; 
		itsMTM.Elt(0) = MTM; 
		itsCash.Elt(0) = Cash; 
		itsPNL.Elt(0) = PNL; 
		itsNotional.Elt(0) = Notional; 
	}
	
	for (int i=1; i<itsNbPath; i++)
	{
		PeriodLen = (cSimulationSchedule.Elt(i) - cSimulationSchedule.Elt(i-1)) / 365. ; 
		
		int SizeDF = GetSizeDFCurve(i); 
		ARM_Vector DatesVector(SizeDF,0.);
		ARM_Vector DFVector(SizeDF,0.); 


		ComputeDiscountFactorCurve(i,DatesVector,DFVector); 

		itsDates[i] = (ARM_Vector*) DatesVector.Clone(); 
		itsDFs[i] = (ARM_Vector*) DFVector.Clone(); 

		if (fabs(itsIsNoRiskyCouponDate.Elt(i))<DB_TOL)
			IsNoRiskyCpnDate = true;
		else 
			IsNoRiskyCpnDate = false;

		// Computing target 
		TargetPV = ComputeTarget(NoRiskyPeriod,
								IsNoRiskyCpnDate,
								InitialValo,
								Target,
								RunningFees,
								DatesVector,
								DFVector); 
			

		// Testing if Cash In or Cash Out
		if ((!itsIsCashIn) && (!itsIsCashOut))
		{
			TimeToRoll -= PeriodLen; 
			ImpliedSpread = itsImpliedSpread.Elt(i) * (1- (RollPeriod - TimeToRoll)/2 * (1-itsSpreadSlope));
		
			if (NumRoll > 0)
					ImpliedSpread = MAX(ImpliedSpread-itsStressRoll, 0.0001);

			Accrued += Notional * Running * PeriodLen; 

			Cash  *= (1 + InstRate * PeriodLen);
			PNL  *= (1 + InstRate * PeriodLen);

			// In Case of Default
			if (itsDefaults.Elt(i)) 
			{
				double Jump = 1; 
				double Rebal ; 
				Jump = (1 - (double)itsDefaults.Elt(i)/(double)(NbAssets-NbDefaultNames)) ;

				Cash -= Notional * (1 - itsRecovery) / (double)(NbAssets-NbDefaultNames) ;
				Rebal = Running * Notional * (1-Jump) * (cSimulationSchedule.Elt(i) - itsRiskyCouponDate.Elt(NumRiskyCpn-1) ) /365.; // (double)(p - CouponDates[NumCoupon-1]) / 262.;
				Cash += Rebal; 
				Accrued -= Rebal; 

				Notional *= Jump; 
			
				NbDefaultNames += itsDefaults.Elt(i); 
			}


			// Non Risky Coupons
			if (itsIsNoRiskyCouponDate.Elt(i))
			{	
				// Coupon paid (Eurib + Target) 
				Cash -= (NonRiskyCpn + Target + RunningFees) * InitialValo * NoRiskyPeriod/365. ; // * (double)CouponPeriod / 262.; 
				PNL += RunningFees * InitialValo * NoRiskyPeriod/365.; //; (double)CouponPeriod / 262.; 

				NumNoRiskyCpn ++;
				// New Coupon period 
				
				if (NumNoRiskyCpn < itsNoRiskyCouponDate.GetSize()) 
					NoRiskyPeriod = (itsNoRiskyCouponDate.Elt(NumNoRiskyCpn) - itsNoRiskyCouponDate.Elt(NumNoRiskyCpn-1)) ; 

				// Setting EURIB 3M 
				NonRiskyCpn = SetEURIB3M(DatesVector,DFVector);  

			}

			InstRate = itsInstRate.Elt(i); 
			Duration = (1 - exp(-(InstRate + ImpliedSpread) * (TimeToRoll + RefMat))) / (InstRate + ImpliedSpread);

	
			// Risky Coupons
			if (itsIsRiskyCouponDate.Elt(i))
			{	
				// Coupon received on CDS on index
				Cash += Running * Notional * RiskyPeriod/365. ; //(double)CouponPeriod / 262.;

				NumRiskyCpn ++;
				// New Coupon period 
				
				if (NumRiskyCpn < itsRiskyCouponDate.GetSize()) 
					RiskyPeriod = (itsRiskyCouponDate.Elt(NumRiskyCpn) - itsRiskyCouponDate.Elt(NumRiskyCpn-1)) ; 

				Accrued = 0.; 

			}

			MTM = (Running - ImpliedSpread) * Duration * Notional + Accrued ;


			// Roll Date 
			if (itsIsRollDate.Elt(i))
			{
				NumRoll ++; 
				
				if (NumRoll < itsRollDate.GetSize())
					RollPeriod = (itsRollDate.Elt(NumRoll) - itsRollDate.Elt(NumRoll-1)) / 365.  ;  
				else 
					RollPeriod = 0.5 ; //262. * PeriodLen ; 
				
				TimeToRoll = RollPeriod ; 

				double NewImpliedSpread = itsImpliedSpread.Elt(i) ;
				NewImpliedSpread = MAX(NewImpliedSpread-itsStressRoll, 0.0001);
				double NewDuration = (1 - exp(-(InstRate + NewImpliedSpread) * (TimeToRoll + RefMat))) / (InstRate + NewImpliedSpread);  
				if ( NumRoll == 1) 
					Cash += (Running - (ImpliedSpread + BidMid + itsStressRoll)) * Duration * Notional ; 
				else
					Cash += (Running - (ImpliedSpread + BidMid)) * Duration * Notional ; 
				Running = NewImpliedSpread ;
				BidMid = MIN(BidMidCap, MAX(0.005 * Running,BidMidFloor));

				Cash -= (Running - (NewImpliedSpread - BidMid)) * NewDuration * Notional ; 

				MTM = (Running - NewImpliedSpread) * NewDuration * Notional + Accrued ;
				NbDefaultNames = 0; 

				Duration = NewDuration ; 
				ImpliedSpread = NewImpliedSpread; 

				}

			V = MTM + Cash; 

			// Testing if Cash In 
			if ((V-TargetPV) > BidMid * Notional * Duration) 
				{
					itsIsCashIn = true ; 
					
					itsNbCashIn++; 
					DayCashIn = i; 
					itsAverageCashInDay += DayCashIn; 

					itsSimulDayTarget.Elt(NumSimu) = cSimulationSchedule.Elt(i); 
					PNL += V - TargetPV - BidMid * Notional * Duration; 
					
					itsSimulDayTarget.Elt(NumSimu) = i; 

				}
				
			// Testing if Cash Out
			if (V < Desactivation) 
				{
					itsIsCashOut = true ;

					itsNbCashOut++; 
					DayCashOut = i;
					itsAverageCashOutDay += DayCashOut; 

					V -= BidMid * Notional * Duration;
					Cash = V; 
					MTM = 0.; 
					if (!(V>0.))
					{
						PNL += V ; 
						V =0.;
					}
				}
	
			// Testing if there is need to rebalance 
			ExpoCap = MIN(ExpoOnV * V, ExpoOnV0 * InitialValo); 						
			Ntheo = MIN(ExpoCap , Alpha * (TargetPV-V) + Beta); 
			TestRebal = Notional / Ntheo; 
			IsRebalance = ((TestRebal < 0.75) || (TestRebal > 1.25)); 
			
				if ((IsRebalance) && (!itsIsCashIn) && (!itsIsCashOut))  
				{
					double NRebal = Ntheo - Notional; 
					int Sgn = (TestRebal > 1 ) ? 1 : -1; 
					double SoultRebal = (ImpliedSpread + Sgn * BidMid - Running) * Duration * NRebal; 
					double AccrRebal = Running * NRebal * (i - itsRiskyCouponDate.Elt(NumRiskyCpn-1)) /365. ; // adding Year Fraction 
					
					SoultRebal -= AccrRebal; 

					Notional = Ntheo ; 
					
					Accrued += AccrRebal; 
					MTM = (Running - ImpliedSpread) * Duration * Notional + Accrued ;
					
					Cash = Cash + SoultRebal ; 
					V = MTM + Cash; 
					CountRebalance ++; 
					
				}
			
			// For logs 
			if ((itsIsDetailedSimulation) && itsNumDetailedSimulation == NumSimu+1) 
			{
				itsRiskyValue.Elt(i) = V; 
				itsTargetValue.Elt(i) = TargetPV; 
				itsMTM.Elt(i) = MTM; 
				itsCash.Elt(i) = Cash; 
				itsPNL.Elt(i) = PNL; 
				itsNotional.Elt(i) = Notional; 
			}

		}
		else // Cash In or Cash Out 
		{
			// Cash In Case
			if (itsIsCashIn)
			{
				V = TargetPV; 
				PNL *= (1 + InstRate * PeriodLen); 

				if (itsIsNoRiskyCouponDate.Elt(i))
				{
					PNL += RunningFees * InitialValo * NoRiskyPeriod / 365.; // adding Year Fraction 
					NumNoRiskyCpn ++;
					// New Coupon period 
				
					if (NumNoRiskyCpn < itsNoRiskyCouponDate.GetSize()) 
						NoRiskyPeriod = (itsNoRiskyCouponDate.Elt(NumNoRiskyCpn) - itsNoRiskyCouponDate.Elt(NumNoRiskyCpn-1)) ; 

					NonRiskyCpn = SetEURIB3M(DatesVector,DFVector); // EUR3M or INFLATION ... 
				}


			}
			
			// Cash Out Case
			if (itsIsCashOut)
			{
				if ((V>0.))			
				{
					V *= (1 + InstRate * PeriodLen); 
					PNL *= (1 + InstRate * PeriodLen); 
				}
			}
			
			// For logs 
			if ((itsIsDetailedSimulation) && itsNumDetailedSimulation == NumSimu+1) 
			{
				itsRiskyValue.Elt(i) = V; 
				itsTargetValue.Elt(i) = TargetPV; 
				itsMTM.Elt(i) = MTM; 
				itsCash.Elt(i) = Cash; 
				itsPNL.Elt(i) = PNL; 
				itsNotional.Elt(i) = Notional; 
			}

			// Always compute inst rate for next period
			InstRate = itsInstRate.Elt(i); 

		}
	}

	if ((itsIsDetailedSimulation) && itsNumDetailedSimulation == NumSimu+1) 
	{
		unsetVectorFlg("RISKY_VALUE"); 			
		adoptVectorValue("RISKY_VALUE", (ARM_Vector*)itsRiskyValue.Clone());

		unsetVectorFlg("TARGET"); 			
		adoptVectorValue("TARGET", (ARM_Vector*)itsTargetValue.Clone());
		
		unsetVectorFlg("MTM"); 			
		adoptVectorValue("MTM", (ARM_Vector*)itsMTM.Clone());
		
		unsetVectorFlg("CASH"); 			
		adoptVectorValue("CASH", (ARM_Vector*)itsCash.Clone());
		
		unsetVectorFlg("PNL"); 			
		adoptVectorValue("PNL", (ARM_Vector*)itsPNL.Clone());
		
		unsetVectorFlg("NOTIONAL"); 			
		adoptVectorValue("NOTIONAL", (ARM_Vector*)itsNotional.Clone());

		unsetVectorFlg("IMPLIED"); 			
		adoptVectorValue("IMPLIED", (ARM_Vector*)itsImpliedSpread.Clone());

		unsetVectorFlg("INSTRATE"); 			
		adoptVectorValue("INSTRATE", (ARM_Vector*)itsInstRate.Clone());

		unsetVectorFlg("DEFAULT"); 			
		adoptVectorValue("DEFAULT", (ARM_Vector*)itsDefaults.Clone());

		unsetVectorFlg("FORWARD"); 			
		adoptVectorValue("FORWARD", (ARM_Vector*)itsFwdRate.Clone());

		unsetVectorFlg("XT"); 			
		adoptVectorValue("XT", (ARM_Vector*)itsXt.Clone());

	}


	if ((!itsIsCashIn) && (!itsIsCashOut)) 
	{
		itsNbNoInNoOut++; 
		itsAverageFinalVal += V ;
	}
	
	itsSimulPNL.Elt(NumSimu) = PNL; 
	itsValo.Elt(NumSimu) = V; 
	itsCountRebalance += CountRebalance; 

}

// Computes the global outputs of the simulations 
void ICM_Pricer_CPDO::ComputeGlobalOutputs()
{
	itsCashInPerc = itsNbCashIn / itsNbSimulations; 
	itsCashOutPerc = itsNbCashOut / itsNbSimulations; 
	if (itsNbNoInNoOut) 
		itsAverageFinalVal /= itsNbNoInNoOut;
	itsCountRebalance /= itsNbSimulations; 
	
	if (itsNbCashIn)
		itsAverageCashInDay /= itsNbCashIn; 

	if (itsNbCashOut)
		itsAverageCashOutDay /= itsNbCashOut; 

	setValue("CASHIN_PERC",itsCashInPerc);
	setValue("CASHOUT_PERC",itsCashOutPerc);
	setValue("AVRG_FIN_VAL",itsAverageFinalVal);
	setValue("COUNT_REB",itsCountRebalance);
	setValue("CASHIN_DAY",itsAverageCashInDay);
	setValue("CASHOUT_DAY",itsAverageCashOutDay);

}

int ICM_Pricer_CPDO::GetSizeDFCurve(const double& NumPath)
{
	const ARM_Vector& cSimulationSchedule = getVectorValue("SIMULATION_SCHEDULE");

	double NbJours = cSimulationSchedule.Elt(NumPath);
	
	int SizeCurve0 = itsNoRiskyCouponDate.GetSize();
	int SizeCurve = 1; 

	for (int i = 1; i < SizeCurve0 ; i++ )
	{
		double Matu = itsNoRiskyCouponDate.Elt(i);
		if (((Matu - NbJours) > DB_TOL)) 
			SizeCurve ++;
	}

	return SizeCurve; 

}



// Generates the Discount Factor Curve for path NumPath of the simulation 
// This function is called at each path of the simulation 
void ICM_Pricer_CPDO::ComputeDiscountFactorCurve(const double& NumPath, 
												 ARM_Vector& DateVect, 
												 ARM_Vector& DFVect)
{

const ARM_Vector& cSimulationSchedule = getVectorValue("SIMULATION_SCHEDULE");

double NbJours = cSimulationSchedule.Elt(NumPath);
double Beta ; 
double Theta ; 
double Matu ; 
int i,j ; 
int SizeCurve0 = itsNoRiskyCouponDate.GetSize();

	DateVect.Elt(0) = 0;
	DFVect.Elt(0) = 1;

	j=1;
	for (i = 1; i <SizeCurve0 ; i++ )
	{

		Matu = itsNoRiskyCouponDate.Elt(i);

		if (((Matu - NbJours) > 0)) 
		{
			Beta = (1 - exp(- itsa * (Matu - NbJours) / 365.)) / itsa ; 

			Theta =  ( (exp(-2 * itsa * (Matu - NbJours)/ 365.) - exp(- 2 * itsa * Matu/ 365.)) - ( 1 - exp ( -2 * itsa * NbJours/ 365.)) ) / (2*itsa) 
				- 2 * ( (exp(-itsa * (Matu - NbJours)/ 365.) - exp(-itsa * Matu/ 365.)) - (1 - exp (-itsa * NbJours/ 365.)) ) / itsa ;
			
			
			DateVect.Elt(j) = Matu - NbJours; 
			DFVect.Elt(j) = itsDF0.Elt(i) / itsDF.Elt(NumPath) * (exp( - Beta * itsXt.Elt(NumPath))) * exp (- itsIRVol * itsIRVol / (2 * itsa * itsa) * Theta); 
			j++ ;
		}

	}


}


double ICM_Pricer_CPDO::ComputeTarget(const int& FirstPeriodLen,
									  const bool& IsCpnDate,
									  const double& V0, 
									  const double& Target, 
									  const double& RunningFees,
									  const ARM_Vector& VecDate,
									  const ARM_Vector& VecDF)
{

	double TargetPV; 
	int SizeCurve = VecDate.GetSize();

	if (SizeCurve > 1) 
	{
		if (IsCpnDate)
			TargetPV = VecDF.Elt(1) * VecDate.Elt(1)  ; 
		else
			TargetPV = VecDF.Elt(1) * FirstPeriodLen ;
	}
	else
		TargetPV = 0.; 

	for (int i = 2 ; i < SizeCurve ; i++)
	{
		TargetPV += VecDF.Elt(i) * (VecDate.Elt(i) - VecDate.Elt(i-1)) ; 
	}

	TargetPV *= V0 ;
	TargetPV *=  (Target + RunningFees); 
	TargetPV /=  365.; 

	TargetPV += V0 ;  
	return TargetPV; 

}

double ICM_Pricer_CPDO::SetEURIB3M(const ARM_Vector& VecDate,
								const ARM_Vector& VecDF)
{
	double NonRiskyCpn ;
	if (VecDF.GetSize()>1)	
		NonRiskyCpn = -365. / VecDate.Elt(1) * log(VecDF.Elt(1)); 
	
	return NonRiskyCpn; 

}


// Have to be defined - Otherwise class is abstract 

double ICM_Pricer_CPDO::Accrued()
{
	ICMTHROW(ERR_INVALID_ARGUMENT,"Not implemented"); 
}
double ICM_Pricer_CPDO::FeeLegPV () 
{
	ICMTHROW(ERR_INVALID_ARGUMENT,"Not implemented"); 
}

double ICM_Pricer_CPDO::DefLegPV () 
{
	ICMTHROW(ERR_INVALID_ARGUMENT,"Not implemented"); 
}

double ICM_Pricer_CPDO::ComputeDuration(void) 
{
	ICMTHROW(ERR_INVALID_ARGUMENT,"Not implemented"); 
}

double ICM_Pricer_CPDO::ComputeSpread(const double& MtM ) 
{
	ICMTHROW(ERR_INVALID_ARGUMENT,"Not implemented"); 
}

double ICM_Pricer_CPDO::ComputeImpliedVol(const double& Price)
{
	ICMTHROW(ERR_INVALID_ARGUMENT,"Not implemented"); 
}

double ICM_Pricer_CPDO::Compute_Fwd_Spread(const ARM_Date &Mty, const ARM_Date &CDS_ExpiryDate, double& dur) 
{
	ICMTHROW(ERR_INVALID_ARGUMENT,"Not implemented"); 
}


// 	virtual 
double ICM_Pricer_CPDO::ComputeSensitivity(const qSENSITIVITY_TYPE& typesensi, 
		const std::string& plot,
		const std::string& label,
		double epsilon , double epsilonGamma) 
{
	//return ICM_Pricer::ComputeSensitivity(typesensi,plot,label,epsilon, epsilonGamma); 
	ICMTHROW(ERR_INVALID_ARGUMENT,"Not implemented"); 

}

void ICM_Pricer_CPDO::DoPriceVector(const std::string& measure)
{

	if ((measure == "SIMULATION_SCHEDULE") || 
			(measure == "SCHED_RISKY") || 
			(measure == "SCHED_ISRISKY") || 
			(measure == "SCHED_ROLL") || 
			(measure == "SCHED_ISROLL") || 
			(measure == "SCHED_NORISKY") ||
			(measure == "SCHED_ISNORISKY"))
	{
		CreateSimulationSchedule();
	}	
	else if ((measure == "RISKY_VALUE") || 
			(measure == "TARGET") || 
			(measure == "MTM") || 
			(measure == "CASH") || 
			(measure == "PNL") || 
			(measure == "NOTIONAL") || 
			(measure == "IMPLIED") || 
			(measure == "INSTRATE") || 
			(measure == "DEFAULT") || 
			(measure == "FORWARD") || 
			(measure == "XT")) 
	{
		Price("CASHIN_PERC");
	}	
/*	else if (measure == "SCHED_ROLL")
	{
		CreateSimulationSchedule();
		adoptVectorValue(measure, (ARM_Vector*)itsRollDate.Clone());
	}	
	else if (measure == "SCHED_NORISKY")
	{
		CreateSimulationSchedule();
		adoptVectorValue(measure,(ARM_Vector*)itsNoRiskyCouponDate.Clone());
	}	*/
	else
		ICMTHROW(ERR_INVALID_ARGUMENT,"Type of measure does not exist"); 

}

double ICM_Pricer_CPDO::DoPrice(const std::string& measure)
{
	PriceVector("SIMULATION_SCHEDULE");

	if (measure == "CASHIN_PERC")
	{
		Simulate() ;
		return getValue("CASHIN_PERC") ;
	}
	else if (measure == "CASHOUT_PERC")
	{
		Simulate() ;
		return getValue("CASHOUT_PERC") ;
	}		
	else if (measure == "COUNT_REB")
	{
		Simulate() ;
		return getValue("COUNT_REB") ;
	}		
	else if (measure == "CASHIN_DAY")
	{
		Simulate() ;
		return getValue("CASHIN_DAY") ;
	}		
	else if (measure == "CASHOUT_DAY")
	{
		Simulate() ;
		return getValue("CASHOUT_DAY") ;
	}		
	else if (measure == "AVRG_FIN_VAL")
	{
		Simulate() ;
		return getValue("AVRG_FIN_VAL") ;
	}		
	else
		ICMTHROW(ERR_INVALID_ARGUMENT,"Type of measure does not exist"); 


}

