#include <CCString.h>

#include "ARM_result.h"
#include "ARM_interglob.h"

double ARMDATE2XLDATE (const CCString& armdate)
{
	int y, m, d;

	sscanf ((const char*)armdate, "%2d.%2d.%4d", &d, &m, &y);
	
	return (DAT_struct_to_ssdate (y, m, d));
}

CCString XLDATE2ARMDATE (double xldate)
{
	int y, m, d;
	char buf[11];
	CCString tmp;

	long long_xldate = (long)xldate;

	DAT_ssdate_to_struct ((double)long_xldate, &y, &m, &d);

	if(d < 10)
	{
		sprintf (buf, "0%1d.", d);
	}
	else
	{
		sprintf (buf, "%2d.", d);
	}
	tmp.Set (buf);
	
	if(m < 10)
	{
		sprintf (buf, "0%1d.", m);
	}
	else
	{
		sprintf (buf, "%2d.", m);
	}
	tmp += CCString (buf);

	sprintf (buf, "%4d", y);

	tmp += CCString (buf);

	return (tmp);
}

long ARM_ConvPriceYield (const CCString& aPriceYield, ARM_result& result)
{
	CCString tmp = aPriceYield;
	tmp.toUpper ();

	if(tmp == "PRICE")
	{
		return K_PRICE;
	}
	if(tmp == "P")
	{
		return K_PRICE;
	}
	if(tmp == "0")
	{
		return K_PRICE;
	}
	if(tmp == "YIELD")
	{
		return K_YIELD;
	}
	if(tmp == "Y")
	{
		return K_YIELD;
	}
	if(tmp == "1")
	{
		return K_YIELD;
	}
    
	result.setRetCode (ARM_DEFAULT_ERR);
	result.setMsg ("ARM_ERR: Invalid Parameter - Valid are Yield, Price, Y, P, 0, 1");

	return ARM_DEFAULT_ERR;
}



long ARM_ConvIrIndName (const CCString& aIrIndName, ARM_result& result)
{
	CCString tmp = aIrIndName;
	tmp.toUpper ();

	if(tmp == "PIBOR1M")
	{
		return K_PIBOR1M;
	}
	if(tmp == "PIBOR3M")
	{
		return K_PIBOR3M;
	}
	if(tmp == "PIBOR6M")
	{
		return K_PIBOR6M;
	}
	if(tmp == "PIBOR1Y")
	{
		return K_PIBOR1Y;
	}
	if(tmp == "LIBOR1M")
	{
		return K_LIBOR1M;
	}
	if(tmp == "LIBOR3M")
	{
		return K_LIBOR3M;
	}
	if(tmp == "LIBOR6M")
	{
		return K_LIBOR6M;
	}
	if(tmp == "LIBOR1Y")
	{
		return K_LIBOR1Y;
	}
	if(tmp == "EURIBOR1M")
	{
		return K_EURIBOR1M;
	}
	if(tmp == "EURIBOR3M")
	{
		return K_EURIBOR3M;
	}
	if(tmp == "EURIBOR6M")
	{
		return K_EURIBOR6M;
	}
	if(tmp == "EURIBOR1Y")
	{
		return K_EURIBOR1Y;
	}
	if(tmp == "EURIBOR12M")
	{
		return K_EURIBOR1Y;
	}
	
	result.setRetCode (ARM_DEFAULT_ERR);
	result.setMsg ("ARM_ERR: Invalid Type Name - Valid are Euribor1m, Euribor3m, Euribor6m, Euribor1y, Libor1m, Libor3m, Libor6m, Libor1y");

	return ARM_DEFAULT_ERR;
}



long ARM_ConvIrIndNameToFreq (const CCString& aIrIndName, ARM_result& result)
{
	CCString tmp = aIrIndName;
	tmp.toUpper ();

	if(tmp == "PIBOR1M")
	{
		return K_MONTHLY;
	}
	if(tmp == "PIBOR3M")
	{
		return K_QUARTERLY;
	}
	if(tmp == "PIBOR6M")
	{
		return K_SEMIANNUAL;
	}
	if(tmp == "PIBOR1Y")
	{
		return K_ANNUAL;
	}
	if(tmp == "LIBOR1M")
	{
		return K_MONTHLY;
	}
	if(tmp == "LIBOR3M")
	{
		return K_QUARTERLY;
	}
	if(tmp == "LIBOR6M")
	{
		return K_SEMIANNUAL;
	}
	if(tmp == "LIBOR1Y")
	{
		return K_ANNUAL;
	}
	if(tmp == "EURIBOR1M")
	{
		return K_MONTHLY;
	}
	if(tmp == "EURIBOR3M")
	{
		return K_QUARTERLY;
	}
	if(tmp == "EURIBOR6M")
	{
		return K_SEMIANNUAL;
	}
	if(tmp == "EURIBOR1Y")
	{
		return K_ANNUAL;
	}
	if(tmp == "EURIBOR12M")
	{
		return K_ANNUAL;
	}
	
	result.setRetCode (ARM_DEFAULT_ERR);
	result.setMsg ("ARM_ERR: Invalid Type Name - Valid are Euribor1m, Euribor3m, Euribor6m, Euribor1y, Libor1m, Libor3m, Libor6m, Libor1y");

	return ARM_DEFAULT_ERR;
}



long ARM_ConvTmIxName (const CCString& aTmIxName, ARM_result& result)
{ 
	CCString tmp = aTmIxName;
	tmp.toUpper ();

    if(tmp == "T4M")
	{
         return 26;
	}
    if(tmp == "TAM")
	{
         return 28;
	}
    if(tmp == "TMP")
	{
         return 30;
	}
    if(tmp == "T4M_FIXED")
	{
         return 27;
	}

	result.setRetCode (ARM_DEFAULT_ERR);
	result.setMsg ("ARM_ERR: Invalid Type Name - Valid are T4M, TAM, TMP, T4M_FIXED");

	return ARM_DEFAULT_ERR;
}



long ARM_ConvYieldOrVol (const CCString& yieldOrVol, ARM_result& result)
{ 
	CCString tmp = yieldOrVol;
	tmp.toUpper ();

    if(tmp == "Y")
	{
         return 0;
	}
    if(tmp == "V")
	{
         return 1;
	}
    
	result.setRetCode (ARM_DEFAULT_ERR);
	result.setMsg ("ARM_ERR: Invalid Type Name - Valid are Y, V");

	return ARM_DEFAULT_ERR;
}



long ARM_ConvCalcMod (const CCString& calcMod, ARM_result& result)
{ 
	CCString tmp = calcMod;
	tmp.toUpper ();

    if(tmp == "NOR")
	{
         return 0;
	}
    if(tmp == "LOGNOR")
	{
         return 1;
	}
	if(tmp == "SQR")
	{
         return 2;
	}
    
	result.setRetCode (ARM_DEFAULT_ERR);
	result.setMsg ("ARM_ERR: Invalid Type Name - Valid are NOR, LOGNOR, SQR");

	return ARM_DEFAULT_ERR;
}



long ARM_ConvFrequency (const CCString& aFrequency, ARM_result& result)
{ 
	CCString tmp = aFrequency;
	tmp.toUpper ();

	if(tmp == "-1")
	{
		return K_DEF_FREQ;
	}
    if(tmp == "A")
	{
        return K_ANNUAL;
	}
    if(tmp == "S")
	{
        return K_SEMIANNUAL;
	}
    if(tmp == "Q")
	{
        return K_QUARTERLY;
	}
    if(tmp == "B")
	{
        return K_BIMONTHLY;
	}
    if(tmp == "M")
	{
        return K_MONTHLY;
	}
    if(tmp == "W")
	{
		return K_WEEKLY;
	}
    if(tmp == "D")
	{
		return K_DAILY;
	}
    if(tmp == "Z")
	{
        return K_ZEROCOUPON;
	}
    if(tmp == "1")
	{
        return K_ANNUAL;
	}
    if(tmp == "2")
	{
		return K_SEMIANNUAL;
	}
    if(tmp == "4")
	{
        return K_QUARTERLY;
	}
    if(tmp == "6")
	{
        return K_BIMONTHLY;
	}
    if(tmp == "12")
	{
		return K_MONTHLY;
	}
    if(tmp == "52")
	{
        return K_WEEKLY;
	}
    if(tmp == "365")
	{
        return K_DAILY;
	}
    if(tmp == "0")
	{
         return K_ZEROCOUPON;
	}

	result.setRetCode (ARM_DEFAULT_ERR);
	result.setMsg ("ARM_ERR: Invalid Type Name - Valid are A, S, Q, B(imonth), M, D, Z");

	return ARM_DEFAULT_ERR;
}



long ARM_ConvDecompFrequency (const CCString& aFrequency, ARM_result& result)
{
	CCString tmp = aFrequency;
	tmp.toUpper ();

    if(tmp == "C")
	{
		return K_COMP_CONT;
	}
    if(tmp == "P")
	{
        return K_COMP_PROP;
	}
    if(tmp == "A")
	{
        return K_COMP_ANNUAL;
	}
    if(tmp == "S")
	{
        return K_COMP_SEMIANNUAL;
	}
    if(tmp == "Q")
	{
        return K_COMP_QUARTERLY;
	}
    if(tmp == "M")
	{
        return K_COMP_MONTHLY;
	}
    if(tmp == "B")
	{
        return K_COMP_BIMONTHLY;
	}
    if(tmp == "D360")
	{
        return K_COMP_DAILY_360;
	}
    if(tmp == "D365")
	{
        return K_COMP_DAILY_365;
	}
    if(tmp == "0")
	{
        return K_COMP_CONT;
	}
    if(tmp == "-1")
	{
        return K_COMP_PROP;
	}
    if(tmp == "1")
	{
        return K_COMP_ANNUAL;
	}
    if(tmp == "2")
	{
        return K_COMP_SEMIANNUAL;
	}
    if(tmp == "4")
	{
        return K_COMP_QUARTERLY;
	}
    if(tmp == "12")
	{
        return K_COMP_MONTHLY;
	}
    if(tmp == "6")
	{
        return K_COMP_BIMONTHLY;
	}
    if(tmp == "360")
	{
        return K_COMP_DAILY_360;
	}
    if(tmp == "365")
	{
        return K_COMP_DAILY_365;
	}

	result.setRetCode (ARM_DEFAULT_ERR);
	result.setMsg ("ARM_ERR: Invalid Type Name - Valid are A, S, Q, B(imonth), M, D, Z");

	return ARM_DEFAULT_ERR;
}


 
long ARM_ConvIrType (const CCString& aIrType, ARM_result& result)
{
	CCString tmp = aIrType;
	tmp.toUpper ();

    if(tmp == "FIXED")
	{
        return K_FIXED;
	}
    if(tmp == "LIBOR1M")
	{
        return K_LIBOR1M;
	}
    if(tmp == "LIBOR2M")
	{
        return K_LIBOR2M;
	}
    if(tmp == "LIBOR3M")
	{
        return K_LIBOR3M;
	}
    if(tmp == "LIBOR6M")
	{
        return K_LIBOR6M;
	}
    if(tmp == "LIBOR1Y")
	{
        return K_LIBOR1Y;
	}
    if(tmp == "LIBOR12M")
	{
		return K_LIBOR1Y;
	}
    if(tmp == "PIBOR1M")
	{
        return K_PIBOR1M;
	}
    if(tmp == "PIBOR2M")
	{
        return K_PIBOR2M;
	}
    if(tmp == "PIBOR3M")
	{
        return K_PIBOR3M;
	}
    if(tmp == "PIBOR6M")
	{
        return K_PIBOR6M;
	}
    if(tmp == "PIBOR1Y")
	{
        return K_PIBOR1Y;
	}
    if(tmp == "EURIBOR1M")
	{
        return K_EURIBOR1M;
	}
    if(tmp == "EURIBOR2M")
	{
        return K_EURIBOR2M;
	}
	if(tmp == "EURIBOR3M")
	{
        return K_EURIBOR3M;
	}
    if(tmp == "EURIBOR6M")
	{
        return K_EURIBOR6M;
	}
    if(tmp == "EURIBOR1Y")
	{
        return K_EURIBOR1Y;
	}
    if(tmp == "EURIBOR12M")
	{
        return K_EURIBOR1Y;
	}
    if(tmp == "CMT1")
	{
        return K_CMT1;
	}
    if(tmp == "CMT2")
	{
        return K_CMT2;
	}
    if(tmp == "CMT5")
	{
        return K_CMT5;
	}
    if(tmp == "CMT10")
	{
        return K_CMT10;
	}
    if(tmp == "CMT15")
	{
        return K_CMT15;
	}
    if(tmp == "CMT20")
	{
        return K_CMT20;
	}
    if(tmp == "CMT30")
	{
        return K_CMT30;
	}
    if(tmp == "CMS1")
	{
        return K_CMS1;
	}
    if(tmp == "CMS2")
	{
        return K_CMS2;
	}
    if(tmp == "CMS3")
	{
        return K_CMS3;
	}
    if(tmp == "CMS4")
	{
        return K_CMS4;
	}
    if(tmp == "CMS5")
	{
        return K_CMS5;
	}
    if(tmp == "CMS6")
	{
        return K_CMS6;
	}   
    if(tmp == "CMS7")
	{
        return K_CMS7;
	}   
    if(tmp == "CMS8")
	{
        return K_CMS8;
	}   
    if(tmp == "CMS9")
	{
        return K_CMS9;
	}   
    if(tmp == "CMS10")
	{
        return K_CMS10;
	}
    if(tmp == "CMS10")
	{
        return K_CMS10;
	}
    if(tmp == "CMS11")
	{
        return K_CMS11;
	}
    if(tmp == "CMS12")
	{
        return K_CMS12;
	}
    if(tmp == "CMS13")
	{
        return K_CMS13;
	}
    if(tmp == "CMS14")
	{
        return K_CMS14;
	}
    if(tmp == "CMS15")
	{
        return K_CMS15;
	}
    if(tmp == "CMS16")
	{
        return K_CMS16;
	}
    if(tmp == "CMS17")
	{
        return K_CMS17;
	}
    if(tmp == "CMS18")
	{
        return K_CMS18;
	}
    if(tmp == "CMS19")
	{
        return K_CMS19;
	}
    if(tmp == "CMS20")
	{
        return K_CMS20;
	}
    if(tmp == "CMS21")
	{
        return K_CMS21;
	}
   if(tmp == "CMS22")
	{
        return K_CMS22;
	}
   if(tmp == "CMS23")
	{
        return K_CMS23;
	}
   if(tmp == "CMS24")
	{
        return K_CMS24;
	}
   if(tmp == "CMS25")
	{
        return K_CMS25;
	}
   if(tmp == "CMS26")
	{
        return K_CMS26;
	}
   if(tmp == "CMS27")
	{
        return K_CMS27;
	}
   if(tmp == "CMS28")
	{
        return K_CMS28;
	}
   if(tmp == "CMS29")
	{
        return K_CMS29;
	}
    if(tmp == "CMS30")
	{
        return K_CMS30;
	}
    if(tmp == "TEC5")
	{
        return K_TEC5;
	}
    if(tmp == "TEC10")
	{
        return K_TEC10;
	}
    if(tmp == "T4M")
	{
        return K_T4M;
	}
    if(tmp == "T4M_FIXED")
	{
        return K_T4M_FIXED;
	}
    if(tmp == "TAM")
	{
        return K_TAM;
	}
    if(tmp == "TAG")
	{
        return K_TAG;
	}
    if(tmp == "TMP")
	{
        return K_TMP;
	}
	if(tmp == "STD")
	{
		return K_STD;
	}
 
	result.setRetCode (ARM_DEFAULT_ERR);
	result.setMsg ("ARM_ERR: Invalid Type index - Valid are Libor3m Pibor6m CMS10 etc ...");

	return ARM_DEFAULT_ERR;
}



long ARM_ConvDayCount (const CCString& aDayCount, ARM_result& result)
{
	CCString tmp = aDayCount;
	tmp.toUpper ();

	if(tmp == "ACTUAL")
	{
        return KACTUAL_ACTUAL;
	}
    if(tmp == "A365")
	{
        return KACTUAL_365;
	}
    if(tmp == "A360")
	{
        return KACTUAL_360;
	}
    if(tmp == "30/360")
	{
        return K30_360;
	}
    if(tmp == "ACTREAL")
	{
        return KACTUAL_REAL;
	}
    if(tmp == "1")
	{
        return KACTUAL_ACTUAL;
	}
    if(tmp == "2")
	{
        return KACTUAL_365;
	}
    if(tmp == "3")
	{
        return KACTUAL_360;
	}
    if(tmp == "4")
	{
        return K30_360;
	}
    if(tmp == "5")
	{
        return KACTUAL_REAL;
	}
	if(tmp == "6")
	{
        return KACTUAL_FEB29;
	}
	if(tmp == "7")
	{
        return KACTUAL_ISMA;
	}
	if(tmp == "-1")
	{
		return -1;
	}

	result.setRetCode (ARM_DEFAULT_ERR);
	result.setMsg ("ARM_ERR: Invalid DayCount - Valid are ACTUAL: A365, A360, 30/360, ACTREAL");

	return ARM_DEFAULT_ERR;
}



long ARM_ConvCompMeth (const CCString& aCompMeth, ARM_result& result)
{
	CCString tmp = aCompMeth;
	tmp.toUpper ();

	if(tmp == "C")
	{
        return K_COMP_CONT;
	}
    if(tmp == "P")
	{
        return K_COMP_PROP;
	}
    if(tmp == "A")
	{
        return K_COMP_ANNUAL;
	}
    if(tmp == "S")
	{
        return K_COMP_SEMIANNUAL;
	}
    if(tmp == "Q")
	{
        return K_COMP_QUARTERLY;
	}
    if(tmp == "B")
	{
        return K_COMP_BIMONTHLY;
	}
    if(tmp == "M")
	{
        return K_COMP_MONTHLY;
	}
    if(tmp == "D360")
	{
        return K_COMP_DAILY_360;
	}
    if(tmp == "D365")
	{
        return K_COMP_DAILY_365;
	}
    if(tmp == "0")
	{
        return K_COMP_CONT;
	}
    if(tmp == "-1")
	{
        return K_COMP_PROP;
	}
    if(tmp == "1")
	{
        return K_COMP_ANNUAL;
	}
    if(tmp == "2")
	{
        return K_COMP_SEMIANNUAL;
	}
    if(tmp == "4")
	{
        return K_COMP_QUARTERLY;
	}

	result.setRetCode (ARM_DEFAULT_ERR);
	result.setMsg ("ARM_ERR: Invalid Compounding Method - Valid are: 'C'ontinuous, 'P'roportional (Lineaire), 'A', 'S', 'Q', 'B', 'M', 'D360', 'D365'");

	return ARM_DEFAULT_ERR;
}



long ARM_ConvForwardYieldMethod (const CCString& aForwardYieldMeth, ARM_result& result)
{
	CCString tmp = aForwardYieldMeth;
	tmp.toUpper ();

	if(tmp == "MON")
	{
        return -1;
	}
	if(tmp == "-1")
	{
		return -1;
	}
    if(tmp == "CONT")
	{
        return 0;
	}
	if(tmp[0] == '0')
	{
		return 0;
	}
	if(tmp == "ACTU")
	{
        return 1;
	}
	if(tmp[0] == '1')
	{
		return 1;
	}
	
	result.setRetCode (ARM_DEFAULT_ERR);
	result.setMsg ("ARM_ERR: Invalid Forward Yield Method - Valid are: 'MON', '-1', 'CONT', '0', 'ACTU', '1'");

	return ARM_DEFAULT_ERR;
}



long ARM_ConvCallOrPut (const CCString& aCorP, ARM_result& result)
{
	CCString tmp = aCorP;
	tmp.toUpper ();

    if(tmp == "CALL")
	{
        return K_CALL;
	}
    if(tmp == "PUT")
	{
        return K_PUT;
	}
    if(tmp == "C")
	{
        return K_CALL;
	}
    if(tmp == "P")
    {
		return K_PUT;
	}
    if(tmp == "1")
	{
        return K_CALL;
	}
    if(tmp == "-1")
	{
        return K_PUT;
	}
 
	result.setRetCode (ARM_DEFAULT_ERR);
	result.setMsg ("ARM_ERR: Invalid Option Type - Valid are: (Call, C, 1) Or (Put, P, -1)");

	return ARM_DEFAULT_ERR;
}



long ARM_ConvParam (const CCString& param, ARM_result& result)
{
	CCString tmp = param;
	tmp.toUpper ();

    if(tmp == "PRICE")
	{
        return 0;
	}
    if(tmp == "DELTA")
	{
        return 1;
	}
    
	result.setRetCode (ARM_DEFAULT_ERR);
	result.setMsg ("ARM_ERR: Invalid Param - Valid are: PRICE, DELTA");

	return ARM_DEFAULT_ERR;
}



long ARM_ConvMktType (const CCString& aMkt, ARM_result& result)
{
	CCString tmp = aMkt;
	tmp.toUpper ();

    if(tmp[0] == 'F')
	{
        return K_FUT;
	}
    if(tmp[0] == 'M')
	{
        return K_MM;
	}
    if(tmp[0] == 'S')
	{
        return K_SWAP;
	}
    if(tmp[0] == 'B')
	{
        return K_BOND;
	}
    if(tmp[0] == '0')
	{
        return K_FUT;
	}
    if(tmp[0] == '1')
	{
        return K_MM;
	}
    if(tmp[0] == '2')
	{
        return K_SWAP;
	}
    if(tmp[0] == '3')
	{
        return K_BOND;
	}

	result.setRetCode (ARM_DEFAULT_ERR);
	result.setMsg ("ARM_ERR: Invalid Market - Valid are: (FUT,0) - (MM,1) - (SWAP,2) - (BOND,3)");

	return ARM_DEFAULT_ERR;
}



long ARM_ConvCvMethod (const CCString& aCvMeth, ARM_result& result)
{
	CCString tmp = aCvMeth;
	tmp.toUpper ();
	
	if(tmp[0] == 'P')
	{
		return K_PAR;
	}
    if(tmp[0] == 'R')
	{
        return K_RAW;
	}
    if(tmp[0] == '0')
	{
        return K_PAR;
	}
    if(tmp[0] == '1')
	{
       return K_RAW;
	}

	result.setRetCode (ARM_DEFAULT_ERR);
	result.setMsg ("ARM_ERR: Invalid Cv Method - Valid are: (PAR,0) - (RAW,1)");

	return ARM_DEFAULT_ERR;
}


long ARM_ConvInterpMethod (const CCString& aIntMeth, ARM_result& result)
{
	CCString tmp = aIntMeth;
	tmp.toUpper ();

    if(tmp[0] == 'C')
	{
		return K_CONTINUOUS;
	}
    if(tmp[0] == 'L')
	{
        return K_LINEAR;
	}
    if(tmp[0] == '0')
	{
        return K_CONTINUOUS;
	}
    if(tmp[0] == '1')
	{
        return K_LINEAR;
	}
	
	result.setRetCode (ARM_DEFAULT_ERR);
	result.setMsg ("ARM_ERR: Invalid Interpol Method - Valid are: (CONT,0) - (LIN,1)");

	return ARM_DEFAULT_ERR;
}



long ARM_ConvCapOrFloor (const CCString& aCorF, ARM_result& result)
{ 
	CCString tmp = aCorF;
	tmp.toUpper ();

    if(tmp == "CAP")
	{
        return K_CAP;
	}
    if(tmp == "FLOOR")
	{
        return K_FLOOR;
	}
    if(tmp == "C")
	{
        return K_CAP;
	}
    if(tmp == "F")
	{
        return K_FLOOR;
	}
    if(tmp == "1")
	{
        return K_CAP;
	}
    if(tmp == "-1")
	{
        return K_FLOOR;
	}

	result.setRetCode (ARM_DEFAULT_ERR);
	result.setMsg ("ARM_ERR: Invalid CapOrFloor - Valid are : Cap, C, 1 or Floor, F, -1");

	return ARM_DEFAULT_ERR;
}



long ARM_ConvExerciseType (const CCString& aOptType, ARM_result& result)
{
	CCString tmp = aOptType;
	tmp.toUpper ();

    if(tmp == "E")
	{
		return K_EUROPEAN;
	}
    if(tmp == "A")
	{
        return K_AMERICAN;
	}
    if(tmp == "B")
	{
        return K_BERMUDAN;
	}
 
	result.setRetCode (ARM_DEFAULT_ERR);
	result.setMsg ("ARM_ERR: Invalid Exercise Type - Valid are: A{merican}, E{uropean}, B{ermudan}");

	return ARM_DEFAULT_ERR;
}



long ARM_ConvRecOrPay (const CCString& aRorP, ARM_result& result)
{
	CCString tmp = aRorP;
	tmp.toUpper ();

    if(tmp == "R")
	{
        return K_RCV;
	}
    if(tmp == "P")
	{
		return K_PAY;
	}
    if(tmp == "REC")
	{
		return K_RCV;
	}
    if(tmp == "PAY")
	{
		return K_PAY;
	}
    if(tmp == "1")
	{
		return K_RCV;
	}
    if(tmp == "-1")
	{
		return K_PAY;
	}

	result.setRetCode (ARM_DEFAULT_ERR);
	result.setMsg ("ARM_ERR: Invalid RcvOrPay - Valid are REC, R, r, 1 or PAY, P, p, -1");

	return ARM_DEFAULT_ERR;
}



long ARM_ConvFwdRule (const CCString& aFwdRule, ARM_result& result)
{
	CCString tmp = aFwdRule;
	tmp.toUpper ();

    if(tmp == "F")
	{
        return K_FOLLOWING;
	}
    if(tmp == "MF")
	{
        return K_MOD_FOLLOWING;
	}
    if(tmp == "P")
	{
        return K_PREVIOUS;
	}
    if(tmp == "MP")
	{
        return K_MOD_PREVIOUS;
	}
    if(tmp == "1")
	{
        return K_FOLLOWING;
	}
    if(tmp == "2")
	{
        return K_MOD_FOLLOWING;
	}
    if(tmp == "-1")
	{
        return K_PREVIOUS;
	}
    if(tmp == "-2")
	{
        return K_MOD_PREVIOUS;
	}

	result.setRetCode (ARM_DEFAULT_ERR);
	result.setMsg ("ARM_ERR: Invalid Fwd Rule - Valid are: F(ollowing), f, 1 or P(revious), p, -1");

	return ARM_DEFAULT_ERR;
}



long ARM_ConvRule (const CCString& rule, ARM_result& result)
{
	CCString tmp = rule;
	tmp.toUpper ();

    if(tmp == "F")
	{
        return 1;
	}
	if(tmp == "1")
	{
        return 1;
	}
    if(tmp == "P")
	{
        return -1;
	}
	if(tmp == "-1")
	{
        return -1;
	}
    if(tmp == "MF")
	{
        return 2;
	}
	if(tmp == "2")
	{
        return 2;
	}
    if(tmp == "MP")
	{
        return -2;
	}
	if(tmp == "-2")
	{
        return -2;
	}
    
	result.setRetCode (ARM_DEFAULT_ERR);
	result.setMsg ("ARM_ERR: Invalid Rule - Valid are: F, MF, P, MP");

	return ARM_DEFAULT_ERR;
}


 
long ARM_ConvPayResetRule (const CCString& aRule, ARM_result& result)
{
	CCString tmp = aRule;
	tmp.toUpper ();

    if(tmp == "ADV")
	{
        return K_ADVANCE;
	}
    if(tmp == "ADVANCE")
	{
        return K_ADVANCE;
	}
    if(tmp == "ARR")
	{
        return K_ARREARS;
	}
    if(tmp == "ARREARS")
	{
        return K_ARREARS;
	}
    if(tmp == "1")
	{
        return K_ADVANCE;
	}
    if(tmp == "-1")
	{
        return K_ARREARS;
	}

	result.setRetCode (ARM_DEFAULT_ERR);
	result.setMsg ("ARM_ERR: Invalid Payment or Reset Rule - Valid are  ADV, ADVANCE, 1 or ARR, ARREARS, -1");

	return ARM_DEFAULT_ERR;
}


 
long ARM_ConvIntRule (const CCString& aRule, ARM_result& result)
{
	CCString tmp = aRule;
	tmp.toUpper ();

    if(tmp == "ADJ")
	{
		return K_ADJUSTED;
	}
    if(tmp == "A")
	{
        return K_ADJUSTED;
	}
    if(tmp == "NOADJ")
	{
        return K_UNADJUSTED;
	}
    if(tmp == "UNADJ")
	{
        return K_UNADJUSTED;
	}
    if(tmp == "N")
	{
        return K_UNADJUSTED;
	}
    if(tmp == "U")
	{
        return K_UNADJUSTED;
	}
    if(tmp == "1")
	{
        return K_ADJUSTED;
	}
    if(tmp == "0")
	{
        return K_UNADJUSTED;
	}

	result.setRetCode (ARM_DEFAULT_ERR);
	result.setMsg ("ARM_ERR: Invalid Interest Rule - Valid are  ADJ, A, 1 or NOADJ, N, 0");

	return ARM_DEFAULT_ERR;
}



long ARM_ConvStubRule (const CCString& aRule, ARM_result& result)
{
	CCString tmp = aRule;
	tmp.toUpper ();

	if(tmp == "SS")
	{
        return K_SHORTSTART;
	}
    if(tmp == "LS")
	{
        return K_LONGSTART;
	}
    if(tmp == "SE")
	{
        return K_SHORTEND;
	}
    if(tmp == "LE")
	{
        return K_LONGEND;
	}
    if(tmp == "1")
	{
        return K_SHORTSTART;
	}
    if(tmp == "2")
	{
        return K_LONGSTART;
	}
    if(tmp == "3")
	{
        return K_SHORTEND;
	}
    if(tmp == "4")
	{
        return K_LONGEND;
	}
 
	result.setRetCode (ARM_DEFAULT_ERR);
	result.setMsg ("Invalid Stub Rule - Valid are  SS , LS, SE, LE or 1, 2, 3, 4");

	return ARM_DEFAULT_ERR;
}


 
long ARM_ConvObjectClass (const CCString& idUnderlying, ARM_result& result)
{
	CCString tmp = idUnderlying;
	tmp.toUpper ();

	if(tmp == "BOND")
	{
		return 1;
	}
	if(tmp == "PORT")
	{
		return 1;
	}

	result.setRetCode (ARM_DEFAULT_ERR);
	result.setMsg ("Invalid Object Class");

	return ARM_DEFAULT_ERR;
}



long ARM_ConvSvtyParam (const CCString& aParam, ARM_result& result)
{
	CCString tmp = aParam;
	tmp.toUpper ();

	if(tmp == "DELTA")
	{
		return 0;
	}
	if(tmp == "GAMMA")
	{
		return 1;
	}
	if(tmp == "RHO")
	{
		return 2;
	}
	if(tmp == "THETA")
	{
		return 3;
	}
	if(tmp == "VEGA")
	{
		return 4;
	}
	if(tmp == "KAPPA")
	{
		return 5;
	}
	if(tmp[0] == '0')
	{
		return 0;
	}
	if(tmp[0] == '1')
	{
		return 1;
	}
	if(tmp[0] == '2')
	{
		return 2;
	}
	if(tmp[0] == '3')
	{
		return 3;
	}
	if(tmp[0] == '4')
	{
		return 4;
	}
	if(tmp[0] == '5')
	{
		return 5;
	}

	result.setRetCode (ARM_DEFAULT_ERR);
	result.setMsg ("Invalid Parameter");

	return ARM_DEFAULT_ERR;
}



long ARM_ConvCurvSvtyParam (const CCString& aParam, ARM_result& result)
{
	CCString tmp = aParam;
	tmp.toUpper ();

	if(tmp == "TAUX")
	{
		return 10;
	}
	if(tmp == "VOL")
	{
		return 11;
	}
	long buf;
	if(sscanf (tmp, "%ld", &buf) == 1)
	{
		return (buf + 10);
	}

	result.setRetCode (ARM_DEFAULT_ERR);
	result.setMsg ("Invalid Parameter");

	return ARM_DEFAULT_ERR;
}



long ARM_ConvUpDown (const CCString& aUpDown, ARM_result& result)
{
	CCString tmp = aUpDown;
	tmp.toUpper ();

	if(tmp == "UP")
	{
		return 1;
	}
	if(tmp == "DOWN")
	{
		return -1;
	}

	long buf;
	if((sprintf (tmp, "%ld", &buf)) == 1)
	{
		return buf;
	}

	result.setRetCode (ARM_DEFAULT_ERR);
	result.setMsg ("Invalid Parameter - Valid are Up, Down, 1, -1");

	return ARM_DEFAULT_ERR;
}



long ARM_ConvInOut (const CCString& aInOut, ARM_result& result)
{
	CCString tmp = aInOut;
	tmp.toUpper ();

	if(tmp == "IN")
	{
		return 1;
	}
	if(tmp == "OUT")
	{
		return -1;
	}

	long buf;
	if((sprintf (tmp, "%ld", &buf)) == 1)
	{
		return buf;
	}

	result.setRetCode (ARM_DEFAULT_ERR);
	result.setMsg ("Invalid Parameter - Valid are In, 1, Out, -1");

	return ARM_DEFAULT_ERR;
}



long ARM_ConvVolType (const CCString& aVolType, ARM_result& result)
{
	CCString tmp = aVolType;
	tmp.toUpper ();

	if(tmp == "ATM")
	{
		return K_ATMF_VOL;
	}
	if(tmp[0] == 'A')
	{
		return K_ATMF_VOL;
	}
	if(tmp == "SMILE")
	{
		return K_SMILE_VOL;
	}
	if(tmp[0] == 'S')
	{
		return K_SMILE_VOL;
	}

	long buf;
	if((sprintf (tmp, "%ld", &buf)) == 1)
	{
		return buf;
	}

	result.setRetCode (ARM_DEFAULT_ERR);
	result.setMsg ("Invalid Parameter - Valid are ATM, Smile, A, S");

	return ARM_DEFAULT_ERR;
}

long ARM_ConvIrIndIdToFreqId (long IrIndId, ARM_result& result)
{
	if(IrIndId == K_PIBOR1M)
	{
		return K_MONTHLY;
	}
	if(IrIndId == K_PIBOR3M)
	{
		return K_QUARTERLY;
	}
	if(IrIndId == K_PIBOR6M)
	{
		return K_SEMIANNUAL;
	}
	if(IrIndId == K_PIBOR1Y)
	{
		return K_ANNUAL;
	}
	if(IrIndId == K_LIBOR1M)
	{
		return K_MONTHLY;
	}
	if(IrIndId == K_LIBOR3M)
	{
		return K_QUARTERLY;
	}
	if(IrIndId == K_LIBOR6M)
	{
		return K_SEMIANNUAL;
	}
	if(IrIndId == K_LIBOR1Y)
	{
		return K_ANNUAL;
	}
	if(IrIndId == K_EURIBOR1M)
	{
		return K_MONTHLY;
	}
	if(IrIndId == K_EURIBOR3M)
	{
		return K_QUARTERLY;
	}
	if(IrIndId == K_EURIBOR6M)
	{
		return K_SEMIANNUAL;
	}
	if(IrIndId == K_EURIBOR1Y)
	{
		return K_ANNUAL;
	}
	
	return ARM_DEFAULT_ERR;
}

long GetNumObjectId (const CCString& stringObjectId)
{
	if(!stringObjectId)
	{
		return ARM_KO;
	}

	if(stringObjectId == "ARM_ERR")
	{
		return ARM_KO;
	}

	if(stringObjectId == "ARM_FATAL_ERR")
	{
		return ARM_KO;
	}

	if(stringObjectId == "OBSOLETE")
	{
		return ARM_KO;
	}

	if(stringObjectId == "ERROR")
	{
		return ARM_KO;
	}

	CCString curClass = GetStringObjectClass (stringObjectId);

	if(!((curClass == ZERO_CURVE_LIN_CLASS) ||
		(curClass == ZERO_CURVE_SPL_CLASS) ||
		(curClass == ZERO_CURVE_VSK_CLASS) ||
		(curClass == ZERO_CURVE_FLAT_CLASS) ||
		(curClass == ZERO_CURVE_SPLCUB_CLASS) ||
		(curClass == ZERO_CURVE_CUBDIFF_CLASS) ||
		(curClass == ZERO_CURVE_INTSMO_CLASS) ||
        (curClass == VOL_CURVE_LIN_CLASS) ||
		(curClass == VOL_FLAT_CLASS) ||
        (curClass == BOND_CLASS) ||
        (curClass == CCY_CLASS) ||
        (curClass == FOREX_CLASS) ||
        (curClass == PIBFUT_CLASS) ||
        (curClass == IRFUT_CLASS) ||
        (curClass == NNNFUT_CLASS) ||
        (curClass == BDFUT_CLASS) ||
        (curClass == IRINDEX_CLASS) ||
        (curClass == SWAPLEG_CLASS) ||
        (curClass == CMSLEG_CLASS) ||
        (curClass == CMTLEG_CLASS) ||
        (curClass == T4MLEG_CLASS) ||
        (curClass == SWAP_CLASS) ||
        (curClass == CAPFLOOR_CLASS) ||
        (curClass == RNGNOTE_CLASS) ||
        (curClass == REVERSE_CLASS) ||
        (curClass == REVERSECOUPON_CLASS) ||
        (curClass == FLEXIBLECAPFLOOR_CLASS) ||
        (curClass == OPTION_CLASS) ||
        (curClass == BARRIER_CLASS) ||
        (curClass == SWAPTION_CLASS) ||
        (curClass == SWAPTION_CAPFLOOR_CLASS) ||
        (curClass == IASEC_CLASS) ||
        (curClass == BSMODEL_CLASS) ||
        (curClass == YIELD_CURVE_BASIC_CLASS) ||
        (curClass == GYCMODEL_CLASS) ||
        (curClass == GYCLSMODEL_CLASS) ||
        (curClass == HW2FMODEL_CLASS) ||
        (curClass == HW2F_TREE_MODEL_CLASS) ||
        (curClass == DFG_YIELD_CURVE_HWTREE_CLASS) ||
        (curClass == DFBS_CLASS) ||
		(curClass == DFFXBS_CLASS) ||
	    (curClass == DFFXBSPLAIN_CLASS) ||
	    (curClass == DFHWXSIGVAR_CLASS) ||
        (curClass == DFHWSIGVARTREE_CLASS) ||
        (curClass == GAUSSIAN_2YC_MODEL_CLASS) ||
        (curClass == CRR_TREE_CLASS) ||
        (curClass == YIELD_CURVE_HWTREE_CLASS) ||
        (curClass == HWSIGCONST_TREE_CLASS) ||
        (curClass == HWSIGVAR_TREE_CLASS) ||
        (curClass == HWSIGVAR_ANALYTIC_CLASS) ||
        (curClass == HW2FSIGVAR_ANALYTIC_CLASS) ||
        (curClass == BKIRTREE_CLASS) ||
        (curClass == I3DTREE_HWTREE_CLASS) ||
        (curClass == HWFNMONTECARLO_CLASS) ||
        (curClass == HWFNMONTECARLOSV_CLASS) ||
        (curClass == HW2FFNMONTECARLOSV_CLASS) ||
        (curClass == PF_CLASS) ||
        (curClass == STRUCTURE_CLASS) ||
        (curClass == XSTYLE_CLASS) ||
        (curClass == REFVAL_CLASS) ||
		(curClass == STICKY_CLASS) ||
		(curClass == RATCHET_CLASS) ||
        (curClass == IAREFVAL_CLASS)))
	{
		return ARM_KO;
	}

	char buf[50];
    
	const char* aString = (const char *)stringObjectId;

	strcpy (buf, &(aString[5]));

	if((!buf) || (strlen (buf) == 0))
	{
		return ARM_KO;
	}

	return (atol (buf));
}



CCString GetStringObjectClass (const CCString& stringObjectId)
{
	char buf[5];

	strncpy (buf, (const char*)stringObjectId, 4);
	buf[4] = '\0';

	CCString stringObjectClass (buf);

	return (stringObjectClass);
}