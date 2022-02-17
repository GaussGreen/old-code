/* ===========================================================================
   
	 FILENAME:     swp_f_ycinstr.c

     PURPOSE:      A yield curve can be thought of as composed 
	               of sets of instruments of various types
	               (swap, cash, futs);
	                                 here they are.

   =========================================================================== */
#include "math.h"
#include "swp_h_all.h"

static Err interp_future_string(String fname,
			Date *fut_last_trading,
			Date *fut_start, 
			Date *fut_end,
			int *nm, 
			SrtCcyParam *ccy_param)
{
	SrtMonth sm;
	int sy;
	int num_wed = 0;
	char m[4];
	Err err;
	
/* get month */
	strncpy(m,fname,3);
	m[3] = '\0';
	if(err = interp_month(m,&sm)) return serror("bad future month %s",fname);	

/* get year */
 	sy = atoi(fname+3);
	if(sy < 0 || sy > 99) return serror("bad future year %s",fname);	
	if(sy < 70) sy += 2000;
	else sy += 1900;

	if(fname[5] == 'M')
		*nm = 1;
	else if(fname[5] == 'Q'){
		if ( sm != SRT_DEC && sm != SRT_MAR && 
		     sm != SRT_JUN && sm != SRT_SEP )
		{
			return serror("3m Futs start DEC,MAR,JUN or SEP. %s",fname);	
		}
		*nm = 3;
	}
	else
		return serror("future must end with M or Q: %s",fname);

	if(err = ccy_param->fut_func(*nm,ccy_param->spot_lag,
		sy,sm,fut_last_trading,fut_start, fut_end))
		return err;


	return 0;
}	


Err string_to_YCInstr(String iname, YCInstr *instr, double *rate, 
		SrtCcyParam *ccy_param, Date spot)
{
	int ny=0,nm=0,nd=0;
	int num_brk_month = 0;
	Err err;
	if(!isdigit(iname[0]))
	{
		if(err =  interp_future_string(iname,
				  &instr->dp.futdates.settle,
				  &instr->dp.futdates.start,
				  &instr->dp.futdates.end,
				  &nm,ccy_param))
			return err;
		if(nm == 1)
			instr->type = FUT1MINSTR;
		if(nm == 3)instr->type = FUT3MINSTR;      
	}
	else
	{
		if(err = interp_tenor_string(iname,&ny,&nm,&nd))
				return err;
		instr->dp.cashdates.start = spot;
		if(nd>0)
		{
			instr->dp.cashdates.end = add_unit(spot,nd,SRT_DAY,SUCCEEDING);
			instr->dp.cashdates.basis_code = ccy_param->cash_basis_code;
			instr->type = CASHINSTR;
		}
		else 
		if(ny == 0 && nm > 0 && nm <= 12)
		{
			instr->dp.cashdates.end = add_unit(spot,nm,SRT_MONTH,ccy_param->cash_bus_day_conv);
			instr->dp.cashdates.basis_code = ccy_param->cash_basis_code;
			instr->type = CASHINSTR;
		}
		else
		{
			nm = nm + 12 * ny;
			instr->dp.swapdates.basis_code = ccy_param->swap_basis_code;
			instr->dp.swapdates.compd = ccy_param->compd;
			num_brk_month = nm%(12/(ccy_param->compd));
       	    instr->dp.swapdates.end = add_unit(spot,nm,SRT_MONTH,NO_BUSDAY_CONVENTION);
			instr->dp.swapdates.nfp = 
					nm * (ccy_param->compd) / 12;

			if(num_brk_month)
				instr->dp.swapdates.first_full_fixing = 
					add_unit(spot,num_brk_month,SRT_MONTH,
						ccy_param->swap_bus_day_conv);
			else 
				instr->dp.swapdates.first_full_fixing = spot;
			
			instr->type = SWAPINSTR;
			instr->dp.swapdates.spot_lag = ccy_param->spot_lag;
			instr->dp.swapdates.direction = BKWD;
		}
	}
	if(rate)
		instr->rate = *rate;
	
	return NULL;
}	



Err YC_get_Instr(YC_Obj *yc, YCInstr *instr, int i)
{
	if(i<0)return serror("Bad instance index %d",i);
	switch(instr->type){
	  case SWAPINSTR:
		if(i >= YC_field_Dlength(yc,YC_SWAP_RATE))
			return serror("Bad swap instance index %d",i);
		instr->rate = 
			YC_field_Dget(yc,YC_SWAP_RATE,i);
		instr->dp.swapdates.basis_code = 
			(SrtBasisCode)(YC_field_Iget(yc,YC_SWAP_BASIS_CODE,i));

/** would be more aesthetically pleasing if this was a field get ***/
		instr->dp.swapdates.start = 
			YC_Dateget(yc,YC_SPOT);

		instr->dp.swapdates.end = 
			YC_field_Dateget(yc,YC_SWAPEND_DATE,i);
		instr->dp.swapdates.first_full_fixing = 
			YC_field_Dateget(yc,YC_SWAP_FIRST_FULL_FIXING,i);
		instr->dp.swapdates.nfp = 
			(int)(YC_field_Iget(yc,YC_SWAP_NUM_FULL_PERIOD,i));
		instr->dp.swapdates.compd = 
			(SrtCompounding)(YC_field_Iget(yc,YC_SWAP_COMPD,i));
		instr->dp.swapdates.direction = FWD;
		break;
	  case CASHINSTR:
		if(i >= YC_field_Dlength(yc,YC_CASH_RATE))
			return serror("Bad cash instance index %d",i);
		instr->rate = 
			(SrtBasisCode)(YC_field_Dget(yc,YC_CASH_RATE,i));
		instr->dp.cashdates.basis_code = 
			(SrtBasisCode)(YC_field_Iget(yc,YC_CASH_BASIS_CODE,i));

/** would be more aesthetically pleasing if this was a field get ***/
		instr->dp.cashdates.start = 
			YC_Dateget(yc,YC_SPOT);
		instr->dp.cashdates.end = 
			YC_field_Dateget(yc,YC_CASH_DATE,i);
		break;
	  case FUT1MINSTR:
		if(i >= YC_field_Dlength(yc,YC_M1_FUTPRICE))
			return serror("Bad fut1m instance index %d",i);
		instr->rate = 
			YC_field_Dget(yc,YC_M1_FUTPRICE,i);
		instr->dp.futdates.basis_code = 
			(SrtBasisCode)(YC_field_Iget(yc,YC_FUT1M_BASIS_CODE,i));
		instr->dp.futdates.start = 
			YC_field_Dateget(yc,YC_M1_FUTSTARTDATE,i);
		instr->dp.futdates.end = 
			YC_field_Dateget(yc,YC_M1_FUTENDDATE,i);
		instr->dp.futdates.settle = 
			YC_field_Dateget(yc,YC_M1_FUTSETTLEDATE,i);
		break;
	  case FUT3MINSTR:
		if(i >= YC_field_Dlength(yc,YC_M3_FUTPRICE))
			return serror("Bad fut3m instance index %d",i);
		instr->rate = 
			YC_field_Dget(yc,YC_M3_FUTPRICE,i);
		instr->dp.futdates.basis_code = 
			(SrtBasisCode)(YC_field_Iget(yc,YC_FUT3M_BASIS_CODE,i));
		instr->dp.futdates.start = 
			YC_field_Dateget(yc,YC_M3_FUTSTARTDATE,i);
		instr->dp.futdates.end = 
			YC_field_Dateget(yc,YC_M3_FUTENDDATE,i);
		instr->dp.futdates.settle = 
			YC_field_Dateget(yc,YC_M3_FUTSETTLEDATE,i);
		break;
	  default:
		return serror("Unknown Instrument %d",instr->type);
	}
	return NULL;
}

Err YC_set_Instr(YC_Obj *yc, YCInstr *instr, int i)
{
	if(i<0)return serror("Bad instance index %d",i);
	switch(instr->type){
	  case SWAPINSTR:
		if(i >= YC_field_Dlength(yc,YC_SWAP_RATE))
			return serror("Bad swap instance index %d",i);
		YC_field_Dset(yc,YC_SWAP_RATE,i,instr->rate);  
		YC_field_Iset(yc,YC_SWAP_BASIS_CODE,i,
			instr->dp.swapdates.basis_code);
		YC_field_Dateset(yc,YC_SWAPEND_DATE,i,
			instr->dp.swapdates.end);
		YC_field_Dateset(yc,YC_SWAP_FIRST_FULL_FIXING,i,
			instr->dp.swapdates.first_full_fixing);
		YC_field_Iset(yc,YC_SWAP_NUM_FULL_PERIOD,i,
			instr->dp.swapdates.nfp); 
		YC_field_Iset(yc,YC_SWAP_COMPD,i,
			instr->dp.swapdates.compd); 
		if(i > 0 && instr->dp.swapdates.end <= 
			YC_field_Dateget(yc,YC_SWAPEND_DATE,i-1))return
			serror("YC_set_instr swaps must be in order.");
		break;
	  case CASHINSTR:
		if(i >= YC_field_Dlength(yc,YC_CASH_RATE))
			return serror("Bad cash instance index %d",i);
		YC_field_Dset(yc,YC_CASH_RATE,i,instr->rate);
		YC_field_Iset(yc,YC_CASH_BASIS_CODE,i,
			instr->dp.cashdates.basis_code); 
		YC_field_Dateset(yc,YC_CASH_DATE,i,
			instr->dp.cashdates.end);
		if(i > 0 && instr->dp.cashdates.end <= 
			YC_field_Dateget(yc,YC_CASH_DATE,i-1))return
			serror("YC_set_instr cash rates must be in order.");
		break;
	  case FUT1MINSTR:
		if(i >= YC_field_Dlength(yc,YC_M1_FUTPRICE))
			return serror("Bad fut1m instance index %d",i);
		YC_field_Dset(yc,YC_M1_FUTPRICE,i,
			instr->rate);
		YC_field_Iset(yc,YC_FUT1M_BASIS_CODE,i,
			instr->dp.futdates.basis_code);
		YC_field_Dateset(yc,YC_M1_FUTSTARTDATE,i,
			instr->dp.futdates.start);
		YC_field_Dateset(yc,YC_M1_FUTENDDATE,i,
			instr->dp.futdates.end);
		YC_field_Dateset(yc,YC_M1_FUTSETTLEDATE,i,
			instr->dp.futdates.settle);
		if(i > 0 && instr->dp.futdates.start != 
			YC_field_Dateget(yc,YC_M1_FUTENDDATE,i-1))return
			serror("YC_set_instr fut1m dates must fit together.");
		break;
	  case FUT3MINSTR:
		if(i >= YC_field_Dlength(yc,YC_M3_FUTPRICE))
			return serror("Bad fut3m instance index %d",i);
		YC_field_Dset(yc,YC_M3_FUTPRICE,i,
			instr->rate);
		YC_field_Iset(yc,YC_FUT3M_BASIS_CODE,i,
			instr->dp.futdates.basis_code);
		YC_field_Dateset(yc,YC_M3_FUTSTARTDATE,i,
			instr->dp.futdates.start);
		YC_field_Dateset(yc,YC_M3_FUTENDDATE,i,
			instr->dp.futdates.end);
		YC_field_Dateset(yc,YC_M3_FUTSETTLEDATE,i,
			instr->dp.futdates.settle);
		if(i > 0 && instr->dp.futdates.start != 
			YC_field_Dateget(yc,YC_M3_FUTENDDATE,i-1))return
			serror("YC_set_instr fut3m dates must fit together.");
		break;
	  default:
		return serror("Unknown Instrument %d",instr->type);
	}
	return NULL;
}



/** 	return 0 if the two instruments described have
	the same date descriptors.
	(not necessarily the same rates)
**/
 
int YCInstr_comp(YCInstr *in1, YCInstr *in2)
{
	int d,b;
        int result;
	
	if(in1->type != in2->type) return -1;
	switch(in1->type){
	  case SWAPINSTR:
		d = in1->dp.swapdates.direction;
		in1->dp.swapdates.direction = in2->dp.swapdates.direction;
		b = memcmp(&(in1->dp.swapdates),&(in2->dp.swapdates),
			sizeof(SwapDP));  
		in1->dp.swapdates.direction = (SrtSwapDateDir)(d);
		result =  b;
		break;
	  case CASHINSTR:
                result = (memcmp(&(in1->dp.cashdates),&(in2->dp.cashdates),
			sizeof(CashDP)));
		break;  
	  case FUT1MINSTR:
	  case FUT3MINSTR:
                result = (memcmp(&(in1->dp.futdates),&(in2->dp.futdates),
			sizeof(FutDP)));
		break;  
	  default:
		return -1;

	}

	return result;
}  

/* 
	returns the index in yc of the YCInstr in yc of type
	instr->type, or -1 if doesn't exist.
*/
		
int YC_index_Instr(YC_Obj *yc, YCInstr *in1)
{
	int i,numinst;
	YCInstr in2;
	if(in1->type == SWAPINSTR)
		numinst = YC_Iget(yc,YC_NUMSWAP);
	if(in1->type == CASHINSTR)
		numinst = YC_Iget(yc,YC_NUMCASH);
	if(in1->type == FUT1MINSTR)
		numinst = YC_Iget(yc,YC_NUMFUT1M);
	if(in1->type == FUT3MINSTR)
		numinst = YC_Iget(yc,YC_NUMFUT3M);
	in2.type = in1->type;
		
	for(i=0;i<numinst;i++){
		YC_get_Instr(yc,&in2,i);
		if(!YCInstr_comp(in1,&in2)) return i;
	}
	return -1;
}







