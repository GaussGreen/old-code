/***
	Filename    swp_h_cmt_param.h
	author	    O Van Eyseren
	written	    Oct 23 1995

***/

#ifndef SWP_H_CMT_PARAM_H
#define SWP_H_CMT_PARAM_H

#include "swp_h_ycinstr.h"

typedef enum CMTCode { CMT1 = 1, CMT2 = 2, CMT3 = 3, CMT5 = 5, CMT7 = 7, 
			CMT10 = 10, CMT20 = 20, CMT30 = 30 ,  TEC10 = -10,
			LASTCMTCODE}CMTCode;


typedef struct{                            
	CMTCode 	cmt_code;
	long        cmt_mat;
	char		yc_name[SRTBUFSZ];
	char		vc_name[SRTBUFSZ];
	char		mkt_name[SRTBUFSZ];
	Err			(*cmt_getvol)(long,long,double,double*);
	double		flatvol;
	SrtDiffusionType voltype;
	char		swap_ref_rate[SRTBUFSZ];
	BusDayConv 	cmt_bus_day_conv;
	BusDayConv 	cms_bus_day_conv;
	BasisCode 	cmt_basis_code;
	BasisCode 	cms_basis_code;
	BasisCode 	swap_basis_code;
	BasisCode 	bond_basis_code;
	SrtCompounding 	cmt_freq;
	SrtCompounding 	cms_freq;
	SrtCompounding 	swap_compd;
	SrtCompounding 	bond_compd;
	InterpMethod 	interp_method;
	SRT_Boolean 	use_prop_vol_flg;
	SRT_Boolean 	interp_on_spread_flg;
	int		allocated;
	int 		spot_lag;
	} CMT_Param_Struct, CMT_Param, *CMT_Param_Ptr;

typedef struct{
	YCInstrType	type;   /* Will only use CASHINST and SWAPINSTR */
	double 		spread;     /* CMS_CMT spread */
	DateParam 	cmt_dp;     /* A SwapDP for the CMT leg*/ 
	DateParam 	cms_dp;     /* A SwapDP for the CMS leg*/ 
}CMTInstr;
 

CMT_Param *new_CMT_Param();

int free_CMT_Param(CMT_Param *cps);

CMT_Param_Struct *init_CMT_Param(String yc_name,
								 String vc_name,
								 String mkt_name,
								 Err (*GetVol)(long, long, double, double *), /* volatility function for reference swptns */
								 CMTCode CMT_code);

Err interp_vol_calc_method( String str, Message	*val);

Err interp_interp_on_spread( String str, Message *val);

Err interp_cmt_string(String cmt_index_string,  CMTCode *cmt_code);
 
Err CMT_string_set_param(CMT_Param_Struct *cmt_param,
                         String param_name,      
						 String param_val);

#endif
