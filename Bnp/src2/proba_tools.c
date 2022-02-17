#include 	"s2020_olink_dir:uf2020.h"     	/* s2020 openlink functions   */
#include        "srt$d_src_inc:srt_h_all.h"
#include        "srt$d_src_inc:opfnctns.h"
#include	"srt$d_src_inc:utallhdr.h"
#include	"proba_tools.h"
#include	"srt$d_src_inc:srt_h_gen_2020.h"

/******************************************************************************/
/*******************************************************************************
*			Private Static Declarations
*******************************************************************************/

static 	double 	NI_var;
static	double	NI_vol;
static	double	NI_mu;
static	double	NI_mu2;

/******************************************************************************/


/*******************************************************************************
*			Private Function Prototypes
*******************************************************************************/

static 	double 	int_g(
			double 	t
			);
/******************************************************************************/


/*******************************************************************************
**                      External variables for error handing
*******************************************************************************/

#pragma nostandard

EXTERND	int	srt_err_ok                 ;
EXTERND	jmp_buf	srt_err_env                ;
EXTERND	char	sort_2020_err_buf[ 300]    ;
EXTERND char	srt_dbg_alloc_buffer[ 200] ;

#pragma standard

/******************************************************************************/

 
/*******************************************************************************
**                      Error routines a la POISSON
*******************************************************************************/

/* ==========================================================================
   Declare the error routine
   ========================================================================== */

char*	eroutine(
		int	err
		)
{
	if (!err)
	{
		return "";	/* 0 corresponds to no error */
        }
	else
	{
		return ERR_BUF_2020;
	}
}


/* ==========================================================================
   Error initialise routine
   ========================================================================== */

int	ufl_init(
		void
		)
{
	ufl_fpass( desc_info);	/* pass function desscriptor table address */
	ufl_ipass( &lib_info);	/* pass library info structure address     */
	ufl_epass( eroutine);	/* pass error-handler address		   */

	return( 0);
}


/* ==========================================================================
   Error close routine
   ========================================================================== */

int	ufl_close(
		void
		)
{
	return( 0);
}

/******************************************************************************/


/*******************************************************************************
*                                                                           
* FUNCTION     : ??                                         		    
*                                                                           
* PURPOSE      : ??    					    	          
*                                                                           
* DESCRIPTION  : ??						       	    
*									    
* PARAMETERS   : ??                                                         
*                      - ??			                            
*                      - Input/output                             	    
*              : ??                                           	      	    
*                      - ??                            		         
*                      - Input/output                      		    
*              : ??	                                                    
*                      - ??		                                    
*                      - Input/output                               	    
*                                                                           
* RETURNS      : ??		                                            
*                      - ??                                   		    
*                      - Long                                               
*                                                                           
*****************************************************************************/
/************************** CODE ********************************************/
/************************ life_option out ********************************/


/************************************************************************/
/*********************PROBA *********************************************/

	/************** proba max *********************/
double proba_max_fct(double spot,double barrier,double vol, 
double mat,double r_d,double r_f, int flag)
{
double mu,u,d1,d2,price,coeff;
int pl;

if (spot>=barrier) return(0);
price=0.0;
pl=1-2*flag;
mu = r_d - r_f - pl*vol*vol/2;
u=log(barrier/spot);
coeff=exp(2*mu*u/vol/vol);
d1=(u-mu*mat)/vol/sqrt(mat);
d2=(-u-mu*mat)/vol/sqrt(mat);

price = (srt_f_optnrmdis(d1)-coeff*srt_f_optnrmdis(d2));
price *= exp(-r_d*(1-flag)*mat-r_f*flag*mat);

return(price);
}
/************************************************************************/
/************************************************************************/

	/************** proba min *********************/
double proba_min_fct(double spot,double barrier,double vol, 
double mat,double r_d, double r_f , int flag)
{
double mu,l,d1,d2,price,coeff;
int pl;

if (spot<=barrier) return(0);
price=0.0;
pl=1-2*flag;
mu = r_d - r_f - pl*vol*vol/2;
l=log(barrier/spot);
coeff=exp(2*mu*l/vol/vol);
d1=(-l+mu*mat)/vol/sqrt(mat);
d2=(l+mu*mat)/vol/sqrt(mat);

price = (srt_f_optnrmdis(d1)-coeff*srt_f_optnrmdis(d2));
price *= exp(-r_d*(1-flag)*mat-r_f*flag*mat);

return(price);
}
/************************************************************************/
/************************************************************************/

	/************ proba max explod ****************/
double proba_max_explo_fct(double spot,double barrier,double vol, 
double mat,double r_d, double r_f,int flag)
{
double mu,mu1,mu2,u,d1,d2,price,coeff;
int pl;

if (spot>=barrier) return(1);
price=0.0;
pl=1-2*flag;
mu1 = r_d - r_f - pl*vol*vol/2;
mu2 = pl*sqrt(mu1*mu1+2*vol*vol*(r_d*(1-flag)+r_f*flag));
mu=mu1-mu2;
u=log(barrier/spot);
coeff=exp(2*mu2*u/vol/vol);
d1=(u-mu2*mat)/vol/sqrt(mat);
d2=(-u-mu2*mat)/vol/sqrt(mat);

price = 1-(srt_f_optnrmdis(d1)-coeff*srt_f_optnrmdis(d2));
price *= exp(u*mu/vol/vol);

return(price);
}
/************************************************************************/
/************************************************************************/

	/************ proba min explod ****************/
double proba_min_explo_fct(double spot,double barrier,double vol, 
double mat,double r_d, double r_f,int flag)
{
double mu,mu1,mu2,l,d1,d2,price,coeff;
int pl;

if (spot<=barrier) return(1);
price=0.0;
pl=1-2*flag;
mu1 = r_d - r_f - pl*vol*vol/2;
mu2 = pl*sqrt(mu1*mu1+2*vol*vol*(r_d*(1-flag)+r_f*flag));
mu=mu1-mu2;
l=log(barrier/spot);
coeff=exp(2*mu2*l/vol/vol);
d1=(-l+mu2*mat)/vol/sqrt(mat);
d2=(l+mu2*mat)/vol/sqrt(mat);

price = 1-(srt_f_optnrmdis(d1)-coeff*srt_f_optnrmdis(d2));
price *= exp(l*mu/vol/vol);

return(price);
}
/************************************************************************/
/************************************************************************/

         /************Range proba ***********************/

             /*********** Gauss method *************/
double dens_gauss_n(double var,double x1_n, double x2_n, 
		double mu,double vol,double mat)
{
double result,coeff1,d1,coeff2,d2;

coeff1=exp(mu*x1_n/(vol*vol));
d1=(var-x1_n-mu*mat)/(vol*sqrt(mat));
coeff2=exp(mu*x2_n/(vol*vol));
d2=(var-x2_n-mu*mat)/(vol*sqrt(mat));
result=coeff1*srt_f_optnrmdis(d1)-coeff2*srt_f_optnrmdis(d2);
return(result);
}

double proba_gauss_fct(double spot,double b_do, double b_up, 
double vol, double opt_mat,double r_d, double r_f , int nb_term,int flag)
{
double x1_n,x2_n,a,b,k,mu,bound_do,bound_up,price;
double price_k;
int i,pl;

if ((spot>=b_up)||(spot<=b_do)) return(0);
pl=1-2*flag;
a=log(b_up/spot);
b=log(spot/b_do);
mu = r_d - r_f - pl*vol*vol/2;
price=0.0;

bound_up = a;
bound_do = -b;

	for (i = -nb_term;i<=nb_term;i++)
 	{
        x1_n = 2*i*(a+b);
        x2_n = 2*a-x1_n;

        price_k = dens_gauss_n(bound_up,x1_n,x2_n,mu,vol,opt_mat)
  		-dens_gauss_n(bound_do,x1_n,x2_n,mu,vol,opt_mat);

         price += exp(-r_d*(1-flag)*opt_mat-r_f*flag*opt_mat)*price_k;

        }
return(price);
}
/************************************************************************/
/************************************************************************/

/***??????????????????????????????????????????????????????????????***/
/************Range proba : Gauss method+explo***********************/
static 	double 	int_g(
			double 	t
			)
{
	static 	double 	expo;
	static	double	result;
        
	expo   = exp( -pow( NI_var - NI_mu2 * t, 2 ) /
	            ( 2 * NI_vol * NI_vol * t       ) );
	
	
	result = (NI_var+NI_mu*t)/(2*NI_vol*t*sqrt( 2 * SRT_PI * t ) ) * expo;
	
	return( result);
}



double dens_gauss_n_explo(double var,double x1_n, double x2_n, 
		double mu,double mu2,double vol,double mat)
{
double result,coeff1,integ1,coeff2,integ2;
double (*f)(double);

coeff1=exp(mu*x1_n/(vol*vol))*exp((mu-mu2)*(var-x1_n)/(vol*vol));
	NI_var      = var-x1_n;                             
	f	  = int_g;
	integ1	  = sm_qsimp( f, 0.00001, mat, 1.0e-06);
coeff2=exp(mu*x2_n/(vol*vol))*exp((mu-mu2)*(var-x2_n)/(vol*vol));
	NI_var      = var-x2_n;                             
	f	  = int_g;
	integ2	  = sm_qsimp( f, 0.00001, mat, 1.0e-06);
result = coeff1*integ1-coeff2*integ2;
return(result);
}

double proba_gauss_explo_fct(double spot,double b_do, double b_up, 
double vol, double opt_mat,double r_d, double r_f , int nb_term,int flag)
{
double x1_n,x2_n,a,b,k,mu,bound_do,bound_up,price;
double price_k,mu2;
int i,pl;

if ((spot>=b_up)||(spot<=b_do)) return(1);
pl=1-2*flag;
a=log(b_up/spot);
b=log(spot/b_do);
mu = r_d - r_f - pl*vol*vol/2;
mu2 = sqrt(mu*mu+2*vol*vol*(r_d*(1-flag)+r_f*flag));
	NI_vol    = vol;
	NI_mu = mu;
	NI_mu2 = mu2;
price=0.0;

bound_up = a;
bound_do = -b;

	for (i = -nb_term;i<=nb_term;i++)
 	{
        x1_n = 2*i*(a+b);
        x2_n = 2*a-x1_n;

        price_k = dens_gauss_n_explo(bound_up,x1_n,x2_n,mu,mu2,vol,opt_mat)
  		-dens_gauss_n_explo(bound_do,x1_n,x2_n,mu,mu2,vol,opt_mat);

         price += price_k;

        }
return(price);
}

/***??????????????????????????????????????????????????????????????***/

/************************************************************************/
/************************************************************************/
                   /****** fourier method *******/
double dens_four_n(int n,double x,double a, double b,double coeff1, double coeff2)

{
double var,result;

var = (double)(n*SRT_PI*(x+b)/(a+b));
result = exp(x*coeff1);
result /= (coeff1*coeff1+coeff2*coeff2);
result *= (coeff1*sin(var)-coeff2*cos(var));

return(result);
} 

double proba_four_fct(double spot, double b_do, double b_up, 
double vol, double opt_mat,double r_d, double r_f , int nb_term,int flag)

{
double a,b,k,four_n,bound_do,bound_up,price,coeff1,coeff2;
double mu,lam_n,price_k;
int i,pl;

pl=1-2*flag;
a=log(b_up/spot);
b=log(spot/b_do);
mu = r_d - r_f - pl*vol*vol/2;
price=0.0;

if ((spot<=b_do)||(spot>=b_up)) return(0);

bound_up = a;
bound_do = -b;

	 for (i=1;i<=nb_term;i++)
 	 {
	coeff2 = (double)(i*SRT_PI/(a+b));
        four_n = 2*sin(b*coeff2)/(a+b);

	coeff1 = mu/vol/vol;
	lam_n=(mu*coeff1+coeff2*coeff2*vol*vol)/2;
	price_k = four_n*exp(-lam_n*opt_mat)*
	  (dens_four_n(i,bound_up,a,b,coeff1,coeff2)-dens_four_n(i,bound_do,a,b,coeff1,coeff2)); 

 
         price += exp(-r_d*(1-flag)*opt_mat-r_f*flag*opt_mat)*price_k;
         }

return(price);
}
/************************************************************************/
/************************************************************************/

/********** Proba to touch b_yes and not b_no : fourier ***************/

double probaf_yesno_fct(double spot, double b_yes, double b_no, 
double vol, double mat,double r_d, double r_f , int nb_term, int flag)

{
double a,b,price,vn,yn;
double mu,term;
int i,pl;

pl = 1-2*flag;
a = log(b_no/spot);
b = log(spot/b_yes);
mu = r_d - r_f - pl*vol*vol/2;
price = 0.0;

	 for (i=1;i<=nb_term;i++)
 	 {
         vn = -mu*mu/(2*vol*vol)-i*i*SRT_PI*SRT_PI*vol*vol/(2*(a+b)*(a+b));
         yn = exp(vn*mat)-1;
	term = pow(-1,i+1)*(SRT_PI*vol*vol/(a+b)/(a+b))*i*yn/vn*sin(SRT_PI*a/(a+b));
         price += exp(-r_d*(1-flag)*mat-r_f*flag*mat)*exp(-b*mu/vol/vol)*term;
         }

return(price);
}
/************************************************************************/
/************************************************************************/

/********** Proba to touch b_yes and not b_no : fourier+explo ***************/

double probaf_yesno_explo_fct(double spot, double b_yes, double b_no, 
double vol, double mat,double r_d, double r_f , int nb_term)

{
double a,b,price,vn,yn;
double mu,term;
int i,pl;

a=log(b_no/spot);
b=log(spot/b_yes);
mu = r_d - r_f - vol*vol/2;
price=0.0;

	 for (i=1;i<=nb_term;i++)
 	 {
         vn = -r_d-mu*mu/(2*vol*vol)-i*i*SRT_PI*SRT_PI*vol*vol/(2*(a+b)*(a+b));
         yn = exp(vn*mat)-1;
	term = pow(-1,i+1)*(SRT_PI*vol*vol/(a+b)/(a+b))*i*yn/vn*sin(SRT_PI*a/(a+b));
         price += exp(-b*mu/vol/vol)*term;
         }

return(price);
}

/************************************************************************/
/************************************************************************/

/********** Proba to touch b_yes and not b_no : gaussian ***************/
double probag_yesno_fct(double spot, double b_yes, double b_no, 
double vol, double mat,double r_d, double r_f , int nb_term, int flag)

{
double a,b,price,yn,coeff_tot,coeff;
double mu,term,d1,d2;
int i,pl,mult,mult2;


if (b_yes<b_no) 
	{
	if (b_yes>=spot) return(exp(-r_d*(1-flag)*mat-r_f*flag*mat));
	if (b_no<=spot) return(0);
	mult2 = 1;
        }
	else 
	{
	if (b_yes<=spot) return(exp(-r_d*(1-flag)*mat-r_f*flag*mat));
	if (b_no>=spot) return(0);
	mult2 = -1;
	}
pl 	= 1-2*flag; 
a  	= log(b_yes/spot);
b  	= log(spot/b_no);
mu 	= r_d - r_f - pl*vol*vol/2;
price 	= 0.0;

	 for (i = -nb_term;i<=nb_term;i++)
 	 {
         yn = -(2*i+1)*(a+b)+b;
	if (yn>0) mult=1;
	          else mult = -1;
	coeff=exp(-2*mu*yn/vol/vol);
	d1 = mult*(-yn-mu*mat)/vol/sqrt(mat);
	d2 = mult*(-yn+mu*mat)/vol/sqrt(mat);
         coeff_tot = exp(-mu*2*i*(a+b)/vol/vol);
	if (yn==0) term = -0.5;
	else term = mult*(srt_f_optnrmdis(d1)+coeff*srt_f_optnrmdis(d2));
         price += exp(-r_d*(1-flag)*mat-r_f*flag*mat)*term*coeff_tot;
         }

return(mult2*price);
}
/************************************************************************/
/************************************************************************/

/********** Proba to touch b_yes and not b_no : gaussian+explo ***************/
double probag_yesno_explo_fct(double spot, double b_yes, double b_no, 
double vol, double mat,double r_d, double r_f , int nb_term, int flag)

{
double a,b,price,yn,coeff_tot,coeff;
double mu,term,d1,d2,mu2;
int i,mult,mult2,pl;

if (b_yes<b_no) 
	{
	if (b_yes>=spot) return(1);
	if (b_no<=spot) return(0);
	mult2 = 1;
        }
	else 
	{
	if (b_yes<=spot) return(1);
	if (b_no>=spot) return(0);
	mult2 = -1;
	}

pl 	= 1-2*flag; 
a = log(b_yes/spot);
b = log(spot/b_no);
mu = r_d - r_f - pl*vol*vol/2;
mu2 = -sqrt(mu*mu+2*vol*vol*(r_d*(1-flag)+r_f*flag));
price = 0.0;

	 for (i = -nb_term;i<=nb_term;i++)
 	 {
         yn = -(2*i+1)*(a+b)+b;
	if (yn>0) mult=1;
	          else mult = -1;
	coeff=exp(-2*mu2*yn/vol/vol);
	d1 = mult*(-yn-mu2*mat)/vol/sqrt(mat);
	d2 = mult*(-yn+mu2*mat)/vol/sqrt(mat);
         coeff_tot = exp(-mu*2*i*(a+b)/vol/vol)*exp(yn*(mu2-mu)/vol/vol);
	if (yn==0) term = -0.5;
	else term = mult*(srt_f_optnrmdis(d1)+coeff*srt_f_optnrmdis(d2));
         price += term*coeff_tot;
         }

return(mult2*price);
}

/***************************************************************************/
/************************ Interface ***********************************/
/***************************************************************************/


int 	proba_max_addin	(
				int 	argc, 
				UF_ARG 	*argv, 
				UF_ARG 	*ret
				)
{
  double spot,barrier,vol, mat,rd,rf; 
  int flag,dbg=0;	
  double premium;
	
	int cur_index=0;
	Err err;
        double answer;

INIT_DOUBLE(spot);
INIT_DOUBLE(barrier);
INIT_DOUBLE(vol);
INIT_DOUBLE(mat);
INIT_DOUBLE(rd);
INIT_DOUBLE(rf);
INIT_INT(flag);

/***        if(argc > cur_index)
		INIT_INT(dbg);
	if(dbg)	lib$signal(SS$_DEBUG);***/

if(spot>barrier) error_2020(" it is worth ZERO");
answer=proba_max_fct(spot,barrier,vol,mat,rd, rf,flag); 

   if (err = srt_f_wrtdbl2020 ( answer , ret ))              
	{
		error_2020(err);
	}                

      return(0);
}

int 	proba_min_addin	(
				int 	argc, 
				UF_ARG 	*argv, 
				UF_ARG 	*ret
				)
{
  double spot,barrier,vol, mat,rd,rf; 
  int flag,dbg=0;	
  double premium;
	
	int cur_index=0;
	Err err;
        double answer;

INIT_DOUBLE(spot);
INIT_DOUBLE(barrier);
INIT_DOUBLE(vol);
INIT_DOUBLE(mat);
INIT_DOUBLE(rd);
INIT_DOUBLE(rf);
INIT_INT(flag);

/***        if(argc > cur_index)
		INIT_INT(dbg);
	if(dbg)	lib$signal(SS$_DEBUG);    ****/

if(spot<barrier) error_2020(" it is worth ZERO");
answer=proba_min_fct(spot,barrier,vol,mat,rd, rf,flag);

   if (err = srt_f_wrtdbl2020 ( answer , ret ))              
	{
		error_2020(err);
	}                   

      return(0);
}

int 	proba_max_explo_addin	(
				int 	argc, 
				UF_ARG 	*argv, 
				UF_ARG 	*ret
				)
{
  double spot,barrier,vol, mat,rd,rf; 
  int flag,dbg=0;	
  double premium;
	
	int cur_index=0;
	Err err;
        double answer;

INIT_DOUBLE(spot);
INIT_DOUBLE(barrier);
INIT_DOUBLE(vol);
INIT_DOUBLE(mat);
INIT_DOUBLE(rd);
INIT_DOUBLE(rf);
INIT_INT(flag);

/***        if(argc > cur_index)
		INIT_INT(dbg);
	if(dbg)	lib$signal(SS$_DEBUG);***/

if(spot>barrier) error_2020(" it is worth ZERO");
answer=proba_max_explo_fct(spot,barrier,vol,mat,rd, rf,flag);

   if (err = srt_f_wrtdbl2020 ( answer , ret ))              
	{
		error_2020(err);
	}                   

      return(0);
}

int 	proba_min_explo_addin	(
				int 	argc, 
				UF_ARG 	*argv, 
				UF_ARG 	*ret
				)
{
  double spot,barrier,vol, mat,rd,rf; 
  int flag,dbg=0;	
  double premium;
	
	int cur_index=0;
	Err err;
        double answer;

INIT_DOUBLE(spot);
INIT_DOUBLE(barrier);
INIT_DOUBLE(vol);
INIT_DOUBLE(mat);
INIT_DOUBLE(rd);
INIT_DOUBLE(rf);
INIT_INT(flag);

/***        if(argc > cur_index)
		INIT_INT(dbg);
	if(dbg)	lib$signal(SS$_DEBUG); ***/

if(spot<barrier) error_2020(" it is worth ZERO");
answer=proba_min_explo_fct(spot,barrier,vol,mat,rd, rf,flag);

   if (err = srt_f_wrtdbl2020 ( answer , ret ))              
	{
		error_2020(err);
	}                   

      return(0);
}

int 	proba_gauss_addin	(
				int 	argc, 
				UF_ARG 	*argv, 
				UF_ARG 	*ret
				)
{
  double spot,strike, b_up,b_do, vol, mat,rd,rf; 
  int opt_type,num_term,flag,dbg=0;	
  double premium;
	
	int cur_index=0;
	Err err;
        double answer;

INIT_DOUBLE(spot);
INIT_DOUBLE(b_do);
INIT_DOUBLE(b_up);
INIT_DOUBLE(vol);
INIT_DOUBLE(mat);
INIT_DOUBLE(rd);
INIT_DOUBLE(rf);
INIT_INT(num_term);
INIT_INT(flag);

/***        if(argc > cur_index)
		INIT_INT(dbg);
	if(dbg)	lib$signal(SS$_DEBUG); ***/

if(b_do>b_up) error_2020(" upper barrier must be > lower barrier");
answer=proba_gauss_fct(spot,b_do,b_up, vol,mat,rd, rf,num_term,flag);

   if (err = srt_f_wrtdbl2020 ( answer , ret ))              
	{
		error_2020(err);
	}                   

      return(0);
}

int 	proba_gauss_explo_addin	(
				int 	argc, 
				UF_ARG 	*argv, 
				UF_ARG 	*ret
				)
{
  double spot,strike, b_up,b_do, vol, mat,rd,rf; 
  int opt_type,num_term,flag,dbg=0;	
  double premium;
	
	int cur_index=0;
	Err err;
        double answer;

INIT_DOUBLE(spot);
INIT_DOUBLE(b_do);
INIT_DOUBLE(b_up);
INIT_DOUBLE(vol);
INIT_DOUBLE(mat);
INIT_DOUBLE(rd);
INIT_DOUBLE(rf);
INIT_INT(num_term);
INIT_INT(flag);

/***        if(argc > cur_index)
		INIT_INT(dbg);
	if(dbg)	lib$signal(SS$_DEBUG); ***/

if(b_do>b_up) error_2020(" upper barrier must be > lower barrier");
answer=proba_gauss_explo_fct(spot,b_do,b_up, vol,mat,rd, rf,num_term,flag);

   if (err = srt_f_wrtdbl2020 ( answer , ret ))              
	{
		error_2020(err);
	}                   

      return(0);
}


int 	proba_four_addin	(
				int 	argc, 
				UF_ARG 	*argv, 
				UF_ARG 	*ret
				)
{
  double spot,strike, b_up,b_do, vol, mat,rd,rf; 
  int opt_type,num_term,flag,dbg=0;	
  double premium;
	
	int cur_index=0;
	Err err;
        double answer;

INIT_DOUBLE(spot);
INIT_DOUBLE(b_do);
INIT_DOUBLE(b_up);
INIT_DOUBLE(vol);
INIT_DOUBLE(mat);
INIT_DOUBLE(rd);
INIT_DOUBLE(rf);
INIT_INT(num_term);
INIT_INT(flag);

/***        if(argc > cur_index)
		INIT_INT(dbg);
	if(dbg)	lib$signal(SS$_DEBUG); ***/

if(b_do>b_up) error_2020(" upper barrier must be > lower barrier");
answer=proba_four_fct(spot,b_do,b_up, vol,mat,rd, rf,num_term,flag);

   if (err = srt_f_wrtdbl2020 ( answer , ret ))              
	{
		error_2020(err);
	}                   

      return(0);
}

int 	probaf_yesno_addin	(
				int 	argc, 
				UF_ARG 	*argv, 
				UF_ARG 	*ret
				)
{
  double spot,strike, b_yes,b_no, vol, mat,rd,rf; 
  int opt_type,num_term,flag,dbg=0;	
  double premium;
	
	int cur_index=0;
	Err err;
        double answer;

INIT_DOUBLE(spot);
INIT_DOUBLE(b_yes);
INIT_DOUBLE(b_no);
INIT_DOUBLE(vol);
INIT_DOUBLE(mat);
INIT_DOUBLE(rd);
INIT_DOUBLE(rf);
INIT_INT(num_term);
INIT_INT(flag);

/***        if(argc > cur_index)
		INIT_INT(dbg);
	if(dbg)	lib$signal(SS$_DEBUG);  ***/

answer=probaf_yesno_fct(spot,b_yes,b_no, vol,mat,rd, rf,num_term,flag);

   if (err = srt_f_wrtdbl2020 ( answer , ret ))              
	{
		error_2020(err);
	}                   

      return(0);
}

int 	probaf_yesno_explo_addin	(
				int 	argc, 
				UF_ARG 	*argv, 
				UF_ARG 	*ret
				)
{
  double spot,strike, b_yes,b_no, vol, mat,rd,rf; 
  int opt_type,num_term,flag,dbg=0;	
  double premium;
	
	int cur_index=0;
	Err err;
        double answer;

INIT_DOUBLE(spot);
INIT_DOUBLE(b_yes);
INIT_DOUBLE(b_no);
INIT_DOUBLE(vol);
INIT_DOUBLE(mat);
INIT_DOUBLE(rd);
INIT_DOUBLE(rf);
INIT_INT(num_term);

/***        if(argc > cur_index)
		INIT_INT(dbg);
	if(dbg)	lib$signal(SS$_DEBUG);***/

answer=probaf_yesno_explo_fct(spot,b_yes,b_no, vol,mat,rd, rf,num_term);

   if (err = srt_f_wrtdbl2020 ( answer , ret ))              
	{
		error_2020(err);
	}                   

      return(0);
}

int 	probag_yesno_addin	(
				int 	argc, 
				UF_ARG 	*argv, 
				UF_ARG 	*ret
				)
{
  double spot,strike, b_yes,b_no, vol, mat,rd,rf; 
  int opt_type,num_term,flag,dbg=0;	
  double premium;
	
	int cur_index=0;
	Err err;
        double answer;

INIT_DOUBLE(spot);
INIT_DOUBLE(b_yes);
INIT_DOUBLE(b_no);
INIT_DOUBLE(vol);
INIT_DOUBLE(mat);
INIT_DOUBLE(rd);
INIT_DOUBLE(rf);
INIT_INT(num_term);
INIT_INT(flag);

/***        if(argc > cur_index)
		INIT_INT(dbg);
	if(dbg)	lib$signal(SS$_DEBUG); ***/

if (((b_yes>spot)&&(b_no>spot))||((b_yes<spot)&&(b_no<spot))) 
error_2020("spot must be between the 2 barriers");

answer=probag_yesno_fct(spot,b_yes,b_no, vol,mat,rd, rf,num_term,flag);

   if (err = srt_f_wrtdbl2020 ( answer , ret ))              
	{
		error_2020(err);
	}                   

      return(0);
}

int 	probag_yesno_explo_addin	(
				int 	argc, 
				UF_ARG 	*argv, 
				UF_ARG 	*ret
				)
{
  double spot,strike, b_yes,b_no, vol, mat,rd,rf; 
  int opt_type,num_term,flag,dbg=0;	
  double premium;
	
	int cur_index=0;
	Err err;
        double answer;

INIT_DOUBLE(spot);
INIT_DOUBLE(b_yes);
INIT_DOUBLE(b_no);
INIT_DOUBLE(vol);
INIT_DOUBLE(mat);
INIT_DOUBLE(rd);
INIT_DOUBLE(rf);
INIT_INT(num_term);
INIT_INT(flag);

/***        if(argc > cur_index)
		INIT_INT(dbg);
	if(dbg)	lib$signal(SS$_DEBUG);   ***/

if (((b_yes>spot)&&(b_no>spot))||((b_yes<spot)&&(b_no<spot))) 
error_2020("spot must be between the 2 barriers");

answer=probag_yesno_explo_fct(spot,b_yes,b_no,vol,mat,rd,rf,num_term,flag);

   if (err = srt_f_wrtdbl2020 ( answer , ret ))              
	{
		error_2020(err);
	}                   

      return(0);
}

