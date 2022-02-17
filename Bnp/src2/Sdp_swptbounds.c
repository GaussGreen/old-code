/*  Author: AdA */
/*  Date:   22/06/1999 */

/*  Associated header: swptbounds.h */

/*  This wraps the srt_sdp function to solve for the market correl. matrix */
/*  ---------------------------------------------------------------------- */

#include "stdio.h"
#include "math.h"
#include "utallhdr.h"
#include "sdp_sdp.h"
#include "BGM.h"
/*  For DEBUG only *************** */
#ifdef _DEBUG
#include "sdp_sdplib.h"
#endif
/*  ****************************** */

int matu(double *vecs, int dim)
{
	int i,res=0;
	for (i=1;i<=dim;i++) if ((vecs[i]!=0.0)&&(res<1)) 
	{
		res=i;
	}
	return res;
}

Err build_swmat(double *vecs, double **matres, int dim)
{
	int i,j;

	for (i=1;i<=dim;i++) for (j=1;j<=dim;j++) matres[i][j]=vecs[i]*vecs[j];
	return NULL;
}

Err addma(double **addto, double **addit, double scal, int dim)
{
	int i,j;
	for (i=1;i<=dim;i++) for (j=1;j<=dim;j++) addto[i][j]+=scal*addit[i][j];
	return NULL;
}

Err minusma(double **matin, double **matout, int dim)
{
	int i,j;
	for (i=1;i<=dim;i++) for (j=1;j<=dim;j++) matout[i][j]=-matin[i][j];
	return NULL;
}

Err rotate_left(double *vecs, double *vecres, int shift, int dim)
{
	int i;
	for (i=1;i<=dim;i++) vecres[i]=vecs[(((i-1)+shift)%dim)+1];
	return NULL;
}

Err build_instrument(double *vecs, double **matres, int dim)
{
	double **matbuf=dmatrix(1,dim,1,dim);
	double *vecbuf=dvector(1,dim);
	int i,j;
	double matur;

	for (i=1;i<=dim;i++) for (j=1;j<=dim;j++) matres[i][j]=0;
	matur=matu(vecs,dim);
	for (i=1;i<=matur;i++)
	{
		rotate_left(vecs,vecbuf,i-1,dim);
		build_swmat(vecbuf,matbuf,dim);
		addma(matres,matbuf,1/matur,dim);
	}
	free_dmatrix(matbuf,1,dim,1,dim);
	free_dvector(vecbuf,1,dim);
	return NULL;
}
 
Err build_instrument_with_cvgs(double *vecs, double **matres, int dim, double *coverages)
{
	double **matbuf=dmatrix(1,dim,1,dim);
	double *vecbuf=dvector(1,dim);
	int i,j;
	double sumcvgs;
	double matur;
	
	sumcvgs = 0.0;
	for (i=1;i<=dim;i++) for (j=1;j<=dim;j++) matres[i][j]=0;
	matur=matu(vecs,dim);

	for (i=1;i<=matur;i++)
	{
		sumcvgs+=coverages[i];
	}


	for (i=1;i<=matur;i++)
	{
		rotate_left(vecs,vecbuf,i-1,dim);
		build_swmat(vecbuf,matbuf,dim);
		addma(matres,matbuf,coverages[i]/sumcvgs,dim);
	}
	free_dmatrix(matbuf,1,dim,1,dim);
	free_dvector(vecbuf,1,dim);
	return NULL;
}



/*  Tests a particular instrument to see if it is a caplet or not. */
/*  Returns 0 if not a caplet, else returns the caplet maturity. */
int is_caplet(double **mata, int dim)
{
	int i,j,res=1;
	for (i=1;i<=dim;i++) for (j=i+1;j<=dim;j++)
	{
		if (mata[i][j]!=0) res=0;
	}
	if (res==1)
	{
		for (i=1;i<=dim;i++) if (mata[i][i]!=0) res=i;
	}
	return res;
}


/*  Test for R^2 partial order on caplets. */
int is_bigger(double *caplet_a, double *caplet_b)
{
	int res=1;

	if (caplet_a[1]>=caplet_b[1]) res=0; 
	if (caplet_a[2]>=caplet_b[2]) res=0;
	return res;
}


/*  This fudges a caplet variance so that it fits. */
/*  It takes a Caplet which has a too low vol. and interpolates  */
/*  linearly in sigma^2 * T. */
/*  takes instruments int the form (caplet maturity, variance). */
/*  Makes sigma^2 * T linear, so it aligns the three sigma^2 * T. */
Err fudge_var(double *caplet_before, double *caplet_after, double *caplet_to_fudge)
{
	double bufa=0,bufb=0;
	Err err=NULL;

	bufa=(caplet_after[1]-caplet_before[1])/(caplet_after[2]-caplet_before[2]);
	if (bufa<=0) err=serror("Smoothing impossible ...");
	if (err==NULL)
	{
		bufb=caplet_before[1];
		caplet_to_fudge[1]=(bufb+(caplet_to_fudge[2]-caplet_before[2])*bufa);
	}
	return err;
}



/*  This takes the all calibration instruments and vars set to look for basic  */
/*  infeasibility at the Caplet level. */
/*  It tests for the presence of non-increasing sigma^2 * T and fudges the input variance */
/*  in the cases where some Caplets are out of range. These ones are adapted so that sigma^2 * T  */
/*  is closest to linear. */
/*  Succeeds if only isolated caplets are out of range, fails if series of caplets are out of range. */
Err smooth_caplets(double ***const_mat, double *const_val, int dim, int num_const, int print_level)
{
	char buffer[512];
	Err errbuf=NULL,err=NULL;
	int i,j,k=0;
	int pos_before=0,pos_after=num_const, buf_pos_before=0, buf_pos_after=0, is_last=0;
	double **vecpos=dmatrix(0,num_const,1,2);
	int has_fudged=0,fudge_point=0;

	/*  For DEBUG */
	#ifdef _DEBUG
	double *vec_debug=dvector(1,num_const);
	#endif
	/*  --------- */

	if (print_level>=1)
	{
		sprintf(buffer,"Caplets vols: Testing quasi-stationnarity ... \n");
		smessage(buffer);
	}
	/*  Start Scanning for bad caplets. */
	/*  First, initialize the vec_ppos vector. */
	for (i=1;i<=num_const;i++)
	{
		vecpos[i][2]=is_caplet(const_mat[i],dim);
		if (vecpos[i][2]!=0) vecpos[i][1]=vecpos[i][2]*const_val[i];
	}
	/*  ***** For DEBUG ******** */
	#ifdef _DEBUG
	for (k=1;k<=num_const;k++) vec_debug[k]=vecpos[k][1];
	visu_vec(vec_debug,num_const,51);
	#endif
	/*  ************************ */
	/*  Grabs all caplets one by one. */
	for (i=1;i<=num_const;i++)
	{
		pos_before=0;pos_after=num_const+1;
		/*  Then test with all shorter caplets */
		for (j=1;j<=num_const;j++) if (vecpos[j][2]<vecpos[i][2])
		{
			if (!(is_bigger(vecpos[j],vecpos[i])) && (j>pos_before) ) pos_before=j;
		}
		/*  If the caplet failed look for the nearest longer one. */
		if (pos_before!=0)
		{
			for (j=1;j<=num_const;j++) if ((vecpos[i][2]<vecpos[j][2]) && (j<pos_after)) pos_after=j;
			/*  If no longer caplet and grab another nearest shorter one. */
			if (pos_after==num_const+1) 
			{
				pos_after=pos_before; pos_before=0;
				for (j=1;j<=num_const;j++) 
				{
					if ((vecpos[j][2]<vecpos[pos_after][2]) && (vecpos[j][2]>vecpos[pos_before][2])) pos_before=j;
				}
			}
			/*  If not, directly fudge Caplets and admit it... */
			errbuf=fudge_var(vecpos[pos_before],vecpos[pos_after],vecpos[i]);
			fudge_point=i;
			/*  If the fudge of this Caplet is impossible (next Caplet vol. lower */
			/*  than that of the Caplet before) then try to lower the caplet vol. before. */
			if (errbuf!=NULL)
			{
				/*  Grab another caplet immediately before (as above) */
				buf_pos_after=pos_after; pos_after=i; buf_pos_before=pos_before; pos_before=0;
				for (j=1;j<=num_const;j++) 
				{
					if (vecpos[j][2]<vecpos[buf_pos_before][2] && (vecpos[j][2]>vecpos[pos_before][2])) pos_before=j;
				}
				/*  if it has found one, try to fudge... */
				if (vecpos[pos_before][2]!=0)
				{
					errbuf=fudge_var(vecpos[pos_before],vecpos[pos_after],vecpos[buf_pos_before]);
					fudge_point=buf_pos_before;
					if (errbuf!=NULL) err=errbuf;
				}
				else
				{
					/*  If not, Interpolate backward. */
					errbuf=fudge_var(vecpos[i],vecpos[buf_pos_after],vecpos[buf_pos_before]);
					fudge_point=buf_pos_before;
					if (errbuf!=NULL) err=errbuf;
				}
			}
			/*  Make public excuses about this shamefull procedure... */
			if (print_level>=1)
			{
				has_fudged=1;
				smessage("Adjusting Caplet Mat %d Y",(int)(vecpos[fudge_point][2]));
			}
		/*  go to the next caplet. */
		}
	}
	/*  Apply fudge. */
	for (i=1;i<=num_const;i++) if (vecpos[i][2]!=0) const_val[i]=(vecpos[i][1])/vecpos[i][2];
	if (print_level>=1) smessage("\n");
	/*  ***** For DEBUG ******** */
	#ifdef _DEBUG
	for (k=1;k<=num_const;k++) vec_debug[k]=vecpos[k][1];
	visu_vec(vec_debug,num_const,52);
	#endif
	/*  ************************ */
	/*  Free */
	#ifdef _DEBUG
	free_dvector(vec_debug,1,num_const);
	#endif
	free_dmatrix(vecpos,0,num_const,1,2);
	/*  Returns the error. */
	return err;
}


/*  Solves the SDP program calibrating on market normal vars. */
/*  It gets the instuments as one line per instrument (In puts have to be properly transposed) */

Err compute_bounds(double **calib_inst, double *market_vars, int ninst, double *vec_swpt, int dim, int minmax, double *result, int *return_error, double *error, int printlevel, double toler, int niter)
{
	int k;
	Err err=NULL;
	double ***calib_cube=dcube(1,ninst,1,dim,1,dim);
	double **dirmat=dmatrix(1,dim,1,dim);
	double **matbuf=dmatrix(1,dim,1,dim);
	double **matbufb=dmatrix(1,dim,1,dim);

	if (printlevel>=1)
	{
		smessage("Normal BGM: Calibrating & pricing.... \n");
	}
	build_instrument(vec_swpt,dirmat,dim);
	if (minmax==0) minusma(dirmat,dirmat,dim);
	for (k=1;k<=ninst;k++) build_instrument(calib_inst[k],calib_cube[k],dim);
	/*  Fudge Caplet vols. */
	smooth_caplets(calib_cube,market_vars,dim,ninst,printlevel);
	/*  If fudge failed exit, else proceed to SDP. */
	if (err!=NULL) 
	{
		smessage("Hopeless Caplet variances, giving up... \n");
		*result=0;
		*return_error=0;
	}
	else
	{
	err=srt_sdp(calib_cube,market_vars,dirmat,dim,ninst,matbuf,result,error,return_error,toler,printlevel,niter);
	}
	/*  Switch min-max. */
	if (minmax==0) *result=(-(*result));
	/*  Free everything... */
	free_dcube(calib_cube,1,ninst,1,dim,1,dim);
	free_dmatrix(dirmat,1,dim,1,dim);
	free_dmatrix(matbuf,1,dim,1,dim);
	free_dmatrix(matbufb,1,dim,1,dim);
	return err;
}

