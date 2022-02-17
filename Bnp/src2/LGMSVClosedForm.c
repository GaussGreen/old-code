/* ==========================================================================
   FILE_NAME:	LGMSVClosedForm.c

   PURPOSE:		
	
   DATE:		03/21/02
   
   AUTHOR:		P.A.
   ========================================================================== */

#include "LGMSVPDE.h"
#include "opfnctns.h"
#include "lgmsvclosedform.h"


#define MAX_CPN			600
#define PI2	6.28318530717958647692528676655900576839433879875021164194988918462

/* For swaping two variables */
#define LGMSVPDESwap(Type,A,B) {Type C; C=A; A=B; B=C;}

static void prod_comp(double	re1,
					  double	im1,
					  double	re2,
					  double	im2,
					  double	*re3,
					  double	*im3)
{
	*re3 = re1 * re2 - im1 * im2;
	*im3 = re1 * im2 + im1 * re2;
}

static void div_comp( double	re1,
					  double	im1,
					  double	re2,
					  double	im2,
					  double	*re3,
					  double	*im3)
{
static double div;

	div = re2 * re2 + im2 * im2;
	prod_comp(re1, im1, re2, -im2, re3, im3);
	*re3 /= div;
	*im3 /= div;
}

static void sqr_comp(	double	re1,
						double	im1,
						double	*re3,
						double	*im3)
{
static double norm;

	norm = sqrt(re1 * re1 + im1 * im1);
	*re3 = sqrt((re1 + norm) / 2.0);
	*im3 = sqrt((-re1 + norm) / 2.0);
	if (im1 < 0.0)
	{
		*im3 *= -1.0;
	}
}

/* -------------------------------------------------------------------------------------------------------------
	LGMSVSaveDensity	

  -------------------------------------------------------------------------------------------------------------- */				
void LGMSVSaveTFDensity( /* Inputs */
					 int iNbPhi,
					 int iNbft,
					 
					 double ***Density)
{	
int  iNumPhi,iNumft;
FILE *fid;
	
	fid = fopen("C:\\TFDensity.txt","wt");

	/* Save the grids */
	for (iNumPhi = 1; iNumPhi<=iNbPhi; iNumPhi++)
	{		
		for (iNumft = 1; iNumft<iNbft; iNumft=iNumft+2)
		{
			fprintf(fid," %g	",sqrt(Density[1][iNumPhi][iNumft]*Density[1][iNumPhi][iNumft]+Density[1][iNumPhi][iNumft+1]*Density[1][iNumPhi][iNumft+1]));
		}
			
		fprintf(fid,"\n");
	}

	fclose(fid);

}
/* ------------------------------------------------------------------------------------------------------------------
	From Numerical Recipes


   ------------------------------------------------------------------------------------------------------------------ */

#define NR_END 1
#define FREE_ARG char*
#undef nrerror
#define nrerror(ERROR_STRING){nrerror_print(ERROR_STRING);return NULL;}

/* modified not to exit */
static void nrerror_print(char error_text[])
/* Numerical Recipes standard error handler */
{
	FILE *f;
	f = fopen("nrerror.out","a");
	if(f)
	{
		fprintf(f,"Numerical Recipes run-time error...\n");
		fprintf(f,"%s\n",error_text);
		fclose(f);
	}

/*	fprintf(stderr,"...now exiting to system...\n");
	exit(1); */
}


#define SWAP(a,b) tempr=(a);(a)=(b);(b)=tempr

void fourn(double data[], unsigned long nn[], int ndim, int isign)
{
	int idim;
	unsigned long i1,i2,i3,i2rev,i3rev,ip1,ip2,ip3,ifp1,ifp2;
	unsigned long ibit,k1,k2,n,nprev,nrem,ntot;
	double tempi,tempr;
	double theta,wi,wpi,wpr,wr,wtemp;

	for (ntot=1,idim=1;idim<=ndim;idim++)
		ntot *= nn[idim];
	nprev=1;
	for (idim=ndim;idim>=1;idim--) {
		n=nn[idim];
		nrem=ntot/(n*nprev);
		ip1=nprev << 1;
		ip2=ip1*n;
		ip3=ip2*nrem;
		i2rev=1;
		for (i2=1;i2<=ip2;i2+=ip1) {
			if (i2 < i2rev) {
				for (i1=i2;i1<=i2+ip1-2;i1+=2) {
					for (i3=i1;i3<=ip3;i3+=ip2) {
						i3rev=i2rev+i3-i2;
						SWAP(data[i3],data[i3rev]);
						SWAP(data[i3+1],data[i3rev+1]);
					}
				}
			}
			ibit=ip2 >> 1;
			while (ibit >= ip1 && i2rev > ibit) {
				i2rev -= ibit;
				ibit >>= 1;
			}
			i2rev += ibit;
		}
		ifp1=ip1;
		while (ifp1 < ip2) {
			ifp2=ifp1 << 1;
			theta=isign*6.28318530717959/(ifp2/ip1);
			wtemp=sin(0.5*theta);
			wpr = -2.0*wtemp*wtemp;
			wpi=sin(theta);
			wr=1.0;
			wi=0.0;
			for (i3=1;i3<=ifp1;i3+=ip1) {
				for (i1=i3;i1<=i3+ip1-2;i1+=2) {
					for (i2=i1;i2<=ip3;i2+=ifp2) {
						k1=i2;
						k2=k1+ifp1;
						tempr=(double)wr*data[k2]-(double)wi*data[k2+1];
						tempi=(double)wr*data[k2+1]+(double)wi*data[k2];
						data[k2]=data[k1]-tempr;
						data[k2+1]=data[k1+1]-tempi;
						data[k1] += tempr;
						data[k1+1] += tempi;
					}
				}
				wr=(wtemp=wr)*wpr-wi*wpi+wr;
				wi=wi*wpr+wtemp*wpi+wi;
			}
			ifp1=ifp2;
		}
		nprev *= n;
	}
}
#undef SWAP

void rlft3(double ***data, double **speq, unsigned long nn1, unsigned long nn2,
	unsigned long nn3, int isign)
{
	void fourn(double data[], unsigned long nn[], int ndim, int isign);
	unsigned long i1,i2,i3,j1,j2,j3,nn[4],ii3;
	double theta,wi,wpi,wpr,wr,wtemp;
	double c1,c2,h1r,h1i,h2r,h2i;

	c1=0.5;
	c2 = -0.5*isign;
	theta=isign*(6.28318530717959/nn3);
	wtemp=sin(0.5*theta);
	wpr = -2.0*wtemp*wtemp;
	wpi=sin(theta);
	nn[1]=nn1;
	nn[2]=nn2;
	nn[3]=nn3 >> 1;
	if (isign == 1) {
		fourn(&data[1][1][1]-1,nn,3,isign);
		for (i1=1;i1<=nn1;i1++)
			for (i2=1,j2=0;i2<=nn2;i2++) {
				speq[i1][++j2]=data[i1][i2][1];
				speq[i1][++j2]=data[i1][i2][2];
			}
	}
	for (i1=1;i1<=nn1;i1++) {
		j1=(i1 != 1 ? nn1-i1+2 : 1);
		wr=1.0;
		wi=0.0;
		for (ii3=1,i3=1;i3<=(nn3>>2)+1;i3++,ii3+=2) {
			for (i2=1;i2<=nn2;i2++) {
				if (i3 == 1) {
					j2=(i2 != 1 ? ((nn2-i2)<<1)+3 : 1);
					h1r=c1*(data[i1][i2][1]+speq[j1][j2]);
					h1i=c1*(data[i1][i2][2]-speq[j1][j2+1]);
					h2i=c2*(data[i1][i2][1]-speq[j1][j2]);
					h2r= -c2*(data[i1][i2][2]+speq[j1][j2+1]);
					data[i1][i2][1]=h1r+h2r;
					data[i1][i2][2]=h1i+h2i;
					speq[j1][j2]=h1r-h2r;
					speq[j1][j2+1]=h2i-h1i;
				} else {
					j2=(i2 != 1 ? nn2-i2+2 : 1);
					j3=nn3+3-(i3<<1);
					h1r=c1*(data[i1][i2][ii3]+data[j1][j2][j3]);
					h1i=c1*(data[i1][i2][ii3+1]-data[j1][j2][j3+1]);
					h2i=c2*(data[i1][i2][ii3]-data[j1][j2][j3]);
					h2r= -c2*(data[i1][i2][ii3+1]+data[j1][j2][j3+1]);
					data[i1][i2][ii3]=h1r+wr*h2r-wi*h2i;
					data[i1][i2][ii3+1]=h1i+wr*h2i+wi*h2r;
					data[j1][j2][j3]=h1r-wr*h2r+wi*h2i;
					data[j1][j2][j3+1]= -h1i+wr*h2i+wi*h2r;
				}
			}
			wr=(wtemp=wr)*wpr-wi*wpi+wr;
			wi=wi*wpr+wtemp*wpi+wi;
		}
	}
	if (isign == -1)
		fourn(&data[1][1][1]-1,nn,3,isign);
}
/* (C) Copr. 1986-92 Numerical Recipes Software +135[)6=. */

LGMSVSolFunc1 ***LGMSVSolFunc3tensor(long nrl, long nrh, long ncl, long nch, long ndl, long ndh)
/* allocate a double 3tensor with range t[nrl..nrh][ncl..nch][ndl..ndh] */
{
	long i,j,nrow=nrh-nrl+1,ncol=nch-ncl+1,ndep=ndh-ndl+1;
	LGMSVSolFunc1 ***t;

	/* allocate pointers to pointers to rows */
	t=(LGMSVSolFunc1 ***) srt_calloc(nrow+NR_END,sizeof(LGMSVSolFunc1**));
	if (!t) nrerror("allocation failure 1 in LGMSVSolFunc3tensor()");
	t += NR_END;
	t -= nrl;

	/* allocate pointers to rows and set pointers to them */
	t[nrl]=(LGMSVSolFunc1 **) srt_calloc(nrow*ncol+NR_END,sizeof(LGMSVSolFunc1*));
	if (!t[nrl])
	{
		free((FREE_ARG) (t+nrl-NR_END));
		nrerror("allocation failure 2 in LGMSVSolFunc3tensor()");
	}
	t[nrl] += NR_END;
	t[nrl] -= ncl;

	/* allocate rows and set pointers to them */
	t[nrl][ncl]=(LGMSVSolFunc1 *) srt_calloc(nrow*ncol*ndep+NR_END,sizeof(LGMSVSolFunc1));
	if (!t[nrl][ncl])
	{
		free((FREE_ARG) (t[nrl]+ncl-NR_END));
		free((FREE_ARG) (t+nrl-NR_END));
		nrerror("allocation failure 3 in LGMSVSolFunc3tensor()");
	}
	t[nrl][ncl] += NR_END;
	t[nrl][ncl] -= ndl;

	for(j=ncl+1;j<=nch;j++) t[nrl][j]=t[nrl][j-1]+ndep;
	for(i=nrl+1;i<=nrh;i++) {
		t[i]=t[i-1]+ncol;
		t[i][ncl]=t[i-1][ncl]+ncol*ndep;
		for(j=ncl+1;j<=nch;j++) t[i][j]=t[i][j-1]+ndep;
	}

	/* return pointer to array of pointers to rows */
	return t;
}

void free_LGMSVSolFunc3tensor(LGMSVSolFunc1 ***t, long nrl, long nrh, long ncl, long nch,
	long ndl, long ndh)
/* free a double LGMSVSolFunc3tensor allocated by LGMSVSolFunc3tensor() */
{
	free((FREE_ARG) (t[nrl][ncl]+ndl-NR_END));
	free((FREE_ARG) (t[nrl]+ncl-NR_END));
	free((FREE_ARG) (t+nrl-NR_END));
	t = NULL;
}

/* -------------------------------------------------------------------------------------------------------------
	LGMSVFillUsefullMoment	

  -------------------------------------------------------------------------------------------------------------- */
void LGMSVFillUsefullMoment(/* Input */
							int i,
							int j,
							int k,

							/* Output */
							long ***IsUsefullMoment)
{
	IsUsefullMoment[i][j][k]=LGMSV_TRUE;
	if (i>0)
	{
		LGMSVFillUsefullMoment(i-1,j,k+1,IsUsefullMoment);
		LGMSVFillUsefullMoment(i-1,j,k,IsUsefullMoment);
	}
	if (k>0)
		LGMSVFillUsefullMoment(i,j,k-1,IsUsefullMoment);

	if (j>0)
	{	
		LGMSVFillUsefullMoment(i,j-1,k,IsUsefullMoment);
		if (j>1)
			LGMSVFillUsefullMoment(i,j-2,k+1,IsUsefullMoment);
	}	
}


/* -------------------------------------------------------------------------------------------------------------
	LGMSVUpdateMomentValue	

  -------------------------------------------------------------------------------------------------------------- */
void LGMSVUpdateMomentValue(/* Input */
							LGMSVSolFunc1 ***MomentFunc,
							double	dDt,
							int		iNbCumulant,
							long	***IsUsefullMoment,

							double  *lambdaArray,

							/* Output */
							double  ***MomentValue)
{
/* Declaration of local variables */
int i,j,k;
	
	for (i=0;i<=iNbCumulant;i++)
		for (j=0;j<=iNbCumulant;j++)
			for (k=0;k<=iNbCumulant;k++)
				if (IsUsefullMoment[i][j][k]=LGMSV_TRUE)
						MomentValue[i][j][k] = LGMSVFuncValue1(&MomentFunc[i][j][k],dDt,lambdaArray,0);
}
	

/* -------------------------------------------------------------------------------------------------------------
	LGMSVFillDescription	

  -------------------------------------------------------------------------------------------------------------- */
void LGMSVFillDescription(/* Input */
							int		iNbCumulant,

							double dLambdaX,
							double dAlpha,
							double dLambdaEps,
							double dRho,

							long	***IsUsefullMoment,
							double  ***MomentValue,
							double	dSigt,

							/* Output */
							LGMSVSolFunc1 ***MomentFunc)
{
/* Declaration of local variables */
int i,j,k;
double dSigt2, dAlphaRhoSigt, dAlpha2;

	dSigt2 = dSigt*dSigt;
	dAlphaRhoSigt = dAlpha*dRho*dSigt;
	dAlpha2 = dAlpha*dAlpha;

	for (i=0;i<=iNbCumulant;i++)
		for (j=0;j<=iNbCumulant;j++)
			for (k=0;k<=iNbCumulant;k++)
			{
				if (IsUsefullMoment[i][j][k]=LGMSV_TRUE)
				{
					MomentFunc[i][j][k].dXt1 = MomentValue[i][j][k];
					MomentFunc[i][j][k].dLambda = (2.0*i+j)*dLambdaX+k*dLambdaEps;
					MomentFunc[i][j][k].iNbFunction = 5;
					
					if (i>0)
					{
						MomentFunc[i][j][k].bIsft1[0]=0;
						MomentFunc[i][j][k].a[0]=i*dSigt2;
						MomentFunc[i][j][k].pft[0]=&MomentFunc[i-1][j][k+1];

						MomentFunc[i][j][k].bIsft1[1]=0;
						MomentFunc[i][j][k].a[1]=-i*dSigt2;
						MomentFunc[i][j][k].pft[1]=&MomentFunc[i-1][j][k];
					}
					else
					{
						MomentFunc[i][j][k].bIsft1[0]=1;
						MomentFunc[i][j][k].a[0]=0;

						MomentFunc[i][j][k].bIsft1[1]=1;
						MomentFunc[i][j][k].a[1]=0;
					}

					if (k>0)
					{
						MomentFunc[i][j][k].bIsft1[2]=0;
						MomentFunc[i][j][k].a[2]=k*(dLambdaEps+(k-1.0)/2.0*dAlpha2);
						MomentFunc[i][j][k].pft[2]=&MomentFunc[i][j][k-1];
					}
					else
					{
						MomentFunc[i][j][k].bIsft1[2]=1;
						MomentFunc[i][j][k].a[2]=0;
					}

					if (j>0)
					{
						MomentFunc[i][j][k].bIsft1[3]=0;
						MomentFunc[i][j][k].a[3]=j*k*dAlphaRhoSigt;
						MomentFunc[i][j][k].pft[3]=&MomentFunc[i][j-1][k];

						if (j>1)
						{
							MomentFunc[i][j][k].bIsft1[4]=0;
							MomentFunc[i][j][k].a[4]=j*(j-1.0)/2.0*dSigt2;
							MomentFunc[i][j][k].pft[4]=&MomentFunc[i][j-2][k+1];
						}
						else
						{
							MomentFunc[i][j][k].bIsft1[4]=1;
							MomentFunc[i][j][k].a[4]=0;
						}
					}
					else
					{
						MomentFunc[i][j][k].bIsft1[3]=1;
						MomentFunc[i][j][k].a[3]=0;

						MomentFunc[i][j][k].bIsft1[4]=1;
						MomentFunc[i][j][k].a[4]=0;
					}	
				}
			}
}

/* -------------------------------------------------------------------------------------------------------------
	LGMSVMomentCalculation	

  -------------------------------------------------------------------------------------------------------------- */
void LGMSVMomentCalculation(
					 /* Inputs */
					 int	iNbCumulant,	/* Number of cumulant for the approximation of the (f,Phi) density */

					 /* Parameter of diffusion */
					 /* BE CAREFULL : Alpha and LambdaEps for V=Eps^2 */
					 double dLambdaX,
					 
					 double dAlpha,
					 double dLambdaEps,

					 double dRho,

					 double dTExercice,		/* Exercice of the swaption in years from today */
					 double dTStar,			/* Tstar */
					 
					 int	iNbSigTime,		/* Term Structure of g(t) */
					 double *SigTime,		
					 double	*Sig,

					 double *lambdaArray,

					 /* Outputs */
					 long	***IsUsefullMoment,
					 LGMSVSolFunc1 ***MomentFunc,
					 double ***MomentValue,
					 double *pPhitMean)
{
/* Declaration of locals variables */
int i, j, k;
int iNumOrder;
int iNumTime, iIndExercice;

double dDt;
double dInvDerBeta, dCoef;

LGMSVSolFunc1 FuncPhitMean;
double dPhitMean;
double dTstart, dTend;

	/* Definition of the function E(Phit) */
	
	FuncPhitMean.dLambda = 2*dLambdaX;
	FuncPhitMean.iNbFunction = 1;
	FuncPhitMean.bIsft1[0]=1;
	FuncPhitMean.a[0]=1;

	/* Init the useful moment to false */
	for (i=0;i<=iNbCumulant;i++)
		for (j=0;j<=iNbCumulant;j++)
			for (k=0;k<=iNbCumulant;k++)
				IsUsefullMoment[i][j][k]=LGMSV_FALSE;

	/* Fill recursively the usefull moment */
	/* Such that we can calculate all Moment(i,N-i) */
	for (iNumOrder=1; iNumOrder<=iNbCumulant; iNumOrder++)
		for (i=0;i<=iNumOrder;i++)
			LGMSVFillUsefullMoment(i,iNumOrder-i,0,IsUsefullMoment);

	/* Init the value of all the usefull moment at t=0 */
	for (i=0;i<=iNbCumulant;i++)
		for (j=0;j<=iNbCumulant;j++)
			for (k=0;k<=iNbCumulant;k++)
			{
				if ((k>=0) && (i==0) && (j==0))
					MomentValue[i][j][k]=1;
				else
					MomentValue[i][j][k]=0;
			}

	dPhitMean = 0; 
	
	/* Calculation of all the intermediaire */	
	iNumTime = 0;
	while ((iNumTime<iNbSigTime-1)&&(SigTime[iNumTime+1]<(dTExercice + 1.0E-10)))
		iNumTime++;
	iIndExercice = iNumTime ;

	dTend = 0.0;
			
	for (iNumTime=0; iNumTime<=iIndExercice; iNumTime++)
	{
		/* Start = End */
		dTstart = dTend;

		dTend = SigTime[iNumTime];

		if (iNumTime >0)
		{
			/* Update of the value of all the useful moment at time SigTime[iNumTime]*/
			LGMSVUpdateMomentValue(/* Input */
									MomentFunc,
									dDt,
									iNbCumulant,
									IsUsefullMoment,

									lambdaArray,

									/* Output */
									MomentValue);

			/* Evaluation of E(Phit) at t=SigTime[iNumTime] */
			dPhitMean = LGMSVFuncValue1(&FuncPhitMean,dDt,lambdaArray,0);
		}
		
		
		/* Fill the description  of all the usefull moments function at t */
		LGMSVFillDescription(/* Input */
							iNbCumulant,

							dLambdaX,
							dAlpha,
							dLambdaEps,
							dRho,

							IsUsefullMoment,
							MomentValue,
							Sig[iNumTime],

							/* Output */
							MomentFunc);

		/* Fill the description of the E(Phit) function */
		FuncPhitMean.dXt1 = dPhitMean;
		FuncPhitMean.a[0] = Sig[iNumTime]*Sig[iNumTime];

		/* Calculation of dDt */
		dDt = dTend-dTstart;

	}
	
	dDt = dTExercice - dTstart;

	/* Evaluation of all the moment E(PhiBar^i*f^j*V^0) = E(PhiBar^i*Y^j*V^0)*exp(-j*lx*(Tstar-t))         */
	/* 		MomentFunc(i,j,k) = E(PhiBar^i*Y^j*V^k)														   */
	/*			where PhiBar = Phi - E(Phi) and Y=exp(lx*(Tstar-t))*(f(t,Tstar)-f(0,Tstar))				   */
	dInvDerBeta = exp(-dLambdaX*(dTStar-dTExercice));
	dCoef = 1.0;
	
	for (iNumOrder=1; iNumOrder<=iNbCumulant; iNumOrder++)
	{
		dCoef = 1.0;
		for (i=0;i<=iNumOrder;i++)
		{
				MomentValue[iNumOrder-i][i][0] = dCoef*LGMSVFuncValue1(&MomentFunc[iNumOrder-i][i][0],dDt,lambdaArray,0);
				dCoef *= dInvDerBeta;
		}
	}

	/*if (IsUsefullMoment[0][0][2]==LGMSV_TRUE)
		MomentValue[0][0][2] = LGMSVFuncValue1(MomentFunc[0][0][2],dDt,lambdaArray,0); */
	
	/* Evaluation of E(Phit) at t=dTExercice */
	dPhitMean = LGMSVFuncValue1(&FuncPhitMean,dDt,lambdaArray,0);

	*pPhitMean = dPhitMean;
}							

/* -------------------------------------------------------------------------------------------------------------
	LGMSVCalculateDensity	

  -------------------------------------------------------------------------------------------------------------- */
void LGMSVCalculateDensity(/* Input */
							int iNbCumulant,
							int iNbPhi,
							int iNbft,
							
							double dPhiFreqStep,
							double dftFreqStep,
							double ***MomentValue,

							double *Combination,

							/* Outputs */
							double ***Density,
							double **Densityspeq)
{
double d2pi, d2pik_kFact, dpi;
int iIsComplex;
double dSign;
int iNumPhiFreq, iNumOrder, iNumftFreq;
int i,j;
double dPhifreq;
double dCoef;

double ***DensityTF;
double **DensityTFspeq;

	/* Initialisation */
	DensityTF = Density;
	DensityTFspeq = Densityspeq;

	/* Fill the TF of the density */
	dpi = 3.1415926535897932384626433832795;
	d2pi = 2*dpi;
	d2pik_kFact = 1;
	
	/* Init the C(i,k) for k=0 */
	for (i=-1;i<=iNbCumulant;i++)
		Combination[i]=0;
	Combination[0]=1;

	/* Init the density with 1 = order 0 */
	for (iNumPhiFreq=1;iNumPhiFreq<=iNbPhi;iNumPhiFreq++)
	{
		/* Fill frequency 0 to (iNbft/2.0-1) *Dfreq */
		for (iNumftFreq=1;iNumftFreq<=(iNbft/2.0);iNumftFreq++)
		{
			/* Real Part */
			DensityTF[1][iNumPhiFreq][2*(iNumftFreq-1)+1] = 1;
			/* Imaginary Part */
			DensityTF[1][iNumPhiFreq][2*iNumftFreq] = 0;
		}
		/* Fill frequency (iNbft/2.0) *Dfreq */
		
		/* Real Part */
		DensityTFspeq[1][2*(iNumPhiFreq-1)+1] = 1;
		/* Imaginary Part */
		DensityTFspeq[1][2*iNumPhiFreq] = 0;
	}
	
	for (iNumOrder=1;iNumOrder<=iNbCumulant;iNumOrder++)
	{
		/* For Calculation of the coefficient : (2*pi*i)^k/k! */
		d2pik_kFact = d2pik_kFact*d2pi/((double) iNumOrder);
		iIsComplex = (int) fmod(iNumOrder,2.0);
		dSign = 1.0-2.0*((int) fmod(iNumOrder/2.0,2.0));

		/* Calculation of all the C(i,iNumOrder) by Pascal Recursion */
		for (i=iNumOrder;i>=0;i--)
			Combination[i] += Combination[i-1];
		
		for (j=0; j<=iNumOrder; j++)
		{
			for (iNumPhiFreq=1;iNumPhiFreq<=iNbPhi;iNumPhiFreq++)
			{	
				dPhifreq = (iNumPhiFreq-1-iNbPhi*(int)((iNumPhiFreq-1)/(iNbPhi/2.0)))*dPhiFreqStep;
				for (iNumftFreq=1;iNumftFreq<=(iNbft/2.0);iNumftFreq++)
				{
					DensityTF[1][iNumPhiFreq][2*(iNumftFreq-1)+1+iIsComplex] += 
						d2pik_kFact*dSign*Combination[j]*pow(dftFreqStep*(iNumftFreq-1),j)*
						pow(dPhifreq,iNumOrder-j)*MomentValue[iNumOrder-j][j][0]; 
					/*DensityTF[1][iNumPhiFreq][2*(iNumftFreq-1)+1] = exp(-2*dpi*dpi*dPhifreq*dPhifreq)
						*exp(-2*dpi*dpi*dftFreqStep*(iNumftFreq-1)*dftFreqStep*(iNumftFreq-1));
					DensityTF[1][iNumPhiFreq][2*(iNumftFreq-1)+2] = 0; */
				}
				
				/* For a real function TF verify conj(H(f))=H(-f) */
				/* Or numerical TF is periodic, so H(fe/2) should be real */
				if (iIsComplex == LGMSV_FALSE)
				{
					DensityTFspeq[1][2*(iNumPhiFreq-1)+1+iIsComplex] += 
						d2pik_kFact*dSign*Combination[j]*pow(dftFreqStep*iNbft/2.0,j)*
						pow(dPhifreq,iNumOrder-j)*MomentValue[iNumOrder-j][j][0]; 
				}
				/*DensityTFspeq[1][2*(iNumPhiFreq-1)+1] = exp(-2*dpi*dpi*dPhifreq*dPhifreq)*exp(-2*dpi*dpi*dftFreqStep*(iNbft/2.0)*dftFreqStep*(iNbft/2.0));
				DensityTFspeq[1][2*(iNumPhiFreq-1)+2] = 0;*/

			}
		}
	}

	/* Normalisation such that rlft3 gives us the true ifft of TF(P(f_phi,f_ft))) */
	dCoef = dftFreqStep*dPhiFreqStep*2.0;/*(iNbPhi*iNbft/2.0); */
	
	for (iNumPhiFreq=1;iNumPhiFreq<=iNbPhi;iNumPhiFreq++)
	{
		/* frequency 0 to (iNbft/2.0-1) *Dfreq */
		for (iNumftFreq=1;iNumftFreq<=(iNbft/2.0);iNumftFreq++)
		{
			/* Real Part */
			DensityTF[1][iNumPhiFreq][2*(iNumftFreq-1)+1] *= dCoef;
			/* Imaginary Part */
			DensityTF[1][iNumPhiFreq][2*iNumftFreq] *= dCoef;
		}
		/* frequency (iNbft/2.0) *Dfreq */
		
		/* Real Part */
		DensityTFspeq[1][2*(iNumPhiFreq-1)+1] *= dCoef;
		/* Imaginary Part */
		DensityTFspeq[1][2*iNumPhiFreq] *= dCoef;
	}


/*	for (iNumPhiFreq=1;iNumPhiFreq<=iNbPhi;iNumPhiFreq++)
	{
		/* frequency 0 to (iNbft/2.0-1) *Dfreq */
/*		for (iNumftFreq=1;iNumftFreq<=(iNbft/2.0);iNumftFreq++)
		{
			/* Real Part */
/*			DensityTF[1][iNumPhiFreq][2*(iNumftFreq-1)+1] =max(DensityTF[1][iNumPhiFreq][2*(iNumftFreq-1)+1],0);
		}
		/* frequency (iNbft/2.0) *Dfreq */
		
		/* Real Part */
/*		DensityTFspeq[1][2*(iNumPhiFreq-1)+1] = max(DensityTFspeq[1][2*(iNumPhiFreq-1)+1],0);;
	}




	/* Calculation of the density by iFFT */
	rlft3(DensityTF,DensityTFspeq,1,iNbPhi,iNbft,-1);
}							
/* -------------------------------------------------------------------------------------------------------------
	LGMSVCalculateDensityCumulant	

  -------------------------------------------------------------------------------------------------------------- */
void LGMSVCalculateDensityCumulant(/* Input */
							int iNbCumulant,
							int iNbPhi,
							int iNbft,
							
							double dPhiFreqStep,
							double dftFreqStep,
							double ***MomentValue,

							double *Combination,
							double **CumulantRealPart,
							double **CumulantImagPart,
							double **XRealPart,
							double **XImagPart,
							double **TempRealPart,
							double **TempImagPart,

							/* Outputs */
							double ***Density,
							double **Densityspeq)
{
double d2pi, d2pik_kFact, dpi;
int iIsComplex;
double dSign;
int iNumPhiFreq, iNumOrder, iNumftFreq, iNumOrder1;
int i,j,j1, n;
double dPhifreq;
double dCoef;

double ***DensityTF;
double **DensityTFspeq;

double dXReal, dXImag, dCumulantReal, dCumulantImag;
double dftFreq;
double dExpReal, dImag;

	/* Initialisation */
	DensityTF = Density;
	DensityTFspeq = Densityspeq;

	/* Fill the TF of the density */
	dpi = 3.1415926535897932384626433832795;
	d2pi = 2.0*dpi;
	d2pik_kFact = 1;
	
	/* Init the C(i,k) for k=0 */
	for (i=-1;i<=iNbCumulant;i++)
		Combination[i]=0;
	Combination[0]=1;

	
	/* Fill X */
	for (iNumOrder=1;iNumOrder<=iNbCumulant;iNumOrder++)
	{
		/* For Calculation of the coefficient : (2*pi*i)^k/k! */
		iIsComplex = (int) fmod(iNumOrder,2.0);
		dSign = 1.0-2.0*((int) fmod(iNumOrder/2.0,2.0));
		d2pik_kFact = d2pik_kFact*d2pi/((double) (iNumOrder));

		/* Calculation of all the C(i,iNumOrder) by Pascal Recursion */
		for (i=iNumOrder;i>=0;i--)
			Combination[i] += Combination[i-1];
		
		for (j=0; j<=iNumOrder; j++)
		{
			if (iIsComplex == LGMSV_TRUE)
				XImagPart[iNumOrder-j][j]=d2pik_kFact*dSign*Combination[j]*MomentValue[iNumOrder-j][j][0];
			else
				XRealPart[iNumOrder-j][j]=d2pik_kFact*dSign*Combination[j]*MomentValue[iNumOrder-j][j][0];
		}
	
	}

	/* Calculation of the Cumulants */
	
	/* Init to 0 */
	for (iNumOrder=0;iNumOrder<=iNbCumulant;iNumOrder++)
	{
		for (j=0; j<=iNumOrder; j++)
		{
			CumulantRealPart[iNumOrder-j][j]=0;
			CumulantImagPart[iNumOrder-j][j]=0;
		}
	}

	/* Recurence calculation */
	for (n=iNbCumulant; n>=1;n--)
	{
		/* Compute 1-n/(n+1)*Un+1 */
		dCoef = -n/((double)(n+1));
		for (iNumOrder=0;iNumOrder<=iNbCumulant;iNumOrder++)
		{
			for (j=0; j<=iNumOrder; j++)
			{
				CumulantRealPart[iNumOrder-j][j] *=dCoef;
				CumulantImagPart[iNumOrder-j][j] *=dCoef;
			}
		}
		CumulantRealPart[0][0] +=1.0;

		/* Compute x*(1-n/(n+1)*Un+1) */

		for (iNumOrder=0;iNumOrder<=iNbCumulant;iNumOrder++)
		{
			for (j=0; j<=iNumOrder; j++)
			{
				TempRealPart[iNumOrder-j][j]=0;
				TempImagPart[iNumOrder-j][j]=0;
			}
		}
			
		for (iNumOrder=0;iNumOrder<=iNbCumulant;iNumOrder++)
		{
			for (j=0; j<=iNumOrder; j++)
			{
				dXReal=XRealPart[iNumOrder-j][j];
				dXImag=XImagPart[iNumOrder-j][j];

				for (iNumOrder1=0;iNumOrder1<=iNbCumulant-iNumOrder;iNumOrder1++)
				{
					for (j1=0; j1<=iNumOrder1; j1++)
					{
						dCumulantReal=CumulantRealPart[iNumOrder1-j1][j1];
						dCumulantImag=CumulantImagPart[iNumOrder1-j1][j1];
						
						TempRealPart[iNumOrder+iNumOrder1-j-j1][j+j1] += dXReal*dCumulantReal-dXImag*dCumulantImag;	
						TempImagPart[iNumOrder+iNumOrder1-j-j1][j+j1] += dXReal*dCumulantImag+dXImag*dCumulantReal;
					}
				}
			}
		}

		/* Swap the two arrays */
		LGMSVPDESwap(double **,TempRealPart,CumulantRealPart);
		LGMSVPDESwap(double **,TempImagPart,CumulantImagPart);
	}
	
	/* Init the log of the TF density with 0 */
	for (iNumPhiFreq=1;iNumPhiFreq<=iNbPhi;iNumPhiFreq++)
	{
		/* Fill frequency 0 to (iNbft/2.0-1) *Dfreq */
		for (iNumftFreq=1;iNumftFreq<=(iNbft/2.0);iNumftFreq++)
		{
			/* Real Part */
			DensityTF[1][iNumPhiFreq][2*(iNumftFreq-1)+1] = 0;
			/* Imaginary Part */
			DensityTF[1][iNumPhiFreq][2*iNumftFreq] = 0;
		}
		/* Fill frequency (iNbft/2.0) *Dfreq */
		
		/* Real Part */
		DensityTFspeq[1][2*(iNumPhiFreq-1)+1] = 0;
		/* Imaginary Part */
		DensityTFspeq[1][2*iNumPhiFreq] = 0;
	}

	/* Calculate the log of the TF of the density */
	for (iNumOrder=1;iNumOrder<=iNbCumulant;iNumOrder++)
	{
		for (j=0; j<=iNumOrder; j++)
		{
			for (iNumPhiFreq=1;iNumPhiFreq<=iNbPhi;iNumPhiFreq++)
			{	
				dPhifreq = (iNumPhiFreq-1-iNbPhi*(int)((iNumPhiFreq-1)/(iNbPhi/2.0)))*dPhiFreqStep;
				for (iNumftFreq=1;iNumftFreq<=(iNbft/2.0);iNumftFreq++)
				{
					dftFreq=dftFreqStep*(iNumftFreq-1);
					
					/* Real Part */
					DensityTF[1][iNumPhiFreq][2*(iNumftFreq-1)+1] += pow(dftFreq,j)*
						pow(dPhifreq,iNumOrder-j)*CumulantRealPart[iNumOrder-j][j];
					
					/* Imaginary Part */
					DensityTF[1][iNumPhiFreq][2*(iNumftFreq-1)+2] = fmod(DensityTF[1][iNumPhiFreq][2*(iNumftFreq-1)+2]+
						pow(dftFreq,j)*pow(dPhifreq,iNumOrder-j)*CumulantImagPart[iNumOrder-j][j],d2pi); 
				}
				
				/* For a real function TF verify conj(H(f))=H(-f) */
				/* Or numerical TF is periodic, so H(fe/2) should be real */
				DensityTFspeq[1][2*(iNumPhiFreq-1)+1] += pow(dftFreqStep*iNbft/2.0,j)*
						pow(dPhifreq,iNumOrder-j)*CumulantRealPart[iNumOrder-j][j]; 
			}
		}
	}


/*	
	for (iNumPhiFreq=1;iNumPhiFreq<=iNbPhi;iNumPhiFreq++)
	{	
		dPhifreq = (iNumPhiFreq-1-iNbPhi*(int)((iNumPhiFreq-1)/(iNbPhi/2.0)))*dPhiFreqStep;
				
		for (iNumftFreq=1;iNumftFreq<=(iNbft/2.0);iNumftFreq++)
		{
			dftFreq=dftFreqStep*(iNumftFreq-1);
			/* Real Part */
/*			DensityTF[1][iNumPhiFreq][2*(iNumftFreq-1)+1] =CumulantRealPart[i][j]+
					pow(dftFreq,j)*
						pow(dPhifreq,iNumOrder-j)*CumulantRealPart[iNumOrder-j][j];

			/* Imaginary Part */
/*			DensityTF[1][iNumPhiFreq][2*(iNumftFreq-1)+2] = fmod(DensityTF[1][iNumPhiFreq][2*(iNumftFreq-1)+2]+
				pow(dftFreq,j)*pow(dPhifreq,iNumOrder-j)*CumulantImagPart[iNumOrder-j][j],d2pi); 
					



		}
	
	}

*/





















	/* Save the TF density */
	LGMSVSaveTFDensity( /* Inputs */
					iNbPhi,
					 iNbft,
					 
					 DensityTF);

	/* Calculate of the exp(log(TF)) of the density */
	for (iNumPhiFreq=1;iNumPhiFreq<=iNbPhi;iNumPhiFreq++)
	{	
		for (iNumftFreq=1;iNumftFreq<=(iNbft/2.0);iNumftFreq++)
		{	
			dExpReal = exp(DensityTF[1][iNumPhiFreq][2*(iNumftFreq-1)+1]);
			dImag = DensityTF[1][iNumPhiFreq][2*(iNumftFreq-1)+2];

			/* Real Part */
			DensityTF[1][iNumPhiFreq][2*(iNumftFreq-1)+1]=cos(dImag)*dExpReal;
			/* Imaginary Part */
			DensityTF[1][iNumPhiFreq][2*(iNumftFreq-1)+2]=sin(dImag)*dExpReal;
		}
		DensityTFspeq[1][2*(iNumPhiFreq-1)+1]=exp(DensityTFspeq[1][2*(iNumPhiFreq-1)+1]);
	}

	/* Normalisation such that rlft3 gives us the true ifft of TF(P(f_phi,f_ft))) */
	dCoef = dftFreqStep*dPhiFreqStep*2.0;/*(iNbPhi*iNbft/2.0); */
	
	for (iNumPhiFreq=1;iNumPhiFreq<=iNbPhi;iNumPhiFreq++)
	{
		/* frequency 0 to (iNbft/2.0-1) *Dfreq */
		for (iNumftFreq=1;iNumftFreq<=(iNbft/2.0);iNumftFreq++)
		{
			/* Real Part */
			DensityTF[1][iNumPhiFreq][2*(iNumftFreq-1)+1] *= dCoef;
			/* Imaginary Part */
			DensityTF[1][iNumPhiFreq][2*iNumftFreq] *= dCoef;
		}
		/* frequency (iNbft/2.0) *Dfreq */
		
		/* Real Part */
		DensityTFspeq[1][2*(iNumPhiFreq-1)+1] *= dCoef;
		/* Imaginary Part */
		DensityTFspeq[1][2*iNumPhiFreq] *= dCoef;
	}

	/* Calculation of the density by iFFT */
	rlft3(DensityTF,DensityTFspeq,1,iNbPhi,iNbft,-1);
}



/* -------------------------------------------------------------------------------------------------------------
	LGMSVFillPayoff	

  -------------------------------------------------------------------------------------------------------------- */
void LGMSVFillPayoff(/* Input */
						int iNbPhi,
						int iNbft,

						double dLambdaX,

						double dPhitMean,
						
						int iIndexPhiMean,
						double dPhiStep,

						int iIndexft0,
						double dftStep,
						
						double dExTime,	/* Exercice of the swaption in years from today */
						
						double dTStar,			/* Tstar */						

						int	iNbCoupon,		/* Description of the cashflows */
						double *CouponTime,
						double *Coupon,						
						
						/* Outputs */
						double **Payoff)
{
int		iNumDff,iNumPhi,iNumft  ;
double	dDerBetaTexTstar, dBetaTexTstar, dBetaDerBetaTexTstar, dBeta2TexTstar;
double	dBetaTexTi, dBeta2TexTi;
double	dCoef;
int		LimPhi, LimFt;
double	CoefPhi, CoefF, ExpCoefPhi, ExpCoefF, CoefFLim, CoefPhiLim;	

	/* Calculation of beta(Tex,Tstar)*D/DTstar ( beta(Tex,Tstar)) and D/DTstar ( beta(Tex,Tstar)) */
	/* Use for constructing X(Tex) from Phi(Tex) and f(Tex,Tstar)-f(0,Tstar)                      */
	dDerBetaTexTstar =  exp(-dLambdaX * (dTStar - dExTime));
	dBetaTexTstar = (1.0 - dDerBetaTexTstar) / dLambdaX;
	dBetaDerBetaTexTstar = dDerBetaTexTstar * dBetaTexTstar;
	dBeta2TexTstar = dBetaTexTstar * dBetaTexTstar;
		
	LimPhi = iNbPhi - iIndexPhiMean + 1;
	LimFt = iNbft - iIndexft0 + 1;
	
	/* Fill the Payoff to be calculated by expectation */	
	/*	Payoff(Today) = B(Today,Tstar)*Expectation under QTstar [Payoff(Tex)/B(Tex,Tstar) |FToday]					*/
	/*																												*/
	/* Payoff(Tex)/B(Tex,Tstar)*B(Today,Tstar) = Sum of Ci*DF(Tex,Ti)/DF(Tex,Tstar)*DF(Today,Tstar)					*/																		
	/*										= Sum of Ci*DFF(Tex,Ti)/DFF(Tex,Tstar)*DF(Today,Tstar)*					*/
	/*				exp(-(beta(Tex,Ti)^2-beta(Tex,Tstar)^2)*Phi(Ti)-(beta(Tex,Ti)-beta(Tex,Tstar))*X(Ti))			*/
	/* Where (D/DTstar(Beta(Ti,Tstar)))*X(Ti) =																		*/
	/*	f(Ti,Tstar)-f(0,Tstar)-(D/DTstar(Beta(Ti,Tstar))*Beta(Ti,Tstar))(Ti,Tstar)									*/

	/* Initialisation */
	for (iNumPhi=1; iNumPhi<=iNbPhi; iNumPhi++)
	{
		memset(Payoff[iNumPhi], 0, iNbft * sizeof(double));
	}

	for (iNumDff=0;iNumDff<iNbCoupon;iNumDff++)
	{		
		dBetaTexTi = (1.0 - exp(-dLambdaX * (CouponTime[iNumDff] - dExTime))) / dLambdaX;
		dBeta2TexTi = dBetaTexTi * dBetaTexTi;

		CoefF = -(dBetaTexTi - dBetaTexTstar) / dDerBetaTexTstar;
		CoefPhi = -0.5 * (dBeta2TexTi - dBeta2TexTstar) - CoefF * dBetaDerBetaTexTstar;
		ExpCoefF = exp(CoefF * dftStep);
		ExpCoefPhi = exp(CoefPhi * dPhiStep);
		CoefFLim = exp(-iNbft * CoefF * dftStep);
		CoefPhiLim = exp(-iNbPhi * CoefPhi * dPhiStep);

		dCoef = exp(CoefPhi * dPhitMean);

		for (iNumPhi=1; iNumPhi<=iNbPhi; iNumPhi++)
		{
			for (iNumft=1; iNumft<=iNbft; iNumft++)
			{																	
				Payoff[iNumPhi][iNumft] += Coupon[iNumDff] * dCoef;
				
				dCoef *= ExpCoefF;
				if (iNumft == LimFt)
				{					
					dCoef *= CoefFLim;
				}
			}
			
			dCoef *= ExpCoefPhi;			
			if (iNumPhi == LimPhi)
			{
				dCoef *= CoefPhiLim;
			}
		}
	}

	/* Option part : PayOff = Max(PayOff,0) */
	for (iNumPhi=1; iNumPhi<=iNbPhi; iNumPhi++)
		for (iNumft=1; iNumft<=iNbft; iNumft++)
			Payoff[iNumPhi][iNumft] = max(Payoff[iNumPhi][iNumft],0); 
							
}

/* -------------------------------------------------------------------------------------------------------------
	LGMSVSaveDensity	

  -------------------------------------------------------------------------------------------------------------- */				
void LGMSVSaveDensity( /* Inputs */
					 int iNbPhi,
					 int iNbft,

					 int iIndexPhiMean,
					 int iIndexft0,

					 double dPhiStep,
					 double dftStep,
					 double dPhitMean,
					 double	dPhitStd,
					 
					 double ***Density)
{	
int  iNumPhi,iNumft;
int iIndexft,iIndexPhi;
double dfTexTstar, dPhiTex;
double dSum;
double	prod;

FILE *fid;
	
	fid = fopen("C:\\Density.txt","wt");
	fprintf(fid," dPhistep = 	\n");
	fprintf(fid,"%g	\n",dPhiStep);
	fprintf(fid," dftStep = \n");
	fprintf(fid," %g	\n",dftStep);
	fprintf(fid,"PhiMean	%g	PhiStd	%g\n", dPhitMean, dPhitStd);

	fprintf(fid,"(Phi,ft)		");

	for (iNumft=1;iNumft<=iNbft;iNumft++)
	{		
		dfTexTstar = (iNumft - iIndexft0) * dftStep;
		fprintf(fid," %g	",dfTexTstar);
	}
	
	fprintf(fid, "\n");
	fprintf(fid, "\n");

	prod = dPhiStep * dftStep;

	/* Save the grids */
	for (iNumPhi = 1; iNumPhi<=iNbPhi; iNumPhi++)
	{		
		dSum = 0.0;		

		iIndexPhi = iNumPhi - iIndexPhiMean + 1;
		if (iIndexPhi < 1)
		{
			iIndexPhi += iNbPhi;
		}

		dPhiTex = dPhitMean + (iNumPhi - iIndexPhiMean) * dPhiStep; 
		
		fprintf(fid," %g		",dPhiTex);

		for (iNumft = 1; iNumft<=iNbft; iNumft++)
		{ 
			iIndexft = iNumft - iIndexft0 + 1;
			if (iIndexft < 1)
			{
				iIndexft += iNbft;
			}
					   
			fprintf(fid," %g	", Density[1][iIndexPhi][iIndexft] * prod);
			dSum += Density[1][iIndexPhi][iIndexft];
		}
			
		fprintf(fid,"	%g\n", dSum * prod);
	}

	fprintf(fid, "\n");
	fprintf(fid, "		");

	for (iNumft = 1; iNumft<=iNbft; iNumft++)
	{		
		dSum = 0.0;
		iIndexft = iNumft-(iNumft>=iIndexft0)*iNbft+iNbft-(iIndexft0-1);

		for (iNumPhi = 1; iNumPhi<=iNbPhi; iNumPhi++)
		{
			dSum += Density[1][iNumPhi][iIndexft];
		}
			
		fprintf(fid,"%g	", dSum*prod);
	}

	fclose(fid);	

}

/* -------------------------------------------------------------------------------------------------------------
	LGMSVSaveMoment	

  -------------------------------------------------------------------------------------------------------------- */		
void LGMSVSaveMoment( /* Inputs */
					 int iNbCumulant,				 
					 double ***MomentValue)
{	
int  j,i;
FILE *fid;
	
	fid = fopen("C:\\Moment.txt","wt");

	/* Save the moments */
	for (i=0; i<=iNbCumulant; i++)
	{
		for (j=0;j<=iNbCumulant-i;j++)
			fprintf(fid," %g	",MomentValue[i][j][0]);
		fprintf(fid, "\n");
	}

	fclose(fid);

}

/* -------------------------------------------------------------------------------------------------------------
	LGMSVSaveMoment	

  -------------------------------------------------------------------------------------------------------------- */		
void LGMSVSavePayOff( /* Inputs */
					 int iNbPhi,
					 int iNbft,

					 int iIndexPhiMean,
					 int iIndexft0,

					 double dPhiStep,
					 double dftStep,
					 double dPhitMean,

					 double **PayOff)
{	
int  iNumft,iNumPhi;
int iIndexft,iIndexPhi;
double dfTexTstar, dPhiTex;
FILE *fid;
	
	fid = fopen("C:\\PayOff.txt","wt");

	fprintf(fid,"(Phi,ft)		");
	/* Save  */
	for (iNumft=1;iNumft<=iNbft;iNumft++)
	{
		iIndexft = iNumft-(iNumft>=iIndexft0)*iNbft+iNbft-(iIndexft0-1);
		dfTexTstar = (iIndexft-1-((iIndexft-1)>(iNbft-iIndexft0))*iNbft)*dftStep;
		fprintf(fid," %g	",dfTexTstar);
	}
	fprintf(fid, "\n");
	for (iNumPhi=1;iNumPhi<=iNbPhi;iNumPhi++)
	{
		iIndexPhi = iNumPhi-(iNumPhi>=iIndexPhiMean)*iNbPhi+iNbPhi-(iIndexPhiMean-1);
		dPhiTex = dPhitMean+(iIndexPhi-1-((iIndexPhi-1)>(iNbPhi-iIndexPhiMean))*iNbPhi)*dPhiStep; 
		
		fprintf(fid," %g		",dPhiTex);
		for (iNumft=1;iNumft<=iNbft;iNumft++)
		{
			iIndexft = iNumft-(iNumft>=iIndexft0)*iNbft+iNbft-(iIndexft0-1);
			fprintf(fid," %g	",PayOff[iIndexPhi][iIndexft]);
		}
		fprintf(fid, "\n");
	}
	fclose(fid);

}							
/* -------------------------------------------------------------------------------------------------------------
	LGMSVClosedForm	

  -------------------------------------------------------------------------------------------------------------- */
void LGMSVClosedForm(/* Inputs */
					 int	iNbCumulant,	/* Number of cumulant for the approximation of the (f,Phi) density */
					 int	iNbPhi,			/* Number of Phi : Should be a power of two */
					 int	iNbft,			/* Number of ft : Should be a power of two */	

					 /* Parameter of diffusion */
					 /* BE CAREFULL : Alpha and LambdaEps for V=Eps^2 */
					 double dLambdaX,
					 
					 double dAlpha,
					 double dLambdaEps,

					 double dRho,

					 /* Parameter of grids */
					 double	iNbSigmaPhiGridLeft,
					 double iNbSigmaPhiGridRight,


					 double	iNbSigmaftLeft,
					 double	iNbSigmaftRight,


					 long	lExDate,		/* Exercice date of the swaption  */
					 double dExTime,	/* Exercice of the swaption in years from today */
					 
					 double dTStar,			/* Tstar in years from today */
					 
					 int	iNbSigTime,		/* Term Structure of g(t) */
					 double *SigTime,		
					 double	*Sig,

					 int	iNbCoupon,		/* Description of the cashflows */
					 double	*CouponTime,
					 long   *CouponDate,
					 double *Coupon,

					 char	*cYieldCurve,	/* Yield Curve */

					 /* Outputs */
					 double *Price)
{
/* Declaration of locals variables */
Err				err = NULL;
double			lambdaArray[50];

/* For moments calculation */
LGMSVSolFunc1	***MomentFunc = NULL;
double			***MomentValue = NULL;
long			***IsUsefullMoment = NULL;
double			*Combination = NULL;

double			***Density = NULL;
double			**Densityspeq = NULL;
double			**PayOff = NULL;

double			**CumulantRealPart = NULL;
double			**CumulantImagPart = NULL;
double			**XRealPart = NULL;
double			**XImagPart = NULL;
double			**TempRealPart = NULL;
double			**TempImagPart = NULL;

double			*DFF_div_BTTstar = NULL;


double			dPhitMean, dPhitStd, dLogPhiStd, dPhiStep,dPhiMin, dNewPhiMin;
double			dftStd, dftStep, dftMin;
double			dPhiFreqStep, dftFreqStep;
int				iIndexPhiMean ,  iIndexft0;
int				iNumPhi, iNumft;

double			dIntegral;
double			Limit, RealPart, dPhifreq, NbPhiMin;

	/* Initialisation and memory allocation */

	/* At least up to second order moment */
	iNbCumulant = max(iNbCumulant,2);
	
	/* MomentFunc(i,j,k) = E(PhiBar^i*Y^j*V^k) where PhiBar = Phi - E(Phi) a nd Y=exp(lx*(Tstar-t))*(f(t,Tstar)-f(0,Tstar)) */
	MomentFunc =LGMSVSolFunc3tensor(0, iNbCumulant,
									0, iNbCumulant,
									0, iNbCumulant);

	MomentValue =f3tensor(0, iNbCumulant,
						  0, iNbCumulant,
						  0, iNbCumulant);

	IsUsefullMoment =l3tensor(0, iNbCumulant,
							0, iNbCumulant,
							0, iNbCumulant);

	Combination = dvector(-1, iNbCumulant);

	
	Density =f3tensor(1, 1,
						1, iNbPhi,
						1, iNbft);
	
	Densityspeq =dmatrix(1, 1,1, 2*iNbPhi);
	
	PayOff = dmatrix(1, iNbPhi,
					 1, iNbft);

	CumulantRealPart = dmatrix(0, iNbCumulant,0, iNbCumulant);
	CumulantImagPart = dmatrix(0, iNbCumulant,0, iNbCumulant);
	XRealPart = dmatrix(0, iNbCumulant,0, iNbCumulant);
	XImagPart = dmatrix(0, iNbCumulant,0, iNbCumulant);
	TempRealPart =  dmatrix(0, iNbCumulant,0, iNbCumulant);
	TempImagPart =  dmatrix(0, iNbCumulant,0, iNbCumulant);

	DFF_div_BTTstar = dvector(0, iNbCoupon-1);

	memset(lambdaArray, 0, 50 * sizeof(double));

	/* Gestion of allocation errors */
	if (!MomentFunc || !lambdaArray || !MomentValue || !IsUsefullMoment ||!DFF_div_BTTstar
		|| !Combination || !Density || !Densityspeq || !PayOff || !CumulantRealPart || !CumulantImagPart
		|| !XRealPart || !XImagPart || !TempRealPart || !TempImagPart)
	{
		err = "Memory allocation error (1) in LGMSVClosedForm";
		goto FREE_RETURN;
	}

	/* Calculation of all the moments */
	LGMSVMomentCalculation(
					 /* Inputs */
					 iNbCumulant,	/* Number of cumulant for the approximation of the (f,Phi) density */

					 /* Parameter of diffusion */
					 /* BE CAREFULL : Alpha and LambdaEps for V=Eps^2 */
					 dLambdaX,
					 
					 dAlpha,
					 dLambdaEps,

					 dRho,

					 dExTime,			/* Exercice of the swaption in years from today */
					 dTStar,			/* Tstar */
					 
					 iNbSigTime,		/* Term Structure of g(t) */
					 SigTime,		
					 Sig,

					 lambdaArray,

					 /* Outputs */
					 IsUsefullMoment,
					 MomentFunc,
					 MomentValue,
					 &dPhitMean);
	
	/* Save Moments */
	LGMSVSaveMoment( /* Inputs */
					 iNbCumulant,				 
					 MomentValue);

	/* Grid Phi Definition */
	dPhitStd = sqrt(MomentValue[2][0][0]); /* 1 */
	dLogPhiStd = sqrt(log(1+dPhitStd*dPhitStd/dPhitMean/dPhitMean));

/*	dPhiStep = dPhitMean*(exp(iNbSigmaPhiGrid*dLogPhiStd)-exp(-iNbSigmaPhiGrid*dLogPhiStd))/(iNbPhi-1);/*((iNbSigmaPhiGrid+iNbSigmaPhiGrid)*dPhitStd)/(iNbPhi-1);/*dPhitMean*(exp(iNbSigmaPhiGrid*dLogPhiStd)-exp(-iNbSigmaPhiGrid*dLogPhiStd))/(iNbPhi-1); */

/*	dPhiMin = dPhitMean*exp(-iNbSigmaPhiGrid*dLogPhiStd); /*dPhitMean-iNbSigmaPhiGrid*dPhitStd;/*dPhitMean*exp(-iNbSigmaPhiGrid*dLogPhiStd);*/

/*	iIndexPhiMean =(int) floor((dPhitMean-dPhiMin)/dPhiStep)+1;
	dNewPhiMin = dPhitMean-(iIndexPhiMean-1)*dPhiStep; */


	Limit = log(1.0E-03);

	dPhiStep = (iNbSigmaPhiGridLeft +iNbSigmaPhiGridRight)*dPhitStd/(iNbPhi-1);
	dPhifreq = 1.0 / (dPhiStep*iNbPhi) * 3.1416;
	RealPart = 0.0;
	/*
	LGMSVCalculateReClosedForm(dPhifreq, 0, Sig[0], dAlpha, dRho, dExTime, &RealPart);
	*/

	while (RealPart > Limit)
	{
		dPhifreq *= 1.2;
		/*
		LGMSVCalculateReClosedForm(dPhifreq, 0, Sig[0], dAlpha, dRho, dExTime, &RealPart);
		*/
	}

	/*
	dPhiStep = (iNbSigmaPhiGridLeft +iNbSigmaPhiGridRight)*dPhitStd/(iNbPhi-1);
	*/

	dPhiStep = 1.0 / (dPhifreq / 3.1416);

	/*
	dPhiMin = dPhitMean-iNbSigmaPhiGridLeft*dPhitStd;
	*/

	dPhiMin = max(dPhitMean - iNbPhi * iNbSigmaPhiGridLeft / (iNbSigmaPhiGridLeft + iNbSigmaPhiGridRight) * dPhiStep, 0.0);
	NbPhiMin = (iNbSigmaPhiGridLeft + iNbSigmaPhiGridRight) * dPhitStd / dPhiStep;
	NbPhiMin = (int) (exp(((int) (log(NbPhiMin) / log(2.0)) + 1) * log(2)) + 1.0E-16);

	iIndexPhiMean =(int) floor((dPhitMean-dPhiMin)/dPhiStep)+1;
	dNewPhiMin = dPhitMean-(iIndexPhiMean-1)*dPhiStep;
	
	/* Grid ft Definition */
	dftStd = sqrt(MomentValue[0][2][0]); /* 1 */
	dftStep = ((iNbSigmaftLeft+iNbSigmaftRight)*dftStd)/(iNbft-1);
	dftMin = 0.0 - iNbSigmaftLeft*dftStd;
	iIndexft0 =(int)(-dftMin/dftStep+1e-08)+1;
	dftMin = -(iIndexft0-1)*dftStep;

	/* Frequency step in fourier domain */
	dftFreqStep =1.0/(dftStep*iNbft);
	dPhiFreqStep = 1.0/(dPhiStep*iNbPhi);

	/* Calculation of the density (Phibar,ft) */
/*	LGMSVCalculateDensity(/* Input */
/*						iNbCumulant,
						iNbPhi,
						iNbft,

						dPhiFreqStep,
						dftFreqStep,
						MomentValue,

						Combination,

						/* Outputs */
/*						Density,
						Densityspeq); */


	if (dLambdaEps > 1.0E-12)
	{
		/*
		LGMSVCalculateDensityCumulant(
								iNbCumulant,
								iNbPhi,
								iNbft,
								
								dPhiFreqStep,
								dftFreqStep,
								MomentValue,

								Combination,
								CumulantRealPart,
								CumulantImagPart,
								XRealPart,
								XImagPart,
								TempRealPart,
								TempImagPart,								
								Density,
								Densityspeq);
		*/

		/*
		LGMSVCalculateDensityPDE(iNbCumulant,
								iNbPhi,
								iNbft,
								
								dPhiFreqStep,
								dftFreqStep,

								dLambdaX,
								Sig[0],
								dAlpha,
								dLambdaEps,
								dRho,
								dTStar,
								dExTime,
								50,
								Density,
								Densityspeq);
		*/
	}
	else
	{
		/*
		LGMSVCalculateDensityClosedForm(iNbCumulant,
										iNbPhi,
										iNbft,
										
										dPhiFreqStep,
										dftFreqStep,
										Sig[0],
										dAlpha,
										dRho,
										dExTime,
										Density,
										Densityspeq);
		*/
	}

	LGMSVSaveDensity( /* Inputs */
					 iNbPhi,
					 iNbft,
					 iIndexPhiMean,
					 iIndexft0,
					 dPhiStep,
					 dftStep,
					 dPhitMean,
					 dPhitStd,
					 Density);

	/* Fill The Payoff of the swaption */
	LGMSVFillPayoff(/* Input */
						iNbPhi,
						iNbft,
						dLambdaX,
						dPhitMean,						
						iIndexPhiMean,
						dPhiStep,
						iIndexft0,
						dftStep,						
						dExTime,		/* Exercice of the swaption in years from today */
						dTStar,			/* Tstar */					 
						iNbCoupon,		/* Description of the cashflows */						
						CouponTime,
						Coupon,
						/*
						DFF_div_BTTstar,
						*/
						/* Outputs */
						PayOff);

	/* Save the PayOff */
	LGMSVSavePayOff( /* Inputs */
					iNbPhi,
					iNbft,

					 iIndexPhiMean,
					 iIndexft0,

					 dPhiStep,
					 dftStep,
					 dPhitMean,

					 PayOff);

	/* Calculation of the Integral of the payoff times the density */
	dIntegral = 0;
	for (iNumPhi=1;iNumPhi<=iNbPhi;iNumPhi++)
		for (iNumft=1;iNumft<=iNbft;iNumft++)
			/* Calculation of the payoff */
				dIntegral += Density[1][iNumPhi][iNumft]*PayOff[iNumPhi][iNumft];
	dIntegral *= dPhiStep*dftStep;

	/* Return the result */
	*Price = dIntegral;

	
FREE_RETURN:

	/* free memory */
	if (MomentFunc)
		free_LGMSVSolFunc3tensor(MomentFunc, 0, iNbCumulant,0, iNbCumulant,0, iNbCumulant);

	if (MomentValue)
		free_f3tensor(MomentValue, 0, iNbCumulant,0, iNbCumulant,0, iNbCumulant);
	
	if (IsUsefullMoment)
		free_l3tensor(IsUsefullMoment, 0, iNbCumulant,0, iNbCumulant,0, iNbCumulant);

	if (Combination)
		free_dvector(Combination,-1,iNbCumulant);

	if (CumulantRealPart) 
		free_dmatrix(CumulantRealPart,0, iNbCumulant,0, iNbCumulant);

	if (CumulantImagPart)
		free_dmatrix(CumulantImagPart,0, iNbCumulant,0, iNbCumulant);

	if (XRealPart)
		free_dmatrix(XRealPart,0, iNbCumulant,0, iNbCumulant);

	if (XImagPart)
		free_dmatrix(XImagPart,0, iNbCumulant,0, iNbCumulant);

	if (TempRealPart)
		free_dmatrix(TempRealPart,0, iNbCumulant,0, iNbCumulant);
	
	if (TempImagPart)
		free_dmatrix(TempImagPart,0, iNbCumulant,0, iNbCumulant);

	
	if (DFF_div_BTTstar)
		free_dvector(DFF_div_BTTstar,0,iNbCoupon-1);

	if (Density)
		free_f3tensor(Density, 1, 1,1, iNbPhi,1, iNbft);

	if (Densityspeq)
		free_dmatrix(Densityspeq, 1, 1,1, iNbPhi);

	if (PayOff)
		free_dmatrix(PayOff, 1, iNbPhi, 1,iNbft);
}

/* -------------------------------------------------------------------------------------------------------------
	LGMSVCalib	

  -------------------------------------------------------------------------------------------------------------- */
Err LGMSVCalib(
				char			*yc_name,						/*	Name of the yield curve */
				char			*ref_rate_name,					/*	Name of the reference rate */
				char			*swaption_freq,					/*	Frequency and basis of underlying swaptions */
				char			*swaption_basis,

				long			lExDate,
				long			lEndDate,
				double			dStrike,

				double			dTau,
				double			dg,
				double			dAlpha,
				double			dRho,
				double			dLambdaEps,
				long			iNbCumulant,
				long			iNbPhi,
				long			iNbft,

				double			iNbSigmaPhiGridLeft,
				double			iNbSigmaPhiGridRight,
				double			iNbSigmaXGridLeft,
				double			iNbSigmaXGridRight,

				/* Output */
				double			*pSwaptionPrice)
{
	int				i,   ncpn;
	SrtCompounding	ifreq;
	SrtBasisCode	ibasis;
	long			cpn_date[MAX_CPN];
	double			cpn_time[MAX_CPN],
					cpn_cvg[MAX_CPN];
				
	long			theo_date, act_date, temp_date;
	long			today; 


	SrtCurvePtr		yc_ptr;
	Err				err				= NULL;

	double SigTime = 0;

	if (iNbSigmaPhiGridLeft == 0)
		iNbSigmaPhiGridLeft = LGMSV_iNbSigmaExtremePoints;

	if (iNbSigmaPhiGridRight == 0)
		iNbSigmaPhiGridRight = LGMSV_iNbSigmaExtremePoints;

	yc_ptr = lookup_curve (yc_name);
	if (!yc_ptr)
	{
		err = "Yield Curve not found";
		goto FREE_RETURN;
	}
	today = get_today_from_curve (yc_ptr);	

	/*	1.)	Setup the bond schedule and its coupons */
	
	/*	Coupons */

	err = interp_compounding (swaption_freq, &ifreq);
	if (err)
	{
		goto FREE_RETURN;
	}

	err = interp_basis (swaption_basis, &ibasis);
	if (err)
	{
		goto FREE_RETURN;
	}

	theo_date = lEndDate;
	act_date = bus_date_method (theo_date, MODIFIED_SUCCEEDING);
	ncpn = 1;

	while (act_date > lExDate)
	{
		theo_date = add_unit (theo_date, - 12 / ifreq, SRT_MONTH, NO_BUSDAY_CONVENTION);
		act_date = bus_date_method (theo_date, MODIFIED_SUCCEEDING);
		ncpn++;
	}
	ncpn--;

	if (ncpn < 2)
	{
		err = "Not enough coupons";
		goto FREE_RETURN;		
	}

	theo_date = lEndDate;
	act_date = bus_date_method (theo_date, MODIFIED_SUCCEEDING);
	i = ncpn - 1;

	while (i >= 0)
	{
		cpn_time[i] = (act_date - today) * YEARS_IN_DAY; 
		cpn_date[i] = act_date;

		theo_date = add_unit (theo_date, - 12 / ifreq, SRT_MONTH, NO_BUSDAY_CONVENTION);

		temp_date = bus_date_method (theo_date, MODIFIED_SUCCEEDING);
		cpn_cvg[i] = coverage (temp_date, act_date, ibasis);
		act_date = temp_date;

		i--;
	}
	cpn_cvg[0] = 0.0;


	
	/* Calculation of the iNbSigmaXGridRight */
	/* SHOULD BE A SHARED FUNCTION           */
	if (iNbSigmaXGridRight==0)
		if (fabs(dLambdaEps)<1E-06)
		{
			iNbSigmaXGridRight = 7.5+25*max(dRho,0)*dAlpha ;
		}
		else
		{
			iNbSigmaXGridRight = 7.5+25*max(dRho,0)*dAlpha*
					sqrt((1-exp(-2*dLambdaEps*cpn_time[ncpn - 1]))/(2*dLambdaEps*cpn_time[ncpn - 1]));
		}

	if (iNbSigmaXGridLeft == 0)
		iNbSigmaXGridLeft = LGMSV_iNbSigmaXGridLeft;
	

	/* Fill the Coupon */

	/* -K*Lvl */
	for (i=1;i<ncpn;i++)
	{
		cpn_cvg[i] *= -dStrike;
	}
	/* +S*LVL = (if S=SwapCash Rate) Df(Start)-Df(end) */
	cpn_cvg[0] +=1;
	cpn_cvg[ncpn - 1] -=1;


	
	LGMSVClosedForm(/* Inputs */
					 iNbCumulant,	/* Number of cumulant for the approximation of the (f,Phi) density */
					 iNbPhi,			/* Number of Phi : Should be a power of two */
					 iNbft,			/* Number of ft : Should be a power of two */	

					 /* Parameter of diffusion */
					 /* BE CAREFULL : Alpha and LambdaEps for V=Eps^2 */
					 1.0/dTau,
					 
					 2*dAlpha,		/* Alpha of V = Eps^2 */
					 2*dLambdaEps,	/* LambdaEps of V = Eps^2 */

					 dRho,

					 /* Parameter of grids */
					 iNbSigmaPhiGridLeft,
					 iNbSigmaPhiGridRight,


					 iNbSigmaXGridLeft,
					 iNbSigmaXGridRight,


					 lExDate,							/* Exercice date of the swaption  */
					 (lExDate-today)*YEARS_IN_DAY,		/* Exercice of the swaption in years from today */
					 
					 LGMSV_Tstar,			/* Tstar in years from today */
					 
					 1,		/* Term Structure of g(t) */
					 &SigTime,		
					 &dg,

					 ncpn,		/* Description of the cashflows */
					 cpn_time,
					 cpn_date,
					 cpn_cvg,

					 yc_name,	/* Yield Curve */

					 /* Outputs */
					 pSwaptionPrice);


FREE_RETURN:

	if (err)
	{
	
	}

	return err;
}





static void solve_poly2_comp(double	Are,
							double	Aim,
							double	Bre,
							double	Bim,
							double	Cre,
							double	Cim,
							double	*x1re,
							double	*x1im,
							double	*x2re,
							double	*x2im)
{
static double DeltaRe, DeltaIm, DeltaRe2, DeltaIm2;

	prod_comp(Bre, Bim, Bre, Bim, &DeltaRe, &DeltaIm);
	prod_comp(Are, Aim, Cre, Cim, &DeltaRe2, &DeltaIm2);
	DeltaRe -= 4.0 * DeltaRe2;
	DeltaIm -= 4.0 * DeltaIm2;

	sqr_comp(DeltaRe, DeltaIm, &DeltaRe, &DeltaIm);

	*x1re = -Bre + DeltaRe;
	*x1im = -Bim + DeltaIm;
	*x2re = -Bre - DeltaRe;
	*x2im = -Bim - DeltaIm;

	div_comp(*x1re, *x1im, 2.0 * Are, 2.0 * Aim, x1re, x1im);
	div_comp(*x2re, *x2im, 2.0 * Are, 2.0 * Aim, x2re, x2im);
}

static void LGMSVCalculateReClosedFormPDE(	
											/* Parameters */
											double	lambdaX,
											double	sig0,											
											double	alpha,
											double	lambdaEps,
											double	rho,
											double	TStar,			/* Tstar in years from today */

											double	T,
											double	dPhifreq,
											double	dftFreq,

											int		NbPoint,

											double	*ResARe,
											double	*ResAIm,

											double	*ResDRe,
											double	*ResDIm
											)
{
double	Are, Aim, Bre, Bim, Cre, Cim, DRe, DIm, ARe, AIm, DReTemp, CreTemp, CimTemp, BimTemp;
double	alpha2 = alpha * alpha;
double	dt, t, expt;
int		i;

	dt = T / NbPoint;

	Are = -0.5 * alpha2 * dt;
	Aim = 0.0;
	Bre = lambdaEps * dt;
	Bim = -alpha * rho * sig0 * dftFreq * exp(-lambdaX * TStar) * dt;
	Cre = 0.5 * sig0 * sig0  * dftFreq * dftFreq * exp(-2.0 * lambdaX * TStar) * dt;
	Cim = -sig0 * sig0 * dPhifreq * exp(-2.0 * lambdaX * T) * dt;

	ARe = 0.0;
	AIm = 0.0;

	DRe = 0.0;
	DIm = 0.0;

	t = T;

	for (i=0; i<NbPoint; i++)
	{
		ARe += Bre * DRe;
		AIm += Bre * DIm;

		DReTemp = DRe;

		expt = exp(lambdaX * t);

		BimTemp = Bim * expt;
		expt *= expt;
		CreTemp = Cre * expt;
		CimTemp = Cim * expt;

		DRe -= Are * (DRe * DRe - DIm * DIm) + Bre * DRe - BimTemp * DIm + CreTemp;
		DIm -= 2.0  * Are * DReTemp * DIm + Bre * DIm + BimTemp * DRe + CimTemp;

		t -= dt;
	}

	*ResARe = ARe;
	*ResAIm = AIm;

	*ResDRe = DRe;
	*ResDIm = DIm;
}

void LGMSVCalculateReClosedForm(double	dPhifreq,
								double	dftFreq,

								/* Parameters */
								double	sig0,
								double	alpha,
								double	rho,
								double	T,

								double	*CaracRe)
{
double	Are, Aim, Bre, Bim, Cre, Cim;
double	x1Re, x1Im, x2Re, x2Im;
double	TempRe, TempIm, DeltaRe, DeltaIm;
double	alpha2 = alpha * alpha;
double	dpi = 3.1415926535897932384626433832795;

	Are = -0.5 * alpha2;
	Aim = 0.0;
	Bre = 0.0;
	Bim = -alpha * rho * sig0 * dftFreq;
	Cre = 0.5 * sig0 * sig0  * dftFreq * dftFreq;
	Cim = -sig0 * sig0 * dPhifreq;

	solve_poly2_comp(Are, Aim, Bre, Bim, Cre, Cim, &x1Re, &x1Im, &x2Re, &x2Im);

	DeltaRe = Are * (x1Re - x2Re);
	DeltaIm = Are * (x1Im - x2Im);

	if (fabs(DeltaRe) + fabs(DeltaIm) > 1.0E-12)
	{
		TempRe = exp(-DeltaRe * T);
		TempIm = fmod(-DeltaIm * T, PI2);
		
		DeltaRe = TempRe * cos(TempIm);
		DeltaIm = TempRe * sin(TempIm);

		div_comp(x1Re, x1Im, x2Re, x2Im, &TempRe, &TempIm);
		prod_comp(TempRe, TempIm, DeltaRe, DeltaIm, &TempRe, &TempIm);
		div_comp(1.0 - DeltaRe, -DeltaIm, 1.0 - TempRe, -TempIm, &TempRe, &TempIm);
		prod_comp(x1Re, x1Im, TempRe, TempIm, &TempRe, &TempIm);
	}
	else
	{
		div_comp(x1Re * T, x1Im * T, x1Re * T + 2.0 / alpha2, x1Im, &TempRe, &TempIm);
		prod_comp(TempRe, TempIm, x1Re, x1Im, &TempRe, &TempIm);
	}

	*CaracRe = TempRe;
}

void LGMSVCalculateDensityClosedForm(
							/* Input */
							int iNbCumulant,
							int iNbPhi,
							int iNbft,
							
							double dPhiFreqStep,
							double dftFreqStep,

							/* Parameters */
							double	sig0,
							double	alpha,
							double	rho,
							double	T,
							
							/* Outputs */
							double ***Density,
							double **DensitySpeq)
{
double	dpi;
int		iNumPhiFreq, iNumftFreq;
double	dPhifreq;
double	dCoef;

double	dftFreq;
double	dExpReal, dImag;

double	x1Re, x1Im, x2Re, x2Im,
		DeltaRe, DeltaIm,
		TempRe, TempIm;
double	alpha2, rho22, sig0alpha, rhosig0alpha, sig02, sig0alpha2, PhiMean;

double	fact;
double	X, Y;
double	DriftPhi;
int		nbFtHalf;

	/* Initialisation */
	Density = Density;
	DensitySpeq = DensitySpeq;

	/* Fill the TF of the density */
	dpi = 3.1415926535897932384626433832795;	

	rho22 = (1.0 - rho * rho) / 2.0;
	alpha2 = alpha * alpha;
	sig0alpha = sig0 * alpha;
	sig0alpha2 = sig0alpha * sig0alpha;
	rhosig0alpha = rho * sig0alpha;
	sig02 = sig0 * sig0;

	PhiMean = sig02 * T;

	nbFtHalf = (int) (iNbft / 2.0) + 1;

	/* Normalisation such that rlft3 gives us the true ifft of TF(P(f_phi,f_ft))) */
	dCoef = dftFreqStep*dPhiFreqStep*2.0;/*(iNbPhi*iNbft/2.0); */
	
	for (iNumPhiFreq=1;iNumPhiFreq<=iNbPhi;iNumPhiFreq++)
	{	
		dPhifreq = PI2 * (iNumPhiFreq-1-iNbPhi*(int)((iNumPhiFreq-1)/(iNbPhi/2.0)))*dPhiFreqStep;
		Y = sig0alpha2 * dPhifreq;
		DriftPhi = dPhifreq * PhiMean;

		if (dPhifreq > -1.0E-16)
		{
			fact = 1.0;
		}
		else
		{
			fact = -1.0;
		}

		for (iNumftFreq=1; iNumftFreq<=nbFtHalf; iNumftFreq++)
		{
			dftFreq = PI2 * dftFreqStep*(iNumftFreq-1);			
			X = sig0alpha * dftFreq;

			TempRe = X * X * rho22;
			TempIm = sqrt(Y * Y + TempRe * TempRe);

			DeltaRe = sqrt(TempRe + TempIm);
			DeltaIm = -fact * sqrt(-TempRe + TempIm);			

			x1Re = -DeltaRe / alpha2;
			x1Im = -(DeltaIm + rho * X) / alpha2;

			if (fabs(DeltaRe) + fabs(DeltaIm) > 0.0)
			{				
				x2Re = DeltaRe / alpha2;
				x2Im = (DeltaIm - rho * X) / alpha2;

				TempRe = exp(-DeltaRe * T);
				TempIm = fmod(-DeltaIm * T, PI2);
				
				DeltaRe = TempRe * cos(TempIm);
				DeltaIm = TempRe * sin(TempIm);

				div_comp(x1Re, x1Im, x2Re, x2Im, &TempRe, &TempIm);
				prod_comp(TempRe, TempIm, DeltaRe, DeltaIm, &TempRe, &TempIm);
				div_comp(1.0 - DeltaRe, -DeltaIm, 1.0 - TempRe, -TempIm, &TempRe, &TempIm);
				prod_comp(x1Re, x1Im, TempRe, TempIm, &TempRe, &TempIm);				
			}
			else
			{				
				div_comp(x1Re * T, x1Im * T, x1Re * T + 2.0 / alpha2, x1Im, &TempRe, &TempIm);
				prod_comp(TempRe, TempIm, x1Re, x1Im, &TempRe, &TempIm);
			}

			LGMSVCalculateReClosedForm(dPhifreq, dftFreq, sig0, alpha, rho, T, &DeltaRe);
			
			if (iNumftFreq < nbFtHalf)
			{
				dExpReal = dCoef * exp(TempRe);
				dImag = fmod(TempIm - DriftPhi, PI2);

				/* Real Part */
				Density[1][iNumPhiFreq][2*iNumftFreq-1] = dExpReal * cos(dImag);		
				/* Imaginary Part */
				Density[1][iNumPhiFreq][2*iNumftFreq] = dExpReal * sin(dImag);
			}
			else
			{
				DensitySpeq[1][2*iNumPhiFreq-1] = dCoef * exp(TempRe);
			}
		}
	}	
	
	/* Save the TF density */
	LGMSVSaveTFDensity( /* Inputs */
						iNbPhi,
						iNbft,					 
						Density);

	/* Calculation of the density by iFFT */
	rlft3(Density,DensitySpeq,1,iNbPhi,iNbft,-1);
}

void LGMSVCalculateDensityPDE(
							/* Input */
							int iNbCumulant,
							int iNbPhi,
							int iNbft,
							
							double dPhiFreqStep,
							double dftFreqStep,

							/* Parameters */
							double	lambdaX,
							double	sig0,											
							double	alpha,
							double	lambdaEps,
							double	rho,
							double	TStar,			/* Tstar in years from today */

							double	T,

							int		NbPoint,
							
							/* Outputs */
							double ***Density,
							double **DensitySpeq)
{
double	dpi;
int		iNumPhiFreq, iNumftFreq;
double	dPhifreq;
double	dCoef;

double	dftFreq;
double	dExpReal, dImag;

double	ARe, AIm, DRe, DIm,
		TempRe, TempIm;
double	PhiMean;
double	constPhi;

double	DriftPhi;
int		nbFtHalf;


	/* Fill the TF of the density */
	dpi = 3.1415926535897932384626433832795;	

	PhiMean = sig0 * sig0 * (1.0 - exp(-2.0 * lambdaX * T)) / (2.0 * lambdaX);

	constPhi = exp(-2.0 * lambdaX * T);

	nbFtHalf = (int) (iNbft / 2.0) + 1;

	/* Normalisation such that rlft3 gives us the true ifft of TF(P(f_phi,f_ft))) */
	dCoef = dftFreqStep*dPhiFreqStep*2.0;/*(iNbPhi*iNbft/2.0); */
	
	for (iNumPhiFreq=1;iNumPhiFreq<=iNbPhi;iNumPhiFreq++)
	{	
		dPhifreq = PI2 * (iNumPhiFreq-1-iNbPhi*(int)((iNumPhiFreq-1)/(iNbPhi/2.0)))*dPhiFreqStep;
		DriftPhi = dPhifreq * PhiMean;
		
		for (iNumftFreq=1; iNumftFreq<=nbFtHalf; iNumftFreq++)
		{
			dftFreq = PI2 * dftFreqStep*(iNumftFreq-1);						

			LGMSVCalculateReClosedFormPDE(	lambdaX, sig0, alpha, lambdaEps, rho, TStar, T,
											dPhifreq, dftFreq, NbPoint, &ARe, &AIm, &DRe, &DIm);

			TempRe = ARe + DRe;
			TempIm = AIm - DriftPhi + DIm;
			
			if (iNumftFreq < nbFtHalf)
			{
				dExpReal = dCoef * exp(TempRe);
				dImag = fmod(TempIm, PI2);

				/* Real Part */
				Density[1][iNumPhiFreq][2*iNumftFreq-1] = dExpReal * cos(dImag);		
				/* Imaginary Part */
				Density[1][iNumPhiFreq][2*iNumftFreq] = dExpReal * sin(dImag);
			}
			else
			{
				DensitySpeq[1][2*iNumPhiFreq-1] = dCoef * exp(TempRe);
			}
		}
	}	
	
	/* Save the TF density */
	LGMSVSaveTFDensity( /* Inputs */
						iNbPhi,
						iNbft,					 
						Density);

	/* Calculation of the density by iFFT */
	rlft3(Density,DensitySpeq,1,iNbPhi,iNbft,-1);
}