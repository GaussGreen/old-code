#include        <utallhdr.h"
#include        <opHeston.h"
#include        "opfnctns.h"
#include		<math.h"


#define R 0.61803399
#define C (1.0-R)
#define SHFT2(a,b,c) (a)=(b);(b)=(c);
#define SHFT3(a,b,c,d) (a)=(b);(b)=(c);(c)=(d);

Err HestonPrice(	
					double	Forward,
					double* Strike,
					int		nStrikes,
					double	Maturity,
					double	Sigma,
					double	Alpha,
					double  SigmaIfty, 
					double	b,
					double  Gamma,
					double	Rho,
					double	Disc,
					double  UpperBound,
					SrtCallPutType 	call_put, 
					SrtGreekType greek,
					SRT_Boolean isVolInfFix,
					int IntegrType,
					int nSteps,
					double 	*result
					)
{


	double a,CosT1,SinT1,CosT2,SinT2;
	double *x,*w,*x1,*w1,G11,*G12,G21,*G22,*GKK,n1,n2,n3,RePsi1,ImPsi1,RePsi2,ImPsi2;
	double DRePsi1,DRePsi2,DImPsi1,DImPsi2,DDRePsi1,DDRePsi2,DDImPsi1,DDImPsi2,*premium;
	double esp,esp2,ImAppr1,ImAppr2,dStdDev,dApproxBSVol,nStdDev=5.0;
	int i,j,k,nStp=5,iter,iThresh,nOptSteps1,nOptSteps2,nSteps_opt;
	double Int,Int1,NodeLeft,NodeRight,NodeLeft1,NodeRight1,dScaling=1.,dGear;
	Err err = NULL;

//	IntegrType=3;

	G12=dvector(1,nStrikes);
	G22=dvector(1,nStrikes);
	GKK=dvector(1,nStrikes);
	premium=dvector(1,nStrikes); 

	/* here I define the max n. of Std. Dev. */

	dApproxBSVol= (Forward+Gamma)/Forward*Sigma;
	dStdDev = Forward*sqrt(exp(dApproxBSVol*dApproxBSVol*Maturity)-1.);

	/* if Infinity Vol is fixed, then a is set to the initial vol */
	if (isVolInfFix) SigmaIfty=Sigma;
	a=SigmaIfty*SigmaIfty*b;

	esp=exp(-b*Maturity);
	esp2=exp(-Maturity*(b-Alpha*Rho));


	/* Here I redefinee Strikes and Forward (to improve the convergence) by a scaling factor */

//	Forward+=Gamma;
	Forward=(Forward+Gamma)*dScaling;

	for (j=1; j<=nStrikes;j++)
	{
		Strike[j] = (Strike[j]+Gamma)*dScaling;

//		Strike[j] += Gamma;
	}

	switch (greek) {

	case PREMIUM:

		switch (IntegrType){

		case 1:

		x=dvector(1,nStp);
		w=dvector(1,nStp);
		x1=dvector(1,nStp);
		w1=dvector(1,nStp);

		for (j=1;j<=nStrikes;j++) {
			G12[j]=0.0;
			G22[j]=0.0;
		}

			ImAppr2=-(log(Forward)+(esp2-1)*Sigma*Sigma/2./(b-Alpha*Rho)+Maturity*a/2./(b-Alpha*Rho)
				                +a*(esp2-1.)/2/(b-Alpha*Rho)/(b-Alpha*Rho));
			ImAppr1=-(log(Forward)+Sigma*Sigma/b*(1.-esp)-Maturity*a/2./b+a/2./b/b*(1.-esp));


		for (j=1;j<=nStrikes;j++) {

			iter=0;
			NodeLeft=0.0;
			NodeLeft1=0.0;
			do  {

			/* find the next node on the axis */

			NodeRight = HestonFindRightNode(ImAppr1,Strike[j],NodeLeft,isVolInfFix,Maturity,
											a,b,Alpha,Rho,Forward,Sigma,1.0,greek);
			GaussLeg(NodeLeft, NodeRight, x, w, nStp);

			NodeRight1 = HestonFindRightNode(ImAppr2,Strike[j],NodeLeft,isVolInfFix,Maturity,
											a,b,Alpha,Rho,Forward,Sigma,0.0,greek);
			GaussLeg(NodeLeft1, NodeRight1, x1, w1, nStp);

				Int=0.0;
				Int1=0.0;
				for (k=1;k<=nStp;k++) {

					err=PsiFunction (isVolInfFix,Maturity,1.,-x[k],a,b,Alpha,Rho,log(Forward),
						             Sigma*Sigma,greek,&RePsi1,&ImPsi1,
									 &DRePsi1,&DImPsi1,&DDRePsi1,&DDImPsi1);

					G12[j]+=w[k]/x[k]*exp(RePsi1)*sin(ImPsi1+x[k]*log(Strike[j])); 

					err=PsiFunction (isVolInfFix,Maturity,0.,-x1[k],a,b,Alpha,Rho,log(Forward),
									 Sigma*Sigma,greek,&RePsi2,&ImPsi2,
									 &DRePsi2,&DImPsi2,&DDRePsi2,&DDImPsi2);

					G22[j]+=w1[k]/x1[k]*exp(RePsi2)*sin(ImPsi2+x1[k]*log(Strike[j]));	
			
				}

			NodeLeft=NodeRight;
			NodeLeft1=NodeRight1;

			iter++;

			} while ((NodeLeft <= UpperBound) && (iter <= 1000));
		}

			/* computes the remaining two terms that need no integration */

			err = PsiFunction (isVolInfFix,Maturity,1.,0.,a,b,Alpha,Rho,log(Forward),Sigma*Sigma,greek,&RePsi1,&ImPsi1,&DRePsi1,&DImPsi1,&DDRePsi1,&DDImPsi1);
			G11=exp(RePsi1);

			err=PsiFunction (isVolInfFix,Maturity,0.,0.,a,b,Alpha,Rho,log(Forward),Sigma*Sigma,greek,&RePsi2,&ImPsi2,&DRePsi2,&DImPsi2,&DDRePsi2,&DDImPsi2);
			G21=exp(RePsi2);

			/* finally puts all together and evaluates the option */

			for (j=1;j<=nStrikes;j++) {
				
				n1=G11/2.-1./SRT_PI*G12[j];
				n2=G21/2.-1./SRT_PI*G22[j];

				switch (call_put){

					case SRT_CALL:
					result[j]=DMAX(1.e-13,(n1-Strike[j]*n2)*Disc);

					if (result[j]< (Forward-Strike[j])*Disc) result[j]=(Forward-Strike[j])*Disc+1.e-8;

				break;

				default:
					result[j]= DMAX(1.e-13,(n1-Strike[j]*n2)*Disc - (Forward-Strike[j])*Disc);

					if (result[j]< (-Forward+Strike[j])*Disc) result[j]=(-Forward+Strike[j])*Disc+1.e-8;	
				break;

				}
			}

			free_dvector(x,1,nStp);
			free_dvector(w,1,nStp);
			free_dvector(x1,1,nStp);
			free_dvector(w1,1,nStp);

		break;


		case 0:

		x=dvector(1,nSteps);
		w=dvector(1,nSteps);
		GaussLeg(0., UpperBound, x, w, nSteps);

		for (j=1;j<=nStrikes;j++) {
			G12[j]=0.0;
			G22[j]=0.0;
		}

/*
			// computes the approximated Re(Psi)/v
			
			PsiFunction (isVolInfFix,Maturity,1.,-UpperBound,a,b,Alpha,Rho,log(Forward),Sigma*Sigma,greek,&ReUpBound,&ImUpBound,&DRePsi1,&DImPsi1,&DDRePsi1,&DDImPsi1);
			PsiFunction (isVolInfFix,Maturity,1.,0.0,a,b,Alpha,Rho,log(Forward),Sigma*Sigma,greek,&ReLowBound,&ImLowBound,&DRePsi1,&DImPsi1,&DDRePsi1,&DDImPsi1);

			ReAppr1=(ReUpBound-ReLowBound)/(UpperBound-x[2]);
			ReAppr1c=ReLowBound;

			PsiFunction (isVolInfFix,Maturity,0.,-UpperBound,a,b,Alpha,Rho,log(Forward),Sigma*Sigma,greek,&ReUpBound,&ImUpBound,&DRePsi1,&DImPsi1,&DDRePsi1,&DDImPsi1);
			PsiFunction (isVolInfFix,Maturity,0.,0.0,a,b,Alpha,Rho,log(Forward),Sigma*Sigma,greek,&ReLowBound,&ImLowBound,&DRePsi1,&DImPsi1,&DDRePsi1,&DDImPsi1);

			ReAppr2=(ReUpBound-ReLowBound)/(UpperBound-x[2]);
			ReAppr2c=ReLowBound;

			// computes the approximated Im(Psi)/v 

			ImAppr2=-(log(Forward)+(esp2-1)*Sigma*Sigma/2./(b-Alpha*Rho)+Maturity*a/2./(b-Alpha*Rho)
				                +a*(esp2-1.)/2/(b-Alpha*Rho)/(b-Alpha*Rho));
			ImAppr1=-(log(Forward)+Sigma*Sigma/b*(1.-esp)-Maturity*a/2./b+a/2./b/b*(1.-esp));

*/
		for (i=1;i<=nSteps;i++) {

			err=PsiFunction (isVolInfFix,Maturity,1.,-x[i],a,b,Alpha,Rho,log(Forward),Sigma*Sigma,greek,&RePsi1,&ImPsi1,&DRePsi1,&DImPsi1,&DDRePsi1,&DDImPsi1);
			err=PsiFunction (isVolInfFix,Maturity,0.,-x[i],a,b,Alpha,Rho,log(Forward),Sigma*Sigma,greek,&RePsi2,&ImPsi2,&DRePsi2,&DImPsi2,&DDRePsi2,&DDImPsi2);

			// runs over the strikes 

			for (j=1;j<=nStrikes;j++) {


				G12[j]+=w[i]/x[i]*exp(RePsi1)*sin(ImPsi1+x[i]*log(Strike[j]));
				G22[j]+=w[i]/x[i]*exp(RePsi2)*sin(ImPsi2+x[i]*log(Strike[j]));

/*
				G12[j]+=w[i]/x[i]*exp(RePsi1)*sin(ImPsi1+x[i]*log(Strike[j])) -
					    w[i]/x[i]*exp(ReAppr1c+ReAppr1*x[i]) * sin(x[i]*log(Strike[j])+x[i]*ImAppr1);

				G22[j]+=w[i]/x[i]*exp(RePsi2)*sin(ImPsi2+x[i]*log(Strike[j]))-
	 		            w[i]/x[i]*exp(ReAppr2c+ReAppr2*x[i]) * sin(x[i]*log(Strike[j])+x[i]*ImAppr2);
*/

			}
		}

			err = PsiFunction (isVolInfFix,Maturity,1.,0.,a,b,Alpha,Rho,log(Forward),Sigma*Sigma,greek,&RePsi1,&ImPsi1,&DRePsi1,&DImPsi1,&DDRePsi1,&DDImPsi1);
			G11=exp(RePsi1);

			err=PsiFunction (isVolInfFix,Maturity,0.,0.,a,b,Alpha,Rho,log(Forward),Sigma*Sigma,greek,&RePsi2,&ImPsi2,&DRePsi2,&DImPsi2,&DDRePsi2,&DDImPsi2);
			G21=exp(RePsi2);

			for (j=1;j<=nStrikes;j++) {

/*
				ContrVar = exp(ReAppr1c)*atan(-(log(Strike[j])+ImAppr1)/ReAppr1);
				n1=G11/2.-1./SRT_PI*(G12[j]+ContrVar);

				ContrVar = exp(ReAppr2c)*atan(-(log(Strike[j])+ImAppr2)/ReAppr2);
				n2=G21/2.-1./SRT_PI*(G22[j]+ContrVar);
*/				
				n1=G11/2.-1./SRT_PI*G12[j];
				n2=G21/2.-1./SRT_PI*G22[j];

				switch (call_put){

					case SRT_CALL:
					result[j]=DMAX(1.e-13,(n1-Strike[j]*n2)*Disc)/dScaling;

					if (result[j]< (Forward-Strike[j])*Disc) result[j]=(Forward-Strike[j])*Disc/dScaling+1.e-8;

				break;

				default:
					result[j]= DMAX(1.e-13,(n1-Strike[j]*n2)*Disc/dScaling - (Forward-Strike[j])*Disc)/dScaling;

					if (result[j]< (-Forward+Strike[j])*Disc) result[j]=(-Forward+Strike[j])*Disc+1.e-8;	
				break;

				}
			}

			free_dvector(x,1,nSteps);
			free_dvector(w,1,nSteps);

			break;

		case 3:


		for (j=1;j<=nStrikes;j++) {
			G12[j]=0.0;
			G22[j]=0.0;
		}


			ImAppr2=-(log(Forward)+(esp2-1)*Sigma*Sigma/2./(b-Alpha*Rho)+Maturity*a/2./(b-Alpha*Rho)
				                +a*(esp2-1.)/2/(b-Alpha*Rho)/(b-Alpha*Rho));
			ImAppr1=-(log(Forward)+Sigma*Sigma/b*(1.-esp)-Maturity*a/2./b+a/2./b/b*(1.-esp));


	UpperBound=60.+ 20*pow(DMAX(0.0,5.-Maturity),1.4);


	for (j=1;j<=nStrikes;j++) {	

			// computes the approximated n. of oscillations in a unit interval 
	        // and then the required n. integration steps. For very short maturities, 
		    // the precision is forced to increse, as there is little time value


		// Here it branches following the rule: if Strike is larger than nStdDev from the forward 
		// or smaller than -nStdDev, then it simply evaluates the price at +- nStdDev
		

		dGear = 2 + 3./(Maturity*Maturity*Maturity);
		nOptSteps1 = IMAX(IMIN((int)( dGear*UpperBound*(fabs(ImAppr1+log(Strike[j])))/SRT_PI),nSteps),22);
		nOptSteps2 = IMAX(IMIN((int)( dGear*UpperBound*(fabs(ImAppr2+log(Strike[j])))/SRT_PI),nSteps),22);

		nSteps_opt=IMAX(nOptSteps1,nOptSteps2);
		x=dvector(1,nSteps_opt);
		w=dvector(1,nSteps_opt);

		GaussLeg(0., UpperBound, x, w, nSteps_opt);

		Strike[j] = DMIN(Strike[j],Forward + nStdDev*dStdDev);
		Strike[j] = DMAX(Strike[j], DMAX(0.0015,Forward - nStdDev*dStdDev));

		for (i=1;i<=nSteps_opt;i++) {

			err=PsiFunction (isVolInfFix,Maturity,1.,-x[i],a,b,Alpha,Rho,log(Forward),Sigma*Sigma,greek,&RePsi1,&ImPsi1,&DRePsi1,&DImPsi1,&DDRePsi1,&DDImPsi1);
			err=PsiFunction (isVolInfFix,Maturity,0.,-x[i],a,b,Alpha,Rho,log(Forward),Sigma*Sigma,greek,&RePsi2,&ImPsi2,&DRePsi2,&DImPsi2,&DDRePsi2,&DDImPsi2);

				G12[j]+=w[i]/x[i]*exp(RePsi1)*sin(ImPsi1+x[i]*log(Strike[j]));
				G22[j]+=w[i]/x[i]*exp(RePsi2)*sin(ImPsi2+x[i]*log(Strike[j]));
			
		}

			err = PsiFunction (isVolInfFix,Maturity,1.,0.,a,b,Alpha,Rho,log(Forward),Sigma*Sigma,greek,&RePsi1,&ImPsi1,&DRePsi1,&DImPsi1,&DDRePsi1,&DDImPsi1);
			G11=exp(RePsi1);

			err=PsiFunction (isVolInfFix,Maturity,0.,0.,a,b,Alpha,Rho,log(Forward),Sigma*Sigma,greek,&RePsi2,&ImPsi2,&DRePsi2,&DImPsi2,&DDRePsi2,&DDImPsi2);
			G21=exp(RePsi2);

				
				n1=G11/2.-1./SRT_PI*G12[j];
				n2=G21/2.-1./SRT_PI*G22[j];

				switch (call_put){

					case SRT_CALL:
					result[j]=DMAX(1.e-13,(n1-Strike[j]*n2)*Disc)/dScaling;

					if (result[j]< (Forward-Strike[j])*Disc) result[j]=(Forward-Strike[j])*Disc/dScaling+1.e-8;

				break;

				default:
					result[j]= DMAX(1.e-13,(n1-Strike[j]*n2)*Disc/dScaling - (Forward-Strike[j])*Disc)/dScaling;

					if (result[j]< (-Forward+Strike[j])*Disc) result[j]=(-Forward+Strike[j])*Disc+1.e-8;	
				break;

				}

			free_dvector(x,1,nSteps_opt);
			free_dvector(w,1,nSteps_opt);
	}		


			break;


		default:

		/* runs over strikes to set the breakpoint between the two methods */

		iThresh=nStrikes;
		for (j=0;j<nStrikes;j++) {
			if (Strike[j+1] > 0.003*dScaling)  {
				iThresh = j;
				break;
			}
		}

		/* now I integrate with the optimised algorithm for strikes < breakpoint and with the parallel
		   algorithm for strikes  > breakpoint */

		x=dvector(1,nStp);
		w=dvector(1,nStp);
		x1=dvector(1,nStp);
		w1=dvector(1,nStp);

		for (j=1;j<=iThresh;j++) {
			G12[j]=0.0;
			G22[j]=0.0;
		}

			ImAppr2=-(log(Forward)+(esp2-1)*Sigma*Sigma/2./(b-Alpha*Rho)+Maturity*a/2./(b-Alpha*Rho)
				                +a*(esp2-1.)/2/(b-Alpha*Rho)/(b-Alpha*Rho));
			ImAppr1=-(log(Forward)+Sigma*Sigma/b*(1.-esp)-Maturity*a/2./b+a/2./b/b*(1.-esp));


		for (j=1;j<=iThresh;j++) {

			iter=0;
			NodeLeft=0.0;
			NodeLeft1=0.0;
			do  {

			/* find the next node on the axis */

			NodeRight = HestonFindRightNode(ImAppr1,Strike[j],NodeLeft,isVolInfFix,Maturity,
											a,b,Alpha,Rho,Forward,Sigma,1.0,greek);
			GaussLeg(NodeLeft, NodeRight, x, w, nStp);

			NodeRight1 = HestonFindRightNode(ImAppr2,Strike[j],NodeLeft,isVolInfFix,Maturity,
											a,b,Alpha,Rho,Forward,Sigma,0.0,greek);
			GaussLeg(NodeLeft1, NodeRight1, x1, w1, nStp);

				Int=0.0;
				Int1=0.0;
				for (k=1;k<=nStp;k++) {

					err=PsiFunction (isVolInfFix,Maturity,1.,-x[k],a,b,Alpha,Rho,log(Forward),
						             Sigma*Sigma,greek,&RePsi1,&ImPsi1,
									 &DRePsi1,&DImPsi1,&DDRePsi1,&DDImPsi1);

					G12[j]+=w[k]/x[k]*exp(RePsi1)*sin(ImPsi1+x[k]*log(Strike[j])); 

					err=PsiFunction (isVolInfFix,Maturity,0.,-x1[k],a,b,Alpha,Rho,log(Forward),
									 Sigma*Sigma,greek,&RePsi2,&ImPsi2,
									 &DRePsi2,&DImPsi2,&DDRePsi2,&DDImPsi2);

					G22[j]+=w1[k]/x1[k]*exp(RePsi2)*sin(ImPsi2+x1[k]*log(Strike[j]));	
			
				}

			NodeLeft=NodeRight;
			NodeLeft1=NodeRight1;

			iter++;

			} while ((NodeLeft <= UpperBound) && (iter <= 1000));
		}

			free_dvector(x,1,nStp);
			free_dvector(w,1,nStp);
			free_dvector(x1,1,nStp);
			free_dvector(w1,1,nStp);

			/* second part, now I use the parallel integrator for strikes > 30 bps */

			x=dvector(1,nSteps);
			w=dvector(1,nSteps);
			GaussLeg(0., UpperBound, x, w, nSteps);

		for (j=iThresh+1;j<=nStrikes ;j++) {
			G12[j]=0.0;
			G22[j]=0.0;
		}


		for (i=1;i<=nSteps;i++) {

			err=PsiFunction (isVolInfFix,Maturity,1.,-x[i],a,b,Alpha,Rho,log(Forward),Sigma*Sigma,greek,&RePsi1,&ImPsi1,&DRePsi1,&DImPsi1,&DDRePsi1,&DDImPsi1);
			err=PsiFunction (isVolInfFix,Maturity,0.,-x[i],a,b,Alpha,Rho,log(Forward),Sigma*Sigma,greek,&RePsi2,&ImPsi2,&DRePsi2,&DImPsi2,&DDRePsi2,&DDImPsi2);

			// runs over the strikes 

			for (j=iThresh+1;j<=nStrikes ;j++) {


				G12[j]+=w[i]/x[i]*exp(RePsi1)*sin(ImPsi1+x[i]*log(Strike[j]));
				G22[j]+=w[i]/x[i]*exp(RePsi2)*sin(ImPsi2+x[i]*log(Strike[j]));

			}
		}

			err = PsiFunction (isVolInfFix,Maturity,1.,0.,a,b,Alpha,Rho,log(Forward),Sigma*Sigma,greek,&RePsi1,&ImPsi1,&DRePsi1,&DImPsi1,&DDRePsi1,&DDImPsi1);
			G11=exp(RePsi1);

			err=PsiFunction (isVolInfFix,Maturity,0.,0.,a,b,Alpha,Rho,log(Forward),Sigma*Sigma,greek,&RePsi2,&ImPsi2,&DRePsi2,&DImPsi2,&DDRePsi2,&DDImPsi2);
			G21=exp(RePsi2);

			/* Now sums up all together */

			for (j=1;j<=nStrikes;j++) {
				
				n1=G11/2.-1./SRT_PI*G12[j];
				n2=G21/2.-1./SRT_PI*G22[j];

				switch (call_put){

					case SRT_CALL:
					result[j]=DMAX(1.e-13,(n1-Strike[j]*n2)*Disc)/dScaling;

					if (result[j]< (Forward-Strike[j])*Disc) result[j]=(Forward-Strike[j])*Disc/dScaling+1.e-8;

				break;

				default:
					result[j]= DMAX(1.e-13,(n1-Strike[j]*n2)*Disc/dScaling - (Forward-Strike[j])*Disc)/dScaling;

					if (result[j]< (-Forward+Strike[j])*Disc) result[j]=(-Forward+Strike[j])*Disc+1.e-8;	
				break;

				}
			}

			free_dvector(x,1,nSteps);
			free_dvector(w,1,nSteps);

			break;

		}

	break;

	case DELTA:

		x=dvector(1,nSteps);
		w=dvector(1,nSteps);
		GaussLeg(0., UpperBound, x, w, nSteps);

		for (j=1;j<=nStrikes;j++) {
			G12[j]=0.0;
			G22[j]=0.0;
		}

		for (i=1;i<=nSteps;i++) {

			err=PsiFunction (isVolInfFix,Maturity,1.,-x[i],a,b,Alpha,Rho,log(Forward),Sigma*Sigma,greek,&RePsi1,&ImPsi1,&DRePsi1,&DImPsi1,&DDRePsi1,&DDImPsi1);
			err=PsiFunction (isVolInfFix,Maturity,0.,-x[i],a,b,Alpha,Rho,log(Forward),Sigma*Sigma,greek,&RePsi2,&ImPsi2,&DRePsi2,&DImPsi2,&DDRePsi2,&DDImPsi2);

			// runs over the strikes 

			for (j=1;j<=nStrikes;j++) {

				SinT1=sin(ImPsi1+x[i]*log(Strike[j]));
				CosT1=cos(ImPsi1+x[i]*log(Strike[j]));

				SinT2=sin(ImPsi2+x[i]*log(Strike[j]));
				CosT2=cos(ImPsi2+x[i]*log(Strike[j]));

				G12[j]+=w[i]/x[i]*exp(RePsi1)*(SinT1*DRePsi1+CosT1*DImPsi1);
				G22[j]+=w[i]/x[i]*exp(RePsi2)*(SinT2*DRePsi2+CosT2*DImPsi2);

			}
		}

			err=PsiFunction (isVolInfFix,Maturity,1.,0.,a,b,Alpha,Rho,log(Forward),Sigma*Sigma,greek,&RePsi1,&ImPsi1,&DRePsi1,&DImPsi1,&DDRePsi1,&DDImPsi1);
			G11=exp(RePsi1)*DRePsi1;

			err=PsiFunction (isVolInfFix,Maturity,0.,0.,a,b,Alpha,Rho,log(Forward),Sigma*Sigma,greek,&RePsi2,&ImPsi2,&DRePsi2,&DImPsi2,&DDRePsi2,&DDImPsi2);
			G21=exp(RePsi2)*DRePsi2;

			for (j=1;j<=nStrikes;j++) {

				n1=G11/2.-1./SRT_PI*G12[j];
				n2=G21/2.-1./SRT_PI*G22[j];
		
				switch (call_put){

				case SRT_CALL:
					result[j]=DMAX(1.e-13,(n1-Strike[j]*n2)*Disc);

				break;

 				default:
					result[j]= DMAX(1.e-13,(n1-Strike[j]*n2)*Disc) - Disc;
	
				break;

				}
			}

		free_dvector(x,1,nSteps);
		free_dvector(w,1,nSteps);

		break;

	case GAMMA:

		x=dvector(1,nSteps);
		w=dvector(1,nSteps);
		GaussLeg(0., UpperBound, x, w, nSteps);

		for (j=1;j<=nStrikes;j++) {
			G12[j]=0.0;
			G22[j]=0.0;
		}

		for (i=1;i<=nSteps;i++) {

			err=PsiFunction (isVolInfFix,Maturity,1.,-x[i],a,b,Alpha,Rho,log(Forward),Sigma*Sigma,greek,&RePsi1,&ImPsi1,&DRePsi1,&DImPsi1,&DDRePsi1,&DDImPsi1);
			err=PsiFunction (isVolInfFix,Maturity,0.,-x[i],a,b,Alpha,Rho,log(Forward),Sigma*Sigma,greek,&RePsi2,&ImPsi2,&DRePsi2,&DImPsi2,&DDRePsi2,&DDImPsi2);

			// runs over the strikes 

			for (j=1;j<=nStrikes;j++) {

				SinT1=sin(ImPsi1+x[i]*log(Strike[j]));
				CosT1=cos(ImPsi1+x[i]*log(Strike[j]));
 
				SinT2=sin(ImPsi2+x[i]*log(Strike[j]));
				CosT2=cos(ImPsi2+x[i]*log(Strike[j]));

				G12[j]+=w[i]/x[i]*exp(RePsi1)*(SinT1*DRePsi1*DRePsi1+CosT1*DImPsi1*DRePsi1+DDRePsi1*SinT1+
					                           DRePsi1*DImPsi1*CosT1-DImPsi1*DImPsi1*SinT1+DDImPsi1*CosT1);
				G22[j]+=w[i]/x[i]*exp(RePsi2)*(SinT2*DRePsi2*DRePsi2+CosT2*DImPsi2*DRePsi2+DDRePsi2*SinT2+
					                           DRePsi2*DImPsi2*CosT2-DImPsi2*DImPsi2*SinT2+DDImPsi2*CosT2);

			}
		}

			err=PsiFunction (isVolInfFix,Maturity,1.,0.,a,b,Alpha,Rho,log(Forward),Sigma*Sigma,greek,&RePsi1,&ImPsi1,&DRePsi1,&DImPsi1,&DDRePsi1,&DDImPsi1);
			G11=exp(RePsi1)*(DRePsi1*DRePsi1+DDRePsi1);

			err=PsiFunction (isVolInfFix,Maturity,0.,0.,a,b,Alpha,Rho,log(Forward),Sigma*Sigma,greek,&RePsi2,&ImPsi2,&DRePsi2,&DImPsi2,&DDRePsi2,&DDImPsi2);
			G21=exp(RePsi2)*(DRePsi2*DRePsi2+DDRePsi2);


			for (j=1;j<=nStrikes;j++) {

				n1=G11/2.-1./SRT_PI*G12[j];
				n2=G21/2.-1./SRT_PI*G22[j];
		
				result[j]=DMAX(1.e-13,(n1-Strike[j]*n2)*Disc);	

			}

		free_dvector(x,1,nSteps);
		free_dvector(w,1,nSteps);
		break;

	case VEGA:

		x=dvector(1,nSteps);
		w=dvector(1,nSteps);
		GaussLeg(0., UpperBound, x, w, nSteps);

		for (j=1;j<=nStrikes;j++) {
			G12[j]=0.0;
			G22[j]=0.0;
		}

		for (i=1;i<=nSteps;i++) {

			err=PsiFunction (isVolInfFix,Maturity,1.,-x[i],a,b,Alpha,Rho,log(Forward),Sigma*Sigma,greek,&RePsi1,&ImPsi1,&DRePsi1,&DImPsi1,&DDRePsi1,&DDImPsi1);
			err=PsiFunction (isVolInfFix,Maturity,0.,-x[i],a,b,Alpha,Rho,log(Forward),Sigma*Sigma,greek,&RePsi2,&ImPsi2,&DRePsi2,&DImPsi2,&DDRePsi2,&DDImPsi2);

			// runs over the strikes 

			for (j=1;j<=nStrikes;j++) {

				SinT1=sin(ImPsi1+x[i]*log(Strike[j]));
				CosT1=cos(ImPsi1+x[i]*log(Strike[j]));

				SinT2=sin(ImPsi2+x[i]*log(Strike[j]));
				CosT2=cos(ImPsi2+x[i]*log(Strike[j]));

				G12[j]+=w[i]/x[i]*exp(RePsi1)*(SinT1*DRePsi1+CosT1*DImPsi1);
				G22[j]+=w[i]/x[i]*exp(RePsi2)*(SinT2*DRePsi2+CosT2*DImPsi2);

			}
		}


			err=PsiFunction (isVolInfFix,Maturity,1.,0.,a,b,Alpha,Rho,log(Forward),Sigma*Sigma,greek,&RePsi1,&ImPsi1,&DRePsi1,&DImPsi1,&DDRePsi1,&DDImPsi1);
			G11=exp(RePsi1)*DRePsi1;

			err=PsiFunction (isVolInfFix,Maturity,0.,0.,a,b,Alpha,Rho,log(Forward),Sigma*Sigma,greek,&RePsi2,&ImPsi2,&DRePsi2,&DImPsi2,&DDRePsi2,&DDImPsi2);
			G21=exp(RePsi2)*DRePsi2;


			for (j=1;j<=nStrikes;j++) {

				n1=G11/2.-1./SRT_PI*G12[j];
				n2=G21/2.-1./SRT_PI*G22[j];
		

					result[j]=(n1-Strike[j]*n2)*Disc;

			}

	free_dvector(x,1,nSteps);
	free_dvector(w,1,nSteps);
	break;

	case VOLGA:

		x=dvector(1,nSteps);
		w=dvector(1,nSteps);
		GaussLeg(0., UpperBound, x, w, nSteps);

		for (j=1;j<=nStrikes;j++) {
			G12[j]=0.0;
			G22[j]=0.0;
		}

		for (i=1;i<=nSteps;i++) {

			err=PsiFunction (isVolInfFix,Maturity,1.,-x[i],a,b,Alpha,Rho,log(Forward),Sigma*Sigma,greek,&RePsi1,&ImPsi1,&DRePsi1,&DImPsi1,&DDRePsi1,&DDImPsi1);
			err=PsiFunction (isVolInfFix,Maturity,0.,-x[i],a,b,Alpha,Rho,log(Forward),Sigma*Sigma,greek,&RePsi2,&ImPsi2,&DRePsi2,&DImPsi2,&DDRePsi2,&DDImPsi2);

			// runs over the strikes 

			for (j=1;j<=nStrikes;j++) {

				SinT1=sin(ImPsi1+x[i]*log(Strike[j]));
				CosT1=cos(ImPsi1+x[i]*log(Strike[j]));

				SinT2=sin(ImPsi2+x[i]*log(Strike[j]));
				CosT2=cos(ImPsi2+x[i]*log(Strike[j]));

				G12[j]+=w[i]/x[i]*exp(RePsi1)*(SinT1*DRePsi1*DRePsi1+CosT1*DImPsi1*DRePsi1+DDRePsi1*SinT1+
					                           DRePsi1*DImPsi1*CosT1-DImPsi1*DImPsi1*SinT1+DDImPsi1*CosT1);
				G22[j]+=w[i]/x[i]*exp(RePsi2)*(SinT2*DRePsi2*DRePsi2+CosT2*DImPsi2*DRePsi2+DDRePsi2*SinT2+
					                           DRePsi2*DImPsi2*CosT2-DImPsi2*DImPsi2*SinT2+DDImPsi2*CosT2);

			}
		}

			err=PsiFunction (isVolInfFix,Maturity,1.,0.,a,b,Alpha,Rho,log(Forward),Sigma*Sigma,greek,&RePsi1,&ImPsi1,&DRePsi1,&DImPsi1,&DDRePsi1,&DDImPsi1);
			G11=exp(RePsi1)*(DRePsi1*DRePsi1+DDRePsi1);

			err=PsiFunction (isVolInfFix,Maturity,0.,0.,a,b,Alpha,Rho,log(Forward),Sigma*Sigma,greek,&RePsi2,&ImPsi2,&DRePsi2,&DImPsi2,&DDRePsi2,&DDImPsi2);
			G21=exp(RePsi2)*(DRePsi2*DRePsi2+DDRePsi2);


			for (j=1;j<=nStrikes;j++) {

				n1=G11/2.-1./SRT_PI*G12[j];
				n2=G21/2.-1./SRT_PI*G22[j];
		
				result[j]=(n1-Strike[j]*n2)*Disc;	

			}

		free_dvector(x,1,nSteps);
		free_dvector(w,1,nSteps);
		break;

	case VANNA:

		x=dvector(1,nSteps);
		w=dvector(1,nSteps);
		GaussLeg(0., UpperBound, x, w, nSteps);

		for (j=1;j<=nStrikes;j++) {
			G12[j]=0.0;
			G22[j]=0.0;
		}

		for (i=1;i<=nSteps;i++) {

			err=PsiFunction (isVolInfFix,Maturity,1.,-x[i],a,b,Alpha,Rho,log(Forward),Sigma*Sigma,greek,&RePsi1,&ImPsi1,&DRePsi1,&DImPsi1,&DDRePsi1,&DDImPsi1);
			err=PsiFunction (isVolInfFix,Maturity,0.,-x[i],a,b,Alpha,Rho,log(Forward),Sigma*Sigma,greek,&RePsi2,&ImPsi2,&DRePsi2,&DImPsi2,&DDRePsi2,&DDImPsi2);

			// runs over the strikes 

			for (j=1;j<=nStrikes;j++) {

				SinT1=sin(ImPsi1+x[i]*log(Strike[j]));
				CosT1=cos(ImPsi1+x[i]*log(Strike[j]));

				SinT2=sin(ImPsi2+x[i]*log(Strike[j]));
				CosT2=cos(ImPsi2+x[i]*log(Strike[j]));

				G12[j]+=w[i]/x[i]*exp(RePsi1)*(SinT1*DDRePsi1*DRePsi1+CosT1*DRePsi1*DDImPsi1+
					                           DImPsi1*DDRePsi1*CosT1-DDImPsi1*DImPsi1*SinT1);
				G22[j]+=w[i]/x[i]*exp(RePsi2)*(SinT2*DDRePsi2*DRePsi2+CosT2*DRePsi2*DDImPsi2+
					                           DImPsi2*DDRePsi2*CosT2-DDImPsi2*DImPsi2*SinT2);

			}
		}

			err=PsiFunction (isVolInfFix,Maturity,1.,0.,a,b,Alpha,Rho,log(Forward),Sigma*Sigma,greek,&RePsi1,&ImPsi1,&DRePsi1,&DImPsi1,&DDRePsi1,&DDImPsi1);
			G11=exp(RePsi1)*(DRePsi1*DDRePsi1);

			err=PsiFunction (isVolInfFix,Maturity,0.,0.,a,b,Alpha,Rho,log(Forward),Sigma*Sigma,greek,&RePsi2,&ImPsi2,&DRePsi2,&DImPsi2,&DDRePsi2,&DDImPsi2);
			G21=exp(RePsi2)*(DRePsi2*DDRePsi2);


			for (j=1;j<=nStrikes;j++) {

				n1=G11/2.-1./SRT_PI*G12[j];
				n2=G21/2.-1./SRT_PI*G22[j];
		
				result[j]=(n1-Strike[j]*n2)*Disc;	

			}

		free_dvector(x,1,nSteps);
		free_dvector(w,1,nSteps);
		break;

	case THETA:

		x=dvector(1,nSteps);
		w=dvector(1,nSteps);
		GaussLeg(0., UpperBound, x, w, nSteps);

		/* I first compute the derivative of the term containing the DF */

		greek = PREMIUM;
		for (j=1;j<=nStrikes;j++) {
			G12[j]=0.0;
			G22[j]=0.0;
		}

		for (i=1;i<=nSteps;i++) {

			err=PsiFunction (isVolInfFix,Maturity,1.,-x[i],a,b,Alpha,Rho,log(Forward),Sigma*Sigma,greek,&RePsi1,&ImPsi1,&DRePsi1,&DImPsi1,&DDRePsi1,&DDImPsi1);
			err=PsiFunction (isVolInfFix,Maturity,0.,-x[i],a,b,Alpha,Rho,log(Forward),Sigma*Sigma,greek,&RePsi2,&ImPsi2,&DRePsi2,&DImPsi2,&DDRePsi2,&DDImPsi2);

			// runs over the strikes 

			for (j=1;j<=nStrikes;j++) {

				G12[j]+=w[i]/x[i]*exp(RePsi1)*sin(ImPsi1+x[i]*log(Strike[j]));
				G22[j]+=w[i]/x[i]*exp(RePsi2)*sin(ImPsi2+x[i]*log(Strike[j]));

			}
		}

			err=PsiFunction (isVolInfFix,Maturity,1.,0.,a,b,Alpha,Rho,log(Forward),Sigma*Sigma,greek,&RePsi1,&ImPsi1,&DRePsi1,&DImPsi1,&DDRePsi1,&DDImPsi1);
			G11=exp(RePsi1);

			err=PsiFunction (isVolInfFix,Maturity,0.,0.,a,b,Alpha,Rho,log(Forward),Sigma*Sigma,greek,&RePsi2,&ImPsi2,&DRePsi2,&DImPsi2,&DDRePsi2,&DDImPsi2);
			G21=exp(RePsi2);


			for (j=1;j<=nStrikes;j++) {

				n1=G11/2.-1./SRT_PI*G12[j];
				n2=G21/2.-1./SRT_PI*G22[j];
		
				switch (call_put){

					case SRT_CALL:
					result[j]=-(n1-Strike[j]*n2)*Disc*log(Disc)/Maturity;

				break;

				default:
					result[j]= -(n1-Strike[j]*n2)*Disc*log(Disc)/Maturity
						+ (Forward-Strike[j])*Disc*log(Disc)/Maturity;

				break;

				}
			}

		/* I now compute the derivative of the term equivalent to [Fwd*N(d1)-KN(d2)] */

		greek = THETA;
		for (j=1;j<=nStrikes;j++) {
			G12[j]=0.0;
			G22[j]=0.0;
		}

		for (i=1;i<=nSteps;i++) {

			err=PsiFunction (isVolInfFix,Maturity,1.,-x[i],a,b,Alpha,Rho,log(Forward),Sigma*Sigma,greek,&RePsi1,&ImPsi1,&DRePsi1,&DImPsi1,&DDRePsi1,&DDImPsi1);
			err=PsiFunction (isVolInfFix,Maturity,0.,-x[i],a,b,Alpha,Rho,log(Forward),Sigma*Sigma,greek,&RePsi2,&ImPsi2,&DRePsi2,&DImPsi2,&DDRePsi2,&DDImPsi2);

			// runs over the strikes 

			for (j=1;j<=nStrikes;j++) {

				SinT1=sin(ImPsi1+x[i]*log(Strike[j]));
				CosT1=cos(ImPsi1+x[i]*log(Strike[j]));

				SinT2=sin(ImPsi2+x[i]*log(Strike[j]));
				CosT2=cos(ImPsi2+x[i]*log(Strike[j]));

				G12[j]+=w[i]/x[i]*exp(RePsi1)*(SinT1*DRePsi1+CosT1*DImPsi1);
				G22[j]+=w[i]/x[i]*exp(RePsi2)*(SinT2*DRePsi2+CosT2*DImPsi2);

			}
		}

			err=PsiFunction (isVolInfFix,Maturity,1.,0.,a,b,Alpha,Rho,log(Forward),Sigma*Sigma,greek,&RePsi1,&ImPsi1,&DRePsi1,&DImPsi1,&DDRePsi1,&DDImPsi1);
			G11=exp(RePsi1)*DRePsi1;

			err=PsiFunction (isVolInfFix,Maturity,0.,0.,a,b,Alpha,Rho,log(Forward),Sigma*Sigma,greek,&RePsi2,&ImPsi2,&DRePsi2,&DImPsi2,&DDRePsi2,&DDImPsi2);
			G21=exp(RePsi2)*DRePsi2;


			for (j=1;j<=nStrikes;j++) {

				n1=G11/2.-1./SRT_PI*G12[j];
				n2=G21/2.-1./SRT_PI*G22[j];
	
				result[j]+=-(n1-Strike[j]*n2)*Disc;

			}

		free_dvector(x,1,nSteps);
		free_dvector(w,1,nSteps);
		break;

	case DENSITY: 

		x=dvector(1,nSteps);
		w=dvector(1,nSteps);
		GaussLeg(0., UpperBound, x, w, nSteps);

		/* now I split the integration in two parts */

		for (i=1;i<=nSteps;i++) {

			err=PsiFunction (isVolInfFix,Maturity,1.,-x[i],a,b,Alpha,Rho,log(Forward),Sigma*Sigma,greek,&RePsi1,&ImPsi1,&DRePsi1,&DImPsi1,&DDRePsi1,&DDImPsi1);
			err=PsiFunction (isVolInfFix,Maturity,0.,-x[i],a,b,Alpha,Rho,log(Forward),Sigma*Sigma,greek,&RePsi2,&ImPsi2,&DRePsi2,&DImPsi2,&DDRePsi2,&DDImPsi2);

			// runs over the strikes 

			for (j=1;j<=nStrikes;j++) {

				G12[j]+=w[i]*exp(RePsi1)*(-cos(ImPsi1+x[i]*log(Strike[j]))-x[i]*sin(ImPsi1+x[i]*log(Strike[j])))/Strike[j]/Strike[j];
				G22[j]+=w[i]*exp(RePsi2)*(-cos(ImPsi2+x[i]*log(Strike[j]))-x[i]*sin(ImPsi2+x[i]*log(Strike[j])))/Strike[j]/Strike[j];
				GKK[j]+=w[i]*exp(RePsi2)*cos(ImPsi2+x[i]*log(Strike[j]))/Strike[j];
				
			}
		}

			for (j=nStrikes;j>=1;j--) {

				n1=-1./SRT_PI*G12[j];
				n2=-1./SRT_PI*G22[j];
				n3=-1./SRT_PI*GKK[j];

					result[j]=DMAX(0.0,n1-2.*n3-Strike[j]*n2);

			}


		free_dvector(x,1,nSteps);
		free_dvector(w,1,nSteps);
		break;

	case CUMDENSITY: 

		switch (IntegrType){

		case 1:

		x=dvector(1,nStp);
		w=dvector(1,nStp);
		x1=dvector(1,nStp);
		w1=dvector(1,nStp);

		for (j=1;j<=nStrikes;j++) {

			G12[j]=0.0;
			G22[j]=0.0;
			GKK[j]=0.0;
		}

			ImAppr2=-(log(Forward)+(esp2-1)*Sigma*Sigma/2./(b-Alpha*Rho)+Maturity*a/2./(b-Alpha*Rho)
				                +a*(esp2-1.)/2/(b-Alpha*Rho)/(b-Alpha*Rho));
			ImAppr1=-(log(Forward)+Sigma*Sigma/b*(1.-esp)-Maturity*a/2./b+a/2./b/b*(1.-esp));


		for (j=1;j<=nStrikes;j++) {

			iter=0;
			NodeLeft=0.0;
			NodeLeft1=0.0;
			do  {

			/* find the next node on the axis */

			NodeRight = HestonFindRightNode(ImAppr1,Strike[j],NodeLeft,isVolInfFix,Maturity,
											a,b,Alpha,Rho,Forward,Sigma,1.0,greek);
			GaussLeg(NodeLeft, NodeRight, x, w, nStp);

			NodeRight1 = HestonFindRightNode(ImAppr2,Strike[j],NodeLeft,isVolInfFix,Maturity,
											a,b,Alpha,Rho,Forward,Sigma,0.0,greek);
			GaussLeg(NodeLeft1, NodeRight1, x1, w1, nStp);

				Int=0.0;
				Int1=0.0;
				for (k=1;k<=nStp;k++) {


					err=PsiFunction (isVolInfFix,Maturity,1.,-x[k],a,b,Alpha,Rho,log(Forward),
						             Sigma*Sigma,greek,&RePsi1,&ImPsi1,
									 &DRePsi1,&DImPsi1,&DDRePsi1,&DDImPsi1);

					G12[j]+=w[k]*exp(RePsi1)*cos(ImPsi1+x[k]*log(Strike[j]))/Strike[j];
 

					err=PsiFunction (isVolInfFix,Maturity,0.,-x1[k],a,b,Alpha,Rho,log(Forward),
									 Sigma*Sigma,greek,&RePsi2,&ImPsi2,
									 &DRePsi2,&DImPsi2,&DDRePsi2,&DDImPsi2);
	
				    G22[j]+=w1[k]*exp(RePsi2)*cos(ImPsi2+x1[k]*log(Strike[j]))/Strike[j];
					GKK[j]+=w1[k]/x1[k]*exp(RePsi2)*sin(ImPsi2+x1[k]*log(Strike[j]));
			
				}

			NodeLeft=NodeRight;
			NodeLeft1=NodeRight1;

			iter++;

			} while ((NodeLeft <= UpperBound) && (iter <= 1000));
		}

			/* computes the remaining two terms that need no integration */

			err = PsiFunction (isVolInfFix,Maturity,1.,0.,a,b,Alpha,Rho,log(Forward),Sigma*Sigma,greek,&RePsi1,&ImPsi1,&DRePsi1,&DImPsi1,&DDRePsi1,&DDImPsi1);
			G11=exp(RePsi1);

			err=PsiFunction (isVolInfFix,Maturity,0.,0.,a,b,Alpha,Rho,log(Forward),Sigma*Sigma,greek,&RePsi2,&ImPsi2,&DRePsi2,&DImPsi2,&DDRePsi2,&DDImPsi2);
			G21=exp(RePsi2);

			/* finally puts all together and evaluates the option */

			for (j=1;j<=nStrikes;j++) {
				
				n1=-1./SRT_PI*G12[j];
				n2=-1./SRT_PI*G22[j];
				n3=-1./SRT_PI*GKK[j];

				result[j]=DMAX(DMIN(1.0,0.5+(n1-n3-Strike[j]*n2)),0.0);

				
			}

			free_dvector(x,1,nStp);
			free_dvector(w,1,nStp);
			free_dvector(x1,1,nStp);
			free_dvector(w1,1,nStp);

		break;

		default:

		x=dvector(1,nSteps);
		w=dvector(1,nSteps);
		GaussLeg(0., UpperBound, x, w, nSteps);

		for (j=1;j<=nStrikes;j++) {

			G12[j]=0.0;
			G22[j]=0.0;
			GKK[j]=0.0;
		}

		for (i=1;i<=nSteps;i++) {

			err=PsiFunction (isVolInfFix,Maturity,1.,-x[i],a,b,Alpha,Rho,log(Forward),Sigma*Sigma,greek,&RePsi1,&ImPsi1,&DRePsi1,&DImPsi1,&DDRePsi1,&DDImPsi1);
			err=PsiFunction (isVolInfFix,Maturity,0.,-x[i],a,b,Alpha,Rho,log(Forward),Sigma*Sigma,greek,&RePsi2,&ImPsi2,&DRePsi2,&DImPsi2,&DDRePsi2,&DDImPsi2);

			// runs over the strikes 

			for (j=1;j<=nStrikes;j++) {


				G12[j]+=w[i]*exp(RePsi1)*cos(ImPsi1+x[i]*log(Strike[j]))/Strike[j];
				G22[j]+=w[i]*exp(RePsi2)*cos(ImPsi2+x[i]*log(Strike[j]))/Strike[j];
				GKK[j]+=w[i]/x[i]*exp(RePsi2)*sin(ImPsi2+x[i]*log(Strike[j]));
				
			}
		}

			for (j=1;j<=nStrikes;j++) {

				n1=-1./SRT_PI*G12[j];
				n2=-1./SRT_PI*G22[j];
				n3=-1./SRT_PI*GKK[j];

				result[j]=DMAX(DMIN(1.0,0.5+(n1-n3-Strike[j]*n2)),0.0);

			}

		free_dvector(x,1,nSteps);
		free_dvector(w,1,nSteps);

		}

		break;
	}

	Forward-=Gamma;

	for (j=1; j<=nStrikes;j++)
	{
		Strike[j] -= Gamma;
	}

	free_dvector(G12,1,nStrikes);
	free_dvector(G22,1,nStrikes);
	free_dvector(GKK,1,nStrikes);
	free_dvector(premium,1,nStrikes);

	return NULL;
}

/*---------------------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------------------*/

#define NRANSI
#define MAXIT 60
#define UNUSED (-1.11e30)
											
double HestonFindRightNode(
						   double ImAppr,
						   double Strike,
						   double NodeLeft,
	   					   SRT_Boolean isVolInfFix,
						   double Maturity,
						   double a,
						   double b,
						   double Alpha,
						   double Rho,
						   double Forward,
						   double Sigma,
						   double U,
						   SrtGreekType greek
						   )
{	
	double Node=0.0,x1,x2,RePsi,ImPsi,DRePsi,DImPsi,DDRePsi,DDImPsi;
	int j;
	double ans,fh,fl,fm,fnew,s,xh,xl,xm,xnew,xacc=1.e-7;
	Err err=NULL;

	/* finds x1, the left boundary */

	x1 = NodeLeft+fabs(SRT_PI/2./(ImAppr+log(Strike)));
	x2 = NodeLeft+fabs(SRT_PI*3./2./(ImAppr+log(Strike)));

	err=PsiFunction (isVolInfFix,Maturity,U,-x1,a,b,Alpha,Rho,log(Forward),
		             Sigma*Sigma,greek,&RePsi,&ImPsi,
					 &DRePsi,&DImPsi,&DDRePsi,&DDImPsi);
	fl = sin(fabs(ImPsi+x1*log(Strike)));

	err=PsiFunction (isVolInfFix,Maturity,U,-x2,a,b,Alpha,Rho,log(Forward),
		             Sigma*Sigma,greek,&RePsi,&ImPsi,
					 &DRePsi,&DImPsi,&DDRePsi,&DDImPsi);
	fh = sin(fabs(ImPsi+x2*log(Strike)));

	if ((fl > 0.0 && fh < 0.0) || (fl < 0.0 && fh > 0.0)) {
		xl=x1;
		xh=x2;
		ans=UNUSED;
		for (j=1;j<=MAXIT;j++) {
			xm=0.5*(xl+xh);

			err=PsiFunction (isVolInfFix,Maturity,U,-xm,a,b,Alpha,Rho,log(Forward),
							 Sigma*Sigma,greek,&RePsi,&ImPsi,
							 &DRePsi,&DImPsi,&DDRePsi,&DDImPsi);
			fm = sin(fabs(ImPsi+xm*log(Strike)));

			s=sqrt(fm*fm-fl*fh);
			if (s == 0.0) return ans;
			xnew=xm+(xm-xl)*((fl >= fh ? 1.0 : -1.0)*fm/s);
			if (fabs(xnew-ans) <= xacc) return ans;
			ans=xnew;


			err=PsiFunction (isVolInfFix,Maturity,U,-ans,a,b,Alpha,Rho,log(Forward),
							 Sigma*Sigma,greek,&RePsi,&ImPsi,
							 &DRePsi,&DImPsi,&DDRePsi,&DDImPsi);
			fnew = sin(fabs(ImPsi+ans*log(Strike)));

			if (fnew == 0.0) return ans;
			if (SIGN(fm,fnew) != fm) {
				xl=xm;
				fl=fm;
				xh=ans;
				fh=fnew;
			} else if (SIGN(fl,fnew) != fl) {
				xh=ans;
				fh=fnew;
			} else if (SIGN(fh,fnew) != fh) {
				xl=ans;
				fl=fnew;
			} else {
				return 0.0;
			}
			if (fabs(xh-xl) <= xacc) return ans;
		}
		return 0.0;
		 // nrerror("zriddr exceed maximum iterations");
	}
	else {
		if (fl == 0.0) return x1;
		if (fh == 0.0) return x2;
		// nrerror("root must be bracketed in zriddr.");
	}
	
	return 0.0;
}	
	
#undef MAXIT
#undef UNUSED
#undef NRANSI

/*---------------------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------------------*/

Err PsiFunction (
					SRT_Boolean isVolInfFix,
					double T,
					double U,
					double V,
					double a,
					double b,
					double alpha,
					double rho,
					double x1,
					double x2,
					SrtGreekType greek,
					double *RePsi,
					double *ImPsi,
					double *DRePsi,
					double *DImPsi,
					double *DDRePsi,
					double *DDImPsi
					)

{

	double	Phi1,Phi2,Phi3,Phi4,Psi1,Psi2,Psi3,Psi4,NumRe,NumIm,DenRe,DenIm, Chi, Den, DDen,
			PhiG, PsiG, PhiG2, PsiG2, Phi3G, Psi3G, SinT, Sin2T, CosT, Cos2T, U2=U*U, U3=U*U*U, 
			V2=V*V, V3=V*V*V, r=0.0, NumRe2,NumIm2, Phi2G,Psi2G,Esp,Esp2,DNumRe,DNumIm,DLogArg;
	double	alpha2=alpha*alpha, alpharho=alpha*rho,alpharho2=alpharho*alpharho,
			Talpha2=T/alpha2,ralpha2=r*alpha2,aalpha2=a/alpha2,
			Ualpharho=U*alpharho,Valpharho=V*alpharho,
			PowG,ATanG, VPhiG,VPsiG, Den2ReDen2Im,
			SinTEsp,CosTEsp,Sin2TEsp2,Cos2TEsp2;

	Err err = NULL;

	//  Defines the Re and Im part of the square of the function gamma 

	PhiG2=		b*b + alpha2*(U - U2 + V2) - 
				2*b*Ualpharho + alpharho2*(U2 - V2);

	PsiG2=	V*(alpha2*(1 - 2*U) - 
			2*(b*alpharho - U*alpharho2));

	//  Defines the Re and Im part of the function gamma 

	PowG=pow(PhiG2*PhiG2+PsiG2*PsiG2,0.25);
	ATanG=atan(PsiG2/PhiG2)/2.0;
	PhiG=PowG*cos(ATanG);
	PsiG=PowG*sin(ATanG);
	VPhiG=V*PhiG;
	VPsiG=V*PsiG;
	
	Phi2G=PhiG*PhiG;
	Psi2G=PsiG*PsiG;
	Phi3G=Phi2G*PhiG;
	Psi3G=Psi2G*PsiG;

	SinT=sin(T*PsiG);
	Sin2T=SinT*SinT;
	CosT=cos(T*PsiG);
	Cos2T=CosT*CosT;
	Esp=exp(-T*PhiG);
	Esp2=Esp*Esp;
	SinTEsp=SinT*Esp;
	CosTEsp=CosT*Esp;
	Sin2TEsp2=Sin2T*Esp2;
	Cos2TEsp2=Cos2T*Esp2;

	//  Defines the L1 function 
	Phi1=x1*U;
	Psi1=x1*V;

	//  Defines the L2 function 
	NumRe=
		-x2*((1 - 2*CosTEsp + Cos2TEsp2 + Sin2TEsp2)*(b*(U - U2 + V2) - alpharho*(U2 - U3 + V2 - U*V2)) + 
		(1 - Cos2TEsp2 - Sin2TEsp2)*(PhiG*(U - U2 + V2) + VPsiG*(1 - 2*U)) - 
		2*SinTEsp*(VPhiG*(1 - 2*U) - PsiG*(U - U2 + V2)));


	NumIm=
		-x2*((1 - 2*CosTEsp + Cos2TEsp2 + Sin2TEsp2)*(b*V*(1 - 2*U) + alpharho*(U2*V + V3)) + 
		(1 - Cos2TEsp2 - Sin2TEsp2)*(VPhiG*(1 - 2*U) - PsiG*(U - U2 + V2)) + 
		2*SinTEsp*(PhiG*(U - U2 + V2) + VPsiG*(1 - 2*U)));


	DenRe=
		b - Ualpharho + PhiG - 
		Esp*(CosT*(b - Ualpharho - PhiG) - 
		SinT*(Valpharho + PsiG));

	DenIm=
		-Valpharho + PsiG + 
		Esp*(CosT*(Valpharho + PsiG) + 
		SinT*(b - Ualpharho - PhiG));

	Den2ReDen2Im=1.0/(DenRe*DenRe+DenIm*DenIm);
	Phi2=NumRe*Den2ReDen2Im;
	Psi2=NumIm*Den2ReDen2Im;

	//  Defines the L3 function 

	Phi3=Talpha2*(ralpha2*(U-1)+a*(b-Ualpharho-PhiG));
	Psi3=Talpha2*(ralpha2*V-a*(Valpharho+PsiG));


	//  Defines the L4 function 

	Chi=Phi2G+Psi2G;
	Den=2*Chi;


	NumRe=
		PhiG*(b - Ualpharho) - Phi2G + 2*Chi - Valpharho*PsiG - Psi2G - 
		Esp*(CosT*(PhiG*(b - Ualpharho) - Phi2G - Valpharho*PsiG - Psi2G) - 
		SinT*(Valpharho*PhiG + PsiG*(b - Ualpharho)));

	NumIm=
		-Valpharho*PhiG - PsiG*(b - Ualpharho) + 
		Esp*(CosT*(Valpharho*PhiG + PsiG*(b - Ualpharho)) + 
		SinT*(PhiG*(b - Ualpharho) - Phi2G - Valpharho*PsiG -  Psi2G));

	Phi4=-aalpha2*log((NumRe*NumRe+NumIm*NumIm)/(Den*Den));
	Psi4=-2*aalpha2*atan(NumIm/NumRe);

	// computes Re and Im part of the Psi function 

	*RePsi=Phi1+Phi2+Phi3+Phi4;
	*ImPsi=Psi1+Psi2+Psi3+Psi4;	

	// Now I compute the greeks 

	switch (greek) {


	case PREMIUM:

		*DRePsi=0.0;
		*DImPsi=0.0;
		*DDRePsi=0.0;
		*DDImPsi=0.0;
		break;


	case DELTA:

		x1=exp(x1); // exp(x1) is the fwd. 
		*DRePsi=1./x1*U;
		*DImPsi=1./x1*V;
		*DDRePsi=0.0;
		*DDImPsi=0.0;

		break;

	case GAMMA:

		x1=exp(x1); // exp(x1) is the fwd. 
		*DRePsi=1./x1*U;
		*DImPsi=1./x1*V;
		*DDRePsi=-1./x1/x1*U;
		*DDImPsi=-1./x1/x1*V;

		break;

	case VEGA:

//		a=2.*b*sqrt(x2);
		x2=2.*sqrt(x2);

		NumRe=	-(b*U*x2) + b*U2*x2 - b*V2*x2 + U2*x2*alpha*rho - 
					U3*x2*alpha*rho + V2*x2*alpha*rho - 
					U*V2*x2*alpha*rho - U*x2*PhiG + U2*x2*PhiG - 
					V2*x2*PhiG - V*x2*PsiG + 2*U*V*x2*PsiG + 
					(2*b*U*x2*CosT)*Esp - 
					(2*b*U2*x2*CosT)*Esp + 
					(2*b*V2*x2*CosT)*Esp - 
					(2*U2*x2*alpha*rho*CosT)*Esp + 
					(2*U3*x2*alpha*rho*CosT)*Esp - 
					(2*V2*x2*alpha*rho*CosT)*Esp + 
					(2*U*V2*x2*alpha*rho*CosT)*Esp - 
					(b*U*x2*Cos2T)*Esp2 + 
					(b*U2*x2*Cos2T)*Esp2 - 
					(b*V2*x2*Cos2T)*Esp2 + 
					(U2*x2*alpha*rho*Cos2T)*Esp2 - 
					(U3*x2*alpha*rho*Cos2T)*Esp2 + 
					(V2*x2*alpha*rho*Cos2T)*Esp2 - 
					(U*V2*x2*alpha*rho*Cos2T)*Esp2 + 
					(U*x2*PhiG*Cos2T)*Esp2 - 
					(U2*x2*PhiG*Cos2T)*Esp2 + 
					(V2*x2*PhiG*Cos2T)*Esp2 + 
					(V*x2*PsiG*Cos2T)*Esp2 - 
					(2*U*V*x2*PsiG*Cos2T)*Esp2 + 
					(2*V*x2*PhiG*SinT)*Esp - 
					(4*U*V*x2*PhiG*SinT)*Esp - 
					(2*U*x2*PsiG*SinT)*Esp + 
					(2*U2*x2*PsiG*SinT)*Esp - 
					(2*V2*x2*PsiG*SinT)*Esp - 
					(b*U*x2*Sin2T)*Esp2 + 
					(b*U2*x2*Sin2T)*Esp2 - 
					(b*V2*x2*Sin2T)*Esp2 + 
					(U2*x2*alpha*rho*Sin2T)*Esp2 - 
					(U3*x2*alpha*rho*Sin2T)*Esp2 + 
					(V2*x2*alpha*rho*Sin2T)*Esp2 - 
					(U*V2*x2*alpha*rho*Sin2T)*Esp2 + 
					(U*x2*PhiG*Sin2T)*Esp2 - 
					(U2*x2*PhiG*Sin2T)*Esp2 + 
					(V2*x2*PhiG*Sin2T)*Esp2 + 
					(V*x2*PsiG*Sin2T)*Esp2 - 
					(2*U*V*x2*PsiG*Sin2T)*Esp2;

		NumIm=   -(b*V*x2) + 2*b*U*V*x2 - U2*V*x2*alpha*rho -
					V3*x2*alpha*rho - V*x2*PhiG + 2*U*V*x2*PhiG + 
					U*x2*PsiG - U2*x2*PsiG + V2*x2*PsiG + (2*b*V*x2*CosT)*Esp - 
					(4*b*U*V*x2*CosT)*Esp + (2*U2*V*x2*alpha*rho*CosT)*Esp +
					(2*V3*x2*alpha*rho*CosT)*Esp - 
					(b*V*x2*Cos2T)*Esp2 +  
					(2*b*U*V*x2*Cos2T)*Esp2 - 
					(U2*V*x2*alpha*rho*Cos2T)*Esp2 - 
					(V3*x2*alpha*rho*Cos2T)*Esp2 + 
					(V*x2*PhiG*Cos2T)*Esp2 - 
					(2*U*V*x2*PhiG*Cos2T)*Esp2 - 
					(U*x2*PsiG*Cos2T)*Esp2 + 
					(U2*x2*PsiG*Cos2T)*Esp2 - 
					(V2*x2*PsiG*Cos2T)*Esp2 - 
					(2*U*x2*PhiG*SinT)*Esp + 
					(2*U2*x2*PhiG*SinT)*Esp - 
					(2*V2*x2*PhiG*SinT)*Esp - 
					(2*V*x2*PsiG*SinT)*Esp + 
					(4*U*V*x2*PsiG*SinT)*Esp - 
					(b*V*x2*Sin2T)*Esp2 + 
					(2*b*U*V*x2*Sin2T)*Esp2 - 
					(U2*V*x2*alpha*rho*Sin2T)*Esp2 - 
					(V3*x2*alpha*rho*Sin2T)*Esp2 + 
					(V*x2*PhiG*Sin2T)*Esp2 - 
					(2*U*V*x2*PhiG*Sin2T)*Esp2 - 
					(U*x2*PsiG*Sin2T)*Esp2 + 
					(U2*x2*PsiG*Sin2T)*Esp2 - 
					(V2*x2*PsiG*Sin2T)*Esp2;

	DenRe=
		b - U*alpha*rho + PhiG - (b*CosT)*Esp + 
		(U*alpha*rho*CosT)*Esp + (PhiG*CosT)*Esp + 
		(V*alpha*rho*SinT)*Esp + (PsiG*SinT)*Esp;

	DenIm=
		-(V*alpha*rho) + PsiG + (V*alpha*rho*CosT)*Esp + 
		(PsiG*CosT)*Esp + (b*SinT)*Esp - (U*alpha*rho*SinT)*Esp -
		(PhiG*SinT)*Esp;

		*DRePsi=NumRe/(DenRe*DenRe+DenIm*DenIm);
		*DImPsi=NumIm/(DenRe*DenRe+DenIm*DenIm);
		*DDRePsi=0.0;
		*DDImPsi=0.0;

	/*
	This part is needed if one is interested in the case when Sigma(0) = Sigma(infinity)
    */

		if (isVolInfFix) {

			a=b*x2;
			Phi3=T/alpha/alpha*(a*b-a*U*alpha*rho-a*PhiG);
			Psi3=T/alpha/alpha*(-a*V*alpha*rho-a*PsiG);
			Chi=PhiG*PhiG+PsiG*PsiG;
			Den=2*Chi;

			NumRe2=
				b*PhiG - U*alpha*rho*PhiG - Phi2G + 2*Chi - V*alpha*
				rho*PsiG - Psi2G -  (b*PhiG*CosT)*Esp + 
				(U*alpha*rho*PhiG*CosT)*Esp + (Phi2G*CosT)*Esp + (V*alpha*rho*
				PsiG*CosT)*Esp +  (Psi2G*CosT)*Esp + (V*alpha*rho*
				PhiG*SinT)*Esp +  (b*PsiG*SinT)*Esp - (U*alpha*rho*PsiG*SinT)*Esp;

			NumIm2=
				-(V*alpha*rho*PhiG) - b*PsiG + U*alpha*rho*PsiG + (V*
				alpha*rho*PhiG*CosT)*Esp + (b*PsiG*CosT)*Esp - 
				(U*alpha*rho*PsiG*CosT)*Esp + (b*PhiG*SinT)*Esp - 
				(U*alpha*rho*PhiG*SinT)*Esp - (Phi2G*SinT)*Esp - (V*alpha*rho*
				PsiG*SinT)*Esp -  (Psi2G*SinT)*Esp;

			Phi4=-a/alpha/alpha*log(NumRe2*NumRe2/Den/Den+NumIm2*NumIm2/Den/Den);
			Psi4=-2.*a/alpha/alpha*atan(NumIm2/NumRe2);


				*DRePsi+=Phi3+Phi4;
				*DImPsi+=Psi3+Psi4;
				*DDRePsi=0.0;
				*DDImPsi=0.0;
		}

		break;

	case VANNA:

		x1=exp(x1); // exp(x1) is the fwd. 
		*DRePsi=1./x1*U;
		*DImPsi=1./x1*V;				
		

		x2=2.*sqrt(x2);

		NumRe=	-(b*U*x2) + b*U2*x2 - b*V2*x2 + U2*x2*alpha*rho - 
					U3*x2*alpha*rho + V2*x2*alpha*rho - 
					U*V2*x2*alpha*rho - U*x2*PhiG + U2*x2*PhiG - 
					V2*x2*PhiG - V*x2*PsiG + 2*U*V*x2*PsiG + 
					(2*b*U*x2*CosT)*Esp - 
					(2*b*U2*x2*CosT)*Esp + 
					(2*b*V2*x2*CosT)*Esp - 
					(2*U2*x2*alpha*rho*CosT)*Esp + 
					(2*U3*x2*alpha*rho*CosT)*Esp - 
					(2*V2*x2*alpha*rho*CosT)*Esp + 
					(2*U*V2*x2*alpha*rho*CosT)*Esp - 
					(b*U*x2*Cos2T)*Esp2 + 
					(b*U2*x2*Cos2T)*Esp2 - 
					(b*V2*x2*Cos2T)*Esp2 + 
					(U2*x2*alpha*rho*Cos2T)*Esp2 - 
					(U3*x2*alpha*rho*Cos2T)*Esp2 + 
					(V2*x2*alpha*rho*Cos2T)*Esp2 - 
					(U*V2*x2*alpha*rho*Cos2T)*Esp2 + 
					(U*x2*PhiG*Cos2T)*Esp2 - 
					(U2*x2*PhiG*Cos2T)*Esp2 + 
					(V2*x2*PhiG*Cos2T)*Esp2 + 
					(V*x2*PsiG*Cos2T)*Esp2 - 
					(2*U*V*x2*PsiG*Cos2T)*Esp2 + 
					(2*V*x2*PhiG*SinT)*Esp - 
					(4*U*V*x2*PhiG*SinT)*Esp - 
					(2*U*x2*PsiG*SinT)*Esp + 
					(2*U2*x2*PsiG*SinT)*Esp - 
					(2*V2*x2*PsiG*SinT)*Esp - 
					(b*U*x2*Sin2T)*Esp2 + 
					(b*U2*x2*Sin2T)*Esp2 - 
					(b*V2*x2*Sin2T)*Esp2 + 
					(U2*x2*alpha*rho*Sin2T)*Esp2 - 
					(U3*x2*alpha*rho*Sin2T)*Esp2 + 
					(V2*x2*alpha*rho*Sin2T)*Esp2 - 
					(U*V2*x2*alpha*rho*Sin2T)*Esp2 + 
					(U*x2*PhiG*Sin2T)*Esp2 - 
					(U2*x2*PhiG*Sin2T)*Esp2 + 
					(V2*x2*PhiG*Sin2T)*Esp2 + 
					(V*x2*PsiG*Sin2T)*Esp2 - 
					(2*U*V*x2*PsiG*Sin2T)*Esp2;

		NumIm=   -(b*V*x2) + 2*b*U*V*x2 - U2*V*x2*alpha*rho -
					V3*x2*alpha*rho - V*x2*PhiG + 2*U*V*x2*PhiG + 
					U*x2*PsiG - U2*x2*PsiG + V2*x2*PsiG + (2*b*V*x2*CosT)*Esp - 
					(4*b*U*V*x2*CosT)*Esp + (2*U2*V*x2*alpha*rho*CosT)*Esp +
					(2*V3*x2*alpha*rho*CosT)*Esp - 
					(b*V*x2*Cos2T)*Esp2 +  
					(2*b*U*V*x2*Cos2T)*Esp2 - 
					(U2*V*x2*alpha*rho*Cos2T)*Esp2 - 
					(V3*x2*alpha*rho*Cos2T)*Esp2 + 
					(V*x2*PhiG*Cos2T)*Esp2 - 
					(2*U*V*x2*PhiG*Cos2T)*Esp2 - 
					(U*x2*PsiG*Cos2T)*Esp2 + 
					(U2*x2*PsiG*Cos2T)*Esp2 - 
					(V2*x2*PsiG*Cos2T)*Esp2 - 
					(2*U*x2*PhiG*SinT)*Esp + 
					(2*U2*x2*PhiG*SinT)*Esp - 
					(2*V2*x2*PhiG*SinT)*Esp - 
					(2*V*x2*PsiG*SinT)*Esp + 
					(4*U*V*x2*PsiG*SinT)*Esp - 
					(b*V*x2*Sin2T)*Esp2 + 
					(2*b*U*V*x2*Sin2T)*Esp2 - 
					(U2*V*x2*alpha*rho*Sin2T)*Esp2 - 
					(V3*x2*alpha*rho*Sin2T)*Esp2 + 
					(V*x2*PhiG*Sin2T)*Esp2 - 
					(2*U*V*x2*PhiG*Sin2T)*Esp2 - 
					(U*x2*PsiG*Sin2T)*Esp2 + 
					(U2*x2*PsiG*Sin2T)*Esp2 - 
					(V2*x2*PsiG*Sin2T)*Esp2;

	DenRe=
		b - U*alpha*rho + PhiG - (b*CosT)*Esp + 
		(U*alpha*rho*CosT)*Esp + (PhiG*CosT)*Esp + 
		(V*alpha*rho*SinT)*Esp + (PsiG*SinT)*Esp;

	DenIm=
		-(V*alpha*rho) + PsiG + (V*alpha*rho*CosT)*Esp + 
		(PsiG*CosT)*Esp + (b*SinT)*Esp - (U*alpha*rho*SinT)*Esp -
		(PhiG*SinT)*Esp;

		*DDRePsi=NumRe/(DenRe*DenRe+DenIm*DenIm);
		*DDImPsi=NumIm/(DenRe*DenRe+DenIm*DenIm);

	/*

	This part is needed if one is interested in the case when Sigma(0) = Sigma(infinity)
	*/
		if (isVolInfFix) {

			a=b*x2;
			Phi3=T/alpha/alpha*(a*b-a*U*alpha*rho-a*PhiG);
			Psi3=T/alpha/alpha*(-a*V*alpha*rho-a*PsiG);
			Chi=PhiG*PhiG+PsiG*PsiG;
			Den=2*Chi;

			NumRe2=
				b*PhiG - U*alpha*rho*PhiG - Phi2G + 2*Chi - V*alpha*
				rho*PsiG - Psi2G -  (b*PhiG*CosT)*Esp + 
				(U*alpha*rho*PhiG*CosT)*Esp + (Phi2G*CosT)*Esp + (V*alpha*rho*
				PsiG*CosT)*Esp +  (Psi2G*CosT)*Esp + (V*alpha*rho*
				PhiG*SinT)*Esp +  (b*PsiG*SinT)*Esp - (U*alpha*rho*PsiG*SinT)*Esp;

			NumIm2=
				-(V*alpha*rho*PhiG) - b*PsiG + U*alpha*rho*PsiG + (V*
				alpha*rho*PhiG*CosT)*Esp + (b*PsiG*CosT)*Esp - 
				(U*alpha*rho*PsiG*CosT)*Esp + (b*PhiG*SinT)*Esp - 
				(U*alpha*rho*PhiG*SinT)*Esp - (Phi2G*SinT)*Esp - (V*alpha*rho*
				PsiG*SinT)*Esp -  (Psi2G*SinT)*Esp;

			Phi4=-a/alpha/alpha*log(NumRe2*NumRe2/Den/Den+NumIm2*NumIm2/Den/Den);
			Psi4=-2.*a/alpha/alpha*atan(NumIm2/NumRe2);

				*DDRePsi+=Phi3+Phi4;
				*DDImPsi+=Psi3+Psi4;
		}

		break;

	case VOLGA:


		x2=2.*sqrt(x2);

		NumRe=	-(b*U*x2) + b*U2*x2 - b*V2*x2 + U2*x2*alpha*rho - 
					U3*x2*alpha*rho + V2*x2*alpha*rho - 
					U*V2*x2*alpha*rho - U*x2*PhiG + U2*x2*PhiG - 
					V2*x2*PhiG - V*x2*PsiG + 2*U*V*x2*PsiG + 
					(2*b*U*x2*CosT)*Esp - 
					(2*b*U2*x2*CosT)*Esp + 
					(2*b*V2*x2*CosT)*Esp - 
					(2*U2*x2*alpha*rho*CosT)*Esp + 
					(2*U3*x2*alpha*rho*CosT)*Esp - 
					(2*V2*x2*alpha*rho*CosT)*Esp + 
					(2*U*V2*x2*alpha*rho*CosT)*Esp - 
					(b*U*x2*Cos2T)*Esp2 + 
					(b*U2*x2*Cos2T)*Esp2 - 
					(b*V2*x2*Cos2T)*Esp2 + 
					(U2*x2*alpha*rho*Cos2T)*Esp2 - 
					(U3*x2*alpha*rho*Cos2T)*Esp2 + 
					(V2*x2*alpha*rho*Cos2T)*Esp2 - 
					(U*V2*x2*alpha*rho*Cos2T)*Esp2 + 
					(U*x2*PhiG*Cos2T)*Esp2 - 
					(U2*x2*PhiG*Cos2T)*Esp2 + 
					(V2*x2*PhiG*Cos2T)*Esp2 + 
					(V*x2*PsiG*Cos2T)*Esp2 - 
					(2*U*V*x2*PsiG*Cos2T)*Esp2 + 
					(2*V*x2*PhiG*SinT)*Esp - 
					(4*U*V*x2*PhiG*SinT)*Esp - 
					(2*U*x2*PsiG*SinT)*Esp + 
					(2*U2*x2*PsiG*SinT)*Esp - 
					(2*V2*x2*PsiG*SinT)*Esp - 
					(b*U*x2*Sin2T)*Esp2 + 
					(b*U2*x2*Sin2T)*Esp2 - 
					(b*V2*x2*Sin2T)*Esp2 + 
					(U2*x2*alpha*rho*Sin2T)*Esp2 - 
					(U3*x2*alpha*rho*Sin2T)*Esp2 + 
					(V2*x2*alpha*rho*Sin2T)*Esp2 - 
					(U*V2*x2*alpha*rho*Sin2T)*Esp2 + 
					(U*x2*PhiG*Sin2T)*Esp2 - 
					(U2*x2*PhiG*Sin2T)*Esp2 + 
					(V2*x2*PhiG*Sin2T)*Esp2 + 
					(V*x2*PsiG*Sin2T)*Esp2 - 
					(2*U*V*x2*PsiG*Sin2T)*Esp2;

		NumIm=   -(b*V*x2) + 2*b*U*V*x2 - U2*V*x2*alpha*rho -
					V3*x2*alpha*rho - V*x2*PhiG + 2*U*V*x2*PhiG + 
					U*x2*PsiG - U2*x2*PsiG + V2*x2*PsiG + (2*b*V*x2*CosT)*Esp - 
					(4*b*U*V*x2*CosT)*Esp + (2*U2*V*x2*alpha*rho*CosT)*Esp +
					(2*V3*x2*alpha*rho*CosT)*Esp - 
					(b*V*x2*Cos2T)*Esp2 +  
					(2*b*U*V*x2*Cos2T)*Esp2 - 
					(U2*V*x2*alpha*rho*Cos2T)*Esp2 - 
					(V3*x2*alpha*rho*Cos2T)*Esp2 + 
					(V*x2*PhiG*Cos2T)*Esp2 - 
					(2*U*V*x2*PhiG*Cos2T)*Esp2 - 
					(U*x2*PsiG*Cos2T)*Esp2 + 
					(U2*x2*PsiG*Cos2T)*Esp2 - 
					(V2*x2*PsiG*Cos2T)*Esp2 - 
					(2*U*x2*PhiG*SinT)*Esp + 
					(2*U2*x2*PhiG*SinT)*Esp - 
					(2*V2*x2*PhiG*SinT)*Esp - 
					(2*V*x2*PsiG*SinT)*Esp + 
					(4*U*V*x2*PsiG*SinT)*Esp - 
					(b*V*x2*Sin2T)*Esp2 + 
					(2*b*U*V*x2*Sin2T)*Esp2 - 
					(U2*V*x2*alpha*rho*Sin2T)*Esp2 - 
					(V3*x2*alpha*rho*Sin2T)*Esp2 + 
					(V*x2*PhiG*Sin2T)*Esp2 - 
					(2*U*V*x2*PhiG*Sin2T)*Esp2 - 
					(U*x2*PsiG*Sin2T)*Esp2 + 
					(U2*x2*PsiG*Sin2T)*Esp2 - 
					(V2*x2*PsiG*Sin2T)*Esp2;

	DenRe=
		b - U*alpha*rho + PhiG - (b*CosT)*Esp + 
		(U*alpha*rho*CosT)*Esp + (PhiG*CosT)*Esp + 
		(V*alpha*rho*SinT)*Esp + (PsiG*SinT)*Esp;

	DenIm=
		-(V*alpha*rho) + PsiG + (V*alpha*rho*CosT)*Esp + 
		(PsiG*CosT)*Esp + (b*SinT)*Esp - (U*alpha*rho*SinT)*Esp -
		(PhiG*SinT)*Esp;

		*DRePsi=NumRe/(DenRe*DenRe+DenIm*DenIm);
		*DImPsi=NumIm/(DenRe*DenRe+DenIm*DenIm);

	/*
	This part is needed if one is interested in the case when Sigma(0) = Sigma(infinity)
	*/
		if (isVolInfFix) {

			a=b*x2;
			Phi3=T/alpha/alpha*(a*b-a*U*alpha*rho-a*PhiG);
			Psi3=T/alpha/alpha*(-a*V*alpha*rho-a*PsiG);

			Chi=PhiG*PhiG+PsiG*PsiG;
			Den=2*Chi;

			NumRe2=
				b*PhiG - U*alpha*rho*PhiG - Phi2G + 2*Chi - V*alpha*
				rho*PsiG - Psi2G -  (b*PhiG*CosT)*Esp + 
				(U*alpha*rho*PhiG*CosT)*Esp + (Phi2G*CosT)*Esp + (V*alpha*rho*
				PsiG*CosT)*Esp +  (Psi2G*CosT)*Esp + (V*alpha*rho*
				PhiG*SinT)*Esp +  (b*PsiG*SinT)*Esp - (U*alpha*rho*PsiG*SinT)*Esp;

			NumIm2=
				-(V*alpha*rho*PhiG) - b*PsiG + U*alpha*rho*PsiG + (V*
				alpha*rho*PhiG*CosT)*Esp + (b*PsiG*CosT)*Esp - 
				(U*alpha*rho*PsiG*CosT)*Esp + (b*PhiG*SinT)*Esp - 
				(U*alpha*rho*PhiG*SinT)*Esp - (Phi2G*SinT)*Esp - (V*alpha*rho*
				PsiG*SinT)*Esp -  (Psi2G*SinT)*Esp;

			Phi4=-a/alpha/alpha*log(NumRe2*NumRe2/Den/Den+NumIm2*NumIm2/Den/Den);
			Psi4=-2.*a/alpha/alpha*atan(NumIm2/NumRe2);

		*DRePsi+=Phi3+Phi4;
		*DImPsi+=Psi3+Psi4;

		}


		x2=2.;
		NumRe=	-(b*U*x2) + b*U2*x2 - b*V2*x2 + U2*x2*alpha*rho - 
					U3*x2*alpha*rho + V2*x2*alpha*rho - 
					U*V2*x2*alpha*rho - U*x2*PhiG + U2*x2*PhiG - 
					V2*x2*PhiG - V*x2*PsiG + 2*U*V*x2*PsiG + 
					(2*b*U*x2*CosT)*Esp - 
					(2*b*U2*x2*CosT)*Esp + 
					(2*b*V2*x2*CosT)*Esp - 
					(2*U2*x2*alpha*rho*CosT)*Esp + 
					(2*U3*x2*alpha*rho*CosT)*Esp - 
					(2*V2*x2*alpha*rho*CosT)*Esp + 
					(2*U*V2*x2*alpha*rho*CosT)*Esp - 
					(b*U*x2*Cos2T)*Esp2 + 
					(b*U2*x2*Cos2T)*Esp2 - 
					(b*V2*x2*Cos2T)*Esp2 + 
					(U2*x2*alpha*rho*Cos2T)*Esp2 - 
					(U3*x2*alpha*rho*Cos2T)*Esp2 + 
					(V2*x2*alpha*rho*Cos2T)*Esp2 - 
					(U*V2*x2*alpha*rho*Cos2T)*Esp2 + 
					(U*x2*PhiG*Cos2T)*Esp2 - 
					(U2*x2*PhiG*Cos2T)*Esp2 + 
					(V2*x2*PhiG*Cos2T)*Esp2 + 
					(V*x2*PsiG*Cos2T)*Esp2 - 
					(2*U*V*x2*PsiG*Cos2T)*Esp2 + 
					(2*V*x2*PhiG*SinT)*Esp - 
					(4*U*V*x2*PhiG*SinT)*Esp - 
					(2*U*x2*PsiG*SinT)*Esp + 
					(2*U2*x2*PsiG*SinT)*Esp - 
					(2*V2*x2*PsiG*SinT)*Esp - 
					(b*U*x2*Sin2T)*Esp2 + 
					(b*U2*x2*Sin2T)*Esp2 - 
					(b*V2*x2*Sin2T)*Esp2 + 
					(U2*x2*alpha*rho*Sin2T)*Esp2 - 
					(U3*x2*alpha*rho*Sin2T)*Esp2 + 
					(V2*x2*alpha*rho*Sin2T)*Esp2 - 
					(U*V2*x2*alpha*rho*Sin2T)*Esp2 + 
					(U*x2*PhiG*Sin2T)*Esp2 - 
					(U2*x2*PhiG*Sin2T)*Esp2 + 
					(V2*x2*PhiG*Sin2T)*Esp2 + 
					(V*x2*PsiG*Sin2T)*Esp2 - 
					(2*U*V*x2*PsiG*Sin2T)*Esp2;

		NumIm=   -(b*V*x2) + 2*b*U*V*x2 - U2*V*x2*alpha*rho -
					V3*x2*alpha*rho - V*x2*PhiG + 2*U*V*x2*PhiG + 
					U*x2*PsiG - U2*x2*PsiG + V2*x2*PsiG + (2*b*V*x2*CosT)*Esp - 
					(4*b*U*V*x2*CosT)*Esp + (2*U2*V*x2*alpha*rho*CosT)*Esp +
					(2*V3*x2*alpha*rho*CosT)*Esp - 
					(b*V*x2*Cos2T)*Esp2 +  
					(2*b*U*V*x2*Cos2T)*Esp2 - 
					(U2*V*x2*alpha*rho*Cos2T)*Esp2 - 
					(V3*x2*alpha*rho*Cos2T)*Esp2 + 
					(V*x2*PhiG*Cos2T)*Esp2 - 
					(2*U*V*x2*PhiG*Cos2T)*Esp2 - 
					(U*x2*PsiG*Cos2T)*Esp2 + 
					(U2*x2*PsiG*Cos2T)*Esp2 - 
					(V2*x2*PsiG*Cos2T)*Esp2 - 
					(2*U*x2*PhiG*SinT)*Esp + 
					(2*U2*x2*PhiG*SinT)*Esp - 
					(2*V2*x2*PhiG*SinT)*Esp - 
					(2*V*x2*PsiG*SinT)*Esp + 
					(4*U*V*x2*PsiG*SinT)*Esp - 
					(b*V*x2*Sin2T)*Esp2 + 
					(2*b*U*V*x2*Sin2T)*Esp2 - 
					(U2*V*x2*alpha*rho*Sin2T)*Esp2 - 
					(V3*x2*alpha*rho*Sin2T)*Esp2 + 
					(V*x2*PhiG*Sin2T)*Esp2 - 
					(2*U*V*x2*PhiG*Sin2T)*Esp2 - 
					(U*x2*PsiG*Sin2T)*Esp2 + 
					(U2*x2*PsiG*Sin2T)*Esp2 - 
					(V2*x2*PsiG*Sin2T)*Esp2;

	DenRe=
		b - U*alpha*rho + PhiG - (b*CosT)*Esp + 
		(U*alpha*rho*CosT)*Esp + (PhiG*CosT)*Esp + 
		(V*alpha*rho*SinT)*Esp + (PsiG*SinT)*Esp;

	DenIm=
		-(V*alpha*rho) + PsiG + (V*alpha*rho*CosT)*Esp + 
		(PsiG*CosT)*Esp + (b*SinT)*Esp - (U*alpha*rho*SinT)*Esp -
		(PhiG*SinT)*Esp;

		*DDRePsi=NumRe/(DenRe*DenRe+DenIm*DenIm);
		*DDImPsi=NumIm/(DenRe*DenRe+DenIm*DenIm);

	/*
	This part is needed if one is interested in the case when Sigma(0) = Sigma(infinity)
	*/
		if (isVolInfFix) {

			a=2.*b;
			Phi3=T/alpha/alpha*(a*b-a*U*alpha*rho-a*PhiG);
			Psi3=T/alpha/alpha*(-a*V*alpha*rho-a*PsiG);
			Chi=PhiG*PhiG+PsiG*PsiG;
			Den=2*Chi;

			NumRe2=
				b*PhiG - U*alpha*rho*PhiG - Phi2G + 2*Chi - V*alpha*
				rho*PsiG - Psi2G -  (b*PhiG*CosT)*Esp + 
				(U*alpha*rho*PhiG*CosT)*Esp + (Phi2G*CosT)*Esp + (V*alpha*rho*
				PsiG*CosT)*Esp +  (Psi2G*CosT)*Esp + (V*alpha*rho*
				PhiG*SinT)*Esp +  (b*PsiG*SinT)*Esp - (U*alpha*rho*PsiG*SinT)*Esp;

			NumIm2=
				-(V*alpha*rho*PhiG) - b*PsiG + U*alpha*rho*PsiG + (V*
				alpha*rho*PhiG*CosT)*Esp + (b*PsiG*CosT)*Esp - 
				(U*alpha*rho*PsiG*CosT)*Esp + (b*PhiG*SinT)*Esp - 
				(U*alpha*rho*PhiG*SinT)*Esp - (Phi2G*SinT)*Esp - (V*alpha*rho*
				PsiG*SinT)*Esp -  (Psi2G*SinT)*Esp;

			Phi4=-a/alpha/alpha*log(NumRe2*NumRe2/Den/Den+NumIm2*NumIm2/Den/Den);
			Psi4=-2.*a/alpha/alpha*atan(NumIm2/NumRe2);

		
			*DDRePsi+=Phi3+Phi4;
			*DDImPsi+=Psi3+Psi4;
		}

		break;

		case THETA:

		/* L2 term contribution */

		NumRe=	-(b*U*x2) + b*U2*x2 - b*V2*x2 + U2*x2*alpha*rho - 
					U3*x2*alpha*rho + V2*x2*alpha*rho - 
					U*V2*x2*alpha*rho - U*x2*PhiG + U2*x2*PhiG - 
					V2*x2*PhiG - V*x2*PsiG + 2*U*V*x2*PsiG + 
					(2*b*U*x2*CosT)*Esp - 
					(2*b*U2*x2*CosT)*Esp + 
					(2*b*V2*x2*CosT)*Esp - 
					(2*U2*x2*alpha*rho*CosT)*Esp + 
					(2*U3*x2*alpha*rho*CosT)*Esp - 
					(2*V2*x2*alpha*rho*CosT)*Esp + 
					(2*U*V2*x2*alpha*rho*CosT)*Esp - 
					(b*U*x2*Cos2T)*Esp2 + 
					(b*U2*x2*Cos2T)*Esp2 - 
					(b*V2*x2*Cos2T)*Esp2 + 
					(U2*x2*alpha*rho*Cos2T)*Esp2 - 
					(U3*x2*alpha*rho*Cos2T)*Esp2 + 
					(V2*x2*alpha*rho*Cos2T)*Esp2 - 
					(U*V2*x2*alpha*rho*Cos2T)*Esp2 + 
					(U*x2*PhiG*Cos2T)*Esp2 - 
					(U2*x2*PhiG*Cos2T)*Esp2 + 
					(V2*x2*PhiG*Cos2T)*Esp2 + 
					(V*x2*PsiG*Cos2T)*Esp2 - 
					(2*U*V*x2*PsiG*Cos2T)*Esp2 + 
					(2*V*x2*PhiG*SinT)*Esp - 
					(4*U*V*x2*PhiG*SinT)*Esp - 
					(2*U*x2*PsiG*SinT)*Esp + 
					(2*U2*x2*PsiG*SinT)*Esp - 
					(2*V2*x2*PsiG*SinT)*Esp - 
					(b*U*x2*Sin2T)*Esp2 + 
					(b*U2*x2*Sin2T)*Esp2 - 
					(b*V2*x2*Sin2T)*Esp2 + 
					(U2*x2*alpha*rho*Sin2T)*Esp2 - 
					(U3*x2*alpha*rho*Sin2T)*Esp2 + 
					(V2*x2*alpha*rho*Sin2T)*Esp2 - 
					(U*V2*x2*alpha*rho*Sin2T)*Esp2 + 
					(U*x2*PhiG*Sin2T)*Esp2 - 
					(U2*x2*PhiG*Sin2T)*Esp2 + 
					(V2*x2*PhiG*Sin2T)*Esp2 + 
					(V*x2*PsiG*Sin2T)*Esp2 - 
					(2*U*V*x2*PsiG*Sin2T)*Esp2;

		NumIm=   -(b*V*x2) + 2*b*U*V*x2 - U2*V*x2*alpha*rho -
					V3*x2*alpha*rho - V*x2*PhiG + 2*U*V*x2*PhiG + 
					U*x2*PsiG - U2*x2*PsiG + V2*x2*PsiG + (2*b*V*x2*CosT)*Esp - 
					(4*b*U*V*x2*CosT)*Esp + (2*U2*V*x2*alpha*rho*CosT)*Esp +
					(2*V3*x2*alpha*rho*CosT)*Esp - 
					(b*V*x2*Cos2T)*Esp2 +  
					(2*b*U*V*x2*Cos2T)*Esp2 - 
					(U2*V*x2*alpha*rho*Cos2T)*Esp2 - 
					(V3*x2*alpha*rho*Cos2T)*Esp2 + 
					(V*x2*PhiG*Cos2T)*Esp2 - 
					(2*U*V*x2*PhiG*Cos2T)*Esp2 - 
					(U*x2*PsiG*Cos2T)*Esp2 + 
					(U2*x2*PsiG*Cos2T)*Esp2 - 
					(V2*x2*PsiG*Cos2T)*Esp2 - 
					(2*U*x2*PhiG*SinT)*Esp + 
					(2*U2*x2*PhiG*SinT)*Esp - 
					(2*V2*x2*PhiG*SinT)*Esp - 
					(2*V*x2*PsiG*SinT)*Esp + 
					(4*U*V*x2*PsiG*SinT)*Esp - 
					(b*V*x2*Sin2T)*Esp2 + 
					(2*b*U*V*x2*Sin2T)*Esp2 - 
					(U2*V*x2*alpha*rho*Sin2T)*Esp2 - 
					(V3*x2*alpha*rho*Sin2T)*Esp2 + 
					(V*x2*PhiG*Sin2T)*Esp2 - 
					(2*U*V*x2*PhiG*Sin2T)*Esp2 - 
					(U*x2*PsiG*Sin2T)*Esp2 + 
					(U2*x2*PsiG*Sin2T)*Esp2 - 
					(V2*x2*PsiG*Sin2T)*Esp2;

		DNumRe=	(2*x2*(PhiG*(b*U - b*U2 + b*V2 - U2*alpha*
				rho + U3*alpha*rho - V2*alpha*rho + 
				U*V2*alpha*rho - U*PhiG + U2*PhiG - 
				V2*PhiG - V*PsiG + 2*U*V*PsiG)*Esp2 + 

				(b*(-U + U2 - V2)*PhiG + 
				U2*alpha*rho*PhiG - U3*alpha*rho*PhiG + 
				V2*alpha*rho*PhiG - U*V2*alpha*rho*
				PhiG + V*PhiG*PsiG - 2*U*V*PhiG*PsiG - U*Psi2G + 
				U2*Psi2G - V2*Psi2G)*Esp*CosT + 

				Esp*(-(V*Phi2G) + 2*U*V*Phi2G - b*U*
				PsiG + b*U2*PsiG - b*V2*PsiG + 
				U2*alpha*rho*PsiG - U3*alpha*rho*
				PsiG + V2*alpha*rho*PsiG - 
				U*V2*alpha*rho*PsiG + U*PhiG*PsiG - U2*
				PhiG*PsiG + V2*PhiG*PsiG)*SinT)); /* OK*/

		DNumIm=	(2*x2*(PhiG*(b*V - 2*b*U*V + U2*V*alpha*rho + V3*
				alpha*rho - V*PhiG + 2*U*V*PhiG + U*PsiG - 
				U2*PsiG + V2*PsiG)*Esp2 + 
				
				Esp*(b*(-1 + 2*U)*V*PhiG - U2*V*alpha*rho*PhiG - 
				V3*alpha*rho*PhiG - U*PhiG*PsiG + 
				U2*PhiG*PsiG - V2*PhiG*PsiG - 
				V*Psi2G + 2*U*V*Psi2G)*CosT - 

				Esp*(-(U*Phi2G) + U2*Phi2G - 
				V2*Phi2G + b*V*PsiG - 2*b*U*V*PsiG + 
				U2*V*alpha*rho*PsiG + V3*alpha*rho*
				PsiG - V*PhiG*PsiG + 2*U*V*PhiG*PsiG)*SinT)); /* OK */

		DDen=	(-2*(PhiG*(b*b - 2*b*U*alpha*rho + U2*
				alpha*alpha*rho*rho + V2*alpha*alpha*rho*rho - 
				2*b*PhiG + 2*U*alpha*rho*PhiG + Phi2G + 2*V*
				alpha*rho*PsiG + Psi2G)*Esp2 + 

				Esp*(-(b*b*PhiG) + 2*b*U*alpha*rho*
				PhiG - U2*alpha*alpha*rho*rho*PhiG - 
				V2*alpha*alpha*rho*rho*PhiG + 
				Phi3G - 2*V*alpha*rho*PhiG*PsiG - 2*b*Psi2G + 
				2*U*alpha*rho*Psi2G + PhiG*Psi2G)*
				CosT + 

				Esp*(2*V*alpha*rho*Phi2G - b*b*
				PsiG + 2*b*U*alpha*rho*PsiG - 
				U2*alpha*alpha*rho*rho*PsiG - V2*
				alpha*alpha*rho*rho*PsiG + 2*b*PhiG*PsiG - 
				2*U*alpha*rho*PhiG*PsiG + Phi2G*PsiG + 
				Psi3G)*SinT));            /* OK */

	DenRe=
		b - U*alpha*rho + PhiG - (b*CosT)*Esp + 
		(U*alpha*rho*CosT)*Esp + (PhiG*CosT)*Esp + 
		(V*alpha*rho*SinT)*Esp + (PsiG*SinT)*Esp;

	DenIm=
		-(V*alpha*rho) + PsiG + (V*alpha*rho*CosT)*Esp + 
		(PsiG*CosT)*Esp + (b*SinT)*Esp - (U*alpha*rho*SinT)*Esp -
		(PhiG*SinT)*Esp;

	Den=DenRe*DenRe+DenIm*DenIm;

		Phi2=(DNumRe*Den-DDen*NumRe)/Den/Den;  /* OK */
		Psi2=(DNumIm*Den-DDen*NumIm)/Den/Den;  /* OK */

		/* L3 term contribution */

		Phi3=1./alpha/alpha*(r*alpha*alpha*(U-1)+a*b-a*U*alpha*rho-a*PhiG);
		Psi3=1./alpha/alpha*(r*V*alpha*alpha-a*V*alpha*rho-a*PsiG);


		/* L4 term contribution */

	Chi=PhiG*PhiG+PsiG*PsiG;
	Den=2*Chi;

	NumRe=
		b*PhiG - U*alpha*rho*PhiG - Phi2G + 2*Chi - V*alpha*
		rho*PsiG - Psi2G -  (b*PhiG*CosT)*Esp + 
		(U*alpha*rho*PhiG*CosT)*Esp + (Phi2G*CosT)*Esp + (V*alpha*rho*
		PsiG*CosT)*Esp +  (Psi2G*CosT)*Esp + (V*alpha*rho*
		PhiG*SinT)*Esp +  (b*PsiG*SinT)*Esp - (U*alpha*rho*PsiG*SinT)*Esp; /* OK */

	NumIm=
		-(V*alpha*rho*PhiG) - b*PsiG + U*alpha*rho*PsiG + (V*
		alpha*rho*PhiG*CosT)*Esp + (b*PsiG*CosT)*Esp - 
		(U*alpha*rho*PsiG*CosT)*Esp + (b*PhiG*SinT)*Esp - 
		(U*alpha*rho*PhiG*SinT)*Esp - (Phi2G*SinT)*Esp - (V*alpha*rho*
		PsiG*SinT)*Esp -  (Psi2G*SinT)*Esp;  /* OK */
		
	DLogArg=-(PhiG*(b*b - 2*b*U*alpha*rho + U2*alpha*alpha*
			rho*rho + V2*alpha*alpha*rho*rho - 2*b*PhiG + 
			2*U*alpha*rho*PhiG + Phi2G + 2*V*alpha*rho*
			PsiG + Psi2G)*Esp2 + 
			Esp*(-(b*b*PhiG) + 2*b*U*alpha*rho*PhiG -
			U2*alpha*alpha*rho*rho*PhiG - 
			V2*alpha*alpha*rho*rho*PhiG + Phi3G
			- 2*V*alpha*rho*PhiG*PsiG - 2*b*Psi2G + 
			2*U*alpha*rho*Psi2G + PhiG*Psi2G)*
			CosT +			
			Esp*(2*V*alpha*rho*Phi2G - b*b*
			PsiG + 2*b*U*alpha*rho*PsiG - 
			U2*alpha*alpha*rho*rho*PsiG - V2*
			alpha*alpha*rho*rho*PsiG + 2*b*PhiG*PsiG - 
			2*U*alpha*rho*PhiG*PsiG + Phi2G*PsiG + 
			Psi3G)*SinT)/
			(2.*(Phi2G + Psi2G)); /* OK */

		Phi4=-a/alpha/alpha/(NumRe*NumRe/Den/Den+NumIm*NumIm/Den/Den)*DLogArg;
	
		DNumRe=-((Phi2G + Psi2G)*((-b + U*alpha*rho + PhiG)*CosT + (V*alpha*rho + PsiG)*SinT))*Esp; /*OK */
		DNumIm=((Phi2G + Psi2G)*(-((V*alpha*rho + PsiG)*CosT) + (-b + U*alpha*rho + PhiG)*SinT))*Esp; /* OK */
		Psi4=(DNumIm*NumRe-DNumRe*NumIm)/NumRe/NumRe;
		Psi4=-2*a/alpha/alpha*Psi4*1./(1+NumIm*NumIm/NumRe/NumRe);

		/* computes the output */

		*DRePsi=Phi2+Phi3+Phi4;
		*DImPsi=Psi2+Psi3+Psi4;
		*DDRePsi=0.0;
		*DDImPsi=0.0;

		break;

	case DENSITY:

		*DRePsi=0.0;
		*DImPsi=0.0;
		*DDRePsi=0.0;
		*DDImPsi=0.0;

	break;

	}

	return err;
}	
	
/*-----------------------------------------------------------------------------------------------*/								
/*-----------------------------------------------------------------------------------------------*/	
/*-----------------------------------------------------------------------------------------------*/	

Err HestonATMVol(	
					double	Forward,
					double	Maturity,
					double	Sigma,
					double	Alpha,
					double  SigmaIfty,
					double	b,
					double  Gamma,
					double	Rho,
					double	Disc,
					double  UpperBound,
					int		nSteps,
					SRT_Boolean isVolInfFix,
					double 	*result)
{
	Err err=NULL;
	int nStrikes = 1;
	double *Strike;
	double *Price;
	double vol;
	int IntegrType=3;

	Strike = dvector(0,1);
	Price = dvector(0,1);
	Strike[1] = Forward;

	if (isVolInfFix) SigmaIfty=Sigma;

	err = HestonPrice(	
						Forward,
						Strike,
						nStrikes,
						Maturity,
						Sigma,
						Alpha,
						SigmaIfty,
						b,
						Gamma,
						Rho,
						Disc,
						UpperBound,
						SRT_CALL, 
						PREMIUM,
						isVolInfFix,
						IntegrType,
						nSteps,
					 	Price);

	err = srt_f_optimpvol	(Price[1],Forward,Forward,Maturity,1,SRT_CALL,SRT_LOGNORMAL,&vol);
	*result = vol;
	free_dvector(Strike,0,1);
	free_dvector(Price,0,1);

	return err;
}

/*-----------------------------------------------------------------------------------------------*/								
/*-----------------------------------------------------------------------------------------------*/	
/*-----------------------------------------------------------------------------------------------*/	
/* Calibration of the sigma parameter in Heston to the ATM vol (shift is given) */

Err HestonCalibrateSigmaInit(
							 double	Forward,
							 double	Maturity,
							 double	ATMVol,
							 double	AlphaHeston,
							 double SigmaIftyHeston,
							 double	LambdaHeston,
							 double ShiftHeston,
							 double	RhoHeston,
							 double	UpperBound,
							 int	nSteps,
							 SRT_Boolean isVolInfFix,
							 double *result)
{
	
/*	double	SigmaMin,SigmaMax,SigmaMid;

	double	VolMin,VolMax,VolMid; */
	int		n_iter_max = 50;
	int		current_iter;
	Err		err= NULL;
	double	SigmaGuess,VolGuess;
	double	Sigma1,Vol1,Sigma2,Vol2,Vol1Shift;
	double	errorv;
	double	prec=1.e-7;
	double	shift;
	double	der;

	/* First Guess */
	SigmaGuess = ATMVol*Forward/(Forward+ShiftHeston);
	if (isVolInfFix) SigmaIftyHeston = SigmaGuess;

	err = HestonATMVol( Forward,
						Maturity,
						SigmaGuess,
						AlphaHeston,
						SigmaIftyHeston,
						LambdaHeston,
						ShiftHeston,
						RhoHeston,
						1,
						UpperBound,
						nSteps,
						isVolInfFix,
						&VolGuess);
	if (err)
	{
		return err;
	}

	errorv = ATMVol - VolGuess;
	if (fabs(errorv)<prec)
	{
		*result = VolGuess;
		return err;
	}

	/* Calculates a second point with the derivative */
	Sigma1 = SigmaGuess;
	Vol1 = VolGuess;
	shift = 0.0005;

	if (errorv < 0)		shift *= -1;
	
	if (isVolInfFix) SigmaIftyHeston = Sigma1+shift;

	err = HestonATMVol(Forward,
						Maturity,
						Sigma1+shift,
						AlphaHeston,
						SigmaIftyHeston,
						LambdaHeston,
						ShiftHeston,
						RhoHeston,
						1,
						UpperBound,
						nSteps,
						isVolInfFix,
						&Vol1Shift);

	der = (Vol1Shift - Vol1)/shift;

	Sigma2 = Sigma1 + errorv/der; 

	current_iter = 1;

	/* Newton */
	while (current_iter <= n_iter_max) 
	{

	if (isVolInfFix) SigmaIftyHeston = Sigma2;

		err = HestonATMVol(Forward,
						Maturity,
						Sigma2,
						AlphaHeston,
						SigmaIftyHeston,
						LambdaHeston,
						ShiftHeston,
						RhoHeston,
						1,
						UpperBound,
						nSteps,
						isVolInfFix,
						&Vol2);

		if (err) return err;
	
		errorv = ATMVol - Vol2;
		if (fabs(errorv)<prec)
		{
			*result = Sigma2;
			return err;
		}
		der = (Vol2 - Vol1)/(Sigma2 - Sigma1);

		Sigma1 = Sigma2;
		Vol1 = Vol2;
		Sigma2 = Sigma1 +errorv/der;
		current_iter++;
	}
	
	*result = Sigma2;
	return err;


/* Dichotomy */
/* 

	SigmaMax = 2*ATMVol*Forward/(Forward+ShiftHeston);
	SigmaMin = 0.25*ATMVol*Forward/(Forward+ShiftHeston);
	SigmaMid = 0.5*(SigmaMax+SigmaMin);
	err = HestonATMVol(	
						Forward,
						Maturity,
						SigmaMax,
						AlphaHeston,
						LambdaHeston,
						ShiftHeston,
						RhoHeston,
						1,
						UpperBound,
						nSteps,
						&VolMax);
	err = HestonATMVol(	
						Forward,
						Maturity,
						SigmaMin,
						AlphaHeston,
						LambdaHeston,
						ShiftHeston,
						RhoHeston,
						1,
						UpperBound,
						nSteps,
						&VolMin);
	if (VolMin>ATMVol || VolMax < ATMVol)
	{
		err = "Dichotomy failed";
		return err;
	}
	else
	while (SigmaMax - SigmaMin > prec & current_iter < n_iter_max)
	{
		err = HestonATMVol(	
						Forward,
						Maturity,
						SigmaMid,
						AlphaHeston,
						LambdaHeston,
						ShiftHeston,
						RhoHeston,
						1,
						UpperBound,
						nSteps,
						&VolMid);
		if (VolMid < ATMVol)
		{
			SigmaMin=SigmaMid;
			SigmaMid = 0.5*(SigmaMin + SigmaMax);
		}
		else
		{
			SigmaMax = SigmaMid;
			SigmaMid = 0.5*(SigmaMin + SigmaMax);
		}
		current_iter+=1;
	}
	
	if (SigmaMax - SigmaMin>prec)
		err = "Dichotomy failed";
	else
		*result = SigmaMid;
	return err; */
}

/*-----------------------------------------------------------------------------------------------*/								
/*-----------------------------------------------------------------------------------------------*/	
/*-----------------------------------------------------------------------------------------------*/	


Err CalibrateHestonToSABR(
						  double	Forward,
						  double	Maturity,
						  double	SigmaBeta,
						  double	AlphaSABR,
						  double	BetaSABR,
						  double	RhoSABR,
						  double	*SigmaHeston,
						  double	*AlphaHeston,
						  double    *SigmaIftyHeston,
						  double	*LambdaHeston,
						  double	*ShiftHeston,
						  double	*RhoHeston,
						  double	UpperBound,
						  int		nSteps,
						  int		nStdDev,
						  SRT_Boolean isVolInfFix,
						  SrtCalibrationType CalibrType
						)
{
	Err err=NULL;
	double	ATMVol;
	double	slope,volshift;
	double	KMin,KMax,volKMin,volKMax;
	double	alphamin,alphamax,alphamid,precalpha;
	double	sigma0;
	double	*Strikes;
	double	*Prices;
	double	*Vols;
	double	conv,convSABR,dApproxStdev;
	double  x,y,z;
	int n_iter_max=100,n_iter_max2=20,i,IntegrType=0;
	int current_iter=0,nK=9,trial_iter=0;
	double f1,f2,x0,x1,x2,x3,*StrikeVec,xMINold,xMINnew,MINnew,tol=1.e-5,ax,bx,cx;
	double atm_BSSABR_vol,atm_BSHESTON_vol,gap;

	Strikes = dvector(1,2);
	Prices = dvector(1,2);
	Vols = dvector(1,2);

	/* Calculates the ATM Vol */
	err = srt_f_optsarbvol(
							Forward,
							Forward,
							Maturity,
							SigmaBeta,
							AlphaSABR,
							BetaSABR,
							RhoSABR,
							SRT_BETAVOL,
							SRT_LOGNORMAL,
							&ATMVol);

	/* 
	KMax = Forward*(1+nStdDev*ATMVol*sqrt(Maturity));
	if (nStdDev*Forward - KMax > 0.01)
	{
		KMin = nStdDev*Forward - KMax;
	}
	else
	{
		KMin = 0.01;
	}
	*/

	dApproxStdev = SigmaBeta / pow(Forward, 1.0 - BetaSABR) * sqrt(Maturity);

	KMin = Forward * exp( - 0.5 * dApproxStdev * dApproxStdev - nStdDev * dApproxStdev); // modifier....
	KMax = Forward * exp( - 0.5 * dApproxStdev * dApproxStdev + nStdDev * dApproxStdev);	

	Strikes[1] = KMin;
	Strikes[2] = KMax;


	/* Calculates the vols far from the money */
	err = srt_f_optsarbvol(Forward,KMin,Maturity,SigmaBeta,AlphaSABR,BetaSABR,RhoSABR,SRT_BETAVOL,SRT_LOGNORMAL,&volKMin);
	err = srt_f_optsarbvol(Forward,KMax,Maturity,SigmaBeta,AlphaSABR,BetaSABR,RhoSABR,SRT_BETAVOL,SRT_LOGNORMAL,&volKMax);

	x = volKMax + volKMin;
	y = x*0.5;
	z = y - ATMVol;
	 convSABR = z;
	
							
	/* Calculates the shift */
	slope = (BetaSABR-1)*ATMVol/Forward*(1-0.5/(1+(1-BetaSABR)*(1-BetaSABR)/24*ATMVol*ATMVol*Maturity));
	volshift = 2*slope*Forward+ATMVol;
	*ShiftHeston = ATMVol*Forward/volshift*(1-ATMVol*ATMVol*Maturity/24)/(1-volshift*volshift*Maturity/24)-Forward;


	/* Calculates the rho */
	*RhoHeston = RhoSABR;

	/* So far, the mean-reversion remains equal to the initial guess ... */
	/* ... *LambdaHeston = *LambdaHeston; */

	/* Calculates the strikes up to 2 std */
/*
	KMax = Forward*(1+2*ATMVol*sqrt(Maturity));
	if (2*Forward - KMax > 0.01)
	{
		KMin = 2*Forward - KMax;
	}
	else
	{
		KMin = 0.01;
	}

	Strikes[1] = KMin;
	Strikes[2] = KMax;

*/
	/******************************************************************************************************/
	/* Calibrates the alpha */
    /******************************************************************************************************/
	
	switch (CalibrType) {
	
	case MATCH_CONV:

	alphamin = 0.001;
	alphamax = max(2.,AlphaSABR);
	precalpha = 0.000001;	

	/* Dichotomy on the alpha */
	while ((alphamax-alphamin>precalpha) && (current_iter<n_iter_max))
	{
		alphamid = 0.5*(alphamin+alphamax);

		err = HestonCalibrateSigmaInit(
							 	Forward,
							 	Maturity,
							 	ATMVol,
							 	alphamid,
								*SigmaIftyHeston,
							 	*LambdaHeston,
								*ShiftHeston,
							 	*RhoHeston,
							 	UpperBound,
							 	nSteps,
								isVolInfFix,
								&sigma0);

		if (err) return err;

		err = HestonPrice(	
						Forward,
						Strikes,
						2,
						Maturity,
						sigma0,
						alphamid,
						*SigmaIftyHeston,
						*LambdaHeston,
						*ShiftHeston,
						*RhoHeston,
						1,
						UpperBound,
					 	SRT_CALL, 
						PREMIUM,
						isVolInfFix,
						IntegrType,
						nSteps,
					 	Prices
					);

		if (err) return err;

		err = srt_f_optimpvol	(Prices[1],Forward,Strikes[1],Maturity,1,SRT_CALL,SRT_LOGNORMAL,&(Vols[1]));
		err = srt_f_optimpvol	(Prices[2],Forward,Strikes[2],Maturity,1,SRT_CALL,SRT_LOGNORMAL,&(Vols[2]));

		conv = ((Vols[1]+Vols[2])*0.5-ATMVol);

		if (conv<convSABR)
			alphamin = alphamid;
		else
			alphamax = alphamid;
		current_iter+=1;

	}

	*AlphaHeston = alphamid;

	err = HestonCalibrateSigmaInit(
						 	Forward,
						 	Maturity,
						 	ATMVol,
						 	alphamid,
							*SigmaIftyHeston,
						 	*LambdaHeston,
							*ShiftHeston,
						 	*RhoHeston,
						 	UpperBound,
						 	nSteps,
							isVolInfFix,
							&sigma0);

	*SigmaHeston = sigma0;

	break;
	case CHI2_MIN:

	StrikeVec=dvector(1,nK);

	for (i=1;i<=nK;i++){
	
		StrikeVec[i]=KMin+(i-1.)*(KMax-KMin)/(nK-1.);
	}

	ax = 0.000001;
	cx = max(2.0,AlphaSABR);
	bx=0.1; // (ax+cx)/16.; //(ax+cx)/8.;
	precalpha = 0.0000000001;

		err = HestonCalibrateSigmaInit(
							 	Forward,
							 	Maturity,
							 	ATMVol,
							 	bx,
								*SigmaIftyHeston,
							 	*LambdaHeston,
								*ShiftHeston,
							 	*RhoHeston,
							 	UpperBound,
							 	nSteps,
								isVolInfFix,
								&sigma0);   /* 0.09 this is the first guess */
	    xMINnew=bx;

		do {


			xMINold=xMINnew;
			x0=ax;
			x3=cx;

			if (fabs(cx-bx) > fabs(bx-ax)) {
				x1=bx;
				x2=bx+C*(cx-bx);
			} else {
				x2=bx;
				x1=bx-C*(bx-ax);
			}

			err = srt_f_optsarbvol(Forward,Forward,Maturity,SigmaBeta,AlphaSABR,BetaSABR,RhoSABR,
								  SRT_BETAVOL,SRT_LOGNORMAL,&atm_BSSABR_vol);

			/**********************  computes the Chi2 at x1 ***********************/

			err = HestonATMVol(Forward,Maturity,sigma0,x1,*SigmaIftyHeston,*LambdaHeston,*ShiftHeston,
					            *RhoHeston,1.0,UpperBound,nSteps,isVolInfFix,&atm_BSHESTON_vol);

			gap = atm_BSSABR_vol-atm_BSHESTON_vol;

			f1=SmileChi2(Forward,Maturity,sigma0,x1,*SigmaIftyHeston,*LambdaHeston,*ShiftHeston,*RhoHeston,
				         UpperBound,nSteps,nK,StrikeVec,SigmaBeta,AlphaSABR,BetaSABR,RhoSABR,gap,isVolInfFix);


			/**********************  computes the Chi2 at x1 **********************/

			err = HestonATMVol(Forward,Maturity,sigma0,x2,*SigmaIftyHeston,*LambdaHeston,*ShiftHeston,
					            *RhoHeston,1.0,UpperBound,nSteps,isVolInfFix,&atm_BSHESTON_vol);

			gap = atm_BSSABR_vol-atm_BSHESTON_vol;

			f2=SmileChi2(Forward,Maturity,sigma0,x2,*SigmaIftyHeston,*LambdaHeston,*ShiftHeston,*RhoHeston,
				         UpperBound,nSteps,nK,StrikeVec,SigmaBeta,AlphaSABR,BetaSABR,RhoSABR,gap,isVolInfFix);


			while (fabs(x3-x0) > tol*(fabs(x1)+fabs(x2))) {

				if (f2 < f1) {

					SHFT3(x0,x1,x2,R*x1+C*x3);
					f1=f2;

					err = HestonATMVol(Forward,Maturity,sigma0,x2,*SigmaIftyHeston,*LambdaHeston,*ShiftHeston,
					            *RhoHeston,1.0,UpperBound,nSteps,isVolInfFix,&atm_BSHESTON_vol);

					gap = atm_BSSABR_vol-atm_BSHESTON_vol;

					f2=SmileChi2(Forward,Maturity,sigma0,x2,*SigmaIftyHeston,*LambdaHeston,*ShiftHeston,*RhoHeston,
				         UpperBound,nSteps,nK,StrikeVec,SigmaBeta,AlphaSABR,BetaSABR,RhoSABR,gap,isVolInfFix);

				} else {

					SHFT3(x3,x2,x1,R*x2+C*x0);
					f2=f1;

					err = HestonATMVol(Forward,Maturity,sigma0,x1,*SigmaIftyHeston,*LambdaHeston,*ShiftHeston,
					            *RhoHeston,1.0,UpperBound,nSteps,isVolInfFix,&atm_BSHESTON_vol);

					gap = atm_BSSABR_vol-atm_BSHESTON_vol;

					f1=SmileChi2(Forward,Maturity,sigma0,x1,*SigmaIftyHeston,*LambdaHeston,*ShiftHeston,*RhoHeston,
				         UpperBound,nSteps,nK,StrikeVec,SigmaBeta,AlphaSABR,BetaSABR,RhoSABR,gap,isVolInfFix);
				}
			}

			if (f1 < f2) {
				xMINnew=x1;
				MINnew=f1;
			} else {
				xMINnew=x2;
				MINnew=f2;
			}
			

		err = HestonCalibrateSigmaInit(
							 	Forward,
							 	Maturity,
							 	ATMVol,
							 	xMINnew,
								*SigmaIftyHeston,
							 	*LambdaHeston,
								*ShiftHeston,
							 	*RhoHeston,
							 	UpperBound,
							 	nSteps,
								isVolInfFix,
								&sigma0);

		if (sigma0 > 0.15) /* > 0.2 */ {
			/* sigma0=0.09; */

		bx=bx+pow(-1.0,trial_iter)*0.01*(trial_iter+1);
		err = HestonCalibrateSigmaInit(
							 	Forward,
							 	Maturity,
							 	ATMVol,
							 	bx,
								*SigmaIftyHeston,
							 	*LambdaHeston,
								*ShiftHeston,
							 	*RhoHeston,
							 	UpperBound,
							 	nSteps,
								isVolInfFix,
								&sigma0);

			current_iter=-1;
			trial_iter++;
			if (bx < 0.0) current_iter = n_iter_max2;
		}

		current_iter+=1;

		} 	while ((xMINnew-xMINold > precalpha) && (current_iter < n_iter_max2+1));

		*AlphaHeston = xMINnew;
		err = HestonCalibrateSigmaInit(
							 	Forward,
							 	Maturity,
							 	ATMVol,
							 	xMINnew,
								*SigmaIftyHeston,
							 	*LambdaHeston,
								*ShiftHeston,
							 	*RhoHeston,
							 	UpperBound,
							 	nSteps,
								isVolInfFix,
								&sigma0);
		*SigmaHeston = sigma0;

		free_dvector(StrikeVec,1,11);
		break;
	}

	free_dvector(Strikes,1,2);
	free_dvector(Prices,1,2);
	free_dvector(Vols,1,2);
	return err;

}

/*-----------------------------------------------------------------------------------------------*/								
/*            This function computes the chi2 error on the smile								 */
/*-----------------------------------------------------------------------------------------------*/

double  SmileChi2(   double	Forward,
				     double	Maturity,
					 double	Sigma0,
					 double	AlphaHeston,
					 double SigmaIftyHeston,
					 double	LambdaHeston,
					 double ShiftHeston,
					 double	RhoHeston,
					 double	UpperBound,
					 int	nSteps,
					 int	nK,
					 double *StrikeVec,
					 double SigmaBeta,
					 double AlphaSABR,
					 double BetaSABR,
					 double RhoSABR,
					 double gap,
					 SRT_Boolean isVolInfFix
					 )
{
		double BS_SABRVol,*BS_HESTONVol,Chi2=0.0;
		Err err = NULL;
		int i;

		BS_HESTONVol = dvector(1,nK);

		err = HestonVol( Forward,
						 StrikeVec,
						 nK,
						 Maturity,
						 Sigma0,
						 SigmaIftyHeston,
						 AlphaHeston,
						 ShiftHeston,
						 RhoHeston,
						 LambdaHeston,
						 HESTON_TO_LOG,			
						 UpperBound,
						 nSteps,
						 isVolInfFix,
						 BS_HESTONVol
						 );

		for (i=1;i<=nK;i++){

			err = srt_f_optsarbvol(	Forward,
									StrikeVec[i],
									Maturity,
									SigmaBeta,
									AlphaSABR,
									BetaSABR,
									RhoSABR,
									SRT_BETAVOL,
									SRT_LOGNORMAL,
									&BS_SABRVol);

		Chi2+=(BS_HESTONVol[i]+gap-BS_SABRVol)*(BS_HESTONVol[i]+gap-BS_SABRVol);
		}

		free_dvector(BS_HESTONVol,1,nK);
		return Chi2;
}

/*-----------------------------------------------------------------------------------------------*/								
/*            This function returns the shift on the Fwd. lognormal dynamics that                */
/*            gives the same skew as an equivalent Beta for the SABR model	                     */
/*-----------------------------------------------------------------------------------------------*/

Err  HestonGetShift(	
					double 		Forward,
					double 		Maturity,
					double 		SigmaBetaSABR,
					double		AlphaSABR,
					double		BetaSABR,
					double		RhoSABR,
					double		ATMVol,
					double		*shift)
{

	double slope,volshift;
	Err err=NULL;

	if (SigmaBetaSABR == 0.0) {

	}  else {

	/* Calculates the ATM Vol */
	err = srt_f_optsarbvol(
							Forward,
							Forward,
							Maturity,
							SigmaBetaSABR,
							AlphaSABR,
							BetaSABR,
							RhoSABR,
							SRT_BETAVOL,
							SRT_LOGNORMAL,
							&ATMVol);

		if (err) return err;

	}

	/* Calculates the shift */
	slope = (BetaSABR-1)*ATMVol/Forward*(1-0.5/(1+(1-BetaSABR)*(1-BetaSABR)/24*ATMVol*ATMVol*Maturity));
	volshift = 2*slope*Forward+ATMVol;
	*shift = ATMVol*Forward/volshift*(1.-ATMVol*ATMVol*Maturity/24.)/(1.-volshift*volshift*Maturity/24.)-Forward;

		
	return err;

}

/*-----------------------------------------------------------------------------------------------*/	
/*        This function converts among different types of volatilities by matching the premiums  */
/*        of the corresponding calls/puts. It is the equivalent of the add-in SABRVol adapted    */
/*        for the Heston case                                                                    */
/*-----------------------------------------------------------------------------------------------*/	

#define MAXIT 100
Err  HestonVol(	
				double Forward,
				double *Strikes,
				int	   nStrikes,
				double Maturity,
				double Vol,
				double VolInfty,
				double Alpha,
				double Gamma,
				double Rho,
				double MeanRev,
				SrtVolConversion  TypeConversion,			
				double UpperBound,
				int nSteps,
				SRT_Boolean isVolInfFix,
				double *NewVols
			 )
{

	Err err = NULL;
	int i,j,IntegrType=3; /* sets up the integration type to the "mixed" one */
	double df,dx,dxold,f,fh,fl,x1=0.0001,x2=2.0;
	double temp,xh,xl,rts,xacc=1.e-7;
	double BSPremium,/*HestonPremium,*/ImplVol;
	double *result=NULL,*result_low=NULL,*result_high=NULL,*Strike=NULL;

	result=dvector(1,nStrikes);
	result_low=dvector(1,nStrikes);
	result_high=dvector(1,nStrikes);
	Strike=dvector(1,1);

	switch (TypeConversion){

		case LOG_TO_LOG:

		break;

		case LOG_TO_HESTON:

			/* computes premium at the left boundary */
			err =  HestonPrice(	
								Forward,
								Strikes,
								1,
								Maturity,
								x1,
								Alpha,
								VolInfty, 
								MeanRev,
								Gamma,
								Rho,
								1.0,
								UpperBound,
								SRT_CALL, 
								PREMIUM,
								isVolInfFix,
								IntegrType,
								nSteps,
								result_low
							  );

			/* computes premium at the right boundary */
			err =  HestonPrice(	
								Forward,
								Strikes,
								1,
								Maturity,
								x2,
								Alpha,
								VolInfty, 
								MeanRev,
								Gamma,
								Rho,
								1.0,
								UpperBound,
								SRT_CALL, 
								PREMIUM,
								isVolInfFix,
								IntegrType,
								nSteps,
								result_high
							  );


			for (i=1;i<=nStrikes;i++) {

				BSPremium = srt_f_optblksch(	
											Forward,
											Strikes[i],
											Vol,
											Maturity, 
											1.0, 
											SRT_CALL, 
											PREMIUM
									);


				fl=result_low[i]-BSPremium;
				fh=result_high[i]-BSPremium;

				if ((fl > 0.0 && fh > 0.0) || (fl < 0.0 && fh < 0.0))  /* Throw exception: not bracketed */ ;

				if (fl == 0.0) {
					NewVols[i]=x1;
					goto EndLoop;
					// return err;
				}
				if (fh == 0.0) {
					NewVols[i]=x2;
					goto EndLoop;
					// return err;
				}
				if (fl < 0.0) {
					xl=x1;
					xh=x2;
				} else {
					xh=x1;
					xl=x2;
				}

				rts=0.5*(x1+x2);
				dxold=fabs(x2-x1);
				dx=dxold;
				
				Strike[1]=Strikes[i];

				err =  HestonPrice(	
									Forward,
									Strike,
									1,
									Maturity,
									rts,
									Alpha,
									VolInfty, 
									MeanRev,
									Gamma,
									Rho,
									1.0,
									UpperBound,
									SRT_CALL, 
									PREMIUM,
									isVolInfFix,
									IntegrType,
									nSteps,
									result
								  );

				f=result[1]-BSPremium;

				err =  HestonPrice(	
									Forward,
									Strike,
									1,
									Maturity,
									rts,
									Alpha,
									VolInfty, 
									MeanRev,
									Gamma,
									Rho,
									1.0,
									UpperBound,
									SRT_CALL, 
									VEGA,
									isVolInfFix,
									IntegrType,
									nSteps,
									result
								  );
				df=result[1];

					for (j=1;j<=MAXIT;j++) {
						if ((((rts-xh)*df-f)*((rts-xl)*df-f) >= 0.0)
							|| (fabs(2.0*f) > fabs(dxold*df))) {
							dxold=dx;
							dx=0.5*(xh-xl);
							rts=xl+dx;				

							if (xl == rts) {
								NewVols[i]=rts;
								goto EndLoop;
								// return err;
							}

						} else {
							dxold=dx;
							dx=f/df;
							temp=rts;
							rts -= dx;		

							if (temp == rts) {
								NewVols[i]=rts;
								goto EndLoop;
								// return err;
							}

						}
						if (fabs(dx) < xacc) {
							NewVols[i]=rts;
							goto EndLoop;
							// return err;
						}	

						err =  HestonPrice(	
											Forward,
											Strike,
											1,
											Maturity,
											rts,
											Alpha,
											VolInfty, 
											MeanRev,
											Gamma,
											Rho,
											1.0,
											UpperBound,
											SRT_CALL, 
											PREMIUM,
											isVolInfFix,
											IntegrType,
											nSteps,
											result
										  );

						f=result[1]-BSPremium;	

						if (f < 0.0)
							xl=rts;
						else
							xh=rts;
					}		

EndLoop: ;
		}			

		  return err;


		break;

		case HESTON_TO_LOG:

			err =  HestonPrice(	
								Forward,
								Strikes,
								nStrikes,
								Maturity,
								Vol,
								Alpha,
								VolInfty, 
								MeanRev,
								Gamma,
								Rho,
								1.0,
								UpperBound,
								SRT_CALL, 
								PREMIUM,
								isVolInfFix,
								IntegrType,
								nSteps,
								result
							  );

			for (j=1;j<=nStrikes;j++) {

				err = srt_f_optimpvol(result[j],Forward, Strikes[j], Maturity, 1.0, 
											  SRT_CALL, SRT_LOGNORMAL, &ImplVol);

				NewVols[j]=ImplVol;
			}

		break;

		case HESTON_TO_HESTON:

		break;

		case HESTON_TO_NORMAL:

			err =  HestonPrice(	
								Forward,
								Strikes,
								nStrikes,
								Maturity,
								Vol,
								Alpha,
								VolInfty, 
								MeanRev,
								Gamma,
								Rho,
								1.0,
								UpperBound,
								SRT_CALL, 
								PREMIUM,
								isVolInfFix,
								IntegrType,
								nSteps,
								result
							  );

			for (j=1;j<=nStrikes;j++) {

				err = srt_f_optimpvol(result[j],Forward, Strikes[j], Maturity, 1.0, 
											  SRT_CALL, SRT_NORMAL, &ImplVol);

				NewVols[j]=ImplVol;
			}


		break;

		case NORMAL_TO_HESTON:

		
			err =  HestonPrice(	
								Forward,
								Strikes,
								1,
								Maturity,
								x1,
								Alpha,
								VolInfty, 
								MeanRev,
								Gamma,
								Rho,
								1.0,
								UpperBound,
								SRT_CALL, 
								PREMIUM,
								isVolInfFix,
								IntegrType,
								nSteps,
								result_low
							  );

			err =  HestonPrice(	
								Forward,
								Strikes,
								1,
								Maturity,
								x2,
								Alpha,
								VolInfty, 
								MeanRev,
								Gamma,
								Rho,
								1.0,
								UpperBound,
								SRT_CALL, 
								PREMIUM,
								isVolInfFix,
								IntegrType,
								nSteps,
								result_high
							  );


			for (i=1;i<=nStrikes;i++) {

				BSPremium = srt_f_optblknrm(	
											Forward,
											Strikes[i],
											Vol,
											Maturity, 
											1.0, 
											SRT_CALL, 
											PREMIUM
									);

				fl=result_low[i]-BSPremium;
				fh=result_high[i]-BSPremium;

				if ((fl > 0.0 && fh > 0.0) || (fl < 0.0 && fh < 0.0))   ;

				if (fl == 0.0) {
					NewVols[i]=x1;
					goto EndLoop2;
					// return err;
				}
				if (fh == 0.0) {
					NewVols[i]=x2;
					goto EndLoop2;
					// return err;
				}
				if (fl < 0.0) {
					xl=x1;
					xh=x2;
				} else {
					xh=x1;
					xl=x2;
				}

				rts=0.5*(x1+x2);
				dxold=fabs(x2-x1);
				dx=dxold;
				
				Strike[1]=Strikes[i];

				err =  HestonPrice(	
									Forward,
									Strike,
									1,
									Maturity,
									rts,
									Alpha,
									VolInfty, 
									MeanRev,
									Gamma,
									Rho,
									1.0,
									UpperBound,
									SRT_CALL, 
									PREMIUM,
									isVolInfFix,
									IntegrType,
									nSteps,
									result
								  );

				f=result[1]-BSPremium;

				err =  HestonPrice(	
									Forward,
									Strike,
									1,
									Maturity,
									rts,
									Alpha,
									VolInfty, 
									MeanRev,
									Gamma,
									Rho,
									1.0,
									UpperBound,
									SRT_CALL, 
									VEGA,
									isVolInfFix,
									IntegrType,
									nSteps,
									result
								  );
				df=result[1];

					for (j=1;j<=MAXIT;j++) {
						if ((((rts-xh)*df-f)*((rts-xl)*df-f) >= 0.0)
							|| (fabs(2.0*f) > fabs(dxold*df))) {
							dxold=dx;
							dx=0.5*(xh-xl);
							rts=xl+dx;				

							if (xl == rts) {
								NewVols[i]=rts;
								goto EndLoop;
								// return err;
							}

						} else {
							dxold=dx;
							dx=f/df;
							temp=rts;
							rts -= dx;		

							if (temp == rts) {
								NewVols[i]=rts;
								goto EndLoop2;
								// return err;
							}

						}
						if (fabs(dx) < xacc) {
							NewVols[i]=rts;
							goto EndLoop2;
							// return err;
						}	

						err =  HestonPrice(	
											Forward,
											Strike,
											1,
											Maturity,
											rts,
											Alpha,
											VolInfty, 
											MeanRev,
											Gamma,
											Rho,
											1.0,
											UpperBound,
											SRT_CALL, 
											PREMIUM,
											isVolInfFix,
											IntegrType,
											nSteps,
											result
										  );

						f=result[1]-BSPremium;	

						if (f < 0.0)
							xl=rts;
						else
							xh=rts;
					}		

EndLoop2: ;
		}			

		  return err;

		break;

		case NORMAL_TO_LOG:


			for (j=1;j<=nStrikes;j++) {

				BSPremium = srt_f_optblknrm(Forward,Strikes[j],Vol,Maturity,1.0,SRT_CALL,PREMIUM);

				err = srt_f_optimpvol(BSPremium,Forward, Strikes[j], Maturity, 1.0, 
											  SRT_CALL, SRT_LOGNORMAL, &ImplVol);

				NewVols[j]=ImplVol;
			}

		break;

		case LOG_TO_NORMAL:


			for (j=1;j<=nStrikes;j++) {

				BSPremium = srt_f_optblksch(Forward,Strikes[j],Vol,Maturity,1.0,SRT_CALL,PREMIUM);

				err = srt_f_optimpvol(BSPremium,Forward, Strikes[j], Maturity, 1.0, 
											  SRT_CALL, SRT_NORMAL, &ImplVol);

				NewVols[j]=ImplVol;
			}

		break;

	}

	free_dvector(result,1,nStrikes);
	free_dvector(result_low,1,nStrikes);
	free_dvector(result_high,1,nStrikes);
	free_dvector(Strike,1,1);

	return err;
}

#undef MAXIT
#undef C
#undef R
#undef SHFT2
#undef SHFT3
