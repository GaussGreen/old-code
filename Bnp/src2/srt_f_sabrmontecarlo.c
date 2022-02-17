#include		<math.h"
#include        "utallhdr.h"
#include        <OPFNCTNS.H>
#include		<num_h_sobol.h"
#include		<srt_h_all.h"
#include		<num_h_abs.h"
#include		<srt_h_resetable.h"



Err monte_carlo_sabr(
						double **rndm_mat,
						double forw,
						double vovol,
						double beta,
						double rho,
						double sigma,
						double num_paths,
						double maturity,
						double num_steps,
						int	   path_num,
						double *Forward,
						int modeltype,
						int sampletype,
						int *pathnogood)


{
double g[2];
double a1;
double a2;
double Eta;
double delta_t;
double sqrt_delta_t;
double delta;
int		 step;
Err		err=NULL;
long     idum;
double rand1,rand2;
double statevar,statevar0,statevar1;
idum = -path_num;
uniform(&idum);

*pathnogood = 0;



/* use a1 and a2 as symmetrical correlated values for the Brownians */
	
	a1 = sqrt((1+rho)/2);
	a2 = sqrt((1-rho)/2);

	step = 1;
	delta_t = maturity/num_steps;
	sqrt_delta_t = sqrt(delta_t);
	Eta = sigma;
	Forward[0] = forw;
	Forward[1] = forw;
	if (beta!=1)
	{
		statevar = pow(forw, 2-2*beta)/pow(sigma*(1-beta),2);
		delta = (1-2*beta)/(1-beta);
	}

	


	while (step<num_steps)

	{



			if (sampletype == 0)/*Randsam*/
			{
				rand1= inv_cumnorm_fast(uniform(&idum));
				rand2 = inv_cumnorm_fast(uniform(&idum));
				g[0] = a1*rand1+a2*rand2;
				g[1] = a1*rand1-a2*rand2;
			}

			else
			{
				g[0] = a1*rndm_mat[1][step]+a2*rndm_mat[2][step];
				g[1] = a1*rndm_mat[1][step]-a2*rndm_mat[2][step];
			}

			switch (modeltype)
	
			{
			case 0:
			/* Euler Discretization*/
 
			Forward[0]+=Eta*pow(Forward[0],beta)*g[0]*sqrt_delta_t;
			Forward[1]+=sigma*pow(Forward[1],beta)*g[0]*sqrt_delta_t;
			break;

			case 1:

			/*bessel-square discretization*/
			statevar = pow(Forward[0],2-2*beta)/pow(Eta*(1-beta),2);
			statevar = statevar+delta*delta_t+2*sqrt(statevar)*g[0]*sqrt_delta_t;
			Forward[0] = pow(statevar*(1-beta)*(1-beta)*Eta*Eta,1/(2-2*beta));
			
			break;

			case 2:

			/* Bessel discretization*/
  			statevar0 = pow(Forward[0],1-beta)/(1-beta);
			statevar0 = statevar0+Eta*g[0]*sqrt_delta_t-beta*Eta*Eta*delta_t/(2*(1-beta)*statevar0);
			statevar1 = pow(Forward[1],1-beta)/(1-beta);
			statevar1 = statevar1+sigma*g[0]*sqrt_delta_t-beta*sigma*sigma*delta_t/(2*(1-beta)*statevar1);
			Forward[0] = pow((1-beta)*statevar0,1/(1-beta));
			Forward[1] = pow((1-beta)*statevar1,1/(1-beta));
			break;

			case 3:

			/* Mihlstein discretization */		

				Forward[0] += pow(Forward[0],beta)*(Eta*g[0]*sqrt_delta_t)+0.5*Eta*Eta*beta*pow(Forward[0],2*beta-1)*(g[0]*g[0]-1)*delta_t;
				Forward[1] += pow(Forward[1],beta)*(sigma*g[0]*sqrt_delta_t)+0.5*sigma*sigma*beta*pow(Forward[1],2*beta-1)*(g[0]*g[0]-1)*delta_t;
				//+jumps[step]/sigma+(U2*lambda2-U1*lambda1)*delta_t);
			break;
			}

			Eta = Eta*exp(vovol*g[1]*sqrt_delta_t-0.5*vovol*vovol*delta_t);
			
			if (Forward[0] < 0)
			{
/*				*pathnogood = 1;*/
				Forward[0] = 0;
			}

		step++;
	}	


/* free memory*/


  return NULL;
}





Err srt_f_sabrmontecarlo(
					double forw,
					double vovol,
					double beta,
					double rho,
					double num_paths,
					double num_steps,
					double sigma,
					double maturity,
					double *strike,
					int    num_strikes,
					int modeltype, 
					int sampletype,/*0: randsam, 1: abs, 2:sobol, 3: SPECTRUNC*/
					double *impvol
					)
					
{
double ***cube=NULL;
int cur_path;
double *Forward;
double **price=NULL;
double *times_at_steps;
int i;
long seed;
double *blkbeta;
int pathnogood;
Err		err=NULL;
double actnumpaths = num_paths;
double *correl = NULL;
double **stdev= NULL;
double correlation = 0;

/*Memory allocation*/
price = dmatrix(1,2,1,num_strikes);
stdev = dmatrix(1,2,1,num_strikes);
correl = dvector(1,num_strikes);
cube = f3tensor(1,(long) num_paths,1,2,1,(long) num_steps); 
Forward = dvector(1,2);
blkbeta = dvector(1,num_strikes);
times_at_steps = dvector(1,(long) num_steps);
/* Choice of a method */

for (i = 1;i<num_steps+1;i++)
{
	times_at_steps[i] = i*maturity/num_steps;
}

seed = -123456789;

/* Spectrunc*/
if (sampletype == 3)
{
	if (err = SpecTruncCube(cube,times_at_steps,0.0001,1,(long) num_paths,1,2,1,(long) num_steps,&seed))
	{
		free_dmatrix(price,1,2,1,num_strikes);
		free_f3tensor(cube,1,(long) num_paths,1,2,1,(long) num_steps);
		free_dvector(Forward,1,2);
		free_dvector(blkbeta,1,num_strikes);
		smessage("error in sobol cube");
		return err;
	}
}

/*Sobol*/
if (sampletype==2)
{
	if (err = sobol_init(1,(long) num_paths,1,2,1,(long) num_steps))
	{
		free_dmatrix(price,1,2,1,num_strikes);
		free_f3tensor(cube,1,(long) num_paths,1,2,1,(long) num_steps);
		free_dvector(Forward,1,2);
		free_dvector(blkbeta,1,num_strikes);
		smessage("error in sobol cube");
		return err;
	}



	if (err = sobol_cube(cube,1,(long) num_paths,1,2,1,(long) num_steps))
	{
		free_dmatrix(price,1,2,1,num_strikes);
		free_f3tensor(cube,1,(long) num_paths,1,2,1,(long) num_steps);
		free_dvector(Forward,1,2);
		free_dvector(blkbeta,1,num_strikes);
		smessage("error in sobol cube");
		return err;
	}

}

/*Abs*/
else if (sampletype == 1)/*abs*/
{
	if (err = ABSCube(cube,1,(long) num_paths,1,2,1,(long) num_steps, &seed))
	{
		free_dmatrix(price,1,2,1,num_strikes);
		free_dvector(Forward,1,2);
		free_f3tensor(cube,1,(long) num_paths,1,2,1,(long) num_steps);
		free_dvector(blkbeta,1,num_strikes);
		smessage("error in abs cube");
		return err;
	}

}
 

		for (i = 1; i<=num_strikes; i++)
		{
			blkbeta[i] = srt_f_optblkschbeta(forw,strike[i],sigma,maturity,beta, 1.0,SRT_CALL, SRT_PREMIUM);
		}
		


for (cur_path = 1; cur_path<= num_paths;cur_path++)
{
	
		if (err=monte_carlo_sabr(cube[cur_path],forw,vovol,beta,rho,sigma,num_paths,maturity,num_steps,cur_path+1, Forward, modeltype, sampletype,&pathnogood))
		{
			free_dmatrix(price,1,2,1,num_strikes);
			free_f3tensor(cube,1,(long) num_paths,1,2,1,(long) num_steps);
			free_dvector(Forward,1,2);
			free_dvector(blkbeta,1,num_strikes);
			smessage("error in monte carlo");
			return err;
		}
	
	if (pathnogood == 0)
	{
		for (i = 1; i<=num_strikes; i++)
		{
			price[1][i] += DMAX(Forward[0]-strike[i],0);
			stdev[1][i] += DMAX(Forward[0]-strike[i],0)*DMAX(Forward[0]-strike[i],0);
			correl[i]   += DMAX(Forward[0]-strike[i],0)*DMAX(Forward[1]-strike[i],0);
			price[2][i] += DMAX(Forward[1]-strike[i],0);
			stdev[2][i] += DMAX(Forward[1]-strike[i],0)*DMAX(Forward[1]-strike[i],0);

		}
	}	
	else
	{
		actnumpaths--;
	}
				
}

		
	
	for (i=1;i<=num_strikes;i++)
	{
		

		correlation = (correl[i]-price[1][i]*price[2][i]/actnumpaths)/sqrt((stdev[1][i]-price[1][i]*price[1][i]/actnumpaths)*(stdev[2][i]-price[2][i]*price[2][i]/actnumpaths));
		
		price[1][i]=price[1][i]/actnumpaths-correlation*
			sqrt(stdev[1][i]-price[1][i]*price[1][i]/actnumpaths)/sqrt(stdev[2][i]-price[2][i]*price[2][i]/actnumpaths)*(price[2][i]/actnumpaths-blkbeta[i]);


/*		price[2][i]=sqrt(stdev[1][i]-price[1][i]*price[1][i]/actnumpaths)/sqrt(actnumpaths*(actnumpaths-1))*(1-correlation*correlation);*/

		if (err=srt_f_optimpvol(price[1][i],forw,strike[i],maturity,1,SRT_CALL,SRT_LOGNORMAL,&(impvol[i])))
		{
			impvol[i] = 0;
			smessage("error in implied vol");
			err = NULL;
		}

		
	}

/*free memory*/
		free_dmatrix(price,1,2,1,num_strikes);
		free_dmatrix(stdev,1,2,1,num_strikes);
		free_dvector(correl,1,num_strikes);
		free_f3tensor(cube,1,(long) num_paths,1,2,1,(long) num_steps);
		free_dvector(Forward,1,2);
		free_dvector(blkbeta,1,num_strikes);
		free_dvector(times_at_steps,1,(long) num_steps);
		return NULL;

}








