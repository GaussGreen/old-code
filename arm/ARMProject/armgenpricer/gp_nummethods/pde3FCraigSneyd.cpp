/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 *	\file pde3FCraigSneyd.cpp
 *	\author  K Belkheir
 *	\version 1.0
 *	\date July 2006
 *
 */

#include "gpnummethods/pde3Fnumericalschemes.h"
#include "gpinfra/pricingmodel.h"
#include "gpbase/vectormanip.h"
#include "gpbase/gpvector.h"
#include "gpbase/typedef.h"
#include "gpinfra/pricingstates.h"
#include "gpbase/gpmatrixlinalg.h"
#include "gpnummethods/argconvdefault.h"



CC_BEGIN_NAMESPACE( ARM )


// Karim Functions

////////////////////////////////////////////////////
///	Class  : NONE
///	Routine: Tridiag
///	Returns: void
///	Action : Solves for a vector u[1..n] the 
///   tridiagonal linear set given by equation (2.4.1).
///   a[1..n], b[1..n], c[1..n], and r[1..n] are input
///   vectors and are not modified.
////////////////////////////////////////////////////
void tridiag(ARM_GP_Vector &gam, ARM_GP_Vector &a, ARM_GP_Vector &b, ARM_GP_Vector &c, ARM_GP_Vector &r, ARM_GP_Vector &u, size_t n)
{
	unsigned long j;
	double bet;
	if (b[0] == 0.0) throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"b[1] == 0 in tridiagonalsolve" );;
	//If this happens then you should rewrite your equations as a set of order N - 1, w ith u2
	//trivially eliminated.
	u[0]=r[0]/(bet=b[0]);
	//cout << "u["<<0<<"]= "<<u[0]<<endl;
	for (j=2;j<=n;j++) 
	{ //Decomposition and forward substitution.
		gam[j-1]=c[j-2]/bet;
		bet=b[j-1]-a[j-1]*gam[j-1];
		if (bet == 0.0) throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"bet== 0 in tridiagonalsolve" );;
		u[j-1]=(r[j-1]-a[j-1]*u[j-2])/bet; 
	}
	for (j=(n-1);j>=1;j--)
	{
		u[j-1] -= gam[j]*u[j];//Backsubstitution.
	}
}

////////////////////////////////////////////////////
///	Class  : NONE
///	Routine: Tridiag2
///	Returns: void
///	Action : Solves for a vector u[1..n] the 
///   tridiagonal linear set given by equation (2.4.1).
///   a[1..n], b[1..n], c[1..n], and r[1..n] are input
///   vectors and are not modified.
////////////////////////////////////////////////////
void tridiag2(ARM_GP_Vector &gam, ARM_GP_Vector &a, ARM_GP_Vector &b, ARM_GP_Vector &c, ARM_GP_Vector &r, ARM_GP_Vector &u, size_t n)
{
	ARM_GP_Vector::iterator gamit	= gam.begin();
	ARM_GP_Vector::iterator ait	= a.begin();
	ARM_GP_Vector::iterator bit	= b.begin();
	ARM_GP_Vector::iterator cit	= c.begin();
	ARM_GP_Vector::iterator rit	= r.begin();
	ARM_GP_Vector::iterator uit	= u.begin();
	unsigned long j;
	double bet;
	if ((*bit) == 0.0) throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"b[1] == 0 in tridiagonalsolve" );;
	//If this happens then you should rewrite your equations as a set of order N - 1, w ith u2
	//trivially eliminated.
	(*uit)=(*rit)/(bet=*bit);
	gamit++;
	bit++;
	ait++;
	rit++;
	uit++;
	for (j=2;j<=n;j++) 
	{ //Decomposition and forward substitution.
		(*gamit)=(*cit)/bet;
		bet=(*bit)-(*ait)*(*gamit);
		if (bet == 0.0) throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"bet== 0 in tridiagonalsolve" );;
		(*uit)= ((*rit)-(*ait)*(*(uit-1)))/bet;

		gamit++;
		ait++;
		bit++;
		cit++;
		rit++;
		uit++;
	}
	for (j=(n-1);j>=1;j--)
	{
		(*(uit-1)) -= (*gamit)*(*uit);//Backsubstitution.
		uit--;
		gamit--;
	}
}
////////////////////////////////////////////////////
///	Class  : NONE
///	Routine: PermutXY
///	Returns: void
///	Action : Stock in the vector V the values of the
///   vector U read first in the y direction, then 
///   in the x direction and finally in the z one
////////////////////////////////////////////////////
void PermutXY(std::vector<double> &U, int nx, int ny, int nz, std::vector<double> &V )
{
	int i,j,k,I,I_per,K,K_per;
	unsigned int lengthU =U.size();
	
	//increments directionnels
	int e2=nx+1;
	int e3=(nx+1)*(ny+1);

	int e2_per=ny+1;
	int e3_per=(ny+1)*(nx+1);

	for(K=0;K<lengthU;K++)
	{
		I=K%e3;
		i=I%e2;
		j=I/e2;
		k=K/e3;
		I_per=e2_per*i+j;
		K_per=e3_per*k+I_per;
		V[K_per]=U[K];
	}
}

////////////////////////////////////////////////////
///	Class  : NONE
///	Routine: PermutXZ
///	Returns: void
///	Action : Stock in the vector V the values of the
///   vector U read first in the z direction, then 
///   in the y direction and finally in the x one
////////////////////////////////////////////////////
void PermutXZ(std::vector<double> &U, int nx, int ny, int nz, std::vector<double> &V)
{
	int i,j,k,I,I_per,K,K_per;
	unsigned int lengthU =U.size();
	
	//increments directionnels
	int e2=nx+1;
	int e3=(nx+1)*(ny+1);

	int e2_per=nz+1;
	int e3_per=(nz+1)*(ny+1);

	for(K=0;K<lengthU;K++)
	{
		I=K%e3;
		i=I%e2;
		j=I/e2;
		k=K/e3;
		I_per=e2_per*j+k;
		K_per=e3_per*i+I_per;
		V[K_per]=U[K];
	}
}

////////////////////////////////////////////////////
///	Class  : NONE
///	Routine: InvPermut
///	Returns: void
///	Action : Cancell the effect of the 2 functions
///   above
////////////////////////////////////////////////////
void InvPermut(ARM_GP_Vector &U, int nx, int ny, int nz, ARM_GP_Vector &V)
{
	int i,j,k,Izxy,Ixyz,Kzxy,Kxyz;
	unsigned int lengthU =U.size();
	
	//increments directionnels;
	int e2zxy=nz+1;
	int e3zxy=(nz+1)*(nx+1);

	int e2xyz=nx+1;
	int e3xyz=(nx+1)*(ny+1);

	for(Kzxy=0;Kzxy<lengthU;Kzxy++)
	{
		Izxy=Kzxy%e3zxy;
		i=Izxy%e2zxy;
		j=Izxy/e2zxy;
		k=Kzxy/e3zxy;
		Ixyz=e2xyz*k+j;
		Kxyz=e3xyz*i+Ixyz;
		V[Kxyz]=U[Kzxy];
	}
}

////////////////////////////////////////////////////
///	Class  : NONE
///	Routine: InvPermut2
///	Returns: void
///	Action : Cancell the effect of the 2 functions
///   above
////////////////////////////////////////////////////
void InvPermut2(ARM_GP_Vector &U, int nx, int ny, int nz, ARM_GP_Vector &V)
{
	int i,j,k,Izxy,Ixyz,Kzxy,Kxyz;
	unsigned int lengthU =U.size();
	
	//increments directionnels;
	int e2zxy=nz+1;
	int e3zxy=(nz+1)*(nx+1);

	int e2xyz=nx+1;
	int e3xyz=(nx+1)*(ny+1);

	ARM_GP_Vector::iterator Uit;
	ARM_GP_Vector::iterator Vit;

	for(Kzxy=0;Kzxy<lengthU;Kzxy++)
	{
		Izxy=Kzxy%e3zxy;
		i=Izxy%e2zxy;
		j=Izxy/e2zxy;
		k=Kzxy/e3zxy;
		Ixyz=e2xyz*k+j;
		Kxyz=e3xyz*i+Ixyz;
		Uit=U.begin()+Kzxy;
		Vit=V.begin()+Kxyz;
		(*Vit)=(*Uit);
	}
}

////////////////////////////////////////////////////
///	Class  : ARM_PDE3FCraigSneydNumericalScheme 
///	Routine: Reconstruct
///	Returns: void
///	Action : Build the vector newU which is the 
///		extension of the vector U_tilde to the 
///		borders according to the choice of the 
///		type of "LimitCondition":
///		the knowledge of oldU is necessary only 
///		for the Boundary condition "Belkheir"
////////////////////////////////////////////////////
void ARM_PDE3FCraigSneydNumericalScheme::Reconstruct(ARM_GP_Vector &U_tilde, ARM_GP_Vector &oldU, ARM_GP_Vector &newU )
{
	unsigned int lengthU=oldU.size();
	unsigned int lengthUtilde=U_tilde.size();

	int i,i_tilde,j,j_tilde,k,k_tilde,I,I_tilde,K,K_tilde;
	
	size_t Nx=GetNX()-1;
	size_t Ny=GetNY()-1;
	size_t Nz=GetNZ()-1;

	size_t Nx_tilde=Nx-2;
	size_t Ny_tilde=Ny-2;
	size_t Nz_tilde=Nz-2;
	
	//increments directionnels
	int e1=1;
	int e2=Nx+1;
	int e3=(Nx+1)*(Ny+1);
	int e1_tilde=1.0;
	int e2_tilde=Nx_tilde+1;
	int e3_tilde=(Nx_tilde+1)*(Ny_tilde+1);

	if ( lengthUtilde!=(Nx-1)*(Ny-1)*(Nz-1) ) 
	{
		throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"Utilde de mauvaise taille");
	}

	int ind=0;
	if ( (itsBC==VonNeumann)||(itsBC==Belkheir) )
	{
		ind=1;
	}
	int flag=0;
	switch(ind)
	{
		case 0:
		
			for(K=0;K<lengthU;K++)
			{
				I=K%e3;
				i=I%e2;
				j=I/e2;
				k=K/e3;
				if((i%Nx!=0)&&(j%Ny!=0)&&(k%Nz!=0))//points interrieurs
				{
					i_tilde=i-1;
					j_tilde=j-1;
					I_tilde=e2_tilde*j_tilde+i_tilde;
					k_tilde=k-1;
					K_tilde=e3_tilde*k_tilde+I_tilde;
					newU[K]=U_tilde[K_tilde];
				}
				else {newU[K]=0;}
			}
		break;
		
		case 1:
			if (itsBC==Belkheir){flag=1;}
			for(K=0;K<lengthU;K++)
			{
				I=K%e3;
				i=I%e2;
				j=I/e2;
				k=K/e3;
				//points interrieurs
				if((i%Nx!=0)&&(j%Ny!=0)&&(k%Nz!=0))
				{
					i_tilde=i-1;
					j_tilde=j-1;
					I_tilde=e2_tilde*j_tilde+i_tilde;
					k_tilde=k-1;
					K_tilde=e3_tilde*k_tilde+I_tilde;
					newU[K]=U_tilde[K_tilde];
				}
				//les 6 faces interrieures
				else if ((j==0)&&(i%Nx!=0)&&(k%Nz!=0))
				{
					i_tilde=i-1;
					j_tilde=0;
					I_tilde=e2_tilde*j_tilde+i_tilde;
					k_tilde=k-1;
					K_tilde=e3_tilde*k_tilde+I_tilde;
					newU[K]=U_tilde[K_tilde]		+	flag*(oldU[K]-oldU[K+e2]);
				}
				else if ((i==0)&&(j%Ny!=0)&&(k%Nz!=0))
				{
					i_tilde=0;
					j_tilde=j-1;
					I_tilde=e2_tilde*j_tilde+i_tilde;
					k_tilde=k-1;
					K_tilde=e3_tilde*k_tilde+I_tilde;
					newU[K]=U_tilde[K_tilde]		+	flag*(oldU[K]-oldU[K+e1]);
				}
				else if ((j==Ny)&&(i%Nx!=0)&&(k%Nz!=0))
				{
					i_tilde=i-1;
					j_tilde=j-2;
					I_tilde=e2_tilde*j_tilde+i_tilde;
					k_tilde=k-1;
					K_tilde=e3_tilde*k_tilde+I_tilde;
					newU[K]=U_tilde[K_tilde]		+	flag*(oldU[K]-oldU[K-e2]);
				}
				else if ((i==Nx)&&(j%Ny!=0)&&(k%Nz!=0))
				{
					i_tilde=i-2;
					j_tilde=j-1;
					I_tilde=e2_tilde*j_tilde+i_tilde;
					k_tilde=k-1;
					K_tilde=e3_tilde*k_tilde+I_tilde;
					newU[K]=U_tilde[K_tilde]		+	flag*(oldU[K]-oldU[K-e1]);
				}
				else if ((k==0)&&(j%Ny!=0)&&(i%Nx!=0))
				{
					i_tilde=i-1;
					j_tilde=j-1;
					I_tilde=e2_tilde*j_tilde+i_tilde;
					k_tilde=0;
					K_tilde=e3_tilde*k_tilde+I_tilde;
					newU[K]=U_tilde[K_tilde]		+	flag*(oldU[K]-oldU[K+e3]);
				}
				else if ((k==Nz)&&(j%Ny!=0)&&(i%Nx!=0))
				{
					i_tilde=i-1;
					j_tilde=j-1;
					I_tilde=e2_tilde*j_tilde+i_tilde;
					k_tilde=k-2;
					K_tilde=e3_tilde*k_tilde+I_tilde;
					newU[K]=U_tilde[K_tilde]		+	flag*(oldU[K]-oldU[K-e3]);
				}
				//les 12 arretes interrieures
				
				//4 arretes Basses
				else if ((j==0)&&(i%Nx!=0)&&(k==0))
				{
					i_tilde=i-1;
					j_tilde=0;
					I_tilde=e2_tilde*j_tilde+i_tilde;
					k_tilde=0;
					K_tilde=e3_tilde*k_tilde+I_tilde;
					newU[K]=U_tilde[K_tilde]		+	flag*(oldU[K]-oldU[K+e2+e3]);
				}
				else if ((i==0)&&(j%Ny!=0)&&(k==0))
				{
					i_tilde=0;
					j_tilde=j-1;
					I_tilde=e2_tilde*j_tilde+i_tilde;
					k_tilde=0;
					K_tilde=e3_tilde*k_tilde+I_tilde;
					newU[K]=U_tilde[K_tilde]		+	flag*(oldU[K]-oldU[K+e1+e3]);
				}
				else if ((j==Ny)&&(i%Nx!=0)&&(k==0))
				{
					i_tilde=i-1;
					j_tilde=j-2;
					I_tilde=e2_tilde*j_tilde+i_tilde;
					k_tilde=0;
					K_tilde=e3_tilde*k_tilde+I_tilde;
					newU[K]=U_tilde[K_tilde]		+	flag*(oldU[K]-oldU[K-e2+e3]);
				}
				else if ((i==Nx)&&(j%Ny!=0)&&(k==0))
				{
					i_tilde=i-2;
					j_tilde=j-1;
					I_tilde=e2_tilde*j_tilde+i_tilde;
					k_tilde=0;
					K_tilde=e3_tilde*k_tilde+I_tilde;
					newU[K]=U_tilde[K_tilde]		+	flag*(oldU[K]-oldU[K-e1+e3]);
				}

				//4 arretes intermediaires
				else if ((i==0)&&(j==0)&&(k%Nz!=0))
				{
					i_tilde=0;
					j_tilde=0;
					I_tilde=e2_tilde*j_tilde+i_tilde;
					k_tilde=k-1;
					K_tilde=e3_tilde*k_tilde+I_tilde;
					newU[K]=U_tilde[K_tilde]		+	flag*(oldU[K]-oldU[K+e1+e2]);
				}
				else if ((i==0)&&(j==Ny)&&(k%Nz!=0))
				{
					i_tilde=0;
					j_tilde=j-2;
					I_tilde=e2_tilde*j_tilde+i_tilde;
					k_tilde=k-1;
					K_tilde=e3_tilde*k_tilde+I_tilde;
					newU[K]=U_tilde[K_tilde]		+	flag*(oldU[K]-oldU[K+e1-e2]);
				}
				else if ((i==Nx)&&(j==Ny)&&(k%Nz!=0))
				{
					i_tilde=i-2;
					j_tilde=j-2;
					I_tilde=e2_tilde*j_tilde+i_tilde;
					k_tilde=k-1;
					K_tilde=e3_tilde*k_tilde+I_tilde;
					newU[K]=U_tilde[K_tilde]		+	flag*(oldU[K]-oldU[K-e1-e2]);
				}
				else if ((i==Nx)&&(j==0)&&(k%Nz!=0))
				{
					i_tilde=i-2;
					j_tilde=0;
					I_tilde=e2_tilde*j_tilde+i_tilde;
					k_tilde=k-1;
					K_tilde=e3_tilde*k_tilde+I_tilde;
					newU[K]=U_tilde[K_tilde]		+	flag*(oldU[K]-oldU[K-e1+e2]);
				}
				//4 arretes Hautes
				else if ((j==0)&&(i%Nx!=0)&&(k==Nz))
				{
					i_tilde=i-1;
					j_tilde=0;
					I_tilde=e2_tilde*j_tilde+i_tilde;
					k_tilde=k-2;
					K_tilde=e3_tilde*k_tilde+I_tilde;
					newU[K]=U_tilde[K_tilde]		+	flag*(oldU[K]-oldU[K+e2-e3]);
				}
				else if ((i==0)&&(j%Ny!=0)&&(k==Nz))
				{
					i_tilde=0;
					j_tilde=j-1;
					I_tilde=e2_tilde*j_tilde+i_tilde;
					k_tilde=k-2;
					K_tilde=e3_tilde*k_tilde+I_tilde;
					newU[K]=U_tilde[K_tilde]		+	flag*(oldU[K]-oldU[K+e1-e3]);
				}
				else if ((j==Ny)&&(i%Nx!=0)&&(k==Nz))
				{
					i_tilde=i-1;
					j_tilde=j-2;
					I_tilde=e2_tilde*j_tilde+i_tilde;
					k_tilde=k-2;
					K_tilde=e3_tilde*k_tilde+I_tilde;
					newU[K]=U_tilde[K_tilde]		+	flag*(oldU[K]-oldU[K-e2-e3]);
				}
				else if ((i==Nx)&&(j%Ny!=0)&&(k==Nz))
				{
					i_tilde=i-2;
					j_tilde=j-1;
					I_tilde=e2_tilde*j_tilde+i_tilde;
					k_tilde=k-2;
					K_tilde=e3_tilde*k_tilde+I_tilde;
					newU[K]=U_tilde[K_tilde]		+	flag*(oldU[K]-oldU[K-e1-e3]);
				}
			}	
			//les Huits sommets
			newU[0]=U_tilde[0]	+	flag*(oldU[0]-oldU[0+e1+e2+e3]);
			newU[Nx]=U_tilde[Nx_tilde]		+	flag*(oldU[Nx]-oldU[Nx -e1+e2+e3]);
			newU[Ny*e2]=U_tilde[Ny_tilde*e2_tilde]		+	flag*(oldU[Ny*e2]-oldU[Ny*e2 +e1-e2+e3]);
			newU[Nx+Ny*e2]=U_tilde[Nx_tilde+Ny_tilde*e2_tilde]		+	flag*(oldU[Nx+Ny*e2]-oldU[Nx+Ny*e2 -e1-e2+e3]);
					
			newU[Nz*e3]=U_tilde[Nz_tilde*e3_tilde]	+	flag*(oldU[Nz*e3]-oldU[Nz*e3 +e1+e2-e3]);
			newU[Nx+Nz*e3]=U_tilde[Nx_tilde+Nz_tilde*e3_tilde]		+	flag*(oldU[Nx+Nz*e3]-oldU[Nx+Nz*e3 -e1+e2-e3]);
			newU[Ny*e2+Nz*e3]=U_tilde[Ny_tilde*e2_tilde+Nz_tilde*e3_tilde]		+	flag*(oldU[Ny*e2+Nz*e3]-oldU[Ny*e2+Nz*e3 +e1-e2-e3]);
			newU[Nx+Ny*e2+Nz*e3]=U_tilde[Nx_tilde+Ny_tilde*e2_tilde+Nz_tilde*e3_tilde]		+	flag*(oldU[Nx+Ny*e2+Nz*e3]-oldU[Nx+Ny*e2+Nz*e3 -e1-e2-e3]);
		break;
		
		default:
			ARM_THROW( ERR_INVALID_ARGUMENT, "Unknown Boundary Condition");
		}
}


////////////////////////////////////////////////////
///	Class  : ARM_PDE3FCraigSneydNumericalScheme 
///	Routine: Reconstruct2
///	Returns: void
///	Action : Fill the boudary of the PDE based on the
/// conditions (Dirichlet, Von Neumann)
////////////////////////////////////////////////////
void ARM_PDE3FCraigSneydNumericalScheme::Reconstruct2(ARM_GP_Vector &U)
{	
	for (ARM_GP_T_Vector<size_t>::iterator it1 = itsReconstructFrom.begin(), it2 = itsReconstructTo.begin();
	it1!=itsReconstructFrom.end(); ++it1,++it2)
		(*(U.begin()+(*it2))) = (*(U.begin()+(*it1)));
}



// GP Functions

////////////////////////////////////////////////////
///	Class  : ARM_PDE3FCraigSneydNumericalScheme
///	Routine: constructor
///	Returns:
///	Action : 
////////////////////////////////////////////////////

ARM_PDE3FCraigSneydNumericalScheme::ARM_PDE3FCraigSneydNumericalScheme(
		size_t NX, 
		size_t NY, 
		size_t NZ,
		double theta1, 
		double theta2, 
		double theta3, 
		BoundConditionType BC, 
		double lambda,
		GridType gridType,
		const ARM_GP_Matrix& gridData,
		const ARM_GP_Matrix& schedulerData
		)
:
ARM_PDE3FNumericalScheme(NX,NY,NZ,gridType,gridData,schedulerData),
itsTheta1(theta1),
itsTheta2(theta2),
itsTheta3(theta3),
itsBC(BC),
itsLambda(lambda),

itsPrevO(NULL),
itsPrevPx(NULL),
itsPrevPy(NULL),
itsPrevPz(NULL),
itsPrevQxx(NULL),
itsPrevQyy(NULL),
itsPrevQzz(NULL),
itsPrevQxy(NULL),
itsPrevQyz(NULL),
itsPrevQzx(NULL),

itsO(NULL),
itsPx(NULL),
itsPy(NULL),
itsPz(NULL),
itsQxx(NULL),
itsQyy(NULL),
itsQzz(NULL),
itsQxy(NULL),
itsQyz(NULL),
itsQzx(NULL),

itsA1(NULL),
itsB1(NULL),
itsC1(NULL),

itsA2(NULL),
itsB2(NULL),
itsC2(NULL),

itsA3(NULL),
itsB3(NULL),
itsC3(NULL),

itsD(NULL),
itsGam(NULL),
itsU_new(NULL),
itsUtilde(NULL),
itsUtilde_temp(NULL),

itsIsBound(NULL),
itsReconstructFrom(NULL),
itsReconstructTo(NULL)
{
}

////////////////////////////////////////////////////
///	Class  : ARM_PDE3FCraigSneydNumericalScheme
///	Routine: copy constructor
///	Returns:
///	Action : 
////////////////////////////////////////////////////

ARM_PDE3FCraigSneydNumericalScheme::ARM_PDE3FCraigSneydNumericalScheme(const ARM_PDE3FCraigSneydNumericalScheme& rhs)
:
ARM_PDE3FNumericalScheme(rhs),
itsTheta1(rhs.itsTheta1),
itsTheta2(rhs.itsTheta2),
itsTheta3(rhs.itsTheta3),
itsBC(rhs.itsBC),
itsLambda(rhs.itsLambda),

itsPrevO(rhs.itsPrevO),
itsPrevPx(rhs.itsPrevPx),
itsPrevPy(rhs.itsPrevPy),
itsPrevPz(rhs.itsPrevPz),
itsPrevQxx(rhs.itsPrevQxx),
itsPrevQyy(rhs.itsPrevQyy),
itsPrevQzz(rhs.itsPrevQzz),
itsPrevQxy(rhs.itsPrevQxy),
itsPrevQyz(rhs.itsPrevQyz),
itsPrevQzx(rhs.itsPrevQzx),

itsO(rhs.itsO),
itsPx(rhs.itsPx),
itsPy(rhs.itsPy),
itsPz(rhs.itsPz),
itsQxx(rhs.itsQxx),
itsQyy(rhs.itsQyy),
itsQzz(rhs.itsQzz),
itsQxy(rhs.itsQxy),
itsQyz(rhs.itsQyz),
itsQzx(rhs.itsQzx),

itsA1(rhs.itsA1),
itsB1(rhs.itsB1),
itsC1(rhs.itsC1),

itsA2(rhs.itsA2),
itsB2(rhs.itsB2),
itsC2(rhs.itsC2),

itsA3(rhs.itsA3),
itsB3(rhs.itsB3),
itsC3(rhs.itsC3),

itsD(rhs.itsD),
itsGam(rhs.itsGam),
itsU_new(rhs.itsU_new),
itsUtilde(rhs.itsUtilde),
itsUtilde_temp(rhs.itsUtilde_temp),

itsIsBound(rhs.itsIsBound),
itsReconstructFrom(rhs.itsReconstructFrom),
itsReconstructTo(rhs.itsReconstructTo)
{
}

////////////////////////////////////////////////////
///	Class  : ARM_PDE3FCraigSneydNumericalScheme
///	Routine: destructor
///	Returns:
///	Action : 
////////////////////////////////////////////////////

ARM_PDE3FCraigSneydNumericalScheme::~ARM_PDE3FCraigSneydNumericalScheme()
{
}

////////////////////////////////////////////////////
///	Class  : ARM_PDE3FNumericalScheme
///	Routine: Init
///	Returns: ARM_PricingStatesPtr
///	Action : Initialiation of the NumericalScheme
///  builds nummethodstates, discretization schemes
///  and the coefficients of the PDE 
////////////////////////////////////////////////////

ARM_PricingStatesPtr ARM_PDE3FCraigSneydNumericalScheme::Init( const ARM_PricingModel& model )
{
	ARM_PricingStatesPtr tempPricingStates;
	tempPricingStates=Init2(model);//before it was Init1
	return tempPricingStates;
}

////////////////////////////////////////////////////
///	Class  : ARM_PDE3FNumericalScheme
///	Routine: Init1
///	Returns: ARM_PricingStatesPtr
///	Action : Initialiation of the NumericalScheme
///  builds nummethodstates, discretization schemes
///  and the coefficients of the PDE 
///			OLD VERSION WHICH IS SLOW
////////////////////////////////////////////////////


ARM_PricingStatesPtr ARM_PDE3FCraigSneydNumericalScheme::Init1( const ARM_PricingModel& model )
{
	//execution of the Init of ARM_PDE3FNumericalScheme
	ARM_PricingStatesPtr tempPricingStates = ARM_PDE3FNumericalScheme::Init( model );

	//final date
	ARM_GP_VectorPtr timeSteps( new ARM_GP_Vector(*model.GetNumMethod()->GetTimeSteps()));

	//Initialisation of the coeffs of the PDE
	size_t lengthU=GetXGrid()->size();
	itsQxx		= ARM_GP_VectorPtr( new ARM_GP_Vector(lengthU) );
	itsQyy		= ARM_GP_VectorPtr( new ARM_GP_Vector(lengthU) );
	itsQzz		= ARM_GP_VectorPtr( new ARM_GP_Vector(lengthU) );
	itsQxy		= ARM_GP_VectorPtr( new ARM_GP_Vector(lengthU) );
	itsQyz		= ARM_GP_VectorPtr( new ARM_GP_Vector(lengthU) );
	itsQzx		= ARM_GP_VectorPtr( new ARM_GP_Vector(lengthU) );
	itsPx		= ARM_GP_VectorPtr( new ARM_GP_Vector(lengthU) );
	itsPy		= ARM_GP_VectorPtr( new ARM_GP_Vector(lengthU) );
	itsPz		= ARM_GP_VectorPtr( new ARM_GP_Vector(lengthU) );
	itsO		= ARM_GP_VectorPtr( new ARM_GP_Vector(lengthU) );
	itsPrevQxx  = ARM_GP_VectorPtr( new ARM_GP_Vector(lengthU) );
	itsPrevQyy  = ARM_GP_VectorPtr( new ARM_GP_Vector(lengthU) );
	itsPrevQzz  = ARM_GP_VectorPtr( new ARM_GP_Vector(lengthU) );
	itsPrevQxy  = ARM_GP_VectorPtr( new ARM_GP_Vector(lengthU) );
	itsPrevQyz  = ARM_GP_VectorPtr( new ARM_GP_Vector(lengthU) );
	itsPrevQzx  = ARM_GP_VectorPtr( new ARM_GP_Vector(lengthU) );
	itsPrevPx   = ARM_GP_VectorPtr( new ARM_GP_Vector(lengthU) );
	itsPrevPy	= ARM_GP_VectorPtr( new ARM_GP_Vector(lengthU) );
	itsPrevPz	= ARM_GP_VectorPtr( new ARM_GP_Vector(lengthU) );
	itsPrevO	= ARM_GP_VectorPtr( new ARM_GP_Vector(lengthU) );
	
	//Setting of the coefficients
	model.UpdatePDE3DCoeffs(
	timeSteps->size()-2,
	tempPricingStates,
	itsQxx,
	itsQyy,
	itsQzz,
	itsQxy,
	itsQyz,
	itsQzx,
	itsPx,
	itsPy,
	itsPz,
	itsO,
	false,
	true);

	//initialisation of the others vectors of the class
	//tailles
	size_t Nx=GetNX()-1;
	size_t Ny=GetNY()-1;
	size_t Nz=GetNZ()-1;	
	size_t lengthUtilde=(Nx-1)*(Ny-1)*(Nz-1);
	
	itsA1 = ARM_GP_VectorPtr( new ARM_GP_Vector( lengthUtilde ) );
	itsB1 = ARM_GP_VectorPtr( new ARM_GP_Vector( lengthUtilde ) );
	itsC1 = ARM_GP_VectorPtr( new ARM_GP_Vector( lengthUtilde ) );

	itsA2 = ARM_GP_VectorPtr( new ARM_GP_Vector( lengthUtilde ) );
	itsB2 = ARM_GP_VectorPtr( new ARM_GP_Vector( lengthUtilde ) );
	itsC2 = ARM_GP_VectorPtr( new ARM_GP_Vector( lengthUtilde ) );
	
	itsA3 = ARM_GP_VectorPtr( new ARM_GP_Vector( lengthUtilde ) );
	itsB3 = ARM_GP_VectorPtr( new ARM_GP_Vector( lengthUtilde ) );
	itsC3 = ARM_GP_VectorPtr( new ARM_GP_Vector( lengthUtilde ) );
	
	itsD = ARM_GP_VectorPtr( new ARM_GP_Vector( lengthUtilde ) );
	itsGam = ARM_GP_VectorPtr( new ARM_GP_Vector( lengthUtilde ) );
	itsUtilde = ARM_GP_VectorPtr( new ARM_GP_Vector( lengthUtilde ) );
	itsUtilde_temp = ARM_GP_VectorPtr( new ARM_GP_Vector( lengthUtilde ) );
	itsU_new = ARM_GP_VectorPtr( new ARM_GP_Vector( GetXGrid()->size() ) );
	
	return tempPricingStates; 
}


////////////////////////////////////////////////////
///	Class  : ARM_PDE3FNumericalScheme
///	Routine: Init2
///	Returns: ARM_PricingStatesPtr
///	Action : Initialiation of the NumericalScheme
///  builds nummethodstates, discretization schemes
///  and the coefficients of the PDE 
////////////////////////////////////////////////////

ARM_PricingStatesPtr ARM_PDE3FCraigSneydNumericalScheme::Init2( const ARM_PricingModel& model )
{
	//execution of the Init of ARM_PDE3FNumericalScheme
	ARM_PricingStatesPtr tempPricingStates = ARM_PDE3FNumericalScheme::Init( model, itsLambda );

	//final date
	ARM_GP_VectorPtr timeSteps( new ARM_GP_Vector(*model.GetNumMethod()->GetTimeSteps() ));

	//Initialisation of the coeffs of the PDE
	size_t lengthU=GetXGrid()->size();
	itsQxx		= ARM_GP_VectorPtr( new ARM_GP_Vector(lengthU) );
	itsQyy		= ARM_GP_VectorPtr( new ARM_GP_Vector(lengthU) );
	itsQzz		= ARM_GP_VectorPtr( new ARM_GP_Vector(lengthU) );
	itsQxy		= ARM_GP_VectorPtr( new ARM_GP_Vector(lengthU) );
	itsQyz		= ARM_GP_VectorPtr( new ARM_GP_Vector(lengthU) );
	itsQzx		= ARM_GP_VectorPtr( new ARM_GP_Vector(lengthU) );
	itsPx		= ARM_GP_VectorPtr( new ARM_GP_Vector(lengthU) );
	itsPy		= ARM_GP_VectorPtr( new ARM_GP_Vector(lengthU) );
	itsPz		= ARM_GP_VectorPtr( new ARM_GP_Vector(lengthU) );
	itsO		= ARM_GP_VectorPtr( new ARM_GP_Vector(lengthU) );
	itsPrevQxx  = ARM_GP_VectorPtr( new ARM_GP_Vector(lengthU) );
	itsPrevQyy  = ARM_GP_VectorPtr( new ARM_GP_Vector(lengthU) );
	itsPrevQzz  = ARM_GP_VectorPtr( new ARM_GP_Vector(lengthU) );
	itsPrevQxy  = ARM_GP_VectorPtr( new ARM_GP_Vector(lengthU) );
	itsPrevQyz  = ARM_GP_VectorPtr( new ARM_GP_Vector(lengthU) );
	itsPrevQzx  = ARM_GP_VectorPtr( new ARM_GP_Vector(lengthU) );
	itsPrevPx   = ARM_GP_VectorPtr( new ARM_GP_Vector(lengthU) );
	itsPrevPy	= ARM_GP_VectorPtr( new ARM_GP_Vector(lengthU) );
	itsPrevPz	= ARM_GP_VectorPtr( new ARM_GP_Vector(lengthU) );
	itsPrevO	= ARM_GP_VectorPtr( new ARM_GP_Vector(lengthU) );
	
	//setting of the model states
	model.PdeStatesToModelStates(tempPricingStates, timeSteps->size()-1, itsLambda);

	//Setting of the coefficients
	model.UpdatePDE3DCoeffs(
	timeSteps->size()-1,
	tempPricingStates,
	itsQxx,
	itsQyy,
	itsQzz,
	itsQxy,
	itsQyz,
	itsQzx,
	itsPx,
	itsPy,
	itsPz,
	itsO,
	itsLambda,
	true);
	

	//initialisation of the others vectors of the class
	//tailles
	
	itsA1 = ARM_GP_VectorPtr( new ARM_GP_Vector( lengthU, 0.0 ) );
	itsB1 = ARM_GP_VectorPtr( new ARM_GP_Vector( lengthU, 1.0 ) );
	itsC1 = ARM_GP_VectorPtr( new ARM_GP_Vector( lengthU, 0.0 ) );

	itsA2 = ARM_GP_VectorPtr( new ARM_GP_Vector( lengthU, 0.0 ) );
	itsB2 = ARM_GP_VectorPtr( new ARM_GP_Vector( lengthU, 1.0 ) );
	itsC2 = ARM_GP_VectorPtr( new ARM_GP_Vector( lengthU, 0.0 ) );
	
	itsA3 = ARM_GP_VectorPtr( new ARM_GP_Vector( lengthU, 0.0 ) );
	itsB3 = ARM_GP_VectorPtr( new ARM_GP_Vector( lengthU, 1.0 ) );
	itsC3 = ARM_GP_VectorPtr( new ARM_GP_Vector( lengthU, 0.0 ) );
	
	itsD = ARM_GP_VectorPtr( new ARM_GP_Vector( lengthU ) );
	itsGam = ARM_GP_VectorPtr( new ARM_GP_Vector( lengthU ) );
	itsU_temp = ARM_GP_VectorPtr( new ARM_GP_Vector( lengthU, 0.0 ) );
	itsU_new = ARM_GP_VectorPtr( new ARM_GP_Vector( lengthU ) );

	itsRxx		= ARM_GP_VectorPtr( new ARM_GP_Vector(lengthU) );
	itsRyy		= ARM_GP_VectorPtr( new ARM_GP_Vector(lengthU) );
	itsRzz		= ARM_GP_VectorPtr( new ARM_GP_Vector(lengthU) );
	itsRxy		= ARM_GP_VectorPtr( new ARM_GP_Vector(lengthU) );
	itsRyz		= ARM_GP_VectorPtr( new ARM_GP_Vector(lengthU) );
	itsRzx		= ARM_GP_VectorPtr( new ARM_GP_Vector(lengthU) );
	itsSx		= ARM_GP_VectorPtr( new ARM_GP_Vector(lengthU) );
	itsSy		= ARM_GP_VectorPtr( new ARM_GP_Vector(lengthU) );
	itsSz		= ARM_GP_VectorPtr( new ARM_GP_Vector(lengthU) );

	itsIsBound			=	ARM_GP_T_Vector<bool>(lengthU) ;
	InitReconstructAndGeometry();

	return tempPricingStates; 
}


////////////////////////////////////////////////////
///	Class  : ARM_PDE3FNumericalScheme
///	Routine: InitGeometry
///	Returns: void
///	Action : Initialize the geometry coefficient
/// of the scheme
////////////////////////////////////////////////////

void ARM_PDE3FCraigSneydNumericalScheme::InitReconstructAndGeometry()
{
	size_t Nx=GetNX()-1;
	size_t Ny=GetNY()-1;
	size_t Nz=GetNZ()-1;
	
	size_t lengthU=(*itsQxx).size();
	
	//incréments directionnels
	size_t e1=1;
	size_t e2=Nx+1;
	size_t e3=(Nx+1)*(Ny+1);

	size_t K,I;
	size_t i,j,k;

	ARM_GP_Vector::iterator deltaXit;
	ARM_GP_Vector::iterator deltaYit;
	ARM_GP_Vector::iterator deltaZit;
	ARM_GP_Vector::iterator Rxxit	= itsRxx->begin();
	ARM_GP_Vector::iterator Ryyit	= itsRyy->begin();
	ARM_GP_Vector::iterator Rzzit	= itsRzz->begin();
	ARM_GP_Vector::iterator Rxyit	= itsRxy->begin();
	ARM_GP_Vector::iterator Ryzit	= itsRyz->begin();
	ARM_GP_Vector::iterator Rzxit	= itsRzx->begin();
	ARM_GP_Vector::iterator Sxit	= itsSx->begin();
	ARM_GP_Vector::iterator Syit	= itsSy->begin();
	ARM_GP_Vector::iterator Szit	= itsSz->begin();
	ARM_GP_T_Vector<bool>::iterator	IsBoundit = itsIsBound.begin();

	for(K=0;K<lengthU;K++)
	{
		I=K%e3;
		i=I%e2;
		j=I/e2;
		k=K/e3;
				
		//indices qui vont intervenir dans les differences finies
			
		if( (i%Nx!=0)&&(j%Ny!=0)&&(k%Nz!=0) )//interrieur du domaine
		{
			//constantes du schéma

			deltaXit = GetDeltaX()->begin() + i;
			deltaYit = GetDeltaY()->begin() + j;
			deltaZit = GetDeltaZ()->begin() + k;
		
			(*Rxxit)=2/( pow( *(deltaXit-1),2 ) + pow( *(deltaXit),2 ) );
			(*Ryyit)=2/( pow( *(deltaYit-1),2 ) + pow( *(deltaYit),2 ) );
			(*Rzzit)=2/( pow( *(deltaZit-1),2 ) + pow( *(deltaZit),2 ) );
		
			(*Rxyit)=4/( ( *(deltaXit-1) + *(deltaXit) ) * ( *(deltaYit-1) + *(deltaYit) ) );
			(*Ryzit)=4/( ( *(deltaYit-1) + *(deltaYit) ) * ( *(deltaZit-1) + *(deltaZit) ) );
			(*Rzxit)=4/( ( *(deltaZit-1) + *(deltaZit) ) * ( *(deltaXit-1) + *(deltaXit) ) );
			
			(*Sxit)= 2/( *(deltaXit-1) + *(deltaXit) );
			(*Syit)= 2/( *(deltaYit-1) + *(deltaYit) );
			(*Szit)= 2/( *(deltaZit-1) + *(deltaZit) );

			(*IsBoundit) = false;
		}
		else
		{
			(*IsBoundit) = true;

			//les 6 faces interrieures
			if ((j==0)&&(i%Nx!=0)&&(k%Nz!=0))
			{
				itsReconstructFrom.push_back(K+e2);
				itsReconstructTo.push_back(K);

			}
			else if ((i==0)&&(j%Ny!=0)&&(k%Nz!=0))
			{
				itsReconstructFrom.push_back(K+e1);
				itsReconstructTo.push_back(K);
			}
			else if ((j==Ny)&&(i%Nx!=0)&&(k%Nz!=0))
			{
				itsReconstructFrom.push_back(K-e2);
				itsReconstructTo.push_back(K);
			}
			else if ((i==Nx)&&(j%Ny!=0)&&(k%Nz!=0))
			{
				itsReconstructFrom.push_back(K-e1);
				itsReconstructTo.push_back(K);
			}
			else if ((k==0)&&(j%Ny!=0)&&(i%Nx!=0))
			{
				itsReconstructFrom.push_back(K+e3);
				itsReconstructTo.push_back(K);
			}
			else if ((k==Nz)&&(j%Ny!=0)&&(i%Nx!=0))
			{
				itsReconstructFrom.push_back(K-e3);
				itsReconstructTo.push_back(K);
			}
			//les 12 arretes interrieures
			
			//4 arretes Basses
			else if ((j==0)&&(i%Nx!=0)&&(k==0))
			{
				itsReconstructFrom.push_back(K+e2+e3);
				itsReconstructTo.push_back(K);
			}
			else if ((i==0)&&(j%Ny!=0)&&(k==0))
			{
				itsReconstructFrom.push_back(K+e1+e3);
				itsReconstructTo.push_back(K);
			}
			else if ((j==Ny)&&(i%Nx!=0)&&(k==0))
			{
				itsReconstructFrom.push_back(K-e2+e3);
				itsReconstructTo.push_back(K);
			}
			else if ((i==Nx)&&(j%Ny!=0)&&(k==0))
			{
				itsReconstructFrom.push_back(K-e1+e3);
				itsReconstructTo.push_back(K);
			}

			//4 arretes intermediaires
			else if ((i==0)&&(j==0)&&(k%Nz!=0))
			{
				itsReconstructFrom.push_back(K+e1+e2);
				itsReconstructTo.push_back(K);
			}
			else if ((i==0)&&(j==Ny)&&(k%Nz!=0))
			{
				itsReconstructFrom.push_back(K+e1-e2);
				itsReconstructTo.push_back(K);
			}
			else if ((i==Nx)&&(j==Ny)&&(k%Nz!=0))
			{
				itsReconstructFrom.push_back(K-e1-e2);
				itsReconstructTo.push_back(K);
			}
			else if ((i==Nx)&&(j==0)&&(k%Nz!=0))
			{
				itsReconstructFrom.push_back(K-e1+e2);
				itsReconstructTo.push_back(K);
			}
			//4 arretes Hautes
			else if ((j==0)&&(i%Nx!=0)&&(k==Nz))
			{
				itsReconstructFrom.push_back(K+e2-e3);
				itsReconstructTo.push_back(K);
			}
			else if ((i==0)&&(j%Ny!=0)&&(k==Nz))
			{
				itsReconstructFrom.push_back(K+e1-e3);
				itsReconstructTo.push_back(K);
			}
			else if ((j==Ny)&&(i%Nx!=0)&&(k==Nz))
			{
				itsReconstructFrom.push_back(K-e2-e3);
				itsReconstructTo.push_back(K);
			}
			else if ((i==Nx)&&(j%Ny!=0)&&(k==Nz))
			{
				itsReconstructFrom.push_back(K-e1-e3);
				itsReconstructTo.push_back(K);
			}
			//8 sommets
			else if (K==0)
			{
				itsReconstructTo.push_back(0);
				itsReconstructFrom.push_back(e1+e2+e3);
			}
			else if (K==Nx)
			{
				itsReconstructTo.push_back(Nx);
				itsReconstructFrom.push_back(Nx-e1+e2+e3);
			}
			else if (K==Ny*e2)
			{
				itsReconstructTo.push_back(Ny*e2);
				itsReconstructFrom.push_back(Ny*e2+e1-e2+e3);
			}
			else if (K==Nx+Ny*e2)
			{
				itsReconstructTo.push_back(Nx+Ny*e2);
				itsReconstructFrom.push_back(Nx+Ny*e2-e1-e2+e3);
			}
			else if (K==Nz*e3)
			{
				itsReconstructTo.push_back(Nz*e3);
				itsReconstructFrom.push_back(Nz*e3+e1+e2-e3);
			}
			else if (K==Nx+Nz*e3)
			{
				itsReconstructTo.push_back(Nx+Nz*e3);
				itsReconstructFrom.push_back(Nx+Nz*e3-e1+e2-e3);
			}
			else if (K==Ny*e2+Nz*e3)
			{
				itsReconstructTo.push_back(Ny*e2+Nz*e3);
				itsReconstructFrom.push_back(Ny*e2+Nz*e3+e1-e2-e3);
			}
			else 
			{
				itsReconstructTo.push_back(Nx+Ny*e2+Nz*e3);
				itsReconstructFrom.push_back(Nx+Ny*e2+Nz*e3-e1-e2-e3);
			}
		}

		Rxxit++;
		Ryyit++;
		Rzzit++;
		Rxyit++;
		Ryzit++;
		Rzxit++;
		Sxit++;
		Syit++;
		Szit++;
		IsBoundit++;

	}
}

////////////////////////////////////////////////////
///	Class  : ARM_PDE3FCraigSneydNumericalScheme 
///	Routine: Induct
///	Returns: void
///	Action : induct nummethodstates/prices from one date to
///  another
////////////////////////////////////////////////////

void ARM_PDE3FCraigSneydNumericalScheme::Induct( const ARM_PricingModel& model, ARM_PricingStatesPtr& states, int toTimeIdx)
{
	Induct2(model, states, toTimeIdx);//before it was Induct1
}

////////////////////////////////////////////////////
///	Class  : ARM_PDE3FCraigSneydNumericalScheme 
///	Routine: Induct1
///	Returns: void
///	Action : induct nummethodstates/prices from one date to
///  another
			///OLD VERSION which was SLOW
////////////////////////////////////////////////////

/*
void ARM_PDE3FCraigSneydNumericalScheme::Induct1( const ARM_PricingModel& model, ARM_PricingStatesPtr& states, int toTimeIdx)
{
	size_t i;
	size_t PayoffsSize = states->GetPayoffsSize();
	size_t PayoffStatesSize = states->GetPayoffStatesSize();

	// Il faut garder les anciennes valeurs !!!!!!
	(*itsPrevQxx) = (*itsQxx);
	(*itsPrevQyy) = (*itsQyy);
	(*itsPrevQzz) = (*itsQzz);
	(*itsPrevQxy) = (*itsQxy);
	(*itsPrevQyz) = (*itsQyz);
	(*itsPrevQzx) = (*itsQzx);
	(*itsPrevPx) = (*itsPx);
	(*itsPrevPy) = (*itsPy);
	(*itsPrevPz) = (*itsPz);
	(*itsPrevO) = (*itsO);

	//Calculation of the coefficients of the PDE at the next date
	double time = (*GetTimeSteps())[toTimeIdx];
	GetModel()->UpdatePDE3DCoeffs(
		toTimeIdx,
		states,
		itsQxx,
		itsQyy,
		itsQzz,
		itsQxy,
		itsQyz,
		itsQzx,
		itsPx,
		itsPy,
		itsPz,
		itsO,
		false,
		false);


	/// PayoffStatesSize
	//flag to avoid the calculation of the intermediate payoff
	bool isOtherPayoffs;
	//Initialisation of all the coefficients of the tridiag systemes
	InitABC(toTimeIdx);
	//loop
	for( i=0 ; i<PayoffStatesSize ; ++i)
	{
		ARM_PayoffStates& payoffStates = states->GetPayoffStates(i);
		isOtherPayoffs = states->GetPayoffStates(i).GetOtherPayoffsFlag();

		if(isOtherPayoffs)
		{	
			ARM_GP_Matrix& intermediatePayoffs = const_cast<ARM_GP_Matrix&> (payoffStates.GetIntermediatePayoffs() );
			UpdateIntermediatePayoffs( intermediatePayoffs, toTimeIdx );
		}
		ARM_GP_Vector vec = payoffStates.GetPayoffs();
		UpdatePayoffs( vec, toTimeIdx  );
		payoffStates.SetPayoffs( vec );
	}

	/// Payoffs
	for( i=0 ; i<PayoffsSize ; ++i)
	{
		ARM_VectorPtr vec = states->GetPayoffVec(i);
		UpdatePayoffs( *vec, toTimeIdx  );

		for( int j=0; j<vec->size() ; ++j)
			states->SetPayoff( j, i, vec->Elt(j) );

	}
}
*/

////////////////////////////////////////////////////
///	Class  : ARM_PDE3FCraigSneydNumericalScheme 
///	Routine: Induct2
///	Returns: void
///	Action : induct nummethodstates/prices from one date to
///  another
////////////////////////////////////////////////////

void ARM_PDE3FCraigSneydNumericalScheme::Induct2( const ARM_PricingModel& model, ARM_PricingStatesPtr& states, int toTimeIdx)
{
	size_t i;
	size_t PayoffsSize = states->GetPayoffsSize();
	size_t PayoffStatesSize = states->GetPayoffStatesSize();

	// Il faut garder les anciennes valeurs !!!!!!
	(*itsPrevQxx) = (*itsQxx);
	(*itsPrevQyy) = (*itsQyy);
	(*itsPrevQzz) = (*itsQzz);
	(*itsPrevQxy) = (*itsQxy);
	(*itsPrevQyz) = (*itsQyz);
	(*itsPrevQzx) = (*itsQzx);
	(*itsPrevPx) = (*itsPx);
	(*itsPrevPy) = (*itsPy);
	(*itsPrevPz) = (*itsPz);
	(*itsPrevO) = (*itsO);

	//we update the modelstate (must be done in principle after the induct)
	model.PdeStatesToModelStates(states, toTimeIdx, itsLambda);

	//Calculation of the coefficients of the PDE at the next date
	double time = (*GetTimeSteps())[toTimeIdx];
	GetModel()->UpdatePDE3DCoeffs(
		toTimeIdx,
		states,
		itsQxx,
		itsQyy,
		itsQzz,
		itsQxy,
		itsQyz,
		itsQzx,
		itsPx,
		itsPy,
		itsPz,
		itsO,
		itsLambda,
		false);

	double deltaT=( (*GetTimeSteps())[toTimeIdx+1]-(*GetTimeSteps())[toTimeIdx] )/K_YEAR_LEN;
	// Backward To Forward induction;
	deltaT=-deltaT;

	// We first introduce the deltaT (not so ugly ....)

	(*itsRxx) *= deltaT;
	(*itsRyy) *= deltaT;
	(*itsRzz) *= deltaT;
	(*itsRxy) *= deltaT;
	(*itsRyz) *= deltaT;
	(*itsRzx) *= deltaT;
	(*itsSx) *= deltaT;
	(*itsSy) *= deltaT;
	(*itsSz) *= deltaT;


	/// PayoffStatesSize
	//flag to avoid the calculation of the intermediate payoff
	bool isOtherPayoffs;
	//Initialisation of all the coefficients of the tridiag systemes
	InitABC2(toTimeIdx);
	//loop
	for( i=0 ; i<PayoffStatesSize ; ++i)
	{
		ARM_PayoffStates& payoffStates = states->GetPayoffStates(i);
		isOtherPayoffs = payoffStates.GetOtherPayoffsFlag();

		if(isOtherPayoffs)
		{	
			UpdateIntermediatePayoffs2( payoffStates, toTimeIdx );
		}
		std::vector<double> vec = payoffStates.GetPayoffs();
		UpdatePayoffs2( ARM_GP_Vector(vec), toTimeIdx  );
		payoffStates.SetPayoffs( vec );
	}

	// Then we remove the deltaT
	(*itsRxx) /= deltaT;
	(*itsRyy) /= deltaT;
	(*itsRzz) /= deltaT;
	(*itsRxy) /= deltaT;
	(*itsRyz) /= deltaT;
	(*itsRzx) /= deltaT;
	(*itsSx) /= deltaT;
	(*itsSy) /= deltaT;
	(*itsSz) /= deltaT;
}

////////////////////////////////////////////////////
///	Class  : ARM_PDE3FCraigSneydNumericalScheme 
///	Routine: InitABC
///	Returns: void
///	Action : Initialize the A,B & C coefficients
////////////////////////////////////////////////////

void ARM_PDE3FCraigSneydNumericalScheme::InitABC(int toTimeIdx)
{
	double deltaT=( (*GetTimeSteps())[toTimeIdx+1]-(*GetTimeSteps())[toTimeIdx] )/K_YEAR_LEN;
	double rxx;
	double ryy;
	double rzz;
	double rxy;
	double ryz;
	double rzx;
	double sx;
	double sy;
	double sz;

	int i,i_tilde,j,j_tilde,k,k_tilde,I,I_tilde,K,K_tilde;

	size_t Nx=GetNX()-1;
	size_t Ny=GetNY()-1;
	size_t Nz=GetNZ()-1;

	size_t Nx_tilde=Nx-2;
	size_t Ny_tilde=Ny-2;
	size_t Nz_tilde=Nz-2;
	
	size_t lengthU=(*itsQxx).size();
	size_t lengthUtilde=(Nx-1)*(Ny-1)*(Nz-1);
	
	//incréments directionnels
	int e1,e2,e3;
	int e_tilde2,e_tilde3;
	int e1_0,e2_0,e3_0;
	int e_tilde2_0,e_tilde3_0;
	int e_tilde2_1,e_tilde3_1;
	
	e1_0=1;
	e2_0=Nx+1;
	e3_0=(Nx+1)*(Ny+1);

		//int	e_tilde1_0=1;
	e_tilde2_0=Nx_tilde+1;
	e_tilde3_0=(Nx_tilde+1)*(Ny_tilde+1);

	e_tilde2_1=Ny_tilde+1;
	e_tilde3_1=(Nx_tilde+1)*(Ny_tilde+1);
	
	int GCG,CGG,DCG,CDG,CCG;
	int GGC,CGC,DGC,DCC,DDC,CDC,GDC,GCC,CCC;
	int GCD,CGD,DCD,CDD,CCD;
	int CCC_tilde;

	if ( lengthUtilde!=(Nx-1)*(Ny-1)*(Nz-1) ) 
	{
		throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"Utilde de mauvaise taille");
	}

	deltaT=-deltaT;//pour tenir compte que l'on est Backward et non forward
	
	//flag to avoid the recopy of the function
	//reference=Von Neumann
	double flagB1, flagB2;
	switch(itsBC)
	{
		case Dirichlet:
		{
			flagB1=2;
			flagB2=0;
		}
		break;

		case VonNeumann:
		{
			flagB1=1;
			flagB2=1;
		}
		break;

		case Belkheir:
		{
			flagB1=1;
			flagB2=1;
		}
		break;
	}

	//Premiere etape du schema ADI
		
	e1=e1_0;
	e2=e2_0;
	e3=e3_0;

	e_tilde2=e_tilde2_0;
	e_tilde3=e_tilde3_0;

	//Construction des vecteurs a,b,c,d
	for(K=0;K<lengthU;K++)
	{
		I=K%e3;
		i=I%e2;
		j=I/e2;
		k=K/e3;
				
		//indices qui vont intervenir dans les differences finies
		CCC=K;
		DCC=CCC+e1;
		GCC=CCC-e1;
		CDC=CCC+e2;
		CGC=CCC-e2;
				
		GDC=GCC+e2;
		GGC=GCC-e2;
		DDC=DCC+e2;
		DGC=DCC-e2;
		CCG=CCC-e3;
		DCG=CCG+e1;
		GCG=CCG-e1;
		CDG=CCG+e2;

		CGG=CCG-e2;
		CCD=CCC+e3;
		DCD=CCD+e1;
		GCD=CCD-e1;
		CDD=CCD+e2;
		CGD=CCD-e2;
			
		if( (i%Nx!=0)&&(j%Ny!=0)&&(k%Nz!=0) )//interrieur du domaine
		{
			//constantes du schéma
		
			rxx=2*deltaT/(pow((*GetDeltaX())[i-1],2)+pow((*GetDeltaX())[i],2));
			ryy=2*deltaT/(pow((*GetDeltaY())[j-1],2)+pow((*GetDeltaY())[j],2));
			rzz=2*deltaT/(pow((*GetDeltaZ())[k-1],2)+pow((*GetDeltaZ())[k],2));
		
			rxy=4*deltaT/( ((*GetDeltaX())[i-1]+(*GetDeltaX())[i])*((*GetDeltaY())[j-1]+(*GetDeltaY())[j]) );
			ryz=4*deltaT/( ((*GetDeltaY())[j-1]+(*GetDeltaY())[j])*((*GetDeltaZ())[k-1]+(*GetDeltaZ())[k]) );
			rzx=4*deltaT/( ((*GetDeltaZ())[k-1]+(*GetDeltaZ())[k])*((*GetDeltaX())[i-1]+(*GetDeltaX())[i]) );
			
			sx= 2*deltaT/((*GetDeltaX())[i-1]+(*GetDeltaX())[i]);
			sy= 2*deltaT/((*GetDeltaY())[j-1]+(*GetDeltaY())[j]);
			sz= 2*deltaT/((*GetDeltaZ())[k-1]+(*GetDeltaZ())[k]);
			i_tilde=i-1;
			j_tilde=j-1;
			I_tilde=e_tilde2*j_tilde+i_tilde;
			k_tilde=k-1;
			K_tilde=e_tilde3*k_tilde+I_tilde;
			
			if(i_tilde==0)
			{
				(*itsA1)[K_tilde] = 0.0; 
				(*itsB1)[K_tilde]=1 + flagB1*itsTheta1*rxx*0.5*(*itsQxx)[CCC] + flagB2*itsTheta2*sx*0.5*(*itsPx)[CCC]-itsTheta3*deltaT*(*itsO)[CCC];
				(*itsC1)[K_tilde]= -itsTheta1*rxx*0.5*(*itsQxx)[CCC]-itsTheta2*sx*0.5*(*itsPx)[CCC];
			}
			else if(i_tilde==Nx_tilde)
			{
				(*itsA1)[K_tilde]=-itsTheta1*rxx*0.5*(*itsQxx)[CCC] + itsTheta2*sx*0.5*(*itsPx)[CCC];
				(*itsB1)[K_tilde]=1 + flagB1*itsTheta1*rxx*0.5*(*itsQxx)[CCC] - flagB2*itsTheta2*sx*0.5*(*itsPx)[CCC] - itsTheta3*deltaT*(*itsO)[CCC];
				(*itsC1)[K_tilde]=0;
			}
			else
			{
				(*itsA1)[K_tilde]=-itsTheta1*rxx*0.5*(*itsQxx)[CCC] + itsTheta2*sx*0.5*(*itsPx)[CCC];
				(*itsB1)[K_tilde]=1+itsTheta1*rxx*(*itsQxx)[CCC]-itsTheta3*deltaT*(*itsO)[CCC];
				(*itsC1)[K_tilde]=-itsTheta1*rxx*0.5*(*itsQxx)[CCC]-0.5*itsTheta2*sx*(*itsPx)[CCC];
			}
		}
	}
				
	//Deuxieme etape du schema ADI
			
	//Permutations

	e2=Ny+1;
	e3=(Nx+1)*(Ny+1);

	e_tilde2=e_tilde2_1;
	e_tilde3=e_tilde3_1;

	//Construction des vecteurs a,b,c,d
	if(Ny_tilde!=0)
	{
		for(K=0;K<lengthU;K++)
		{
			I=K%e3;
			i=I%e2;
			j=I/e2;
			k=K/e3;
			
			if( (i%Ny!=0)&&(j%Nx!=0)&&(k%Nz!=0) )//interrieur du domaine
			{
				i_tilde=i-1;
				j_tilde=j-1;
				I_tilde=e_tilde2*j_tilde+i_tilde;
				k_tilde=k-1;
				K_tilde=e_tilde3*k_tilde+I_tilde;

				//indices intervenants dans les differences finies mais dans la base d'origine
				CCC=k*e3_0+i*e2_0+j;
				GCC=CCC-e2_0;
				DCC=CCC+e2_0;
				CCC_tilde=k_tilde*e_tilde3_0+i_tilde*e_tilde2_0+j_tilde;

				//constantes du schéma
				
				ryy=2*deltaT/( pow( (*GetDeltaY())[j-1],2 ) + pow( (*GetDeltaY())[j],2 ) );
				sy= 2*deltaT/( (*GetDeltaY())[j-1] + (*GetDeltaY())[j] );
									
				if (i_tilde==0)
				{
					(*itsA2)[K_tilde]=0;
					(*itsB2)[K_tilde]=1 + flagB1*itsTheta1*ryy*0.5*(*itsQyy)[CCC] + flagB2*itsTheta2*sy*0.5*(*itsPy)[CCC];
					(*itsC2)[K_tilde]= -itsTheta1*ryy*0.5*(*itsQyy)[CCC] - itsTheta2*sy*0.5*(*itsPy)[CCC];
				}
				else if (i_tilde==Ny_tilde)
				{
					(*itsA2)[K_tilde]= -itsTheta1*ryy*0.5*(*itsQyy)[CCC] + itsTheta2*sy*0.5*(*itsPy)[CCC];
					(*itsB2)[K_tilde]=1 + flagB1*itsTheta1*ryy*0.5*(*itsQyy)[CCC] - flagB2*itsTheta2*sy*0.5*(*itsPy)[CCC];
					(*itsC2)[K_tilde]=0;
				}
				else
				{
					(*itsA2)[K_tilde]=-itsTheta1*ryy*0.5*(*itsQyy)[CCC]+itsTheta2*sy*0.5*(*itsPy)[CCC];
					(*itsB2)[K_tilde]=1+itsTheta1*ryy*(*itsQyy)[CCC];
					(*itsC2)[K_tilde]=-itsTheta1*ryy*0.5*(*itsQyy)[CCC]-itsTheta2*sy*0.5*(*itsPy)[CCC];
				}
			}
		}
	}
	else
	{
		for(K=0;K<lengthU;K++)
		{
			I=K%e3;
			i=I%e2;
			j=I/e2;
			k=K/e3;
			
			if( (i%Ny!=0)&&(j%Nx!=0)&&(k%Nz!=0) )//interrieur du domaine
			{
				i_tilde=i-1;
				j_tilde=j-1;
				I_tilde=e_tilde2*j_tilde+i_tilde;
				k_tilde=k-1;
				K_tilde=e_tilde3*k_tilde+I_tilde;

				//indices intervenants dans les differences finies mais dans la base d'origine
				CCC=k*e3_0+i*e2_0+j;
				GCC=CCC-e2_0;
				DCC=CCC+e2_0;
				CCC_tilde=k_tilde*e_tilde3_0+i_tilde*e_tilde2_0+j_tilde;

				//constantes du schéma
				
				ryy=2*deltaT/( pow( (*GetDeltaY())[j-1],2 ) + pow( (*GetDeltaY())[j],2 ) );
				sy= 2*deltaT/( (*GetDeltaY())[j-1] + (*GetDeltaY())[j] );
									
				(*itsA2)[K_tilde]=0;
				(*itsB2)[K_tilde]=1;
				(*itsC2)[K_tilde]=0;
			}
		}

	}

	//Troisieme et derniere etape du schema ADI
			
	//permutations
	e2=Nz+1;
	e3=(Nx+1)*(Nz+1);

	e_tilde2=Nz_tilde+1;
	e_tilde3=(Nx_tilde+1)*(Nz_tilde+1);

	//Construction des vecteurs a,b,c,d
	if (Nz_tilde!=0)
	{
		for(K=0;K<lengthU;K++)
		{
			I=K%e3;
			i=I%e2;
			j=I/e2;
			k=K/e3;

			if( (i%Nz!=0)&&(j%Nx!=0)&&(k%Ny!=0) )//interrieur du domaine
			{
				i_tilde=i-1;
				j_tilde=j-1;
				I_tilde=e_tilde2*j_tilde+i_tilde;
				k_tilde=k-1;
				K_tilde=e_tilde3*k_tilde+I_tilde;

				//indices intervenants dans les differences finies
				CCC=i*e3_0+k*e2_0+j;
				GCC=CCC-e3_0;
				DCC=CCC+e3_0;
				CCC_tilde=i_tilde*e_tilde3_0+j_tilde*e_tilde2_1+k_tilde;

				//constantes du schéma
		
				rzz=2*deltaT/( pow( (*GetDeltaZ())[k-1],2 ) + pow( (*GetDeltaZ())[k],2 ) );
				sz= 2*deltaT/( (*GetDeltaZ())[k-1] + (*GetDeltaZ())[k] );

				if(i_tilde==0)
				{
					(*itsA3)[K_tilde]=0;
					(*itsB3)[K_tilde]=1 + flagB1*itsTheta1*rzz*0.5*(*itsQzz)[CCC] + flagB2*itsTheta2*sz*0.5*(*itsPz)[CCC];
					(*itsC3)[K_tilde]=-itsTheta1*rzz*0.5*(*itsQzz)[CCC]-itsTheta2*sz*0.5*(*itsPz)[CCC];
				}
				else if(i_tilde==Nz_tilde)
				{
					(*itsA3)[K_tilde]= -itsTheta1*rzz*0.5*(*itsQzz)[CCC] + itsTheta2*sz*0.5*(*itsPz)[CCC];
					(*itsB3)[K_tilde]=1 + flagB1*itsTheta1*rzz*0.5*(*itsQzz)[CCC] - flagB2*itsTheta2*sz*0.5*(*itsPz)[CCC];
					(*itsC3)[K_tilde]=0;
				}
				else
				{
					(*itsA3)[K_tilde]=-itsTheta1*rzz*0.5*(*itsQzz)[CCC] + itsTheta2*sz*0.5*(*itsPz)[CCC];
					(*itsB3)[K_tilde]=1 + itsTheta1*rzz*(*itsQzz)[CCC];
					(*itsC3)[K_tilde]=-itsTheta1*rzz*0.5*(*itsQzz)[CCC] - itsTheta2*sz*0.5*(*itsPz)[CCC];
				}
			}
		}
	}
	else
	{
		for(K=0;K<lengthU;K++)
		{
			I=K%e3;
			i=I%e2;
			j=I/e2;
			k=K/e3;

			if( (i%Nz!=0)&&(j%Nx!=0)&&(k%Ny!=0) )//interrieur du domaine
			{
				i_tilde=i-1;
				j_tilde=j-1;
				I_tilde=e_tilde2*j_tilde+i_tilde;
				k_tilde=k-1;
				K_tilde=e_tilde3*k_tilde+I_tilde;

				//indices intervenants dans les differences finies
				CCC=i*e3_0+k*e2_0+j;
				GCC=CCC-e3_0;
				DCC=CCC+e3_0;
				CCC_tilde=i_tilde*e_tilde3_0+j_tilde*e_tilde2_1+k_tilde;

				(*itsA3)[K_tilde]=0;
				(*itsB3)[K_tilde]=1;
				(*itsC3)[K_tilde]=0;
			}
		}
	}					
}


////////////////////////////////////////////////////
///	Class  : ARM_PDE3FCraigSneydNumericalScheme 
///	Routine: InitABC2
///	Returns: void
///	Action : Initialize the A,B & C coefficients
////////////////////////////////////////////////////


void ARM_PDE3FCraigSneydNumericalScheme::InitABC2(int toTimeIdx)
{
	double deltaT=( (*GetTimeSteps())[toTimeIdx+1]-(*GetTimeSteps())[toTimeIdx] )/K_YEAR_LEN;
	// Backward To Forward induction;
	deltaT=-deltaT;

	size_t Nx=GetNX()-1;
	size_t Ny=GetNY()-1;
	size_t Nz=GetNZ()-1;
	
	size_t lengthU=(*itsQxx).size();
	
	//incréments directionnels
	
	size_t e1_xyz=1;
	size_t e2_xyz=Nx+1;
	size_t e3_xyz=(Nx+1)*(Ny+1);

	size_t e1_yxz=1;
	size_t e2_yxz=Ny+1;
	size_t e3_yxz=(Ny+1)*(Nx+1);

	size_t e1_zxy=1;
	size_t e2_zxy=Nz+1;
	size_t e3_zxy=(Nz+1)*(Nx+1);

	// Temporary variable for the permutation
	size_t K_xyz, I_xyz;
	size_t i_xyz, j_xyz, k_xyz;

	if ( lengthU!=(Nx+1)*(Ny+1)*(Nz+1) ) 
	{
		throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"U de mauvaise taille");
	}

	//flag to avoid the recopy of the function
	//reference=Von Neumann
	double flagB1, flagB2;
	switch(itsBC)
	{
		case Dirichlet:
		{
			flagB1=2;
			flagB2=0;
		}
		break;

		case VonNeumann:
		{
			flagB1=1;
			flagB2=1;
		}
		break;

		case Belkheir:
		{
			flagB1=1;
			flagB2=1;
		}
		break;
	}

	ARM_GP_T_Vector<bool>::iterator IsBoundit1 = itsIsBound.begin();
	std::vector<double>::iterator Qxxit = itsQxx->begin();
	std::vector<double>::iterator Qyyit = itsQyy->begin();
	std::vector<double>::iterator Qzzit = itsQzz->begin();
	std::vector<double>::iterator Pxit = itsPx->begin();
	std::vector<double>::iterator Pyit = itsPy->begin();
	std::vector<double>::iterator Pzit = itsPz->begin();
	std::vector<double>::iterator Oit = itsO->begin();
	std::vector<double>::iterator Rxxit = itsRxx->begin();
	std::vector<double>::iterator Ryyit = itsRyy->begin();
	std::vector<double>::iterator Rzzit = itsRzz->begin();
	std::vector<double>::iterator Sxit = itsSx->begin();
	std::vector<double>::iterator Syit = itsSy->begin();
	std::vector<double>::iterator Szit = itsSz->begin();
	std::vector<double>::iterator A1it = itsA1->begin();
	std::vector<double>::iterator B1it = itsB1->begin();
	std::vector<double>::iterator C1it = itsC1->begin();
	std::vector<double>::iterator A2it = itsA2->begin();
	std::vector<double>::iterator B2it = itsB2->begin();
	std::vector<double>::iterator C2it = itsC2->begin();
	std::vector<double>::iterator A3it = itsA3->begin();
	std::vector<double>::iterator B3it = itsB3->begin();
	std::vector<double>::iterator C3it = itsC3->begin();

	size_t dec12=0;
	size_t dec13=0;

	//Construction des vecteurs a,b,c,d
	for(K_xyz=0;K_xyz<lengthU;K_xyz++)
	{		
		if( !(*IsBoundit1) )//interrieur du domaine
		{		
			k_xyz = K_xyz/e3_xyz;
			I_xyz = K_xyz%e3_xyz;
			j_xyz = I_xyz/e2_xyz;
			i_xyz = I_xyz%e2_xyz;
		
			dec12 = k_xyz*e3_yxz+i_xyz*e2_yxz+j_xyz*e1_yxz;
			dec13 = j_xyz*e3_zxy+i_xyz*e2_zxy+k_xyz*e1_zxy;

			// Base d'écriture yxz
			A2it=itsA2->begin()+dec12;
			B2it=itsB2->begin()+dec12;
			C2it=itsC2->begin()+dec12;

			// Base d'écriture zxy
			A3it=itsA3->begin()+dec13;
			B3it=itsB3->begin()+dec13;
			C3it=itsC3->begin()+dec13;

			(*A1it)	=	-itsTheta1*(*Rxxit)*0.5*(*Qxxit) + itsTheta2*(*Sxit)*0.5*(*Pxit);
			(*B1it)	=	1 + itsTheta1*(*Rxxit)*(*Qxxit)  - itsTheta3*deltaT*(*Oit);
			(*C1it)	=	-itsTheta1*(*Rxxit)*0.5*(*Qxxit) - itsTheta2*(*Sxit)*0.5*(*Pxit);

			(*A2it)	=	-itsTheta1*(*Ryyit)*0.5*(*Qyyit) + itsTheta2*(*Syit)*0.5*(*Pyit);
			(*B2it)	=	1 + itsTheta1*(*Ryyit)*(*Qyyit);
			(*C2it)	=	-itsTheta1*(*Ryyit)*0.5*(*Qyyit) - itsTheta2*(*Syit)*0.5*(*Pyit);

			(*A3it)	=	-itsTheta1*(*Rzzit)*0.5*(*Qzzit) + itsTheta2*(*Szit)*0.5*(*Pzit);
			(*B3it)	=	1 + itsTheta1*(*Rzzit)*(*Qzzit);
			(*C3it)	=	-itsTheta1*(*Rzzit)*0.5*(*Qzzit) - itsTheta2*(*Szit)*0.5*(*Pzit);

			if(i_xyz==1)
			{
				(*A1it)	=	0.0; 
				(*B1it)	=	1+ flagB1*itsTheta1*(*Rxxit)*0.5*(*Qxxit) + flagB2*itsTheta2*(*Sxit)*0.5*(*Pxit)-itsTheta3*deltaT*(*Oit);
			}
			if(i_xyz==Nx-1)
			{
				(*B1it)	=	1 + flagB1*itsTheta1*(*Rxxit)*0.5*(*Qxxit) - flagB2*itsTheta2*(*Sxit)*0.5*(*Pxit) - itsTheta3*deltaT*(*Oit);
				(*C1it)	=	0;
			}
			if(j_xyz==1)
			{
				(*A2it)	=	0.0; 
				(*B2it)	=	1+ flagB1*itsTheta1*(*Ryyit)*0.5*(*Qyyit) + flagB2*itsTheta2*(*Syit)*0.5*(*Pyit);
			}
			if(j_xyz==Ny-1)
			{
				(*B2it)	=	1 + flagB1*itsTheta1*(*Ryyit)*0.5*(*Qyyit) - flagB2*itsTheta2*(*Syit)*0.5*(*Pyit);
				(*C2it)	=	0;
			}
			if(k_xyz==1)
			{
				(*A3it)	=	0.0; 
				(*B3it)	=	1+ flagB1*itsTheta1*(*Rzzit)*0.5*(*Qzzit) + flagB2*itsTheta2*(*Szit)*0.5*(*Pzit);
			}
			if(k_xyz==Nz-1)
			{
				(*B3it)	=	1 + flagB1*itsTheta1*(*Rzzit)*0.5*(*Qzzit) - flagB2*itsTheta2*(*Szit)*0.5*(*Pzit);
				(*C3it)	=	0;
			}
			//pour la stabilite de la degenerecence 
			if(Nx==2)
			{
				(*A1it)	=	0.0; 
				(*B1it)	=	1.0;
				(*C1it)	=	0.0;
			}
			if(Ny==2)
			{
				(*A2it)	=	0.0; 
				(*B2it)	=	1.0;
				(*C2it)	=	0.0;
			}
			if(Nz==2)
			{
				(*A3it)	=	0.0; 
				(*B3it)	=	1.0;
				(*C3it)	=	0.0;
			}
		}

		// Base de lecture xyz
		Qxxit++;
		Pxit++;
		Rxxit++;
		Sxit++;
	
		Qyyit++;
		Pyit++;
		Ryyit++;
		Syit++;
		
		Qzzit++;
		Pzit++;
		Rzzit++;
		Szit++;

		Oit++;

		// Base d'écriture xyz
		A1it++;
		B1it++;
		C1it++;

		//incrementation
		IsBoundit1++;
	}
}


////////////////////////////////////////////////////
///	Class  : UpdatePayoffs 
///	Routine: UpdatePayoffs
///	Returns: void
///	Action : Backward a payoff vector
////////////////////////////////////////////////////


void ARM_PDE3FCraigSneydNumericalScheme::UpdatePayoffs(ARM_GP_Vector& U, int toTimeIdx)
{
	double deltaT=( (*GetTimeSteps())[toTimeIdx+1]-(*GetTimeSteps())[toTimeIdx] )/K_YEAR_LEN;
	double rxx;
	double ryy;
	double rzz;
	double rxy;
	double ryz;
	double rzx;
	double sx;
	double sy;
	double sz;

	int i,i_tilde,j,j_tilde,k,k_tilde,I,I_tilde,K,K_tilde;

	size_t Nx=GetNX()-1;
	size_t Ny=GetNY()-1;
	size_t Nz=GetNZ()-1;

	size_t Nx_tilde=Nx-2;
	size_t Ny_tilde=Ny-2;
	size_t Nz_tilde=Nz-2;
	
	unsigned int lengthU=U.size();
	unsigned int lengthUtilde=(Nx-1)*(Ny-1)*(Nz-1);
	
	//incréments directionnels
	int e1,e2,e3;
	int e_tilde2,e_tilde3;
	int e1_0,e2_0,e3_0;
	int e_tilde2_0,e_tilde3_0;
	int e_tilde2_1,e_tilde3_1;
	
	e1_0=1;
	e2_0=Nx+1;
	e3_0=(Nx+1)*(Ny+1);

		int	e_tilde1_0=1;
	e_tilde2_0=Nx_tilde+1;
	e_tilde3_0=(Nx_tilde+1)*(Ny_tilde+1);

	e_tilde2_1=Ny_tilde+1;
	e_tilde3_1=(Nx_tilde+1)*(Ny_tilde+1);
	
	int GCG,CGG,DCG,CDG,CCG;
	int GGC,CGC,DGC,DCC,DDC,CDC,GDC,GCC,CCC;
	int GCD,CGD,DCD,CDD,CCD;
	int CCC_tilde;

	if ( lengthUtilde!=(Nx-1)*(Ny-1)*(Nz-1) ) 
	{
		throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"Utilde de mauvaise taille");
	}

	deltaT=-deltaT;//pour tenir compte que l'on est Backward et non forward
	
	int ind=1;
	if( (itsBC==Dirichlet) || (itsBC==VonNeumann) )
	{
		ind=0;
	}
	switch(ind)
	{
		case 0:
		{
			//Premiere etape du schema ADI
			
			e1=e1_0;
			e2=e2_0;
			e3=e3_0;

			e_tilde2=e_tilde2_0;
			e_tilde3=e_tilde3_0;

			//Construction du vecteur d
			for(K=0;K<lengthU;K++)
			{
				I=K%e3;
				i=I%e2;
				j=I/e2;
				k=K/e3;
				
				//indices qui vont intervenir dans les differences finies
				CCC=K;
				DCC=CCC+e1;
				GCC=CCC-e1;
				CDC=CCC+e2;
				CGC=CCC-e2;
					
				GDC=GCC+e2;
				GGC=GCC-e2;
				DDC=DCC+e2;
				DGC=DCC-e2;

				CCG=CCC-e3;
				DCG=CCG+e1;
				GCG=CCG-e1;
				CDG=CCG+e2;
				CGG=CCG-e2;

				CCD=CCC+e3;
				DCD=CCD+e1;
				GCD=CCD-e1;
				CDD=CCD+e2;
				CGD=CCD-e2;
			
				if( (i%Nx!=0)&&(j%Ny!=0)&&(k%Nz!=0) )//interrieur du domaine
				{
					//constantes du schéma
				
					rxx=2*deltaT/(pow((*GetDeltaX())[i-1],2)+pow((*GetDeltaX())[i],2));
					ryy=2*deltaT/(pow((*GetDeltaY())[j-1],2)+pow((*GetDeltaY())[j],2));
					rzz=2*deltaT/(pow((*GetDeltaZ())[k-1],2)+pow((*GetDeltaZ())[k],2));
				
					rxy=4*deltaT/( ((*GetDeltaX())[i-1]+(*GetDeltaX())[i])*((*GetDeltaY())[j-1]+(*GetDeltaY())[j]) );
					ryz=4*deltaT/( ((*GetDeltaY())[j-1]+(*GetDeltaY())[j])*((*GetDeltaZ())[k-1]+(*GetDeltaZ())[k]) );
					rzx=4*deltaT/( ((*GetDeltaZ())[k-1]+(*GetDeltaZ())[k])*((*GetDeltaX())[i-1]+(*GetDeltaX())[i]) );
				
					sx= 2*deltaT/((*GetDeltaX())[i-1]+(*GetDeltaX())[i]);
					sy= 2*deltaT/((*GetDeltaY())[j-1]+(*GetDeltaY())[j]);
					sz= 2*deltaT/((*GetDeltaZ())[k-1]+(*GetDeltaZ())[k]);

					i_tilde=i-1;
					j_tilde=j-1;
					I_tilde=e_tilde2*j_tilde+i_tilde;
					k_tilde=k-1;
					K_tilde=e_tilde3*k_tilde+I_tilde;
					
					(*itsD)[K_tilde]=U[CCC]
							+(1-itsTheta1)*rxx*0.5*(*itsPrevQxx)[CCC]*( U[GCC] -2*U[CCC] + U[DCC] )
							+ryy*0.5*(*itsPrevQyy)[CCC]*( U[CGC] - 2*U[CCC] + U[CDC] )
							+rzz*0.5*(*itsPrevQzz)[CCC]*( U[CCG] - 2*U[CCC] + U[CCD] )
							+rxy*0.25*(*itsPrevQxy)[CCC]*( U[GGC] - U[DGC] - U[GDC] + U[DDC] )
							+ryz*0.25*(*itsPrevQyz)[CCC]*( U[CGG] - U[CDG] - U[CGD] + U[CDD] )
							+rzx*0.25*(*itsPrevQzx)[CCC]*( U[GCG] - U[DCG] - U[GCD] + U[DCD] )
							+(1-itsTheta2)*sx*0.5*(*itsPrevPx)[CCC]*( -U[GCC] + U[DCC] )
							+sy*0.5*(*itsPrevPy)[CCC]*( -U[CGC] + U[CDC] )
							+sz*0.5*(*itsPrevPz)[CCC]*( -U[CCG] + U[CCD] )
							+(1-itsTheta3)*deltaT*(*itsPrevO)[CCC]*U[CCC];
				}
			}
			
			tridiag(*itsGam,*itsA1,*itsB1,*itsC1,*itsD,*itsUtilde,lengthUtilde);
					
			//Deuxieme etape du schema ADI
			
			//Permutations

			e2=Ny+1;
			e3=(Nx+1)*(Ny+1);

			e_tilde2=e_tilde2_1;
			e_tilde3=e_tilde3_1;

		
			for(K=0;K<lengthU;K++)
			{
				I=K%e3;
				i=I%e2;
				j=I/e2;
				k=K/e3;
			
				if( (i%Ny!=0)&&(j%Nx!=0)&&(k%Nz!=0) )//interrieur du domaine
				{
					i_tilde=i-1;
					j_tilde=j-1;
					I_tilde=e_tilde2*j_tilde+i_tilde;
					k_tilde=k-1;
					K_tilde=e_tilde3*k_tilde+I_tilde;

					//indices intervenants dans les differences finies mais dans la base d'origine
					CCC=k*e3_0+i*e2_0+j;
					GCC=CCC-e2_0;
					DCC=CCC+e2_0;
					CCC_tilde=k_tilde*e_tilde3_0+i_tilde*e_tilde2_0+j_tilde;

					//constantes du schéma
				
					ryy=2*deltaT/( pow( (*GetDeltaY())[j-1],2 ) + pow( (*GetDeltaY())[j],2 ) );
					sy= 2*deltaT/( (*GetDeltaY())[j-1] + (*GetDeltaY())[j] );
				
					(*itsD)[K_tilde]=(*itsUtilde)[CCC_tilde]
							-itsTheta1*ryy*0.5*( (*itsPrevQyy)[CCC]*U[GCC] - 2*(*itsPrevQyy)[CCC]*U[CCC] + (*itsPrevQyy)[CCC]*U[DCC])
							-itsTheta2*sy*0.5*( -(*itsPrevPy)[CCC]*U[GCC] + (*itsPrevPy)[CCC]*U[DCC]);
				}
			}

			tridiag(*itsGam,*itsA2,*itsB2,*itsC2,*itsD,*itsUtilde,lengthUtilde);

			//Troisieme et derniere etape du schema ADI
				
			//permutations
			e2=Nz+1;
			e3=(Nx+1)*(Nz+1);

			e_tilde2=Nz_tilde+1;
			e_tilde3=(Nx_tilde+1)*(Nz_tilde+1);

			//Construction du vecteur d
			for(K=0;K<lengthU;K++)
			{
				I=K%e3;
				i=I%e2;
				j=I/e2;
				k=K/e3;

				if( (i%Nz!=0)&&(j%Nx!=0)&&(k%Ny!=0) )//interrieur du domaine
				{
					i_tilde=i-1;
					j_tilde=j-1;
					I_tilde=e_tilde2*j_tilde+i_tilde;
					k_tilde=k-1;
					K_tilde=e_tilde3*k_tilde+I_tilde;

					//indices intervenants dans les differences finies
					CCC=i*e3_0+k*e2_0+j;
					GCC=CCC-e3_0;
					DCC=CCC+e3_0;
					CCC_tilde=i_tilde*e_tilde3_0+j_tilde*e_tilde2_1+k_tilde;

					//constantes du schéma
					rzz=2*deltaT/( pow( (*GetDeltaZ())[k-1],2 ) + pow( (*GetDeltaZ())[k],2 ) );
					sz= 2*deltaT/( (*GetDeltaZ())[k-1] + (*GetDeltaZ())[k] );
					
					(*itsD)[K_tilde]=(*itsUtilde)[CCC_tilde]
									-itsTheta1*rzz*0.5*(*itsPrevQzz)[CCC]*( U[GCC] -2*U[CCC] + U[DCC] )
									-itsTheta2*sz*0.5*(*itsPrevPz)[CCC]*( -U[GCC] + U[DCC] );
				}
			}
										
  			tridiag(*itsGam,*itsA3,*itsB3,*itsC3,*itsD,*itsUtilde,lengthUtilde);
	
			//permutations finales
			InvPermut(*itsUtilde,Nx_tilde,Ny_tilde,Nz_tilde,*itsUtilde_temp);

			//Construction du vecteur U_temp 
			Reconstruct(*itsUtilde_temp,U,*itsU_new);
		}
		break;

		case 1:
		{
			//Premiere etape du schema ADI
			
			e1=e1_0;
			e2=e2_0;
			e3=e3_0;

			e_tilde2=e_tilde2_0;
			e_tilde3=e_tilde3_0;

			//Construction des vecteurs a,b,c,d
			for(K=0;K<lengthU;K++)
			{
				I=K%e3;
				i=I%e2;
				j=I/e2;
				k=K/e3;
				
				//indices qui vont intervenir dans les differences finies
				CCC=K;
				DCC=CCC+e1;
				GCC=CCC-e1;
				CDC=CCC+e2;
				CGC=CCC-e2;
					
				GDC=GCC+e2;
				GGC=GCC-e2;
				DDC=DCC+e2;
				DGC=DCC-e2;

				CCG=CCC-e3;
				DCG=CCG+e1;
				GCG=CCG-e1;
				CDG=CCG+e2;
				CGG=CCG-e2;

				CCD=CCC+e3;
				DCD=CCD+e1;
				GCD=CCD-e1;
				CDD=CCD+e2;
				CGD=CCD-e2;
			
				if( (i%Nx!=0)&&(j%Ny!=0)&&(k%Nz!=0) )//interrieur du domaine
				{
					//constantes du schéma
				
					rxx=2*deltaT/(pow((*GetDeltaX())[i-1],2)+pow((*GetDeltaX())[i],2));
					ryy=2*deltaT/(pow((*GetDeltaY())[j-1],2)+pow((*GetDeltaY())[j],2));
					rzz=2*deltaT/(pow((*GetDeltaZ())[k-1],2)+pow((*GetDeltaZ())[k],2));
				
					rxy=4*deltaT/( ((*GetDeltaX())[i-1]+(*GetDeltaX())[i])*((*GetDeltaY())[j-1]+(*GetDeltaY())[j]) );
					ryz=4*deltaT/( ((*GetDeltaY())[j-1]+(*GetDeltaY())[j])*((*GetDeltaZ())[k-1]+(*GetDeltaZ())[k]) );
					rzx=4*deltaT/( ((*GetDeltaZ())[k-1]+(*GetDeltaZ())[k])*((*GetDeltaX())[i-1]+(*GetDeltaX())[i]) );
				
					sx= 2*deltaT/((*GetDeltaX())[i-1]+(*GetDeltaX())[i]);
					sy= 2*deltaT/((*GetDeltaY())[j-1]+(*GetDeltaY())[j]);
					sz= 2*deltaT/((*GetDeltaZ())[k-1]+(*GetDeltaZ())[k]);

					i_tilde=i-1;
					j_tilde=j-1;
					I_tilde=e_tilde2*j_tilde+i_tilde;
					k_tilde=k-1;
					K_tilde=e_tilde3*k_tilde+I_tilde;
					
					if(i_tilde==0)
					{
						(*itsD)[K_tilde]=U[CCC]
								+(1-itsTheta1)*rxx*0.5*( (*itsPrevQxx)[CCC]*U[GCC] -2*(*itsPrevQxx)[CCC]*U[CCC] + (*itsPrevQxx)[CCC]*U[DCC] )
								+ryy*0.5*( (*itsPrevQyy)[CCC]*U[CGC] - 2*(*itsPrevQyy)[CCC]*U[CCC] + (*itsPrevQyy)[CCC]*U[CDC] )
								+rzz*0.5*( (*itsPrevQzz)[CCC]*U[CCG] - 2*(*itsPrevQzz)[CCC]*U[CCC] + (*itsPrevQzz)[CCC]*U[CCD] )
								+rxy*0.25*( (*itsPrevQxy)[CCC]*U[GGC] - (*itsPrevQxy)[CCC]*U[DGC] - (*itsPrevQxy)[CCC]*U[GDC] + (*itsPrevQxy)[CCC]*U[DDC] )
								+ryz*0.25*( (*itsPrevQyz)[CCC]*U[CGG] - (*itsPrevQyz)[CCC]*U[CDG] - (*itsPrevQyz)[CCC]*U[CGD] + (*itsPrevQyz)[CCC]*U[CDD] )
								+rzx*0.25*( (*itsPrevQzx)[CCC]*U[GCG] - (*itsPrevQzx)[CCC]*U[DCG] - (*itsPrevQzx)[CCC]*U[GCD] + (*itsPrevQzx)[CCC]*U[DCD] )
								+(1-itsTheta2)*sx*0.5*( -(*itsPrevPx)[CCC]*U[GCC] + (*itsPrevPx)[CCC]*U[DCC])
								+sy*0.5*( -(*itsPrevPy)[CCC]*U[CGC] + (*itsPrevPy)[CCC]*U[CDC] )
								+sz*0.5*( -(*itsPrevPz)[CCC]*U[CCG]+(*itsPrevPz)[CCC]*U[CCD] )
								+(1-itsTheta3)*deltaT*(*itsPrevO)[CCC]*U[CCC]
									+itsTheta1*rxx*0.5*(*itsPrevQxx)[CCC]*(U[CCC]-U[GCC])//
									+itsTheta2*sx*0.5*(*itsPrevPx)[CCC]*(U[CCC]-U[GCC]);//

					}
					else if(i_tilde==Nx_tilde)
					{
						(*itsD)[K_tilde]=U[CCC]
								+(1-itsTheta1)*rxx*0.5*( (*itsPrevQxx)[CCC]*U[GCC] - 2*(*itsPrevQxx)[CCC]*U[CCC] + (*itsPrevQxx)[CCC]*U[DCC] )
								+ryy*0.5*( (*itsPrevQyy)[CCC]*U[CGC] - 2*(*itsPrevQyy)[CCC]*U[CCC] + (*itsPrevQyy)[CCC]*U[CDC] )
								+rzz*0.5*( (*itsPrevQzz)[CCC]*U[CCG] - 2*(*itsPrevQzz)[CCC]*U[CCC] + (*itsPrevQzz)[CCC]*U[CCD] )
								+rxy*0.25*( (*itsPrevQxy)[CCC]*U[GGC] - (*itsPrevQxy)[CCC]*U[DGC] - (*itsPrevQxy)[CCC]*U[GDC] + (*itsPrevQxy)[CCC]*U[DDC] )
								+ryz*0.25*( (*itsPrevQyz)[CCC]*U[CGG] - (*itsPrevQyz)[CCC]*U[CDG] - (*itsPrevQyz)[CCC]*U[CGD] + (*itsPrevQyz)[CCC]*U[CDD] )
								+rzx*0.25*( (*itsPrevQzx)[CCC]*U[GCG] - (*itsPrevQzx)[CCC]*U[DCG] - (*itsPrevQzx)[CCC]*U[GCD] + (*itsPrevQzx)[CCC]*U[DCD] )
								+(1-itsTheta2)*sx*0.5*( -(*itsPrevPx)[CCC]*U[GCC] + (*itsPrevPx)[CCC]*U[DCC] )
								+sy*0.5*( -(*itsPrevPy)[CCC]*U[CGC] + (*itsPrevPy)[CCC]*U[CDC] )
								+sz*0.5*( -(*itsPrevPz)[CCC]*U[CCG] + (*itsPrevPz)[CCC]*U[CCD] )
								+(1-itsTheta3)*deltaT*(*itsPrevO)[CCC]*U[CCC]
									+itsTheta1*rxx*0.5*(*itsPrevQxx)[CCC]*(U[CCC]-U[DCC])
									-itsTheta2*sx*0.5*(*itsPrevPx)[CCC]*(U[CCC]-U[DCC]);

					}
					else
					{
						(*itsD)[K_tilde]=U[CCC]
							+(1-itsTheta1)*rxx*0.5*( (*itsPrevQxx)[CCC]*U[GCC] -2*(*itsPrevQxx)[CCC]*U[CCC] + (*itsPrevQxx)[CCC]*U[DCC] )
							+ryy*0.5*( (*itsPrevQyy)[CCC]*U[CGC] - 2*(*itsPrevQyy)[CCC]*U[CCC] + (*itsPrevQyy)[CCC]*U[CDC] )
							+rzz*0.5*( (*itsPrevQzz)[CCC]*U[CCG] - 2*(*itsPrevQzz)[CCC]*U[CCC] + (*itsPrevQzz)[CCC]*U[CCD] )
							+rxy*0.25*( (*itsPrevQxy)[CCC]*U[GGC] - (*itsPrevQxy)[CCC]*U[DGC] - (*itsPrevQxy)[CCC]*U[GDC] + (*itsPrevQxy)[CCC]*U[DDC] )
							+ryz*0.25*( (*itsPrevQyz)[CCC]*U[CGG] - (*itsPrevQyz)[CCC]*U[CDG] - (*itsPrevQyz)[CCC]*U[CGD] + (*itsPrevQyz)[CCC]*U[CDD] )
							+rzx*0.25*( (*itsPrevQzx)[CCC]*U[GCG] - (*itsPrevQzx)[CCC]*U[DCG] - (*itsPrevQzx)[CCC]*U[GCD] + (*itsPrevQzx)[CCC]*U[DCD] )
							+(1-itsTheta2)*sx*0.5*( -(*itsPrevPx)[CCC]*U[GCC] + (*itsPrevPx)[CCC]*U[DCC] )
							+sy*0.5*( -(*itsPrevPy)[CCC]*U[CGC] + (*itsPrevPy)[CCC]*U[CDC] )
							+sz*0.5*( -(*itsPrevPz)[CCC]*U[CCG] + (*itsPrevPz)[CCC]*U[CCD] )
							+(1-itsTheta3)*deltaT*(*itsPrevO)[CCC]*U[CCC];
					}
				}
			}
			tridiag(*itsGam,*itsA1,*itsB1,*itsC1,*itsD,*itsUtilde,lengthUtilde);
					
			//Deuxieme etape du schema ADI
			
			//Permutations

			e2=Ny+1;
			e3=(Nx+1)*(Ny+1);

			e_tilde2=e_tilde2_1;
			e_tilde3=e_tilde3_1;

			//Construction du vecteur d
			if(Ny_tilde!=0)
			{
				for(K=0;K<lengthU;K++)
				{
					I=K%e3;
					i=I%e2;
					j=I/e2;
					k=K/e3;
			
					if( (i%Ny!=0)&&(j%Nx!=0)&&(k%Nz!=0) )//interrieur du domaine
					{
						i_tilde=i-1;
						j_tilde=j-1;
						I_tilde=e_tilde2*j_tilde+i_tilde;
						k_tilde=k-1;
						K_tilde=e_tilde3*k_tilde+I_tilde;

						//indices intervenants dans les differences finies mais dans la base d'origine
						CCC=k*e3_0+i*e2_0+j;
						GCC=CCC-e2_0;
						DCC=CCC+e2_0;
						CCC_tilde=k_tilde*e_tilde3_0+i_tilde*e_tilde2_0+j_tilde;

						//constantes du schéma
						ryy=2*deltaT/( pow( (*GetDeltaY())[j-1],2 ) + pow( (*GetDeltaY())[j],2 ) );
						sy= 2*deltaT/( (*GetDeltaY())[j-1] + (*GetDeltaY())[j] );
						
						
						if (i_tilde==0)
						{
							(*itsD)[K_tilde]=(*itsUtilde)[CCC_tilde]
									-itsTheta1*ryy*0.5*( (*itsPrevQyy)[CCC]*U[GCC] - 2*(*itsPrevQyy)[CCC]*U[CCC] + (*itsPrevQyy)[CCC]*U[DCC])
									-itsTheta2*sy*0.5*( -(*itsPrevPy)[CCC]*U[GCC] +(*itsPrevPy)[CCC]*U[DCC] )
										+itsTheta1*ryy*0.5*(*itsPrevQyy)[CCC]*( U[CCC]-U[GCC] )
										+itsTheta2*sy*0.5*(*itsPrevQyy)[CCC]*( U[CCC]-U[GCC] );

						}
						else if (i_tilde==Ny_tilde)
						{
							(*itsD)[K_tilde]=(*itsUtilde)[CCC_tilde]
									-itsTheta1*ryy*0.5*( (*itsPrevQyy)[CCC]*U[GCC] - 2*(*itsPrevQyy)[CCC]*U[CCC] + (*itsPrevQyy)[CCC]*U[DCC] )
									-itsTheta2*sy*0.5*( -(*itsPrevPy)[CCC]*U[GCC] + (*itsPrevPy)[CCC]*U[DCC] )
										+itsTheta1*ryy*0.5*(*itsPrevQyy)[CCC]*( U[CCC]-U[DCC] )
										-itsTheta2*sy*0.5*(*itsPrevQyy)[CCC]*( U[CCC]-U[DCC] );
						}
						else
						{
							(*itsD)[K_tilde]=(*itsUtilde)[CCC_tilde]
									-itsTheta1*ryy*0.5*( (*itsPrevQyy)[CCC]*U[GCC] - 2*(*itsPrevQyy)[CCC]*U[CCC] + (*itsPrevQyy)[CCC]*U[DCC] )
									-itsTheta2*sy*0.5*( -(*itsPrevPy)[CCC]*U[GCC] + (*itsPrevPy)[CCC]*U[DCC] );
						}
					}
				}
			}
			else
			{
				for(K=0;K<lengthU;K++)
				{
					I=K%e3;
					i=I%e2;
					j=I/e2;
					k=K/e3;
			
					if( (i%Ny!=0)&&(j%Nx!=0)&&(k%Nz!=0) )//interrieur du domaine
					{
						i_tilde=i-1;
						j_tilde=j-1;
						I_tilde=e_tilde2*j_tilde+i_tilde;
						k_tilde=k-1;
						K_tilde=e_tilde3*k_tilde+I_tilde;

						//indices intervenants dans les differences finies mais dans la base d'origine
						CCC=k*e3_0+i*e2_0+j;
						GCC=CCC-e2_0;
						DCC=CCC+e2_0;
						CCC_tilde=k_tilde*e_tilde3_0+i_tilde*e_tilde2_0+j_tilde;

						//constantes du schéma
						ryy=2*deltaT/( pow( (*GetDeltaY())[j-1],2 ) + pow( (*GetDeltaY())[j],2 ) );
						sy= 2*deltaT/( (*GetDeltaY())[j-1] + (*GetDeltaY())[j] );
									
						(*itsD)[K_tilde]=(*itsUtilde)[CCC_tilde]
								-itsTheta1*ryy*0.5*( (*itsPrevQyy)[CCC]*U[GCC] - 2*(*itsPrevQyy)[CCC]*U[CCC] + (*itsPrevQyy)[CCC]*U[DCC])
								-itsTheta2*sy*0.5*( -(*itsPrevPy)[CCC]*U[GCC] + (*itsPrevPy)[CCC]*U[DCC]);
		
					}
				}

			}

			//Thomas
			tridiag(*itsGam,*itsA2,*itsB2,*itsC2,*itsD,*itsUtilde,lengthUtilde);

			//Troisieme et derniere etape du schema ADI
				
			//permutations
			e2=Nz+1;
			e3=(Nx+1)*(Nz+1);

			e_tilde2=Nz_tilde+1;
			e_tilde3=(Nx_tilde+1)*(Nz_tilde+1);

			//Construction du vecteur d
			if (Nz_tilde!=0)
			{
				for(K=0;K<lengthU;K++)
				{
					I=K%e3;
					i=I%e2;
					j=I/e2;
					k=K/e3;

					if( (i%Nz!=0)&&(j%Nx!=0)&&(k%Ny!=0) )//interrieur du domaine
					{
						i_tilde=i-1;
						j_tilde=j-1;
						I_tilde=e_tilde2*j_tilde+i_tilde;
						k_tilde=k-1;
						K_tilde=e_tilde3*k_tilde+I_tilde;

						//indices intervenants dans les differences finies
						CCC=i*e3_0+k*e2_0+j;
						GCC=CCC-e3_0;
						DCC=CCC+e3_0;
						CCC_tilde=i_tilde*e_tilde3_0+j_tilde*e_tilde2_1+k_tilde;

						//constantes du schéma
						rzz=2*deltaT/( pow( (*GetDeltaZ())[k-1],2 ) + pow( (*GetDeltaZ())[k],2 ) );
						sz= 2*deltaT/( (*GetDeltaZ())[k-1] + (*GetDeltaZ())[k] );

						if(i_tilde==0)
						{
							(*itsD)[K_tilde]=(*itsUtilde)[CCC_tilde]
									-itsTheta1*rzz*0.5*( (*itsPrevQzz)[CCC]*U[GCC] - 2*(*itsPrevQzz)[CCC]*U[CCC] +(*itsPrevQzz)[CCC]*U[DCC] )
									-itsTheta2*sz*0.5*( -(*itsPrevPz)[CCC]*U[GCC] + (*itsPrevPz)[CCC]*U[DCC] )
										+ itsTheta1*rzz*0.5*(*itsQzz)[CCC]*(U[CCC]-U[GCC])
										+ itsTheta2*sz*0.5*(*itsPz)[CCC]*(U[CCC]-U[GCC]);

						}
						else if(i_tilde==Nz_tilde)
						{
							(*itsD)[K_tilde]=(*itsUtilde)[CCC_tilde]
									-itsTheta1*rzz*0.5*( (*itsPrevQzz)[CCC]*U[GCC] - 2*(*itsPrevQzz)[CCC]*U[CCC] + (*itsPrevQzz)[CCC]*U[DCC] )
									-itsTheta2*sz*0.5*( -(*itsPrevPz)[CCC]*U[GCC] + (*itsPrevPz)[CCC]*U[DCC] )
										+ itsTheta1*rzz*0.5*(*itsQzz)[CCC]*(U[CCC]-U[DCC])
										- itsTheta2*sz*0.5*(*itsPz)[CCC]*(U[CCC]-U[DCC]);
						}
						else
						{
							(*itsD)[K_tilde]=(*itsUtilde)[CCC_tilde]
									-itsTheta1*rzz*0.5*( (*itsPrevQzz)[CCC]*U[GCC] - 2*(*itsPrevQzz)[CCC]*U[CCC] + (*itsPrevQzz)[CCC]*U[DCC] )
									-itsTheta2*sz*0.5*( -(*itsPrevPz)[CCC]*U[GCC] + (*itsPrevPz)[CCC]*U[DCC] );
						}
					}
				}
			}
			else
			{
				for(K=0;K<lengthU;K++)
				{
					I=K%e3;
					i=I%e2;
					j=I/e2;
					k=K/e3;

					if( (i%Nz!=0)&&(j%Nx!=0)&&(k%Ny!=0) )//interrieur du domaine
					{
						i_tilde=i-1;
						j_tilde=j-1;
						I_tilde=e_tilde2*j_tilde+i_tilde;
						k_tilde=k-1;
						K_tilde=e_tilde3*k_tilde+I_tilde;

						//indices intervenants dans les differences finies
						CCC=i*e3_0+k*e2_0+j;
						GCC=CCC-e3_0;
						DCC=CCC+e3_0;
						CCC_tilde=i_tilde*e_tilde3_0+j_tilde*e_tilde2_1+k_tilde;

						//constantes du schéma
						rzz=2*deltaT/( pow( (*GetDeltaZ())[k-1],2 ) + pow( (*GetDeltaZ())[k],2 ) );
						sz= 2*deltaT/( (*GetDeltaZ())[k-1] + (*GetDeltaZ())[k] );

						(*itsD)[K_tilde]=(*itsUtilde)[CCC_tilde]
									-itsTheta1*rzz*0.5*( (*itsPrevQzz)[CCC]*U[GCC] - 2*(*itsPrevQzz)[CCC]*U[CCC] + (*itsPrevQzz)[CCC]*U[DCC])
									-itsTheta2*sz*0.5*( -(*itsPrevPz)[CCC]*U[GCC] + (*itsPrevPz)[CCC]*U[DCC]);
		
					}
				}
			}
				
			//Thomas
			tridiag(*itsGam,*itsA3,*itsB3,*itsC3,*itsD,*itsUtilde,lengthUtilde);
	
			//permutations finales
			InvPermut(*itsUtilde,Nx_tilde,Ny_tilde,Nz_tilde,*itsUtilde_temp);

			//Construction du vecteur U_temp 
			Reconstruct(*itsUtilde_temp,U,*itsU_new);
		}
		break;
	}
	U = *itsU_new;
}


////////////////////////////////////////////////////
///	Class  : UpdatePayoffs2
///	Routine: UpdatePayoffs2
///	Returns: void
///	Action : Backward a payoff vector

			///OLD VERSION which was SLOW
////////////////////////////////////////////////////


void ARM_PDE3FCraigSneydNumericalScheme::UpdatePayoffs2(ARM_GP_Vector& U, int toTimeIdx)
{
	double deltaT=( (*GetTimeSteps())[toTimeIdx+1]-(*GetTimeSteps())[toTimeIdx] )/K_YEAR_LEN;

	size_t Nx=GetNX()-1;
	size_t Ny=GetNY()-1;
	size_t Nz=GetNZ()-1;
	
	size_t lengthU=(*itsQxx).size();
	
	//incréments directionnels
	
	size_t e1_xyz=1;
	size_t e2_xyz=Nx+1;
	size_t e3_xyz=(Nx+1)*(Ny+1);

	size_t e1_yxz=1;
	size_t e2_yxz=Ny+1;
	size_t e3_yxz=(Ny+1)*(Nx+1);

	size_t e1_zxy=1;
	size_t e2_zxy=Nz+1;
	size_t e3_zxy=(Nz+1)*(Nx+1);

	// Temporary variable for the permutation
	size_t K_xyz;
	size_t K_yxz, I_yxz;
	size_t i_yxz, j_yxz, k_yxz;
	size_t K_zxy, I_zxy;
	size_t i_zxy, j_zxy, k_zxy;

	if ( lengthU!=(Nx+1)*(Ny+1)*(Nz+1) ) 
	{
		throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"U de mauvaise taille");
	}

	//flag to avoid the recopy of the function
	//reference=Von Neumann
	double flagB1, flagB2;
	switch(itsBC)
	{
		case Dirichlet:
		{
			flagB1=2;
			flagB2=0;
		}
		break;

		case VonNeumann:
		{
			flagB1=1;
			flagB2=1;
		}
		break;

		case Belkheir:
		{
			flagB1=1;
			flagB2=1;
		}
		break;
	}

	deltaT=-deltaT;//pour tenir compte que l'on est Backward et non forward

	ARM_GP_T_Vector<bool>::iterator IsBoundit1 = itsIsBound.begin();
	ARM_GP_T_Vector<bool>::iterator IsBoundit2;
	ARM_GP_T_Vector<bool>::iterator IsBoundit3;
	std::vector<double>::iterator PrevQxxit = itsPrevQxx->begin();
	std::vector<double>::iterator PrevQyyit = itsPrevQyy->begin();
	std::vector<double>::iterator PrevQzzit = itsPrevQzz->begin();
	std::vector<double>::iterator PrevQxyit = itsPrevQxy->begin();
	std::vector<double>::iterator PrevQyzit = itsPrevQyz->begin();
	std::vector<double>::iterator PrevQzxit = itsPrevQzx->begin();
	std::vector<double>::iterator PrevPxit = itsPrevPx->begin();
	std::vector<double>::iterator PrevPyit = itsPrevPy->begin();
	std::vector<double>::iterator PrevPzit = itsPrevPz->begin();
	std::vector<double>::iterator PrevOit = itsPrevO->begin();
	std::vector<double>::iterator Rxxit = itsRxx->begin();
	std::vector<double>::iterator Ryyit = itsRyy->begin();
	std::vector<double>::iterator Rzzit= itsRzz->begin();
	std::vector<double>::iterator Rxyit = itsRxy->begin();
	std::vector<double>::iterator Ryzit = itsRyz->begin();
	std::vector<double>::iterator Rzxit= itsRzx->begin();
	std::vector<double>::iterator Sxit = itsSx->begin();
	std::vector<double>::iterator Syit = itsSy->begin();
	std::vector<double>::iterator Szit= itsSz->begin();

	std::vector<double>::iterator Dit= itsD->begin();


	std::vector<double>::iterator UCCCit=U.begin();
	std::vector<double>::iterator U_tempCCCit;
	std::vector<double>::iterator UDCCit;
	std::vector<double>::iterator UGCCit;
	std::vector<double>::iterator UCDCit;
	std::vector<double>::iterator UCGCit;
	
	std::vector<double>::iterator UGDCit;
	std::vector<double>::iterator UGGCit;
	std::vector<double>::iterator UDDCit;
	std::vector<double>::iterator UDGCit;

	std::vector<double>::iterator UCCGit;
	std::vector<double>::iterator UDCGit;
	std::vector<double>::iterator UGCGit;
	std::vector<double>::iterator UCDGit;
	std::vector<double>::iterator UCGGit;

	std::vector<double>::iterator UCCDit;
	std::vector<double>::iterator UDCDit;
	std::vector<double>::iterator UGCDit;
	std::vector<double>::iterator UCDDit;
	std::vector<double>::iterator UCGDit;
	

	size_t dec = 0;
	size_t dec_temp = 0;

	for(K_xyz=0;K_xyz<lengthU;K_xyz++)
	{
		if( !(*IsBoundit1) )//interrieur du domaine
		{
			//indices qui vont intervenir dans les differences finies
			UCCCit=U.begin()+K_xyz;
			UDCCit=UCCCit+e1_xyz;
			UGCCit=UCCCit-e1_xyz;
			UCDCit=UCCCit+e2_xyz;
			UCGCit=UCCCit-e2_xyz;
				
			UGDCit=UGCCit+e2_xyz;
			UGGCit=UGCCit-e2_xyz;
			UDDCit=UDCCit+e2_xyz;
			UDGCit=UDCCit-e2_xyz;

			UCCGit=UCCCit-e3_xyz;
			UDCGit=UCCGit+e1_xyz;
			UGCGit=UCCGit-e1_xyz;
			UCDGit=UCCGit+e2_xyz;
			UCGGit=UCCGit-e2_xyz;

			UCCDit=UCCCit+e3_xyz;
			UDCDit=UCCDit+e1_xyz;
			UGCDit=UCCDit-e1_xyz;
			UCDDit=UCCDit+e2_xyz;
			UCGDit=UCCDit-e2_xyz;
			
			(*Dit)=(*UCCCit)
					+(1-itsTheta1)*(*Rxxit)*0.5*(*PrevQxxit)*( (*UGCCit) -2*(*UCCCit) + (*UDCCit) )
					+(*Ryyit)*0.5*(*PrevQyyit)*( (*UCGCit) - 2*(*UCCCit) + (*UCDCit) )
					+(*Rzzit)*0.5*(*PrevQzzit)*( (*UCCGit) - 2*(*UCCCit) + (*UCCDit) )
					+(*Rxyit)*0.25*(*PrevQxyit)*( (*UGGCit) - (*UDGCit) - (*UGDCit) + (*UDDCit) )
					+(*Ryzit)*0.25*(*PrevQyzit)*( (*UCGGit) - (*UCDGit) - (*UCGDit) + (*UCDDit) )
					+(*Rzxit)*0.25*(*PrevQzxit)*( (*UGCGit) - (*UDCGit) - (*UGCDit) + (*UDCDit) )
					+(1-itsTheta2)*(*Sxit)*0.5*(*PrevPxit)*( -(*UGCCit) + (*UDCCit) )
					+(*Syit)*0.5*(*PrevPyit)*( -(*UCGCit) + (*UCDCit) )
					+(*Szit)*0.5*(*PrevPzit)*( -(*UCCGit) + (*UCCDit) )
					+(1-itsTheta3)*deltaT*(*PrevOit)*(*UCCCit);
		}
		else
		{
			(*Dit)=0.0;
		}

		// Base de lecture xyz
		IsBoundit1++;
		PrevQxxit++;
		PrevQyyit++;
		PrevQzzit++;
		PrevPxit++;
		PrevPyit++;
		PrevPzit++;
		PrevOit++;
		Rxxit++;
		Ryyit++;
		Rzzit++;
		Rxyit++;
		Ryzit++;
		Rzxit++;
		Sxit++;
		Syit++;
		Szit++;

		// Base d'ecriture xyz
		Dit++;
	}

	tridiag2(*itsGam,*itsA1,*itsB1,*itsC1,*itsD,*itsU_temp,lengthU);

				
	//Deuxieme etape du schema ADI
	
	//Permutations xyz -> yxz

	Dit = itsD->begin();

	for(K_yxz=0;K_yxz<lengthU;K_yxz++)
	{
		k_yxz = K_yxz/e3_yxz;
		I_yxz = K_yxz%e3_yxz;
		j_yxz = I_yxz/e2_yxz;
		i_yxz = I_yxz%e2_yxz;

		dec = k_yxz*e3_xyz+i_yxz*e2_xyz+j_yxz*e1_xyz;

		IsBoundit2 = itsIsBound.begin() + dec;

		if( !(*IsBoundit2) )//interrieur du domaine
		{
			//indices qui vont intervenir dans les differences finies
			
			// Base de lecture xyz
			U_tempCCCit=itsU_temp->begin()+dec;
			UCCCit=U.begin()+dec;
			UDCCit=UCCCit+e2_xyz;
			UGCCit=UCCCit-e2_xyz;
			PrevQyyit = itsPrevQyy->begin()+dec;
			PrevPyit = itsPrevPy->begin()+dec;
			Ryyit = itsRyy->begin()+dec;
			Syit = itsSy->begin()+dec;
		
			(*Dit)=(*U_tempCCCit)
					-itsTheta1*(*Ryyit)*0.5*(*PrevQyyit)*( (*UGCCit) - 2*(*UCCCit) + (*UDCCit) )
					-itsTheta2*(*Syit)*0.5*(*PrevPyit)*( -(*UGCCit) + (*UDCCit) );
		}
		else
		{
			(*Dit)=0.0;
		}
		// Base d'ecriture yxz
		Dit++;
	}

	tridiag2(*itsGam,*itsA2,*itsB2,*itsC2,*itsD,*itsU_temp,lengthU);

	//Troisieme etape du schema ADI
	
	//Permutations yxz -> zxy

	Dit = itsD->begin();

	for(K_zxy=0;K_zxy<lengthU;K_zxy++)
	{
		k_zxy = K_zxy/e3_zxy;
		I_zxy = K_zxy%e3_zxy;
		j_zxy = I_zxy/e2_zxy;
		i_zxy = I_zxy%e2_zxy;


		// Base de lecture xyz
		dec = i_zxy*e3_xyz+k_zxy*e2_xyz+j_zxy*e1_xyz;

		IsBoundit3 = itsIsBound.begin() + dec;

		if( !(*IsBoundit3) )//interrieur du domaine
		{
			//indices qui vont intervenir dans les differences finies
			
			// Base de lecture yxz
			// A vérifier !!!!!!
			dec_temp = i_zxy*e3_yxz + j_zxy*e2_yxz + k_zxy*e1_yxz;

			U_tempCCCit=itsU_temp->begin()+dec_temp;

			// Base de lecture xyz
			UCCCit=U.begin()+dec;
			UDCCit=UCCCit+e3_xyz;
			UGCCit=UCCCit-e3_xyz;
			PrevQzzit = itsPrevQzz->begin()+dec;
			PrevPzit = itsPrevPz->begin()+dec;
			Rzzit = itsRzz->begin()+dec;
			Szit = itsSz->begin()+dec;
		
			(*Dit)=(*U_tempCCCit)
					-itsTheta1*(*Rzzit)*0.5*(*PrevQzzit)*( (*UGCCit) - 2*(*UCCCit) + (*UDCCit) )
					-itsTheta2*(*Szit)*0.5*(*PrevPzit)*( -(*UGCCit) + (*UDCCit) );
		}
		else
		{
			(*Dit)=0.0;
		}
		// Base d'ecriture zxy
		Dit++;
	}

	tridiag2(*itsGam,*itsA3,*itsB3,*itsC3,*itsD,*itsU_temp,lengthU);

	//permutations finales
	InvPermut2(*itsU_temp,Nx,Ny,Nz,*itsU_new);

	//Construction du vecteur U_new 
	Reconstruct2(*itsU_new);
	
	U = *itsU_new;
}

////////////////////////////////////////////////////
///	Class  : ARM_PDE3FCraigSneydNumericalScheme 
///	Routine: UpdateIntermediatePayoffs
///	Returns: void
///	Action : Backward a collection of vectors stocked
///    in a matrix M

			///OLD VERSION which was SLOW
////////////////////////////////////////////////////

/*
void ARM_PDE3FCraigSneydNumericalScheme::UpdateIntermediatePayoffs(ARM_GP_Matrix& vec, int toTimeIdx)
{
	std::vector<double> vec2(vec.cols());

	for( size_t j=0 ; j<vec.rows() ; ++j )
	{
		for( size_t i=0 ; i<vec.cols() ; ++i)
			vec2[i] = vec(j,i);

		UpdatePayoffs( vec2, toTimeIdx );
	
		for( i=0 ; i<vec.cols() ; ++i)
			vec(j,i) = vec2[i];
	}
}
*/

////////////////////////////////////////////////////
///	Class  : ARM_PDE3FCraigSneydNumericalScheme 
///	Routine: UpdateIntermediatePayoffs2
///	Returns: void
///	Action : Backward a collection of vectors stocked
///    in a matrix M
////////////////////////////////////////////////////

void ARM_PDE3FCraigSneydNumericalScheme::UpdateIntermediatePayoffs2(ARM_PayoffStates& payoffStates,  int toTimeIdx)
{
	ARM_GP_Matrix& vec = const_cast<ARM_GP_Matrix&> (payoffStates.GetIntermediatePayoffs() );
//	ARM_GP_VectorPtr vec2(new std::vector<double>(vec.cols()));
	std::vector<double>::iterator vec2it;

	for( size_t j=0 ; j<vec.rows() ; ++j )
	{
		ARM_GP_VectorPtr vec2(new ARM_GP_Vector(vec.cols()));
		vec2it = vec2->begin();
		for( size_t i=0 ; i<vec.cols() ; ++i)
		{
			*vec2it = vec(j,i);
			vec2it++;
		}

		UpdatePayoffs2( *vec2, toTimeIdx );
		payoffStates.SetPayoffSnapshot(j,vec2);

		
		vec2it = vec2->begin();
		for( size_t i=0 ; i<vec.cols() ; ++i)
		{
			vec(j,i) = *vec2it;
			vec2it++;
		}
	}
}

////////////////////////////////////////////////////
///	Class  : ARM_PDE3FCraigSneydNumericalScheme 
///	Routine: toString
///	Returns: string
///	Action : As usual for a toString method
////////////////////////////////////////////////////
string ARM_PDE3FCraigSneydNumericalScheme::toString(const string& indent, const string& nextIndent) const
{
	CC_Ostringstream os;
	os << indent << " 3F Craig Sneyd Numerical Scheme " << CC_NS(std,endl);
	os << indent << "     XGrid Discretization : " << GetNX() << " Points" << CC_NS(std,endl);
	os << indent << "     YGrid Discretization : " << GetNY() << " Points" << CC_NS(std,endl);
	os << indent << "     ZGrid Discretization : " << GetNZ() << " Points" << CC_NS(std,endl);
	os << indent << "     Operture in StdDev mode : " << getStdDevNb() << " Standard Deviations" << CC_NS(std,endl);
	if( (GetXGrid()->size()>0)&&(GetYGrid()->size()>0)&&(GetZGrid()->size()>0) )
	{
		os << indent << "     Xmax  : " << (*GetXGrid())[GetXGrid()->size()-1] << CC_NS(std,endl);
		os << indent << "     Ymax  : " << (*GetYGrid())[GetYGrid()->size()-1] << CC_NS(std,endl);
		os << indent << "     Zmax  : " << (*GetZGrid())[GetZGrid()->size()-1] << CC_NS(std,endl);
	}
	else
	{
		os << indent << "     Xmax  not initilalized " << CC_NS(std,endl);
		os << indent << "     Ymax  not initilalized " << CC_NS(std,endl);
		os << indent << "     Zmax  not initilalized " << CC_NS(std,endl);
	}
	if( (GetDeltaX()->size()>0)&&(GetDeltaY()->size()>0)&&(GetDeltaZ()->size()>0) )
	{
		os << indent << "     DeltaX int  : " << (*GetDeltaX())[(GetDeltaX()->size()-1)/2] << "     DeltaX ext  : " << (*GetDeltaX())[0] << CC_NS(std,endl);
		os << indent << "     DeltaY int  : " << (*GetDeltaY())[(GetDeltaY()->size()-1)/2] << "     DeltaY ext  : " << (*GetDeltaY())[0] << CC_NS(std,endl);
		os << indent << "     DeltaZ int  : " << (*GetDeltaZ())[(GetDeltaZ()->size()-1)/2] << "     DeltaZ ext  : " << (*GetDeltaZ())[0] << CC_NS(std,endl);
	}
	else
	{
		os << indent << "     DeltaX  not initilalized " << CC_NS(std,endl);
		os << indent << "     DeltaY  not initilalized " << CC_NS(std,endl);
		os << indent << "     DeltaZ  not initilalized " << CC_NS(std,endl);
	}
	os << indent << "     Theta1	       : " << itsTheta1 << CC_NS(std,endl);
	os << indent << "     Theta2	       : " << itsTheta2 << CC_NS(std,endl);
	os << indent << "     Theta3	       : " << itsTheta3 << CC_NS(std,endl);
	os << indent << "     Boundary Condition : " << ARM_ArgConvReverse_PDEBoundConditionType.GetString(itsBC) << CC_NS(std,endl);
	os << indent << "     Lambda : " << itsLambda << CC_NS(std,endl);
	return os.str();
}

CC_END_NAMESPACE()


