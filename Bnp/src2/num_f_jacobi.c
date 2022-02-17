/* ========================================================================

   MODULE:	    num_f_jacobi.c

   LIBRARY:	    UTL_LIB

   FUNCTION:	jacobi_diagonalisation

   DESCRIPTION:	From Numerical Recipes to diagonalise a symmetrix matrix
   From NUMERICAL RECIPES IN C
   
   Computes all eigenvalues and eigenvectors of a real symmetric matrix
   sym_mat[0..n-1][0..n-1]. 
   On output, elements of matrix above the diagonal are destroyed !!!!! 
   The vector eigen_val[0..n-1] returns the eigenvalues of sym_mata, 
   sorted by decreasing order.  
   The matrix eigen_vec[0..n-1][0..n-1] is a matrix whose columns
   contain, on output, the normalized eigenvectors of sym_mat.  

  ========================================================================== */

#include "utallhdr.h" 
#include "math.h"           


#define NRANSI
#define ROTATE(a,i,j,k,l) g=a[i][j];h=a[k][l];a[i][j]=g-s*(h+g*tau);\
	a[k][l]=h+s*(g-h*tau);


static Err sort_eigen_values(double eigen_val[], double **eigen_vec, int n)	
/* From NUMERICAL RECIPES IN C
   Given the dimension of the matrix and its eigen_values,
   sorts the eigen values in the eigen_val vector, 
   and replaces the eigen vectors in the right order in **eigen_vec*/
{
	int i,j,k;
	double p;
	Err err=NULL;

	for (i = 0; i<n-1;i++)
	{
		p=eigen_val[k=i];
		for (j = i+1;j<n;j++)
			if (eigen_val[j] >=p) 
				p=eigen_val[k=j];
		if (k !=i)
		{
			eigen_val[k] = eigen_val[i];
			eigen_val[i] = p;
 			for (j=0;j<n;j++)
			{
				p=eigen_vec[j][i];
				eigen_vec[j][i]=eigen_vec[j][k];
				eigen_vec[j][k]=p;
			}
		}
	}
	return err;


}

static Err sort_eigen_values2(double *eigen_val, double **eigen_vec, int n)	
/* From NUMERICAL RECIPES IN C
   Given the dimension of the matrix and its eigen_values,
   sorts the eigen values in the eigen_val vector, 
   and replaces the eigen vectors in the right order in **eigen_vec*/
{
	int i,j,k;
	double p;
	Err err=NULL;

	for (i = 0; i<n-1;i++)
	{
		p=eigen_val[k=i];
		for (j = i+1;j<n;j++)
			if (eigen_val[j] >=p) 
				p=eigen_val[k=j];
		if (k !=i)
		{
			eigen_val[k] = eigen_val[i];
			eigen_val[i] = p;
 			for (j=0;j<n;j++)
			{
				p=eigen_vec[j][i];
				eigen_vec[j][i]=eigen_vec[j][k];
				eigen_vec[j][k]=p;
			}
		}
	}
	return err;


}

/* ---------------------------------------------------------------------------------- */
/* ---------------------------------------------------------------------------------- 
   Computes the eigenvalues  and eigenvectors of a real symmetric matrix 
   sym_mat[0..n-1][0..n-1] On output, elements of sym_mat above the diagonal are 
   destroyed. eigen_val[0..n-1] returns the eigen values of sym_mat. 
   eigen_vec[0..n-1][0..n-1] is a matrix that stores all the corresponding renormalised 
   eigen vectors of sym_mat  ( eigen_vec[...][i] corresponds to eigen_val[i] )
   nrot returns the number of Jacobi rotations that were required
   ---------------------------------------------------------------------------------- */
Err jacobi_diagonalisation(
				double  **sym_mat, 
				int     n, 
				double  eigen_val[], 
				double  **eigen_vec, 
				int     *nrot)	

{
int     i,iq,ip,j;
double  tresh,theta,tau,t,sum,s,h,g,c,*b,*z;
Err     err = NULL;

	b=dvector(0,n-1);
	z=dvector(0,n-1);
	for (ip=0;ip<n;ip++) 
	{
		for (iq=0;iq<n;iq++) 
			eigen_vec[ip][iq]=0.0;
		eigen_vec[ip][ip]=1.0;
	}
	for (ip=0;ip<n;ip++) 
	{
		b[ip]=sym_mat[ip][ip];
		eigen_val[ip]=sym_mat[ip][ip];
		z[ip]=0.0;
	}
	*nrot = 0;
	for (i=1;i<=50;i++) 
	{
		sum=0.0;
		for (ip=0;ip<n-1;ip++) 
		{
			for (iq=ip+1;iq<n;iq++)
				sum += fabs(sym_mat[ip][iq]);
		}
		if (sum == 0.0) 
		{
			free_dvector(z,0,n-1);
			free_dvector(b,0,n-1);
			err = sort_eigen_values(eigen_val, eigen_vec, n); 
			return err;
		}
		if (i < 4)        /* On the first three sweeps */
			tresh=0.2*sum/(n*n);
		else
			tresh=0.0;
		for (ip=0;ip<n-1;ip++) 
		{
			for (iq=ip+1;iq<n;iq++) 
			{
				g=100.0*fabs(sym_mat[ip][iq]);
				if ( (i > 4) && 
				(double)(fabs(eigen_val[ip])+g) 
					== (double)fabs(eigen_val[ip])
				&& (double)(fabs(eigen_val[iq])+g) 
					== (double)fabs(eigen_val[iq]))
					sym_mat[ip][iq]=0.0;
				else 
				if (fabs(sym_mat[ip][iq]) > tresh) 
				{
					h=eigen_val[iq]
						-eigen_val[ip];
					if ((double)(fabs(h)+g) 
							== (double)fabs(h))
						t=(sym_mat[ip][iq])/h;
					else 
					{
						theta=0.5*h
							/(sym_mat[ip][iq]);
						t= 1.0/ (fabs(theta)
							+sqrt(1.0+theta*theta));
						if (theta < 0.0) 
								t = -t;
					}
					c=1.0/sqrt(1+t*t);
					s=t*c;
					tau=s/(1.0+c);
					h=t*sym_mat[ip][iq];
					z[ip] -= h;
					z[iq] += h;
					eigen_val[ip] -= h;
					eigen_val[iq] += h;
					sym_mat[ip][iq]=0.0;
					for (j=0;j<=ip-1;j++) 
					{
						ROTATE(sym_mat,j,ip,j,iq);
					}
					for (j=ip+1;j<=iq-1;j++) 
					{
						ROTATE(sym_mat,ip,j,j,iq);
					}
					for (j=iq+1;j<n;j++)
					{
						ROTATE(sym_mat,ip,j,iq,j);
					}
					for (j=0;j<n;j++)
					{
						ROTATE(eigen_vec,
							j,ip,j,iq);
					}
				}
			}
		}
		for (ip=0;ip<n;ip++) 
		{ 
			b[ip] += z[ip];
			eigen_val[ip]=b[ip];
			z[ip]=0.0;
		}
	}
	sort_eigen_values(eigen_val, eigen_vec, n); 
	return err;
}


Err jacobi_diagonalisation2(
				double  **sym_mat, 
				int     n, 
				double  *eigen_val, 
				double  **eigen_vec, 
				int     *nrot)	

{
int     i,iq,ip,j;
double  tresh,theta,tau,t,sum,s,h,g,c,*b,*z;
Err     err = NULL;

	b=dvector(0,n-1);
	z=dvector(0,n-1);
	for (ip=0;ip<n;ip++) 
	{
		for (iq=0;iq<n;iq++) 
			eigen_vec[ip][iq]=0.0;
		eigen_vec[ip][ip]=1.0;
	}
	for (ip=0;ip<n;ip++) 
	{
		b[ip]=sym_mat[ip][ip];
		eigen_val[ip]=sym_mat[ip][ip];
		z[ip]=0.0;
	}
	*nrot = 0;
	for (i=1;i<=50;i++) 
	{
		sum=0.0;
		for (ip=0;ip<n-1;ip++) 
		{
			for (iq=ip+1;iq<n;iq++)
				sum += fabs(sym_mat[ip][iq]);
		}
		if (sum == 0.0) 
		{
			free_dvector(z,0,n-1);
			free_dvector(b,0,n-1);
			err = sort_eigen_values2(eigen_val, eigen_vec, n); 
			return err;
		}
		if (i < 4)        /* On the first three sweeps */
			tresh=0.2*sum/(n*n);
		else
			tresh=0.0;
		for (ip=0;ip<n-1;ip++) 
		{
			for (iq=ip+1;iq<n;iq++) 
			{
				g=100.0*fabs(sym_mat[ip][iq]);
				if ( (i > 4) && 
				(double)(fabs(eigen_val[ip])+g) 
					== (double)fabs(eigen_val[ip])
				&& (double)(fabs(eigen_val[iq])+g) 
					== (double)fabs(eigen_val[iq]))
					sym_mat[ip][iq]=0.0;
				else 
				if (fabs(sym_mat[ip][iq]) > tresh) 
				{
					h=eigen_val[iq]
						-eigen_val[ip];
					if ((double)(fabs(h)+g) 
							== (double)fabs(h))
						t=(sym_mat[ip][iq])/h;
					else 
					{
						theta=0.5*h
							/(sym_mat[ip][iq]);
						t= 1.0/ (fabs(theta)
							+sqrt(1.0+theta*theta));
						if (theta < 0.0) 
								t = -t;
					}
					c=1.0/sqrt(1+t*t);
					s=t*c;
					tau=s/(1.0+c);
					h=t*sym_mat[ip][iq];
					z[ip] -= h;
					z[iq] += h;
					eigen_val[ip] -= h;
					eigen_val[iq] += h;
					sym_mat[ip][iq]=0.0;
					for (j=0;j<=ip-1;j++) 
					{
						ROTATE(sym_mat,j,ip,j,iq);
					}
					for (j=ip+1;j<=iq-1;j++) 
					{
						ROTATE(sym_mat,ip,j,j,iq);
					}
					for (j=iq+1;j<n;j++)
					{
						ROTATE(sym_mat,ip,j,iq,j);
					}
					for (j=0;j<n;j++)
					{
						ROTATE(eigen_vec,
							j,ip,j,iq);
					}
				}
			}
		}
		for (ip=0;ip<n;ip++) 
		{ 
			b[ip] += z[ip];
			eigen_val[ip]=b[ip];
			z[ip]=0.0;
		}
	}
	sort_eigen_values2(eigen_val, eigen_vec, n); 
	return err;
}

#undef ROTATE
#undef NRANSI
