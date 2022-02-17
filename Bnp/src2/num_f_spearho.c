
// S. Galluccio: 21 November 2000
// Created Source File

#include "math.h"
#include "num_h_allhdr.h"
#define NRANSI

//  Computes the Spearman rho correlation coeeficient given two sets of data,
// that is data1 and data2 of same length n 

double spear_rho(	double *data1,  // input
					double *data2,  // input 
					unsigned long n)			//input
{

	unsigned long i,*indx1,*indx2,*irank1,*irank2;
	double avg1=0.0,avg2=0.0,var1=0.0,var2=0.0,rhocorr=0.0;

	indx1=lvector(1,n);
	indx2=lvector(1,n);
	irank1=lvector(1,n);
	irank2=lvector(1,n);

//////////////// ranks vector1 ////////////////////////////////////////

	index_data(	n, 
				data1, 
				indx1); 
	rank_data(	n, 
				indx1, // input: the indexed table 
				irank1);  // output: the rank table 

//////////////// ranks vector2 ////////////////////////////////////////

	index_data(	n, 
				data2, 
				indx2);	
	rank_data(	n, 
				indx2, // input: the indexed table 
				irank2);  // output: the rank table 

////////////////////////////////// compute the rank correlation matrix
		
			for (i=1;i<=n;i++) {
				avg1+=(irank1[i]-avg1)/i;
				var1+=(irank1[i]*irank1[i]-var1)/i;

				avg2+=(irank2[i]-avg2)/i;
				var2+=(irank2[i]*irank2[i]-var2)/i;

				rhocorr+=(irank1[i]*irank2[i]-rhocorr)/i;					
			}			
		
			rhocorr=(rhocorr-avg1*avg2)/sqrt(var1-avg1*avg1)/sqrt(var2-avg2*avg2);

	//////////////// free the memory //////////////////////////////////

	free_lvector(indx1,1,n);
	free_lvector(indx2,1,n);
	free_lvector(irank1,1,n);
	free_lvector(irank2,1,n);

	return rhocorr;
}
#undef NRANSI

///////////////////////////////////////////////////////////////////////

void rank_vector(	int n, 
					double *w, 
					double *s)
{
	int j=1,ji,jt;
	double t,rank;

	*s=0.0;
	while (j < n) {
		if (w[j+1] != w[j]) {
			w[j]=j;
			++j;
		} else {
			for (jt=j+1;jt<=n && w[jt]==w[j];jt++);
			rank=0.5*(j+jt-1);
			for (ji=j;ji<=(jt-1);ji++) w[ji]=rank;
			t=jt-j;
			*s += t*t*t-t;
			j=jt;
		}
	}
	if (j == n) w[n]=n;
}

