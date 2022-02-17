#include "edginc/config.hpp"
#include "edginc/regression.h"
#include "edginc/Nrfns.hpp"
#include <iostream>
#include <numeric>
#include <vector>
#include <math.h>

DRLIB_BEGIN_NAMESPACE
//#include "alib/macros.h"

//some math functions                                        
static void NRfit(const vector<vector<double> > & x, const vector<double> & y, vector<double> & coeffs, double & chisq, vector<double> & errors);
static void NRsvdfit(const vector<vector<double> > & x, const vector<double> & y, vector<double> & coeffs, double & chisq, vector<double> & errors);
void svbksb(vector<vector<double> > & u, vector<double> & w, vector<vector<double> > & v, int m, int n, vector<double> & b, vector<double> & x);
void svdcmp(vector<vector<double> > & a, int m, int n, vector<double> & w, vector<vector<double> > & v);
static double pythag(double a, double b);
void gaussj(vector<vector<double> > & a, int n, vector<vector<double> > & b, int m);

void PlainRegression::DoTheRegression(vector<vector<double> > & basisValues, vector<double> & dealValues, vector<double> & coeffs)
{
  vector<double> errors;
  double chisq;
  NRfit( basisValues, dealValues, coeffs, chisq, errors);

  for(size_t i=0; i < coeffs.size() ; ++i)
	  clog << "PlainRegression: coeff[" << i << "]= " << coeffs[i] << endl;
  clog << "PlainRegression: chisq= " << chisq << endl;
}

void SVDRegression::DoTheRegression(vector<vector<double> > & basisValues, vector<double> & dealValues, vector<double> & coeffs)
{
  vector<double> errors;
  double chisq;
  NRsvdfit( basisValues, dealValues, coeffs, chisq, errors);
  clog << "SVDRegression: number of samples= " << basisValues.size() << endl;
  for(size_t i=0; i < basisValues.size(); ++i)
  {
	  clog << "[" << i << "] size= " << basisValues[i].size();
	  for(size_t j = 0; j < basisValues[i].size(); ++j)
		  clog << " " << basisValues[i][j];
	  clog << " => " << dealValues[i] << endl;
  }
  clog << "Results of regression: " << endl;
  for(size_t i=0; i < coeffs.size() ; ++i)
	  clog << "SVDRegression: coeff[" << i << "]= " << coeffs[i] << endl;
  clog << "SVDRegression: chisq= " << chisq << endl;
}

static void NRsvdfit(const vector<vector<double> > & f, const vector<double> & y, vector<double> & coeffs, double & chi, vector<double> & errors)
/* Given a set of data points x[1..ndata],y[1..ndata] with individual standard deviations
sig[1..ndata], use chi2 minimization to determine the coe.cients a[1..ma] of the fitting
function y = sum(afunci(x)). Here we solve the fitting equations using singular
value decomposition of the ndata by ma matrix. Arrays u[1..ndata][1..ma],
v[1..ma][1..ma], and w[1..ma] provide workspace on input; on output they define the
singular value decomposition, and can be used to obtain the covariance matrix. The program
returns values for the ma fit parameters a and chi */
{
	int i, j;
	double wmax,thresh,sum;

	int ndat= y.size();
	if (ndat==0) return;
	int ma = f[0].size();	// ma = nr of basis functions

	vector<double> b(ndat);
	vector<double> afunc(ma);
	vector<double> a(ma);
	vector<vector<double> > u(ndat, vector<double> (ma));
	vector<vector<double> > v(ma, vector<double> (ma));
	vector<double> w(ma);

	for (i=0; i<ndat; i++) {
		afunc = f[i];
		for (j=0; j<ma; j++)	u[i][j] = afunc[j];
		b[i] = y[i];
	}

	svdcmp(u,ndat,ma,w,v);		// Singular value decomposition.
	wmax = 0.;
	for (j=0; j<ma; j++)
		if (w[j] > wmax) wmax=w[j];
	thresh = TOL*wmax;
	for (j=0; j<ma; j++)
		if (w[j] < thresh) w[j] = 0.;
	svbksb(u,w,v,ndat,ma,b,a);

	coeffs.clear();
	for (int l=0; l<ma; l++) coeffs.push_back(a[l]);	// This line must be added!!!
	double chisq = 0;
	errors.clear();
	for (i=0; i<ndat; i++) {
		afunc = f[i];
		for (sum=0.,j=0; j<ma; j++) sum += coeffs[j]*afunc[j];
		chisq += SQR(y[i]-sum);
		errors.push_back(y[i]-sum);
	}
	chi = sqrt(chisq/ndat);
}

void svbksb(vector<vector<double> > & u, vector<double> & w, vector<vector<double> > & v, int m, int n, vector<double> & b, vector<double> & x)
/* Solves A·X = B for a vector X, where A is specified by the arrays u[1..m][1..n], w[1..n],
v[1..n][1..n] as returned by svdcmp. m and n are the dimensions of a, and will be equal for
square matrices. b[1..m] is the input right-hand side. x[1..n] is the output solution vector.
No input quantities are destroyed, so the routine may be called sequentially with different b’s */
{
	int jj,j,i;
	double s;
	
	vector<double> tmp(n);

	for (j=0; j<n; j++) {
		s = 0.;
		if (w[j]) {
			for (i=0; i<m; i++) s += u[i][j]*b[i];
			s /= w[j];
		}
		tmp[j] = s;
	}
	for (j=0; j<n; j++) {
		s = 0.;
		for (jj=0; jj<n; jj++) s += v[j][jj]*tmp[jj];
		x[j] = s;
	}
}

void svdcmp(vector<vector<double> > & a, int m, int n, vector<double> & w, vector<vector<double> > & v)
/* Given a matrix a[1..m][1..n], this routine computes its singular value decomposition,
A =U·W·V^T. The matrix U replaces a on output. The diagonal matrix of singular values W is output
as a vector w[1..n]. The matrix V (not the transpose V^T ) is output as v[1..n][1..n] */
{
//	double pythag(double a, double b);
	int flag,i,its,j,jj,k,l,nm;
	double anorm,c,f,g,h,s,scale,x,y,z;

	vector<double> rv1(n);
	g = scale = anorm=0.0;
	for (i=0; i<n; i++) {
		l = i+1;
		rv1[i] = scale*g;
		g = s = scale = 0.;
		if (i < m) {				//			??? i<= m or i<m
			for (k=i; k<m; k++) scale += fabs(a[k][i]);
			if (scale) {
				for (k=i; k<m; k++) {
					a[k][i] /= scale;
					s += a[k][i]*a[k][i];
				}
				f = a[i][i];
				g = -SIGN(sqrt(s),f);
				h = f*g-s;
				a[i][i] = f-g;
				for (j=l; j<n; j++) {
					for (s=0.,k=i; k<m; k++) s += a[k][i]*a[k][j];
					f = s/h;
					for (k=i; k<m; k++) a[k][j] += f*a[k][i];
				}
				for (k=i; k<m; k++) a[k][i] *= scale;
			}
		}
		w[i] = scale*g;
		g = s = scale = 0.;
		if (i < m && i != n-1) {					// ???
			for (k=l; k<n; k++) scale += fabs(a[i][k]);
			if (scale) {
				for (k=l; k<n; k++) {
					a[i][k] /= scale;
					s += a[i][k]*a[i][k];
				}
				f = a[i][l];
				g = -SIGN(sqrt(s),f);
				h = f*g-s;
				a[i][l] = f-g;
				for (k=l; k<n; k++) rv1[k] = a[i][k]/h;
				for (j=l; j<m; j++) {
					for (s=0.,k=l; k<n; k++) s += a[j][k]*a[i][k];
					for (k=l; k<n; k++) a[j][k] += s*rv1[k];
				}
				for (k=l; k<n; k++) a[i][k] *= scale;
			}
		}
		anorm = MAX(anorm,(fabs(w[i])+fabs(rv1[i])));
	}
	for (i=n-1; i>=0; i--) {
		if (i < n-1) {
			if (g) {
				for (j=l; j<n; j++)
					v[j][i] = (a[i][j]/a[i][l])/g;
				for (j=l; j<n; j++) {
					for (s=0.,k=l; k<n; k++) s += a[i][k]*v[k][j];
					for (k=l; k<n; k++) v[k][j] += s*v[k][i];
				}
			}
			for (j=l; j<n; j++) v[i][j] = v[j][i] = 0.;
		}
		v[i][i] = 1.;
		g = rv1[i];
		l = i;
	}
	for (i=MIN(m,n)-1; i>=0; i--) {
		l = i+1;
		g = w[i];
		for (j=l; j<n; j++) a[i][j] = 0.;
		if (g) {
			g = 1./g;
			for (j=l; j<n; j++) {
				for (s=0.,k=l; k<m; k++) s += a[k][i]*a[k][j];
				f = (s/a[i][i])*g;
				for (k=i; k<m; k++) a[k][j] += f*a[k][i];
			}
			for (j=i; j<m; j++) a[j][i] *= g;
		}
		else
			for (j=i; j<m; j++) a[j][i] = 0.;
		++a[i][i];
	}
	for (k=n-1; k>=0; k--) {
		for (its=1; its<=30; its++) {
			flag = 1;
			for (l=k; l>=0; l--) {
				nm = l-1;
				if ((double)(fabs(rv1[l])+anorm) == anorm) {
					flag=0;
					break;
				}
				if ((double)(fabs(w[nm])+anorm) == anorm) break;
			}
			if (flag) {
				c = 0.;
				s = 1.;
				for (i=l; i<=k; i++) {
					f = s*rv1[i];
					rv1[i] = c*rv1[i];
					if ((double)(fabs(f)+anorm) == anorm) break;
					g = w[i];
					h = pythag(f,g);
					w[i] = h;
					h = 1./h;
					c = g*h;
					s = -f*h;
					for (j=0; j<m; j++) {
						y = a[j][nm];
						z = a[j][i];
						a[j][nm] = y*c+z*s;
						a[j][i] = z*c-y*s;
					}
				}
			}
			z = w[k];
			if (l == k) {
				if (z < 0.) {
					w[k] = -z;
					for (j=0; j<n; j++) v[j][k] = -v[j][k];
				}
				break;
			}
			if (its == 100) cout<<"no convergence in 30 svdcmp iterations"<<endl;
			x = w[l];
			nm = k-1;
			y = w[nm];
			g = rv1[nm];
			h = rv1[k];
			f = ((y-z)*(y+z)+(g-h)*(g+h))/(2.*h*y);
			g = pythag(f,1.);
			f = ((x-z)*(x+z)+h*((y/(f+SIGN(g,f)))-h))/x;
			c = s = 1.;
			for (j=l; j<=nm; j++) {
				i = j+1;
				g = rv1[i];
				y = w[i];
				h = s*g;
				g = c*g;
				z = pythag(f,h);
				rv1[j] = z;
				c = f/z;
				s = h/z;
				f = x*c+g*s;
				g = g*c-x*s;
				h = y*s;
				y *= c;
				for (jj=0; jj<n; jj++) {
					x = v[jj][j];
					z = v[jj][i];
					v[jj][j] = x*c+z*s;
					v[jj][i] = z*c-x*s;
				}
				z = pythag(f,h);
				w[j] = z;
				if (z) {
					z = 1./z;
					c = f*z;
					s = h*z;
				}
				f = c*g+s*y;
				x = c*y-s*g;
				for (jj=0; jj<m; jj++) {
					y = a[jj][j];
					z = a[jj][i];
					a[jj][j] = y*c+z*s;
					a[jj][i] = z*c-y*s;
				}
			}
			rv1[l] = 0.;
			rv1[k] = f;
			w[k] = x;
		}
	}
}

static double pythag(double a, double b) {
	double absa = fabs(a);
	double absb = fabs(b);
	if (absa > absb) return absa*sqrt(1.+SQR(absb/absa));
	else return (absb == 0. ? 0. : absb*sqrt(1.+SQR(absa/absb)));
}

/*===================================================================================*/

static void NRfit(const vector<vector<double> > & f, const vector<double> & y, vector<double> & coeffs, double & chi, vector<double> & errors)
/* Given a set of data points x[1..ndat], y[1..ndat] with individual standard deviations
sig[1..ndat], use chi-2 minimization to fit for some or all of the coeficients a[1..ma] of
a function that depends linearly on a, y = sum(a_i*afunc_i(x)). The input array ia[1..ma]
indicates by nonzero entries those components of a that should be .tted for, and by zero entries
those components that should be held fixed at their input values. The program returns values
for a[1..ma], chi-2 = chisq, and the covariance matrix covar[1..ma][1..ma]. (Parameters
held fixed will return zero covariances.) The user supplies a routine funcs(x,afunc,ma) that
returns the ma basis functions evaluated at x = x in the array afunc[1..ma] */
{
	int i,j,k;
	double sum;

	int ndat= y.size();
	if (ndat==0) return;
	int ma = f[0].size();	// ma = nr of basis functions

	vector<vector<double> > beta(ma, vector<double>(1,0));
	vector<double> afunc(ma);
	vector<vector<double> > covar(ma, vector<double>(ma,0));

	for (i=0; i<ndat; i++) {
		afunc = f[i];
		for (j=0; j<ma; j++) {
			for (k=0; k<=j; k++)	covar[j][k] += afunc[j]*afunc[k];
			beta[j][0] += y[i]*afunc[j];
		}
	}
	for (j=0; j<ma; j++)
		for (k=j+1; k<ma; k++)	covar[j][k] = covar[k][j];

	gaussj(covar,ma,beta,1);

	coeffs.clear();
	for (int l=0; l<ma; l++) coeffs.push_back(beta[l][0]);
	double chisq = 0;
	errors.clear();
	for (i=0; i<ndat; i++) {
		afunc = f[i];
		for (sum=0.,j=0; j<ma; j++) sum += coeffs[j]*afunc[j];
		chisq += SQR(y[i]-sum);
		errors.push_back(y[i]-sum);
	}
	chi = sqrt(chisq/ndat);
}

#define SWAP(a,b) {swap=(a);(a)=(b);(b)=swap;}

void gaussj(vector<vector<double> > & a, int n, vector<vector<double> > & b, int m)
/* Linear equation solution by Gauss-Jordan elimination, equation (2.1.1) above. a[1..n][1..n]
is the input matrix. b[1..n][1..m] is input containing the m right-hand side vectors. On
output, a is replaced by its matrix inverse, and b is replaced by the corresponding set of solution vectors */
{
	int i,icol,irow,j,k,l,ll;
	double big,dum,pivinv,swap;
	
	vector<long> indxc(n); // The integer arrays ipiv, indxr, and indxc are used for bookkeeping on the pivoting
	vector<long> indxr(n);
	vector<long>  ipiv(n);

	for (j=0; j<n; j++) ipiv[j]=0;
	for (i=0; i<n; i++) {
		// This is the main loop over the columns to be reduced.
		big=0.0;
		for (j=0; j<n; j++) {
		// This is the outer loop of the search for a pivot element.
			if (ipiv[j] != 1) {
				for (k=0; k<n; k++) {
					if (ipiv[k] == 0) {
						if (fabs(a[j][k]) >= big) {
							big=fabs(a[j][k]);
							irow=j;
							icol=k;
						}
					}
				}
			}	
		}
		++(ipiv[icol]);
/* We now have the pivot element, so we interchange rows, if needed, to put the pivot
element on the diagonal. The columns are not physically interchanged, only relabeled:
indxc[i], the column of the ith pivot element, is the ith column that is reduced, while
indxr[i] is the row in which that pivot element was originally located. If indxr[i] .=
indxc[i] there is an implied column interchange. With this form of bookkeeping, the
solution b’s will end up in the correct order, and the inverse matrix will be scrambled
by columns */
		if (irow != icol) {
			for (l=0; l<n; l++) SWAP(a[irow][l],a[icol][l])
			for (l=0; l<m; l++) SWAP(b[irow][l],b[icol][l])
		}
		indxr[i]=irow; // We are now ready to divide the pivot row by the pivot element, located at irow and icol.
		indxc[i]=icol;
		if (a[icol][icol] == 0.0) cout<<"gaussj: Singular Matrix"<<"\t"<<icol<<"\t"<<endl;
		pivinv=1.0/a[icol][icol];
		a[icol][icol]=1.0;
		for (l=0; l<n; l++) a[icol][l] *= pivinv;
		for (l=0; l<m; l++) b[icol][l] *= pivinv;
		for (ll=0; ll<n; ll++)
			// Next, we reduce the rows except for the pivot one, of course.
			if (ll != icol) {
				dum=a[ll][icol];
				a[ll][icol]=0.0;
				for (l=0; l<n;l++) a[ll][l] -= a[icol][l]*dum;
				for (l=0; l<m;l++) b[ll][l] -= b[icol][l]*dum;
			}
	}
/* This is the end of the main loop over columns of the reduction. It only remains to unscramble
the solution in view of the column interchanges. We do this by interchanging pairs of columns in
the reverse order that the permutation was built up */
	for (l=n-1;l>=0;l--) {
		if (indxr[l] != indxc[l])
			for (k=0; k<n; k++)
				SWAP(a[k][indxr[l]],a[k][indxc[l]]);
	}
	// And we are done.
}

DRLIB_END_NAMESPACE
