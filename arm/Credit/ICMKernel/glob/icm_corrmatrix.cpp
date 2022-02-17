
#include "ARMKernel\glob\firsttoinc.h"

#include "ICMKernel\glob\icm_corrmatrix.h"
#include "ICMKernel\mod\modelmulticurves.h"
#include "ICMKernel\util\icm_cholesky.h"
#include "ICMKernel\util\icm_utils.h"

#include "ICMKernel/glob/icm_maths.h" 
#include <nage04.h>

static void __stdcall confun(long n, long ncnlin, long needc[], double x[], double conf[], double conjac[], Nag_Comm *comm);
static void __stdcall objfun(long m, long n, double x[], double f[], double fjac[], long tdfjac, Nag_Comm *comm);

//	--------------------------------------------------------------------------------------------
void 
ICM_CorrMatrix::Set(const ARM_Date& AsOf,
		 const string& name, 
		 const std::vector<std::string>& labels,
		 const ICM_QMatrix<double>&matrix) 
{
	ARM_Vector  betas (labels.size(),0.);

	ICM_Beta_Correlation::Set(AsOf,name, betas,labels,(ARM_IRIndex*)0,(ARM_IRIndex*)0);
	itsMatrix=matrix; 
}

// ----------------------------
//	Copy of members data
// ----------------------------
void ICM_CorrMatrix::BitwiseCopy(const ARM_Object* src)
{
	ICM_CorrMatrix* CorrMatrix = (ICM_CorrMatrix *) src;
 
	itsMatrix = CorrMatrix->itsMatrix ;
	
	itsFixedCorrelation=CorrMatrix->itsFixedCorrelation ;
	itsIsComputedBetasVector=CorrMatrix->itsIsComputedBetasVector ;
}

// -------------
//	Copy Method 
// -------------
void ICM_CorrMatrix::Copy(const ARM_Object* src)
{
	ICM_Beta_Correlation::Copy(src);

	BitwiseCopy(src);
}

// --------------
//	Clone Method
// --------------
ARM_Object* ICM_CorrMatrix::Clone(void)
{
 ICM_CorrMatrix* theClone = new ICM_CorrMatrix();

 theClone->Copy(this);

 return(theClone);
}
// --------------------------------------------------------------------------------------------------------
// virtual 
ICM_Correlation* ICM_CorrMatrix::GenerateShiftCorrel(const std::string& label,
											 qSENSITIVITY_TYPE typesensi,
											 double epsilon )
{
	ICM_CorrMatrix* correl = (ICM_CorrMatrix*) Clone(); 

	if (itsFixedCorrelation != CREDIT_DEFAULT_VALUE)
	{
		correl->SetFixedCorrelation(GetFixedCorrelation()+epsilon);
		return correl;
	}

	if ((typesensi != ICMBETA_TYPE) && (typesensi != ICMCORRELATION_TYPE)) return (correl);

	int i=0,j=0;
	int nolabel = GetLabelNo(label);

	ICM_QMatrix<double>& Matrix = correl->GetMatrix();
	double defcorr = 0.999;

	if (nolabel == CREDIT_DEFAULT_VALUE) 
	{
	for (i=0; i<itsMatrix.Getnbrows(); i++)
		for (j=i; j<itsMatrix.Getnbcols(); j++)
		{
			double value = MIN(Matrix.Getvalue(i,j)+epsilon,defcorr);
			value = MAX(value,-defcorr);

			Matrix.SetValue(i,j,value);
			Matrix.SetValue(j,i,value);
		}

	for (i=0; i<itsMatrix.Getnbrows(); i++)
		Matrix.SetValue(i,i,1.);

	}
	else
	{
		for (j=0; j<itsMatrix.Getnbrows(); j++)
		{
			double value = MIN(Matrix.Getvalue(nolabel,j)+epsilon,defcorr);
			value = MAX(value,-defcorr);

			Matrix.SetValue(nolabel,j,value);
			Matrix.SetValue(j,nolabel,value);
		}

		Matrix.SetValue(nolabel,nolabel,1.);
	}
	
	correl->ResetBetas();

	return (correl);
}

//	--------------------------------------------------------------------------------------------
void 
ICM_CorrMatrix::SetMatrix(int size, double value)
{
 	itsMatrix.Resize(size,size); 
	itsMatrix.Fill(value); 
	for (int i=0; i<size; i++)
		itsMatrix.SetValue(i,i,1.);
}

// --------------------------------------------------------------------
// Returns Correlation Matrix for a labels vector 
// --------------------------------------------------------------------
inline void 
ICM_CorrMatrix::PopulateCorrelationMatrix(ICM_QMatrix<double>& mat,
										  const std::vector<std::string>& labels, 
										  int size) 
{
//	QUANTIFYER("ICM_CorrMatrix::PopulateCorrelationMatrix"); 
	if (size==0) 
		ICMTHROW(ERR_INVALID_MODEL,"ERROR : PopulateCorrelationMatrix : size is 0 "); 

	mat.Resize(size,size); 
	
	unsigned int j,k; 
	for ( j=0; j<size; j++)
	{
		mat(j,j)=1. ;
		for ( k=j ; k<size; k++)
			mat(k,j)=mat(j,k)= GetCorrelationMatrix(labels[j],labels[k]);
	}

}
// --------------------------------------------------------------------
// Set Correlation Matrix
// --------------------------------------------------------------------
void ICM_CorrMatrix::SetCorrelation(const std::string& label1, const std::string& label2, double value)
{
	int i;
	int row=-1, col=-1;
	int size=GetSize();
	
	if (label1==label2)
		return;

	for ( i=0; i<size; i++)
	{	
		if ( GetLabel(i)==label1)  row = i;
		else if (GetLabel(i)==label2)  col = i;
		if ((row>=0) && (col>=0)) break;
	}	

	if ((row>=0) && (col>=0))
	{	
		itsMatrix.SetValue(row,col,value);
		itsMatrix.SetValue(col,row,value);
		return;
	}
	ICMTHROW(ERR_INVALID_ARGUMENT,"SetCorrelation: can't find "<<label1<<" "<<label2); 

}

// --------------------------------------------------------------------
// Barycentric Adjustment for correlation matrix
// --------------------------------------------------------------------
void ICM_CorrMatrix::ModifyCorrMatrixForBary(double beta, int UP)
{
	double alpha = 1.-beta;
	int size = itsMatrix.Getnbcols();
	double tmp=0.;

	for (int i=0; i<size; i++)
	{
		itsMatrix.SetValue(i,i,1.);

		for (int j=0; j<i; j++)
		{
			tmp = itsMatrix.Getvalue(i,j);

			if (UP)
				tmp = alpha*tmp + (1-alpha);
			else
				tmp = alpha*tmp - (1-alpha);

			itsMatrix.SetValue(i,j,tmp);
			itsMatrix.SetValue(j,i,tmp);
		}
	}

}


// --------------------------------------------------------------------
// Extract Correlation Matrix 
// --------------------------------------------------------------------
ICM_CorrMatrix* ICM_CorrMatrix::ExtractCorrs(const std::vector<std::string>&labels )
{
	int j,k;
	ICM_QMatrix<double> newmatrix(labels.size(),labels.size()); 

 
	for ( j=0; j<labels.size(); j++)
		for ( k=j; k<labels.size(); k++)
		  newmatrix(j,k)  = GetCorrelation(labels[j].c_str(),labels[k].c_str());				

	for ( j=0; j<labels.size(); j++)
		for ( k=j; k<labels.size(); k++)
		  newmatrix(k,j)= newmatrix(j,k);		



	ICM_CorrMatrix* corrmatrix = new ICM_CorrMatrix(GetAsOfDate(),GetStructName(),labels,newmatrix);

 	return (corrmatrix);
}


// --------------------------------------------------------------------
// View Method
// --------------------------------------------------------------------
void ICM_CorrMatrix::View(char* id, FILE* ficOut)
{
	FILE* fOut;
	char  fOutName[200];
	int i = 0;

	if ( ficOut == NULL )
	{
		ARM_GetViewFile(ARM_VIEW_FILE, id, fOutName);
		(void) unlink(fOutName);
		fOut = fopen(fOutName, "w"); 
	}
	else
	{
		fOut = ficOut;
	} 

	fprintf(fOut, "\n ======> Correlation Matrix :\n\n");
		
	if (!itsMatrix.IsEmpty())
	{
	int size = GetSize();
	int k =0;

		for (i = 0; i<size; i++)
		{
			fprintf(fOut, "%s\t", GetLabel(i).c_str()); 

			for (k = 0; k<size; k++)
			{
				fprintf(fOut, "\t %f ", itsMatrix.Getvalue(k,i)); 
			}

			fprintf(fOut, "\n");
		}

	}

	ICM_Beta_Correlation::View(id,ficOut);

	if ( ficOut == NULL )
	{
		fclose(fOut);
	}
}


// --------------------------------------------------------------------
// Returns Correlation for 2 issuers
// --------------------------------------------------------------------
double ICM_CorrMatrix::GetCorrelationMatrix(const std::string& label1, const std::string& label2)  
{
	int i;
	int row=-1,col=-1;

	if (label1==label2)
		return 1.0;

	for ( i=0; i<GetSize(); i++)
	{	
		if (!strcmp(GetLabel(i).c_str(),label1.c_str()))  row = i;
		else if (!strcmp(GetLabel(i).c_str(),label2.c_str()))  col = i;
		if ((row>=0) && (col>=0)) break;
	}	

	if ((row>=0) && (col>=0))
		return itsMatrix.Getvalue(row,col);

	ICMTHROW(ERR_INVALID_MODEL,"GetCorrelationMatrix : unknown "<<label1<<" or "<<label2); 

	return 0.;
}

// --------------------------------------------------------------------
// Computes Betas vector for a list of issuers
// --------------------------------------------------------------------
ARM_Vector* ICM_CorrMatrix::ComputeBetas(int nbissuers,
										 const std::vector<std::string>& labels,
										 const ARM_Vector& nominals,
										 const ARM_Date &Maturity,
										 ARM_Model* Original)
{
	
	if (itsIsComputedBetasVector) return GetBetas(labels,nbissuers);
	
	{
		ICM_ModelMultiCurves* model = (ICM_ModelMultiCurves*)Original;
		ARM_Vector RecoveryRates(nbissuers); 
		ARM_Vector DefProbAtMatu(nbissuers); 
		for(unsigned int i=0;i<nbissuers;i++) 
		{
			RecoveryRates[i]= model->GetRecoveryRate(labels[i]) ; 
			DefProbAtMatu[i]= model->GetDefaultCurve(labels[i])->DefaultProba(Maturity); 
		}
		return ComputeBetas(nbissuers,labels,
				nominals,
				RecoveryRates,
				DefProbAtMatu) ; 
	}

}

// --------------------------------------------------------------------
// Computes Betas vector for a list of issuers
// --------------------------------------------------------------------
ARM_Vector *ICM_CorrMatrix::ComputeBetas(int nbissuers,
										 const std::vector<std::string>& labels,
										 const ARM_Vector& nominals,
										 const ARM_Vector& RecoveryRates,
										 const ARM_Vector& DefProbAtMatu) 
{
		if (itsIsComputedBetasVector) return GetBetas(labels,nbissuers);
	
	int il = 0;
	ARM_Vector * Betas = new ARM_Vector(nbissuers,0.) ;	
	unsigned int i , j , k;
	if (itsFixedCorrelation != CREDIT_DEFAULT_VALUE)
	{
		double aux = sqrt(itsFixedCorrelation);
		for (k=0; k<nbissuers; k++)
			Betas->Elt(k)=aux;

		SetBetas(*Betas,labels);
		itsIsComputedBetasVector = true;
		return (Betas);
	}

	//	Here we extract the Symetric correl matrix according to labels list 
	// double** corrmatrix = GetCorrelationMatrix(labels,nbissuers);
	ICM_QMatrix<double>	theCorrMatrix; 
	
	PopulateCorrelationMatrix(theCorrMatrix,labels,nbissuers); 
	ICM_QMatrix<double>	theSortedMatrix(nbissuers,nbissuers,(double) 0.0); 
	ARM_Vector vNoLabelsSorted(nbissuers);
	SortCorrelationMatrix(theCorrMatrix,labels,theSortedMatrix, vNoLabelsSorted);
	

	ARM_Vector Losses(nbissuers) ;
	for (i = 0; i < nbissuers; i++)		Losses[i] = (1. - RecoveryRates[(int)vNoLabelsSorted[i]])  * nominals[(int)vNoLabelsSorted[i]];

	ARM_Vector InvProb(nbissuers) ;

	for (i = 0; i < nbissuers; i++)
		InvProb[i] = NAG_deviates_normal( DefProbAtMatu[(int)vNoLabelsSorted[i]]) ;
	
	double averageloss = 0. ; 

	ICM_QMatrix<double> MatLosses(nbissuers,nbissuers) ;
	
	for (i = 0; i < nbissuers; i++)
	{
		double value; 
		for (unsigned int j = 0; j < i; j++)
		{
			value=Losses[i] * Losses[j]  ;
			MatLosses(i,j)	= value  ;
			averageloss += value*value ;
		}
	}
	

	averageloss = sqrt(averageloss) ;

	//  We duplicate the initial correlation matrix in an 
	//  appropriate structure. 
// JLA	ICM_QMatrix<double> InitialCorrelations(nbissuers,nbissuers) ;
	ICM_QMatrix<double> InitialCorrelations(theSortedMatrix) ;

	ICM_QMatrix<double> JointProb(nbissuers,nbissuers) ;


	double stdevLosses = 0. ;

	for (i = 0; i < nbissuers; i++)
	{
		JointProb(i,i) = Losses[i] * Losses[i] * InvProb[i] ;			
		
		MatLosses(i,i) = 1 ;

		for (j = 0; j < i; j++)
		{
			JointProb(i,j)  = Losses[i] * Losses[j] * 
				NAG_bivariate_normal_dist(InvProb[i], InvProb[j], InitialCorrelations(i,j)) ; 
			
			MatLosses(i,j)  = 	MatLosses(i,j) / averageloss ;
			MatLosses(j,i)  = 	MatLosses(i,j) ;
			JointProb(j,i)  = JointProb(i,j) ;
			stdevLosses     += JointProb(i,j) ;				
		}

	}

	ARM_Vector InitBetasProv(nbissuers) ; 

	for (i = 0; i < nbissuers; i++)
			InitBetasProv[i] = 0.5 ;

	StrOpt SO ;
	SO.ptMeanLosses   = &averageloss ;
	SO.ptStdLosses    = &stdevLosses ;
	SO.ptLosses       = &Losses ;
	SO.ptInvProb      = &InvProb ;
	SO.ptInitialBetas = &InitBetasProv ; 
	SO.ptWeights      = &MatLosses ;
	SO.ptCorrels      = &InitialCorrelations ;	
			
	long N = 20 ;
	double Epsilon = 0.2 ;				// Condition: 0 = < Epsilon < 1

	double * vecteur = MatrixCorrel_FactorModel(&SO, N, Epsilon) ;
	
	for (j = 0; j < nbissuers; j++)
				Betas->Elt((int)vNoLabelsSorted[j]) = vecteur[j];

	delete[] vecteur; 
	SetBetas(*Betas,labels);
	itsIsComputedBetasVector = true;
	return Betas;
}


// --------------------------------------------------------------------
// Return Cholesky Matrix
// --------------------------------------------------------------------

ICM_QMatrix<double> ICM_CorrMatrix::ComputeCholeskyMatrix()
{
//	QUANTIFYER("ICM_CorrMatrix::ComputeCholeskyMatrix"); 
	int size = GetSize();

	ICM_QMatrix<double> CholeskyMatrix(size,size,0.);

    

	for (int i = 0; i < size; i++)
	{

		for (int j = 0; j < size; j++)
		{

		CholeskyMatrix.SetValue(i,j,itsMatrix.Getvalue(i,j));

		}
	}
	double tmp=0.;

	cholesky Obj=cholesky();

	Obj.init(CholeskyMatrix,size);

	Obj.choldc(size);


		for ( i = 0; i < size; i++)
		{
		CholeskyMatrix.SetValue(i,i,Obj.m_eigenvalues[i]);

			for (int j = i+1; j < size; j++)
			{

			CholeskyMatrix.SetValue(i,j,0.);
			}


		}


	for ( i = 0; i < size; i++)
	{

		for (int j = 0; j < i; j++)
		{
			tmp=Obj.m_matrix(i,j);
			CholeskyMatrix.SetValue(i,j,Obj.m_matrix(i,j));

		tmp=CholeskyMatrix.Getvalue(i,j);

		}
	}

    	
   

	return CholeskyMatrix;
}
// sort on the sum of row vector and an ascii sur of names
void ICM_CorrMatrix::SortCorrelationMatrix(ICM_QMatrix<double> theCorrMatrix, 
										   const std::vector<std::string>& labels,
								ICM_QMatrix<double>& theCorrMatrixSorted, ARM_Vector& vNoLabelsSorted)
{
	int i=0, c=0;
	int ascii = 0, asTmp =0; 
	ARM_Vector vtemp;
	bool nameSort = false;
	
	if ( theCorrMatrix.GetSize() == 0){
		ICMTHROW(ERR_INVALID_DATA,"ICM_CorrMatrix::SortCorrelationMatrix No Data in the matrix");
	}
	vtemp.Resize(theCorrMatrix.GetSize());
	vector<string> vsLabels(theCorrMatrix.GetSize());
	
	
	for (i=0; i< theCorrMatrix.GetSize(); i++){
		ARM_Vector* pTmpVector = theCorrMatrix.RowAsVector(i);
		
		double sum = 0;
		for (int il=0;il<pTmpVector->GetSize();il++) {sum +=pTmpVector->Elt(il);}

		vtemp[i] = sum;
		vsLabels[i]= string(labels[i]); // set labels 
		delete pTmpVector;

	}
	// sort on sum values
	//ARM_Vector v2temp = vtemp.Sort(vNoLabelsSorted, K_DECREASING);
	ARM_Vector v2temp = SortWithConstraints(vtemp,vsLabels, vNoLabelsSorted, K_DECREASING);

	theCorrMatrixSorted.Permut(theCorrMatrix, vNoLabelsSorted);
}

double DistMatrix(int size, double * V, ICM_QMatrix<double> & M)
{
	double sum = 0. ;

	for (unsigned int i = 0; i < size; i++)
	{
		for (unsigned int j = 0; j < i; j++)
		{
			//JLA : unused double m = M(i,j) ;
			//JLA:	unused double v = V[i]*V[j] ;
			double dif = (M(i,j) - V[i]*V[j]) ;
			sum += dif * dif ;
		}
	}
	
	return sqrt(sum) ;

}


static void __stdcall confun(long n, long ncnlin, long needc[], double x[], double conf[], double conjac[], Nag_Comm *comm)
{

	comm->flag = 0 ;
	needc[0] = 1 ;
	conf[0]  = 0.; 
	
	StrOpt * ptSO = (StrOpt *)comm->p ;

	for (unsigned int i = 0; i < ptSO->ptCorrels->Getnbrows(); i++)
	{				
		for (unsigned int j = 0; j < i; j++)
		{
			double cour = (ptSO->ptLosses->Elt(i)) * (ptSO->ptLosses->Elt(j)) * 
				NAG_bivariate_normal_dist(ptSO->ptInvProb->Elt(i), ptSO->ptInvProb->Elt(j), x[i]*x[j]) ; 			

			conf[0]     += cour;				
		}
	}

	conf[0] /= (*ptSO->ptMeanLosses); 
} 

static void __stdcall objfun(long m, long n, double x[], double f[], double fjac[], long tdfjac, Nag_Comm *comm)
{

#define FAKHER_FJAC(I,J) fjac[(I)*tdfjac +(J)]
	
	/* Local variables */
	
	int i, j, k, jmax, jmin;
	
	/* Function to evaluate the objective subfunctions
	/* and their 1st derivatives.*/
	comm->flag = 2 ;

	i = 0 ;
	j = 1 ;
	
	jmax = n-1 ;
	jmin = 1;
	StrOpt * ptSO = (StrOpt *)comm->p ;
	
	for (k = 0; k < m; ++k)
	{
		/* Evaluate objective subfunction f(k) */
		
		if (j <= jmax)
		{
			f[k] =  x[i]*x[j] * ptSO->ptWeights->Getvalue(j,i) ;
			j++ ;
		}
	
		else
		{
			i = i+1 ;
			jmin = jmin + 1 ;
			j = jmin ;
			f[k] =  x[i]*x[j] * ptSO->ptWeights->Getvalue(j,i) ;
			j++ ;
		}
	}
	
	
	// Computation of FJAC , Hessian Matrix, i.e FJAC(i,j) = df[i]/dx[j]
	//		  _									      _
	//		 |x[1]   x[0]   0    0    0   0   0    0   |
	//		 |x[2]    0   x[0]   0    0   0   0    0   |
	//		 |x[3]    0     0   x[0]      0   0    0   |
	//		 |	.	  .	    .    .    .   .   .    .   |
	//		 |	.	  .	    .    . 	  .	  .	  .    .   |
	//		 |	.	  .	    .    .	  .	  .	  .    .   |
	//		 |x[n]    0     0    0    0   0   0  x[0]  |
	//       | 0    x[2]  x[1]   0    0   0   0    0   |
	//		 | 0    x[3]    0   x[1]  0   0   0    0   |
	//FJAC = |	.	  .	    .	 .	  .	  .	  .    .   |
	//		 |	.	  .	    .	 .	  .	  .	  .    .   |
	//		 |	.	  .	    .	 .	  .	  .	  .    .   |
	//		 | 0    x[n]    0    0    0   0   0  x[1]  |
	//		 | 0     0     x[3] x[2]  0   0   0    0   |  
	//		 | 0     0     x[4]  0   x[2] .   0    0   |
	//		 |	.	  .	    .	 .	  .	  .	  .    .   |
	//		 |	.	  .	    .	 .	  .	  .	  .    .   |
	//		 |	.	  .	    .	 .	  .	  .	  .    .   |
	//		 |	.	  .	    .	 .	  .	 x[n] .  x[n-2]|
	//		 | 0     0      0    0    0   0 x[n] x[n-1]|
	//		  -									      -

	int N_iter = n-1 ;
	int N_rows = 0 ;
	k = 0 ;
	

	while (N_iter >0)
	{
		j = 0 ;
		for (i = N_rows + 0; i < N_rows + N_iter; ++i) 
		{
			FAKHER_FJAC(i,k) = x[j+1] * ptSO->ptWeights->Getvalue(k,j+1) ;
			FAKHER_FJAC(i,k+j+1) = x[k] * ptSO->ptWeights->Getvalue(k,j+1) ;
			j++;
		}

		k = k + 1 ;
		N_rows += N_iter ;
		N_iter-- ;
	}

} 


double* MatrixCorrel_FactorModel (StrOpt * ptSO, long MAX_ITER, double OPTIM_TOL)
{
//	QUANTIFYER("MatrixCorrel_FactorModel") ;

	const unsigned int n      = ptSO->ptCorrels->Getnbrows() ;           
	const unsigned int m      = n*(n-1)/2 ;
	const unsigned int nclin  = 0 ;
	const unsigned int ncnlin = 1 ;      
	const unsigned int tda    = n ;
	const unsigned int tdfjac = n ; 

	
	/* Local variables */

	double * a = 0 ;
	double * f    = new double[m] ;
   	double * y    = new double[m] ; 
	double * bl   = new double[n+ncnlin] ; 
	double * bu   = new double[n+ncnlin] ;
	double * x    = new double[n] ;	
	double * fjac = new double[m*tdfjac] ;
	
	double * res    = new double[n] ;	
	double objf ;

	NagError fail; 
	INIT_FAIL(fail); 
	

		/* the y vector of the objective */
	unsigned int jmax, jmin ;
	unsigned int i = 0 ;
	unsigned int j = 1 ;
	
	jmax = n - 1 ;
	jmin = 1 ;
	
	for (unsigned int k = 0; k < m; ++k)
	{
		/* Evaluate objective subfunction f(k) */
		
		if (j <= jmax)
		{
			y[k] = ptSO->ptWeights->Getvalue(j,i) * ptSO->ptCorrels->Getvalue(j,i) ;
			j++;
		}
	
		else
		{
			i = i+1 ;
			jmin = jmin + 1;
			j = jmin;
			y[k] = ptSO->ptWeights->Getvalue(j,i) * ptSO->ptCorrels->Getvalue(j,i) ;
			j++;
		}
	}

		/* lower bounds */
		for (i=0; i<n; i++)
			bl[i] = -0.99;
				
		/* upper bounds */

		for (i=0; i<n; i++)
			bu[i] = 0.99;
		
		bl[n] = *(ptSO->ptStdLosses)/(*ptSO->ptMeanLosses) ; 
		bu[n] = *(ptSO->ptStdLosses)/(*ptSO->ptMeanLosses) ; 
		
		/* the initial point x */

		for (i=0; i<n; i++)
			x[i] = (ptSO->ptInitialBetas->Elt(i)) ;  
							
		/* Solve the problem */


		Nag_Comm comm ; 
		comm.p = (Pointer)(ptSO) ;

		Nag_E04_Opt options ;
		ICM_NAGSILENT(options); 
		options.max_iter = MAX_ITER ;
		options.optim_tol = OPTIM_TOL ;
		options.con_deriv = FALSE ;

			e04unc(m, n, nclin, ncnlin, a, tda, bl, bu, y, objfun,
			   confun, x, &objf, f, fjac, tdfjac, &options,
			   &comm, &fail) ;

		double criteria = DistMatrix(n, x, *(ptSO->ptCorrels)) ;

		unsigned int cmpt = 1 ;

		while ((criteria > 0.001) && (cmpt < 15))
		{

			e04unc(m, n, nclin, ncnlin, a, tda, bl, bu, y, objfun,
			   confun, x, &objf, f, fjac, tdfjac, &options,
			   &comm, &fail) ;
			criteria = DistMatrix(n, x, *(ptSO->ptCorrels)) ;
			cmpt++ ;
		}

		ICM_NAGFREE(options); 
		
		for (i=0; i<n; i++)
			res[i] = x[i] ; 
		
		delete [] x ; delete [] fjac ; delete [] f ; delete [] y ; delete [] bl ; delete [] bu ;

		
		return (res) ;


	if (res)
		delete [] res ; 
	res = NULL ;					

	return (NULL) ;

}




