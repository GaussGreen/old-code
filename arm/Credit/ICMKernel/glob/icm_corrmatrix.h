
#ifndef _ICM_CORRELATION_MATRIX_H_
#define _ICM_CORRELATION_MATRIX_H_

#include <set>
#include <numeric>

#include "ICMKernel\glob\icm_betas_correlation.h"
#include "ICMKernel\glob\icm_enums.h"

struct StrOpt
{
	double * ptMeanLosses ;
	double * ptStdLosses ;
	ARM_Vector  * ptLosses ;
	ARM_Vector  * ptInvProb ;
	ARM_Vector  * ptInitialBetas ;
	ICM_QMatrix<double> * ptWeights ;
	ICM_QMatrix<double> * ptCorrels ;
} ;

double DistMatrix(int size, double * V, ICM_QMatrix<double> & M);
double* MatrixCorrel_FactorModel (StrOpt * ptSO, long MAX_ITER, double OPTIM_TOL);

/*********************************************************************************/
/*! \class  ICM_CorrMatrix icm_corrmatrix.h "icm_corrmatrix.h"
 *  \author D Pouponneau
 *	\version 1.0
 *	\date   April 2004
 *	\brief Creates a correlation matrix object 
/***********************************************************************************/

class ICM_CorrMatrix : public ICM_Beta_Correlation
{

private:
	
	ICM_QMatrix<double> itsMatrix;
	double	itsFixedCorrelation;
	bool	itsIsComputedBetasVector;

protected:
	inline void Init(void)
	{
		ICM_Beta_Correlation::Init(); 
		SetName(ICM_CORRMATRIX);
		// itsMatrix = NULL;
		itsFixedCorrelation = CREDIT_DEFAULT_VALUE;
		itsIsComputedBetasVector = false;
	}
private:
	 void Set(const ARM_Date& AsOf,
		 const string& name, 
		 const std::vector<std::string>& labels,
		 const ICM_QMatrix<double>&matrix) ;

public: 

	ICM_CorrMatrix() {Init();}		

	ICM_CorrMatrix(const ARM_Date& AsOf,
					const string& name, 
					const std::vector<std::string>& labels,
					const ICM_QMatrix<double>&matrix) 
	{
		Init();

		Set(AsOf,name,labels,matrix);
	};

	~ICM_CorrMatrix() {}


	// ----------------------------
	//	Copy of members data
	// ----------------------------
	void BitwiseCopy(const ARM_Object* src) ;


	// -------------
	//	Copy Method 
	// -------------
	void Copy(const ARM_Object* src) ;

	// --------------
	//	Clone Method
	// --------------
	ARM_Object* Clone(void) ;



	inline int GetSize() const { return itsMatrix.Getnbcols();}

	inline void SetMatrix(int size, double value) ;
	

	inline void SetMatrix(double** matrix, int size) ;

	void SetMatrix(ICM_QMatrix<double>* matrix) ;
	void SetMatrix(const ICM_QMatrix<double>&  matrix) ;

	inline const ICM_QMatrix<double>& GetMatrix() const {return itsMatrix;}
	inline ICM_QMatrix<double>& GetMatrix()  {return itsMatrix;}


	double GetCorrelationMatrix(const std::string& label1, const std::string& label2);
	void	PopulateCorrelationMatrix(ICM_QMatrix<double>&,const std::vector<std::string>&labels,int size); 

	ICM_CorrMatrix* ExtractCorrs(const std::vector<std::string>& );
	void SetCorrelation(const std::string& label1, const std::string& label2, double value);
	void ModifyCorrMatrixForBary(double beta, int UP);


	virtual double GetCorrelation(const std::string&  issuer1,
								  const std::string&  issuer2,
								  double maturity = CREDIT_DEFAULT_VALUE,
								  double strike = CREDIT_DEFAULT_VALUE,
								  double actualYF = CREDIT_DEFAULT_VALUE)
	{
		return GetCorrelationMatrix( issuer1,issuer2);
	}

	virtual double GetBeta(const std::string& issuer,
						   double maturity = CREDIT_DEFAULT_VALUE,
						   double strike = CREDIT_DEFAULT_VALUE,
						   double actualYF = CREDIT_DEFAULT_VALUE)
	{

		double Beta = ICM_Beta_Correlation::GetBeta(issuer);

		return (Beta);
	}


	virtual ARM_Vector* ComputeBetas(int nbissuers,
							 const std::vector<std::string>& labels,
							 const ARM_Vector& nominals, 
							 const ARM_Date & Maturity,
							 ARM_Model* model = NULL);

	// JLA : 
	ARM_Vector* ComputeBetas(int nbissuers,
							 const std::vector<std::string>& labels,
							 const ARM_Vector& Nominals, 
							 const ARM_Vector& RecoveryRates,
							 const ARM_Vector& DefProbAtMatu ) ; 


							 

	virtual ICM_QMatrix<double> ComputeCholeskyMatrix();

	void View(char* id, FILE* ficOut);

	virtual void ResetBetas() 
	{	
		itsIsComputedBetasVector = false; 
	}

	virtual ICM_Correlation* GenerateShiftCorrel(const std::string&  label,
												 qSENSITIVITY_TYPE typesensi,
												 double epsilon ) ;
	virtual ICM_Correlation* GenerateShiftBetas(const std::string& label,
												qSENSITIVITY_TYPE typesensi,
												double epsilon  ) 
	{
		//Par défaut on fait une sensi au beta
		ICM_Correlation* Correlation = GenerateShiftCorrel(label,
														   typesensi,
														   epsilon);

		return (Correlation);
	}

	inline void SetFixedCorrelation(double fixcorr) { itsFixedCorrelation = fixcorr; }
	inline double GetFixedCorrelation(void) { return itsFixedCorrelation; }
	void SortCorrelationMatrix(ICM_QMatrix<double> theCorrMatrix, 
							const std::vector<std::string>& labels,
								ICM_QMatrix<double>& theCorrMatrixSorted, ARM_Vector& vNoLabelsSorted);
	inline double Average()
	{
		if (itsFixedCorrelation != CREDIT_DEFAULT_VALUE) return itsFixedCorrelation;
		
		unsigned long  size = itsMatrix.Getnbrows();
		double Avg = 0.;

		for (int il1=0;il1<size; il1++)
		{for (int il2=0;il2<il1; il2++)
		{ Avg+=itsMatrix(il1,il2);}}

		Avg/=(size)*(size-1)/2;

		return (Avg);
	}

};

#endif