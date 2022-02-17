
#include <numeric>
#include "firsttoinc.h"
#include "expt.h"   // for the exceptions

#include <cmath>
#include <complex>


#include "gpclosedforms/qgm_skewcalibration.h"
#include "gpclosedforms/optimization1.h"

#include "nag.h"
#include "nage04.h"

#include "gpbase/port.h"
#include "gpbase/removenagwarning.h"

using namespace std;

CC_BEGIN_NAMESPACE(ARM)


QGM_ParameterSet::QGM_ParameterSet(Optimization_Result_Set* r)
{
	_X = new std::vector<double>(*r->OptimalParamSet);
	_objective=r->OptimalObjective;
}

QGM_ParameterSet* QGM_CalibrateSkew(std::vector<double>* tagetVect,
	  std::vector<double>* weights,
      std::vector<double>* precisions,
	  std::vector<double>* InitVector,
	  std::vector<double>* LBoundVector,
	  std::vector<double>* UBoundVector,
	  int algorithm)
{
	int sizeVector = tagetVect->size();

	class objectiveFuntion : public Optimization_ObjectiveFuntion
	{
    private:
        std::vector<double> itsPresicions;
        std::vector<double> itsWieghts;
        std::vector<double> itsTagetVect;
	  
	public:
		void NAG_CALL operator() (Integer m, Integer n, 
			double x[], /// input
			double f[],	/// output (f(x))
			double fjac[],  /// output  (Df(x,i))
			Integer tdfjac, Nag_Comm *comm)
			
		{
			int i;
            double sumWeight = std::accumulate(itsWieghts.begin(),itsWieghts.end(),0.0);//itsWieghts.sum();
			for (int l =0; l<m*n; ++l)
				fjac[l] = 0;
			for(i=0;i<m;i++)
			{
				f[i] = 0;
                double precision = itsPresicions[i];
                double sqrtWeight = sqrt(itsWieghts[i]/sumWeight);
                double a = sqrtWeight/precision;
				for(int j=i; j<m; ++j)
				{               
					f[i]+=x[j];
					fjac[m*i+j] = 1.0/(m-i)*a;
				}
				f[i]=f[i]/double((m-i))*a;
				double target = itsTagetVect[i];
                f[i]-= target*a;
			}
		}
		objectiveFuntion(std::vector<double>& precisions, std::vector<double>& weights,vector<double>& tagetVect)
        :itsPresicions(precisions),
        itsWieghts(weights),
        itsTagetVect(tagetVect)
        {}	
	};
	objectiveFuntion func(*precisions,*weights,*tagetVect);

	Optimization_Result_Set* result= OptimizeWithDerivatives(&std::vector<double>(sizeVector,0.0),
		&std::vector<double>(sizeVector,0.0),
		&std::vector<double>(sizeVector,1.0),
		&func,
		InitVector,
		LBoundVector,
		UBoundVector,
		algorithm,
		false,
		//"C:\\temp\\QGMnag.txt",
		"",
		0.0001,  // Tolerance
		50);	// max_iter

	QGM_ParameterSet* setptr= new QGM_ParameterSet(result);
	delete result;
	return setptr;


}

CC_END_NAMESPACE()

/*---------------------------------------------------------------------------*/
/*---- End of file ----*/
