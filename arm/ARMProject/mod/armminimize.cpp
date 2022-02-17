#include "armminimize.h"

#include "armfrmmodelmixture.h"

// Take care ! We can't use Nag on WIN32 because header files
// are not available.

#ifdef WIN32
	#include "ournag.h"
	#include <nag_stdlib.h>
	#include <nage04.h>
	#include <nagx02.h>
#endif

#include <math.h>
#include <vector>
#include <algorithm>

ARM_Minimize::ARM_Minimize()
{
	Init();
}

ARM_Minimize::ARM_Minimize(
ARM_FRMModelMixture* model,
const ARM_Vector& lowerBound,
const ARM_Vector& upperBound)
{
	Init();

	itsModel = model;

	itsLowerBound = lowerBound;
	itsUpperBound = upperBound;
}

ARM_Minimize::~ARM_Minimize(void)
{
}

void ARM_Minimize::Init(void)
{
	itsModel = NULL;
}

void ARM_Minimize::BitwiseCopy(const ARM_Object* MinimizeNd)
{
	ARM_Minimize* src = (ARM_Minimize*) MinimizeNd;

	itsModel = src->itsModel;

	itsLowerBound = src->itsLowerBound;
	itsUpperBound = src->itsUpperBound;
}


void ARM_Minimize::Copy(const ARM_Object* MinimizeNd)
{
	ARM_Calibration::Copy(MinimizeNd);
  
	BitwiseCopy(MinimizeNd);
}


ARM_Object* ARM_Minimize::Clone(void)
{
  ARM_Minimize* theClone = new ARM_Minimize();
  
  theClone->Copy(this);
  
  return(theClone);
}

double uniforme()
{
	int rand_val=0;
	while (rand_val == 0)
	{
		rand_val = rand();
	}
	return rand_val/(1.0+RAND_MAX);
}

void ARM_Minimize::StochasticMinimisation(ARM_Vector& x, int index)
{
	int n = itsModel->GetN();

	double min = 0.0;
	double val = 0.0;

	for (int i = 0; i < 100; ++i)
	{
		std::vector<double> volParam(n);
		std::vector<double> spreadParam(n);

		for (int j = 0; j < n; ++j)
		{
			volParam[j] = itsLowerBound[j]+(itsUpperBound[j]-itsLowerBound[j])*uniforme();
		}

		std::sort(volParam.begin(),volParam.end());

		for (j = 0; j < n-1; ++j)
		{
			spreadParam[j] = itsLowerBound[n+j]+(itsUpperBound[n+j]-itsLowerBound[n+j])*uniforme();
		}

		ARM_Vector xVector(2*n-1);

		for (j = 0; j < n; ++j)
		{
			xVector[j] = volParam[j];
		}

		for (j = 0; j < n-1; ++j)
		{
			xVector[n+j] = spreadParam[j];
		}

		if (i == 0)
		{
			min = itsModel->SimpleFuncToMinimize(xVector,index);
			x = xVector;
		}
		else
		{
			val = itsModel->SimpleFuncToMinimize(xVector,index);
			if (min > val)
			{
				min = val;
				x = xVector;
			}
		}
	}
}

void ARM_Minimize::sortVolParam(ARM_Vector& x)
{	
	int n = itsModel->GetN();

	std::vector<double> volParam(n);

	for (int i = 0; i < n; ++i)
	{
		volParam[i] = x[i];
	}

	std::sort(volParam.begin(),volParam.end());

	for (i = 0; i < n; ++i)
	{
		x[i] = volParam[i];
	}
}

bool ARM_Minimize::testValidSolution(const ARM_Vector& x)
{
	int n = itsModel->GetN();

	for (int i = 0; i < n; ++i)
	{
		if (fabs(x[i]-itsLowerBound[i]) < 1e-2)
		{
			return false;
		}

		if (fabs(x[i]-itsUpperBound[i]) < 1e-2)
		{
			return false;
		}
	}

	for (i = 0; i < n-1; ++i)
	{
		if (fabs(x[n+i]-itsLowerBound[n+i]) < 1e-4)
		{
			return false;
		}

		if (fabs(x[n+i]-itsUpperBound[n+i]) < 1e-4)
		{
			return false;
		}
	}
	
	return true;
}

#if WIN32

struct ModelAndIndex
{
	int index;
	ARM_FRMModelMixture* model;
};

void NAG_CALL lsqNd(Integer m, double x[], double* objf,double g[],
                    Nag_Comm *comm)
{
	ModelAndIndex* modelAndIndex = (ModelAndIndex*) comm->p;

	ARM_Vector xVector(m);

	int i;

	for (i = 0; i < m; ++i)
	{
		xVector[i] = x[i];
	}

	*objf = modelAndIndex->model->SimpleFuncToMinimize(xVector,modelAndIndex->index);
}

#endif

void ARM_Minimize::Calibrate(ARM_Vector& x, int index)
{
#if WIN32
	Nag_E04_Opt options;

	static NagError fail;
	Nag_BoundType bound;
	Nag_Comm comm;

	double* var;
	double* g;
	double objf;
	double* bl;
	double* bu;

	int n = 2*itsModel->GetN()-1;
	
	var = new double [n];
	g = new double [n];
	bl = new double [n];
	bu = new double [n];


	bound = Nag_Bounds;

	int i,j;

	for (i = 0; i < n; ++i)
	{
		bl[i] = itsLowerBound[i];
		bu[i] = itsUpperBound[i];

		var[i] = x[i];
	}

	e04xxc(&options);

	// Décommenter les lignes s

	//strcpy(options.outfile,"c:\\log\\naglog.txt");

	//options.print_level = Nag_Soln_Iter_Full;

	fail.print = FALSE;

	struct ModelAndIndex modelAndIndex;

	modelAndIndex.model = itsModel;
	modelAndIndex.index = index;

	comm.p = &modelAndIndex;

	for (i = 0; i < 5; ++i)
	{
		StochasticMinimisation(x,index);

		for (j = 0; j < n; ++j)
		{
			var[j] = x[j];
		}

		e04jbc(n, &lsqNd, bound, bl, bu, var, &objf, g, &options, &comm, &fail);

		for (j = 0; j < n; ++j)
		{
			x[j] = var[j];
		}

		if (testValidSolution(x))
		{
			break;
		}
	}

	e04xzc(&options,"all",&fail);

	delete [] var;
	delete [] g;
	delete [] bl;
	delete [] bu;
#endif
}
