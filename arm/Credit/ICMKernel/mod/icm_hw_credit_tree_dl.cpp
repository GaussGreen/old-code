#include "firsttoinc.h"

#include "icm_hw_credit_tree_dl.h"
#include <iomanip>

// -------------------------------------------------------------------
// Diffusion des probabilités dans l'arbre
// -------------------------------------------------------------------
void ICM_HWTree_intensity::DiffuseTree(bool corp)
{
	DiffuseTree(0,0,corp);
}
// -------------------------------------------------------------------
// Diffusion des probabilités à partir d'un noeud
// -------------------------------------------------------------------
void ICM_HWTree_intensity::DiffuseTree(int i0,int j0,bool corp)
{
	int i,j;

	for (i=0;i<itsTimeStep.size();i++)
	{for (j=0;j<=i;j++)
	{GetTree().data(i,j).proba=0.;}}

	for (i=i0;i<itsTimeStep.size();i++)
	{
		double sumproba = 0.;
		double EP_sumproba = 0.;

		for (j=j0;j<=i-i0+j0;j++)
		{
			double Jsize=0.,Intensity=0.,proba=0.,deltat=0.;
			if ((i>i0) && (j>j0))
			{
				Jsize=GetTree().data(i-1,j-1).jumpsize;
				Intensity=GetTree().data(i-1,j-1).lambda;
				deltat =GetTree().Time(i)-GetTree().Time(i-1);
			}
			else if (i>i0)
			{
				Jsize=GetTree().data(i-1,0).jumpsize;
				Intensity=GetTree().data(i-1,0).lambda;
				deltat =GetTree().Time(i)-GetTree().Time(i-1);
			}

			if ((i==i0) && (j==j0))
				{proba = 1.;}
			else if ((i==j) && (i>j0))
			{
				if (corp)
					proba = GetTree().data(i-1,j-1).proba*Intensity*deltat;
				else
					proba = GetTree().data(i-1,j-1).EP_i_j*Intensity*deltat;
			}
			else if (j==0)
			{
				if (corp)
					proba = GetTree().data(i-1,j).proba*(1.-Intensity*deltat);
				else
					proba = GetTree().data(i-1,j).EP_i_j*(1.-Intensity*deltat);
			}
			else
			{
				if (corp)
					proba = GetTree().data(i-1,j-1).proba*Intensity*deltat +
						GetTree().data(i-1,j).proba*(1.-Intensity*deltat);
				else
					proba = GetTree().data(i-1,j-1).EP_i_j*Intensity*deltat +
						GetTree().data(i-1,j).EP_i_j*(1.-Intensity*deltat);
			}

			if (corp)
				GetTree().data(i,j).proba=proba;
			else
				GetTree().data(i,j).EP_i_j=proba;

			sumproba+= proba;
		}

		if (CHECK_EQUAL(sumproba,1.)==false)
		{
			for (j=j0;j<=i-i0+j0;j++)
			{
				if (corp)
				GetTree().data(i,j).proba=GetTree().data(i,j).proba/sumproba;
				else
				GetTree().data(i,j).EP_i_j=GetTree().data(i,j).EP_i_j/sumproba;
			}
		}
	}

}

double ICM_HWTree_intensity::SumJumps(int& i)
{
	int j;
	double esp=0.,V=0.;

	if (i>=itsTimeStep.size())
		return 0.;

	for (j=0;j<=i;j++)
	{
		esp += GetTree().data(i,j).jumpsize * GetTree().data(i,j).proba;
	}

	return (esp);

}

double ICM_HWTree_intensity::CorpEsp(int& i,bool corp)
{
	int j;
	double esp=0.,V=0.;

	if (i>=itsTimeStep.size())
		return 0.;

	for (j=0;j<=i;j++)
	{
		if (corp) {
		V=GetTree().data(i,j).theta+GetTree().data(i,j).jumpsize;
		esp += V * GetTree().data(i,j).proba;}
		else
			{
		V=GetTree().data(i,j).EP_i_j;
		esp += V * GetTree().data(i,j).proba;}
	}

	return (esp);

}

void ICM_HWTree_intensity::SetTheta(int& period,double value)
{
	int j;

	if (period>=itsTimeStep.size())
		return ;

	for (j=0;j<=period;j++)
	{GetTree().data(period,j).theta=value;}

}

double ICM_HWTree_intensity::SumProba(int& period)
{
	int j;
	double value=0.;

	if (period>=itsTimeStep.size())
		return 0.;

	for (j=0;j<=period;j++)
	{value+=GetTree().data(period,j).proba;}

	return value;
}


//---------------------------------------------------------------------------
// View method
//---------------------------------------------------------------------------
void ICM_HWTree_intensity::View(char* id, FILE* ficOut)
{
	unsigned int i=0,j=0;
	FILE* fOut;
	char  fOutName[200];

	if ( ficOut == NULL )
	{
		ARM_GetViewFile(ARM_VIEW_FILE, id, fOutName);
		(void) unlink(fOutName);
		fOut = fopen(fOutName, "w"); 
	}
	else	fOut = ficOut;
	


	fprintf(fOut,"\n");
	fprintf(fOut,"\n");

	std::stringstream sstr ;
	sstr<<"\t\t\t ----------------- Hull White Binomial Tree ----------------- \n\n"<<std::endl ;
		
	if (view_all)
	{
	sstr<<"-- Dumping Tree : (lambda,theta,jumpsize,proba) "<<std::endl ;
	for(i=0;i<itsHWTree.depth();i++)
	{
		sstr<<(double)itsTimeStep[i]<<"\t"; 
		for(j=0;j<=i;j++) 
			sstr<<"(" 
			<<std::setw(4)<<std::setw(4)<<itsHWTree.data(i,j).lambda
			<<","<<std::setw(4)<<itsHWTree.data(i,j).theta
			<<","<<std::setw(4)<<itsHWTree.data(i,j).jumpsize
			<<","<<std::setw(4)<<itsHWTree.data(i,j).proba
			<<")"<<"\t"; 

		sstr<<std::endl; 
	}
	sstr<<"-- End of Dumping Tree\n"<<std::endl ;
	}

	if (view_theta_only)
	{
	sstr<<"-- Dumping Tree : (theta) "<<std::endl ;
	for(i=0;i<itsHWTree.depth();i++)
	{
		sstr<<(double)itsTimeStep[i]<<"\t"; 
		sstr<<std::setw(4)<<itsHWTree.data(i,0).theta; 
		sstr<<std::endl; 
	}
	sstr<<"-- End of Dumping Tree\n"<<std::endl ;
	}

	if (view_jumpsize_only)
	{
	sstr<<"-- Dumping Tree : (JumpSize) "<<std::endl ;
	{
		for(j=0;j<=(itsHWTree.depth()-1);j++)
		{
		sstr<<(double)j<<"\t"; 
		sstr<<std::setw(4)<<itsHWTree.data(itsHWTree.depth()-1,j).jumpsize; 
		sstr<<std::endl; 
		}
	}
	sstr<<"-- End of Dumping Tree\n"<<std::endl ;
	}

	if (view_lambda_only)
	{
	sstr<<"-- Dumping Tree : (lambda) "<<std::endl ;
	for(i=0;i<itsHWTree.depth();i++)
	{
		sstr<<(double)itsTimeStep[i]<<"\t"; 
		for(j=0;j<=i;j++) 
		{sstr<<std::setw(4)<<itsHWTree.data(i,j).lambda<<"\t";}

		sstr<<std::endl; 
	}
	sstr<<"-- End of Dumping Tree\n"<<std::endl ;
	}

	if (view_proba_only)
	{
	sstr<<"-- Dumping Tree : node probability "<<std::endl ;
	for(i=0;i<itsHWTree.depth();i++)
	{
		sstr<<(double)itsTimeStep[i]<<"\t"; 
		for(j=0;j<=i;j++) 
		{sstr<<std::setw(4)<<itsHWTree.data(i,j).proba<<"\t";} 
		sstr<<std::endl; 
	}
	sstr<<"-- End of Dumping Tree\n"<<std::endl ;
	}

	sstr<<std::endl; 

	fprintf(fOut,"%s\n",sstr.str().c_str()); 
	fprintf(fOut, "\n");

	if ( ficOut == NULL )fclose(fOut);
}

double ICM_HWTree_intensity::EC_Node_until_maturity(int i,int j,vector<double>& Epdef)
{
	double esp=0.,V=0.;

	if (i>=itsTimeStep.size())
		return 0.;

	int size=itsTimeStep.size();
	int sizemax = size;
	int i2=0,j2=0;
	double theta = 0.;
	double jumpsize = 0.;

	DiffuseTree(i,j,true);
	
	for (i2=i;i2<sizemax;i2++)
	{
		esp=0.;
		for (j2=j;j2<=i2-i+j;j2++)
		{
			V=(theta=GetTree().data(i2,j2).theta)+
					(jumpsize=GetTree().data(i2,j2).jumpsize);

			if ((i2==i)&&(j2==j))
				{esp += V;}
			else
				{esp += V * GetTree().data(i2,j2).proba;}
		}
		Epdef[i2]=esp;
	}

	return (esp);

}
