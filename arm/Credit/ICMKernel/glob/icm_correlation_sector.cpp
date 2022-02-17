/* -----------------------------------------------------------------
   -----------------------------------------------------------------

	FILE:		ICM_Correlation_Sector.CPP
	PROJECT:	MOD
	
	DESCRIPTION:	Containor for a Sector Correlation

  -----------------------------------------------------------------

 	CREATION:	February, 2006

	LAST MODIF:	February, 2006
  -----------------------------------------------------------------
   
	ICM Library

		version 1.0
		developped by Laurent Jacquel

  -----------------------------------------------------------------
  ----------------------------------------------------------------- */


#include "ICMKernel\glob\icm_correlation_sector.h"


void	ICM_Correlation_Sector::View(char* id, FILE* ficOut)
{
	FILE* fOut;
	char  fOutName[200];
	int i = 0;

	if (ficOut == NULL)
	{
		ARM_GetViewFile(ARM_VIEW_FILE, id, fOutName);

		(void) unlink(fOutName);

		fOut = fopen(fOutName, "w"); 
	}
	else
	{
		fOut = ficOut;
	} 

	fprintf(fOut, "\n ======> Sectorial Correlation:\n\n");
	
	// PARAMETERS
	fprintf(fOut, "Nb Sectors:\t%u\n", its_Nb_Sectors);
	fprintf(fOut, "Nb Names:\t%u\n", GetLabels().size());

	fprintf(fOut, "Correlation Type:\t");
	switch (its_CorrelationType)
	{
		case TFCT_FULL:
			fprintf(fOut, "\t%s", "Diff Inter - Diff Intra\n");
			fprintf(fOut, "\n\t\tBetas\t\tLambdas\n");
			for (i=0; i<its_Nb_Sectors; i++)
				fprintf(fOut, "\t\t%f\t\t%f\n", its_Betas[i], its_Lambdas[i]);

			break;

		case TFCT_SAME_INTER_DIFF_INTRA:
			fprintf(fOut, "\t%s", "Same Inter - Diff Intra\n");
			fprintf(fOut, "Inter Sector Correlation:\t%f\n", its_Betas[0]*its_Betas[0]*its_Lambdas[0]*its_Lambdas[0]);
	
			fprintf(fOut, "\n\t\tBetas\n");
			for (i=0; i<its_Nb_Sectors; i++)
				fprintf(fOut, "\t\t%f\n", its_Betas[i]);
			break;

		case TFCT_SAME_INTER_SAME_INTRA:
			fprintf(fOut, "\t%s", "Same Inter - Same Intra\n");
			fprintf(fOut, "Inter Sector Correlation:\t%f\n", its_Betas[0]*its_Betas[0]*its_Lambdas[0]*its_Lambdas[0]);
			fprintf(fOut, "Intra Sector Correlation:\t%f\n", its_Betas[0]*its_Betas[0]);
			break;

		fprintf(fOut, "\n");

	}

	if (!GetLabels().empty())
	{
		fprintf(fOut, "\nLabels\t\tSector Membership\t\n");

		for (i=0; i<GetLabels().size(); i++)
		{
			fprintf(fOut, "%s\t\tSector Id: %u\t\t\n", GetLabel(i).c_str(), its_Sector_Membership[i]);
		}
		fprintf(fOut, "\n");
		fprintf(fOut, "\nSector\t\tBetas\t\tLambdas\n");

		if (its_Betas.size() == its_Nb_Sectors)
		{
			if ( its_Lambdas.size() == its_Nb_Sectors)
			{
				for (i=0; i<its_Nb_Sectors; i++)
					fprintf(fOut, "Sector Id: %u\t%lf\t%lf\n", i, its_Betas[i], its_Lambdas[i]);
			}
			else 
			{
				for (i=0; i<its_Nb_Sectors; i++)
					fprintf(fOut, "Sector Id: %u\t%lf\t\t\n", i, its_Betas[i]);
			}
		}
		if (its_Betas_Down.size() == its_Nb_Sectors)
		{
			if ( its_Lambdas_Down.size() == its_Nb_Sectors)
			{
				for (i=0; i<its_Nb_Sectors; i++)
					fprintf(fOut, "Sector Id: %u\t%lf\t%lf\n", i, its_Betas_Down[i], its_Lambdas_Down[i]);
			}
			else 
			{
				for (i=0; i<its_Nb_Sectors; i++)
					fprintf(fOut, "Sector Id: %u\t%lf\t\t\n", i, its_Betas_Down[i]);
			}
		}

	}

	if (ficOut == NULL)
	{
		fclose(fOut);
	}
}