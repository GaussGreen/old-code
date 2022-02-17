
#include "ICMKernel\mod\icm_meta_model.h"


void ICM_Meta_Model::Init(void)
{
	SetName(ICM_META_MODEL);

	itsPricerType.clear();
	itsModels = NULL;  
	itsNbModel = 0;
}


void ICM_Meta_Model::Set(ARM_Model** models,
						 const int& nbofmodels,
						 const vector<ARM_CLASS_NAME>& PricerType)
{

	int i =0;

	if (itsModels)
	{
		delete[] itsModels;
		itsModels = NULL;
	}

	itsModels = new ARM_Model*[nbofmodels];
	itsNbModel = nbofmodels;

	for (i=0;i<itsNbModel;i++)
		itsModels[i] = (models)[i];

	SetStartDate((models)[0]->GetStartDate());

	SetPricerType(PricerType);

}

// ----------------------------
//	Copy of members data
// ----------------------------
void ICM_Meta_Model::BitwiseCopy(const ARM_Object* src)
{
	int i = 0;
    ICM_Meta_Model* pf = (ICM_Meta_Model*) src;

	if (pf->itsModels)
	{
		if (itsModels)
			delete[] itsModels;
		itsModels = NULL;

		itsModels = new ARM_Model*[pf->itsNbModel];

		for (i = 0; i<pf->itsNbModel; i++)
				itsModels[i] = (pf->itsModels)[i];

	}

	SetPricerType(pf->itsPricerType);


	itsNbModel = pf->itsNbModel;

}

// -------------
//	Copy Method 
// -------------
void ICM_Meta_Model::Copy(const ARM_Object* src)
{
    ARM_Model::Copy(src);

    BitwiseCopy(src);
}

// --------------
//	Clone Method
// --------------
ARM_Object* ICM_Meta_Model::Clone(void)
{
ICM_Meta_Model* theClone = new ICM_Meta_Model();

theClone->Copy(this);
 
return(theClone);
}



