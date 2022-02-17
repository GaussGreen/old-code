
#include <ARM\libarm_local\firstToBeIncluded.h>

#include <ICMKernel\crv\icm_constant_piecewise.h>
#include <ICMKernel\glob\icm_corrmatrix.h>


#include <atlbase.h>
#include <ARM\libarm_frometk\XMLTools.h>
#include <ARM\libarm_frometk\arm_local_parsexml_util.h>
#include "VariantTools.h"

#include <ICMKernel\inst\icm_cds.h>
#include <ICMKernel\inst\icm_mez.h>
#include <ICMKernel\inst\icm_nthtd.h>

#include <ARM\libarm_frometk\arm_local_parsexml.h>

#include <ARM\libarm_frometk\arm_local_etoolkit_for_icm.h>
#include <ARM\libarm_frometk\arm_local_parsexml_for_icm.h>


#include <ARM\libarm_frometk\PaserManagerUtilities.h>

#include <ARM\libarm_frometk\arm_local_paesexml_calypso.h> 

#include <libCCatl\CCatl.h>


extern ICM_Cds* ARMLOCAL_ParseCalypsoCDS(const string& chaineXML ,
										 const ARM_Date& date,
										 const string modeltype,
										 string& BookName,
										 string& calypsoId);

extern ICM_Nthtd* ARMLOCAL_ParseCalypsoNTD(const string& chaineXML,
										   const ARM_Date& date,
										   const string modeltype,
										   string& BookName,
										   string& calypsoId);
extern ICM_Mez* ARMLOCAL_ParseCalypsoCDO(const string& chaineXML,
										   const ARM_Date& date,
										   const string modeltype,
										   string& BookName,
										   string& calypsoId);
extern void ConversionAmounts(int nbnames ,
							  ARM_Vector* init_notionals ,
							  double sub_init,
							  double end_init,
							  double& sub_out,
							  double& tranche_out,
							  ARM_Vector* out_notionals);



