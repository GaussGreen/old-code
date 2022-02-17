#include "edginc/config.hpp"
#include "edginc/VolTools.hpp"
#include "edginc/Addin.hpp"
#include "edginc/SwapMaturityVolRequest.hpp"
#include "edginc/VolProcessedBS.hpp"


DRLIB_BEGIN_NAMESPACE

// Excel interface.
void VolTools::load(CClassSP& clazz) {
    clazz->setPublic();                         // make visible to EAS/spreadsheet
    clazz->setDescription("My description");
    REGISTER(VolTools, clazz);
    SUPERCLASS(CObject);                        // Base class
    FIELD(selectedExpiries,"Expiries");
    FIELD(selectedTenors,"Tenors");
    FIELD(irvol,"Swap vol object");
    EMPTY_SHELL_METHOD(defaultConstructor);

    Addin::registerObjectMethod(
            "GetSwapVol", Addin::UTILITIES,
            "Help message",
            false, Addin::expandSimple,
			&VolTools::getSwapVol);    

}

CClassConstSP const VolTools::TYPE = CClass::registerClassLoadMethod("VolTools", typeid(VolTools), load);

// constructor
VolTools::VolTools(ExpiryArraySP Exp, ExpiryArraySP Tenor, IRVolSP ir, CClassConstSP const &type):
CObject(type),
selectedExpiries(Exp),
selectedTenors(Tenor),
irvol(ir)
{
	validatePop2Object();
}

void VolTools::validatePop2Object() 
{
}

// Populate the market swaption vol grid. This is a public function.
CDoubleMatrixSP VolTools::getSwapVol()
{
    try
    {
        int row_size = selectedExpiries->size(); 
        int col_size = selectedTenors->size(); 
		DateTime baseDate = irvol->getBaseDate();

        CDoubleMatrixSP swapVol(new CDoubleMatrix(col_size, row_size));

        for(int i = 0; i < row_size; i++)
		{
			ExpirySP Exp = (*selectedExpiries)[i];
			DateTime expiryDate = Exp->toDate(baseDate);
            for(int j = 0; j < col_size; j++)
			{
				ExpirySP Ten = (*selectedTenors)[j];
				SwapMaturityVolRequest	volReqSwapMaturity(Ten.get());
				CVolProcessedSP			volProcessed2(irvol->getProcessedVol(&volReqSwapMaturity, 0));
				CVolProcessedBSSP		ptrVolProcessedBS(CVolProcessedBSSP::dynamicCast(volProcessed2));
				DateTimeArray			dtExpiryDates;
				dtExpiryDates.push_back(baseDate);
				dtExpiryDates.push_back(expiryDate);
				CDoubleArray            vols(dtExpiryDates.size() - 1);
		        ptrVolProcessedBS->CalcVol(dtExpiryDates, CVolProcessedBS::forward, vols);
                (*swapVol)[j][i] = vols[0];
			}
		}

        return swapVol;
    }
    catch(exception& e) { throw ModelException(e, "VolTools::getSwapVol"); }
}


// Finally we need this function to force the compiler to register our class
bool VolToolsLoad(void) {
    return (VolTools::TYPE != 0);
}

DRLIB_END_NAMESPACE

