
extern int OptionPricing();
extern int CBPricing();
extern int GeneralizedPEPSLikePricing();
extern int AttachedWarrantCBPricing();

extern int CalibrationVolFlatHRTimeComponent();

int main()
{
  OptionPricing();
  
  CBPricing();

  GeneralizedPEPSLikePricing();

  AttachedWarrantCBPricing();

  CalibrationVolFlatHRTimeComponent();

  return 0;
}
