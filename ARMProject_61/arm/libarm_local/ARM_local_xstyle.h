#ifndef ARM_LOCAL_XSTYLE_H
#define ARM_LOCAL_XSTYLE_H



extern long ARMLOCAL_BERMUDANXSTYLE(VECTOR<double>& xDates,
                                    VECTOR<double>& expiryDates,
									ARM_result& result,
									long objId = -1);

extern long ARMLOCAL_EUROPEANXSTYLE(double xdate,
									ARM_result& result,
									long objId = -1);

extern long ARMLOCAL_AMERICANXSTYLE(double xStartDate,
									double xEndDate,
									ARM_result& result,
									long objId = -1);

#endif /* ARM_LOCAL_XSTYLE_H */