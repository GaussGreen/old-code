#ifndef EDR_SRMCONSTANTS_HPP
#define EDR_SRMCONSTANTS_HPP

DRLIB_BEGIN_NAMESPACE

namespace SRMConstants {
    
    MCARLO_DLL extern const double SRM_TINY;   // 1E-10;
    MCARLO_DLL extern const double QCUTOFF;    // 1E-4 Normal model for |q|<QCUTOFF
    MCARLO_DLL extern const double EXP_CUTOFF; // 7.0E+2

    MCARLO_DLL extern const int MAX2QITER;     // 20;

    MCARLO_DLL extern const double MAX_DIFFUSION_STEP_SIZE; // = 1.0;
}

DRLIB_END_NAMESPACE

#endif // EDR_SRMCONSTANTS_HPP
