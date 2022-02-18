///* ==========================================================================
//   declarations for all functions
//   ========================================================================== */
//int proba_max_addin(int argc, UF_ARG* argv, UF_ARG* ret);
//
//int proba_min_addin(int argc, UF_ARG* argv, UF_ARG* ret);
//
//int proba_max_explo_addin(int argc, UF_ARG* argv, UF_ARG* ret);
//
//int proba_min_explo_addin(int argc, UF_ARG* argv, UF_ARG* ret);
//
//int proba_gauss_addin(int argc, UF_ARG* argv, UF_ARG* ret);
//
//int proba_gauss_explo_addin(int argc, UF_ARG* argv, UF_ARG* ret);
//
///**int 	proba_four_addin(
//                        int 	argc,
//                        UF_ARG 	*argv,
//                        UF_ARG 	*ret
//                        );
//
//int 	probaf_yesno_addin(
//                        int 	argc,
//                        UF_ARG 	*argv,
//                        UF_ARG 	*ret
//                        );
//
//int 	probaf_yesno_explo_addin(
//                        int 	argc,
//                        UF_ARG 	*argv,
//                        UF_ARG 	*ret
//                        ); **/
//
//int probag_yesno_addin(int argc, UF_ARG* argv, UF_ARG* ret);
//
//int probag_yesno_explo_addin(int argc, UF_ARG* argv, UF_ARG* ret);
//
///* ==========================================================================
//   Library info structure
//   ========================================================================== */
//
//static UF_INFO lib_info = {"PROBA_TOOLS", "V1.00", __DATE__};
//
///* ==========================================================================
//   Function descriptor table
//   ========================================================================== */
//
//static UF_DESC desc_info[] = {
//    "proba_max",
//    "DDDDDDI",
//    "D",
//    "spot,barrier,vol,mat,rd,rf,dom/for",
//    proba_max_addin,
//
//    "proba_min",
//    "DDDDDDI",
//    "D",
//    "spot,barrier,vol,mat,rd,rf,dom/for",
//    proba_min_addin,
//
//    "proba_max_explo",
//    "DDDDDDI",
//    "D",
//    "spot,barrier,vol,mat,rd,rf,dom/for",
//    proba_max_explo_addin,
//
//    "proba_min_explo",
//    "DDDDDDI",
//    "D",
//    "spot,barrier,vol,mat,rd,rf,dom/for",
//    proba_min_explo_addin,
//
//    "proba_range_gauss",
//    "DDDDDDDII",
//    "D",
//    "spot,bar_do,bar_up,vol,mat,rd,rf,nb_term,dom/for",
//    proba_gauss_addin,
//
//    "proba_range_gauss_explo",
//    "DDDDDDDII",
//    "D",
//    "spot,bar_do,bar_up,vol,mat,rd,rf,nb_term,dom/for",
//    proba_gauss_explo_addin,
//
//    /***	"proba_range_four",
//                    "DDDDDDDII",
//                    "D",
//                    "spot,bar_do,bar_up,vol,mat,rd,rf,nb_term,dom/for",
//            proba_four_addin,
//
//            "probaf_yesno",
//                    "DDDDDDDII",
//                    "D",
//                    "spot,bar_yes,bar_no,vol,mat,rd,rf,nb_term,dom/for",
//            probaf_yesno_addin,
//
//            "probaf_yesno_explo",
//                    "DDDDDDDI",
//                    "D",
//                    "spot,bar_yes,bar_no,vol,mat,rd,rf,nb_term",
//            probaf_yesno_explo_addin,  ****/
//
//    "probag_yesno",
//    "DDDDDDDII",
//    "D",
//    "spot,bar_yes,bar_no,vol,mat,rd,rf,nb_term,dom/for",
//    probag_yesno_addin,
//
//    "probag_yesno_explo",
//    "DDDDDDDII",
//    "D",
//    "spot,bar_yes,bar_no,vol,mat,rd,rf,nb_term,dom/for",
//    probag_yesno_explo_addin,
//
//    "",
//    "",
//    "",
//    "",
//    0
//
//};
