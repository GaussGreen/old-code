#include "imsl_inc.h"

#if defined( _MSC_VER )
#pragma warning( once : 4101 4102 4244 4305 )
#endif

/*
#define Mac2ms
*/

#include <signal.h>

typedef void    (*Signal)();
/* whether to use signals or not (it's expensive) */
#define xUSE_SIGNALS

#ifdef USE_SIGNALS
static void     PROTO(l_handle_signal,(int));
#if defined(COMPUTER_VAX) || defined(COMPUTER_ALFAC_IEEE)
#define SIGABRT SIGIOT
#endif


#define   N_STACK   20

#if defined(IMSL_MACHINE_80X86) || defined(IMSL_MACHINE_NT)
#define N_CATCH		2
static Mint lv_catch[N_CATCH] = {
  SIGILL, SIGSEGV
};
#else
#define N_CATCH		3
static Mint lv_catch[N_CATCH] = {
#if defined(COMPUTER_HP93C) || defined(COMPUTER_HP98C) || defined(COMPUTER_HP97C)
  SIGILL, SIGSEGV, _SIGBUS
#elif defined(Mac2ms)
  SIGILL, SIGSEGV
#else
  SIGILL, SIGSEGV, SIGBUS
#endif
};
#endif

static Signal   lv_signal_stack[N_STACK][N_CATCH];
static Mint     lv_sp = 0;
#endif

void imsl_signal()
{
  return;
}

#ifdef ANSI
void imsl_set_signal(Mint on)
#else
void imsl_set_signal(on)
    Mint     on;
#endif
{
#ifdef USE_SIGNALS
    Mint     *psig;

    if (on) {
        for (psig=lv_catch;  psig < lv_catch+N_CATCH;  psig++)
            lv_signal_stack[lv_sp][*psig]  = signal(*psig,  l_handle_signal);
        lv_sp++;
        imsl_error_param.signal_set = 0;
    } else {
        lv_sp--;
        for (psig=lv_catch;  psig < lv_catch+N_CATCH;  psig++)
            signal(*psig, lv_signal_stack[lv_sp][*psig]);
     }
#endif
}


#ifdef USE_SIGNALS
#ifdef ANSI
static void l_handle_signal(int sig)
#else
static void l_handle_signal(sig)
    Mint	sig;
#endif
{
    Mint        error_type;

    imsl_error_param.signal_set = 1;
    imsl_e1sti (1, (Mint)sig);
    error_type = (Mint)IMSL_TERMINAL;
    switch (sig) {
        case SIGINT:    imsl_ermes(error_type, IMSL_SIGNAL_INT);	break;
        case SIGILL:    imsl_ermes(error_type, IMSL_SIGNAL_ILL);	break;
        case SIGABRT:   imsl_ermes(error_type, IMSL_SIGNAL_ABRT);	break;
        case SIGFPE:    imsl_ermes(error_type, IMSL_SIGNAL_FPE);	break;
#if !(defined(IMSL_MACHINE_80X86) || defined(IMSL_MACHINE_NT))
#if defined(COMPUTER_HP93C ) || defined(COMPUTER_HP98C) || defined(COMPUTER_HP97C)
        case _SIGBUS:    imsl_ermes(error_type, IMSL_SIGNAL_BUS);	break;
#elif defined(Mac2ms)
#else
        case SIGBUS:    imsl_ermes(error_type, IMSL_SIGNAL_BUS);	break;
#endif
#endif
        case SIGSEGV:   imsl_ermes(error_type, IMSL_SIGNAL_SGEV);	break;
    }
    if (lv_sp > 0) {
	longjmp(imsl_error_param.jmpbuf_environ[--imsl_got_environ], 1);
    }
}
#endif

