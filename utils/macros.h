#pragma once

#define NLATMAX           100

#define BOOST_CELL_THRESH    10
#define BOOST_NUMPART_THRESH 250
#define BOOST_BIN_REF        1

//what fraction of particles have to be sorted
//to switch from quicksort to heapsort
#define FRACTION_SORTED_REQD_TO_HEAP_SORT   0.6

#define ADD_DIFF_TIME(t0,t1)            ((t1.tv_sec - t0.tv_sec) + 1e-6*(t1.tv_usec - t0.tv_usec))
#define REALTIME_ELAPSED_NS(t0, t1)     ((t1.tv_sec - t0.tv_sec)*1000000000.0 + (t1.tv_nsec - t0.tv_nsec))

#define ALIGNMENT                32

#define STRINGIFY(x)   #x
#define STR(x) STRINGIFY(x)


#define MIN(X, Y)                        (( (X) < (Y)) ? (X):(Y))
#define MAX(X, Y)                        (( (X) > (Y)) ? (X):(Y))


#define ASSIGN_CELL_TIMINGS(thread_timings, nx1, nx2, timediff, tid, first_cellid, second_cellid) \
    {                                                                   \
        thread_timings->N1 = nx1;                                       \
        thread_timings->N2 = nx2;                                       \
        thread_timings->time_in_ns = timediff;                          \
        thread_timings->tid = tid;                                      \
        thread_timings->first_cellindex = first_cellid;                 \
        thread_timings->second_cellindex = second_cellid;               \
    }

/* Taken from http://stackoverflow.com/questions/19403233/compile-time-struct-size-check-error-out-if-odd
   which is in turn taken from the linux kernel */
/* #define BUILD_BUG_OR_ZERO(e) (sizeof(struct{ int:-!!(e);})) */
/* #define ENSURE_STRUCT_SIZE(e, size)  BUILD_BUG_OR_ZERO(sizeof(e) != size) */
/* However, the previous one gives me an unused-value warning and I do not want
   to turn that compiler warning off. Therefore, this version, which results in
   an unused local typedef warning is used. I turn off the corresponding warning
   in common.mk (-Wno-unused-local-typedefs) via CFLAGS
*/
#define BUILD_BUG_OR_ZERO(cond, msg) typedef volatile char assertion_on_##msg[( !!(cond) )*2-1 ]
#define ENSURE_STRUCT_SIZE(e, size)                 BUILD_BUG_OR_ZERO(sizeof(e) == size, sizeof_struct_config_options)

/* Macro Constants */
//Just to output some colors
#define ANSI_COLOR_RED     "\x1b[31m"
#define ANSI_COLOR_GREEN   "\x1b[32m"
#define ANSI_COLOR_YELLOW  "\x1b[33m"
#define ANSI_COLOR_CYAN    "\x1b[36m"
#define ANSI_COLOR_BLUE    "\x1b[34m"
#define ANSI_COLOR_MAGENTA "\x1b[35m"
#define ANSI_COLOR_RESET   "\x1b[0m"


#define PI_UNICODE    "\u03C0"
#define XI_UNICODE    "\u03BE"
/* #define PIMAX_UNICODE PI_UNICODE"\u2098""\u2090""\u2093" */
/* #define RP_UNICODE    "\u209A" */
#define PIMAX_UNICODE "pimax"
#define RP_UNICODE    "rp"
#define THETA_UNICODE "\u03B8"
#define OMEGA_UNICODE "\u03C9"
#define MU_UNICODE    "\u03BC"

#define PI_SAFE    "pi"
#define XI_SAFE    "xi"
#define PIMAX_SAFE "pimax"
#define RP_SAFE "rp"
#define THETA_SAFE "theta"
#define OMEGA      "omega"
#define MU_SAFE    "mu"


#ifdef USE_UNICODE
#define PI_CHAR PI_UNICODE
#define XI_CHAR XI_UNICODE
#define PIMAX_CHAR PIMAX_UNICODE
#define RP_CHAR  RP_UNICODE
#define MU_CHAR  MU_UNICODE
#define THETA_CHAR THETA_UNICODE
#define OMEGA_CHAR OMEGA_UNICODE
#define UNICODE_WARNING  "\n\
If you see unrendered characters, then your current terminal (likely) does not \n\
support UTF-8. Tun `xterm -u8` to get an UTF-8 compliant terminal (and hopefully) \n\
a font-set that supports Greek characters. For details on how to display unicode \n\
characters on xterm, see:\n\n\
  http://unix.stackexchange.com/questions/196152/xterm-not-displaying-unicode\n\n\
If none of the fixes work, disable ``USE_UNICODE`` option in ``common.mk`` in \n\
the ROOT DIRECTORY of ``Corrfunc`` and re-install the entire packge.\n"
#else
#define PI_CHAR PI_SAFE
#define XI_CHAR XI_SAFE
#define MU_CHAR MU_SAFE
#define PIMAX_CHAR PIMAX_SAFE
#define RP_CHAR    RP_SAFE
#define THETA_CHAR THETA_SAFE
#define OMEGA_CHAR OMEGA_SAFE
#define UNICODE_WARNING "\n"
#endif

/* Function-like macros */
#ifdef NDEBUG
#define XASSERT(EXP, ...)                                do{} while(0)
#else
#define XASSERT(EXP, ...)                                               \
     do { if (!(EXP)) {                                                 \
             fprintf(stderr,"Error in file: %s\tfunc: %s\tline: %d with expression `"#EXP"'\n", __FILE__, __FUNCTION__, __LINE__); \
             fprintf(stderr,__VA_ARGS__);                               \
             fprintf(stderr,ANSI_COLOR_BLUE "Hopefully, input validation. Otherwise, bug in code: please file an issue on GitHub: https://github.com/manodeep/Corrfunc/issues"ANSI_COLOR_RESET"\n"); \
             return EXIT_FAILURE;                                       \
         }                                                              \
     } while (0)
#endif

#ifdef NDEBUG
#define XPRINT(EXP, ...)                                do{} while(0)
#else
#define XPRINT(EXP, ...)                                               \
     do { if (!(EXP)) {                                                 \
             fprintf(stderr,"Error in file: %s\tfunc: %s\tline: %d with expression `"#EXP"'\n", __FILE__, __FUNCTION__, __LINE__); \
             fprintf(stderr,__VA_ARGS__);                               \
             fprintf(stderr,ANSI_COLOR_BLUE "Hopefully, input validation. Otherwise, bug in code: please file an issue on GitHub: https://github.com/manodeep/Corrfunc/issues"ANSI_COLOR_RESET"\n"); \
         }                                                              \
     } while (0)
#endif


#ifdef NDEBUG
#define XRETURN(EXP, VAL, ...)                                do{} while(0)
#else
#define XRETURN(EXP, VAL, ...)                                           \
     do { if (!(EXP)) {                                                 \
             fprintf(stderr,"Error in file: %s\tfunc: %s\tline: %d with expression `"#EXP"'\n", __FILE__, __FUNCTION__, __LINE__); \
             fprintf(stderr,__VA_ARGS__);                               \
             fprintf(stderr,ANSI_COLOR_BLUE "Hopefully, input validation. Otherwise, bug in code: please file an issue on GitHub: https://github.com/manodeep/Corrfunc/issues"ANSI_COLOR_RESET"\n"); \
             return VAL;                                                \
         }                                                              \
     } while (0)
#endif

#define SETUP_INTERRUPT_HANDLERS(handler_name)                          \
     const int interrupt_signals[] = {SIGTERM, SIGINT, SIGHUP};         \
     const size_t nsig = sizeof(interrupt_signals)/sizeof(interrupt_signals[0]); \
     typedef void (* sig_handlers)(int);                                \
     sig_handlers previous_handlers[nsig];                              \
     for(size_t i=0;i<nsig;i++) {                                       \
         int signo = interrupt_signals[i];                              \
         sig_handlers prev = signal(signo, handler_name);               \
         if (prev == SIG_ERR) {                                         \
             fprintf(stderr,"Can not handle signal = %d\n", signo);     \
         }                                                              \
         previous_handlers[i] = prev;                                   \
     }

#define RESET_INTERRUPT_HANDLERS()              \
     for(size_t i=0;i<nsig;i++) {                                       \
         int signo = interrupt_signals[i];                              \
         sig_handlers prev = previous_handlers[i];                      \
         if(prev == SIG_IGN || prev == SIG_ERR) continue;               \
         if(signal(signo, prev) == SIG_ERR) {                           \
             fprintf(stderr,"Could not reset signal handler to default for signal = %d\n", signo); \
         }                                                              \
     }
