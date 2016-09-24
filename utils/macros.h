#pragma once

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

#define PI_SAFE    "pi"
#define XI_SAFE    "xi"
#define PIMAX_SAFE "pimax"
#define RP_SAFE "rp"
#define THETA_SAFE "theta"
#define OMEGA      "omega"


#ifdef USE_UNICODE
#define PI_CHAR PI_UNICODE
#define XI_CHAR XI_UNICODE
#define PIMAX_CHAR PIMAX_UNICODE
#define RP_CHAR  RP_UNICODE
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
             fprintf(stderr,ANSI_COLOR_BLUE "Hopefully, input validation. Otherwise, bug in code: please email Manodeep Sinha <manodeep@gmail.com>"ANSI_COLOR_RESET"\n"); \
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
             fprintf(stderr,ANSI_COLOR_BLUE "Hopefully, input validation. Otherwise, bug in code: please email Manodeep Sinha <manodeep@gmail.com>"ANSI_COLOR_RESET"\n"); \
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
             fprintf(stderr,ANSI_COLOR_BLUE "Hopefully, input validation. Otherwise, bug in code: please email Manodeep Sinha <manodeep@gmail.com>"ANSI_COLOR_RESET"\n"); \
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
     
