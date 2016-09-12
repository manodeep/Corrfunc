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
     
