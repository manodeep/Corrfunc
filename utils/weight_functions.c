#include "defs.h"

typedef enum {
  NONE=-42, /* default */
  PAIR_PRODUCT=0,
  NUM_WEIGHT_TYPE 
} weight_method_t; // type of weighting to apply

//////////////////////////////////
// Weighting functions
//////////////////////////////////

/*
 * The pair weight is the product of the particle weights
 */
static inline DOUBLE pair_product(const pair_struct *pair){
    return pair->weights0[0]*pair->weights1[0];
}


//////////////////////////////////
// Utility functions
//////////////////////////////////

/* Maps a name to weighting method
   `method` will be set on return.
 */
static inline int get_weight_method_by_name(const char *name, weight_method_t *method){
    if(name == NULL || strcmp(name, "")){
        *method = NONE;
        return EXIT_SUCCESS;
    }
    if(strcmp(name, "pair_product") == 0 || strcmp(name, "p") == 0){
        *method = PAIR_PRODUCT;
        return EXIT_SUCCESS;
    }
        
    return EXIT_FAILURE;
}

/* Gives the number of weight arrays required by the given weighting method
 */
static inline int get_num_weights_by_method(const weight_method_t method){
    switch(method){
        case PAIR_PRODUCT:
            return 1;
        default:
        case NONE:
            return 0;
    }
}

/* Gives a pointer to the weight function for the given weighting method
 */
static inline weight_func_t get_weight_func_by_method(const weight_method_t method){
    switch(method){
        case PAIR_PRODUCT:
            return &pair_product;
        default:
        case NONE:
            return NULL;
    }
}
