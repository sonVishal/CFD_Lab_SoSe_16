#include "common.h"

//TODO: (TKS) Where should these swap functions be placed?
//          * Make a common function file
void swapComponentFields(t_component *c, int numComp){
    for (int i = 0; i < numComp; ++i) {
        double* swap = c[i].collideField;
        c[i].collideField = c[i].streamField;
        c[i].streamField = swap;
    }
}
