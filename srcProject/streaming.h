#ifndef _STREAMING_H_
#define _STREAMING_H_
#include <stdio.h>
#include <assert.h>
#include "LBDefinitions.h"

//Wrapper around doStreaming to stream numComp components
void streamComponents(t_component* c, int numComp, int *flagField, int *xlength);

/** carries out the streaming step and writes the respective distribution functions from
 *  collideField to streamField.
 */
void doStreaming(double *streamField, double *collideField, int *flagField,int *xlength);

#endif
