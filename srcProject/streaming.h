#ifndef _STREAMING_H_
#define _STREAMING_H_
#include <stdio.h>
#include <assert.h>
#include "LBDefinitions.h"

/** carries out the streaming step and writes the respective distribution functions from
 *  collideField to streamField.
 */
void doStreaming(double *collideField, double *streamField,int *flagField,int *xlength);

#endif
