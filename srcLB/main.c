#ifndef _MAIN_C_
#define _MAIN_C_

#include "boundary.h"
#include "collision.h"
#include "initLB.h"
#include "streaming.h"
#include "visualLB.h"

//TODO: (TKS) Make LBDefinitions work
#include "LBDefinitions.h"

int main(int argc, char *argv[]){

    double *collideField=NULL;
    double *streamField=NULL;

    /* TODO: (DL) in the worksheet there is "flagField" and "flagfield" -
     * for the moment I think they are the same... */
    int *flagField=NULL;
    int xlength;
    double tau;
    double velocityWall[3];
    int timesteps;
    int timestepsPerPlotting;

    readParameters(&xlength, &tau, velocityWall, &timesteps, &timestepsPerPlotting,argc, &argv[1]);

    /*Initializing pointers*/
    size_t totalsize = (xlength+2)*(xlength+2)*(xlength+2);
    collideField  = (double *)  malloc(Q*totalsize * sizeof( double ));
    streamField   = (double *)  malloc(Q*totalsize * sizeof( double ));

    /* TODO: (DL) we could save memory and use 'short' instead of int */
    /* calloc: only required to set boundary values */
    flagField     = (int *)  calloc(totalsize, sizeof( int ));

    initialiseFields(collideField, streamField, flagField, xlength);

    for(int t = 0; t < timesteps; t++){
	    double *swap=NULL;
	    doStreaming(collideField,streamField,flagField,xlength);

	    swap = collideField;
	    collideField = streamField;
	    streamField = swap;

	    //doCollision(collideField,flagField,&tau,xlength);
	    //treatBoundary(collideField,flagField,velocityWall,xlength);

	    if (t%timestepsPerPlotting==0){
	        //writeVtkOutput(collideField,flagField,argv,t,xlength);
	    }
    }

    free(streamField);
    free(collideField);
    free(flagField);

    return 0;
}

#endif
