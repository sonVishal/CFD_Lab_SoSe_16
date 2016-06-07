
#include "helperMPI.h"

//TODO: (TKS) 
//Reduce bufferSize[6]--> bufferSize[3]?
//          * Can do two iterations of the loop at the same time.
//          * Or complete loop unroll

void communicate(double** sendBuffer, double**readBuffer, double* collideField, int* xlength, int* bufferSize){
    //Run extract, swap, inject for all sides.
    //
    for (int direction = LEFT; direction < BACK; ++direction) {
        extract(sendBuffer, collideField, xlength, bufferSize, direction);
        swap(sendBuffer, readBuffer, bufferSize, direction);
        inject(readBuffer, collideField, xlength, bufferSize, direction);
    }

}

//Copy distributions needed to sendbuffer.
void extract( double** sendBuffer, double* collideField, int* xlength, int* bufferSize, int direction){
}

//Send distributions and wait to receive.
void swap(double** sendBuffer, double**readBuffer, int* bufferSize, int direction){

}

//Copy read buffer to ghost layer
void inject(double** readBuffer, double* collideField, int* xlength, int* bufferSize, int direction){
}
