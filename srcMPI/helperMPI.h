
#ifndef _HELPERMPI_H_
#define _HELPERMPI_H_


typedef struct{
    int startInner, endInner;
    int startOuter, endOuter;
    int fixedValue;

} t_iterPara;

void communicate(double** sendBuffer, double**readBuffer, double* collideField, const t_procData *procData);
void extract( double** sendBuffer, double* collideField, const t_iterPara *iterPara, const t_procData *procData,
              int direction, int* index);
void swap(double** sendBuffer, double** readBuffer, const t_procData *procData, int direction);
void inject(double** readBuffer, double* collideField, t_iterPara *iterPara, const t_procData *procData, 
            int direction, int *index);

void p_setCommIterationParameters(t_iterPara *iterPara, const t_procData *procData, const int direction);
void p_assignIndices(int direction, int *index);

#endif
