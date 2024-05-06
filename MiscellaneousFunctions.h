#ifndef MISCELLANEOUS_FUNCTIONS_H
#define MISCELLANEOUS_FUNCTIONS_H

double moment(int *x, float *w, int length, int order);
double centralMoment(int *x, float *w, int length, int order);
double correlation(int *x1, int *x2, float *w, int length);

#endif