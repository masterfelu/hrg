#include "MiscellaneousFunctions.h"

#include <cmath>

double moment(int *x, float *w, int length, int order){
	double sum = 0;
	for (int i = 0; i < length; ++i)
	{
		sum += pow(x[i],order)*w[i];
	}
	return sum;
}

double centralMoment(int *x, float *w, int length, int order){
	double mean = moment(x,w,length,1);
	double sum = 0;
	for (int i = 0; i < length; ++i)
	{
		sum += pow((x[i]-mean),order)*w[i];
	}
	return sum;
}

double correlation(int *x1, int *x2, float *w, int length){
	double mean1 = moment(x1,w,length,1);
	double mean2 = moment(x2,w,length,1);
	double sum = 0;
	for (int i = 0; i < length; ++i)
	{
		sum += (x1[i]-mean1)*(x2[i]-mean2)*w[i];
	}
	return sum;
}
