#ifndef STATS_H
#define STATS_H

#include <stdlib.h>
#include "dist.h"

typedef struct {
    double *data;
    int n;
} DataArray;

DataArray *data_create(int n);
void data_free(DataArray *d);
double data_mean(DataArray *d);
double data_variance(DataArray *d);
double data_std(DataArray *d);
double data_sum(DataArray *d);

double chisq_test(DataArray *observed, DataArray *expected);
double ttest_one(DataArray *data, double mu);
double ttest_two(DataArray *data1, DataArray *data2);
double ztest(DataArray *data, double mu, double sigma);

#endif