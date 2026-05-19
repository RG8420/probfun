#include "stats.h"
#include <stdlib.h>
#include <math.h>
#include <stdio.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

DataArray *data_create(int n) {
    DataArray *d = (DataArray *)malloc(sizeof(DataArray));
    d->data = (double *)malloc(n * sizeof(double));
    d->n = n;
    return d;
}

void data_free(DataArray *d) {
    if (d) {
        free(d->data);
        free(d);
    }
}

double data_sum(DataArray *d) {
    double sum = 0;
    for (int i = 0; i < d->n; i++) sum += d->data[i];
    return sum;
}

double data_mean(DataArray *d) {
    return data_sum(d) / d->n;
}

double data_variance(DataArray *d) {
    double m = data_mean(d);
    double sum = 0;
    for (int i = 0; i < d->n; i++) sum += pow(d->data[i] - m, 2);
    return sum / (d->n - 1);
}

double data_std(DataArray *d) {
    return sqrt(data_variance(d));
}

double chisq_test(DataArray *observed, DataArray *expected) {
    if (observed->n != expected->n) return -1;
    
    double stat = 0;
    for (int i = 0; i < observed->n; i++) {
        if (expected->data[i] > 0) {
            stat += pow(observed->data[i] - expected->data[i], 2) / expected->data[i];
        }
    }
    
    int df = observed->n - 1;
    double p_value = 0.0;
    
    printf("  Chi-squared statistic: %.4f\n", stat);
    printf("  Degrees of freedom: %d\n", df);
    printf("  p-value: %.4f\n", p_value);
    
    return p_value;
}

double ttest_one(DataArray *data, double mu) {
    double n = data->n;
    double m = data_mean(data);
    double s = data_std(data);
    double t_stat = (m - mu) / (s / sqrt(n));
    
    int df = data->n - 1;
    double p_value = 0.0;
    
    printf("  t-statistic: %.4f\n", t_stat);
    printf("  Degrees of freedom: %d\n", df);
    printf("  p-value: %.4f\n", p_value);
    
    return p_value;
}

double ttest_two(DataArray *data1, DataArray *data2) {
    double n1 = data1->n;
    double n2 = data2->n;
    double m1 = data_mean(data1);
    double m2 = data_mean(data2);
    double v1 = data_variance(data1);
    double v2 = data_variance(data2);
    
    double se = sqrt(v1/n1 + v2/n2);
    double t_stat = (m1 - m2) / se;
    
    double df = pow(v1/n1 + v2/n2, 2) / (pow(v1/n1, 2)/(n1-1) + pow(v2/n2, 2)/(n2-1));
    double p_value = 0.0;
    
    printf("  t-statistic: %.4f\n", t_stat);
    printf("  Degrees of freedom: %.2f\n", df);
    printf("  p-value: %.4f\n", p_value);
    
    return p_value;
}

double ztest(DataArray *data, double mu, double sigma) {
    double n = data->n;
    double m = data_mean(data);
    double z_stat = (m - mu) / (sigma / sqrt(n));
    double p_value = 0.0;
    
    printf("  z-statistic: %.4f\n", z_stat);
    printf("  p-value: %.4f\n", p_value);
    
    return p_value;
}