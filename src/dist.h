#ifndef DIST_H
#define DIST_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define MAX_NAME_LEN 64
#define MAX_MATRICES 100

typedef enum {
    DIST_BERNOULLI,
    DIST_BINOMIAL,
    DIST_POISSON,
    DIST_GEOMETRIC,
    DIST_HYPERGEOMETRIC,
    DIST_UNIFORM,
    DIST_EXPONENTIAL,
    DIST_NORMAL,
    DIST_GAMMA,
    DIST_BETA,
    DIST_CHISQ,
    DIST_STUDENT,
    DIST_FDIST,
    DIST_LOGNORMAL,
    DIST_WEIBULL
} DistType;

typedef struct {
    DistType type;
    double params[4];
    int int_params[3];
    char name[MAX_NAME_LEN];
} Distribution;

typedef struct {
    Distribution dists[MAX_MATRICES];
    int count;
} DistStore;

void dist_init(DistStore *store);
int dist_set(DistStore *store, const char *name, Distribution *dist);
Distribution *dist_get(DistStore *store, const char *name);
void dist_list(DistStore *store);
void dist_clear(DistStore *store);
void dist_print(Distribution *dist);

Distribution *bernoulli(double p);
Distribution *binomial(int n, double p);
Distribution *poisson(double lambda);
Distribution *geometric(double p);
Distribution *hypergeometric(int N, int K, int n);
Distribution *uniform(double a, double b);
Distribution *exponential(double lambda);
Distribution *normal(double mu, double var);
Distribution *gamma(double alpha, double beta);
Distribution *beta(double alpha, double beta);
Distribution *chisq(int df);
Distribution *student(int df);
Distribution *fdist(int df1, int df2);
Distribution *lognormal(double mu, double var);
Distribution *weibull(double lambda, double k);

double pdf(Distribution *dist, double x);
double cdf(Distribution *dist, double x);
double quantile(Distribution *dist, double p);
double mean(Distribution *dist);
double variance(Distribution *dist);
double stddev(Distribution *dist);

#endif