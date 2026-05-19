#ifndef RANDOM_H
#define RANDOM_H

#include <stdlib.h>
#include <time.h>
#include "dist.h"

void seed_random(unsigned int seed);
double rand_uniform(void);
double rand_normal(double mu, double var);
double rand_exponential(double lambda);
double rand_gamma(double alpha, double beta);
double rand_beta(double alpha, double beta);
double rand_chisq(int df);
double rand_student(int df);
double rand_fdist(int df1, int df2);
double rand_lognormal(double mu, double var);
double rand_weibull(double lambda, double k);
double rand_bernoulli(double p);
int rand_binomial(int n, double p);
int rand_poisson(double lambda);
int rand_geometric(double p);
int rand_hypergeometric(int N, int K, int n);

double *rand_samples(Distribution *dist, int n);

#endif