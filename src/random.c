#include "random.h"
#include "dist.h"
#include <stdlib.h>
#include <time.h>
#include <math.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

static unsigned int current_seed = 0;

void seed_random(unsigned int seed) {
    current_seed = seed;
    srand(seed);
}

double rand_uniform(void) {
    return (double)rand() / RAND_MAX;
}

double rand_normal(double mu, double var) {
    double u1 = rand_uniform();
    double u2 = rand_uniform();
    double z = sqrt(-2.0 * log(u1)) * cos(2.0 * M_PI * u2);
    return mu + sqrt(var) * z;
}

double rand_exponential(double lambda) {
    return -log(1 - rand_uniform()) / lambda;
}

static double gamma_variate(double alpha, double beta) {
    if (alpha < 1.0) {
        return gamma_variate(1.0 + alpha, beta) * pow(rand_uniform(), 1.0 / alpha);
    }
    double d = alpha - 1.0 / 3.0;
    double c = 1.0 / sqrt(9.0 * d);
    while (1) {
        double x = rand_normal(0, 1);
        double v = 1.0 + c * x;
        while (v <= 0) {
            x = rand_normal(0, 1);
            v = 1.0 + c * x;
        }
        v = v * v * v;
        double u = rand_uniform();
        if (u < 1 - 0.0331 * x * x * x * x) return d * v * beta;
        if (log(u) < 0.5 * x * x + d * (1 - v + log(v))) return d * v * beta;
    }
}

double rand_gamma(double alpha, double beta) {
    return gamma_variate(alpha, beta) / beta;
}

static double beta_variate(double alpha, double beta) {
    double x = gamma_variate(alpha, 1.0);
    double y = gamma_variate(beta, 1.0);
    return x / (x + y);
}

double rand_beta(double alpha, double beta) {
    return beta_variate(alpha, beta);
}

double rand_chisq(int df) {
    return 2.0 * gamma_variate(df / 2.0, 0.5);
}

double rand_student(int df) {
    double z = rand_normal(0, 1);
    double v = rand_chisq(df) / df;
    return z / sqrt(v);
}

double rand_fdist(int df1, int df2) {
    double x = rand_chisq(df1) / df1;
    double y = rand_chisq(df2) / df2;
    return x / y;
}

double rand_lognormal(double mu, double var) {
    return exp(rand_normal(mu, var));
}

double rand_weibull(double lambda, double k) {
    return lambda * pow(-log(1 - rand_uniform()), 1.0 / k);
}

double rand_bernoulli(double p) {
    return (rand_uniform() < p) ? 1.0 : 0.0;
}

int rand_binomial(int n, double p) {
    int count = 0;
    for (int i = 0; i < n; i++) {
        if (rand_uniform() < p) count++;
    }
    return count;
}

int rand_poisson(double lambda) {
    double L = exp(-lambda);
    int k = 0;
    double p = 1.0;
    do {
        k++;
        p *= rand_uniform();
    } while (p > L && k < 1000);
    return k - 1;
}

int rand_geometric(double p) {
    return (int)ceil(log(1 - rand_uniform()) / log(1 - p));
}

int rand_hypergeometric(int N, int K, int n) {
    int successes = 0;
    for (int i = 0; i < n; i++) {
        if ((double)K / (N - i) > rand_uniform()) {
            successes++;
            K--;
        }
        N--;
    }
    return successes;
}

double *rand_samples(Distribution *dist, int n) {
    double *samples = (double *)malloc(n * sizeof(double));
    for (int i = 0; i < n; i++) {
        switch (dist->type) {
            case DIST_BERNOULLI:
                samples[i] = rand_bernoulli(dist->params[0]);
                break;
            case DIST_BINOMIAL:
                samples[i] = rand_binomial(dist->int_params[0], dist->params[0]);
                break;
            case DIST_POISSON:
                samples[i] = rand_poisson(dist->params[0]);
                break;
            case DIST_GEOMETRIC:
                samples[i] = rand_geometric(dist->params[0]);
                break;
            case DIST_HYPERGEOMETRIC:
                samples[i] = rand_hypergeometric(dist->int_params[0], dist->int_params[1], dist->int_params[2]);
                break;
            case DIST_UNIFORM:
                samples[i] = dist->params[0] + rand_uniform() * (dist->params[1] - dist->params[0]);
                break;
            case DIST_EXPONENTIAL:
                samples[i] = rand_exponential(dist->params[0]);
                break;
            case DIST_NORMAL:
                samples[i] = rand_normal(dist->params[0], dist->params[1]);
                break;
            case DIST_GAMMA:
                samples[i] = rand_gamma(dist->params[0], dist->params[1]);
                break;
            case DIST_BETA:
                samples[i] = rand_beta(dist->params[0], dist->params[1]);
                break;
            case DIST_CHISQ:
                samples[i] = rand_chisq(dist->int_params[0]);
                break;
            case DIST_STUDENT:
                samples[i] = rand_student(dist->int_params[0]);
                break;
            case DIST_FDIST:
                samples[i] = rand_fdist(dist->int_params[0], dist->int_params[1]);
                break;
            case DIST_LOGNORMAL:
                samples[i] = rand_lognormal(dist->params[0], dist->params[1]);
                break;
            case DIST_WEIBULL:
                samples[i] = rand_weibull(dist->params[0], dist->params[1]);
                break;
            default:
                samples[i] = 0;
        }
    }
    return samples;
}