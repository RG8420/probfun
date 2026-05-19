#include "dist.h"
#include "random.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

static double beta_cdf(double x, double a, double b);
static double incomplete_beta(double a, double b, double x);
static double beta_cf(double a, double b, double x);
static double log_gamma(double x);

void dist_init(DistStore *store) {
    store->count = 0;
}

int dist_set(DistStore *store, const char *name, Distribution *dist) {
    for (int i = 0; i < store->count; i++) {
        if (strcmp(store->dists[i].name, name) == 0) {
            store->dists[i] = *dist;
            return 1;
        }
    }
    if (store->count >= MAX_MATRICES) return 0;
    store->dists[store->count] = *dist;
    strncpy(store->dists[store->count].name, name, MAX_NAME_LEN - 1);
    store->dists[store->count].name[MAX_NAME_LEN - 1] = '\0';
    store->count++;
    return 1;
}

Distribution *dist_get(DistStore *store, const char *name) {
    for (int i = 0; i < store->count; i++) {
        if (strcmp(store->dists[i].name, name) == 0) {
            return &store->dists[i];
        }
    }
    return NULL;
}

void dist_list(DistStore *store) {
    if (store->count == 0) {
        printf("  No distributions stored.\n");
        return;
    }
    for (int i = 0; i < store->count; i++) {
        printf("  %s: ", store->dists[i].name);
        dist_print(&store->dists[i]);
    }
}

void dist_clear(DistStore *store) {
    store->count = 0;
}

void dist_print(Distribution *dist) {
    switch (dist->type) {
        case DIST_BERNOULLI:
            printf("Bernoulli(p=%.2f)\n", dist->params[0]);
            break;
        case DIST_BINOMIAL:
            printf("Binomial(n=%d, p=%.2f)\n", dist->int_params[0], dist->params[0]);
            break;
        case DIST_POISSON:
            printf("Poisson(lambda=%.2f)\n", dist->params[0]);
            break;
        case DIST_GEOMETRIC:
            printf("Geometric(p=%.2f)\n", dist->params[0]);
            break;
        case DIST_HYPERGEOMETRIC:
            printf("Hypergeometric(N=%d, K=%d, n=%d)\n", 
                   dist->int_params[0], dist->int_params[1], dist->int_params[2]);
            break;
        case DIST_UNIFORM:
            printf("Uniform(a=%.2f, b=%.2f)\n", dist->params[0], dist->params[1]);
            break;
        case DIST_EXPONENTIAL:
            printf("Exponential(lambda=%.2f)\n", dist->params[0]);
            break;
        case DIST_NORMAL:
            printf("Normal(mu=%.2f, var=%.2f)\n", dist->params[0], dist->params[1]);
            break;
        case DIST_GAMMA:
            printf("Gamma(alpha=%.2f, beta=%.2f)\n", dist->params[0], dist->params[1]);
            break;
        case DIST_BETA:
            printf("Beta(alpha=%.2f, beta=%.2f)\n", dist->params[0], dist->params[1]);
            break;
        case DIST_CHISQ:
            printf("Chi-squared(df=%d)\n", dist->int_params[0]);
            break;
        case DIST_STUDENT:
            printf("Student(df=%d)\n", dist->int_params[0]);
            break;
        case DIST_FDIST:
            printf("F(df1=%d, df2=%d)\n", dist->int_params[0], dist->int_params[1]);
            break;
        case DIST_LOGNORMAL:
            printf("Log-normal(mu=%.2f, var=%.2f)\n", dist->params[0], dist->params[1]);
            break;
        case DIST_WEIBULL:
            printf("Weibull(lambda=%.2f, k=%.2f)\n", dist->params[0], dist->params[1]);
            break;
        default:
            printf("Unknown\n");
    }
}

Distribution *bernoulli(double p) {
    Distribution *d = (Distribution *)malloc(sizeof(Distribution));
    d->type = DIST_BERNOULLI;
    d->params[0] = p;
    return d;
}

Distribution *binomial(int n, double p) {
    Distribution *d = (Distribution *)malloc(sizeof(Distribution));
    d->type = DIST_BINOMIAL;
    d->int_params[0] = n;
    d->params[0] = p;
    return d;
}

Distribution *poisson(double lambda) {
    Distribution *d = (Distribution *)malloc(sizeof(Distribution));
    d->type = DIST_POISSON;
    d->params[0] = lambda;
    return d;
}

Distribution *geometric(double p) {
    Distribution *d = (Distribution *)malloc(sizeof(Distribution));
    d->type = DIST_GEOMETRIC;
    d->params[0] = p;
    return d;
}

Distribution *hypergeometric(int N, int K, int n) {
    Distribution *d = (Distribution *)malloc(sizeof(Distribution));
    d->type = DIST_HYPERGEOMETRIC;
    d->int_params[0] = N;
    d->int_params[1] = K;
    d->int_params[2] = n;
    return d;
}

Distribution *uniform(double a, double b) {
    Distribution *d = (Distribution *)malloc(sizeof(Distribution));
    d->type = DIST_UNIFORM;
    d->params[0] = a;
    d->params[1] = b;
    return d;
}

Distribution *exponential(double lambda) {
    Distribution *d = (Distribution *)malloc(sizeof(Distribution));
    d->type = DIST_EXPONENTIAL;
    d->params[0] = lambda;
    return d;
}

Distribution *normal(double mu, double var) {
    Distribution *d = (Distribution *)malloc(sizeof(Distribution));
    d->type = DIST_NORMAL;
    d->params[0] = mu;
    d->params[1] = var;
    return d;
}

Distribution *gamma(double alpha, double beta) {
    Distribution *d = (Distribution *)malloc(sizeof(Distribution));
    d->type = DIST_GAMMA;
    d->params[0] = alpha;
    d->params[1] = beta;
    return d;
}

Distribution *beta(double alpha, double beta) {
    Distribution *d = (Distribution *)malloc(sizeof(Distribution));
    d->type = DIST_BETA;
    d->params[0] = alpha;
    d->params[1] = beta;
    return d;
}

Distribution *chisq(int df) {
    Distribution *d = (Distribution *)malloc(sizeof(Distribution));
    d->type = DIST_CHISQ;
    d->int_params[0] = df;
    return d;
}

Distribution *student(int df) {
    Distribution *d = (Distribution *)malloc(sizeof(Distribution));
    d->type = DIST_STUDENT;
    d->int_params[0] = df;
    return d;
}

Distribution *fdist(int df1, int df2) {
    Distribution *d = (Distribution *)malloc(sizeof(Distribution));
    d->type = DIST_FDIST;
    d->int_params[0] = df1;
    d->int_params[1] = df2;
    return d;
}

Distribution *lognormal(double mu, double var) {
    Distribution *d = (Distribution *)malloc(sizeof(Distribution));
    d->type = DIST_LOGNORMAL;
    d->params[0] = mu;
    d->params[1] = var;
    return d;
}

Distribution *weibull(double lambda, double k) {
    Distribution *d = (Distribution *)malloc(sizeof(Distribution));
    d->type = DIST_WEIBULL;
    d->params[0] = lambda;
    d->params[1] = k;
    return d;
}

static double gamma_func(double x) {
    double g = 7.0;
    double c[] = {0.99999999999980993, 676.5203681218851, -1259.1392167224028,
                  771.32342877765313, -176.61502916214059, 12.507343278686905,
                  -0.13857109526572012, 9.9843695780195716e-6, 1.5056327351493116e-7};
    
    if (x < 0.5) return M_PI / (sin(M_PI * x) * gamma_func(1 - x));
    
    x -= 1.0;
    double a = c[0];
    for (int i = 1; i < 9; i++) a += c[i] / (x + i);
    double t = x + g + 0.5;
    return sqrt(2 * M_PI) * pow(t, x + 0.5) * exp(-t) * a;
}

static double beta_func(double a, double b) {
    return exp(gamma_func(a) + gamma_func(b) - gamma_func(a + b));
}

static double normal_pdf(double x, double mu, double var) {
    double sigma = sqrt(var);
    return exp(-0.5 * pow((x - mu) / sigma, 2)) / (sigma * sqrt(2 * M_PI));
}

double pdf(Distribution *dist, double x) {
    switch (dist->type) {
        case DIST_BERNOULLI:
            return (x == 0 || x == 1) ? dist->params[0] : 0.0;
        case DIST_BINOMIAL: {
            int k = (int)x;
            if (k < 0 || k > dist->int_params[0]) return 0;
            double p = dist->params[0];
            double binom = 1.0;
            for (int i = 0; i < k; i++) binom *= (dist->int_params[0] - i) / (double)(i + 1);
            return binom * pow(p, k) * pow(1 - p, dist->int_params[0] - k);
        }
        case DIST_POISSON: {
            double lambda = dist->params[0];
            double fact = 1.0;
            for (int i = 1; i <= (int)x; i++) fact *= i;
            return pow(lambda, x) * exp(-lambda) / fact;
        }
        case DIST_GEOMETRIC: {
            double p = dist->params[0];
            return p * pow(1 - p, (int)x - 1);
        }
        case DIST_HYPERGEOMETRIC: {
            int N = dist->int_params[0];
            int K = dist->int_params[1];
            int n = dist->int_params[2];
            int k = (int)x;
            if (k < 0 || k > n || k > K || k < n + K - N) return 0;
            double num = 1.0, den = 1.0;
            for (int i = 0; i < k; i++) num *= (K - i);
            for (int i = 0; i < k; i++) den *= (i + 1);
            double comb1 = num / den;
            num = 1.0; den = 1.0;
            for (int i = 0; i < n - k; i++) num *= (N - K - i);
            for (int i = 0; i < n - k; i++) den *= (i + 1);
            double comb2 = num / den;
            num = 1.0; den = 1.0;
            for (int i = 0; i < n; i++) num *= (N - i);
            for (int i = 0; i < n; i++) den *= (i + 1);
            double combN = num / den;
            return (comb1 * comb2) / combN;
        }
        case DIST_UNIFORM:
            return (x >= dist->params[0] && x <= dist->params[1]) ? 
                   1.0 / (dist->params[1] - dist->params[0]) : 0.0;
        case DIST_EXPONENTIAL:
            return dist->params[0] * exp(-dist->params[0] * x);
        case DIST_NORMAL:
            return normal_pdf(x, dist->params[0], dist->params[1]);
        case DIST_GAMMA: {
            double alpha = dist->params[0];
            double beta = dist->params[1];
            return pow(beta, alpha) * pow(x, alpha - 1) * exp(-beta * x) / gamma_func(alpha);
        }
        case DIST_BETA: {
            double alpha = dist->params[0];
            double beta = dist->params[1];
            return pow(x, alpha - 1) * pow(1 - x, beta - 1) / beta_func(alpha, beta);
        }
        case DIST_CHISQ: {
            int df = dist->int_params[0];
            return pow(x, df / 2.0 - 1) * exp(-x / 2.0) / pow(2, df / 2.0) / gamma_func(df / 2.0);
        }
        case DIST_STUDENT: {
            int df = dist->int_params[0];
            return pow(1 + x * x / df, -(df + 1) / 2.0) / sqrt(df * M_PI) / gamma_func((df + 1) / 2.0) * gamma_func(df / 2.0);
        }
        case DIST_FDIST: {
            int df1 = dist->int_params[0];
            int df2 = dist->int_params[1];
            double num = pow(df1 * x, df1) * pow(df2, df2);
            double den = pow(df1 * x + df2, df1 + df2);
            return sqrt(num / den) / x / beta_func(df1 / 2.0, df2 / 2.0);
        }
        case DIST_LOGNORMAL: {
            double mu = dist->params[0];
            double var = dist->params[1];
            double sigma = sqrt(var);
            return exp(-pow(log(x) - mu, 2) / (2 * sigma * sigma)) / (x * sigma * sqrt(2 * M_PI));
        }
        case DIST_WEIBULL: {
            double lambda = dist->params[0];
            double k = dist->params[1];
            return (k / lambda) * pow(x / lambda, k - 1) * exp(-pow(x / lambda, k));
        }
        default:
            return 0.0;
    }
}

static double normal_cdf(double x, double mu, double var) {
    double sigma = sqrt(var);
    return 0.5 * (1 + erf((x - mu) / (sigma * sqrt(2))));
}

static double erf_inv(double x) {
    double a1 = 0.254829592;
    double a2 = -0.284496736;
    double a3 = 1.421413741;
    double a4 = -1.453152027;
    double a5 = 1.061405429;
    double p = 0.3275911;
    
    double t = x;
    double y = x;
    for (int iter = 0; iter < 10; iter++) {
        double err = erf(y) - t;
        double dexp = exp(y * y);
        y = y - err * dexp * (1 + p * y * y);
    }
    return y;
}

double cdf(Distribution *dist, double x) {
    switch (dist->type) {
        case DIST_UNIFORM:
            if (x < dist->params[0]) return 0;
            if (x > dist->params[1]) return 1;
            return (x - dist->params[0]) / (dist->params[1] - dist->params[0]);
        case DIST_EXPONENTIAL:
            return 1 - exp(-dist->params[0] * x);
        case DIST_NORMAL:
            return normal_cdf(x, dist->params[0], dist->params[1]);
        case DIST_LOGNORMAL:
            return normal_cdf(log(x), dist->params[0], dist->params[1]);
        default:
            return 0.0;
    }
}

double quantile(Distribution *dist, double p) {
    if (p <= 0 || p >= 1) return 0;
    
    switch (dist->type) {
        case DIST_NORMAL: {
            double mu = dist->params[0];
            double var = dist->params[1];
            double sigma = sqrt(var);
            double t = p - 0.5;
            double x = mu + sigma * sqrt(2) * erf_inv(t);
            return x;
        }
        case DIST_BETA: {
            double alpha = dist->params[0];
            double beta = dist->params[1];
            double lower = 0.001, upper = 0.999, x;
            for (int i = 0; i < 50; i++) {
                x = (lower + upper) / 2;
                double bx = beta_cdf(x, alpha, beta);
                if (bx < p) lower = x;
                else upper = x;
            }
            return (lower + upper) / 2;
        }
        default:
            return 0.0;
    }
}

static double beta_cdf(double x, double a, double b) {
    if (x <= 0) return 0;
    if (x >= 1) return 1;
    return incomplete_beta(a, b, x);
}

static double incomplete_beta(double a, double b, double x) {
    double bt = (x == 0 || x == 1) ? 0 : 
        exp(log_gamma(a+b) - log_gamma(a) - log_gamma(b) + a*log(x) + b*log(1-x));
    if (x < (a+1)/(a+b+2)) {
        return bt * beta_cf(a, b, x) / a;
    }
    return 1 - bt * beta_cf(b, a, 1-x) / b;
}

static double beta_cf(double a, double b, double x) {
    double m, m2, aa, c = 1, d = 1 - (a+b)*x/(a+1), h = d;
    if (fabs(d) < 1e-30) d = 1e-30;
    d = 1/d;
    h = d;
    for (m = 1; m <= 100; m++) {
        m2 = 2*m;
        aa = m*(b-m)*x/((a+m2-1)*(a+m2));
        d = 1+aa*d;
        if (fabs(d) < 1e-30) d = 1e-30;
        c = 1+aa/c;
        if (fabs(c) < 1e-30) c = 1e-30;
        d = 1/d;
        h *= d*c;
        aa = -(a+m)*(a+b+m)*x/((a+m2)*(a+m2+1));
        d = 1+aa*d;
        if (fabs(d) < 1e-30) d = 1e-30;
        c = 1+aa/c;
        if (fabs(c) < 1e-30) c = 1e-30;
        d = 1/d;
        h *= d*c;
    }
    return h;
}

static double log_gamma(double x) {
    double c[] = {76.18009172947146, -86.50532032941677, 24.01409824083091,
                  -1.231739572450155, 0.1208650973866179e-2, -0.5395239384953e-5};
    double y = x;
    double tmp = x + 5.5;
    tmp -= (x + 0.5) * log(tmp);
    double ser = 1.000000000190015;
    for (int j = 0; j < 6; j++) ser += c[j] / ++y;
    return -tmp + log(2.5066282746310005 * ser / x);
}

double mean(Distribution *dist) {
    switch (dist->type) {
        case DIST_BERNOULLI: return dist->params[0];
        case DIST_BINOMIAL: return dist->int_params[0] * dist->params[0];
        case DIST_POISSON: return dist->params[0];
        case DIST_GEOMETRIC: return 1.0 / dist->params[0];
        case DIST_UNIFORM: return (dist->params[0] + dist->params[1]) / 2.0;
        case DIST_EXPONENTIAL: return 1.0 / dist->params[0];
        case DIST_NORMAL: return dist->params[0];
        case DIST_GAMMA: return dist->params[0] / dist->params[1];
        case DIST_BETA: return dist->params[0] / (dist->params[0] + dist->params[1]);
        case DIST_CHISQ: return dist->int_params[0];
        case DIST_STUDENT: return 0;
        case DIST_FDIST: {
            int df2 = dist->int_params[1];
            return (df2 > 2) ? (double)dist->int_params[0] / (df2 - 2) : 0;
        }
        case DIST_LOGNORMAL: return exp(dist->params[0] + dist->params[1] / 2.0);
        case DIST_WEIBULL: return dist->params[1] * gamma_func(1 + 1.0 / dist->params[0]);
        default: return 0.0;
    }
}

double variance(Distribution *dist) {
    switch (dist->type) {
        case DIST_BERNOULLI: return dist->params[0] * (1 - dist->params[0]);
        case DIST_BINOMIAL: return dist->int_params[0] * dist->params[0] * (1 - dist->params[0]);
        case DIST_POISSON: return dist->params[0];
        case DIST_GEOMETRIC: return (1 - dist->params[0]) / (dist->params[0] * dist->params[0]);
        case DIST_UNIFORM: return pow(dist->params[1] - dist->params[0], 2) / 12.0;
        case DIST_EXPONENTIAL: return 1.0 / (dist->params[0] * dist->params[0]);
        case DIST_NORMAL: return dist->params[1];
        case DIST_GAMMA: return dist->params[0] / (dist->params[1] * dist->params[1]);
        case DIST_BETA: {
            double a = dist->params[0], b = dist->params[1];
            return (a * b) / (pow(a + b, 2) * (a + b + 1));
        }
        case DIST_CHISQ: return 2 * dist->int_params[0];
        case DIST_STUDENT: return (dist->int_params[0] > 2) ? (double)dist->int_params[0] / (dist->int_params[0] - 2) : 0;
        case DIST_FDIST: {
            int df1 = dist->int_params[0], df2 = dist->int_params[1];
            return (df2 > 4) ? 2.0 * df1 * df1 * (df2 + df1 - 2) / (df2 * pow(df2 - 2, 2) * (df2 - 4)) : 0;
        }
        case DIST_LOGNORMAL: {
            double mu = dist->params[0];
            double var = dist->params[1];
            return (exp(var) - 1) * exp(2 * mu + var);
        }
        case DIST_WEIBULL: {
            double lambda = dist->params[0];
            double k = dist->params[1];
            return pow(lambda, 2) * (gamma_func(1 + 2 / k) - pow(gamma_func(1 + 1 / k), 2));
        }
        default: return 0.0;
    }
}

double stddev(Distribution *dist) {
    return sqrt(variance(dist));
}