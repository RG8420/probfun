#ifndef BAYES_H
#define BAYES_H

#include <stdlib.h>
#include <string.h>
#include "dist.h"

typedef struct {
    Distribution *prior;
    int observed_successes;
    int observed_trials;
    Distribution *posterior;
} BayesianModel;

typedef struct {
    BayesianModel models[MAX_MATRICES];
    int count;
} BayesStore;

void bayes_init(BayesStore *store);
int bayes_set_prior(BayesStore *store, const char *name, Distribution *prior);
BayesianModel *bayes_get(BayesStore *store, const char *name);
void bayes_update(BayesianModel *model, int successes, int trials);
Distribution *bayes_posterior(BayesianModel *model);
double bayes_map(BayesianModel *model);
void bayes_credible_interval(BayesianModel *model, double level, double *lower, double *upper);
void bayes_list(BayesStore *store);

#endif