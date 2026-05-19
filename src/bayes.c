#include "bayes.h"
#include "dist.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

void bayes_init(BayesStore *store) {
    store->count = 0;
}

int bayes_set_prior(BayesStore *store, const char *name, Distribution *prior) {
    for (int i = 0; i < store->count; i++) {
        if (strcmp(store->models[i].prior->name, name) == 0) {
            store->models[i].prior = prior;
            store->models[i].observed_successes = 0;
            store->models[i].observed_trials = 0;
            store->models[i].posterior = NULL;
            return 1;
        }
    }
    if (store->count >= MAX_MATRICES) return 0;
    BayesianModel *model = &store->models[store->count];
    model->prior = prior;
    model->observed_successes = 0;
    model->observed_trials = 0;
    model->posterior = NULL;
    store->count++;
    return 1;
}

BayesianModel *bayes_get(BayesStore *store, const char *name) {
    for (int i = 0; i < store->count; i++) {
        if (store->models[i].prior && 
            strcmp(store->models[i].prior->name, name) == 0) {
            return &store->models[i];
        }
    }
    return NULL;
}

void bayes_update(BayesianModel *model, int successes, int trials) {
    model->observed_successes += successes;
    model->observed_trials += trials;
    
    if (model->posterior) free(model->posterior);
    
    model->posterior = (Distribution *)malloc(sizeof(Distribution));
    model->posterior->type = DIST_BETA;
    model->posterior->params[0] = model->prior->params[0] + model->observed_successes;
    model->posterior->params[1] = model->prior->params[1] + model->observed_trials - model->observed_successes;
}

Distribution *bayes_posterior(BayesianModel *model) {
    if (!model->posterior) {
        bayes_update(model, 0, 0);
    }
    return model->posterior;
}

double bayes_map(BayesianModel *model) {
    Distribution *post = bayes_posterior(model);
    if (!post || post->type != DIST_BETA) return 0;
    double alpha = post->params[0];
    double beta = post->params[1];
    return (alpha - 1) / (alpha + beta - 2);
}

void bayes_credible_interval(BayesianModel *model, double level, double *lower, double *upper) {
    Distribution *post = bayes_posterior(model);
    if (!post || post->type != DIST_BETA) {
        *lower = 0;
        *upper = 1;
        return;
    }
    
    double alpha = post->params[0];
    double beta = post->params[1];
    double half_alpha = (1 - level) / 2.0;
    
    *lower = quantile(post, half_alpha);
    *upper = quantile(post, 1 - half_alpha);
}

void bayes_list(BayesStore *store) {
    if (store->count == 0) {
        printf("  No Bayesian models stored.\n");
        return;
    }
    for (int i = 0; i < store->count; i++) {
        BayesianModel *m = &store->models[i];
        if (m->prior) {
            printf("  %s: Prior=Beta(%.2f, %.2f)", m->prior->name, m->prior->params[0], m->prior->params[1]);
            if (m->posterior) {
                printf(", Posterior=Beta(%.2f, %.2f)", m->posterior->params[0], m->posterior->params[1]);
            }
            printf("\n");
        }
    }
}