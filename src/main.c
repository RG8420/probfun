#include "dist.h"
#include "random.h"
#include "combi.h"
#include "bayes.h"
#include "stats.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define MAX_NAME_LEN 64

static DistStore store;
static BayesStore bayes_store;

#define HELP_TEXT \
"===========================================================\n" \
"          PROBABILITY CALCULATOR HELP\n" \
"===========================================================\n" \
"  Input Syntax:\n" \
"    X = normal(0, 1)        # Create distribution\n" \
"    prior A = beta(2, 5)    # Bayesian prior\n" \
"\n" \
"  Combinatorics:\n" \
"    factorial(n)             # n!\n" \
"    permutations(n, r)      # P(n,r)\n" \
"    combinations(n, r)      # C(n,r)\n" \
"\n" \
"  Discrete Distributions:\n" \
"    bernoulli(p), binomial(n,p), poisson(lambda)\n" \
"    geometric(p), hypergeometric(N,K,n)\n" \
"\n" \
"  Continuous Distributions:\n" \
"    uniform(a,b), exponential(lambda)\n" \
"    normal(mu, var), gamma(a,b), beta(a,b)\n" \
"    chisq(df), student(df), f(df1,df2)\n" \
"    lognormal(mu, var), weibull(lambda,k)\n" \
"\n" \
"  Distribution Functions:\n" \
"    pdf(X, x), cdf(X, x), quantile(X, p)\n" \
"    mean(X), variance(X), stddev(X)\n" \
"    rand(X, n)               # n random samples\n" \
"\n" \
"  Bayesian:\n" \
"    prior name = beta(a,b)  # Set prior\n" \
"    bayes name observed(s,t)# Update posterior\n" \
"    map(name)               # MAP estimate\n" \
"    ci(name, level)         # HPD credible interval\n" \
"\n" \
"  Statistical Tests:\n" \
"    chisq_test([obs], [exp]) # Chi-squared test\n" \
"    ttest(data, mu)         # One-sample t-test\n" \
"    ttest2(data1, data2)    # Two-sample t-test\n" \
"    ztest(data, mu, sigma)  # Z-test\n" \
"\n" \
"  Utilities:\n" \
"    list                    # List stored distributions\n" \
"    blist                   # List Bayesian models\n" \
"    clear                   # Clear all\n" \
"    seed(n)                 # Set random seed\n" \
"    help, quit\n" \
"===========================================================\n"

static void print_prompt(void) {
    printf("probcalc> ");
    fflush(stdout);
}

static int parse_assignment(const char *input, char *var_name, char *expr) {
    const char *eq = strchr(input, '=');
    if (!eq) return 0;
    
    if (strncmp(input, "prior ", 6) == 0) {
        const char *name_start = input + 6;
        while (*name_start == ' ') name_start++;
        const char *name_end = eq;
        while (name_end > name_start && *(name_end-1) == ' ') name_end--;
        int len = name_end - name_start;
        if (len > 0 && len < MAX_NAME_LEN) {
            strncpy(var_name, name_start, len);
            var_name[len] = '\0';
        } else {
            return 0;
        }
        strcpy(expr, eq + 1);
        return strlen(var_name) > 0;
    }
    
    int j = 0;
    for (int i = 0; i < eq - input && j < MAX_NAME_LEN - 1; i++) {
        if (input[i] != ' ' && input[i] != '\t') {
            var_name[j++] = input[i];
        }
    }
    var_name[j] = '\0';
    
    strcpy(expr, eq + 1);
    return strlen(var_name) > 0;
}

static void trim(char *s) {
    while (*s == ' ') {
        memmove(s, s+1, strlen(s));
    }
    int len = strlen(s);
    while (len > 0 && (s[len-1] == ' ' || s[len-1] == ')')) {
        s[--len] = '\0';
    }
}

static Distribution *parse_distribution(char *s) {
    trim(s);
    
    if (strncmp(s, "bernoulli(", 10) == 0) {
        double p = atof(s + 10);
        return bernoulli(p);
    }
    if (strncmp(s, "binomial(", 9) == 0) {
        int n = atoi(s + 9);
        char *p = strchr(s, ',');
        double pval = p ? atof(p + 1) : 0.5;
        return binomial(n, pval);
    }
    if (strncmp(s, "poisson(", 8) == 0) {
        double lambda = atof(s + 8);
        return poisson(lambda);
    }
    if (strncmp(s, "geometric(", 9) == 0) {
        double p = atof(s + 9);
        return geometric(p);
    }
    if (strncmp(s, "hypergeometric(", 15) == 0) {
        int N = atoi(s + 15);
        char *p1 = strchr(s + 15, ',');
        int K = p1 ? atoi(p1 + 1) : 0;
        char *p2 = p1 ? strchr(p1 + 1, ',') : NULL;
        int n = p2 ? atoi(p2 + 1) : 0;
        return hypergeometric(N, K, n);
    }
    if (strncmp(s, "uniform(", 8) == 0) {
        double a = atof(s + 8);
        char *p = strchr(s + 8, ',');
        double b = p ? atof(p + 1) : 1.0;
        return uniform(a, b);
    }
    if (strncmp(s, "exponential(", 12) == 0) {
        double lambda = atof(s + 12);
        return exponential(lambda);
    }
    if (strncmp(s, "normal(", 7) == 0) {
        double mu = atof(s + 7);
        char *p = strchr(s + 7, ',');
        double var = p ? atof(p + 1) : 1.0;
        return normal(mu, var);
    }
    if (strncmp(s, "gamma(", 6) == 0) {
        double a = atof(s + 6);
        char *p = strchr(s + 6, ',');
        double b = p ? atof(p + 1) : 1.0;
        return gamma(a, b);
    }
    if (strncmp(s, "beta(", 5) == 0) {
        double a = atof(s + 5);
        char *p = strchr(s + 5, ',');
        double b = p ? atof(p + 1) : 1.0;
        return beta(a, b);
    }
    if (strncmp(s, "chisq(", 6) == 0) {
        int df = atoi(s + 6);
        return chisq(df);
    }
    if (strncmp(s, "student(", 8) == 0) {
        int df = atoi(s + 8);
        return student(df);
    }
    if (strncmp(s, "fdist(", 6) == 0) {
        int df1 = atoi(s + 6);
        char *p = strchr(s + 6, ',');
        int df2 = p ? atoi(p + 1) : 1;
        return fdist(df1, df2);
    }
    if (strncmp(s, "lognormal(", 10) == 0) {
        double mu = atof(s + 10);
        char *p = strchr(s + 10, ',');
        double var = p ? atof(p + 1) : 1.0;
        return lognormal(mu, var);
    }
    if (strncmp(s, "weibull(", 8) == 0) {
        double lambda = atof(s + 8);
        char *p = strchr(s + 8, ',');
        double k = p ? atof(p + 1) : 1.0;
        return weibull(lambda, k);
    }
    
    return NULL;
}

static void handle_command(char *input) {
    char var_name[MAX_NAME_LEN];
    char expr[1024];
    
    if (parse_assignment(input, var_name, expr)) {
        trim(expr);
        
        if (strncmp(input, "prior ", 6) == 0 || strncmp(var_name, "prior ", 6) == 0) {
            Distribution *prior = parse_distribution(expr);
            if (prior) {
                strncpy(prior->name, var_name, MAX_NAME_LEN - 1);
                prior->name[MAX_NAME_LEN - 1] = '\0';
                bayes_set_prior(&bayes_store, var_name, prior);
                printf("Stored prior as %s:\n", var_name);
                dist_print(prior);
            } else {
                printf("Error: Invalid prior distribution\n");
            }
        } else if (strncmp(expr, "bayes ", 6) == 0 || strncmp(expr, "observed(", 8) == 0 || strncmp(expr, "update ", 7) == 0) {
            char *name;
            if (strncmp(expr, "observed(", 8) == 0) {
                name = expr + 8;
            } else if (strncmp(expr, "update ", 7) == 0) {
                name = expr + 7;
            } else {
                name = expr + 6;
            }
            trim(name);
            
            if (strncmp(expr, "observed(", 8) == 0) {
                char *comma1 = strchr(name, ',');
                if (comma1) {
                    *comma1 = '\0';
                    char *success_str = comma1 + 1;
                    char *comma2 = strchr(success_str, ',');
                    if (comma2) {
                        *comma2 = '\0';
                        int successes = atoi(success_str);
                        int trials = atoi(comma2 + 1);
                        BayesianModel *model = bayes_get(&bayes_store, name);
                        if (model) {
                            bayes_update(model, successes, trials);
                            printf("Updated posterior for %s:\n", name);
                            dist_print(model->posterior);
                        } else {
                            printf("Error: No prior found for %s\n", name);
                        }
                    }
                }
            } else {
                char *obs = strchr(name, ' ');
                if (obs) {
                    *obs = '\0';
                    obs++;
                    trim(obs);
                    char *comma = strchr(obs, ',');
                    int successes = atoi(obs);
                    int trials = comma ? atoi(comma + 1) : 0;
                    BayesianModel *model = bayes_get(&bayes_store, name);
                    if (model) {
                        bayes_update(model, successes, trials);
                        printf("Updated posterior for %s:\n", name);
                        dist_print(model->posterior);
                    } else {
                        printf("Error: No prior found for %s\n", name);
                    }
                }
            }
        } else {
            Distribution *dist = parse_distribution(expr);
            if (dist) {
                dist_set(&store, var_name, dist);
                printf("Stored as %s:\n", var_name);
                dist_print(dist);
            } else {
                printf("Error: Invalid distribution syntax\n");
            }
        }
    } else {
        trim(input);
        
        if (strcmp(input, "list") == 0) {
            dist_list(&store);
        } else if (strcmp(input, "blist") == 0) {
            bayes_list(&bayes_store);
        } else if (strcmp(input, "clear") == 0) {
            dist_clear(&store);
            bayes_init(&bayes_store);
            printf("  Cleared all.\n");
        } else if (strcmp(input, "help") == 0) {
            printf("%s", HELP_TEXT);
        } else if (strncmp(input, "seed(", 5) == 0) {
            unsigned int seed = atoi(input + 5);
            seed_random(seed);
            printf("  Random seed set to %u\n", seed);
        } else if (strncmp(input, "observed(", 9) == 0) {
            char *args = input + 9;
            int len = strlen(args);
            if (len > 0 && args[len-1] == ')') args[--len] = '\0';
            char *comma1 = strchr(args, ',');
            if (comma1) {
                *comma1 = '\0';
                char *success_str = comma1 + 1;
                char *comma2 = strchr(success_str, ',');
                if (comma2) {
                    *comma2 = '\0';
                    int successes = atoi(success_str);
                    int trials = atoi(comma2 + 1);
                    BayesianModel *model = bayes_get(&bayes_store, args);
                    if (model) {
                        bayes_update(model, successes, trials);
                        printf("Updated posterior for %s:\n", args);
                        dist_print(model->posterior);
                    } else {
                        printf("Error: No prior found for %s\n", args);
                    }
                }
            }
        } else if (strncmp(input, "pdf(", 4) == 0) {
            char *args = input + 4;
            int len = strlen(args);
            if (len > 0 && args[len-1] == ')') args[--len] = '\0';
            char *comma = strchr(args, ',');
            if (comma) {
                *comma = '\0';
                Distribution *dist = dist_get(&store, args);
                double x = atof(comma + 1);
                if (dist) printf("  pdf = %.6f\n", pdf(dist, x));
                else printf("Error: Unknown distribution\n");
            }
        } else if (strncmp(input, "cdf(", 4) == 0) {
            char *args = input + 4;
            int len = strlen(args);
            if (len > 0 && args[len-1] == ')') args[--len] = '\0';
            char *comma = strchr(args, ',');
            if (comma) {
                *comma = '\0';
                Distribution *dist = dist_get(&store, args);
                double x = atof(comma + 1);
                if (dist) printf("  cdf = %.6f\n", cdf(dist, x));
                else printf("Error: Unknown distribution\n");
            }
        } else if (strncmp(input, "quantile(", 9) == 0) {
            char *args = input + 9;
            int len = strlen(args);
            if (len > 0 && args[len-1] == ')') args[--len] = '\0';
            char *comma = strchr(args, ',');
            if (comma) {
                *comma = '\0';
                Distribution *dist = dist_get(&store, args);
                double p = atof(comma + 1);
                if (dist) printf("  quantile = %.6f\n", quantile(dist, p));
                else printf("Error: Unknown distribution\n");
            }
        } else if (strncmp(input, "rand(", 5) == 0) {
            char *args = input + 5;
            int len = strlen(args);
            if (len > 0 && args[len-1] == ')') args[--len] = '\0';
            char *comma = strchr(args, ',');
            if (comma) {
                *comma = '\0';
                Distribution *dist = dist_get(&store, args);
                int n = atoi(comma + 1);
                if (dist && n > 0) {
                    double *samples = rand_samples(dist, n);
                    printf("  Samples:\n");
                    for (int i = 0; i < n; i++) {
                        printf("  %.4f\n", samples[i]);
                    }
                    free(samples);
                } else {
                    printf("Error: Unknown distribution or invalid n\n");
                }
            }
        } else if (strncmp(input, "mean(", 5) == 0) {
            char *name = input + 5;
            int len = strlen(name);
            if (len > 0 && name[len-1] == ')') name[--len] = '\0';
            Distribution *dist = dist_get(&store, name);
            if (dist) printf("  mean = %.6f\n", mean(dist));
            else printf("Error: Unknown distribution\n");
        } else if (strncmp(input, "variance(", 9) == 0) {
            char *name = input + 9;
            int len = strlen(name);
            if (len > 0 && name[len-1] == ')') name[--len] = '\0';
            Distribution *dist = dist_get(&store, name);
            if (dist) printf("  variance = %.6f\n", variance(dist));
            else printf("Error: Unknown distribution\n");
        } else if (strncmp(input, "stddev(", 7) == 0) {
            char *name = input + 7;
            int len = strlen(name);
            if (len > 0 && name[len-1] == ')') name[--len] = '\0';
            Distribution *dist = dist_get(&store, name);
            if (dist) printf("  stddev = %.6f\n", stddev(dist));
            else printf("Error: Unknown distribution\n");
        } else if (strncmp(input, "rand(", 5) == 0) {
            char *args = input + 5;
            char *comma = strchr(args, ',');
            if (comma) {
                *comma = '\0';
                Distribution *dist = dist_get(&store, args);
                int n = atoi(comma + 1);
                if (dist && n > 0) {
                    double *samples = rand_samples(dist, n);
                    printf("  Samples:\n");
                    for (int i = 0; i < n; i++) {
                        printf("  %.4f\n", samples[i]);
                    }
                    free(samples);
                } else {
                    printf("Error: Unknown distribution or invalid n\n");
                }
            }
        } else if (strncmp(input, "map(", 4) == 0) {
            char *name = input + 4;
            BayesianModel *model = bayes_get(&bayes_store, name);
            if (model) {
                printf("  MAP = %.6f\n", bayes_map(model));
            } else {
                printf("Error: Unknown Bayesian model\n");
            }
        } else if (strncmp(input, "ci(", 3) == 0) {
            char *args = input + 3;
            char *comma = strchr(args, ',');
            if (comma) {
                *comma = '\0';
                BayesianModel *model = bayes_get(&bayes_store, args);
                double level = atof(comma + 1);
                if (model) {
                    double lower, upper;
                    bayes_credible_interval(model, level, &lower, &upper);
                    printf("  HPD %.0f%% CI: [%.4f, %.4f]\n", level * 100, lower, upper);
                } else {
                    printf("Error: Unknown Bayesian model\n");
                }
            }
        } else if (strncmp(input, "factorial(", 9) == 0) {
            char *arg = input + 9;
            int len = strlen(arg);
            if (len > 0 && arg[0] == '(') memmove(arg, arg+1, len--);
            if (len > 0 && arg[len-1] == ')') arg[--len] = '\0';
            int n = atoi(arg);
            printf("  %lld\n", factorial(n));
        } else if (strncmp(input, "permutations(", 12) == 0) {
            char *args = input + 12;
            int len = strlen(args);
            if (len > 0 && args[0] == '(') memmove(args, args+1, len--);
            if (len > 0 && args[len-1] == ')') args[--len] = '\0';
            char *comma = strchr(args, ',');
            if (comma) {
                *comma = '\0';
                int n = atoi(args);
                int r = atoi(comma + 1);
                printf("  %lld\n", permutations(n, r));
            }
        } else if (strncmp(input, "combinations(", 12) == 0) {
            char *args = input + 12;
            int len = strlen(args);
            if (len > 0 && args[0] == '(') memmove(args, args+1, len--);
            if (len > 0 && args[len-1] == ')') args[--len] = '\0';
            char *comma = strchr(args, ',');
            if (comma) {
                *comma = '\0';
                int n = atoi(args);
                int r = atoi(comma + 1);
                printf("  %lld\n", combinations(n, r));
            }
        } else {
            Distribution *dist = dist_get(&store, input);
            if (dist) {
                dist_print(dist);
            } else {
                printf("Error: Unknown command. Type 'help' for options.\n");
            }
        }
    }
}

int main() {
    dist_init(&store);
    bayes_init(&bayes_store);
    seed_random(time(NULL));
    
    printf("=== Probability Calculator (probcalc) ===\n");
    printf("Type 'help' for operations, 'quit' to exit.\n");
    
    char line[1024];
    while (1) {
        print_prompt();
        if (!fgets(line, sizeof(line), stdin)) break;
        line[strcspn(line, "\n")] = 0;
        
        if (strcmp(line, "quit") == 0 || strcmp(line, "exit") == 0) break;
        if (strlen(line) == 0) continue;
        
        handle_command(line);
    }
    
    printf("Goodbye!\n");
    return 0;
}