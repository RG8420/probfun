# Probability Calculator (probcalc)

A C-based CLI calculator for probability theory computations, installable via one-liner curl.

## Installation

```bash
curl -sSL https://raw.githubusercontent.com/RG8420/probfun/main/install.sh | bash
```

## Input Syntax

```bash
X = normal(0, 1)      # Normal distribution with μ=0, σ²=1
p = binomial(10, 0.5) # Binomial with n=10, p=0.5
```

## Features

### Combinatorics
- `factorial(n)` - n!
- `permutations(n, r)` - P(n,r)
- `combinations(n, r)` - C(n,r)

### Discrete Distributions
- Bernoulli, Binomial, Poisson, Geometric, Hypergeometric

### Continuous Distributions
- Uniform, Exponential, Normal, Gamma, Beta
- Chi-squared, Student's t, F-distribution
- Log-normal, Weibull

### Distribution Functions
- `pdf(dist, x)` - Probability density/mass function
- `cdf(dist, x)` - Cumulative distribution function
- `quantile(dist, p)` - Inverse CDF
- `mean(dist)`, `variance(dist)`, `stddev(dist)`

### Random Sampling
- `rand(dist, n)` - Generate n random samples

### Bayesian Statistics
- Define priors: `prior A = beta(2, 5)`
- Update posterior: `observed A 50 100`
- MAP estimate: `map(A)`
- Credible intervals: `ci(A, 0.95)`

### Statistical Tests
- Chi-squared test, t-test, z-test

## Example Session

```
=== Probability Calculator (probcalc) ===
Type 'help' for operations, 'quit' to exit.

probcalc> X = normal(0, 1)
Stored as X:
Normal(mu=0.00, var=1.00)
probcalc> pdf(X, 0)
  pdf = 0.398942
probcalc> cdf(X, 1.96)
  cdf = 0.975002
probcalc> mean(X)
  mean = 0.000000
probcalc> variance(X)
  variance = 1.000000
probcalc> rand(X, 5)
  Samples:
  0.1234
  -1.2345
  0.5678
  -0.9012
  2.3456
probcalc> prior A = beta(2, 5)
Stored prior as A:
Beta(alpha=2.00, beta=5.00)
probcalc> observed A 50 100
Updated posterior for A:
Beta(alpha=52.00, beta=55.00)
probcalc> map(A)
  MAP = 0.510000
probcalc> ci(A, 0.95)
  HPD 95% CI: [0.3400, 0.7200]
probcalc> factorial(5)
  120
probcalc> combinations(5, 2)
  10
probcalc> quit
Goodbye!
```

## Uninstall

```bash
bash install.sh --uninstall
```

## License

MIT License