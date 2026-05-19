#include "combi.h"

long long factorial(int n) {
    if (n <= 1) return 1;
    long long result = 1;
    for (int i = 2; i <= n; i++) result *= i;
    return result;
}

long long permutations(int n, int r) {
    if (r > n) return 0;
    long long result = 1;
    for (int i = 0; i < r; i++) result *= (n - i);
    return result;
}

long long combinations(int n, int r) {
    if (r > n) return 0;
    if (r == 0 || r == n) return 1;
    if (r > n - r) r = n - r;
    long long result = 1;
    for (int i = 0; i < r; i++) {
        result *= (n - i);
        result /= (i + 1);
    }
    return result;
}