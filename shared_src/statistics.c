#include "statistics.h"

double average(double (*f)(), int n) {
    double sum = 0;
    for (int i = 0; i < n; i++) {
        sum += f();
    }
    return sum / n;
}

double mean(double *a, int n) {
    double sum = 0;
    for (int i = 0; i < n; i++) {
        sum += a[i];
    }
    return sum / n;
}

double meansq(double *a, int n) {
    double sum = 0;
    for (int i = 0; i < n; i++) {
        sum += a[i] * a[i];
    }
    return sum / n;
}

double stddev(double *a, int n) {
    return sqrt( (meansq(a, n) - pow(mean(a, n), 2)) * n / (n - 1) );
}
