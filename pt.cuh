#ifndef PT_DEVICE_CUH
#define PT_DEVICE_CUH

#include <cuda_runtime.h>
#include <math_constants.h>

// Device function to compute the continued fraction for the incomplete beta function
__device__ double betacf(double a, double b, double x) {
    const int MAX_ITER = 100;
    const double EPS = 3.0e-7;
    const double FPMIN = 1.0e-30;

    double qab = a + b;
    double qap = a + 1.0;
    double qam = a - 1.0;

    double c = 1.0;
    double d = 1.0 - qab * x / qap;
    if (fabs(d) < FPMIN) d = FPMIN;
    d = 1.0 / d;
    double h = d;

    for (int m = 1; m <= MAX_ITER; m++) {
        int m2 = 2 * m;
        // Even step
        double aa = m * (b - m) * x / ((qam + m2) * (a + m2));
        d = 1.0 + aa * d;
        if (fabs(d) < FPMIN) d = FPMIN;
        c = 1.0 + aa / c;
        if (fabs(c) < FPMIN) c = FPMIN;
        d = 1.0 / d;
        h *= d * c;
        // Odd step
        aa = -(a + m) * (qab + m) * x / ((a + m2) * (qap + m2));
        d = 1.0 + aa * d;
        if (fabs(d) < FPMIN) d = FPMIN;
        c = 1.0 + aa / c;
        if (fabs(c) < FPMIN) c = FPMIN;
        d = 1.0 / d;
        double del = d * c;
        h *= del;
        if (fabs(del - 1.0) < EPS) break;
    }

    return h;
}

// Device function to compute the incomplete beta function
__device__ double betai(double a, double b, double x) {
    if (x < 0.0 || x > 1.0) return 0.0; // Invalid input
    if (x == 0.0 || x == 1.0) return x; // Edge cases

    // Compute ln(Beta(a, b))
    double ln_beta = lgamma(a) + lgamma(b) - lgamma(a + b);

    // Compute front factor
    double front = exp(log(x) * a + log(1.0 - x) * b - ln_beta) / a;

    // Compute continued fraction
    double cf = betacf(a, b, x);

    return front * cf;
}

// Device function to compute the CDF of the Student's t-distribution
__device__ double pt(double t, double df) {
    double x = df / (df + t * t);
    double a = df / 2.0;
    double b = 0.5;

    double ibeta = betai(a, b, x);

    // Compute the CDF
    double cdf = 0.5 * (1.0 + (t > 0 ? 1.0 : -1.0) * (1.0 - ibeta));

    return cdf;
}

#endif // PT_DEVICE_CUH
