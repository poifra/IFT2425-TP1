//------------------------------------------------------
// module  : Tp-IFT2425-I.1.c
// author  : François Poitras & Guillaume Noël-Martel
// date    :
// version : 1.0
// language: C
// note    :
//------------------------------------------------------

//------------------------------------------------
// FICHIERS INCLUS -------------------------------
//------------------------------------------------
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <iostream>
#include <new>
#include <time.h>

//------------------------------------------------
// DEFINITIONS -----------------------------------
//------------------------------------------------
#define CARRE(X) ((X)*(X))
#define CUBE(X)  ((X)*(X)*(X))

const double epsilon = 1E-5;

/*We have the function (where N is the size of y):
f(c; y[i]) =
    sum(i=1, N, y_i^c * ln(y_i)) / sum(i=1, y_i^c)
    - 1/c
    - 1/N*sum(i=1, N, ln(y_i))
We can thus calculate it using three sums over y.*/
double f(double* y, int N, double x)
{
    double sumA = 0, sumB = 0, sumC = 0;
    for (int i = 0; i < N; i++)
    {
        double ln_yi = log(y[i]);
        double pow_yi_x = pow(y[i], x);

        sumA += pow_yi_x * ln_yi;
        sumB += pow_yi_x;
        sumC += ln_yi;
    }
    return (sumA / sumB) - (1 / x) - (1 / (double)N) * sumC;
}
/* When differentiating over c, the third sum disappears; we also have that
there is quotient of functions and so we have (f/g)' = (f'g - fg')/g².
We have that the term (y_i^c)' = ln(y_i)*y_i^c which is coincidently the
same term as the first sum, we thus have that g' = f and the term of f' is
ln(y_i)^2*y_i^c.
We thus have:
sumA = sum(i=1, N, y_i^c * ln(y_i))
sumdA = sum(i=1, N, y_i^c * ln^2(y_i))
sumB = sum(i=1, N, ln(y_i))
sumdB = sumA
f'(c; y[i]) = (sumdA * sumB - sumA * sumdB) / sumdB^2 + 1/c^2 */
double df(double* y, int N, double x)
{
    double sumA = 0, sumdA = 0, sumB = 0, sumdB = 0;
    for (int i = 0; i < N; i++)
    {
        double ln_yi = log(y[i]);
        double pow_yi_x = pow(y[i], x);

        sumA += pow_yi_x * ln_yi;
        sumdA += pow_yi_x * CARRE(ln_yi);
        sumB += pow_yi_x;
    }
    sumdB = sumA;
    return (sumdA * sumB - sumA * sumdB) / CARRE(sumB) + (1 / CARRE(x));
}
double df_epsilon_a(double* y, int N, double x)
{
    return (f(y, N, x + epsilon)
            - f(y, N, x)
        ) / epsilon;
}
double df_epsilon_b(double* y, int N, double x)
{
    return (- f(y, N, x + 2.0 * epsilon)
            + 8.0 * f(y, N, x + epsilon)
            - 8.0 * f(y, N, x - epsilon)
            + f(y, N, x - 2.0 * epsilon)
        ) / (12.0 * epsilon);
}

double approx_newton(double (*df)(double*, int, double),
                          double* y, int N, double x0, double target,
                          double tolerance, int maxIter)
{
    double x = x0;
    for (int i = 0; i < maxIter; i++)
    {
        double val = f(y, N, x);
        if (fabs(target - val) < tolerance) {
            printf("Tolérance atteinte après %i itérations\n", (i+1));
            break;
        }
        x = x - val/(*df)(y, N, x);
    }
    return x;
}

//----------------------------------------------------------
//----------------------------------------------------------
// PROGRAMME PRINCIPAL -------------------------------------
//----------------------------------------------------------
//----------------------------------------------------------
int main(int argc, char** argv)
{
    double y[] = { 0.11, 0.24, 0.27, 0.52, 1.13, 1.54, 1.71, 1.84, 1.92, 2.01 };
    int N = sizeof(y) / sizeof(double);

    double x0 = 0.25, target = 0, tolerance = 10E-6;
    int maxIter = 100;

    printf("Évaluation de c_mv à l'aide de l'algorithme de Newton:\n");

    printf("1.a)\nUtilisant un epsilon simple ...\n");
    double c_mv_a = approx_newton(*df_epsilon_a, y, N, x0, target, tolerance, maxIter);
    printf("c_mv = %f\n", c_mv_a);

    printf("1.b)\nUtilisant un epsilon plus compliqué ...\n");
    double c_mv_b = approx_newton(*df_epsilon_b, y, N, x0, target, tolerance, maxIter);
    printf("c_mv = %f\n", c_mv_b);

    printf("1.c)\nUtilisant la dérivée analytique ...\n");
    double c_mv_c = approx_newton(*df, y, N, x0, target, tolerance, maxIter);
    printf("c_mv = %f\n", c_mv_c);

    return 0;
}