#include "header.h"
#include <cmath>

double u_exact(double r, double z, double t)
{
   return z * sin(t);
}

double lambda_f(int region, double r, double z, double t)
{
   return 1.0;
}

double gamma_f(int region, double r, double z, double t)
{
   return 0.0;
}

double f_f(int region, double r, double z, double t)
{
   return z * cos(t) - z * sin(t);
}

double u_s1(int funcNum, double r, double z, double t)
{
   return u_exact(r, z, t);
}

double theta_s2(int funcNum, double r, double z, double t)
{
   return 0.0;
}

double beta_s3(int funcNum, double r, double z, double t)
{
   return 0.0;
}

double u_beta_s3(int funcNum, double r, double z, double t)
{
   return 0.0;
}

double u0_f(double r, double z)
{
   return 0.0;
}

double v0_f(double r, double z)
{
   return z;
}

double a0_f(double r, double z)
{
   return 0.0;
}


