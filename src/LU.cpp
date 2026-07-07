#include "header.h"
#include <cmath>
#include <iostream>

double scalar(const std::vector<double> &a, const std::vector<double> &b)
{
   double s = 0.0;

   for (size_t i = 0; i < a.size(); i++) {
      s += a[i] * b[i];
   }

   return s;
}

void multA(const matrix &A, const std::vector<double> &x, std::vector<double> &y)
{
   int n = A.n;
   y.assign(n, 0.0);

   for (int i = 0; i < n; i++) {
      y[i] += A.di[i] * x[i];

      for (int k = A.ig[i]; k < A.ig[i + 1]; k++) {
         int j = A.jg[k];

         y[i] += A.ggl[k] * x[j];
         y[j] += A.ggu[k] * x[i];
      }
   }
}

void calcLU(const matrix &A, matrix &LU)
{
   LU.n = A.n;
   LU.ig = A.ig;
   LU.jg = A.jg;

   LU.di.resize(A.n);
   LU.ggl = A.ggl;
   LU.ggu = A.ggu;

   for (int i = 0; i < A.n; i++) {
      double sumDi = 0.0;

      int i0 = LU.ig[i];
      int i1 = LU.ig[i + 1];

      for (int k = i0; k < i1; k++) {
         int j = LU.jg[k];

         double sumL = 0.0;
         double sumU = 0.0;

         int p1 = i0;
         int p2 = LU.ig[j];

         while (p1 < k && p2 < LU.ig[j + 1]) {
            int c1 = LU.jg[p1];
            int c2 = LU.jg[p2];

            if (c1 == c2) {
               sumL += LU.ggl[p1] * LU.ggu[p2];
               sumU += LU.ggl[p2] * LU.ggu[p1];

               p1++;
               p2++;
            }
            else if (c1 < c2) {
               p1++;
            }
            else {
               p2++;
            }
         }

         LU.ggl[k] = (LU.ggl[k] - sumL) / LU.di[j];
         LU.ggu[k] = LU.ggu[k] - sumU;

         sumDi += LU.ggl[k] * LU.ggu[k];
      }

      LU.di[i] = A.di[i] - sumDi;

      if (std::fabs(LU.di[i]) < 1e-15) {
         LU.di[i] = 1e-15;
      }
   }
}

void solveL(const matrix &LU, const std::vector<double> &b, std::vector<double> &y)
{
   y = b;

   for (int i = 0; i < LU.n; i++) {
      double sum = 0.0;

      for (int k = LU.ig[i]; k < LU.ig[i + 1]; k++) {
         sum += LU.ggl[k] * y[LU.jg[k]];
      }

      y[i] -= sum;
   }
}

void solveU(const matrix &LU, const std::vector<double> &y, std::vector<double> &x)
{
   x = y;

   for (int i = LU.n - 1; i >= 0; i--) {
      x[i] /= LU.di[i];

      for (int k = LU.ig[i]; k < LU.ig[i + 1]; k++) {
         int j = LU.jg[k];
         x[j] -= LU.ggu[k] * x[i];
      }
   }
}

void doMSG_LU(const matrix &A, msg_lu &m, std::vector<double> &x, const std::vector<double> &b, int maxIter, double eps)
{
   int n = A.n;

   matrix LU;
   calcLU(A, LU);

   m.r.resize(n);
   m.z.resize(n);
   m.t.resize(n);
   m.tmp.resize(n);
   m.x_tilde.resize(n);
   m.v1.resize(n);
   m.v2.resize(n);

   multA(A, x, m.tmp);

   for (int i = 0; i < n; i++) {
      m.r[i] = b[i] - m.tmp[i];
   }

   solveL(LU, m.r, m.tmp);
   solveU(LU, m.tmp, m.z);

   double rho_prev = scalar(m.r, m.z);

   double norm_b = scalar(b, b);
   if (std::fabs(norm_b) < 1e-20) {
      norm_b = 1.0;
   }

   double residual = std::sqrt(scalar(m.r, m.r) / norm_b);

   //std::cout << "Iter 0: Residual = " << residual << "\n";

   if (residual < eps || std::fabs(rho_prev) < 1e-30) {
      return;
   }

   std::vector<double> &p = m.v1;

   for (int i = 0; i < n; i++) {
      p[i] = m.z[i];
   }

   for (int k = 1; k <= maxIter; k++) {
      multA(A, p, m.t);

      double denom = scalar(p, m.t);

      if (std::fabs(denom) < 1e-30) {
         break;
      }

      double alpha = rho_prev / denom;

      for (int i = 0; i < n; i++) {
         x[i] += alpha * p[i];
         m.r[i] -= alpha * m.t[i];
      }

      residual = std::sqrt(scalar(m.r, m.r) / norm_b);

      //if (k % 10 == 0 || residual < eps) {
      //   std::cout << "Iter " << k << ": Residual = " << residual << "\n";
      //}

      if (residual < eps) {
         break;
      }

      solveL(LU, m.r, m.tmp);
      solveU(LU, m.tmp, m.z);

      double rho_new = scalar(m.r, m.z);

      if (std::fabs(rho_prev) < 1e-30) {
         break;
      }

      double beta = rho_new / rho_prev;

      for (int i = 0; i < n; i++) {
         p[i] = m.z[i] + beta * p[i];
      }

      rho_prev = rho_new;
   }
}