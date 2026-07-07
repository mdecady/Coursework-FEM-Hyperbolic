#include "header.h"
#include <cmath>
#include <algorithm>
#include <iostream>

bool findEdgeMidNode(int n1, int n2, int &mid) // поиск среднего узла ребра
{
   for (const auto &t : triangles) {
      int v0 = t.n[0];
      int v1 = t.n[1];
      int v2 = t.n[2];

      if ((n1 == v0 && n2 == v1) || (n1 == v1 && n2 == v0)) {
         mid = t.n[3];
         return true;
      }

      if ((n1 == v1 && n2 == v2) || (n1 == v2 && n2 == v1)) {
         mid = t.n[4];
         return true;
      }

      if ((n1 == v2 && n2 == v0) || (n1 == v0 && n2 == v2)) {
         mid = t.n[5];
         return true;
      }
   }

   return false;
}

void generatePortrait(matrix &A) // генерация портрета
{
   int n = (int)nodes.size();

   A.n = n;

   std::vector<std::vector<int>> adj(n);

   for (const auto &t : triangles) {
      for (int i = 0; i < 6; i++) {
         for (int j = 0; j < i; j++) {
            int r = t.n[i];
            int c = t.n[j];

            if (r == c) {
               continue;
            }

            if (r < c) {
               std::swap(r, c);
            }

            adj[r].push_back(c);
         }
      }
   }

   A.ig.resize(n + 1);
   A.ig[0] = 0;

   for (int i = 0; i < n; i++) {
      auto &v = adj[i];

      std::sort(v.begin(), v.end());
      v.erase(std::unique(v.begin(), v.end()), v.end());

      A.ig[i + 1] = A.ig[i] + (int)v.size();
   }

   int nnz = A.ig[n];

   A.jg.resize(nnz);
   A.di.assign(n, 0.0);
   A.ggl.assign(nnz, 0.0);
   A.ggu.assign(nnz, 0.0);

   int k = 0;

   for (int i = 0; i < n; i++) {
      for (int c : adj[i]) {
         A.jg[k++] = c;
      }
   }
}

void addElement(matrix &A, std::vector<double> &b, const double locM[6][6], const double locB[6], const int idx[6])
{
   for (int i = 0; i < 6; i++) {
      int row = idx[i];

      A.di[row] += locM[i][i];
      b[row] += locB[i];

      for (int j = 0; j < i; j++) {
         int col = idx[j];

         int r = row;
         int c = col;

         if (r < c) {
            std::swap(r, c);
         }

         int pos = A.ig[r];

         while (pos < A.ig[r + 1] && A.jg[pos] != c) {
            pos++;
         }

         if (pos < A.ig[r + 1]) {
            A.ggl[pos] += locM[i][j];
            A.ggu[pos] += locM[j][i];
         }
      }
   }
}

void getShapeFuncs(double L1, double L2, double L3, double N[6], double dN_dL[6][3])
{
   N[0] = L1 * (2.0 * L1 - 1.0);
   N[1] = L2 * (2.0 * L2 - 1.0);
   N[2] = L3 * (2.0 * L3 - 1.0);
   N[3] = 4.0 * L1 * L2;
   N[4] = 4.0 * L2 * L3;
   N[5] = 4.0 * L3 * L1;

   dN_dL[0][0] = 4.0 * L1 - 1.0;
   dN_dL[0][1] = 0.0;
   dN_dL[0][2] = 0.0;

   dN_dL[1][0] = 0.0;
   dN_dL[1][1] = 4.0 * L2 - 1.0;
   dN_dL[1][2] = 0.0;

   dN_dL[2][0] = 0.0;
   dN_dL[2][1] = 0.0;
   dN_dL[2][2] = 4.0 * L3 - 1.0;

   dN_dL[3][0] = 4.0 * L2;
   dN_dL[3][1] = 4.0 * L1;
   dN_dL[3][2] = 0.0;

   dN_dL[4][0] = 0.0;
   dN_dL[4][1] = 4.0 * L3;
   dN_dL[4][2] = 4.0 * L2;

   dN_dL[5][0] = 4.0 * L3;
   dN_dL[5][1] = 0.0;
   dN_dL[5][2] = 4.0 * L1;
}

void assemble(matrix &G, matrix &M, std::vector<double> &rs, double t)
{
   std::fill(G.di.begin(), G.di.end(), 0.0);
   std::fill(G.ggl.begin(), G.ggl.end(), 0.0);
   std::fill(G.ggu.begin(), G.ggu.end(), 0.0);

   std::fill(M.di.begin(), M.di.end(), 0.0);
   std::fill(M.ggl.begin(), M.ggl.end(), 0.0);
   std::fill(M.ggu.begin(), M.ggu.end(), 0.0);

   std::fill(rs.begin(), rs.end(), 0.0);

   const double qa = 0.797426985353087;
   const double qb = 0.101286507323456;

   const double qc = 0.059715871789770;
   const double qd = 0.470142064105115;

   const double w1 = 0.225;
   const double w2 = 0.125939180544827;
   const double w3 = 0.132394152788527;

   const double L_QUAD[7][3] = {
      {1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0},
      {qa, qb, qb},
      {qb, qa, qb},
      {qb, qb, qa},
      {qc, qd, qd},
      {qd, qc, qd},
      {qd, qd, qc}
   };

   const double W_QUAD[7] = {
      w1, w2, w2, w2, w3, w3, w3
   };

   for (const auto &tri : triangles) {
      double locG[6][6] = { 0.0 };
      double locM[6][6] = { 0.0 };
      double locB[6] = { 0.0 };

      double R[6];
      double Z[6];

      for (int i = 0; i < 6; i++) {
         R[i] = nodes[tri.n[i]].r;
         Z[i] = nodes[tri.n[i]].z;
      }

      double detJ =
         (R[1] - R[0]) * (Z[2] - Z[0])
         -
         (R[2] - R[0]) * (Z[1] - Z[0]);

      double area = 0.5 * std::fabs(detJ); // площадь треугольника. 

      double dLdr[3] = {
         (Z[1] - Z[2]) / detJ,
         (Z[2] - Z[0]) / detJ,
         (Z[0] - Z[1]) / detJ
      };

      double dLdz[3] = {
         (R[2] - R[1]) / detJ,
         (R[0] - R[2]) / detJ,
         (R[1] - R[0]) / detJ
      };

      double lamVertex[3]; 

      for (int i = 0; i < 3; i++) {
         lamVertex[i] = lambda_f(tri.region, R[i], Z[i], t);
      }

      for (int q = 0; q < 7; q++) {
         double N[6];
         double dN_dL[6][3];

         getShapeFuncs(L_QUAD[q][0], L_QUAD[q][1], L_QUAD[q][2], N, dN_dL);

         double r = 0.0;
         double z = 0.0;

         for (int i = 0; i < 3; i++) {
            r += L_QUAD[q][i] * R[i];
            z += L_QUAD[q][i] * Z[i];
         }

         double lam =
            L_QUAD[q][0] * lamVertex[0]
            +
            L_QUAD[q][1] * lamVertex[1]
            +
            L_QUAD[q][2] * lamVertex[2];

         double gam = gamma_f(tri.region, r, z, t);
         double fval = f_f(tri.region, r, z, t);

         double weight = area * W_QUAD[q] * r;

         double dNdr[6];
         double dNdz[6];

         for (int i = 0; i < 6; i++) {
            dNdr[i] =
               dN_dL[i][0] * dLdr[0]
               +
               dN_dL[i][1] * dLdr[1]
               +
               dN_dL[i][2] * dLdr[2];

            dNdz[i] =
               dN_dL[i][0] * dLdz[0]
               +
               dN_dL[i][1] * dLdz[1]
               +
               dN_dL[i][2] * dLdz[2];
         }

         for (int i = 0; i < 6; i++) {
            for (int j = 0; j <= i; j++) {
               double valG =
                  (
                     lam * (dNdr[i] * dNdr[j] + dNdz[i] * dNdz[j])
                     +
                     gam * N[i] * N[j]
                     )
                  * weight;

               double valM = N[i] * N[j] * weight;

               locG[i][j] += valG;
               locM[i][j] += valM;

               if (i != j) {
                  locG[j][i] += valG;
                  locM[j][i] += valM;
               }
            }

            locB[i] += fval * N[i] * weight;
         }
      }

      addElement(G, rs, locG, locB, tri.n);

      double zeroB[6] = { 0.0 };
      addElement(M, rs, locM, zeroB, tri.n);
   }
}

void applyNaturalBC(matrix &A, std::vector<double> &rs, double t)
{
   const double gPoints[3] = {
      0.5 - 0.3872983346207417,
      0.5,
      0.5 + 0.3872983346207417
   };

   const double gWeights[3] = {
      5.0 / 18.0,
      8.0 / 18.0,
      5.0 / 18.0
   };

   auto integrateEdge = [&] (const Edge &e, int type) {
      int mid;

      if (!findEdgeMidNode(e.node1, e.node2, mid)) {
         return;
      }

      int idx[3] = {
         e.node1,
         mid,
         e.node2
      };

      double dr = nodes[e.node2].r - nodes[e.node1].r;
      double dz = nodes[e.node2].z - nodes[e.node1].z;

      double len = std::sqrt(dr * dr + dz * dz);

      double locM[3][3] = { 0.0 };
      double locB[3] = { 0.0 };

      for (int q = 0; q < 3; q++) {
         double s = gPoints[q];
         double w = gWeights[q];

         double phi[3] = {
            (1.0 - s) * (1.0 - 2.0 * s),
            4.0 * s * (1.0 - s),
            s * (2.0 * s - 1.0)
         };

         double r = 0.0;
         double z = 0.0;

         for (int i = 0; i < 3; i++) {
            r += phi[i] * nodes[idx[i]].r;
            z += phi[i] * nodes[idx[i]].z;
         }

         double weight = w * len * r;

         if (type == 2) {
            double theta = theta_s2(e.funcNum, r, z, t);

            for (int i = 0; i < 3; i++) {
               locB[i] += theta * phi[i] * weight;
            }
         }
         else {
            double beta = beta_s3(e.funcNum, r, z, t);
            double ub = u_beta_s3(e.funcNum, r, z, t);

            for (int i = 0; i < 3; i++) {
               for (int j = 0; j < 3; j++) {
                  locM[i][j] += beta * phi[i] * phi[j] * weight;
               }

               locB[i] += beta * ub * phi[i] * weight;
            }
         }
      }

      for (int i = 0; i < 3; i++) {
         int row = idx[i];

         A.di[row] += locM[i][i];
         rs[row] += locB[i];

         for (int j = 0; j < i; j++) {
            int col = idx[j];

            int r_idx = row;
            int c_idx = col;

            if (r_idx < c_idx) {
               std::swap(r_idx, c_idx);
            }

            int pos = A.ig[r_idx];

            while (pos < A.ig[r_idx + 1] && A.jg[pos] != c_idx) {
               pos++;
            }

            if (pos < A.ig[r_idx + 1]) {
               A.ggl[pos] += locM[i][j];
               A.ggu[pos] += locM[j][i];
            }
         }
      }
      };

   for (const auto &e : s2_edges) {
      integrateEdge(e, 2);
   }

   for (const auto &e : s3_edges) {
      integrateEdge(e, 3);
   }
}

void applyDirichletBC(matrix &A, std::vector<double> &b, double t)
{
   for (const auto &bc : s1_nodes) {
      int idx = bc.node;

      double val = u_s1(bc.funcNum, nodes[idx].r, nodes[idx].z, t);

      A.di[idx] = 1.0;
      b[idx] = val;

      for (int k = A.ig[idx]; k < A.ig[idx + 1]; k++) {
         int j = A.jg[k];

         b[j] -= A.ggu[k] * val;

         A.ggl[k] = 0.0;
         A.ggu[k] = 0.0;
      }

      for (int i = idx + 1; i < A.n; i++) {
         for (int k = A.ig[i]; k < A.ig[i + 1]; k++) {
            if (A.jg[k] == idx) {
               b[i] -= A.ggl[k] * val;

               A.ggl[k] = 0.0;
               A.ggu[k] = 0.0;
            }
         }
      }
   }
}