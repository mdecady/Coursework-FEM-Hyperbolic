#include <iostream>
#include <fstream>
#include <iomanip>
#include <vector>
#include <cmath>
#include <string>
#include "header.h"

using namespace std;

std::vector<Node> nodes;
std::vector<Triangle> triangles;
std::vector<S1Node> s1_nodes;
std::vector<Edge> s2_edges;
std::vector<Edge> s3_edges;

void printTable(double t, const vector<double> &q)
{
   int N = (int)nodes.size();

   cout << "  таблица для времени t = " << fixed << setprecision(6) << t << "\n";
   cout << string(95, '=') << "\n";
   cout << setw(6) << "Node"
      << setw(10) << "r"
      << setw(10) << "z"
      << setw(18) << "Exact"
      << setw(18) << "Obtained"
      << setw(12) << "Diff" << "\n";

   double maxErr = 0.0;

   for (int i = 0; i < N; i++) {
      double ex = u_exact(nodes[i].r, nodes[i].z, t);
      double sol = q[i];
      double diff = sol - ex;

      if (std::fabs(diff) > maxErr) {
         maxErr = std::fabs(diff);
      }

      cout << fixed << setprecision(6)
         << setw(6) << i
         << setw(10) << nodes[i].r
         << setw(10) << nodes[i].z
         << setprecision(10)
         << setw(18) << ex
         << setw(18) << sol
         << scientific << setprecision(4)
         << setw(12) << diff
         << fixed << "\n";
   }

   cout << string(95, '=') << "\n";
   cout << "Max Absolute Error at t = "
      << setprecision(6) << t
      << ": " << scientific << maxErr << "\n\n";
}

double calcMaxError(double t, const vector<double> &q)
{
   double maxErr = 0.0;

   for (int i = 0; i < (int)nodes.size(); i++) {
      double ex = u_exact(nodes[i].r, nodes[i].z, t);
      double diff = std::fabs(q[i] - ex);

      if (diff > maxErr) {
         maxErr = diff;
      }
   }

   return maxErr;
}

vector<double> generateTimeGrid(double T, int steps, bool nonUniform, double alpha = 0.35)
{
   vector<double> time(steps + 1, 0.0);

   if (!nonUniform) { //равномерная сетка
      for (int i = 0; i <= steps; i++) {
         time[i] = T * i / steps;
      }

      return time;
   }

   vector<double> dtRaw(steps, 0.0);
   double sum = 0.0;

   for (int i = 0; i < steps; i++) { // неравномерная сетка
      double factor;

      if (i % 2 == 0) {
         factor = 1.0 - alpha;
      }
      else {
         factor = 1.0 + alpha;
      }

      dtRaw[i] = factor;
      sum += dtRaw[i];
   }

   double scale = T / sum;

   for (int i = 0; i < steps; i++) {
      time[i + 1] = time[i] + dtRaw[i] * scale;
   }

   time[steps] = T;

   return time;
}

void readData(int &maxIter, double &eps)
{
   std::ifstream in;

   in.open("nodes.txt");
   int nNodes;

   in >> nNodes;

   nodes.resize(nNodes);

   for (int i = 0; i < nNodes; i++) {
      in >> nodes[i].r >> nodes[i].z;
   }

   in.close();

   in.open("triangles.txt");

   int nTri;

   in >> nTri;

   triangles.resize(nTri);

   for (int i = 0; i < nTri; i++) {
      for (int j = 0; j < 6; j++) {
         in >> triangles[i].n[j];
      }

      in >> triangles[i].region;
   }

   in.close();

   auto readBC = [] (const char *f, auto &vec, bool isS1) {
      std::ifstream fin(f);

      if (!fin) {
         return;
      }

      int c;

      fin >> c;

      vec.resize(c);

      for (int i = 0; i < c; i++) {
         if (isS1) {
            fin >> ((S1Node &)vec[i]).node
               >> ((S1Node &)vec[i]).funcNum;
         }
         else {
            fin >> ((Edge &)vec[i]).node1
               >> ((Edge &)vec[i]).node2
               >> ((Edge &)vec[i]).funcNum;
         }
      }

      fin.close();
      };

   readBC("s1_nodes.txt", s1_nodes, true);
   readBC("s2_edges.txt", s2_edges, false);
   readBC("s3_edges.txt", s3_edges, false);

   in.open("solver.txt");

   if (in) {
      in >> maxIter >> eps;
   }
   else {
      maxIter = 1000;
      eps = 1e-13;
   }

   in.close();
}

void setBc1Vector(vector<double> &q, double t)
{
   for (const auto &bc : s1_nodes) {
      int idx = bc.node;

      q[idx] = u_s1(bc.funcNum, nodes[idx].r, nodes[idx].z, t);
   }
}

void setBc1SystemByValues(matrix &A, vector<double> &b, const vector<double> &val)
{
   for (const auto &bc : s1_nodes) {
      int idx = bc.node;
      double v = val[idx];

      A.di[idx] = 1.0;
      b[idx] = v;

      for (int k = A.ig[idx]; k < A.ig[idx + 1]; k++) {
         int j = A.jg[k];

         b[j] -= A.ggu[k] * v;

         A.ggl[k] = 0.0;
         A.ggu[k] = 0.0;
      }

      for (int i = idx + 1; i < A.n; i++) {
         for (int k = A.ig[i]; k < A.ig[i + 1]; k++) {
            if (A.jg[k] == idx) {
               b[i] -= A.ggl[k] * v;

               A.ggl[k] = 0.0;
               A.ggu[k] = 0.0;
            }
         }
      }
   }
}

void firstLayers(double tau, double chi, double sigma, int maxIter, double eps, vector<double> &q0, vector<double> &q1)
{
   int N = (int)nodes.size();

   q0.assign(N, 0.0);
   q1.assign(N, 0.0);

   vector<double> v0(N, 0.0);
   vector<double> a0_bc(N, 0.0);

   for (int i = 0; i < N; i++) {
      double r = nodes[i].r;
      double z = nodes[i].z;

      q0[i] = u0_f(r, z);
      v0[i] = v0_f(r, z);
      a0_bc[i] = a0_f(r, z);
   }

   setBc1Vector(q0, 0.0);

   matrix G0, M0;

   generatePortrait(G0);
   generatePortrait(M0);

   vector<double> b0(N, 0.0);
   vector<double> Gq0(N, 0.0);
   vector<double> Mv0(N, 0.0);
   vector<double> bAcc(N, 0.0);
   vector<double> a0(N, 0.0);

   assemble(G0, M0, b0, 0.0);

   applyNaturalBC(G0, b0, 0.0);

   multA(G0, q0, Gq0);
   multA(M0, v0, Mv0);

   for (int i = 0; i < N; i++) {
      bAcc[i] = (b0[i] - Gq0[i] - sigma * Mv0[i]) / chi;
   }

   setBc1SystemByValues(M0, bAcc, a0_bc);

   msg_lu msg0;

   doMSG_LU(M0, msg0, a0, bAcc, maxIter, eps);

   for (int i = 0; i < N; i++) {
      q1[i] = q0[i] + tau * v0[i] + 0.5 * tau * tau * a0[i];
   }

   setBc1Vector(q1, tau);
}

int main()
{
   setlocale(0,"");
   int maxIter;
   double eps;

   readData(maxIter, eps);

   double T = 1.0;
   double baseTau = 0.00625;
   bool nonUniform = true;
   double nonUniformAlpha = 0.35;

   double chi = 1.0;
   double sigma = 1.0;

   int steps = (int)std::round(T / baseTau);

   vector<double> timeGrid = generateTimeGrid(T, steps, nonUniform, nonUniformAlpha);

   int N = (int)nodes.size();

   vector<double> q_old(N, 0.0);
   vector<double> q_prev(N, 0.0);
   vector<double> q(N, 0.0);

   vector<double> b(N, 0.0);
   vector<double> rs(N, 0.0);

   vector<double> Mq_old(N, 0.0);
   vector<double> Mq_prev(N, 0.0);

   double tauFirst = timeGrid[1] - timeGrid[0];

   firstLayers(tauFirst, chi, sigma, maxIter, eps, q_old, q_prev);

   printTable(0.0, q_old);

   matrix G, M, A;

   generatePortrait(G);
   generatePortrait(M);
   generatePortrait(A);

   msg_lu msg;

   cout << "Starting time loop..." << endl;

   for (int k = 2; k <= steps; k++) {
      double t = timeGrid[k];

      double dt0 = timeGrid[k] - timeGrid[k - 1];
      double dt1 = timeGrid[k - 1] - timeGrid[k - 2];
      double dt = dt0 + dt1;

      double coefA =
         2.0 * chi / (dt * dt0)
         +
         sigma * (dt + dt0) / (dt * dt0);

      double coefPrev =
         2.0 * chi / (dt1 * dt0)
         +
         sigma * dt / (dt1 * dt0);

      double coefOld =
         2.0 * chi / (dt1 * dt)
         +
         sigma * dt0 / (dt1 * dt);

      assemble(G, M, b, t);

      for (int i = 0; i < N; i++) {
         A.di[i] = G.di[i] + coefA * M.di[i];
      }

      for (int i = 0; i < (int)A.jg.size(); i++) {
         A.ggl[i] = G.ggl[i] + coefA * M.ggl[i];
         A.ggu[i] = G.ggu[i] + coefA * M.ggu[i];
      }

      multA(M, q_prev, Mq_prev);
      multA(M, q_old, Mq_old);

      for (int i = 0; i < N; i++) {
         rs[i] = b[i] + coefPrev * Mq_prev[i] - coefOld * Mq_old[i];
      }

      applyNaturalBC(A, rs, t);
      applyDirichletBC(A, rs, t);

      q.assign(N, 0.0);

      doMSG_LU(A, msg, q, rs, maxIter, eps);

      if (std::fabs(t - 0.25) < 1e-7 ||
         std::fabs(t - 0.50) < 1e-7 ||
         std::fabs(t - 0.75) < 1e-7 ||
         std::fabs(t - 1.00) < 1e-7) {
         printTable(t, q);
      }

      q_old = q_prev;
      q_prev = q;
   }

   double finalError = calcMaxError(T, q_prev);
   return 0;
}