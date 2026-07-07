#ifndef HEADER_H
#define HEADER_H

#include <vector>
#include <string>
using namespace std;

struct matrix {
   int n;
   vector<int> ig, jg;
   vector<double> di, ggl, ggu;
};

struct msg_lu {
   vector<double> r, z, t, tmp, x_tilde, v1, v2;
};

struct Node {
   double r, z;
};

struct Triangle {
   int n[6];
   int region;
};

struct S1Node {
   int node, funcNum;
};

struct Edge {
   int node1, node2, funcNum;
};

extern vector<Node> nodes;
extern vector<Triangle> triangles;
extern vector<S1Node> s1_nodes;
extern vector<Edge> s2_edges;
extern vector<Edge> s3_edges;

// функции теста
double lambda_f(int region, double r, double z, double t);
double gamma_f(int region, double r, double z, double t);
double f_f(int region, double r, double z, double t);

double u_s1(int funcNum, double r, double z, double t);
double theta_s2(int funcNum, double r, double z, double t);
double beta_s3(int funcNum, double r, double z, double t);
double u_beta_s3(int funcNum, double r, double z, double t);

double u_exact(double r, double z, double t);

// начальные услови€
double u0_f(double r, double z);
double v0_f(double r, double z);
double a0_f(double r, double z);

// решатель
void multA(const matrix &A, const vector<double> &x, vector<double> &y);
void doMSG_LU(const matrix &A, msg_lu &m, vector<double> &x, const vector<double> &b, int maxIter, double eps);

// ћ Ё
void generatePortrait(matrix &A);
void assemble(matrix &G, matrix &M, vector<double> &rs, double t);
void applyNaturalBC(matrix &A, vector<double> &rs, double t);
void applyDirichletBC(matrix &A, vector<double> &b, double t);

void readData(int &maxIter, double &eps);

#endif