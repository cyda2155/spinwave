#ifndef _FUNC_H_
#define _FUNC_H_
#include "headfile_otherh.h"

void func_Hamiltonian(const int, const int, const double, complex<double> *, const double*, double *, const double *, double *, const Atom *, const Atom *);
void func_iteration(const int, const int, const int, const int, const double *, const double, double *, const double *, const Atom *, const Atom *, double *, complex<double> **);
//void func_spin(const int, const complex<double> **, const complex<double> **, const complex<double> **, complex<double> **, complex<double> *, complex<double> *, complex<double> *);
//void func_spin_matrix(const int, const int, const Atom *, const Atom *, complex<double> **, complex<double> **, complex<double> **);
void func_d_bond(const int, double *, const Atom *);

#endif
