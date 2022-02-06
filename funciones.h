#ifndef FUNCIONES_H_
#define FUNCIONES_H_
#include <iostream>
#include <cassert>
#include <cstdlib>
#include <cmath>
#include <vector>
#include <eigen3/Eigen/Dense>
#include <random>
#include <fstream>
#include <omp.h>

Eigen::MatrixXd creacion_particulas(int n, int & seed, double l);
std::vector<double> creacion_posiciones(int n, int & seed, double l);
std::vector<double> creacion_velocidades(int n, int & seed);
void paso(std::vector<double> & pos, std::vector<double> & vel, double & t, double delta_t, double radio, double l);
void posicion_siguiente(std::vector<double> & posiciones, std::vector<double> & velocidades, int particula, double delta_tiempo);
void momento_con_pared(std::vector<double> & posiciones, std::vector<double> & velocidades, int particula, double delta_tiempo, double radio, double l);
void momento_con_particula(std::vector<double> & posiciones, std::vector<double> & velocidades, int particula, int particula_2, double delta_tiempo, double radio, double l);
double aleatorio_real(double min, double max, int & seed);
double aleatorio_entero(double min, double max, int & seed);
void print_vector(std::vector<double> data);
void gnuplot_init_trayectorias(double l);
void gnuplot_trayectorias(std::vector<double> & posiciones, double radio);

#endif
