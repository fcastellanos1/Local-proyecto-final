#ifndef FUNCIONES_H_
#define FUNCIONES_H_
#include <iostream>
#include <cstdlib>
#include <cmath>
#include <vector>
#include <random>
#include <fstream>
#include <string>
#include <omp.h>

std::vector<double> creacion_posiciones(int n, int & seed, double l);
std::vector<double> creacion_posiciones_2(int n, int & seed, double l, double radio);
std::vector<double> creacion_velocidades(int n, int & seed);
std::vector<double> creacion_velocidades_2(int n, int & seed, double temperatura, double masa);
void paso(std::vector<double> & pos, std::vector<double> & vel, double & t, double delta_t, double radio, double l);
void paso_paralelo(std::vector<double> & posiciones, std::vector<double> & velocidades, double & tiempo, double delta_tiempo, double radio, double l);
void posicion_siguiente(std::vector<double> & posiciones, std::vector<double> & velocidades, int particula, double delta_tiempo, double l, double radio);
void momento_con_pared(std::vector<double> & posiciones, std::vector<double> & velocidades, int particula, double delta_tiempo, double radio, double l);
void momento_con_particula(std::vector<double> & posiciones, std::vector<double> & velocidades, int particula, int & particula_2, double delta_tiempo, double radio, double l);
void momento_con_particula_paralelo(std::vector<double> & posiciones, std::vector<double> & velocidades, std::vector<double> & copia, int particula, int & particula_2, double delta_tiempo, double radio, double l);
void hacer_distribucion(std::vector<double> & velocidades, double temperatura, double masa, std::string temperatura_string);
void caja_donde_cae(double v, std::vector<int> & contador, double vmax);
double aleatorio_real(double min, double max, int & seed);
double aleatorio_entero(int min, int max, int & seed);
void print_vector(std::vector<double> data);
void gnuplot_init_trayectorias(double l, double temperatura, double dt);
void gnuplot_trayectorias(std::vector<double> & posiciones, double radio);

#endif
