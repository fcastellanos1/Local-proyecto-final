#include "funciones.h"

Eigen::MatrixXd creacion_particulas(int n, int & seed, double l)
{
  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> matriz;
  matriz.resize(n, 4);
  
  for(int ii = 0; ii<n; ii++){
    matriz(ii, 0) = aleatorio_real(0, l, seed); //x
    matriz(ii, 1) = aleatorio_real(0, l, seed); //y    
    matriz(ii, 2) = aleatorio_real(-1, 1, seed); //vx
    int signo_aleatorio = aleatorio_entero(1,2, seed);
    if(signo_aleatorio == 2){signo_aleatorio = -1;}    
    matriz(ii, 3) = signo_aleatorio*std::sqrt(1 - matriz(ii, 2)*matriz(ii, 2)); //vy
  }

  /* Verificación de orden en memoria de matriz:
  for (int i = 0; i < matriz.size(); i++){
    std::cout << *(matriz.data() + i) << "   ";
    if((i+1)%4 == 0){std::cout<<"\n";}
  }
  std::cout<<"\n"<<matriz<<"\n";
  */
  
  return matriz;
}

std::vector<double> creacion_posiciones(int n, int & seed, double l)
{
  std::vector<double> posiciones(2*n, 0.0);
  for(int ii = 0; ii<2*n; ii++){
    posiciones[ii] = aleatorio_real(0, l, seed);
  }
  return posiciones;
}

std::vector<double> creacion_velocidades(int n, int & seed)
{
  std::vector<double> velocidades(2*n, 0.0);
  for(int ii = 0; ii<2*n-1; ii+=2){
    
    velocidades[ii] = aleatorio_real(-1, 1, seed);
    
    int signo_aleatorio = aleatorio_entero(1,2, seed);
    if(signo_aleatorio == 2){signo_aleatorio = -1;}
    
    velocidades[ii+1] = signo_aleatorio*std::sqrt(1 - velocidades[ii]*velocidades[ii]);
  }
  return velocidades;
}

void paso(std::vector<double> & posiciones, std::vector<double> & velocidades, double & tiempo, double delta_tiempo, double radio, double l)
{
  int n = posiciones.size()/2;
  
  for(int particula = 0; particula<n; particula++){
    int index = 2*particula;
    double x = posiciones[index];
    double y = posiciones[index+1];
    
    if(x>radio && y>radio && std::fabs(x-l)>radio && std::fabs(y-l)>radio){
      sin_colision(posiciones, velocidades, particula, delta_tiempo);
    }
    else{
      con_pared(posiciones, velocidades, particula, delta_tiempo, radio, l);
    }
  }
  
  tiempo += delta_tiempo;
}

void sin_colision(std::vector<double> & posiciones, std::vector<double> & velocidades, int particula, double delta_tiempo)
{
  int index = 2*particula;
  posiciones[index] += velocidades[index]*delta_tiempo;
  posiciones[index+1] += velocidades[index+1]*delta_tiempo;
}

void con_pared(std::vector<double> & posiciones, std::vector<double> & velocidades, int particula, double delta_tiempo, double radio, double l)
{
  int index = 2*particula;
  double x = posiciones[index];
  double y = posiciones[index+1];

  if(x<radio || std::fabs(x-l)<radio || x>l){
    velocidades[index] *= -1;
  }
  
  if(y<radio || std::fabs(y-l)<radio || y>l){
    velocidades[index+1] *= -1;
  }

  posiciones[index] += velocidades[index]*delta_tiempo;
  posiciones[index+1] += velocidades[index+1]*delta_tiempo;
}

double aleatorio_real(double min, double max, int & seed)
{
  std::mt19937 gen(seed);
  std::uniform_real_distribution<double> dist(min, max);
  seed++;
  return dist(gen);
}

double aleatorio_entero(double min, double max, int & seed)
{
  std::mt19937 gen(seed);
  std::uniform_int_distribution<int> dist(min, max);
  seed++;
  return dist(gen);
}

void print_vector(std::vector<double> data)
{
  int n = data.size()/2;
  
  for(int ii = 0; ii<2*n; ii++){
    std::cout<<data[ii]<<"\t";
    if((ii+1)%2 == 0){std::cout<<"\n";}
  }
  
}
