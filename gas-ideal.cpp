#include "funciones.h"

int main(int argc, char **argv)
{
  int SEED = 3;
  int N = 10;
  double R = 0.2;
  double L = 10.0;
  int PASOS = 10;
  double TIEMPO = 0.0;
  double DELTA_TIEMPO = 0.1;
  
  std::vector<double> POSICIONES = creacion_posiciones(N, SEED, L);
  std::vector<double> VELOCIDADES = creacion_velocidades(N, SEED);
  
  for(int step = 1; step<= PASOS; step++){
    std::cout<<"t = "<<TIEMPO<<"\n";
    paso(POSICIONES, VELOCIDADES, TIEMPO, DELTA_TIEMPO, R, L);
    print_vector(POSICIONES);
    std::cout<<"\n";
  }
  
  return 0;
}
