#include "funciones.h"

int main(int argc, char **argv)
{
  int SEED = 4;
  int N = std::atoi(argv[1]); //# de partículas
  double TEMPERATURA = std::atof(argv[2]); //kelvin
  double R = 0.1; //metros
  double L = 10.0; //metros
  int PASOS = 1000;
  double TIEMPO = 0.0; //segundos
  double DELTA_TIEMPO = 0.1; //segundos
  
  
  //std::vector<double> POSICIONES = creacion_posiciones(N, SEED, L);
  std::vector<double> POSICIONES = creacion_posiciones_2(N, SEED, L, R);
  std::vector<double> VELOCIDADES = creacion_velocidades(N, SEED);

  gnuplot_init_trayectorias(L);
  
  for(int step = 1; step<= PASOS; step++){
    //std::cout<<"t = "<<TIEMPO<<"\n";
    paso_paralelo(POSICIONES, VELOCIDADES, TIEMPO, DELTA_TIEMPO, R, L);
    //print_vector(POSICIONES);
    gnuplot_trayectorias(POSICIONES, R);
    std::cout<<"\n";
  }
  
  return 0;
}

//Video:
/*
1ra parte = introducción física, con ecuaciones a tratar.
2da parte = implementación computacional.
3ra parte = Dificultades en el proceso.
4ta parte = animación y resultados finales. Confirmación experimental de la ley.
*/
