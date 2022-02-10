#include "funciones.h"

int main(int argc, char **argv)
{
  int SEED = 1;
  int N = 100; //# de partículas
  double TEMPERATURA = std::atof(argv[1]); //kelvin
  double MASA = 0.032/6.02214199e23; //kilogramos
  double R = 0.005; //metros
  double L = 1.0; //metros
  int PASOS = 100;
  double TIEMPO = 0.0; //segundos
  double DELTA_TIEMPO = 0.5e-5; //segundos
  std::string temperatura_string = argv[1];
  
  //std::vector<double> POSICIONES = creacion_posiciones(N, SEED, L);
  std::vector<double> POSICIONES = creacion_posiciones_2(N, SEED, L, R);
  //std::vector<double> VELOCIDADES = creacion_velocidades(N, SEED);
  std::vector<double> VELOCIDADES = creacion_velocidades_2(N, SEED, TEMPERATURA, MASA);
  
  //gnuplot_init_trayectorias_gp(L, temperatura_string);
  gnuplot_init_trayectorias(L);
  
  for(int step = 1; step<= PASOS; step++){
    //std::cout<<"t = "<<TIEMPO<<"\n";
    
    //paso(POSICIONES, VELOCIDADES, TIEMPO, DELTA_TIEMPO, R, L);
    paso_paralelo(POSICIONES, VELOCIDADES, TIEMPO, DELTA_TIEMPO, R, L);
    
    //print_vector(POSICIONES);
    gnuplot_trayectorias(POSICIONES, R);
    //gnuplot_trayectorias_gp(POSICIONES, R);
    std::cout<<"\n";
  }

  //trayectorias.close();
  
  return 0;
}

//Video:
/*
1ra parte = introducción física, con ecuaciones a tratar.
2da parte = implementación computacional.
3ra parte = Dificultades en el proceso.
4ta parte = animación y resultados finales. Confirmación experimental de la ley.
*/
