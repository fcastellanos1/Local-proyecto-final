#include "funciones.h"

int main(int argc, char **argv)
{
  //Test para ver si las funciones aleatorias funcionan correctamente:
  int semilla = 10;
  int n=10000;
  int media_int = 0;
  double media_double = 0.0;
  int min_int = 2;
  int max_int = 13;
  double min_double = 4.1;
  double max_double = 17.3;
  for(int ii=0;ii<n;ii++){
    media_int += aleatorio_entero(min_int, max_int, semilla);
    media_double += aleatorio_real(min_double, max_double, semilla);
  }
  std::cout<<"Media esperada para enteros: "<<(max_int - min_int)/2.0 + min_int<<" media test: "<<(media_int + 0.0)/n<<"\n";
  std::cout<<"Media esperada para reales: "<<(max_double - min_double)/2.0 + min_double<<" media test: "<<media_double/n<<"\n";
  //Resultado: Verdadero. Acción: implementar en código.

    
  //Test para probar si el sanitizer siempre pelea con
  //el "omp parallel for".
  /*
  int n = 100;
  std::vector<double> velocidades(n, 0.0);
  std::vector<double> copia(n, 1.0);
#pragma omp parallel for
  for(int ii=0; ii<n; ii++){
    copia[ii] = velocidades[ii];
  }
  for(int jj=0; jj<n; jj++){
    std::cout<<copia[jj];
    if(jj==n-1){std::cout<<"\n";}
  }
  */
  //Resultado: Verdadero. Acción: no correr con sanitizers.
    
    return 0;
}
