#include "funciones.h"

std::vector<double> creacion_posiciones(int n, int & seed, double l)
{
  std::vector<double> posiciones(2*n, 0.0);
  for(int ii = 0; ii<2*n; ii++){
    posiciones[ii] = aleatorio_real(0, l, seed);
  }
  return posiciones;
}

std::vector<double> creacion_posiciones_2(int n, int & seed, double l, double radio)
{  
  int ver_o_fal = 0;
  std::vector<double> posiciones(2*n, 0.0);
  for(int ii = 0; ii<n; ii++){
    while(ver_o_fal == 0){
      posiciones[2*ii] = aleatorio_real(radio, l-radio, seed);
      posiciones[2*ii+1] = aleatorio_real(radio, l-radio, seed);
      for(int jj = 0; jj<n; jj++){
	if(jj==ii && jj!=n-1){jj++;}
	if(jj==ii && jj==n-1){ver_o_fal=1; break;}
	double distancia = std::sqrt(std::pow(posiciones[2*ii]-posiciones[2*jj],2) + std::pow(posiciones[2*ii+1]-posiciones[2*jj+1],2));
	if(distancia<=2*radio){break;}
	else{if (jj==n-1) ver_o_fal=1;}
      }
    }
    ver_o_fal = 0;
  }
  return posiciones;
}

std::vector<double> creacion_velocidades(int n, int & seed)
{
  std::vector<double> velocidades(2*n, 0.0);
#pragma omp parallel for
  for(int ii = 0; ii<2*n-1; ii+=2){
    
    velocidades[ii] = aleatorio_real(-1, 1, seed);
    
    int signo_aleatorio = aleatorio_entero(1,2, seed);
    if(signo_aleatorio == 2){signo_aleatorio = -1;}
    
    velocidades[ii+1] = signo_aleatorio*std::sqrt(1 - velocidades[ii]*velocidades[ii]);
  }
  return velocidades;
}

std::vector<double> creacion_velocidades_2(int n, int & seed, double temperatura, double masa)
{
  double K_B = 1.380649e-23; //Joules/Kelvin
  double v_media = std::sqrt(2*K_B*temperatura/masa); //Según Teo. Equip. Energía
  
  std::vector<double> velocidades(2*n, 0.0);
#pragma omp parallel for
  for(int ii = 0; ii<2*n-1; ii+=2){
    
    velocidades[ii] = aleatorio_real(-1, 1, seed);
    
    int signo_aleatorio = aleatorio_entero(1,2, seed);
    if(signo_aleatorio == 2){signo_aleatorio = -1;}
    
    velocidades[ii+1] = signo_aleatorio*std::sqrt(1 - velocidades[ii]*velocidades[ii]);

    velocidades[ii] *= v_media;
    velocidades[ii+1] *= v_media;
  }
  return velocidades;
}

void paso(std::vector<double> & posiciones, std::vector<double> & velocidades, double & tiempo, double delta_tiempo, double radio, double l)
{
  int n = posiciones.size()/2;

  for(int particula = 0; particula<n; particula++){
    double x = posiciones[2*particula];
    double y = posiciones[2*particula+1];
    
    if(x<=radio || y<=radio || std::fabs(x-l)<=radio || std::fabs(y-l)<=radio){
      momento_con_pared(posiciones, velocidades, particula, delta_tiempo, radio, l);
    }
    
    int particula_2 = 0;
    while(particula_2 <n){
      if(particula_2 == particula){particula_2++;}
      double x_2 = posiciones[2*particula_2];
      double y_2 = posiciones[2*particula_2+1];
      double distancia = std::sqrt(std::pow(x-x_2, 2)+std::pow(y-y_2, 2));
      if(distancia<=2*radio){momento_con_particula(posiciones, velocidades, particula, particula_2, delta_tiempo, radio, l);}
      particula_2++;
    }  
    
    posicion_siguiente(posiciones, velocidades, particula, delta_tiempo, l, radio);
  }
  
  tiempo += delta_tiempo;
}

void paso_paralelo(std::vector<double> & posiciones, std::vector<double> & velocidades, double & tiempo, double delta_tiempo, double radio, double l)
{
  int n = posiciones.size()/2;
  
  std::vector<double> copia_velocidades(2*n, 0.0);
#pragma omp parallel for
  for(int ii=0;ii<2*n;ii++){
    copia_velocidades[ii]=velocidades[ii];
  }
  
#pragma omp parallel for
  for(int particula=0; particula<n; particula++){
    int particula_2 = 0;
    while(particula_2 <n){
      if(particula_2 == particula){particula_2++;}
      if(particula_2 >=n){break;}
      double x = posiciones[2*particula];
      double y = posiciones[2*particula+1];
      double x_2 = posiciones[2*particula_2];
      double y_2 = posiciones[2*particula_2+1];
      double distancia = std::sqrt(std::pow(x-x_2, 2)+std::pow(y-y_2, 2));
      if(distancia<=2*radio){momento_con_particula_paralelo(posiciones, velocidades, copia_velocidades, particula, particula_2, delta_tiempo, radio, l);}
      particula_2++;
    }
  }
  //velocidades = copia modificada
#pragma omp parallel for
  for(int ii=0; ii<2*n; ii++){
    velocidades[ii] = copia_velocidades[ii];
  }

#pragma omp parallel for
  for(int particula=0; particula<n; particula++){
    posicion_siguiente(posiciones, velocidades, particula, delta_tiempo, l, radio);
  }
  
  tiempo += delta_tiempo;
}

void posicion_siguiente(std::vector<double> & posiciones, std::vector<double> & velocidades, int particula, double delta_tiempo, double l, double radio)
{
  int index = 2*particula;
  posiciones[index] += velocidades[index]*delta_tiempo;
  posiciones[index+1] += velocidades[index+1]*delta_tiempo;
  
  double x = posiciones[index];
  double y = posiciones[index+1];
  if(x<=radio || y<=radio || std::fabs(x-l)<=radio || std::fabs(y-l)<=radio){
    momento_con_pared(posiciones, velocidades, particula, delta_tiempo, radio, l);
  }
  
  //Si se sale, hay reflexión con la pared:
  if(x<radio){posiciones[index] = 2*radio-x;}
  if(x>l-radio){posiciones[index] = 2*(l-radio) - x;}
  if(y<radio){posiciones[index+1] = 2*radio-y;}
  if(y>l-radio){posiciones[index+1] = 2*(l-radio) - y;}
  
}

void momento_con_pared(std::vector<double> & posiciones, std::vector<double> & velocidades, int particula, double delta_tiempo, double radio, double l)
{
  int index = 2*particula;
  double x = posiciones[index];
  double y = posiciones[index+1];
  
  if(x<=radio || std::fabs(x-l)<=radio){
    velocidades[index] *= -1;
  }
  
  if(y<=radio || std::fabs(y-l)<=radio){
    velocidades[index+1] *= -1;
  }
  
}



void momento_con_particula(std::vector<double> & posiciones, std::vector<double> & velocidades, int particula, int & particula_2, double delta_tiempo, double radio, double l)
{
  std::vector<double> diferencia_centros(2, 0.0);
  std::vector<double> diferencia_velocidades(2, 0.0);

  diferencia_centros[0] = posiciones[2*particula]-posiciones[2*particula_2];
  diferencia_centros[1] = posiciones[2*particula+1]-posiciones[2*particula_2+1];

  diferencia_velocidades[0] = velocidades[2*particula]-velocidades[2*particula_2];
  diferencia_velocidades[1] = velocidades[2*particula+1]-velocidades[2*particula_2+1];

  double modulo_dc = std::pow(diferencia_centros[0], 2) + std::pow(diferencia_centros[1], 2);
  double prod_punto = diferencia_centros[0]*diferencia_velocidades[0] + diferencia_centros[1]*diferencia_velocidades[1];
  
  velocidades[2*particula] -= (prod_punto/modulo_dc)*diferencia_centros[0];
  velocidades[2*particula+1] -= (prod_punto/modulo_dc)*diferencia_centros[1];
  
  velocidades[2*particula_2] += (prod_punto/modulo_dc)*diferencia_centros[0];
  velocidades[2*particula_2+1] += (prod_punto/modulo_dc)*diferencia_centros[1];
}

void momento_con_particula_paralelo(std::vector<double> & posiciones, std::vector<double> & velocidades, std::vector<double> & copia, int particula, int & particula_2, double delta_tiempo, double radio, double l)
{
  std::vector<double> diferencia_centros(2, 0.0);
  std::vector<double> diferencia_velocidades(2, 0.0);

  diferencia_centros[0] = posiciones[2*particula]-posiciones[2*particula_2];
  diferencia_centros[1] = posiciones[2*particula+1]-posiciones[2*particula_2+1];

  diferencia_velocidades[0] = velocidades[2*particula]-velocidades[2*particula_2];
  diferencia_velocidades[1] = velocidades[2*particula+1]-velocidades[2*particula_2+1];

  double modulo_dc = std::pow(diferencia_centros[0], 2) + std::pow(diferencia_centros[1], 2);
  double prod_punto = diferencia_centros[0]*diferencia_velocidades[0] + diferencia_centros[1]*diferencia_velocidades[1];
  
  copia[2*particula] -= (prod_punto/modulo_dc)*diferencia_centros[0];
  copia[2*particula+1] -= (prod_punto/modulo_dc)*diferencia_centros[1];
  
  //velocidades[2*particula_2] += (prod_punto/modulo_dc)*diferencia_centros[0];
  //velocidades[2*particula_2+1] += (prod_punto/modulo_dc)*diferencia_centros[1];
}

void hacer_distribucion(std::vector<double> & velocidades, double temperatura, double masa, std::string temperatura_string)
{
  int particiones = 15;
  int n = velocidades.size()/2;
  double K_B = 1.380649e-23; //Joules/Kelvin
  double v_media = std::sqrt(2*K_B*temperatura/masa); //metros/segundos
  double v_aprox_max = v_media*2.5;
  
  std::string nombre_archivo = "distribucion_" + temperatura_string + ".gp";
  std::ofstream distribucion;
  distribucion.open (nombre_archivo);
  //distribucion << "set xrange[0:" << v_aprox_max << "] " << std::endl;
  //distribucion << "set yrange[0:" << n << "] " << std::endl;

  std::vector<int> contador(particiones, 0);
  double min_caja = 0.0;
  double max_caja;
  double dv = v_aprox_max/particiones;
  
  for(int particula = 0; particula<n; particula++){
    double vx = velocidades[2*particula];
    double vy = velocidades[2*particula+1];
    double v = std::sqrt(std::pow(vx,2)+std::pow(vy,2));
    caja_donde_cae(v, contador, v_aprox_max);
  }
  
  //Init gnuplot:
  distribucion << "k=" << K_B << "\n";
  distribucion << "T=" << temperatura << "\n";
  distribucion << "N=" << n << "\n";
  distribucion << "m=" << masa << "\n";
  
  distribucion << "set title 'Distribución de velocidades para T="<<temperatura<<"K'\n";
  distribucion << "set grid\n";
  distribucion << "set term pdf; set out 'distribucion_"<<temperatura<<".pdf'\n";
  distribucion << "set xlabel 'Velocidad[m/s]'; set ylabel 'Número de partículas'\n";
  distribucion << "plot 4*pi*(m/(2*pi*k*T))**(1.5)*x**2*exp(-(m*x**2)/(2*k*T))*N*100 title 'Teórica', '-' w l t 'Experimental'\n";
  //distribucion << "plot '-' w l t 'Experimental'\n";
  
  for(int ii=0; ii<particiones; ii++){
    max_caja = dv*(ii+1);
    //Impresión
    distribucion << min_caja << "\t" << 0 << "\n";
    distribucion << min_caja << "\t" << contador[ii] << "\n";
    distribucion << max_caja << "\t" << contador[ii] << "\n";
    min_caja = max_caja;
    if(ii==particiones-1){distribucion << max_caja << "\t" << 0 <<"\n";}
  }

  distribucion << "e\n";
  
  distribucion.close();
}


void caja_donde_cae(double v, std::vector<int> & contador, double vmax)
{
  int particiones = contador.size();
  double min_caja = 0.0;
  double max_caja;
  double dv = vmax/particiones;
  for(int ii=0; ii<particiones; ii++){
    max_caja = dv*(ii+1);
    if(v >= min_caja && v < max_caja){contador[ii]++;}
    if(ii==particiones-1 && v>= max_caja){contador[ii]++;}
    min_caja = max_caja;
  }
}


double aleatorio_real(double min, double max, int & seed)
{
  std::mt19937 gen(seed);
  std::uniform_real_distribution<double> dist(min, max);
  seed++;
  return dist(gen);
}

double aleatorio_entero(int min, int max, int & seed)
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

void gnuplot_init_trayectorias(double l, double temperatura, double dt)
{
  std::cout << "set terminal gif animate " << std::endl;
  std::cout << "set out 'trayectorias_"<<temperatura<<".gif' " << std::endl;
  std::cout << "set title 'Trayectorias de partículas con T="<<temperatura<< "K'" << std::endl;
  std::cout << "set xlabel 'x[m]'; set ylabel 'y[m]'" << std::endl;
  std::cout << "set label 'Tiempo de frame = "<<dt<<"s' at 0,0.03" << std::endl;
  std::cout << "set xrange[0:" << l << "] " << std::endl;
  std::cout << "set yrange[0:" << l << "] " << std::endl;
}

void gnuplot_trayectorias(std::vector<double> & posiciones, double radio)
{
  std::cout << "plot '-' with circles notitle" << std::endl;
  int n = posiciones.size()/2;
  for(int ii = 0; ii < n; ii++){
    std::cout<<posiciones[2*ii]<<"\t"<<posiciones[2*ii+1]<<"\t"<<radio<< std::endl;
  }
  std::cout << "e" <<std::endl;
}
