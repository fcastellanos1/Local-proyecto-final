N_THREADS=16
N_PARTICLES=50
TEMPERATURA=500

gif: a.out gas-ideal.cpp funciones.cpp funciones.h
	OMP_NUM_THREADS=$(N_THREADS) ./$< $(N_PARTICLES) $(TEMPERATURA) | gnuplot
	geeqie trayectorias.gif &

a.out: gas-ideal.cpp funciones.h funciones.cpp
	g++ -g -fopenmp -std=c++11 $< funciones.cpp

clean:
	rm -rf *.x *.out *.x.* *~ *.pdf *.png *.gif *.txt *.data *#

#-fsanitize=undefined -fsanitize=thread
