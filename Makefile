N_THREADS=16
N_PARTICLES=100
TEMPERATURA=500

gif: gas-ideal.out gas-ideal.cpp funciones.cpp funciones.h
	OMP_NUM_THREADS=$(N_THREADS) ./$< $(N_PARTICLES) $(TEMPERATURA) | gnuplot
	geeqie trayectorias.gif &

test: testing.out testing.cpp funciones.cpp funciones.h
	OMP_NUM_THREADS=$(N_THREADS) ./$<

gas-ideal.out: gas-ideal.cpp funciones.h funciones.cpp
	g++ -g -fopenmp -std=c++11 $< funciones.cpp -o $@

testing.out: testing.cpp funciones.cpp funciones.h
	g++ -g -fopenmp -std=c++11 $< funciones.cpp -o $@

clean:
	rm -rf *.x *.out *.x.* *~ *.pdf *.png *.gif *.txt *.data *#

#-fsanitize=undefined -fsanitize=thread
