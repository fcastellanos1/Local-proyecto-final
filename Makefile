N_THREADS=16
TEMPERATURA=500
#^^^^^^^^^^^^^^Kelvin

all: gas-ideal.out gas-ideal.cpp funciones.cpp funciones.h
	OMP_NUM_THREADS=$(N_THREADS) ./$< $(TEMPERATURA) | gnuplot
	gnuplot distribucion_$(TEMPERATURA).gp
#	geeqie trayectorias_$(TEMPERATURA).gif &
#	xpdf distribucion_$(TEMPERATURA).pdf &

test: testing.out testing.cpp funciones.cpp funciones.h
	OMP_NUM_THREADS=$(N_THREADS) ./$<

gas-ideal.out: gas-ideal.cpp funciones.h funciones.cpp
	g++ -g -fopenmp -std=c++11 $< funciones.cpp -o $@

testing.out: testing.cpp funciones.cpp funciones.h
	g++ -g -fopenmp -std=c++11 $< funciones.cpp -o $@

clean:
	rm -rf *.x *.out *~ *.pdf *.png *.gp *.gif *.data *#

#-fsanitize=undefined -fsanitize=thread
