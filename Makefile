gif: a.out gas-ideal.cpp funciones.cpp funciones.h
	OMP_NUM_THREADS=1 ./$< | gnuplot

a.out: gas-ideal.cpp funciones.h funciones.cpp
	g++ -g -fopenmp -std=c++11 -fsanitize=undefined -fsanitize=thread $< funciones.cpp

clean:
	rm -rf *.x *.out *.x.* *~ *.pdf *.png *.txt *.data *#
