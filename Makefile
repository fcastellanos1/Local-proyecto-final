ejecutar: a.out gas-ideal.cpp funciones.cpp funciones.h
	./$<

a.out: gas-ideal.cpp funciones.h funciones.cpp
	g++ -g -std=c++11 -fsanitize=leak -fsanitize=address -fsanitize=undefined $< funciones.cpp

clean:
	rm -rf *.x *.out *.x.* *~ *.pdf *.png *.txt *.data *#
