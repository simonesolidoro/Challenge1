# Challenge1
Challenge 1 pacs 2024

in questa repo è implementata la template function Gradiente (la cui dichiarazione e definizione sono entrambe contenute nell Header file Gradiente.hpp nella cartella include). alla funzione vengono passate by reference:
- la funzione da minimizzare (definita nel main) e il suo gradiente(possibile scegliere se grad esatto definito in main o grad calcolato con diff.finite centrate)
- una struct di parametri contenete i dati iniziali(x_0, alpha_0) e i parametri per l'implementazione del metodo,
- il vettore soluzione.


I parametri si possono modificare modificando il file data, dal quale vengono presi tramite la funzione GetPot.

la compilazione avviene con il comando: make

la scelta del metodo per la definizione del learning rate puo essere impostata spegificando dalla command line durante l'esecuzione del file main il parametro scelta:

 	./main scelta=1 ---> Exponential decay

	./main scelta=2 ---> Inverse decay

	./main scelta=3 ---> Armijo rule

        scelta di default: Armijo rule  

la scelta del metodo per la definizione del learning rate puo essere impostata specificando dalla command line durante l'esecuzione del file main il parametro gradiente:

	./main gradiente=1 ---> usato gradiente esatto

        ./main gradiente=0 ---> usata aprox con differenze finite centrate e h=0.01 (modifica di h direttamente in definizione di derivfun in Gradiente.hpp)
        di default: gradiente esatto 


ex:
	 make 
    	./main scelta=3 gradiente=1 

