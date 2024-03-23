# Challenge1
Challenge 1 pacs 2024

in questa repo Ã implementata la template function Gradiente (la cui dichiarazione e definizione sono entrambe contenute nell Header file Gradiente.hpp nella cartella include). alla funzione vengono passate by reference:
- la funzione da minimizzare e il suo gradiente (definite nel main)
- una struct di parametri contenete i dati iniziali(x_0, alpha_0) e i 
  parametri per l'implementazione del metodo,
- il vettore soluzione.


I parametri si possono modificare modificando il file data, dal quale vengono presi tramite la funzione GetPot.

la compilazione avviene con in comando: make

la scelta del metodo per la definizione del learning rate puo essere impostata spegificando dalla command line durante l'esecuzione del file main il parametro scelta:
 	./main scelta=1 ---> Exponential decay
	./main scelta=2 ---> Inverse decay
	./main scelta=3 ---> Armijo rule

