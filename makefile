Random Random_tb.c Random.h Random.o:
	gcc -I../Library -g Random_tb.c -o Random_tb -lm

Random.o Random.c Random.h:
	gcc -I../Library -g Random.c -o Random.o