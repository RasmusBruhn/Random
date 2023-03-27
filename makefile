OPT = -O0
WARNINGS = -Wall -Wextra
LIBRARY = -I../Library -L. -L../Library/libs
OPTIONS = -g
LIBS = -lDebug2 -lRandom

CC = gcc
CCFLAGS = $(WARNINGS) $(LIBRARY) $(LIBS) $(OPTIONS) $(OPT)

ifeq ($(OS),Windows_NT)
	EXEEXT = .exe
	DLLEXT = .dll
else
	EXEEXT = 
	DLLEXT = .o
endif

Random_tb$(EXEEXT): Random_tb.c Random.h libRandom$(DLLEXT)
	$(CC) $(CCFLAGS) Random_tb.c -o Random_tb$(EXEEXT)

libRandom$(DLLEXT): Random.h _Random.h Random.c
	$(CC) $(CCFLAGS) -c Random.c -o libRandom$(DLLEXT)

clean:
	rm -f *.dll *.o *.exe *_tb