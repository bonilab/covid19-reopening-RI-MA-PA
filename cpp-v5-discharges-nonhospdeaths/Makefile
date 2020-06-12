# driver: driver.o rkf.o
# 	lcc  driver.o rkf.o -o driver -lm
# makefile
# Tabs *must* be used for the indentations below;
# spaces cause make syntax errors.

.PHONY: clean

CC=g++
#CFLAGS=-fast -xO4 -xdepend -xarch=v8plusa -xprefetch -xvector -xunroll=8 -fsimple=2 -xsafe=mem
#LIBS=-lm -lgsl -lgslcblas
LIBS=-lm

odesim: $(wildcard *.cpp *.h)
	$(CC) -O3 -c -o rkf.o rkf.cpp 							$(LIBS)
	$(CC) -O3 -c -o derivs.o derivs.cpp 						$(LIBS)
	$(CC) -O3 -c -o driver.o driver.cpp 						$(LIBS)
	$(CC) -O3 -c -o parseargs.o parseargs.cpp 					$(LIBS)
	$(CC) -O3 -c -o generate_trajectories.o generate_trajectories.cpp 		$(LIBS)
	$(CC) -O3 -c -o prms.o prms.cpp 						$(LIBS)
	$(CC) -O3 -o $@ driver.o derivs.o rkf.o generate_trajectories.o prms.o parseargs.o	$(LIBS)

clean:
	rm -f odesim *.o *~

