# Compiler flags: all warnings + debugger meta-data
#CC=icc
CC=gcc
# External libraries:
#LIBS = -lgsl -lgslcblas
LIBS = -lgsl -lgslcblas -lm

# This default rule compiles the executable program
instab:instab_o
	./instab.o	

instab_o:instab.c
	$(CC) -o instab.o instab.c $(LIBS)

prof:generate
	./generator.o
generate:nu_anu_profile_generator.c
	$(CC) -o generator.o nu_anu_profile_generator.c -lm

clean:
	rm -r *.o

#Edited by Manu.
