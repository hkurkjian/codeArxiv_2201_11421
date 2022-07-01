MOD=-I../libperso
LIB_PERSO=-L../libperso -lperso

COMP=ifort
BALISES=-qopenmp -check all -backtrace -O3

COMP=gfortran
BALISES=-fopenmp -fcheck=all -fbacktrace -O3

%.eps : %.gp
	gnuplot $<

%.o : %.f90
	$(COMP) $(MOD) $(BALISES) -c $< -o $@

carte: dspecBCS.o carte.o
	$(COMP) $(BALISES) $^ $(LIB_PERSO) -o carte

coupe: dspecBCS.o coupe.o
	$(COMP) $(BALISES) $^ $(LIB_PERSO) -o coupe

spectreBCS: dspecBCS.o spectreBCS.o
	$(COMP) $(BALISES) $^ $(LIB_PERSO) -o spectreBCS

minplasma: dspecBCS.o minplasma.o
	$(COMP) $(BALISES) $^ $(LIB_PERSO) -o minplasma

propre : 
	rm -f *.o 
	rm -f *.mod
	rm -f *.aux
	rm -f *.syntec.gz
	rm -f spectre
	rm -f *.exe

sup : propre
	rm spectreBCS
	rm minplasma
	rm coupe
	rm carte
