# librairies de PRIMME
LIBP = -L./primme/ -lprimme
# includes de PRIMME
INCP = -I./primme/PRIMMESRC/COMMONSRC/
# librairies de ARPACK
LIBA = -L./arpack-ng/build/lib -larpack
# includes de ARPACK
INCA = -I./arpack-ng/ICB
# libraires Jadamilu
LIBJ = -L./JADAMILU-main/lib/INT64YGNU -ljadamilu -lamd -lsuitesparseconfig ./JADAMILU-main/lib/INT64YGNU/mc64d.o ./JADAMILU-main/lib/INT64YGNU/mc21d.o ./JADAMILU-main/lib/INT64YGNU/mc64s.o ./JADAMILU-main/lib/INT64YGNU/mc21s.o -lgfortran -lgomp
# toutes les librairies
LIB = $(LIBP) $(LIBA) $(LIBJ) -lm -llapack -lblas

COPT = -O3 -Wall

default: main

clean: 
	rm *.o 
	rm main
	
main: main.o prob.o time.o interface_primme.o print_csr.o norme_residu.o eulerp_vec.o arpack.o CtoFortran.o
	cc $(COPT) $^ -o $@ $(LIB)

main.o: main.c prob.h time.h interface_primme.h print_csr.h norme_residu.h eulerp_vec.h arpack.h CtoFortran.h
	cc $(COPT) -c $< -o $@ $(INCP) $(INCA)

%.o: %.c %.h
	cc $(COPT) -c $< -o $@ $(INCP) $(INCA)
