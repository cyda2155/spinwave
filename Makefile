ifndef cc
	cc = icpc #-xAVX
else
	ifeq ($(cc),gcc)
		cc = icpc  #-xAVX
	endif
endif
CFLAGS = -fopenmp #for GDB: -g
LDFLAGS = -fopenmp

OBJECTS= Read_Input_*.o LapackCPP.o  func_*.o LCAO.o

#LIBES= -lmkl_lapack95_lp64 -lmkl_em64t -lguide 
#LIBES= -lmkl_avx -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_lapack95_lp64 -lmkl_core -liomp5

LIBES=  -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_lapack95_lp64 -lmkl_core -liomp5 
OPTIMIZE= -O3 
INCLUDEPATH= -I/opt/intel/mkl/include
LIBPATH= -L/opt/intel/mkl/lib/intel64
#dependent relations
main : $(OBJECTS)
	$(cc) $(LDFLAGS)  $^ -o $@ $(INCLUDEPATH) $(LIBPATH) $(LIBES) $(OPTIMIZE)
%.o:%.cpp
	$(cc) $(CFLAGS) -c $^ $(INCLUDEPATH) $(LIBPATH) $(OPTIMIZE)
#make clean
.PHONY : cleanall cleanobj cleanmod
clean : cleanobj cleanmod
	$(RM) main
cleanobj : 
	$(RM) *.o
cleanmod :
	$(RM) *.mod
