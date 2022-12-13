FC	=ifort -O3 -c
#FC	= gfortran -O3 -c -fcheck=all, -ffpe-trap=invalid,zero,overflow -g -Wall -Wextra -pedantic
#FC	 = f77 -O3 -c
#FC	 = g77  -c
linkF=ifort -o
#linkF=gfortran  -o
#linkF=f77 -static -o
#librF= -llapack -lblas

checkfiles= pseudogapBFKM-NESS.f readin.f utilities_psgBFKM.f selfenergies4.f \
	      dscal.f broyden.f 
obj= pseudogapBFKM-NESS.o readin.o utilities_psgBFKM.o selfenergies4.o \
	      dscal.o broyden.o
.f.o:
	${FC} $<

nca:	$(obj)
	$(linkF) $@ $(obj) $(librF)
