FC	=ifort -O3 -c
#FC	= gfortran -O3 -c
#FC	 = f77 -O3 -c
#FC	 = g77  -c
linkF=ifort -o
#linkF=gfortran  -o
#linkF=f77 -static -o
#librF= -llapack -lblas

checkfiles= pseudogapBFKM-NESS.f readin.f utilities_psgBFKM.f selfenergies4.f nullstelle.f\
	    nca_utilities.f  dscal.f broyden.f 
obj= pseudogapBFKM-NESS.o readin.o utilities_psgBFKM.o selfenergies4.o nullstelle.o\
	    nca_utilities.o  dscal.o broyden.o
.f.o:
	${FC} $<

nca:	$(obj)
	$(linkF) $@ $(obj) $(librF)
