run:
	gfortran cgm_complex.f90 -g -fcheck=all -o cgm.exe
	cgm.exe

test:
	gfortran test.f90 -g -fcheck=all -o test.exe
	test.exe

clean:
	del *.mod *.exe *.txt