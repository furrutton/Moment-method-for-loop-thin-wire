target:
	gfortran -g -fcheck=all cgm_fft.f90 -o debug.exe
	debug.exe
clean:
	del *.mod *.txt *.csv *.exe