target:
	# mpifort test.f90 -o test.exe
	mpif90 test.f90 -o test.exe
run:
	# mpirun -report-uri - -np 2 ./test.exe
	mpirun -n 8 ./test
clean:
	rm *.exe *.mod
