comp := ifort
#opt := -fast # does additional stuff beyond -O3, but needs testing!!!
opt := -O3

# TODO: make use of the Makefile.AlphaProject, but also port Windows stuff there!!!

# MS Windows
ifeq (${OS}, Windows_NT)
  ## see also https://software.intel.com/en-us/compiler_winapp_f (2014-12-03)
  opt := ${opt} -static -Qopenmp-link:static -Qmkl -Qlocation,link,"${VCINSTALLDIR}/bin"
  obj := .obj
  exe := .exe
else
# Linux or Mac OSX
  obj := .o
  exe :=
  opt := ${opt} -mkl -static-intel -openmp-link=static
  uname := ${shell uname}
  # Linux only
  ifeq ($(uname), Linux)
    opt := ${opt} -static -static-libgcc -static-libstdc++
  endif
endif

AlphaBayes: Makefile AlphaBayes.f90 Global${obj} ReadParam${obj} ReadData${obj} InitiateSeed${obj} PearsnR4${obj} \
	RidgeRegression${obj} BayesA${obj} MarkerEffectPostProcessing${obj} gasdev${obj} random_gamma${obj} momentR4${obj} ran1${obj} random_order${obj}
	$(comp) $(opt) AlphaBayes.f90 Global${obj} ReadParam${obj} ReadData${obj} InitiateSeed${obj} PearsnR4${obj} \
		RidgeRegression${obj} BayesA${obj} MarkerEffectPostProcessing${obj} gasdev${obj} random_gamma${obj} momentR4${obj} ran1${obj} random_order${obj} \
		-o AlphaBayes

Global${obj}: Makefile Global.f90
	$(comp) -c $(opt) -o Global${obj} Global.f90

ReadParam${obj}: Makefile ReadParam.f90
	$(comp) -c $(opt) -o ReadParam${obj} ReadParam.f90

ReadData${obj}: Makefile ReadData.f90
	$(comp) -c $(opt) -o ReadData${obj} ReadData.f90

InitiateSeed${obj}: Makefile InitiateSeed.f90
	$(comp) -c $(opt) -o InitiateSeed${obj} InitiateSeed.f90

PearsnR4${obj}: Makefile PearsnR4.f90
	$(comp) -c $(opt) -o PearsnR4${obj} PearsnR4.f90

RidgeRegression${obj}: Makefile RidgeRegression.f90
	$(comp) -c $(opt) -o RidgeRegression${obj} RidgeRegression.f90

BayesA${obj}: Makefile BayesA.f90
	$(comp) -c $(opt) -o BayesA${obj} BayesA.f90

MarkerEffectPostProcessing${obj}: Makefile MarkerEffectPostProcessing.f90
	$(comp) -c $(opt) -o MarkerEffectPostProcessing${obj} MarkerEffectPostProcessing.f90

gasdev${obj}: Makefile gasdev.f90
	$(comp) -c $(opt) -o gasdev${obj} gasdev.f90

random_gamma${obj}: Makefile random_gamma.f90
	$(comp) -c $(opt) -o random_gamma${obj} random_gamma.f90

momentR4${obj}: Makefile momentR4.f90
	$(comp) -c $(opt) -o momentR4${obj} momentR4.f90

ran1${obj}: Makefile ran1.f90
	$(comp) -c $(opt) -o ran1${obj} ran1.f90

random_order${obj}: Makefile random_order.f90
	$(comp) -c $(opt) -o random_order${obj} random_order.f90

clean:
	rm -f *${obj} *.mod

cleanall: clean
	rm -f AlphaBayes${exe}
