comp=ifort
opt= -zero -nowarn
#opt= -lpthread -nbs -zero -nowarn
optfor90= -mkl -O3 $(opt)
#optfor90= -g -traceback $(opt)
uname := ${shell uname}
ifneq ($(uname), Linux)
  optfor90 := $(optfor90) -i-static
endif

AlphaBayes: AlphaBayes.f90 Global.o ReadParam.o ReadData.o InitiateSeed.o PearsnR4.o \
	RidgeRegression.o BayesA.o MarkerEffectPostProcessing.o gasdev.o random_gamma.o momentR4.o ran1.o random_order.o
	$(comp) $(optfor90) AlphaBayes.f90 Global.o ReadParam.o ReadData.o InitiateSeed.o PearsnR4.o \
	RidgeRegression.o BayesA.o MarkerEffectPostProcessing.o gasdev.o random_gamma.o momentR4.o ran1.o random_order.o \
	  -o AlphaBayes

Global.o: Global.f90
	$(comp) -c $(optfor90) -o Global.o Global.f90

ReadParam.o: ReadParam.f90
	$(comp) -c $(optfor90) -o ReadParam.o ReadParam.f90

ReadData.o: ReadData.f90
	$(comp) -c $(optfor90) -o ReadData.o ReadData.f90

InitiateSeed.o: InitiateSeed.f90
	$(comp) -c $(optfor90) -o InitiateSeed.o InitiateSeed.f90

PearsnR4.o: PearsnR4.f90
	$(comp) -c $(optfor90) -o PearsnR4.o PearsnR4.f90

RidgeRegression.o: RidgeRegression.f90
	$(comp) -c $(optfor90) -o RidgeRegression.o RidgeRegression.f90

BayesA.o: BayesA.f90
	$(comp) -c $(optfor90) -o BayesA.o BayesA.f90
	
MarkerEffectPostProcessing.o: MarkerEffectPostProcessing.f90	
	$(comp) -c $(optfor90) -o MarkerEffectPostProcessing.o MarkerEffectPostProcessing.f90
	
gasdev.o: gasdev.f90
	$(comp) -c $(optfor90) -o gasdev.o gasdev.f90

random_gamma.o: random_gamma.f90
	$(comp) -c $(optfor90) -o random_gamma.o random_gamma.f90

momentR4.o: momentR4.f90
	$(comp) -c $(optfor90) -o momentR4.o momentR4.f90

ran1.o: ran1.f90
	$(comp) -c $(optfor90) -o ran1.o ran1.f90

random_order.o: random_order.f90
	$(comp) -c $(optfor90) -o random_order.o random_order.f90

clean:
	rm -f *.o *.mod *.lst
