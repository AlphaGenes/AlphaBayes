ifort Global.f90 -c -O3 -i-static -O3 -m64
ifort ReadData.f90 -c -O3 -i-static -O3 -m64
ifort ReadParam.f90 -c -O3 -i-static -O3 -m64
ifort InitiateSeed.f90 -c -O3 -i-static -O3 -m64
ifort PearsnR4.f90 -c -O3 -i-static -O3 -m64
ifort RidgeRegression.f90 -c -O3 -i-static -O3 -m64
ifort gasdev.f90 -c -O3 -i-static -O3 -m64
ifort momentR4.f90 -c -O3 -i-static -O3 -m64
ifort ran1.f90 -c -O3 -i-static -O3 -m64
ifort random_order.f90 -c -O3 -i-static -O3 -m64
ifort AlphaBayes.f90 -c -O3 -i-static -O3 -m64

ifort *.o -mkl -o AlphaBayes
