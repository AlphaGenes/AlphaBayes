comp=ifort
opt=  -zero -nowarn
#opt=  -lpthread -nbs -zero -nowarn
optfor90= -O2  $(opt)   
#optfor90= -g  -traceback $(opt)   
optfor90t= -O2   $(opt)  

weibull: 	  parinclu.o parchange.o check.o lbfgs.o  \
	  kind.o sparse2.o sparse.o sparssub.o fspak90.o fspakwei.o fspaksub.o  dicfs.o \
	  opti.o predict.o rwsave.o various.o  weibmain.o init.o \
	  backsol.o deriv.o deriv2.o 
	$(comp) $(optfor90)  parinclu.o parchange.o weibmain.o  \
	  backsol.o check.o init.o lbfgs.o  \
	  kind.o sparse.o sparse2.o sparssub.o fspak90.o fspakwei.o fspaksub.o dicfs.o \
	  deriv.o deriv2.o opti.o predict.o rwsave.o various.o   
	 mv a.out weibull.exe


weibmain.o:	 weibmain.f
	$(comp) -c $(optfor90) -oweibmain.o weibmain.f 

parinclu.o:	parinclu.f
	$(comp) -c $(optfor90) -o parinclu.o  parinclu.f

parchange.o:	parchange.f
	$(comp) -c $(optfor90) -o parchange.o parchange.f

init.o:	 	init.f 
	$(comp) -c $(optfor90) init.f 

backsol.o:	 backsol.f
	$(comp) -c $(optfor90) backsol.f 

deriv.o:	deriv.f sparse.o 
	$(comp) -c $(optfor90) deriv.f

deriv2.o:	deriv2.f sparse.o
	$(comp) -c  $(optfor90) deriv2.f 

fspakwei.o:	fspakwei.f 
	$(comp) -c $(optfor90) fspakwei.f

dicfs.o:	dicfs.f 
	$(comp) -c $(optfor90) dicfs.f

sparse2.o:	sparse2.f  
	$(comp) -c $(optfor90) sparse2.f

sparssub.o:	sparssub.f
	$(comp) $(optfor90) -c sparssub.f

sparse.o:	sparse.f90
	$(comp) $(optfor90) -c sparse.f90

kind.o:	kind.f90 
	$(comp) -c $(optfor90)  kind.f90

fspak90.o:	fspak90.f90 
	$(comp) -c $(optfor90)  fspak90.f90

fspaksub.o:	fspaksub.f
	$(comp) $(optfor90) -c fspaksub.f

check.o:	check.f
	$(comp) -c $(optfor90) check.f

lbfgs.o:	lbfgs.f  
	$(comp) -c $(optfor90) lbfgs.f

opti.o:		opti.f 
	$(comp) -c $(optfor90) opti.f 

predict.o:	predict.f 
	$(comp) -c $(optfor90) predict.f 

rwsave.o:	rwsave.f 
	$(comp) -c $(optfor90) rwsave.f 

various.o:	various.f 
	$(comp) -c $(optfor90) various.f

clean:
	rm -f *.o *.mod *.lst

