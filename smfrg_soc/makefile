#================================================
# system parameters
F90 = gfortran
LIBS = -llapack -lblas
OPTS = -O3
#================================================
# define all obj files
OBJS = analysis.o bubble.o antisymmetry.o checksymmetry.o contact.o \
       definejob.o edgestate.o phchannel.o \
       eigenstate.o eigentable.o green.o groupindex.o initV.o \
       Interpolation.o logmesh.o  MyLAPACK.o pointgroup.o runfrg.o \
       setkmesh.o setqmesh.o setMform.o setrmesh.o\
       setup.o singularflow.o smfrg_main.o supercell.o \
       model_cuprate.o model_cuprate2l.o model_lee.o
MOD = smfrg_main.o 
#================================================
# link all to generate exe file
smfrg_soc.out: $(OBJS)
	$(F90) $(OBJS) $(LIBS) -o smfrg_soc.out 
#================================================
clean:
	rm -f *.o *.mod *.x *.a
#================================================
#generate every obj and module files
smfrg_main.o : smfrg_main.f90
	$(F90) $(FCOPTS) -c smfrg_main.f90
antisymmetry.o : antisymmetry.f90 $(MOD)
	$(F90) $(FCOPTS) -c antisymmetry.f90  
analysis.o : analysis.f90 $(MOD) 
	$(F90) $(FCOPTS) -c analysis.f90 
phchannel.o : phchannel.f90 $(MOD)
	$(F90) $(FCOPTS) -c phchannel.f90
bubble.o : bubble.f90 $(MOD)
	$(F90) $(FCOPTS) -c bubble.f90 
checksymmetry.o : checksymmetry.f90 $(MOD)
	$(F90) $(FCOPTS) -c checksymmetry.f90 
contact.o : contact.f90 $(MOD)
	$(F90) $(FCOPTS) -c contact.f90 
definejob.o : definejob.f90 $(MOD)
	$(F90) $(FCOPTS) -c definejob.f90 
edgestate.o : edgestate.f90 $(MOD)
	$(F90) $(FCOPTS) -c edgestate.f90 
eigenstate.o : eigenstate.f90 $(MOD)
	$(F90) $(FCOPTS) -c eigenstate.f90 
eigentable.o : eigentable.f90 $(MOD)
	$(F90) $(FCOPTS) -c eigentable.f90 
green.o : green.f90 $(MOD)
	$(F90) $(FCOPTS) -c green.f90
groupindex.o : groupindex.f90 $(MOD)
	$(F90) $(FCOPTS) -c groupindex.f90 
initV.o : initV.f90 $(MOD)
	$(F90) $(FCOPTS) -c initV.f90 
Interpolation.o : Interpolation.f90 $(MOD)
	$(F90) $(FCOPTS) -c Interpolation.f90
logmesh.o : logmesh.f90 $(MOD)
	$(F90) $(FCOPTS) -c logmesh.f90 
MyLAPACK.o : MyLAPACK.f90
	$(F90) $(FCOPTS) -c MyLAPACK.f90 
pointgroup.o : pointgroup.f90 $(MOD)
	$(F90) $(FCOPTS) -c pointgroup.f90 
runfrg.o : runfrg.f90 $(MOD)
	$(F90) $(FCOPTS) -c runfrg.f90 
setkmesh.o : setkmesh.f90 $(MOD)
	$(F90) $(FCOPTS) -c setkmesh.f90 
setqmesh.o : setqmesh.f90 $(MOD)
	$(F90) $(FCOPTS) -c setqmesh.f90 
setrmesh.o : setrmesh.f90 $(MOD)
	$(F90) $(FCOPTS) -c setrmesh.f90 
setMform.o : setMform.f90 $(MOD)
	$(F90) $(FCOPTS) -c setMform.f90 
setup.o : setup.f90 $(MOD)
	$(F90) $(FCOPTS) -c setup.f90 
singularflow.o : singularflow.f90 $(MOD)
	$(F90) $(FCOPTS) -c singularflow.f90 
supercell.o : supercell.f90 $(MOD)
	$(F90) $(FCOPTS) -c supercell.f90 
model_cuprate.o : model_cuprate.f90 $(MOD)
	$(F90) $(FCOPTS) -c model_cuprate.f90
model_cuprate2l.o : model_cuprate2l.f90 $(MOD)
	$(F90) $(FCOPTS) -c model_cuprate2l.f90
model_lee.o : model_lee.f90 $(MOD)
	$(F90) $(FCOPTS) -c model_lee.f90

























