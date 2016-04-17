#================================================
# system parameters
F90 = gfortran
LIBS = -llapack -lblas
#================================================
# define all obj files
OBJS =  analysis.o \
        bcs_k.o bcseigenmodes.o bubble.o \
        checksymmetry.o contact.o \
        definejob.o \
        eigenstate.o eigentable.o electronphonon.o \
        green.o groupindex.o \
        initv.o interpolation.o \
        logmesh.o looputilities.o \
        meanfield.o MyLAPACK.o \
        pcdeff.o pointgroup.o \
        rpa.o runfrg.o \
        setkmesh.o setqmesh.o setrmesh.o \
        setoform.o setlform.o setmform.o \
        setup.o singularflow.o  \
        smfrg_main.o supercell.o symmetrizepcd.o\
				tmatrix.o \
        model_3bandcuo2.o \
        model_bis2.o model_bis2_1band.o \
				model_compass.o \
				model_dagotto.o \
				model_fs1f.o model_fs1flee.o model_fs2f.o \
				model_hujp.o \
        model_kagome.o \
				model_kfs1l10d.o model_kfs1l5d.o model_kfs1lav.o \
				model_kfs2l20d.o model_kfs2vl.o model_kfsprimitive.o \
				model_kuroki.o \
				model_lhq.o \
        model_openladder.o \
				model_raghu.o model_rice.o \
				model_sq1a1d.o model_sq2a2d.o \
				model_sro1orb.o model_sro2orb.o model_sro3orb.o \
				model_wangfa.o model_afa1l10d.o model_afa1l5d.o \
				model_afa2l5d.o model_hujp2.o model_bfa1l10d.o model_bfa1l5d.o 
MOD  =  smfrg_main.o #$(LIBS) 
#================================================
# link all to generate exe file
smfrg_seq.out: $(OBJS)
	$(F90) $(OBJS) $(LIBS) -o smfrg_seq.out 
#================================================
clean:
	rm -f *.o *.mod *.x 
#generate every obj and module files
smfrg_main.o : smfrg_main.f90
	$(F90) -c smfrg_main.f90
analysis.o : analysis.f90 $(MOD) 
	$(F90) -c analysis.f90 
bcs_k.o : bcs_k.f90 $(MOD)
	$(F90) -c bcs_k.f90
bcseigenmodes.o : bcseigenmodes.f90 $(MOD)
	$(F90) -c bcseigenmodes.f90
bubble.o : bubble.f90 $(MOD)
	$(F90) -c bubble.f90 
checksymmetry.o : checksymmetry.f90 $(MOD)
	$(F90) -c checksymmetry.f90 
contact.o : contact.f90 $(MOD)
	$(F90) -c contact.f90 
definejob.o : definejob.f90 $(MOD)
	$(F90) -c definejob.f90 
eigenstate.o : eigenstate.f90 $(MOD)
	$(F90) -c eigenstate.f90 
eigentable.o : eigentable.f90 $(MOD)
	$(F90) -c eigentable.f90 
electronphonon.o : electronphonon.f90 $(MOD)
	$(F90) -c electronphonon.f90
green.o : green.f90 $(MOD)
	$(F90) -c green.f90
groupindex.o : groupindex.f90 $(MOD)
	$(F90) -c groupindex.f90 
initv.o : initv.f90 $(MOD)
	$(F90) -c initv.f90 
interpolation.o : interpolation.f90 $(MOD)
	$(F90) -c interpolation.f90
logmesh.o : logmesh.f90 $(MOD)
	$(F90) -c logmesh.f90 
looputilities.o : looputilities.f90 $(MOD)
	$(F90) -c looputilities.f90 
MyLAPACK.o : MyLAPACK.f90
	$(F90) -c MyLAPACK.f90
meanfield.o : meanfield.f90
	$(F90) -c meanfield.f90
pcdeff.o : pcdeff.f90 $(MOD) 
	$(F90) -c pcdeff.f90
pointgroup.o : pointgroup.f90 $(MOD)
	$(F90) -c pointgroup.f90
rpa.o : rpa.f90 $(MOD)
	$(F90) -c rpa.f90
runfrg.o : runfrg.f90 $(MOD)
	$(F90) -c runfrg.f90 
setkmesh.o : setkmesh.f90 $(MOD)
	$(F90) -c setkmesh.f90 
setqmesh.o : setqmesh.f90 $(MOD)
	$(F90) -c setqmesh.f90 
setrmesh.o : setrmesh.f90 $(MOD)
	$(F90) -c setrmesh.f90 
setoform.o : setoform.f90 $(MOD)
	$(F90) -c setoform.f90
setlform.o : setlform.f90 $(MOD)
	$(F90) -c setlform.f90
setmform.o : setmform.f90 $(MOD)
	$(F90) -c setmform.f90 
setup.o : setup.f90 $(MOD)
	$(F90) -c setup.f90 
singularflow.o : singularflow.f90 $(MOD)
	$(F90) -c singularflow.f90 
slaveboson.o : slaveboson.f90 $(MOD)
	$(F90) -c slaveboson.f90 
supercell.o : supercell.f90 $(MOD)
	$(F90) -c supercell.f90
symmetrizepcd.o : symmetrizepcd.f90 $(MOD)
	$(F90) -c symmetrizepcd.f90
tmatrix.o : tmatrix.f90 $(MOD)
	$(F90) -c tmatrix.f90
model_3bandcuo2.o : model_3bandcuo2.f90 $(MOD)
	$(F90) -c model_3bandcuo2.f90
model_bis2.o : model_bis2.f90 $(MOD)
	$(F90) -c model_bis2.f90
model_bis2_1band.o : model_bis2_1band.f90 $(MOD)
	$(F90) -c model_bis2_1band.f90
model_compass.o : model_compass.f90 $(MOD)  
	$(F90) -c model_compass.f90
model_dagotto.o : model_dagotto.f90 $(MOD)
	$(F90) -c model_dagotto.f90
model_fs1f.o : model_fs1f.f90 $(MOD)
	$(F90) -c model_fs1f.f90
model_fs1flee.o : model_fs1flee.f90 $(MOD)
	$(F90) -c model_fs1flee.f90
model_fs2f.o : model_fs2f.f90 $(MOD)
	$(F90) -c model_fs2f.f90
model_hujp.o : model_hujp.f90 $(MOD)
	$(F90) -c model_hujp.f90
model_kagome.o : model_kagome.f90 $(MOD)
	$(F90) -c model_kagome.f90
model_kfs1l10d.o : model_kfs1l10d.f90 $(MOD)
	$(F90) -c model_kfs1l10d.f90
model_kfs1l5d.o : model_kfs1l5d.f90 $(MOD)
	$(F90) -c model_kfs1l5d.f90
model_kfs1lav.o : model_kfs1lav.f90 $(MOD)
	$(F90) -c model_kfs1lav.f90
model_kfs2l20d.o : model_kfs2l20d.f90 $(MOD)
	$(F90) -c model_kfs2l20d.f90
model_kfs2vl.o : model_kfs2vl.f90 $(MOD)
	$(F90) -c model_kfs2vl.f90
model_kfsprimitive.o : model_kfsprimitive.f90 $(MOD)
	$(F90) -c model_kfsprimitive.f90
model_kuroki.o : model_kuroki.f90 $(MOD)
	$(F90) -c model_kuroki.f90 
model_lhq.o : model_lhq.f90 $(MOD)
	$(F90) -c model_lhq.f90
model_openladder.o : model_openladder.f90 $(MOD)
	$(F90) -c model_openladder.f90
model_raghu.o : model_raghu.f90 $(MOD)
	$(F90) -c model_raghu.f90 
model_rectangle.o : model_rectangle.f90 $(MOD)
	$(F90) -c model_rectangle.f90
model_rice.o : model_rice.f90 $(MOD)
	$(F90) -c model_rice.f90
model_sq1a1d.o : model_sq1a1d.f90 $(MOD)
	$(F90) -c model_sq1a1d.f90
model_sq2a2d.o : model_sq2a2d.f90 $(MOD)
	$(F90) -c model_sq2a2d.f90
model_sro1orb.o : model_sro1orb.f90 $(MOD)
	$(F90) -c model_sro1orb.f90
model_sro2orb.o : model_sro2orb.f90 $(MOD)
	$(F90) -c model_sro2orb.f90
model_sro3orb.o : model_sro3orb.f90 $(MOD)
	$(F90) -c model_sro3orb.f90
model_wangfa.o : model_wangfa.f90 $(MOD)
	$(F90) -c model_wangfa.f90
model_afa1l10d.o : model_afa1l10d.f90 $(MOD)
	$(F90) -c model_afa1l10d.f90
model_afa1l5d.o : model_afa1l5d.f90 $(MOD)
	$(F90) -c model_afa1l5d.f90
model_afa2l5d.o : model_afa2l5d.f90 $(MOD)
	$(F90) -c model_afa2l5d.f90
model_hujp2.o : model_hujp2.f90 $(MOD)
	$(F90) -c model_hujp2.f90
model_bfa1l10d.o : model_bfa1l10d.f90 $(MOD)
	$(F90) -c model_bfa1l10d.f90 
model_bfa1l5d.o : model_bfa1l5d.f90 $(MOD)
	$(F90) -c model_bfa1l5d.f90 
