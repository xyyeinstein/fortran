module standard_derived_types

  type patchmesh
    integer ntri     
		integer nk                               
    real*8 kf(2),lf                          
    real*8 angle                             
    real*8 kc(2,2),km(2),kb(2,6)             
    real*8, allocatable, dimension (:,:) :: kv        
		real*8, allocatable, dimension (:) :: ek,fk,dS,dL 
  end type patchmesh

  type bandmesh
    integer norb                             
    integer nativeband                       
    integer npatch                           
    character*2 pocket                       
		real*8 mushift                              
		real*8 a(3),b(3),c(3)                        
		real*8 ka(3),kb(3),kc(3)                     
    type (patchmesh), allocatable, dimension (:) :: patch  
  end type bandmesh

  type generation
    integer np                                    
    real*8 w                                      
		real*8 area                                   
    real*8 dp                                   
    integer, allocatable, dimension (:) :: info   
    real*8, allocatable, dimension (:,:) :: p     
	real*8, allocatable, dimension (:,:,:) :: tri     
  end type generation

  type meshconfig 
    integer ng                                         
    integer npoints                                    
    type (generation), allocatable, dimension (:) :: g 
  end type meshconfig

  type Mformconfig
    integer atom,nbatom        
    integer Lform              
    integer Oform              
    character*3 form           
  end type Mformconfig

  type Lformconfig
    character*3 form                             
    integer nr,ndim                              
		real*8, allocatable, dimension (:,:) :: r    
		real*8, allocatable, dimension (:) :: f      
  end type Lformconfig

  type Oformconfig
    character*3 form                                
    character*3 parity                              
    integer npair                                   
    integer, allocatable, dimension (:,:) :: pair   
	  real*8, allocatable, dimension (:) :: f         
  end type Oformconfig

  type modelblock
    integer norb,natom                                
		real*8, allocatable, dimension (:,:) :: ra        
		real*8, dimension (3) :: a,b,c                    
		real*8, dimension (2) :: ka,kb                    
    character*3 group                                 
		integer ng                                         
    character*3, allocatable, dimension (:) :: orbits     
		real*8 U,Vnn,Jnn,JH,Uab,mu,filling                            
		logical zerotz                                       
    character*10 sample                                   
    complex*16, allocatable, dimension (:,:,:,:,:) :: t  
    real*8 t1,t2,t3,split
    logical propergauge            
	  logical periodiclayers         
  	logical ladder                 
    integer nleg                   
  	logical ephcoupled              
		real*8 Veph                    
    integer nest                   
	  real*8 qnest(2,12)             
    integer nband,nkf                                             
    integer, allocatable :: nativeband(:),npatch(:)               
		character*2, allocatable :: pocket(:)                         
	  real*8, allocatable :: kf(:,:),lf(:),pdos(:),mushift(:)       
    integer, allocatable :: band(:)                               
	  real*8 Udd,Upp,Vpd,Vpp,Vdd    
    integer n4logqmesh
  end type modelblock

  type leadingchannel
    real*8 S
    integer iq,is
    real*8 q(2)
    character*3 form
  end type leadingchannel

  
  type outputblock
		real*8 Vmax(4)                        
		real*8 q(2,3)                         
		real*8 dope                           
    type (leadingchannel) :: LeadingPair(10)  
  end type outputblock

  type orderparameter
    integer atom,nbatom
		real*8 bond(2)
    complex*16,allocatable,dimension(:,:) ::  fo
  end type

  type mfblock
    logical slaveboson,dochi,dosc,domag,breakpoint,finished,tunemu
    integer nambu,nbond1atom,nbond,nl,nstep
    real*8 Qsdw(2),nhole,errmin,eta,temp,errmin4mu
    real*8,allocatable,dimension(:,:) :: bonds
    type (orderparameter), allocatable, dimension(:) :: chi,delta,mag
    real*8, allocatable, dimension(:,:,:) :: J
  end type

end module standard_derived_types


module workplace
    use standard_derived_types

    real*8 pi,twopi                                 
    complex*16 one                                  

    integer norb,natom,nambu,norb1atom                              

  integer nLform                                                          
  character*3, allocatable, dimension (:) :: Lformtable                   
  real*8, allocatable, dimension (:,:) :: LformBasis                      
  type (LFormConfig), allocatable, target, dimension (:) :: Lformfunc     

  integer nOform                                                          
  type (OFormconfig), allocatable, target, dimension (:) :: Oformfunc     

  integer ndim                                                            
  logical SkipInterOrbitPair,skipOddOrbPair                               
  type (Mformconfig), allocatable, target, dimension (:) :: Mformfunc     
  type (modelblock), target :: model    
  type (outputblock), target :: output

  integer nq                                      
  real*8, allocatable, dimension (:,:) :: qv      
  type (meshconfig) :: qmesh                      

  integer nk                                      
  real*8, allocatable, dimension (:,:) :: kv      
  type (meshconfig) :: kmesh                      
  logical appendkmesh,appendfs                           
  character*20 kmeshfile                          

  integer nR,Lcontact                             
  integer, allocatable, dimension (:,:) :: R      

  integer nkx,nky                                 
  real*8, allocatable, dimension (:,:,:) :: ek    
  real*8, allocatable, dimension (:,:) :: BdGek
  complex*16,allocatable,dimension(:,:,:) :: BdGAk

  character*3 group,guess                         
  integer ng                                      
  logical usegroup                                
	integer, allocatable :: indexgroup(:,:,:)        
  real*8, allocatable :: Xg(:,:,:)                
  integer, allocatable :: indexqg(:,:),indexkg(:,:)            
  integer, allocatable :: indexRg(:,:,:,:,:)       
  integer nw                     
  real*8 wmax,wmin,wir,wirx0,diverge   

  complex*16, allocatable, dimension (:,:,:) :: P,C,D      
  complex*16, allocatable, dimension (:,:) :: Ploc,Cloc,Dloc

  logical DataReady              
  logical AppendBreakPoint       
  logical QuickSearch            
  logical outputXq                


  logical useprojectors,Appendprojectors  
  character*20 projectorfile              
  integer nP2C,nP2D,nC2P,nC2D,nD2C,nD2P   
  integer, allocatable, dimension (:,:) :: infoP2C,infoP2D,infoC2P,infoC2D,infoD2C,infoD2P  
	complex*16, allocatable, dimension (:) :: P2C,P2D,C2P,C2D,D2C,D2P                         

  logical BCSflow                         
  logical AppendX0,useX0                  
  logical uniformkmesh,uniformqmesh
  
  logical doMF
  type (mfblock) :: MF

end module workplace


!=======================================
program main
use workplace

call definejob()
call setupworkplace()
call runfrg()

end
