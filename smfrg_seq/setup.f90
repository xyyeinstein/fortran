subroutine setupworkplace()
  use workplace
  implicit none

  real*8, dimension (3) :: da,db
  real*8 w

  call dualvector(model%b,model%c,model%a,da); model%ka=da(1:2)*2
  call dualvector(model%c,model%a,model%b,db); model%kb=db(1:2)*2

  print*,'setting up meshes...'
  call meshkvec(); print*,'kmesh done.'
  call meshqvec(); print*,'qmesh done.'
  call setRmesh(); print*,'rmesh done.'
  
  if(usegroup)then; call setgroupindex(); print*,'group index done.'; end if
  
  !if(useX0.and.(.not.appendX0))return

  if(useprojectors)then; call projectors(); print*,'projectors initialized'; end if
  call initV(); print*,'Vertex initialized'
  if(dataready)then 
    open(10,file='Vout.dat',form='UNFORMATTED')
    read(10)w,P,C,D
    close(10)
  end if
   
  return
end subroutine setupworkplace

