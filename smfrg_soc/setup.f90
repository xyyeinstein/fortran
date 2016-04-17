subroutine setupworkplace()
  use workplace
  implicit none

  !call eigentable(); print*, 'ek saved to file'
  call meshqvec(); print*,'qmesh done.'
  call meshkvec(); print*,'kmesh done.'
  call setRmesh(); print*,'rmesh done.'
  
  if(usegroup)then; call setgroupindex(); print*,'group index done.'; end if

  if(.not.dataready)call initV(); print*,'Vertex initialized'
  if(useprojectors)then; call projectors(); print*,'projectors initialized'; end if

  return
end subroutine setupworkplace

