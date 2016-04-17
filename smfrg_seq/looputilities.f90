function Qforbidden(Q1,Q2,Q3,Q4)
  logical Qforbidden
  real*8 Q1,Q2,Q3,Q4,Q
  Q=Q1+Q2-Q3-Q4
  Q=abs( mod(Q,2.d0) )
  Qforbidden=(Q>1.e-5)
  return
end function Qforbidden

function Qconserved(Q1,Q2,Q3,Q4)
  logical Qconserved
  real*8 Q1,Q2,Q3,Q4
  logical Qforbidden
  Qconserved=.not.Qforbidden(Q1,Q2,Q3,Q4)
  return
end function Qconserved

