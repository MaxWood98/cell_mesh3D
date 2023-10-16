!Cell Mesh 3D Linear Algebra Module
!Max Wood - mw16116@bristol.ac.uk
!Univeristy of Bristol - Department of Aerospace Engineering

!Version 1.0
!Updated 11-10-2023

!Module
module cellmesh3d_linalg_mod
use cellmesh3d_data_mod
contains 


!4x4 Determinant Function ===========================
function determinant_4(M) result(det)
implicit none 

!Result 
real(dp) :: det 

!Variables - Import
real(dp) :: M(4,4)

!Variables - Local 
real(dp) :: det1,det2,det3,det4
real(dp) :: M3(3,3)

!Evaluate sub-matrix determinants 
M3(1,:) = (/M(2,2),M(2,3),M(2,4)/)
M3(2,:) = (/M(3,2),M(3,3),M(3,4)/)
M3(3,:) = (/M(4,2),M(4,3),M(4,4)/)
det1 = determinant_3(M3)
M3(1,:) = (/M(2,1),M(2,3),M(2,4)/)
M3(2,:) = (/M(3,1),M(3,3),M(3,4)/)
M3(3,:) = (/M(4,1),M(4,3),M(4,4)/)
det2 = determinant_3(M3)
M3(1,:) = (/M(2,1),M(2,2),M(2,4)/)
M3(2,:) = (/M(3,1),M(3,2),M(3,4)/)
M3(3,:) = (/M(4,1),M(4,2),M(4,4)/)
det3 = determinant_3(M3)
M3(1,:) = (/M(2,1),M(2,2),M(2,3)/)
M3(2,:) = (/M(3,1),M(3,2),M(3,3)/)
M3(3,:) = (/M(4,1),M(4,2),M(4,3)/)
det4 = determinant_3(M3)

!Determinant 
det = M(1,1)*det1 - M(1,2)*det2 + M(1,3)*det3 - M(1,4)*det4
return 
end function determinant_4




!3x3 Determinant Function ===========================
function determinant_3(M) result(det)
implicit none 

!Result 
real(dp) :: det 

!Variables - Import
real(dp) :: M(3,3)

!Determinant 
det = M(1,1)*(M(2,2)*M(3,3) - M(2,3)*M(3,2)) - M(1,2)*(M(2,1)*M(3,3) - M(2,3)*M(3,1)) + M(1,3)*(M(2,1)*M(3,2) - M(2,2)*M(3,1))
return 
end function determinant_3


end module cellmesh3d_linalg_mod