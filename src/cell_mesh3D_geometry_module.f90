!2D and 3D Geometry Routine Module
!Max Wood - mw16116@bristol.ac.uk
!Univeristy of Bristol - Department of Aerospace Engineering

!Version 10.0
!Updated 24-10-2023

!Geometry subroutines module
module cellmesh3d_geometry_mod
use cellmesh3d_linalg_mod
use ieee_arithmetic!, only: ieee_value,IEEE_QUIET_NAN
contains


!Set comparison precision function ===========================
function set_cmp_prc(reflen,min_precision_val) result(cval)
implicit none 

!Variables - Import 
real(dp) :: reflen,min_precision_val,cval

!Set comparision zero bound to the maximum of the reference length and the round off error minimum precision bound 
cval = max(reflen,min_precision_val)
return 
end function set_cmp_prc




!Cell volumes subroutine ===========================
subroutine get_cell_volumes(Cvol,volume_mesh)
implicit none 

!Variables - Import
real(dp), dimension(:), allocatable :: Cvol
type(vol_mesh_data) :: volume_mesh

!Variables - Local 
integer(in) :: ff
real(dp) :: Aface,Vface
real(dp) :: Nf(3),Nfn(3),Xpf(3)

!Initialise
if (allocated(Cvol)) then 
    deallocate(Cvol)
end if
allocate(Cvol(volume_mesh%ncell))
Cvol(:) = 0.0d0

!Accumulate volumes 
do ff=1,volume_mesh%nface

    !Normal vector of this face
    Nf = newell_normal(volume_mesh%faces(ff)%nvtx,volume_mesh%faces(ff)%vertices,volume_mesh%vtx)

    !Area of this face
    Aface = norm2(Nf)

    !Construct volume contribution of this face
    if (Aface .NE. 0.0d0) then 
        Nfn = Nf/Aface
    else
        !print *, 'zface',ff
        Nfn = 0.0d0 
    end if 
    Xpf(:) = volume_mesh%vtx(volume_mesh%faces(ff)%vertices(1),:)
    Vface = arbitrary_face_vol(Aface,Xpf,Nfn)

    !Accumulate to adjacent cells 
    if (volume_mesh%faces(ff)%cleft .GT. 0) then 
        Cvol(volume_mesh%faces(ff)%cleft) = Cvol(volume_mesh%faces(ff)%cleft) + Vface
    end if 
    if (volume_mesh%faces(ff)%cright .GT. 0) then 
        Cvol(volume_mesh%faces(ff)%cright) = Cvol(volume_mesh%faces(ff)%cright) - Vface
    end if 
end do
return 
end subroutine get_cell_volumes




!Ray-Plane Intersection Subroutine (Moller-Trombore Method) ======================
subroutine ray_plane_intersection(Xi,Yi,Zi,u,v,ti,O,D,V1,V2,V3)
implicit none
!Intersection of a ray with a plane defined by three points
!Assumes the ray will intersect the plane
!If the ray is paralell to the plane then M = inf

!Variables - Import
real(dp) :: V1(3),V2(3),V3(3)
real(dp) :: O(3),D(3)
real(dp) :: Xi,Yi,Zi
real(dp) :: u,v,ti

!Variables - Local 
real(dp) :: M,Mden
real(dp) :: E1(3),E2(3),T(3),P(3),Q(3)

!Construct parameters
E1(:) = V2(:) - V1(:)
E2(:) = V3(:) - V1(:)
T(:) = O(:) - V1(:)
P(1) = D(2)*E2(3) - D(3)*E2(2)
P(2) = -(D(1)*E2(3) - D(3)*E2(1))
P(3) = D(1)*E2(2) - D(2)*E2(1)
Q(1) = T(2)*E1(3) - T(3)*E1(2)
Q(2) = -(T(1)*E1(3) - T(3)*E1(1))
Q(3) = T(1)*E1(2) - T(2)*E1(1)

!Identify intersection position on plane
Mden = P(1)*E1(1) + P(2)*E1(2) + P(3)*E1(3)
if (Mden .NE. 0.0d0) then 
    M = 1.0d0/(P(1)*E1(1) + P(2)*E1(2) + P(3)*E1(3))
    u = M*(P(1)*T(1) + P(2)*T(2) + P(3)*T(3))
    v = M*(Q(1)*D(1) + Q(2)*D(2) + Q(3)*D(3))
    ti = M*(Q(1)*E2(1) + Q(2)*E2(2) + Q(3)*E2(3))
    Xi = (1.0d0 - u - v)*V1(1) + u*V2(1) + v*V3(1)
    Yi = (1.0d0 - u - v)*V1(2) + u*V2(2) + v*V3(2)
    Zi = (1.0d0 - u - v)*V1(3) + u*V2(3) + v*V3(3)
else
    u = ieee_value(1.0d0,IEEE_QUIET_NAN)
    v = ieee_value(1.0d0,IEEE_QUIET_NAN)
    ti = ieee_value(1.0d0,IEEE_QUIET_NAN)
    Xi = ieee_value(1.0d0,IEEE_QUIET_NAN)
    Yi = ieee_value(1.0d0,IEEE_QUIET_NAN)
    Zi = ieee_value(1.0d0,IEEE_QUIET_NAN)
end if 
return
end subroutine ray_plane_intersection




!Ray-Triangle Intersection Subroutine (watertight WOOP method) ======================
subroutine ray_triangle_intersection(intri,iscoplanar,xi,yi,zi,u,v,w,t,Us,Vs,Ws,org,dir,vt1,vt2,vt3)
implicit none 

!Variables - Import
integer(in) :: intri,iscoplanar
real(dp) :: xi,yi,zi,u,v,w,t,Us,Vs,Ws
real(dp) :: org(3),dir(3),vt1(3),vt2(3),vt3(3)

!Variables - Local 
integer(in) :: kx,ky,kz,kxt,kyt
real(dp) :: Ax,Ay,Az,Bx,By,Bz,Cx,Cy,Cz
real(dp) :: Sx,Sy,Sz,Ts,Det,rcpDet
real(dp) :: A(3),B(3),C(3)

!Find maximal ray direction and perform coordinate swap to make this z
kz = maxloc(abs(dir),1)
kx = kz + 1
if (kx == 4) then 
    kx = 1
end if 
ky = kx + 1
if (ky == 4) then 
    ky = 1
end if 

!Perform winding direction preserving swap of x and y coordinates
if (dir(kz) .LT. 0) then 
    kxt = kx
    kyt = ky
    kx = kyt
    ky = kxt
end if 

!Find shear constants 
Sx = dir(kx)/dir(kz)
Sy = dir(ky)/dir(kz)
Sz = 1.0d0/dir(kz)

!Triangle vertices relative to the ray origin 
A = vt1 - org
B = vt2 - org
C = vt3 - org

!Perform shear and scale of triangle vertices 
Ax = A(kx) - Sx*A(kz)
Ay = A(ky) - Sy*A(kz)
Bx = B(kx) - Sx*B(kz)
By = B(ky) - Sy*B(kz)
Cx = C(kx) - Sx*C(kz)
Cy = C(ky) - Sy*C(kz)

!Find sclaed barycentric coordinates of the intersection location 
Us = Cx*By - Cy*Bx
Vs = Ax*Cy - Ay*Cx
Ws = Bx*Ay - By*Ax

!Calculate determinant 
Det = Us + Vs + Ws

!Check for triangle containement 
if (((Us .LT. 0.0d0) .OR. (Vs .LT. 0.0d0) .OR. (Ws .LT. 0.0d0)) .AND. &
((Us .GT. 0.0d0) .OR. (Vs .GT. 0.0d0) .OR. (Ws .GT. 0.0d0))) then 
    intri = 0
else
    intri = 1
end if 

!Flag if coplanar 
if (Det == 0.0d0) then 
    iscoplanar = 1
else
    iscoplanar = 0
end if 

!Calculate actual properties if not coplanar
if (iscoplanar == 0) then  

    !Find scaled z coordinates for hit distance 
    Az = Sz*A(kz)
    Bz = Sz*B(kz)
    Cz = Sz*C(kz)
    Ts = Us*Az + Vs*Bz + Ws*Cz
    
    !Normalise to actual barycentric coordinates and hit distance 
    rcpDet = 1.0d0/Det
    u = Us*rcpDet
    v = Vs*rcpDet
    w = Ws*rcpDet
    t = Ts*rcpDet
    
    !Find actual hit location 
    xi = u*vt1(1) + v*vt2(1) + w*vt3(1)
    yi = u*vt1(2) + v*vt2(2) + w*vt3(2)
    zi = u*vt1(3) + v*vt2(3) + w*vt3(3)
else
    
    !Set nan return 
    u = ieee_value(1.0d0,IEEE_QUIET_NAN)
    v = ieee_value(1.0d0,IEEE_QUIET_NAN)
    w = ieee_value(1.0d0,IEEE_QUIET_NAN)
    t = ieee_value(1.0d0,IEEE_QUIET_NAN)
    xi = ieee_value(1.0d0,IEEE_QUIET_NAN)
    yi = ieee_value(1.0d0,IEEE_QUIET_NAN)
    zi = ieee_value(1.0d0,IEEE_QUIET_NAN)
end if 
return 
end subroutine ray_triangle_intersection




!Line-triangle intersection bool function ===========================
function line_tri_intersect_bool(vl1,vl2,vt1,vt2,vt3) result(int_bool)
implicit none 

!Variables - Import
integer(in) :: int_bool
real(dp) :: vl1(3),vl2(3),vt1(3),vt2(3),vt3(3)

!Variables - Local 
character(len=2) :: int_location_on_tri
integer(in) :: in_triangle,in_line 

!Classify intersection state with the triangle 
int_location_on_tri = line_tri_intersect_state(vl1,vl2,vt1,vt2,vt3)
if (int_location_on_tri .NE. 'ot') then 
    in_triangle = 1
else
    in_triangle = 0
end if 

!Check if within line 
in_line = line_tri_span_bool(vl1,vl2,vt1,vt2,vt3)

!Set insersection state 
if ((in_triangle == 1) .AND. (in_line == 1)) then 
    int_bool = 1
else
    int_bool = 0 
end if 
return 
end function line_tri_intersect_bool




!Line-triangle intersection state function ===========================
function line_tri_intersect_state(vl1,vl2,vt1,vt2,vt3) result(vtx_loc)
implicit none 

!Variables - Import
character(len=2) :: vtx_loc
real(dp) :: vl1(3),vl2(3),vt1(3),vt2(3),vt3(3)

!Variables - Local 
integer(in) :: intri,iscoplanar
real(dp) :: xi,yi,zi,u,v,w,t,Us,Vs,Ws
real(dp) :: org(3),dir(3)

!Set ray parameters 
org(:) = vl1(:)
dir(:) = vl2(:) - vl1(:)

!Find intersection
call ray_triangle_intersection(intri,iscoplanar,xi,yi,zi,u,v,w,t,Us,Vs,Ws,org,dir,vt1,vt2,vt3)

!Classify location 
if ((intri == 1) .AND. (iscoplanar == 0))then 
    vtx_loc = 'in'
else
    vtx_loc = 'ot'
end if 
return 
end function line_tri_intersect_state




!Line triangle span bool test function ===========================
function line_tri_span_bool(vl1,vl2,vt1,vt2,vt3) result(span_bool)
implicit none 

!Variables - Import
integer(in) :: span_bool
real(dp) :: vl1(3),vl2(3),vt1(3),vt2(3),vt3(3)

!Variables - Local 
character(len=1) span_state

!Find the span state of the line and triangle 
span_state = line_tri_span_state(vl1,vl2,vt1,vt2,vt3)

!Set the state
if ((span_state == 'c') .OR. (span_state == 't')) then !crosses or touches the plane of the triangle 
    span_bool = 1
else !is co-planer and co-incident with or does not cross the plane of the triangle 
    span_bool = 0
end if 
return 
end function line_tri_span_bool




!Line triangle span state function ===========================
function line_tri_span_state(vl1,vl2,vt1,vt2,vt3) result(span_state)
implicit none 

!Variables - Import
character(len=1) span_state
real(dp) :: vl1(3),vl2(3),vt1(3),vt2(3),vt3(3)

!Variables - Local 
integer(in) :: Nvzero
real(dp) :: vol1,vol2

!Volumes of the two line-triangle tetraheda
vol1 = tetrahedra_volume_x6(vt1,vt2,vt3,vl1)
vol2 = tetrahedra_volume_x6(vt1,vt2,vt3,vl2)

!Find number of zero volumes 
Nvzero = 0 
if (abs(vol1) == 0.0d0) then 
    Nvzero = Nvzero + 1
end if
if (abs(vol2) == 0.0d0) then 
    Nvzero = Nvzero + 1
end if

!Set span state 
if (sign(1.0d0,vol1) .NE. sign(1.0d0,vol2)) then !proper cross of plane 
    span_state = 'c'
else
    if (Nvzero == 0) then !volumes are non-zero and the same sign so the line does not cross the plane 
        span_state = 'n'
    elseif (Nvzero == 1) then !one volume is zero so touches the plane 
        span_state = 't'
    else !both volumes are zero so the line is coplanar with the triangle 
        span_state = 'p'
    end if 
end if 
return 
end function line_tri_span_state




!Line triangle intersect type function ===========================
function line_tri_intersect_type(vl1,vl2,vt1,vt2,vt3) result(int_type)
implicit none 

!Variables - Import
character(len=2) int_type
real(dp) :: vl1(3),vl2(3),vt1(3),vt2(3),vt3(3)

!Variables - Local 
integer(in) :: Nvzero
integer(in) :: iszero(2)
real(dp) :: vol1,vol2

!Volumes of the two line-triangle tetraheda
vol1 = tetrahedra_volume_x6(vt1,vt2,vt3,vl1)
vol2 = tetrahedra_volume_x6(vt1,vt2,vt3,vl2)

!Find number of zero volumes 
Nvzero = 0 
iszero(:) = 0 
if (abs(vol1) == 0.0d0) then 
    Nvzero = Nvzero + 1
    iszero(1) = 1
end if
if (abs(vol2) == 0.0d0) then 
    Nvzero = Nvzero + 1
    iszero(2) = 1
end if

!Set span state (entry if vol1 <= 0 and vol2 >= 0)
if (Nvzero == 0) then !both volumes are non zero 
    if (sign(1.0d0,vol1) .NE. sign(1.0d0,vol2)) then !proper cross of plane 
        if (sign(1.0d0,vol1) .LT. 0.0d0) then !in proper
            int_type = 'ip'
        else !out proper
            int_type = 'op'
        end if 
    else !no intersect
        int_type = 'nc'
    end if 
else !improper cross 
    if (Nvzero == 1) then !one volume is zero so touches the plane 
        if (iszero(1) == 1) then !vol1 = 0 
            if (sign(1.0d0,vol2) .GT. 0.0d0) then !in touch
                int_type = 'it'
            else !out touch
                int_type = 'ot'
            end if 
        else !vol2 = 0 
            if (sign(1.0d0,vol1) .LT. 0.0d0) then !in touch
                int_type = 'it'
            else !out touch
                int_type = 'ot'
            end if 
        end if 
    else !both volumes are zero so the line is coplanar with the triangle 
        int_type = 'cp'
    end if 
end if
return 
end function line_tri_intersect_type




!Line-triangle intersection Function ===========================
function line_tri_intersect(vl1,vl2,vt1,vt2,vt3) result(vint)
implicit none 

!Variables - Import
real(dp) :: vint(3)
real(dp) :: vl1(3),vl2(3),vt1(3),vt2(3),vt3(3)

!Variables - Local 
integer(in) :: intri,iscoplanar
real(dp) :: xi,yi,zi,u,v,w,t,Us,Vs,Ws
real(dp) :: org(3),dir(3)

!Set ray parameters 
org(:) = vl1(:)
dir(:) = vl2(:) - vl1(:)

!Find intersection
call ray_triangle_intersection(intri,iscoplanar,xi,yi,zi,u,v,w,t,Us,Vs,Ws,org,dir,vt1,vt2,vt3)

!Store
vint(1) = xi
vint(2) = yi
vint(3) = zi

!Bound within line 
if (t .GT. 1.0d0) then 
    vint(:) = vl2(:)
elseif (t .LT. 0.0d0) then 
    vint(:) = vl1(:)
end if 
return 
end function line_tri_intersect




!Barycentric Line-triangle intersection Function ===========================
function baryc_line_tri_intersect(vl1,vl2,vt1,vt2,vt3) result(bcint)
implicit none 

!Variables - Import
real(dp) :: bcint(3)
real(dp) :: vl1(3),vl2(3),vt1(3),vt2(3),vt3(3)

!Variables - Local 
integer(in) :: intri,iscoplanar
real(dp) :: xi,yi,zi,u,v,w,t,Us,Vs,Ws
real(dp) :: org(3),dir(3)

!Set ray parameters 
org(:) = vl1(:)
dir(:) = vl2(:) - vl1(:)

!Find intersection
call ray_triangle_intersection(intri,iscoplanar,xi,yi,zi,u,v,w,t,Us,Vs,Ws,org,dir,vt1,vt2,vt3)

!Store
bcint(1) = u
bcint(2) = v
bcint(3) = w
return 
end function baryc_line_tri_intersect




!Function to classify a vertices location WRT a triangle using barycentric coordinates ===========================
function vtx_bary_tri_loc(u,v,w,ztol) result(vtx_loc)
implicit none 

!Variables - Import
character(len=2) :: vtx_loc
real(dp) :: u,v,w,ztol

!Variables - Local 
real(dp) :: ub,lb

!Set upper and lower bounds 
ub = 1.0d0 + ztol
lb = -ztol

!Determine location (bias to vertices when nearby for safety)
if ((u <= ub) .AND. (v <= ub) .AND. (w <= ub) .AND. (u >= lb) .AND. (v >= lb) .AND. (w >= lb)) then !within bounds of internal
    if ((u >= ztol) .AND. (v <= 2.0d0*ztol) .AND. (w <= 2.0d0*ztol)) then !vertex 1 (in tollerance)
        vtx_loc = 'v1'
    elseif ((u <= 2.0d0*ztol) .AND. (v >= ztol) .AND. (w <= 2.0d0*ztol)) then !vertex 2 (in tollerance)
        vtx_loc = 'v2'
    elseif ((u <= 2.0d0*ztol) .AND. (v <= 2.0d0*ztol) .AND. (w >= ztol)) then !vertex 3 (in tollerance)
        vtx_loc = 'v3'
    elseif ((u >= ztol) .AND. (v >= ztol) .AND. (w <= ztol)) then !edge (v1-v2) (in tollerance)
        vtx_loc = 'e1'
    elseif ((u <= ztol) .AND. (v >= ztol) .AND. (w >= ztol)) then !edge (v2-v3) (in tollerance)
        vtx_loc = 'e2'
    elseif ((u >= ztol) .AND. (v <= ztol) .AND. (w >= ztol)) then !edge (v3-v1) (in tollerance)
        vtx_loc = 'e3'
    elseif ((u >= ztol) .AND. (v >= ztol) .AND. (w >= ztol)) then !inside triangle (in tollerance)
    ! elseif ((u >= 0.0d0) .AND. (v >= 0.0d0) .AND. (w >= 0.0d0)) then !inside triangle (exact)
        vtx_loc = 'in'
    elseif ((u <= -ztol) .OR. (v <= -ztol) .OR. (w <= -ztol)) then !outside triangle (in tollerance)
    ! elseif ((u < 0.0d0) .OR. (v < 0.0d0) .OR. (w < 0.0d0)) then !outside triangle (exact)
        vtx_loc = 'ot'
    else
        vtx_loc = 'uc'
    endif 
else 
    vtx_loc = 'ot'
end if 
return 
end function vtx_bary_tri_loc




!Find barycentric coordinates function ===========================
function get_barycentric_coordinates(vp,vt1,vt2,vt3) result(barycp)
implicit none 

!Variables - Import
real(dp) :: vp(3),vt1(3),vt2(3),vt3(3),barycp(3)

!Variables - Local 
real(dp) :: M0,M1,M2,M3,D
real(dp) :: A(3),B(3),C(3)

!Base triangle vectors 
A = vt2 - vt1
B = vt3 - vt1
C = vp - vt1

!Matrix entries
M0 = dot_product(A,A)
M1 = dot_product(B,A)
M2 = M1
M3 = dot_product(B,B)

!Determinant 
D = M0*M3 - M1*M2

!Coordinates (u,v,w)
barycp(1) = (dot_product(C,A)*M3 - M1*dot_product(C,B))/D
barycp(2) = (M0*dot_product(C,B) - dot_product(C,A)*M2)/D
barycp(3) = 1.0d0 - barycp(1) - barycp(2)
return 
end function get_barycentric_coordinates




!Project to barycentric coordinates function ===========================
function project_to_barycentric(vp,vt1,vt2,vt3) result(bcint)
implicit none 

!Variables - Import
real(dp) :: bcint(3)
real(dp) :: vp(3),vt1(3),vt2(3),vt3(3)

!Variables - Local 
integer(in) :: intri,iscoplanar
real(dp) :: xi,yi,zi,u,v,w,t,Us,Vs,Ws
real(dp) :: org(3),dir(3),ntri(3)

!Triangle normal 
ntri = crossp3(vt2-vt1,vt3-vt1)

!Set ray parameters 
org(:) = vp(:)
dir(:) = ntri(:)

!Find intersection
call ray_triangle_intersection(intri,iscoplanar,xi,yi,zi,u,v,w,t,Us,Vs,Ws,org,dir,vt1,vt2,vt3)

!Store
bcint(1) = u
bcint(2) = v
bcint(3) = w
return 
end function project_to_barycentric




!Minimum distance point to triangle function ===========================
function min_dist_point_to_tri(vp,vt1,vt2,vt3) result(vid)
implicit none 

!Variables - Import
real(dp) :: vp(3),vt1(3),vt2(3),vt3(3),vid(4)

!Variables - Local 
integer(in) :: within_tri,mindedge,intri,iscoplanar
real(dp) :: xi,yi,zi,u,v,w,t,Us,Vs,Ws
real(dp) :: O(3),D(3),Nt(3),vid_e1(4),vid_e2(4),vid_e3(4),edists(3)

!Triangle normal 
Nt = crossp3(vt2-vt1,vt3-vt1)

!Set ray parameters 
O(:) = vp(:)
D(:) = Nt(:)

!Find intersection with triangle plane 
call ray_triangle_intersection(intri,iscoplanar,xi,yi,zi,u,v,w,t,Us,Vs,Ws,O,D,vt1,vt2,vt3)
within_tri = intri

!Cases on triangle containment 
if (within_tri == 1) then !return point and distance to point 
    vid(1) = Xi 
    vid(2) = Yi 
    vid(3) = Zi 
    vid(4) = norm2(vid(1:3)-vp(:))
else !check for closest point on edges 
    vid_e1 = min_dist_point_to_edge(vt1,vt2,vp)
    vid_e2 = min_dist_point_to_edge(vt2,vt3,vp)
    vid_e3 = min_dist_point_to_edge(vt3,vt1,vp)
    edists(1) = vid_e1(4)
    edists(2) = vid_e2(4)
    edists(3) = vid_e3(4)
    mindedge = minloc(edists,1)
    if (mindedge == 1) then 
        vid(:) = vid_e1(:)
    elseif (mindedge == 2) then 
        vid(:) = vid_e2(:)
    elseif (mindedge == 3) then 
        vid(:) = vid_e3(:)
    end if 
end if 
return 
end function min_dist_point_to_tri




!Minimum distant point to edge function ===========================
function min_dist_point_to_edge(ve1,ve2,vp) result(vid)
implicit none 

!Variables - Import
real(dp) :: vp(3),ve1(3),ve2(3),vid(4)

!Variables - Local 
real(dp) :: t,f
real(dp) :: edir(3)

!Evaluate closest point 
edir(:) = ve2(:) - ve1(:)
edir(:) = edir(:)/norm2(edir(:))
t = dot_product(vp - ve1,edir)
vid(1:3) = ve1(:) + t*edir(:)
f = dot_product(vid(1:3) - ve1(:),edir(:))/norm2(ve2(:) - ve1(:))
if (f .GT. 1.0d0) then 
    vid(1:3) = ve2(:)
elseif (f .LT. 0.0d0) then 
    vid(1:3) = ve1(:)
end if 
vid(4) = norm2(vid(1:3) - vp(:))
return 
end function min_dist_point_to_edge




!Line-square intersection bool Function===========================
function line_square_intersect_bool(vl1,vl2,vs1,vs2,vs3,vs4) result(int_bool)
implicit none 

!Variables - Import
integer(in) :: int_bool
real(dp) :: vl1(3),vl2(3),vs1(3),vs2(3),vs3(3),vs4(3)

!Variables - Local 
integer(in) :: int_tri1,int_tri2

!Test triangle 1
int_tri1 = line_tri_intersect_bool(vl1,vl2,vs1,vs2,vs3)

!Test triangle 2
int_tri2 = line_tri_intersect_bool(vl1,vl2,vs1,vs3,vs4)

!Set intersection state
if ((int_tri1  == 1) .OR. (int_tri2 == 1)) then 
    int_bool = 1
else
    int_bool = 0 
end if 
return 
end function line_square_intersect_bool




!Triangle - Cube Intersection Checking Function ===========================
function tri_cube_intersect_bool(vc1,vc2,vc3,vc4,vc5,vc6,vc7,vc8,vt1,vt2,vt3) result(int_type)
implicit none 

!Variables - Import
integer(in) :: int_type
real(dp) :: vc1(3),vc2(3),vc3(3),vc4(3),vc5(3),vc6(3),vc7(3),vc8(3),vt1(3),vt2(3),vt3(3)

!Variables - Local 
integer(in) :: int_e1,int_e2,int_e3,int_e4,int_e5,int_e6,int_e7,int_e8,int_e9,int_e10,int_e11,int_e12
integer(in) :: int_t_f1,int_t_f2,int_t_f3,int_t_f4,int_t_f5,int_t_f6
integer(in) :: tri_e1_cint,tri_e2_cint,tri_e3_cint,tv1_cont,tv2_cont,tv3_cont,cube_int_tri
real(dp) :: xmin,xmax,ymin,ymax,zmin,zmax

!Test intersection of each edge of the triangle with the cube 
!e1: vt1 -> vt2
int_t_f1 = line_square_intersect_bool(vt1,vt2,vc1,vc2,vc3,vc4)
int_t_f2 = line_square_intersect_bool(vt1,vt2,vc1,vc5,vc6,vc2)
int_t_f3 = line_square_intersect_bool(vt1,vt2,vc2,vc6,vc7,vc3)
int_t_f4 = line_square_intersect_bool(vt1,vt2,vc4,vc8,vc7,vc3)
int_t_f5 = line_square_intersect_bool(vt1,vt2,vc1,vc5,vc8,vc4)
int_t_f6 = line_square_intersect_bool(vt1,vt2,vc5,vc6,vc7,vc8)
if (max(int_t_f1,int_t_f2,int_t_f3,int_t_f4,int_t_f5,int_t_f6) == 1) then 
    tri_e1_cint = 1
else
    tri_e1_cint = 0 
end if 

!e2: vt2 -> vt3
int_t_f1 = line_square_intersect_bool(vt2,vt3,vc1,vc2,vc3,vc4)
int_t_f2 = line_square_intersect_bool(vt2,vt3,vc1,vc5,vc6,vc2)
int_t_f3 = line_square_intersect_bool(vt2,vt3,vc2,vc6,vc7,vc3)
int_t_f4 = line_square_intersect_bool(vt2,vt3,vc4,vc8,vc7,vc3)
int_t_f5 = line_square_intersect_bool(vt2,vt3,vc1,vc5,vc8,vc4)
int_t_f6 = line_square_intersect_bool(vt2,vt3,vc5,vc6,vc7,vc8)
if (max(int_t_f1,int_t_f2,int_t_f3,int_t_f4,int_t_f5,int_t_f6) == 1) then 
    tri_e2_cint = 1
else
    tri_e2_cint = 0 
end if 

!e3: vt3 -> vt1
int_t_f1 = line_square_intersect_bool(vt3,vt1,vc1,vc2,vc3,vc4)
int_t_f2 = line_square_intersect_bool(vt3,vt1,vc1,vc5,vc6,vc2)
int_t_f3 = line_square_intersect_bool(vt3,vt1,vc2,vc6,vc7,vc3)
int_t_f4 = line_square_intersect_bool(vt3,vt1,vc4,vc8,vc7,vc3)
int_t_f5 = line_square_intersect_bool(vt3,vt1,vc1,vc5,vc8,vc4)
int_t_f6 = line_square_intersect_bool(vt3,vt1,vc5,vc6,vc7,vc8)
if (max(int_t_f1,int_t_f2,int_t_f3,int_t_f4,int_t_f5,int_t_f6) == 1) then 
    tri_e3_cint = 1
else
    tri_e3_cint = 0 
end if 

!Test intersection of each edge of the cube with the triangle 
int_e1 = line_tri_intersect_bool(vc1,vc2,vt1,vt2,vt3)
int_e2 = line_tri_intersect_bool(vc2,vc3,vt1,vt2,vt3)
int_e3 = line_tri_intersect_bool(vc3,vc4,vt1,vt2,vt3)
int_e4 = line_tri_intersect_bool(vc4,vc1,vt1,vt2,vt3)
int_e5 = line_tri_intersect_bool(vc1,vc5,vt1,vt2,vt3)
int_e6 = line_tri_intersect_bool(vc2,vc6,vt1,vt2,vt3)
int_e7 = line_tri_intersect_bool(vc3,vc7,vt1,vt2,vt3)
int_e8 = line_tri_intersect_bool(vc4,vc8,vt1,vt2,vt3)
int_e9 = line_tri_intersect_bool(vc5,vc6,vt1,vt2,vt3)
int_e10 = line_tri_intersect_bool(vc6,vc7,vt1,vt2,vt3)
int_e11 = line_tri_intersect_bool(vc7,vc8,vt1,vt2,vt3)
int_e12 = line_tri_intersect_bool(vc8,vc5,vt1,vt2,vt3)
if (max(int_e1,int_e2,int_e3,int_e4,int_e5,int_e6,int_e7,int_e8,int_e9,int_e10,int_e11,int_e12) == 1) then 
    cube_int_tri = 1 
else
    cube_int_tri = 0 
end if 

!Test containment of each vertex of the triangle within the cube 
xmin = min(vc1(1),vc2(1),vc3(1),vc4(1),vc5(1),vc6(1),vc7(1),vc8(1))
xmax = max(vc1(1),vc2(1),vc3(1),vc4(1),vc5(1),vc6(1),vc7(1),vc8(1)) 
ymin = min(vc1(2),vc2(2),vc3(2),vc4(2),vc5(2),vc6(2),vc7(2),vc8(2)) 
ymax = max(vc1(2),vc2(2),vc3(2),vc4(2),vc5(2),vc6(2),vc7(2),vc8(2))  
zmin = min(vc1(3),vc2(3),vc3(3),vc4(3),vc5(3),vc6(3),vc7(3),vc8(3))  
zmax = max(vc1(3),vc2(3),vc3(3),vc4(3),vc5(3),vc6(3),vc7(3),vc8(3))   
tv1_cont = 0
if ((vt1(1) .GE. xmin) .AND. (vt1(1) .LE. xmax)) then 
    if ((vt1(2) .GE. ymin) .AND. (vt1(2) .LE. ymax)) then 
        if ((vt1(3) .GE. zmin) .AND. (vt1(3) .LE. zmax)) then 
            tv1_cont = 1
        end if 
    end if 
end if 
tv2_cont = 0
if ((vt2(1) .GE. xmin) .AND. (vt2(1) .LE. xmax)) then 
    if ((vt2(2) .GE. ymin) .AND. (vt2(2) .LE. ymax)) then 
        if ((vt2(3) .GE. zmin) .AND. (vt2(3) .LE. zmax)) then 
            tv2_cont = 1
        end if 
    end if 
end if 
tv3_cont = 0
if ((vt3(1) .GE. xmin) .AND. (vt3(1) .LE. xmax)) then 
    if ((vt3(2) .GE. ymin) .AND. (vt3(2) .LE. ymax)) then 
        if ((vt3(3) .GE. zmin) .AND. (vt3(3) .LE. zmax)) then 
            tv3_cont = 1
        end if 
    end if 
end if 

!Set intersection status 
if (max(tri_e1_cint,tri_e2_cint,tri_e3_cint,tv1_cont,tv2_cont,tv3_cont,cube_int_tri) == 1) then 
    int_type = 1
else
    int_type = 0 
end if 
return 
end function tri_cube_intersect_bool




!Triangle - Cube Edge Intersection Checking Function ===========================
function tri_cubeedge_intersect_bool(vc1,vc2,vc3,vc4,vc5,vc6,vc7,vc8,vt1,vt2,vt3) result(int_type)
implicit none 

!Variables - Import
integer(in) :: int_type
real(dp) :: vc1(3),vc2(3),vc3(3),vc4(3),vc5(3),vc6(3),vc7(3),vc8(3),vt1(3),vt2(3),vt3(3)

!Variables - Local 
integer(in) :: int_e1,int_e2,int_e3,int_e4,int_e5,int_e6,int_e7,int_e8,int_e9,int_e10,int_e11,int_e12

!Test intersection of each edge of the cube with the triangle 
int_e1 = line_tri_intersect_bool(vc1,vc2,vt1,vt2,vt3)
int_e2 = line_tri_intersect_bool(vc2,vc3,vt1,vt2,vt3)
int_e3 = line_tri_intersect_bool(vc3,vc4,vt1,vt2,vt3)
int_e4 = line_tri_intersect_bool(vc4,vc1,vt1,vt2,vt3)
int_e5 = line_tri_intersect_bool(vc1,vc5,vt1,vt2,vt3)
int_e6 = line_tri_intersect_bool(vc2,vc6,vt1,vt2,vt3)
int_e7 = line_tri_intersect_bool(vc3,vc7,vt1,vt2,vt3)
int_e8 = line_tri_intersect_bool(vc4,vc8,vt1,vt2,vt3)
int_e9 = line_tri_intersect_bool(vc5,vc6,vt1,vt2,vt3)
int_e10 = line_tri_intersect_bool(vc6,vc7,vt1,vt2,vt3)
int_e11 = line_tri_intersect_bool(vc7,vc8,vt1,vt2,vt3)
int_e12 = line_tri_intersect_bool(vc8,vc5,vt1,vt2,vt3)

!Set intersection status 
if (max(int_e1,int_e2,int_e3,int_e4,int_e5,int_e6,int_e7,int_e8,int_e9,int_e10,int_e11,int_e12) == 1) then 
    int_type = 1 
else
    int_type = 0 
end if 
return 
end function tri_cubeedge_intersect_bool




!3D Cross Product Subroutine ===========================
function crossp3(a,b) result(n)
implicit none

!Variables - Import
real(dp) :: n(3),a(3),b(3)

!Calculate
n(1) = a(2)*b(3) - a(3)*b(2)
n(2) = -(a(1)*b(3) - a(3)*b(1))
n(3) = a(1)*b(2) - a(2)*b(1)
return
end function crossp3
    


!Signed volume of tetrahedra function ===========================
function tetrahedra_volume_x6(V1,V2,V3,V4) result(Volx6) 
implicit none
!V1-V3 define triangle base plane, V4 is off plane point

!Variables - Import
real(dp) :: Volx6
real(dp) :: V1(3),V2(3),V3(3),V4(3)

!Variables - Local 
real(dp) :: m11,m12,m13,m21,m22,m23,m31,m32,m33

!Calculate signed volume
m11 = V1(1) - V4(1)
m12 = V1(2) - V4(2)
m13 = V1(3) - V4(3)
m21 = V2(1) - V4(1)
m22 = V2(2) - V4(2)
m23 = V2(3) - V4(3)
m31 = V3(1) - V4(1)
m32 = V3(2) - V4(2)
m33 = V3(3) - V4(3)
Volx6 = (m11*(m22*m33 - m23*m32) - m12*(m21*m33 - m23*m31) + m13*(m21*m32 - m22*m31))
return
end function tetrahedra_volume_x6




!Newell normal function ===========================
function newell_normal(Nvtxf,face_vtx,vertices) result(N)
implicit none 

!Variables - Import
integer(in) :: Nvtxf
real(dp) :: N(3)
integer(in), dimension(:) :: face_vtx
real(dp), dimension(:,:) :: vertices
real(dp) :: vtxC(3),vtxN(3)

!Variables - Local 
integer(in) :: vv,vc,vn 

!Initialise 
N(:) = 0.0d0 

!Accumulate normal vector to the face
do vv=1,Nvtxf
    vc = vv 
    vn = mod(vv,Nvtxf) + 1
    vc = face_vtx(vc)
    vn = face_vtx(vn)
    vtxC(:) = vertices(vc,:)
    vtxN(:) = vertices(vn,:)
    N(1) = N(1) - 0.5d0*(vtxN(3) + vtxC(3))*(vtxN(2) - vtxC(2))
    N(2) = N(2) - 0.5d0*(vtxN(1) + vtxC(1))*(vtxN(3) - vtxC(3))
    N(3) = N(3) - 0.5d0*(vtxN(2) + vtxC(2))*(vtxN(1) - vtxC(1))
end do 
return 
end function newell_normal




!Arbitrary face volume function ===========================
function arbitrary_face_vol(Aface,Xpf,Nfn) result(Vface)
implicit none 

!Variables - Import
real(dp) :: Vface,Aface
real(dp) :: Xpf(3),Nfn(3)

!Evaluate signed volume of face
Vface = (1.0d0/3.0d0)*Aface*dot_product(Xpf,Nfn)
return 
end function arbitrary_face_vol




!Angle between vectors function ===========================
function ang_vec2vec(vec1,vec2) result(ang)
implicit none 
    
!Variables - Import
real(dp) :: ang,vec1(3),vec2(3)

!Variables - Local 
real(dp) :: dotval

!Find angle 
dotval = dot_product(vec1,vec2)
dotval = dotval/(norm2(vec1)*norm2(vec2))
if (dotval .GT. 1.0d0) then 
    dotval = 1.0d0 
elseif (dotval .LT. -1.0d0) then 
    dotval = -1.0d0 
end if 
ang = abs(acos(dotval))
return 
end function ang_vec2vec




!Wendland C2 function ===========================
function wendlandc2(d,Rs) result(W)
implicit none 

!Variables - Import
real(dp) :: d,Rs,W

!Set function value 
d = d/Rs 
if (d .GE. 1.0d0) then 
    W = 0.0d0 
else
    W = ((1.0d0 - d)**4)*(4.0d0*d + 1.0d0)
end if 
return 
end function wendlandc2


end module cellmesh3d_geometry_mod