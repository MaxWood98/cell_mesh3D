!2D and 3D Geometry Routine Module
!Max Wood - mw16116@bristol.ac.uk
!Univeristy of Bristol - Department of Aerospace Engineering

!Version 8.1
!Updated 31-07-2023

!Geometry subroutines module
module cellmesh3d_geometry_mod
! use cellmesh3d_data_mod
use cellmesh3d_adtree_mod
use ieee_arithmetic, only: ieee_value,IEEE_QUIET_NAN
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




!Edge - surface transition check function (dot product) ===========================
function check_edge_surf_transitiondp(vint,ve1,ve2,surface_mesh,surface_adtree,nselected,node_select,cm3dopt) result(trans_state)
implicit none 
   
!System data
type(tree_data) :: surface_adtree
type(surface_data) :: surface_mesh
type(cm3d_options) :: cm3dopt

!Variables - Import 
integer(in) :: trans_state,nselected
integer(in) :: node_select(surface_adtree%nnode)
real(dp) :: ve1(3),ve2(3),vint(3)

!Variables - Local 
integer(in) :: nn,kk,seg_f,intfound,ftgt,tri_idx
real(dp) :: distC,dist_ref
real(dp) :: vf1(3),vf2(3),vf3(3),tri_normal(3),segdir(3),vid(4),vid_f(4)

!Set reference length 
dist_ref = 2.0d0*cm3dopt%far_field_bound

!Segment direction 
segdir(:) = ve2(:) - ve1(:)

!Find minimum distance to the surface and the corresponding surface triangle
distC = dist_ref
intfound = 0 
do nn=1,nselected
    do kk=1,surface_adtree%tree(node_select(nn))%nentry
        ftgt = surface_adtree%tree(node_select(nn))%entry(kk)
        vf1(:) = surface_mesh%vertices(surface_mesh%connectivity(ftgt,1),:)
        vf2(:) = surface_mesh%vertices(surface_mesh%connectivity(ftgt,2),:)
        vf3(:) = surface_mesh%vertices(surface_mesh%connectivity(ftgt,3),:)
        vid = min_dist_point_to_tri(vint,vf1,vf2,vf3)
        if ((vid(4) .LT. distC) .OR. (intfound == 0)) then 
            distC = vid(4)
            vid_f(:) = vid(:)
            seg_f = ftgt 
            intfound = 1
        end if 
    end do 
end do 
tri_idx = seg_f

!Check against this triangle 
vf1(:) = surface_mesh%vertices(surface_mesh%connectivity(tri_idx,1),:)
vf2(:) = surface_mesh%vertices(surface_mesh%connectivity(tri_idx,2),:)
vf3(:) = surface_mesh%vertices(surface_mesh%connectivity(tri_idx,3),:)
tri_normal = crossp3(vf2-vf1,vf3-vf1)
if (dot_product(tri_normal,segdir) .GT. 0.0d0) then !Exit 
    trans_state = -1
elseif (dot_product(tri_normal,segdir) .LT. 0.0d0) then !Entry 
    trans_state = 1
else !No transition
    trans_state = 0 
end if
return 
end function check_edge_surf_transitiondp




!Edge - surface transition check function (perturbation) ===========================
function check_edge_surf_transition(vs1,vs2,surface_mesh,surface_adtree,nselected,node_select,cm3dopt) result(trans_state)
implicit none 

!System data
type(tree_data) :: surface_adtree
type(surface_data) :: surface_mesh
type(cm3d_options) :: cm3dopt

!Variables - Import 
integer(in) :: trans_state,nselected
integer(in) :: node_select(surface_adtree%nnode)
real(dp) :: vs1(3),vs2(3)

!Variables - Local 
integer(in) :: tri_s1,tri_s2,side_vs1,side_vs2
real(dp) :: dist_ref

!Set reference length 
dist_ref = 2.0d0*cm3dopt%far_field_bound

!Vs1 find closest surface segment and the side of this segment 
call get_vtx_closest_tri_side(tri_s1,side_vs1,vs1,dist_ref,surface_adtree,nselected,node_select,surface_mesh) 

!Vs2 find closest surface segment and the side of this segment 
call get_vtx_closest_tri_side(tri_s2,side_vs2,vs2,dist_ref,surface_adtree,nselected,node_select,surface_mesh) 

!Set transition status 
if (side_vs1 .NE. side_vs2) then !transition
    if (side_vs1 .GT. 0) then !transition in 
        trans_state = 1
    else !transition out 
        trans_state = -1
    end if 
else !No transition -> check with dot product of edge direction 
    trans_state = 0 
end if
return 
end function check_edge_surf_transition




!Check vertex closest segment side subroutine ===========================
subroutine get_vtx_closest_tri_side(tri_idx,triside,vtx_tgt,dist_ref,surface_adtree,nselected,node_select,surface_mesh) 
implicit none 

!System data
type(tree_data) :: surface_adtree
type(surface_data) :: surface_mesh

!Variables - Import 
integer(in) :: tri_idx,triside,nselected
integer(in) :: node_select(surface_adtree%nnode)
real(dp) :: dist_ref
real(dp) :: vtx_tgt(3)

!Variables - Local 
integer(in) :: nn,kk,seg_f,intfound,ftgt 
real(dp) :: distC
real(dp) :: vf1(3),vf2(3),vf3(3),vid(4),vid_f(4)

!Find minimum distance to the surface and the corresponding surface triangle
distC = dist_ref
intfound = 0 
do nn=1,nselected
    do kk=1,surface_adtree%tree(node_select(nn))%nentry
        ftgt = surface_adtree%tree(node_select(nn))%entry(kk)
        vf1(:) = surface_mesh%vertices(surface_mesh%connectivity(ftgt,1),:)
        vf2(:) = surface_mesh%vertices(surface_mesh%connectivity(ftgt,2),:)
        vf3(:) = surface_mesh%vertices(surface_mesh%connectivity(ftgt,3),:)
        vid = min_dist_point_to_tri(vtx_tgt,vf1,vf2,vf3)
        if ((vid(4) .LT. distC) .OR. (intfound == 0)) then 
            distC = vid(4)
            vid_f(:) = vid(:)
            seg_f = ftgt 
            intfound = 1
        end if 
    end do 
end do 
tri_idx = seg_f

!Get triangle side for this vertex (1 = in direction of normal | -1 = opposite direction of normal | 0 = exactly on triangle)
vf1(:) = surface_mesh%vertices(surface_mesh%connectivity(tri_idx,1),:)
vf2(:) = surface_mesh%vertices(surface_mesh%connectivity(tri_idx,2),:)
vf3(:) = surface_mesh%vertices(surface_mesh%connectivity(tri_idx,3),:)
triside = check_tri_side(vtx_tgt,vf1,vf2,vf3) 
return 
end subroutine get_vtx_closest_tri_side




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
        print *, 'zface',ff
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




!==========================================================================
!2D Routines ==============================================================
!==========================================================================




!Minimum distance from a point to a line segment function ===========================
function mindist_point2line(vp,v1,v2) result(vid)
implicit none 

!Variables - Import
real(dp) :: vid(3) !xi | yi | dist 
real(dp) :: vp(2),v1(2),v2(2)

!Variables - Local 
real(dp) :: dx,dy,nx,ny
real(dp) :: vn1(2),vn2(2)

!Find line segment direction 
dx = v2(1) - v1(1)
dy = v2(2) - v1(2)

!Check for zero length segment 
if (dx*dx + dy*dy == 0.0d0) then 

    !Return intersection point is vi = v1 = v2
    vid(1:2) = v1(:)
else

    !Segment normal 
    nx = -dy
    ny = dx

    !Build normal direction vertices 
    vn1(1) = vp(1) + nx 
    vn1(2) = vp(2) + ny 
    vn2(1) = vp(1) - nx 
    vn2(2) = vp(2) - ny

    !Find intersection on base line segment 
    vid(1:2) = line_line_intersection_loc_inl1(v1,v2,vn1,vn2) 
end if

!Find distance from this intersect to vp 
vid(3) = sqrt((vp(1) - vid(1))**2 + (vp(2) - vid(2))**2)
return 
end function mindist_point2line




!Minimum distance from a point to a line segment function (unbounded) ===========================
function mindist_point2line_ub(vp,v1,v2) result(vid)
implicit none 

!Variables - Import
real(dp) :: vid(3) !xi | yi | dist 
real(dp) :: vp(2),v1(2),v2(2)

!Variables - Local 
real(dp) :: dx,dy,nx,ny
real(dp) :: vn1(2),vn2(2)

!Find line segment direction 
dx = v2(1) - v1(1)
dy = v2(2) - v1(2)

!Check for zero length segment 
if (dx*dx + dy*dy == 0.0d0) then 

    !Return intersection point is vi = v1 = v2
    vid(1:2) = v1(:)
    vid(3) = 0.0d0 
else

    !Segment normal 
    nx = -dy
    ny = dx

    !Build normal direction vertices 
    vn1(1) = vp(1) + nx 
    vn1(2) = vp(2) + ny 
    vn2(1) = vp(1) - nx 
    vn2(2) = vp(2) - ny

    !Find intersection on base line segment 
    vid(1:2) = line_line_intersection_loc_ub(v1,v2,vn1,vn2) 
end if

!Find distance from this intersect to vp 
vid(3) = sqrt((vp(1) - vid(1))**2 + (vp(2) - vid(2))**2)
return 
end function mindist_point2line_ub




!Line-Line intersection location calculation (bounded) -> L1 = v1 v2 | L2 = v3 v4 ===========================
function line_line_intersection_loc_inl1(v1,v2,v3,v4) result(vi)
implicit none 

!Variables - Import
real(dp) :: v1(2),v2(2),v3(2),v4(2)

!Variables - Local
real(dp) :: Dval,t

!Result
real(dp) :: vi(2)

!Intersection denominator 
Dval = (v1(1) - v2(1))*(v3(2) - v4(2)) - (v1(2) - v2(2))*(v3(1) - v4(1))

!Test paralell
if (Dval .NE. 0.0d0) then !If not parallel -> find intersection location
    t = ((v1(1) - v3(1))*(v3(2) - v4(2)) - (v1(2) - v3(2))*(v3(1) - v4(1)))/Dval !L1 parameter
    if (t .GT. 1.0d0) then 
        t = 1.0d0 
    elseif (t .LT. 0.0d0) then 
        t = 0.0d0 
    end if
    vi(1) = (v1(1) + t*(v2(1) - v1(1)))
    vi(2) = (v1(2) + t*(v2(2) - v1(2)))
else !Take as not intersecting as parallel  
    vi(:) = ieee_value(1.0d0,IEEE_QUIET_NAN)
end if
return 
end function line_line_intersection_loc_inl1




!Line-Line intersection location calculation (unbounded) -> L1 = v1 v2 | L2 = v3 v4 ===========================
function line_line_intersection_loc_ub(v1,v2,v3,v4) result(vi)
implicit none 

!Variables - Import
real(dp) :: v1(2),v2(2),v3(2),v4(2)

!Variables - Local
real(dp) :: Dval,t

!Result
real(dp) :: vi(2)

!Intersection denominator 
Dval = (v1(1) - v2(1))*(v3(2) - v4(2)) - (v1(2) - v2(2))*(v3(1) - v4(1))

!Test paralell
if (Dval .NE. 0.0d0) then !If not parallel -> find intersection location
    t = ((v1(1) - v3(1))*(v3(2) - v4(2)) - (v1(2) - v3(2))*(v3(1) - v4(1)))/Dval !L1 parameter
    vi(1) = (v1(1) + t*(v2(1) - v1(1)))
    vi(2) = (v1(2) + t*(v2(2) - v1(2)))
else !Take as not intersecting as parallel  
    vi(:) = ieee_value(1.0d0,IEEE_QUIET_NAN)
end if
return 
end function line_line_intersection_loc_ub




!Paralell line-line intersection -> L1 = v1 v2 | L2 = v3 v4 ===========================
subroutine line_line_intersection_parallel(nintersect,vi,v1,v2,v3,v4) 
implicit none 

!Variables - Import
integer(in) :: nintersect
real(dp) :: v1(2),v2(2),v3(2),v4(2),vi(2,2)

!Variables - Local
integer(in) :: v1_in_l2,v2_in_l2,v3_in_l1,v4_in_l1
real(dp) :: L1,L2,v1frac_l2,v2frac_l2,v3frac_l1,v4frac_l1

!Reset 
nintersect = 0 
vi(:,:) = 0.0d0 

!Segment lengths
L1 = norm2(v2(:) - v1(:))
L2 = norm2(v4(:) - v3(:))

!L1 fractions on L2
v1frac_l2 = dot_product(v1(:) - v3(:),v4(:) - v3(:))/(L2*L2) 
v2frac_l2 = dot_product(v2(:) - v3(:),v4(:) - v3(:))/(L2*L2)

!L2 fractions on L1
v3frac_l1 = dot_product(v3(:) - v1(:),v2(:) - v1(:))/(L1*L1) 
v4frac_l1 = dot_product(v4(:) - v1(:),v2(:) - v1(:))/(L1*L1) 

!Vertex end containments 
if ((v1frac_l2 .GE. 0.0d0) .AND. (v1frac_l2 .LE. 1.0d0)) then 
    v1_in_l2 = 1
else
    v1_in_l2 = 0 
end if 
if ((v2frac_l2 .GE. 0.0d0) .AND. (v2frac_l2 .LE. 1.0d0)) then 
    v2_in_l2 = 1
else
    v2_in_l2 = 0 
end if 
if ((v3frac_l1 .GE. 0.0d0) .AND. (v3frac_l1 .LE. 1.0d0)) then 
    v3_in_l1 = 1
else
    v3_in_l1 = 0 
end if 
if ((v4frac_l1 .GE. 0.0d0) .AND. (v4frac_l1 .LE. 1.0d0)) then 
    v4_in_l1 = 1
else
    v4_in_l1 = 0 
end if 

!If there is no overlap return no intersections 
if (max(v1_in_l2,v2_in_l2,v3_in_l1,v4_in_l1) == 0) then 
    nintersect = 0 
    return 
end if 

!If edge L1 is fully contained within L2 then return error case to remove l1
if ((v1_in_l2 == 1) .AND. (v2_in_l2 == 1)) then 
    nintersect = -1
    ! print *, 'l1 in l2'
    return 
end if

!If edge l2 is fully contained within l1 then return each end of l2 as intersections 
if ((v3_in_l1 == 1) .AND. (v4_in_l1 == 1)) then 
    nintersect = 2 
    vi(1,:) = v3(:)
    vi(2,:) = v4(:)
end if 

!If one vertex of l2 is contained within l1 then return this as a single intersection
if (v3_in_l1 .NE. v4_in_l1) then 
    nintersect = 1 
    if (v3_in_l1 == 1) then 
        vi(1,:) = v3(:)
    elseif (v4_in_l1 == 1) then 
        vi(1,:) = v4(:)
    end if 
end if 
return 
end subroutine line_line_intersection_parallel




!Line parallel factor calculation -> L1 = v1 v2 | L2 = v3 v4 ===========================
function get_line_line_parallel_factor(v1,v2,v3,v4) result(pfact)
implicit none 

!Variables - Import
real(dp) :: v1(2),v2(2),v3(2),v4(2)
real(dp) :: pfact

!Variables - Local 
real(dp) :: L1,L2

!Segment lengths
L1 = norm2(v2(:) - v1(:))
L2 = norm2(v4(:) - v3(:))

!Find parallel factor 
pfact = dot_product((v2(:) - v1(:)),(v3(:) - v4(:)))
if ((L1 .GT. 0.0d0) .AND. (L2 .GT. 0.0d0)) then !valid edge
    pfact = abs(pfact)/(L1*L2)
else !one edge is of zero length hence parallelism is meaningless
    pfact = 0.0d0 
end if 
return 
end function get_line_line_parallel_factor




!Line line normal distance function -> L1 = v1 v2 | L2 = v3 v4 ===========================
function line_line_normdist(v1,v2,v3,v4) result(ndist)
implicit none 

!Variables - Import
real(dp) :: v1(2),v2(2),v3(2),v4(2)
real(dp) :: ndist

!Variables - Local 
real(dp) :: dist_v1_l2,dist_v2_l2
real(dp) :: vi(3)

!v1 to L2 
vi = mindist_point2line(v1,v3,v4)
dist_v1_l2 = vi(3)

!v2 to L2 
vi = mindist_point2line(v2,v3,v4)
dist_v2_l2 = vi(3)

!Select final distance
ndist = min(dist_v1_l2,dist_v2_l2)
return 
end function line_line_normdist




!Segment - Rectangle Intersection Checking Function ===========================
function seg_rectangle_intersect_bool(vp1,vp2,vp3,vp4,A,B) result(int_type)
implicit none 

!Variables - Import
integer(in) :: int_type
real(dp) :: vp1(2),vp2(2),vp3(2),vp4(2),A(2),B(2)

!Variables - Local 
integer(in) :: int_e1,int_e2,int_e3,int_e4,in_A,in_B,int_init
real(dp) :: xmin,xmax,ymin,ymax

!Test intersection of line segment A->B against each edge of the rectangle
int_init = seg_seg_intersect_bool(A,B,vp1,vp2)
if (int_init .NE. 0) then 
    int_e1 = 1
else 
    int_e1 = 0
end if
int_init = seg_seg_intersect_bool(A,B,vp2,vp3)
if (int_init .NE. 0) then 
    int_e2 = 1
else 
    int_e2 = 0
end if
int_init = seg_seg_intersect_bool(A,B,vp3,vp4)
if (int_init .NE. 0) then 
    int_e3 = 1
else 
    int_e3 = 0
end if
int_init = seg_seg_intersect_bool(A,B,vp4,vp1)
if (int_init .NE. 0) then 
    int_e4 = 1
else 
    int_e4 = 0
end if

!Test containment of each segment end 
in_A = 0
in_B = 0 
xmin = min(vp1(1),vp2(1),vp3(1),vp4(1))
xmax = max(vp1(1),vp2(1),vp3(1),vp4(1))
ymin = min(vp1(2),vp2(2),vp3(2),vp4(2))
ymax = max(vp1(2),vp2(2),vp3(2),vp4(2))
if ((A(1) .GE. xmin) .AND. (A(1) .LE. xmax)) then 
    if ((A(2) .GE. ymin) .AND. (A(2) .LE. ymax)) then 
        in_A = 1
    end if    
end if
if ((B(1) .GE. xmin) .AND. (B(1) .LE. xmax)) then 
    if ((B(2) .GE. ymin) .AND. (B(2) .LE. ymax)) then 
        in_B = 1
    end if    
end if

!Set containment status 
if ((int_e1 == 1) .OR. (int_e2 == 1) .OR. (int_e3 == 1) .OR. (int_e4 == 1) .OR. (in_A == 1) .OR. (in_B == 1)) then 
    int_type = 1
else
    int_type = 0 
end if
return 
end function seg_rectangle_intersect_bool




!Segment - Segment Intersection Checking Function ===========================
function seg_seg_intersect_bool(A,B,C,D) result(int_type) !Segments A->B || C->D
implicit none 

!Variables - Import
integer(in) :: int_type
real(dp) :: A(2),B(2),C(2),D(2)

!Variables - Local 
integer(in) :: intA,intB,intC,intD
real(dp) :: cp_abc,cp_abd,cp_cda,cp_cdb
real(dp) :: L1,L2,Afcd,Bfcd,Cfab,Dfab

!Find segment cross products
cp_abc = cp2d(A,B,C)
cp_abd = cp2d(A,B,D)
cp_cda = cp2d(C,D,A)
cp_cdb = cp2d(C,D,B)

!Classify intersection (0 = none | 1 = proper | 2 = vertex-line touch | 3 = verex-vertex touch | 4 = colinear with overlap)
int_type = 0 
if ((cp_abc*cp_abd .LT. 0.0d0) .AND. (cp_cda*cp_cdb .LT. 0.0d0)) then !Proper intersection (both cross within each others length)
    int_type = 1
elseif ((cp_abc == 0.0d0) .AND. (cp_abd == 0.0d0) .AND. (cp_cda == 0.0d0) .AND. (cp_cdb == 0.0d0)) then !Co-linear
    L1 = norm2(A(:)-B(:))
    L2 = norm2(C(:)-D(:))
    if ((L1 == 0.0d0) .OR. (L2 == 0.0d0)) then !one segment is zero length so ignore this check 
        int_type = 0
    else
        Afcd = dot_product(A(:)-D(:),C(:)-D(:))/(L2**2)
        Bfcd = dot_product(B(:)-D(:),C(:)-D(:))/(L2**2)
        Cfab = dot_product(C(:)-B(:),A(:)-B(:))/(L1**2)
        Dfab = dot_product(D(:)-B(:),A(:)-B(:))/(L1**2)
        if ((Afcd .GT. 0.0d0) .AND. (Afcd .LT. 1.0d0)) then !Contained within segment 
            intA = 4
        elseif ((Afcd == 0.0d0) .OR. (Afcd == 1.0d0)) then !On segment vertex
            intA = 3
        else !Outside segment 
            intA = 0
        end if 
        if ((Bfcd .GT. 0.0d0) .AND. (Bfcd .LT. 1.0d0)) then !Contained within segment 
            intB = 4
        elseif ((Bfcd == 0.0d0) .OR. (Bfcd == 1.0d0)) then !On segment vertex
            intB = 3
        else !Outside segment 
            intB = 0
        end if 
        if ((Cfab .GT. 0.0d0) .AND. (Cfab .LT. 1.0d0)) then !Contained within segment 
            intC = 4
        elseif ((Cfab == 0.0d0) .OR. (Cfab == 1.0d0)) then !On segment vertex
            intC = 3
        else !Outside segment 
            intC = 0
        end if 
        if ((Dfab .GT. 0.0d0) .AND. (Dfab .LT. 1.0d0)) then !Contained within segment 
            intD = 4
        elseif ((Dfab == 0.0d0) .OR. (Dfab == 1.0d0)) then !On segment vertex
            intD = 3
        else !Outside segment 
            intD = 0
        end if 
        if ((intA == 4) .OR. (intB == 4) .OR. (intC == 4) .OR. (intD == 4)) then !Segment - Segment overlap 
            int_type = 4
        elseif ((intA == 3) .OR. (intB == 3) .OR. (intC == 3) .OR. (intD == 3)) then !Vertex - Vertex touch 
            int_type = 3
        else !No intersection
            int_type = 0 
        end if 
    end if 
elseif (cp_abc == 0.0d0) then !C lies on AB
    if (cp_cda*cp_cdb .LT. 0.0d0) then !Within AB
        int_type = 2
    elseif (cp_cda*cp_cdb .GT. 0.0d0) then !Outside AB
        int_type = 0 
    elseif (cp_cda*cp_cdb == 0.0d0) then !On vertex A or B
        int_type = 3
    end if 
elseif (cp_abd == 0.0d0) then !D lies on AB
    if (cp_cda*cp_cdb .LT. 0.0d0) then !Within AB
        int_type = 2
    elseif (cp_cda*cp_cdb .GT. 0.0d0) then !Outside AB
        int_type = 0 
    elseif (cp_cda*cp_cdb == 0.0d0) then !On vertex A or B
        int_type = 3
    end if
elseif (cp_cda == 0.0d0) then !A lies on CD
    if (cp_abc*cp_abd .LT. 0.0d0) then !Within CD
        int_type = 2
    elseif (cp_abc*cp_abd .GT. 0.0d0) then !Outside CD  
        int_type = 0 
    elseif (cp_abc*cp_abd == 0.0d0) then !On vertex C or D
        int_type = 3
    end if
elseif (cp_cdb == 0.0d0) then !B lies on CD
    if (cp_abc*cp_abd .LT. 0.0d0) then !Within CD
        int_type = 2
    elseif (cp_abc*cp_abd .GT. 0.0d0) then !Outside CD  
        int_type = 0 
    elseif (cp_abc*cp_abd == 0.0d0) then !On vertex C or D
        int_type = 3
    end if
else !No intersection within either lines length
    int_type = 0 
end if  
return 
end function seg_seg_intersect_bool




!Side of segment checking function ===========================
function check_segment_side(vl1,vl2,vp) result(side)
implicit none 

!Variables - Import
integer(in) :: side
real(dp) :: vl1(2),vl2(2),vp(2)

!Variables - Local 
real(dp) :: cp_vl1_vl2_vp

!Find segment cross product
cp_vl1_vl2_vp = cp2d(vl1,vl2,vp)

!Clasify side of segment
if (cp_vl1_vl2_vp .GT. 0.0d0) then !'Right' normal direction side
    side = 1
elseif (cp_vl1_vl2_vp .LT. 0.0d0) then !'Left' normal opposide direction side
    side = -1 
else !Lies on the segment vl1->vl2
    side = 0 
end if
return 
end function check_segment_side




!Anticlockwise andgle function ===========================
function acw_ang(A,B,C) result(ang) !vectors (A->B) and (A->C)
implicit none 

!Variables - Import
real(dp) :: ang
real(dp) :: A(2),B(2),C(2)

!Variables - Local 
real(dp) :: dotp,crossp 
real(dp) :: vec1N(2),vec2N(2)

!Normalised direction vectors 
vec1N(1) = B(1) - A(1)
vec1N(2) = B(2) - A(2)
vec1N = vec1N/norm2(vec1N)
vec2N(1) = C(1) - A(1)
vec2N(2) = C(2) - A(2)
vec2N = vec2N/norm2(vec2N)

!Cross and dot products
dotp = dot_product(Vec1N,Vec2N)
crossp = cp2dv(vec1N,vec2N)

!Return angle 
ang = atan2(crossp,dotp)
if (ang .LT. 0.0d0) then 
    ang = ang + 8.0d0*atan(1.0d0)
end if 
return 
end function acw_ang




!2D Cross product functions ===========================
function cp2d(A,B,C) result(val) !vectors (A->B) and (A->C)
implicit none 

!Variables - Import
real(dp) :: val
real(dp) :: A(2),B(2),C(2)

!Evaluate
val = (B(1) - A(1))*(C(2) - A(2)) - (B(2) - A(2))*(C(1) - A(1))
return 
end function cp2d

function cp2dv(vec1N,vec2N) result(val) !vectors (vec1N = (A->B)) and (vec2N = (A->C))
implicit none 

!Variables - Import
real(dp) :: val
real(dp) :: vec1N(2),vec2N(2)

!Evaluate
val = vec1N(1)*vec2N(2) - vec1N(2)*vec2N(1)
return 
end function cp2dv
    



!Line segment area contribution function ===========================
function Asegment(v1,v2) result(Aseg) !+ve area for CCW oriented shapes
implicit none 

!Variables - Import
real(dp) :: v1(2),v2(2)

!Result
real(dp) :: Aseg

!Segment area
Aseg = 0.5d0*(v1(1)*v2(2) - v2(1)*v1(2))
return 
end function Asegment




!==========================================================================
!3D Routines ==============================================================
!==========================================================================




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




!Line-triangle intersection bool Function ===========================
function line_tri_intersect_bool(vl1,vl2,vt1,vt2,vt3) result(int_bool)
implicit none 

!Variables - Import
integer(in) :: int_bool
real(dp) :: vl1(3),vl2(3),vt1(3),vt2(3),vt3(3)

!Variables - Local 
integer(in) :: within_tri_and_line
real(dp) :: Xi,Yi,Zi,u,v,ti,baryval
real(dp) :: O(3),D(3),Nt(3)

!Set ray parameters 
O(:) = vl1(:)
D(:) = vl2(:) - vl1(:)

!Triangle normal 
Nt = crossp3(vt2-vt1,vt3-vt1)

!Check intersection
if (dot_product(D,Nt) .NE. 0.0d0) then !Intersection 

    !Find intersection
    call ray_plane_intersection(Xi,Yi,Zi,u,v,ti,O,D,vt1,vt2,vt3)
    baryval = u + v 
    if (isnan(baryval)) then !error returned
        int_bool = 0 
        return 
    end if 

    !Check if within the triangle and line
    within_tri_and_line = 0 
    if ((ti .GE. 0.0d0) .AND. (ti .LE. 1.0d0)) then 
        if ((u .GE. 0.0d0) .AND. (u .LE. 1.0d0)) then 
            if ((v .GE. 0.0d0) .AND. (v .LE. 1.0d0)) then 
                if ((baryval .GE. 0.0d0) .AND.(baryval .LE. 1.0d0)) then 
                    within_tri_and_line = 1
                end if 
            end if 
        end if 
    end if 
    if (within_tri_and_line == 1) then 
        int_bool = 1
    else
        int_bool = 0 
    end if 
else !Ray is parallel to plane so no intersection
    int_bool = 0 
end if 
return 
end function line_tri_intersect_bool




!Line-triangle intersection Function ===========================
function line_tri_intersect(vl1,vl2,vt1,vt2,vt3) result(vint)
implicit none 

!Variables - Import
real(dp) :: vint(3)
real(dp) :: vl1(3),vl2(3),vt1(3),vt2(3),vt3(3)

!Variables - Local 
real(dp) :: Xi,Yi,Zi,u,v,ti
real(dp) :: O(3),D(3),Nt(3)

!Set ray parameters 
O(:) = vl1(:)
D(:) = vl2(:) - vl1(:)

!Triangle normal 
Nt = crossp3(vt2-vt1,vt3-vt1)

!Check intersection
if (dot_product(D,Nt) .NE. 0.0d0) then !Intersection 

    !Find intersection
    call ray_plane_intersection(Xi,Yi,Zi,u,v,ti,O,D,vt1,vt2,vt3)

    !Store
    vint(1) = Xi
    vint(2) = Yi
    vint(3) = Zi

    !Bound within line 
    if (ti .GT. 1.0d0) then 
        vint(:) = vl2(:)
    elseif (ti .LT. 0.0d0) then 
        vint(:) = vl1(:)
    end if 
else !Ray is parallel to plane so no intersection
    vint(:) = ieee_value(1.0d0,IEEE_QUIET_NAN)
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
real(dp) :: Xi,Yi,Zi,u,v,ti
real(dp) :: O(3),D(3),Nt(3)

!Set ray parameters 
O(:) = vl1(:)
D(:) = vl2(:) - vl1(:)

!Triangle normal 
Nt = crossp3(vt2-vt1,vt3-vt1)

!Check intersection
if (dot_product(D,Nt) .NE. 0.0d0) then !Intersection 

    !Find intersection
    call ray_plane_intersection(Xi,Yi,Zi,u,v,ti,O,D,vt1,vt2,vt3)

    !Store
    bcint(1) = u
    bcint(2) = v
    bcint(3) = 1.0d0 - u - v
else !Ray is parallel to plane so no intersection
    bcint(:) = ieee_value(1.0d0,IEEE_QUIET_NAN)
end if 
return 
end function baryc_line_tri_intersect




!Function to classify a vertices location WRT a triangle using barycentric coordinates ===========================
function vtx_bary_tri_loc(u,v,w,ztol) result(vtx_loc)
implicit none 

!Variables - Import
character(len=2) :: vtx_loc
real(dp) :: u,v,w,ztol

!Determine location
vtx_loc = 'na'
if ((u <= -ztol) .OR. (v <= -ztol) .OR. (w <= -ztol)) then !outside triangle (in tollerance)
    vtx_loc = 'ot'
elseif ((u >= ztol) .AND. (v >= ztol) .AND. (w >= ztol)) then !inside triangle (in tollerance)
    vtx_loc = 'in'
elseif ((u >= ztol) .AND. (v <= ztol) .AND. (w >= ztol)) then !edge (v1-v2) (in tollerance)
    vtx_loc = 'e1'
elseif ((u >= ztol) .AND. (v >= ztol) .AND. (w <= ztol)) then !edge (v2-v3) (in tollerance)
    vtx_loc = 'e2'
elseif ((u <= ztol) .AND. (v >= ztol) .AND. (w >= ztol)) then !edge (v3-v1) (in tollerance)
    vtx_loc = 'e3'
elseif ((u <= ztol) .AND. (v <= ztol) .AND. (w >= ztol)) then !vertex 1 (in tollerance)
    vtx_loc = 'v1'
elseif ((u >= ztol) .AND. (v <= ztol) .AND. (w <= ztol)) then !vertex 2 (in tollerance)
    vtx_loc = 'v2'
elseif ((u <= ztol) .AND. (v >= ztol) .AND. (w <= ztol)) then !vertex 3 (in tollerance)
    vtx_loc = 'v3'
else
    vtx_loc = 'uc'
end if 
return 
end function vtx_bary_tri_loc




!Cartesian to barycentric coordinate projection on plane ===========================
function cart2baryc_onplane(vp,vt1,vt2,vt3) result(vpbc)
implicit none 

!Variables - Import
real(dp) :: vp(3),vt1(3),vt2(3),vt3(3),vpbc(3)

!Variables - Local 
real(dp) :: Xi,Yi,Zi,u,v,w,ti
real(dp) :: O(3),D(3),Nt(3)

!Triangle normal 
Nt = crossp3(vt2-vt1,vt3-vt1)

!Set ray parameters 
O(:) = vp(:)
D(:) = Nt(:)

!Find location on the triangle plane in barycentric coordinates
call ray_plane_intersection(Xi,Yi,Zi,u,v,ti,O,D,vt1,vt2,vt3)
w = 1.0d0 - u - v

!Return barycentric coordinates
vpbc(1) = u  
vpbc(2) = v 
vpbc(3) = w 
return 
end function cart2baryc_onplane




!Minimum distance point to triangle function ===========================
function min_dist_point_to_tri(vp,vt1,vt2,vt3) result(vid)
implicit none 

!Variables - Import
real(dp) :: vp(3),vt1(3),vt2(3),vt3(3),vid(4)

!Variables - Local 
integer(in) :: within_tri,mindedge
real(dp) :: Xi,Yi,Zi,u,v,ti,baryval
real(dp) :: O(3),D(3),Nt(3),vid_e1(4),vid_e2(4),vid_e3(4),edists(3)

!Triangle normal 
Nt = crossp3(vt2-vt1,vt3-vt1)

!Set ray parameters 
O(:) = vp(:)
D(:) = Nt(:)

!Find intersection with triangle plane 
call ray_plane_intersection(Xi,Yi,Zi,u,v,ti,O,D,vt1,vt2,vt3)
baryval = u + v 
within_tri = 0 
if ((u .GE. 0.0d0) .AND. (u .LE. 1.0d0)) then 
    if ((v .GE. 0.0d0) .AND. (v .LE. 1.0d0)) then 
        if ((baryval .GE. 0.0d0) .AND.(baryval .LE. 1.0d0)) then 
            within_tri = 1
        end if 
    end if 
end if 

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




!Check side of triangle function ===========================
function check_tri_side(vp,vt1,vt2,vt3) result(tri_side)
implicit none 

!Variables - Import
integer(in) :: tri_side
real(dp) :: vp(3),vt1(3),vt2(3),vt3(3)

!Variables - Local 
real(dp) :: tetvol 

!Find tetrahedra signed volume 
tetvol = tetrahedra_volume_x6(vt1,vt2,vt3,vp)

!Clasify side of triangle
if (tetvol .GT. 0.0d0) then !'Right' normal direction side
    tri_side = -1
elseif (tetvol .LT. 0.0d0) then !'Left' normal opposide direction side
    tri_side = 1 
else !Lies on the segment vl1->vl2
    tri_side = 0 
end if
return 
end function check_tri_side




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
!f = norm2(vid(1:3) - ve1(:))/norm2(ve2(:) - ve1(:))
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
function tetrahedra_volume_x6(V1,V2,V3,V4) result(Volx6) !6*volume to accuratly extract sign
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
Volx6 = (m11*(m22*m33 - m23*m32) - m12*(m21*m33 - m23*m31) + m13*(m21*m32 - m22*m31)) !For accurate sign of volume
!Vol = Volx6*(1.0d0/6.0d0) !Actual volume
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
!vec1 / vec2 must be unit vectors
    
!Variables - Import
real(dp) :: ang,vec1(3),vec2(3)

!Variables - Local 
real(dp) :: dotval,pi

!Define pi 
pi = 3.141592653589793d0 

!Find angle 
dotval = dot_product(vec1,vec2)
if (dotval .GT. 1.0d0) then 
    dotval = 1.0d0 
elseif (dotval .LT. -1.0d0) then 
    dotval = -1.0d0 
end if 
ang = abs(acos(dotval)*(180.0d0/pi))
return 
end function ang_vec2vec


end module cellmesh3d_geometry_mod