!Cell Mesh 3D Surface Mesh Module
!Max Wood - mw16116@bristol.ac.uk
!Univeristy of Bristol - Department of Aerospace Engineering

!Version 1.0
!Updated 03-08-2023

!Module
module cellmesh3d_surface_mod
use cellmesh3d_geometry_mod
contains 


!Orient surface mesh subroutine ===========================
subroutine orient_surface(surface_mesh,cm3dopt)
implicit none 

!Variables - Import
type(cm3d_options) :: cm3dopt
type(surface_data) :: surface_mesh

!Variables - Local 
integer(in) :: ff
integer(in) :: ftemp(3)
real(dp) :: Vface,vol_total 
real(dp) :: vt1(3),vt2(3),vt3(3),v0(3)

!Evaluate volume of objects tetrahedra_volume_x6(V1,V2,V3,V4) 
v0(:) = 0.0d0 
vol_total = 0.0d0 
do ff=1,surface_mesh%nfcs
    vt1(:) = surface_mesh%vertices(surface_mesh%connectivity(ff,1),:)
    vt2(:) = surface_mesh%vertices(surface_mesh%connectivity(ff,2),:)
    vt3(:) = surface_mesh%vertices(surface_mesh%connectivity(ff,3),:)
    Vface = tetrahedra_volume_x6(vt1,vt2,vt3,v0) 
    vol_total = vol_total + Vface
end do 

!Flip if negative volume 
if (vol_total .LT. 0.0d0) then 
    do ff=1,surface_mesh%nfcs
        ftemp(:) = surface_mesh%connectivity(ff,:)
        surface_mesh%connectivity(ff,1) = ftemp(3)
        surface_mesh%connectivity(ff,2) = ftemp(2)
        surface_mesh%connectivity(ff,3) = ftemp(1)
    end do 
    if (cm3dopt%dispt == 1) then
        write(*,'(A)') '    {surface mesh re-oriented for positive volume}'
    end if
end if 
return 
end subroutine orient_surface




!Evaluate surface radius of curvature subroutine ===========================
subroutine evaluate_surf_rcurv(surface_mesh,cm3dopt)
use ieee_arithmetic
implicit none 

!Variables - Import
type(cm3d_options) :: cm3dopt
type(surface_data) :: surface_mesh

!Variables - Local
integer(in) :: ff,ee,aa,vv
integer(in) :: f1,f2,v1,v2,v3,Nedge,fv0,etgt,ftgt,ev1,ev2,vi,vj
integer(in) :: edges(2*surface_mesh%nfcs,2)
integer(in), dimension(:,:), allocatable :: V2E,V2F,F2E,E2F
real(dp) :: cot_ang1,cot_ang2,Abrc
real(dp) :: Aface(surface_mesh%nfcs),vt1(3),vt2(3),lapP(3)
real(dp), dimension(:,:), allocatable :: cot_weight

!Build surface mesh connectivity 
allocate(surface_mesh%valence(surface_mesh%nvtx))
call get_valence(surface_mesh%valence,surface_mesh%maxValence,surface_mesh%connectivity,surface_mesh%nfcs,surface_mesh%nvtx)
call construct_edges(Nedge,edges,surface_mesh%nvtx,surface_mesh%nfcs,surface_mesh%connectivity)
call get_connectivity(V2E,V2F,F2E,E2F,surface_mesh%valence,surface_mesh%maxvalence,Nedge,surface_mesh%nfcs,surface_mesh%nvtx,&
edges,surface_mesh%connectivity)

!Build cotangent weights 
allocate(cot_weight(Nedge,2))
cot_weight(:,:) = 0.0d0 
do ee=1,Nedge 

    !Faces on this edge
    f1 = E2F(ee,1)
    f2 = E2F(ee,2)

    !Edge vertices 
    ev1 = edges(ee,1)
    ev2 = edges(ee,2)

    !Evaluate for each face
    if (f1 .GT. 0) then 

        !Find other vertex in this face
        fv0 = 0
        do aa=1,3
            if ((surface_mesh%connectivity(f1,aa) .NE. ev1) .AND. (surface_mesh%connectivity(f1,aa) .NE. ev2)) then 
                fv0 = surface_mesh%connectivity(f1,aa)
                exit
            end if 
        end do 

        !Find cot weight f1 
        vt1 = surface_mesh%vertices(ev1,:) - surface_mesh%vertices(fv0,:)
        vt2 = surface_mesh%vertices(ev2,:) - surface_mesh%vertices(fv0,:)
        cot_ang1 = dot_product(vt1,vt2)/norm2(crossp3(vt1,vt2)) 
        if (isnan(cot_ang1)) then 
            cot_ang1 = 0.0d0 
        end if 
    else
        cot_ang1 = 0.0d0 
    end if 
    if (f2 .GT. 0) then 

        !Find other vertex in this face
        fv0 = 0
        do aa=1,3
            if ((surface_mesh%connectivity(f2,aa) .NE. ev1) .AND. (surface_mesh%connectivity(f2,aa) .NE. ev2)) then 
                fv0 = surface_mesh%connectivity(f2,aa)
                exit
            end if 
        end do 

        !Find cot weight f2 
        vt1 = surface_mesh%vertices(ev1,:) - surface_mesh%vertices(fv0,:)
        vt2 = surface_mesh%vertices(ev2,:) - surface_mesh%vertices(fv0,:)
        cot_ang2 = dot_product(vt1,vt2)/norm2(crossp3(vt1,vt2)) 
        if (isnan(cot_ang2)) then 
            cot_ang2 = 0.0d0 
        end if 
    else
        cot_ang2 = 0.0d0 
    end if 

    !Store
    if ((f1 .GT. 0) .AND. (f2 .GT. 0)) then 
        cot_weight(ee,1) = cot_ang1
        cot_weight(ee,2) = cot_ang2
    end if 
end do 

!Find area of each face
do ff=1,surface_mesh%nfcs

    !%Vertices
    v1 = surface_mesh%connectivity(ff,1)
    v2 = surface_mesh%connectivity(ff,2)
    v3 = surface_mesh%connectivity(ff,3)

    !Vectors
    vt1(:) = surface_mesh%vertices(v2,:) - surface_mesh%vertices(v1,:)
    vt2(:) = surface_mesh%vertices(v3,:) - surface_mesh%vertices(v1,:)

    !Area
    Aface(ff) = 0.5d0*norm2(crossp3(vt1,vt2))
end do 

!Find curvature at each surface vertex
allocate(surface_mesh%vtx_rcurv(surface_mesh%nvtx))
surface_mesh%vtx_rcurv(:) = 0.0d0 
do vv=1,surface_mesh%nvtx 

    !Accumulate barycentric area
    Abrc = 0.0d0 
    do aa=1,surface_mesh%maxValence
        ftgt = V2F(vv,aa)
        if (ftgt .GT. 0) then 
            Abrc = Abrc + Aface(ftgt)
        end if 
    end do 
    Abrc = Abrc/3.0d0

    !Accumulate curvature
    lapP(:) = 0
    do aa=1,surface_mesh%valence(vv)
        etgt = V2E(vv,aa)
        if (etgt .GT. 0) then  

            !Vertices 
            if (edges(etgt,1) == vv) then 
                vi = edges(etgt,1)
                vj = edges(etgt,2)
            else
                vi = edges(etgt,2)
                vj = edges(etgt,1)
            end if 

            !Accumulate
            lapP(:) = lapP(:) + (1.0d0/(2.0d0*Abrc))*(cot_weight(etgt,1) + cot_weight(etgt,2))*&
            (surface_mesh%vertices(vj,:) - surface_mesh%vertices(vi,:))
        end if 
    end do 

    !Store mean radius of curvature 
    surface_mesh%vtx_rcurv(vv) = (1.0d0/(0.5d0*norm2(lapP)*cm3dopt%surfRcurvM + 1e-12_dp))
end do 
! open(11,file='io/vrcurve.dat')
!     do vv=1,surface_mesh%nvtx 
!         write(11,*) surface_mesh%vtx_rcurv(vv) 
!     end do 
! close(11)

!Average curvature to faces
allocate(surface_mesh%face_rcurv(surface_mesh%nfcs))
allocate(surface_mesh%face_maxcurv(surface_mesh%nfcs))
surface_mesh%face_rcurv(:) = 0.0d0 
surface_mesh%face_maxcurv(:) = 0.0d0 
do ff=1,surface_mesh%nfcs
    surface_mesh%face_rcurv(ff) = (surface_mesh%vtx_rcurv(surface_mesh%connectivity(ff,1)) + &
    surface_mesh%vtx_rcurv(surface_mesh%connectivity(ff,2)) + &
    surface_mesh%vtx_rcurv(surface_mesh%connectivity(ff,3)))/3.0d0
    surface_mesh%face_maxcurv(ff) = max(surface_mesh%vtx_rcurv(surface_mesh%connectivity(ff,1)),&
    surface_mesh%vtx_rcurv(surface_mesh%connectivity(ff,2)),&
    surface_mesh%vtx_rcurv(surface_mesh%connectivity(ff,3)))
end do 
return 
end subroutine evaluate_surf_rcurv




!Subroutine to evaluate valence in a mesh ===========================
subroutine get_valence(valence,maxValence,faces,nface,nvtx)
implicit none 

!Variables - Import
integer(in) :: maxValence,nface,nvtx
integer(in), dimension(:) :: valence
integer(in), dimension(:,:) :: faces

!Variables - Local
integer(in) :: ff,ee,vv
integer(in) :: ev1,ev2,ubValence,evalid
integer(in), dimension(:,:), allocatable :: vconnect

!Initialse valence array
valence(:) = 0 

!Upper bound of maximum valence 
valence(:) = 0 
do ff=1,nface
    do ee=1,3

        !Edge end vertices
        ev1 = ee
        ev2 = mod(ee,3) + 1
        ev1 = faces(ff,ev1) 
        ev2 = faces(ff,ev2) 

        !Accumulate valence 
        valence(ev1) = valence(ev1) + 1
        valence(ev2) = valence(ev2) + 1
    end do 
end do 
ubValence = 2*maxval(valence)

!Construct actual valence of each vertex
allocate(vconnect(nvtx,ubValence))
vconnect(:,:) = 0 
valence(:) = 0
do ff=1,nface
    do ee=1,3

        !Edge end vertices
        ev1 = ee
        ev2 = mod(ee,3) + 1
        ev1 = faces(ff,ev1) 
        ev2 = faces(ff,ev2) 

        !Check against vconnect 
        evalid = 1
        do vv=1,ubValence
            if (vconnect(ev1,vv) == ev2) then 
                evalid = 0
                exit
            end if 
            if (vconnect(ev2,vv) == ev1) then 
                evalid = 0
                exit 
            end if 
        end do 

        !Add valence if new edge
        if (evalid == 1) then 
            
            !Increment valence on each vertex
            valence(ev1) = valence(ev1) + 1
            valence(ev2) = valence(ev2) + 1

            !Update vconnect
            do vv=1,ubValence
                if (vconnect(ev1,vv) == 0) then 
                    vconnect(ev1,vv) = ev2 
                    exit 
                end if 
            end do 
            do vv=1,ubValence
                if (vconnect(ev2,vv) == 0) then 
                    vconnect(ev2,vv) = ev1 
                    exit 
                end if 
            end do 
        end if 
    end do 
end do 

!Set maximum valence
maxValence = maxval(valence(:))
return 
end subroutine get_valence




!Subroutine to construct edges from triangular faces in a mesh ===========================
subroutine construct_edges(Nedge,edges,Nvtx,Nface,faces)
implicit none 

!Variables - Import
integer(in) :: Nedge,Nvtx,Nface
integer(in), dimension(:,:) :: edges,faces

!Variables - Local 
integer(in) :: ff,ee,vv 
integer(in) :: maxvalence,ev1,ev2,evalid
integer(in) :: valence(Nvtx)
integer(in), dimension(:,:), allocatable :: vconnect

!Set maxvalence
valence(:) = 0 
do ff=1,Nface
    do ee=1,3

        !Edge end vertices
        ev1 = ee
        ev2 = mod(ee,3) + 1
        ev1 = faces(ff,ev1) 
        ev2 = faces(ff,ev2) 

        !Accumulate valence 
        valence(ev1) = valence(ev1) + 1
        valence(ev2) = valence(ev2) + 1
    end do 
end do 
maxvalence = 2*maxval(valence)

!Allocate vconnect
allocate(vconnect(Nvtx,maxvalence))
vconnect(:,:) = 0 

!Build edges 
Nedge = 0 
do ff=1,Nface
    do ee=1,3

        !Edge end vertices
        ev1 = ee
        ev2 = mod(ee,3) + 1
        ev1 = faces(ff,ev1) 
        ev2 = faces(ff,ev2) 

        !Check against vconnect 
        evalid = 1 
        do vv=1,maxvalence
            if (vconnect(ev1,vv) == ev2) then 
                evalid = 0
                exit 
            end if 
            if (vconnect(ev2,vv) == ev1) then 
                evalid = 0
                exit 
            end if 
        end do 

        !Add if valid 
        if (evalid == 1) then 

            !Increment edge count 
            Nedge = Nedge + 1

            !Add edge 
            edges(Nedge,1) = ev1
            edges(Nedge,2) = ev2

            !Add edge to connection structure     
            do vv=1,maxvalence
                if (vconnect(ev1,vv) == 0) then 
                    vconnect(ev1,vv) = ev2 
                    exit 
                end if 
            end do 
            do vv=1,maxvalence
                if (vconnect(ev2,vv) == 0) then 
                    vconnect(ev2,vv) = ev1 
                    exit 
                end if 
            end do 
        end if 
    end do 
end do 
return 
end subroutine construct_edges




!Subroutine to construct connectivity ===========================
subroutine get_connectivity(V2E,V2F,F2E,E2F,valence,maxvalence,Nedge,Nface,Nvtx,edges,faces)
implicit none 

!Variables - Import
integer(in) :: Nedge,Nface,Nvtx
integer(in), dimension(:,:) :: edges,faces
integer(in), dimension(:) :: valence
integer(in), dimension(:,:), allocatable :: V2E,V2F,F2E,E2F

!Variables - Local 
integer(in) :: ff,ee,vv
integer(in) :: ev1,ev2,edgc,vtxc,maxvalence

!Build V2F
allocate(V2F(Nvtx,maxvalence))
V2F(:,:) = 0 
do ff=1,Nface
    do ee=1,3
        vtxc = faces(ff,ee)
        do vv=1,maxvalence
            if (V2F(vtxc,vv) == 0) then 
                V2F(vtxc,vv) = ff
                exit
            end if 
        end do 
    end do 
end do 

!Build V2E
allocate(V2E(Nvtx,maxvalence))
V2E(:,:) = 0 
do ee=1,Nedge
    
    !Edge end vertices
    ev1 = edges(ee,1) 
    ev2 = edges(ee,2) 

    !Add to V2E
    do vv=1,maxvalence
        if (V2E(ev1,vv) == 0) then 
            V2E(ev1,vv) = ee
            exit
        end if 
    end do 
    do vv=1,maxvalence
        if (V2E(ev2,vv) == 0) then 
            V2E(ev2,vv) = ee
            exit
        end if 
    end do
end do 

!Build F2E (ordered)
allocate(F2E(Nface,3))
F2E(:,:) = 0 
do ff=1,Nface
    do ee=1,3

        !Edge end vertices
        ev1 = ee
        ev2 = mod(ee,3) + 1
        ev1 = faces(ff,ev1) 
        ev2 = faces(ff,ev2) 

        !Find edge joining ev1-ev2 in this face from V2E of ev1
        do vv=1,valence(ev1)
            edgc = V2E(ev1,vv)
            if ((edges(edgc,1) == ev2) .OR. (edges(edgc,2) == ev2)) then 
                F2E(ff,ee) = edgc
                exit
            end if 
        end do 
    end do 
end do 

!Build E2F
allocate(E2F(Nedge,2))
E2F(:,:) = 0 
do ff=1,Nface
    do ee=1,3
        edgc = F2E(ff,ee)
        if ((E2F(edgc,1) .NE. ff) .AND. (E2F(edgc,2) .NE. ff)) then 
            if (E2F(edgc,1) == 0) then 
                E2F(edgc,1) = ff
            elseif (E2F(edgc,2) == 0) then 
                E2F(edgc,2) = ff
            end if 
        end if 
    end do 
end do 
return 
end subroutine get_connectivity


end module cellmesh3d_surface_mod