!Cell Mesh 3D Surface Mesh Module
!Max Wood - mw16116@bristol.ac.uk
!Univeristy of Bristol - Department of Aerospace Engineering

!Version 2.0
!Updated 30-10-2023

!Module
module cellmesh3d_surface_mod
use cellmesh3d_geometry_mod
use cellmesh3d_connectivity_mod
contains 


!Preprocess surface mesh ===========================
subroutine preprocess_surface_mesh(surface_mesh,cm3dopt)
implicit none 

!Variables - Import
type(cm3d_options) :: cm3dopt
type(surface_data) :: surface_mesh

!Orient surface mesh for positive object volume (this ensures normal vector convention is correct)
call orient_surface(surface_mesh,cm3dopt)

!Build surface mesh connectivity 
call get_tri_valence(surface_mesh%valence,surface_mesh%maxValence,surface_mesh%connectivity,surface_mesh%nfcs,surface_mesh%nvtx)
call construct_edges(surface_mesh%nedge,surface_mesh%edges,surface_mesh%valence,surface_mesh%nvtx,surface_mesh%nfcs,&
surface_mesh%connectivity)
call get_connectivity(surface_mesh%V2E,surface_mesh%V2F,surface_mesh%F2E,surface_mesh%E2F,surface_mesh%valence,&
surface_mesh%maxvalence,surface_mesh%nedge,surface_mesh%nfcs,surface_mesh%nvtx,surface_mesh%edges,surface_mesh%connectivity)
call construct_surface_normals(surface_mesh)

!Evaluate surface mesh face curvature 
call evaluate_surf_rcurv(surface_mesh,cm3dopt)
return 
end subroutine preprocess_surface_mesh




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
integer(in) :: f1,f2,v1,v2,v3,fv0,etgt,ftgt,ev1,ev2,vi,vj
real(dp) :: pi,cot_ang1,cot_ang2,Abrc,angtot,kl
real(dp) :: Aface(surface_mesh%nfcs),vt1(3),vt2(3),lapP(3)
real(dp), dimension(:,:), allocatable :: cot_weight

!Define pi
pi = 4.0d0*atan(1.0d0)

!Build cotangent weights 
allocate(cot_weight(surface_mesh%nedge,2))
cot_weight(:,:) = 0.0d0 
do ee=1,surface_mesh%nedge 

    !Faces on this edge
    f1 = surface_mesh%E2F(ee,1)
    f2 = surface_mesh%E2F(ee,2)

    !Edge vertices 
    ev1 = surface_mesh%edges(ee,1)
    ev2 = surface_mesh%edges(ee,2)

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
allocate(surface_mesh%face_area(surface_mesh%nfcs))
do ff=1,surface_mesh%nfcs

    !Vertices
    v1 = surface_mesh%connectivity(ff,1)
    v2 = surface_mesh%connectivity(ff,2)
    v3 = surface_mesh%connectivity(ff,3)

    !Vectors
    vt1(:) = surface_mesh%vertices(v2,:) - surface_mesh%vertices(v1,:)
    vt2(:) = surface_mesh%vertices(v3,:) - surface_mesh%vertices(v1,:)

    !Area
    Aface(ff) = 0.5d0*norm2(crossp3(vt1,vt2))
    surface_mesh%face_area(ff) = Aface(ff)
end do 

!Find curvature at each surface vertex
allocate(surface_mesh%vtx_meancurv(surface_mesh%nvtx))
allocate(surface_mesh%vtx_gausscurv(surface_mesh%nvtx))
surface_mesh%vtx_meancurv(:) = 0.0d0 
surface_mesh%vtx_gausscurv(:) = 0.0d0 
do vv=1,surface_mesh%nvtx 

    !Accumulate barycentric area
    Abrc = 0.0d0 
    do aa=1,surface_mesh%maxValence
        ftgt = surface_mesh%V2F(vv,aa)
        if (ftgt .GT. 0) then 
            Abrc = Abrc + Aface(ftgt)
        end if 
    end do 
    Abrc = Abrc/3.0d0

    !Accumulate mean curvature
    lapP(:) = 0
    do aa=1,surface_mesh%valence(vv)
        etgt = surface_mesh%V2E(vv,aa)
        if (etgt .GT. 0) then  

            !Vertices 
            if (surface_mesh%edges(etgt,1) == vv) then 
                vi = surface_mesh%edges(etgt,1)
                vj = surface_mesh%edges(etgt,2)
            else
                vi = surface_mesh%edges(etgt,2)
                vj = surface_mesh%edges(etgt,1)
            end if 

            !Accumulate
            lapP(:) = lapP(:) + (1.0d0/(2.0d0*Abrc))*(cot_weight(etgt,1) + cot_weight(etgt,2))*&
            (surface_mesh%vertices(vj,:) - surface_mesh%vertices(vi,:))
        end if 
    end do 
    surface_mesh%vtx_meancurv(vv) = 0.5d0*norm2(lapP)*sign(1.0d0,dot_product(surface_mesh%vtx_normal(vv,:),-lapP(:)))

    !Acculumate gaussian curvature 
    angtot = 0.0d0 
    do ff=1,surface_mesh%valence(vv)
        ftgt = surface_mesh%V2F(vv,ff)
        if (ftgt .GT. 0) then  

            !Find other vertices in this face
            vi = 0
            do aa=1,3
                if (surface_mesh%connectivity(ftgt,aa) .NE. vv) then 
                    vi = surface_mesh%connectivity(ftgt,aa)
                    exit
                end if 
            end do 
            vj = 0
            do aa=1,3
                if ((surface_mesh%connectivity(ftgt,aa) .NE. vv) .AND. (surface_mesh%connectivity(ftgt,aa) .NE. vi)) then 
                    vj = surface_mesh%connectivity(ftgt,aa)
                    exit
                end if 
            end do 

            !Vectors
            vt1(:) = surface_mesh%vertices(vi,:) - surface_mesh%vertices(vv,:)
            vt2(:) = surface_mesh%vertices(vj,:) - surface_mesh%vertices(vv,:)

            !Find angle vi-vv-vj 
            angtot = angtot + ang_vec2vec(vt1,vt2)
        end if 
    end do 
    surface_mesh%vtx_gausscurv(vv) = (2.0d0*pi - angtot)/Abrc
end do 

!Find principal curvatures 
allocate(surface_mesh%vtx_k1(surface_mesh%nvtx))
allocate(surface_mesh%vtx_k2(surface_mesh%nvtx))
do vv=1,surface_mesh%nvtx 
    surface_mesh%vtx_k1(vv) = surface_mesh%vtx_meancurv(vv) + &
    sqrt(abs(surface_mesh%vtx_meancurv(vv)**2 - surface_mesh%vtx_gausscurv(vv)))
    surface_mesh%vtx_k2(vv) = surface_mesh%vtx_meancurv(vv) - &
    sqrt(abs(surface_mesh%vtx_meancurv(vv)**2 - surface_mesh%vtx_gausscurv(vv)))
end do 

!Construct radius of curvature using the average of the magnitudes of the principal curvatures 
allocate(surface_mesh%vtx_rcurv(surface_mesh%nvtx))
surface_mesh%vtx_rcurv(:) = 0.0d0 
do vv=1,surface_mesh%nvtx 
    kl = 0.5d0*(abs(surface_mesh%vtx_k1(vv)) + abs(surface_mesh%vtx_k2(vv)))
    surface_mesh%vtx_rcurv(vv) = (1.0d0/(kl*cm3dopt%surfRcurvM + 1e-12_dp))
end do 

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

!Debug export curvatures 
! open(11,file='io/vrcurve.dat')
!     do vv=1,surface_mesh%nvtx 
!         ! write(11,*) surface_mesh%vtx_rcurv(vv) 
!         ! write(11,*) surface_mesh%vtx_meancurv(vv) 
!         ! write(11,*) surface_mesh%vtx_gausscurv(vv) 
!         write(11,*) surface_mesh%vtx_k1(vv),surface_mesh%vtx_k2(vv),surface_mesh%vtx_meancurv(vv),&
!         surface_mesh%vtx_gausscurv(vv),surface_mesh%vtx_rcurv(vv)
!     end do 
! close(11)
return 
end subroutine evaluate_surf_rcurv




!Construct surface normal vectors ===========================
subroutine construct_surface_normals(surface_mesh)
implicit none 

!Variables - Import 
type(surface_data) :: surface_mesh

!Variables - Local 
integer(in) :: ff,vv
integer(in) :: v1,v2,v3
real(dp) :: normag
real(dp) :: vt1(3),vt2(3)

!Construct face normals 
allocate(surface_mesh%face_normal(surface_mesh%nfcs,3))
do ff=1,surface_mesh%nfcs

    !Vertices
    v1 = surface_mesh%connectivity(ff,1)
    v2 = surface_mesh%connectivity(ff,2)
    v3 = surface_mesh%connectivity(ff,3)

    !Vectors
    vt1(:) = surface_mesh%vertices(v2,:) - surface_mesh%vertices(v1,:)
    vt2(:) = surface_mesh%vertices(v3,:) - surface_mesh%vertices(v1,:)

    !Normal
    surface_mesh%face_normal(ff,:) = crossp3(vt1,vt2)
    normag = norm2(surface_mesh%face_normal(ff,:))
    if (normag .NE. 0.0d0) then 
        surface_mesh%face_normal(ff,:) = surface_mesh%face_normal(ff,:)/normag
    else
        surface_mesh%face_normal(ff,:) = 0.0d0 
        print *, '** zero area surface face ',ff
    end if 
end do

!Construct vertex normals 
allocate(surface_mesh%vtx_normal(surface_mesh%nvtx,3))
surface_mesh%vtx_normal(:,:) = 0.0d0 
do ff=1,surface_mesh%nfcs

    !Vertices
    v1 = surface_mesh%connectivity(ff,1)
    v2 = surface_mesh%connectivity(ff,2)
    v3 = surface_mesh%connectivity(ff,3)

    !Accumulate
    surface_mesh%vtx_normal(v1,:) = surface_mesh%vtx_normal(v1,:) + surface_mesh%face_normal(ff,:)
    surface_mesh%vtx_normal(v2,:) = surface_mesh%vtx_normal(v2,:) + surface_mesh%face_normal(ff,:)
    surface_mesh%vtx_normal(v3,:) = surface_mesh%vtx_normal(v3,:) + surface_mesh%face_normal(ff,:)
end do 
do vv=1,surface_mesh%nvtx
    normag = norm2(surface_mesh%vtx_normal(vv,:))
    if (normag .NE. 0.0d0) then 
        surface_mesh%vtx_normal(vv,:) = surface_mesh%vtx_normal(vv,:)/normag
    else
        surface_mesh%vtx_normal(vv,:) = 0.0d0 
    end if 
end do 
return 
end subroutine construct_surface_normals


end module cellmesh3d_surface_mod