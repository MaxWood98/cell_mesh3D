!cell_mesh3d gradient coupling module
!Max Wood - mw16116@bristol.ac.uk
!Univeristy of Bristol - Department of Aerospace Engineering

!Version 1.1
!Updated 07-11-2023

!Module
module cellmesh3d_gradient_coupling_mod
use cellmesh3d_linalg_mod
use cellmesh3d_adtree_mod
use cellmesh3d_geometry_mod
use cellmesh3d_connectivity_mod
contains 

!Gradient coupling matrix construction subroutine ===========================
subroutine construct_surfvol_grad_coupling(volume_mesh,surface_mesh,Ndim,node_minDIVsize,global_target_pad,cm3dopt)
implicit none 

!Variables - Import
integer(in) :: Ndim,node_minDIVsize
real(dp) :: global_target_pad
type(surface_data) :: surface_mesh
type(vol_mesh_data) :: volume_mesh
type(cm3d_options) :: cm3dopt

!Variables - Local 
integer(in) :: ii,ff,vv,nn,kk,aa,pp
integer(in) :: Nfsurf,ftgt,etgt,vtgt,vtgta,nselected,maxValence,Npts,Npts_ss
integer(in) :: fminD,vtx2,vtx3,nvsurf_sel,nvsurf_selN,nadd,vsins,vselect
integer(in) :: face_surf(volume_mesh%nface),vsurf_gnn(cm3dopt%glink_nnn),smoothinterp(cm3dopt%glink_nsmooth+1)
integer(in) :: vsurf_select(volume_mesh%nvtx),vsurf_dtagged(volume_mesh%nvtx)
integer(in) :: vsurf_tag(volume_mesh%nvtx),vsurf_tag_sm(surface_mesh%nvtx)
integer(in), dimension(:), allocatable :: node_select,valence
integer(in), dimension(:,:), allocatable :: V2F
real(dp) :: lpad,zxmin,zxmax,zymin,zymax,zzmin,zzmax,dist_c,dist_n,vs_maxd,Rs
real(dp) :: vf1(3),vf2(3),vf3(3),vid(4),vsurf_dist(volume_mesh%nvtx)
real(dp) :: R(cm3dopt%glink_nnn,cm3dopt%glink_nnn),Ri(cm3dopt%glink_nnn,cm3dopt%glink_nnn)
real(dp) :: surf2volRBF(cm3dopt%glink_nnn),smooth_s(cm3dopt%glink_nsmooth+1)
real(dp), dimension(:,:), allocatable :: tvtx
type(face_data), dimension(:), allocatable :: faces_surf
type(tree_data) :: sv_adtree

!Allocate surface interpolation structure 
allocate(volume_mesh%vsinterp(surface_mesh%nvtx))

!Extract surface mesh within volume mesh 
Nfsurf = 0 
face_surf(:) = 0
do ff=1,volume_mesh%nface
    if ((volume_mesh%faces(ff)%cleft == -1) .OR. (volume_mesh%faces(ff)%cright == -1)) then 
        Nfsurf = Nfsurf + 1
        face_surf(Nfsurf) = ff
    end if 
end do 
allocate(faces_surf(Nfsurf))
do ff=1,Nfsurf
    ftgt = face_surf(ff)
    faces_surf(ff)%nvtx = volume_mesh%faces(ftgt)%nvtx
    faces_surf(ff)%cleft = volume_mesh%faces(ftgt)%cleft
    faces_surf(ff)%cright = volume_mesh%faces(ftgt)%cright
    allocate(faces_surf(ff)%vertices(volume_mesh%faces(ftgt)%nvtx))
    allocate(faces_surf(ff)%edges(volume_mesh%faces(ftgt)%nvtx))
    faces_surf(ff)%vertices(:) = volume_mesh%faces(ftgt)%vertices(:)
    faces_surf(ff)%edges(:) = 0 
end do 

!Build V2F on this surface 
call get_arb_valence(valence,maxValence,faces_surf,Nfsurf,volume_mesh%nvtx)
call build_V2F(V2F,maxValence,faces_surf,Nfsurf,volume_mesh%nvtx)

!Build adtree on this volume surface mesh 
if (cm3dopt%dispt == 1) then
    write(*,'(A)') '    {constructing AD-tree on surface within volume mesh}'
end if
allocate(tvtx(Nfsurf,6))
do ff=1,Nfsurf
    tvtx(ff,1) = minval(volume_mesh%vtx(faces_surf(ff)%vertices(:),1)) !xmin
    tvtx(ff,2) = minval(volume_mesh%vtx(faces_surf(ff)%vertices(:),2)) !ymin
    tvtx(ff,3) = minval(volume_mesh%vtx(faces_surf(ff)%vertices(:),3)) !zmin
    tvtx(ff,4) = maxval(volume_mesh%vtx(faces_surf(ff)%vertices(:),1)) !xmax
    tvtx(ff,5) = maxval(volume_mesh%vtx(faces_surf(ff)%vertices(:),2)) !ymax
    tvtx(ff,6) = maxval(volume_mesh%vtx(faces_surf(ff)%vertices(:),3)) !zmax
end do
call build_ADtree(sv_adtree,Ndim,cm3dopt%max_depth,node_minDIVsize,tvtx,global_target_pad,cm3dopt%dispt)
nselected = 0
allocate(node_select(sv_adtree%nnode))
node_select(:) = 0 

!Display
if (cm3dopt%dispt == 1) then
    write(*,'(A)') '    {constructing volume to surface interpolations}'
end if

!Build local interpolation for each surface point
zxmin = 0.0d0 
zxmax = 0.0d0 
zymin = 0.0d0 
zymax = 0.0d0 
zzmin = 0.0d0 
zzmax = 0.0d0
R(:,:) = 0.0d0
Ri(:,:) = 0.0d0
vsurf_tag(:) = 0 
vsurf_select(:) = 0 
vsurf_dtagged(:) = 0
vsurf_dist(:) = 0.0d0  
do vv=1,surface_mesh%nvtx

    !Find the distance to the farthest connected surface vertex from this vertex
    dist_c = 0.0d0 
    do aa=1,surface_mesh%valence(vv)
        etgt = surface_mesh%V2E(vv,aa)
        dist_n = norm2(surface_mesh%vertices(surface_mesh%edges(etgt,2),:) - surface_mesh%vertices(surface_mesh%edges(etgt,1),:))
        if (dist_n .GT. dist_c) then 
            dist_c = dist_n
        end if 
    end do 
    
    !Set search length as 4x the distance to the farthes other surface vertex
    lpad = 4.0d0*dist_c

    !Find the closest volume mesh surface face to the current surface mesh vertex
    zxmin = surface_mesh%vertices(vv,1) - lpad !tgt bounding box -> xmin
    zxmax = surface_mesh%vertices(vv,1) + lpad !tgt bounding box -> xmax
    zymin = surface_mesh%vertices(vv,2) - lpad !tgt bounding box -> ymin
    zymax = surface_mesh%vertices(vv,2) + lpad !tgt bounding box -> ymax
    zzmin = surface_mesh%vertices(vv,3) - lpad !tgt bounding box -> zmin
    zzmax = surface_mesh%vertices(vv,3) + lpad !tgt bounding box -> zmax
    call search_ADtree(nselected,node_select,sv_adtree,zxmin,zxmax,zymin,zymax,zzmin,zzmax)
    fminD = 0
    dist_c = 2.0d0*cm3dopt%far_field_bound
    do nn=1,nselected
        do kk=1,sv_adtree%tree(node_select(nn))%nentry

            !Target face
            ftgt = sv_adtree%tree(node_select(nn))%entry(kk)

            !Check triangulation of this face 
            vf1(:) = volume_mesh%vtx(faces_surf(ftgt)%vertices(1),:)
            do aa=2,faces_surf(ftgt)%nvtx-1

                !Triangle on this face 
                vtx2 = aa 
                vtx3 = aa + 1
                vtx2 = faces_surf(ftgt)%vertices(vtx2)
                vtx3 = faces_surf(ftgt)%vertices(vtx3)
                vf2(:) = volume_mesh%vtx(vtx2,:)
                vf3(:) = volume_mesh%vtx(vtx3,:)

                !Test distance from surface point to this face
                vid = min_dist_point_to_tri(surface_mesh%vertices(vv,:),vf1,vf2,vf3)
                dist_n = vid(4)
                if (dist_n .LE. dist_c) then 
                    fminD = ftgt
                    dist_c = dist_n
                end if
            end do 
        end do 
    end do 

    !If non zero target face (else if zero then this face is outside the mesh domain so can be ignored)
    if (fminD .NE. 0) then !build interpolation structure 

        !Add base face to surface tag structure 
        nvsurf_sel = 0 
        do aa=1,faces_surf(fminD)%nvtx
            nvsurf_sel = nvsurf_sel + 1
            vsurf_select(nvsurf_sel) = faces_surf(fminD)%vertices(aa)
            vsurf_dist(nvsurf_sel) = norm2(volume_mesh%vtx(faces_surf(fminD)%vertices(aa),:) - surface_mesh%vertices(vv,:))
            vsurf_tag(faces_surf(fminD)%vertices(aa)) = 1
        end do 
        vs_maxd = maxval(vsurf_dist(1:nvsurf_sel))

        !Find glink_nnn nearest neighbours superset
        do ff=1,volume_mesh%nvtx !Flood to add valid adjacent new points 
            nadd = 0 
            nvsurf_selN = nvsurf_sel
            do pp=1,nvsurf_sel

                !Target vertex
                vtgt = vsurf_select(pp)

                !Each face on this vertex 
                do aa=1,valence(vtgt)

                    !Target face 
                    ftgt = V2F(vtgt,aa)

                    !All vertices on this face 
                    if (ftgt .GT. 0) then 
                        do kk=1,faces_surf(ftgt)%nvtx
                            vtgta = faces_surf(ftgt)%vertices(kk) 
                            if (vsurf_tag(vtgta) == 0) then !new vertex
                                if (nvsurf_sel .LE. cm3dopt%glink_nnn) then !If not enough points yet found 
                                    nadd = nadd + 1
                                    nvsurf_selN = nvsurf_selN + 1
                                    vsurf_select(nvsurf_selN) = vtgta
                                    vsurf_dist(nvsurf_selN) = norm2(volume_mesh%vtx(vtgta,:) - surface_mesh%vertices(vv,:))
                                    vsurf_tag(vtgta) = 1
                                    if (vsurf_dist(nvsurf_selN) .GT. vs_maxd) then !Update maximum distance if required
                                        vs_maxd = vsurf_dist(nvsurf_selN)
                                    end if
                                elseif (norm2(volume_mesh%vtx(vtgta,:) - surface_mesh%vertices(vv,:)) .LE. vs_maxd) then !If enough points have been found but this point is closer than the current maximum distance point
                                    nadd = nadd + 1
                                    nvsurf_selN = nvsurf_selN + 1
                                    vsurf_select(nvsurf_selN) = vtgta
                                    vsurf_dist(nvsurf_selN) = norm2(volume_mesh%vtx(vtgta,:) - surface_mesh%vertices(vv,:))
                                    vsurf_tag(vtgta) = 1
                                    if (vsurf_dist(nvsurf_selN) .GT. vs_maxd) then !Update maximum distance if required
                                        vs_maxd = vsurf_dist(nvsurf_selN)
                                    end if
                                end if 
                            end if 
                        end do 
                    end if 
                end do 
            end do 

            !Set new point count 
            nvsurf_sel = nvsurf_selN

            !Exit if no new points added 
            if (nadd == 0) then 
                exit 
            end if 
        end do 
        vsurf_tag(vsurf_select(1:nvsurf_sel)) = 0 !Reset vertex tags 

        !Find actual glink_nnn nearest neighbours
        vsins = 0 
        vsurf_gnn(:) = 0 
        vsurf_dtagged(1:nvsurf_sel) = 0 
        do pp=1,nvsurf_sel
            vselect = minloc(vsurf_dist(1:nvsurf_sel),1,mask = vsurf_dtagged(1:nvsurf_sel) == 0)
            vsurf_dtagged(vselect) = 1
            vtgt = vsurf_select(vselect)
            vsins = vsins + 1
            vsurf_gnn(vsins) = vtgt 
            if (vsins == cm3dopt%glink_nnn) then !Exit when all points are found 
                exit 
            end if 
        end do 

        !Reset arrays 
        Ri(:,:) = 0.0d0
        surf2volRBF(:) = 0.0d0 

        !Cases
        if (vsins .LT. cm3dopt%glink_nnn) then !If too few point have been found 

            !Build local dependance matrix and set support radius 
            call build_RBF_influence(R,Rs,vsins,vsurf_gnn,volume_mesh%vtx,cm3dopt)

            !Invert dependance matrix 
            call matinv(R(1:vsins,1:vsins),Ri(1:vsins,1:vsins),vsins)

            !Build RBF dependance of surface point on all selected volume points 
            do ii=1,vsins
                dist_c = norm2(volume_mesh%vtx(vsurf_gnn(ii),:) - surface_mesh%vertices(vv,:))
                surf2volRBF(ii) = wendlandc2(dist_c,Rs)
            end do 

            !Set number of points 
            Npts = vsins
        else !If full number of points found
            if (cm3dopt%glink_nnn .GT. 1) then !More than one surface point requested

                !Build local dependance matrix and set support radius 
                call build_RBF_influence(R,Rs,cm3dopt%glink_nnn,vsurf_gnn,volume_mesh%vtx,cm3dopt)

                !Invert dependance matrix 
                call matinv(R,Ri,cm3dopt%glink_nnn)

                !Build RBF dependance of surface point on all selected volume points 
                do ii=1,cm3dopt%glink_nnn
                    dist_c = norm2(volume_mesh%vtx(vsurf_gnn(ii),:) - surface_mesh%vertices(vv,:))
                    surf2volRBF(ii) = wendlandc2(dist_c,Rs)
                end do 

                !Set number of points 
                Npts = cm3dopt%glink_nnn 
            else !One surface point requested
                surf2volRBF(1) = 1.0d0 
                Ri(1,1) = 1.0d0 
                Npts = 1
            end if 
        end if 

        !Store interpolation structure 
        volume_mesh%vsinterp(vv)%npnts_vi = Npts
        allocate(volume_mesh%vsinterp(vv)%vol_pnts(cm3dopt%glink_nnn))
        allocate(volume_mesh%vsinterp(vv)%surf2volRBF(cm3dopt%glink_nnn))
        allocate(volume_mesh%vsinterp(vv)%Ri(cm3dopt%glink_nnn,cm3dopt%glink_nnn))
        volume_mesh%vsinterp(vv)%vol_pnts(:) = 0.0d0 
        volume_mesh%vsinterp(vv)%surf2volRBF(:) = 0.0d0 
        volume_mesh%vsinterp(vv)%Ri(:,:) = 0.0d0 
        volume_mesh%vsinterp(vv)%vol_pnts(1:Npts) = vsurf_gnn(1:Npts)
        volume_mesh%vsinterp(vv)%surf2volRBF(1:Npts) = surf2volRBF(1:Npts)
        volume_mesh%vsinterp(vv)%Ri(1:Npts,1:Npts) = Ri(1:Npts,1:Npts)
    else !set unused interpolation structure 
        volume_mesh%vsinterp(vv)%npnts_vi = 0 
        allocate(volume_mesh%vsinterp(vv)%vol_pnts(cm3dopt%glink_nnn))
        allocate(volume_mesh%vsinterp(vv)%surf2volRBF(cm3dopt%glink_nnn))
        allocate(volume_mesh%vsinterp(vv)%Ri(cm3dopt%glink_nnn,cm3dopt%glink_nnn))
        volume_mesh%vsinterp(vv)%vol_pnts(:) = 0
        volume_mesh%vsinterp(vv)%surf2volRBF(:) = 0
        volume_mesh%vsinterp(vv)%Ri(:,:) = 0.0d0 
    end if 
end do 

!Display
if (cm3dopt%dispt == 1) then
    write(*,'(A)') '    {constructing surface interpolation smoothing}'
end if

!Add linked vertices to smooth each surface vertex 
vsurf_tag_sm(:) = 0 
do vv=1,surface_mesh%nvtx
    if (volume_mesh%vsinterp(vv)%npnts_vi .GT. 0) then !build smoothing structure

        !Add base vertex surface tag structure 
        nvsurf_sel = 1 
        vsurf_select(nvsurf_sel) = vv 
        vsurf_dist(nvsurf_sel) = 0.0d0 
        vsurf_tag_sm(vv) = 1
        vs_maxd = 0.0d0 

        !Find the superset of cm3dopt%glink_nsmooth+1 vertices near surface vertex vv
        do ff=1,surface_mesh%nvtx !Flood to add valid adjacent new points 
            nadd = 0 
            nvsurf_selN = nvsurf_sel
            do pp=1,nvsurf_sel

                !Target vertex
                vtgt = vsurf_select(pp)

                !Each face on this vertex 
                do aa=1,surface_mesh%valence(vtgt)

                    !Target face 
                    ftgt = surface_mesh%V2F(vtgt,aa)

                    !All vertices on this face 
                    do kk=1,3
                        vtgta = surface_mesh%connectivity(ftgt,kk)
                        if (vsurf_tag_sm(vtgta) == 0) then !new vertex
                            if (nvsurf_sel .LE. cm3dopt%glink_nsmooth+1) then !If not enough points yet found 
                                nadd = nadd + 1
                                nvsurf_selN = nvsurf_selN + 1
                                vsurf_select(nvsurf_selN) = vtgta
                                vsurf_dist(nvsurf_selN) = norm2(surface_mesh%vertices(vtgta,:) - surface_mesh%vertices(vv,:))
                                vsurf_tag_sm(vtgta) = 1
                                if (vsurf_dist(nvsurf_selN) .GT. vs_maxd) then !Update maximum distance if required
                                    vs_maxd = vsurf_dist(nvsurf_selN)
                                end if
                            elseif (norm2(surface_mesh%vertices(vtgta,:) - surface_mesh%vertices(vv,:)) .LE. vs_maxd) then !If enough points have been found but this point is closer than the current maximum distance point
                                nadd = nadd + 1
                                nvsurf_selN = nvsurf_selN + 1
                                vsurf_select(nvsurf_selN) = vtgta
                                vsurf_dist(nvsurf_selN) = norm2(surface_mesh%vertices(vtgta,:) - surface_mesh%vertices(vv,:))
                                vsurf_tag_sm(vtgta) = 1
                                if (vsurf_dist(nvsurf_selN) .GT. vs_maxd) then !Update maximum distance if required
                                    vs_maxd = vsurf_dist(nvsurf_selN)
                                end if
                            end if 
                        end if 
                    end do 
                end do 
            end do 

            !Set new point count 
            nvsurf_sel = nvsurf_selN

            !Exit if no new points added 
            if (nadd == 0) then 
                exit 
            end if 
        end do 
        vsurf_tag_sm(vsurf_select(1:nvsurf_sel)) = 0 !Reset vertex tags 

        !Find actual cm3dopt%glink_nsmooth+1 nearest neighbours
        vsins = 0 
        smoothinterp(:) = 0 
        vsurf_dtagged(1:nvsurf_sel) = 0 
        do pp=1,nvsurf_sel
            vselect = minloc(vsurf_dist(1:nvsurf_sel),1,mask = vsurf_dtagged(1:nvsurf_sel) == 0)
            vsurf_dtagged(vselect) = 1
            vtgt = vsurf_select(vselect)
            vsins = vsins + 1
            smoothinterp(vsins) = vtgt 
            if (vsins == cm3dopt%glink_nsmooth+1) then !Exit when all points are found 
                exit 
            end if 
        end do 

        !Set the number of points found for use in surface smoothing 
        Npts_ss = vsins
        volume_mesh%vsinterp(vv)%npnts_ss = vsins

        !Find the distance to each surface smoothing point 
        smooth_s(:) = 0.0d0 
        do aa=1,Npts_ss
            smooth_s(aa) = norm2(surface_mesh%vertices(smoothinterp(aa),:) - surface_mesh%vertices(vv,:))
        end do 

        !Allocate and store smoothing points
        allocate(volume_mesh%vsinterp(vv)%surf_smooth_pnts(cm3dopt%glink_nsmooth+1))
        allocate(volume_mesh%vsinterp(vv)%surf_smoothRBF(cm3dopt%glink_nsmooth+1))
        volume_mesh%vsinterp(vv)%surf_smoothRBF(:) = 0.0d0 
        volume_mesh%vsinterp(vv)%surf_smooth_pnts(:) = smoothinterp(:)

        !Evaluate and store RBF weighting for each point on the surface
        if (cm3dopt%glink_nsmooth+1 .GE. 1) then 
            Rs = 1.1d0*maxval(smooth_s(:))
            do aa=1,Npts_ss
                volume_mesh%vsinterp(vv)%surf_smoothRBF(aa) = wendlandc2(smooth_s(aa),Rs)
            end do 
        else
            volume_mesh%vsinterp(vv)%surf_smoothRBF(1) = 1.0d0 
        end if 
    else !set unused smoothing structure 
        allocate(volume_mesh%vsinterp(vv)%surf_smooth_pnts(cm3dopt%glink_nsmooth+1))
        allocate(volume_mesh%vsinterp(vv)%surf_smoothRBF(cm3dopt%glink_nsmooth+1))
        volume_mesh%vsinterp(vv)%surf_smooth_pnts(:) = 0
        volume_mesh%vsinterp(vv)%surf_smoothRBF(:) = 0 
    end if 
end do 
return 
end subroutine construct_surfvol_grad_coupling




!Gradient projection subroutine ===========================
subroutine project_gradients(gradient_surf,gradient_vol,volume_mesh,surface_mesh,Ndim,node_minDIVsize,global_target_pad,cm3dopt)
implicit none 

!Variables - Import
integer(in) :: Ndim,node_minDIVsize
real(dp) :: global_target_pad
real(dp), dimension(:,:) :: gradient_vol
real(dp), dimension(:,:), allocatable :: gradient_surf
type(surface_data) :: surface_mesh
type(vol_mesh_data) :: volume_mesh
type(cm3d_options) :: cm3dopt

!Variables - Local 
integer(in) :: ii,ff,vv,nn,kk,aa,pp
integer(in) :: Nfsurf,ftgt,etgt,vtgt,vtgta,nselected,maxValence,Npts,Npts_ss
integer(in) :: fminD,vtx2,vtx3,nvsurf_sel,nvsurf_selN,nadd,vsins,vselect
integer(in) :: face_surf(volume_mesh%nface),vsurf_gnn(cm3dopt%glink_nnn),smoothinterp(cm3dopt%glink_nsmooth+1)
integer(in) :: vsurf_select(volume_mesh%nvtx),vsurf_dtagged(volume_mesh%nvtx)
integer(in) :: vsurf_tag(volume_mesh%nvtx),vsurf_tag_sm(surface_mesh%nvtx)
integer(in), dimension(:), allocatable :: node_select,valence
integer(in), dimension(:,:), allocatable :: V2F
real(dp) :: lpad,zxmin,zxmax,zymin,zymax,zzmin,zzmax,dist_c,dist_n,vs_maxd,Rs
real(dp) :: vf1(3),vf2(3),vf3(3),vid(4),vsurf_dist(volume_mesh%nvtx)
real(dp) :: R(cm3dopt%glink_nnn,cm3dopt%glink_nnn),Ri(cm3dopt%glink_nnn,cm3dopt%glink_nnn)
real(dp) :: surf2volRBF(cm3dopt%glink_nnn),smooth_s(cm3dopt%glink_nsmooth+1),surf_smoothRBF(cm3dopt%glink_nsmooth+1)
real(dp), dimension(:,:), allocatable :: tvtx,gradient_surf_temp
type(face_data), dimension(:), allocatable :: faces_surf
type(tree_data) :: sv_adtree

!Allocate and initialise surface gradient structures 
allocate(gradient_surf(surface_mesh%nvtx,3))
allocate(gradient_surf_temp(surface_mesh%nvtx,3))
gradient_surf(:,:) = 0.0d0 
gradient_surf_temp(:,:) = 0.0d0 

!Allocate surface interpolation structure 
allocate(volume_mesh%vsinterp(surface_mesh%nvtx))

!Extract surface mesh within volume mesh 
Nfsurf = 0 
face_surf(:) = 0
do ff=1,volume_mesh%nface
    if ((volume_mesh%faces(ff)%cleft == -1) .OR. (volume_mesh%faces(ff)%cright == -1)) then 
        Nfsurf = Nfsurf + 1
        face_surf(Nfsurf) = ff
    end if 
end do 
allocate(faces_surf(Nfsurf))
do ff=1,Nfsurf
    ftgt = face_surf(ff)
    faces_surf(ff)%nvtx = volume_mesh%faces(ftgt)%nvtx
    faces_surf(ff)%cleft = volume_mesh%faces(ftgt)%cleft
    faces_surf(ff)%cright = volume_mesh%faces(ftgt)%cright
    allocate(faces_surf(ff)%vertices(volume_mesh%faces(ftgt)%nvtx))
    allocate(faces_surf(ff)%edges(volume_mesh%faces(ftgt)%nvtx))
    faces_surf(ff)%vertices(:) = volume_mesh%faces(ftgt)%vertices(:)
    faces_surf(ff)%edges(:) = 0 
end do 

!Build V2F on this surface 
call get_arb_valence(valence,maxValence,faces_surf,Nfsurf,volume_mesh%nvtx)
call build_V2F(V2F,maxValence,faces_surf,Nfsurf,volume_mesh%nvtx)

!Build adtree on this volume surface mesh 
if (cm3dopt%dispt == 1) then
    write(*,'(A)') '    {constructing AD-tree on surface within volume mesh}'
end if
allocate(tvtx(Nfsurf,6))
do ff=1,Nfsurf
    tvtx(ff,1) = minval(volume_mesh%vtx(faces_surf(ff)%vertices(:),1)) !xmin
    tvtx(ff,2) = minval(volume_mesh%vtx(faces_surf(ff)%vertices(:),2)) !ymin
    tvtx(ff,3) = minval(volume_mesh%vtx(faces_surf(ff)%vertices(:),3)) !zmin
    tvtx(ff,4) = maxval(volume_mesh%vtx(faces_surf(ff)%vertices(:),1)) !xmax
    tvtx(ff,5) = maxval(volume_mesh%vtx(faces_surf(ff)%vertices(:),2)) !ymax
    tvtx(ff,6) = maxval(volume_mesh%vtx(faces_surf(ff)%vertices(:),3)) !zmax
end do
call build_ADtree(sv_adtree,Ndim,cm3dopt%max_depth,node_minDIVsize,tvtx,global_target_pad,cm3dopt%dispt)
nselected = 0
allocate(node_select(sv_adtree%nnode))
node_select(:) = 0 

!Display
if (cm3dopt%dispt == 1) then
    write(*,'(A)') '    {interpolating gradients to the surface}'
end if

!Build local interpolation for each surface point
zxmin = 0.0d0 
zxmax = 0.0d0 
zymin = 0.0d0 
zymax = 0.0d0 
zzmin = 0.0d0 
zzmax = 0.0d0
R(:,:) = 0.0d0
Ri(:,:) = 0.0d0
vsurf_tag(:) = 0 
vsurf_select(:) = 0 
vsurf_dtagged(:) = 0
vsurf_dist(:) = 0.0d0  
do vv=1,surface_mesh%nvtx

    !Find the distance to the farthest connected surface vertex from this vertex
    dist_c = 0.0d0 
    do aa=1,surface_mesh%valence(vv)
        etgt = surface_mesh%V2E(vv,aa)
        dist_n = norm2(surface_mesh%vertices(surface_mesh%edges(etgt,2),:) - surface_mesh%vertices(surface_mesh%edges(etgt,1),:))
        if (dist_n .GT. dist_c) then 
            dist_c = dist_n
        end if 
    end do 
    
    !Set search length as 4x the distance to the farthest other surface vertex
    lpad = 4.0d0*dist_c

    !Find the closest volume mesh surface face to the current surface mesh vertex
    zxmin = surface_mesh%vertices(vv,1) - lpad !tgt bounding box -> xmin
    zxmax = surface_mesh%vertices(vv,1) + lpad !tgt bounding box -> xmax
    zymin = surface_mesh%vertices(vv,2) - lpad !tgt bounding box -> ymin
    zymax = surface_mesh%vertices(vv,2) + lpad !tgt bounding box -> ymax
    zzmin = surface_mesh%vertices(vv,3) - lpad !tgt bounding box -> zmin
    zzmax = surface_mesh%vertices(vv,3) + lpad !tgt bounding box -> zmax
    call search_ADtree(nselected,node_select,sv_adtree,zxmin,zxmax,zymin,zymax,zzmin,zzmax)
    fminD = 0
    dist_c = 2.0d0*cm3dopt%far_field_bound
    do nn=1,nselected
        do kk=1,sv_adtree%tree(node_select(nn))%nentry

            !Target face
            ftgt = sv_adtree%tree(node_select(nn))%entry(kk)

            !Check triangulation of this face 
            vf1(:) = volume_mesh%vtx(faces_surf(ftgt)%vertices(1),:)
            do aa=2,faces_surf(ftgt)%nvtx-1

                !Triangle on this face 
                vtx2 = aa 
                vtx3 = aa + 1
                vtx2 = faces_surf(ftgt)%vertices(vtx2)
                vtx3 = faces_surf(ftgt)%vertices(vtx3)
                vf2(:) = volume_mesh%vtx(vtx2,:)
                vf3(:) = volume_mesh%vtx(vtx3,:)

                !Test distance from surface point to this face
                vid = min_dist_point_to_tri(surface_mesh%vertices(vv,:),vf1,vf2,vf3)
                dist_n = vid(4)
                if (dist_n .LE. dist_c) then 
                    fminD = ftgt
                    dist_c = dist_n
                end if
            end do 
        end do 
    end do 

    !If non zero target face (else if zero then this face is outside the mesh domain so can be ignored)
    if (fminD .NE. 0) then !build interpolation structure 

        !Add base face to surface tag structure 
        nvsurf_sel = 0 
        do aa=1,faces_surf(fminD)%nvtx
            nvsurf_sel = nvsurf_sel + 1
            vsurf_select(nvsurf_sel) = faces_surf(fminD)%vertices(aa)
            vsurf_dist(nvsurf_sel) = norm2(volume_mesh%vtx(faces_surf(fminD)%vertices(aa),:) - surface_mesh%vertices(vv,:))
            vsurf_tag(faces_surf(fminD)%vertices(aa)) = 1
        end do 
        vs_maxd = maxval(vsurf_dist(1:nvsurf_sel))

        !Find glink_nnn nearest neighbours superset
        do ff=1,volume_mesh%nvtx !Flood to add valid adjacent new points 
            nadd = 0 
            nvsurf_selN = nvsurf_sel
            do pp=1,nvsurf_sel

                !Target vertex
                vtgt = vsurf_select(pp)

                !Each face on this vertex 
                do aa=1,valence(vtgt)

                    !Target face 
                    ftgt = V2F(vtgt,aa)

                    !All vertices on this face 
                    if (ftgt .GT. 0) then 
                        do kk=1,faces_surf(ftgt)%nvtx
                            vtgta = faces_surf(ftgt)%vertices(kk) 
                            if (vsurf_tag(vtgta) == 0) then !new vertex
                                if (nvsurf_sel .LE. cm3dopt%glink_nnn) then !If not enough points yet found 
                                    nadd = nadd + 1
                                    nvsurf_selN = nvsurf_selN + 1
                                    vsurf_select(nvsurf_selN) = vtgta
                                    vsurf_dist(nvsurf_selN) = norm2(volume_mesh%vtx(vtgta,:) - surface_mesh%vertices(vv,:))
                                    vsurf_tag(vtgta) = 1
                                    if (vsurf_dist(nvsurf_selN) .GT. vs_maxd) then !Update maximum distance if required
                                        vs_maxd = vsurf_dist(nvsurf_selN)
                                    end if
                                elseif (norm2(volume_mesh%vtx(vtgta,:) - surface_mesh%vertices(vv,:)) .LE. vs_maxd) then !If enough points have been found but this point is closer than the current maximum distance point
                                    nadd = nadd + 1
                                    nvsurf_selN = nvsurf_selN + 1
                                    vsurf_select(nvsurf_selN) = vtgta
                                    vsurf_dist(nvsurf_selN) = norm2(volume_mesh%vtx(vtgta,:) - surface_mesh%vertices(vv,:))
                                    vsurf_tag(vtgta) = 1
                                    if (vsurf_dist(nvsurf_selN) .GT. vs_maxd) then !Update maximum distance if required
                                        vs_maxd = vsurf_dist(nvsurf_selN)
                                    end if
                                end if 
                            end if 
                        end do 
                    end if 
                end do 
            end do 

            !Set new point count 
            nvsurf_sel = nvsurf_selN

            !Exit if no new points added 
            if (nadd == 0) then 
                exit 
            end if 
        end do 
        vsurf_tag(vsurf_select(1:nvsurf_sel)) = 0 !Reset vertex tags 

        !Find actual glink_nnn nearest neighbours
        vsins = 0 
        vsurf_gnn(:) = 0 
        vsurf_dtagged(1:nvsurf_sel) = 0 
        do pp=1,nvsurf_sel
            vselect = minloc(vsurf_dist(1:nvsurf_sel),1,mask = vsurf_dtagged(1:nvsurf_sel) == 0)
            vsurf_dtagged(vselect) = 1
            vtgt = vsurf_select(vselect)
            vsins = vsins + 1
            vsurf_gnn(vsins) = vtgt 
            if (vsins == cm3dopt%glink_nnn) then !Exit when all points are found 
                exit 
            end if 
        end do 

        !Reset arrays 
        Ri(:,:) = 0.0d0
        surf2volRBF(:) = 0.0d0 

        !Cases
        if (vsins .LT. cm3dopt%glink_nnn) then !If too few point have been found 

            !Build local dependance matrix and set support radius 
            call build_RBF_influence(R,Rs,vsins,vsurf_gnn,volume_mesh%vtx,cm3dopt)

            !Invert dependance matrix 
            call matinv(R(1:vsins,1:vsins),Ri(1:vsins,1:vsins),vsins)

            !Build RBF dependance of surface point on all selected volume points 
            do ii=1,vsins
                dist_c = norm2(volume_mesh%vtx(vsurf_gnn(ii),:) - surface_mesh%vertices(vv,:))
                surf2volRBF(ii) = wendlandc2(dist_c,Rs)
            end do 

            !Set number of points 
            Npts = vsins
        else !If full number of points found
            if (cm3dopt%glink_nnn .GT. 1) then !More than one surface point requested

                !Build local dependance matrix and set support radius 
                call build_RBF_influence(R,Rs,cm3dopt%glink_nnn,vsurf_gnn,volume_mesh%vtx,cm3dopt)

                !Invert dependance matrix 
                call matinv(R,Ri,cm3dopt%glink_nnn)

                !Build RBF dependance of surface point on all selected volume points 
                do ii=1,cm3dopt%glink_nnn
                    dist_c = norm2(volume_mesh%vtx(vsurf_gnn(ii),:) - surface_mesh%vertices(vv,:))
                    surf2volRBF(ii) = wendlandc2(dist_c,Rs)
                end do 

                !Set number of points 
                Npts = cm3dopt%glink_nnn 
            else !One surface point requested
                surf2volRBF(1) = 1.0d0 
                Ri(1,1) = 1.0d0 
                Npts = 1
            end if 
        end if 

        !Project to this surface vertex 
        gradient_surf(vv,:) = project_vtx_vgrad_2_surfgrad(gradient_vol,surf2volRBF(1:Npts),Ri,vsurf_gnn(1:Npts),Npts)
        volume_mesh%vsinterp(vv)%npnts_vi = Npts
    else !no face found so ignore projection
        !do nothing 
    end if 
end do 

!Return if no smoothing 
if (cm3dopt%glink_nsmooth == 0) then 
    return 
end if

!Store surface gradients in the temporary structure
gradient_surf_temp(:,:) = gradient_surf(:,:)

!Display
if (cm3dopt%dispt == 1) then
    write(*,'(A)') '    {smoothing surface gradients}'
end if

!Add linked vertices to smooth each surface vertex 
vsurf_tag_sm(:) = 0
do vv=1,surface_mesh%nvtx
    if (volume_mesh%vsinterp(vv)%npnts_vi .GT. 0) then !smooth this vertex

        !Add base vertex surface tag structure 
        nvsurf_sel = 1 
        vsurf_select(nvsurf_sel) = vv 
        vsurf_dist(nvsurf_sel) = 0.0d0 
        vsurf_tag_sm(vv) = 1
        vs_maxd = 0.0d0 

        !Find the superset of cm3dopt%glink_nsmooth+1 vertices near surface vertex vv
        do ff=1,surface_mesh%nvtx !Flood to add valid adjacent new points 
            nadd = 0 
            nvsurf_selN = nvsurf_sel
            do pp=1,nvsurf_sel

                !Target vertex
                vtgt = vsurf_select(pp)

                !Each face on this vertex 
                do aa=1,surface_mesh%valence(vtgt)

                    !Target face 
                    ftgt = surface_mesh%V2F(vtgt,aa)

                    !All vertices on this face 
                    do kk=1,3
                        vtgta = surface_mesh%connectivity(ftgt,kk)
                        if (vsurf_tag_sm(vtgta) == 0) then !new vertex
                            if (nvsurf_sel .LE. cm3dopt%glink_nsmooth+1) then !If not enough points yet found 
                                nadd = nadd + 1
                                nvsurf_selN = nvsurf_selN + 1
                                vsurf_select(nvsurf_selN) = vtgta
                                vsurf_dist(nvsurf_selN) = norm2(surface_mesh%vertices(vtgta,:) - surface_mesh%vertices(vv,:))
                                vsurf_tag_sm(vtgta) = 1
                                if (vsurf_dist(nvsurf_selN) .GT. vs_maxd) then !Update maximum distance if required
                                    vs_maxd = vsurf_dist(nvsurf_selN)
                                end if
                            elseif (norm2(surface_mesh%vertices(vtgta,:) - surface_mesh%vertices(vv,:)) .LE. vs_maxd) then !If enough points have been found but this point is closer than the current maximum distance point
                                nadd = nadd + 1
                                nvsurf_selN = nvsurf_selN + 1
                                vsurf_select(nvsurf_selN) = vtgta
                                vsurf_dist(nvsurf_selN) = norm2(surface_mesh%vertices(vtgta,:) - surface_mesh%vertices(vv,:))
                                vsurf_tag_sm(vtgta) = 1
                                if (vsurf_dist(nvsurf_selN) .GT. vs_maxd) then !Update maximum distance if required
                                    vs_maxd = vsurf_dist(nvsurf_selN)
                                end if
                            end if 
                        end if 
                    end do 
                end do 
            end do 

            !Set new point count 
            nvsurf_sel = nvsurf_selN

            !Exit if no new points added 
            if (nadd == 0) then 
                exit 
            end if 
        end do 
        vsurf_tag_sm(vsurf_select(1:nvsurf_sel)) = 0 !Reset vertex tags 

        !Find actual cm3dopt%glink_nsmooth+1 nearest neighbours
        vsins = 0 
        smoothinterp(:) = 0 
        vsurf_dtagged(1:nvsurf_sel) = 0 
        do pp=1,nvsurf_sel
            vselect = minloc(vsurf_dist(1:nvsurf_sel),1,mask = vsurf_dtagged(1:nvsurf_sel) == 0)
            vsurf_dtagged(vselect) = 1
            vtgt = vsurf_select(vselect)
            vsins = vsins + 1
            smoothinterp(vsins) = vtgt 
            if (vsins == cm3dopt%glink_nsmooth+1) then !Exit when all points are found 
                exit 
            end if 
        end do 

        !Set the number of points found for use in surface smoothing 
        Npts_ss = vsins
        volume_mesh%vsinterp(vv)%npnts_ss = vsins

        !Find the distance to each surface smoothing point 
        smooth_s(:) = 0.0d0 
        do aa=1,Npts_ss
            smooth_s(aa) = norm2(surface_mesh%vertices(smoothinterp(aa),:) - surface_mesh%vertices(vv,:))
        end do 

        !Build RBF values at smoothing points 
        surf_smoothRBF(:) = 0.0d0 
        if (cm3dopt%glink_nsmooth+1 .GE. 1) then 
            Rs = 1.1d0*maxval(smooth_s(:))
            do aa=1,Npts_ss
                surf_smoothRBF(aa) = wendlandc2(smooth_s(aa),Rs)
            end do 
        else
            surf_smoothRBF(1) = 1.0d0 
        end if 

        !Apply smoothing to this vertex 
        gradient_surf(vv,:) = smooth_surface_gradient(gradient_surf_temp,surf_smoothRBF,smoothinterp,Npts_ss)
    else !no interpolation here so ignore smoothing 
        !do nothing 
    end if 
end do 
return 
end subroutine project_gradients




!Subroutine to build RBF influence ===========================
subroutine build_RBF_influence(R,R_sup,Npoint,point_list,vertices,cm3dopt)
implicit none 

!Variables - Import
integer(in) :: Npoint 
integer(in), dimension(:) :: point_list
real(dp) :: R_sup
real(dp), dimension(:,:) :: R,vertices 
type(cm3d_options) :: cm3dopt

!Variables - Local
integer(in) :: ii,jj 
integer(in) :: v1,v2
real(dp) :: mindist(Npoint)

!Populate distance matrix
R(:,:) = 0.0d0 
mindist(:) = ieee_value(1.0d0,IEEE_POSITIVE_INF)
do ii=1,Npoint
    v1 = point_list(ii)
    do jj=1,Npoint
        v2 = point_list(jj)
        R(ii,jj) = norm2(vertices(v2,:) - vertices(v1,:)) 
        if ((R(ii,jj) .GT. 0.0d0) .AND. (R(ii,jj) .LT. mindist(ii))) then 
            mindist(ii) = R(ii,jj)
        end if 
    end do 
end do 

!Set support radius 
R_sup = 50.0d0*maxval(R(1:Npoint,1:Npoint))

!Evaluate RBF values for each distance entry
do ii=1,Npoint
    do jj=1,Npoint
        R(ii,jj) = wendlandc2(R(ii,jj),R_sup)
    end do
end do 

!Relax interpolation at interior points 
if ((Npoint .GT. 2) .AND. (cm3dopt%RBF_relax .NE. 0.0d0)) then 
    do ii=2,Npoint-1
        R(ii,ii) = R(ii,ii) + cm3dopt%RBF_relax/(mindist(ii)*R_sup)
    end do 
end if 
return 
end subroutine build_RBF_influence




!Smooth gradient at vertex function ===========================
function smooth_surface_gradient(gradient_surf_unsmooth,smoothRBF,surfpoints,Nsurfpoint) result(sgrad_smooth)
implicit none 

!Result
real(dp) :: sgrad_smooth(3)

!Variables - Import
integer(in) :: Nsurfpoint
integer(in), dimension(:) :: surfpoints
real(dp), dimension(:) :: smoothRBF
real(dp), dimension(:,:) :: gradient_surf_unsmooth

!Variables - Local
integer(in) :: ii 
real(dp) :: weight_sum

!Find smoothed gradient at this vertex
weight_sum = 0.0d0 
sgrad_smooth(:) = 0.0d0 
do ii=1,Nsurfpoint
    weight_sum = weight_sum + smoothRBF(ii)
    sgrad_smooth(:) = sgrad_smooth(:) + gradient_surf_unsmooth(surfpoints(ii),:)*smoothRBF(ii)
end do 
sgrad_smooth(:) = sgrad_smooth(:)/weight_sum
return 
end function smooth_surface_gradient




!Project gradient to vertex function ===========================
function project_vtx_vgrad_2_surfgrad(gradient_vol,surf2volRBF,Ri,volpoints,Nvolp) result(vsgrad)
implicit none 

!Result
real(dp) :: vsgrad(3)

!Variables - Import
integer(in) :: Nvolp
integer(in), dimension(:) :: volpoints
real(dp), dimension(:) :: surf2volRBF
real(dp), dimension(:,:) :: Ri,gradient_vol

!Variables - Local
integer(in) :: ii 
real(dp) :: gamma_xyz(Nvolp,3),gradient_vol_loc(Nvolp,3)

!Populate local gradient vector
do ii=1,Nvolp
    gradient_vol_loc(ii,:) = gradient_vol(volpoints(ii),:)
end do 

!Solve gamma coefficients 
gamma_xyz(:,1) = matmul(Ri(1:Nvolp,1:Nvolp),gradient_vol_loc(:,1))
gamma_xyz(:,2) = matmul(Ri(1:Nvolp,1:Nvolp),gradient_vol_loc(:,2))
gamma_xyz(:,3) = matmul(Ri(1:Nvolp,1:Nvolp),gradient_vol_loc(:,3))

!Construct interpolated gradient 
vsgrad(1) = dot_product(gamma_xyz(:,1),surf2volRBF(:))
vsgrad(2) = dot_product(gamma_xyz(:,2),surf2volRBF(:))
vsgrad(3) = dot_product(gamma_xyz(:,3),surf2volRBF(:))
return 
end function project_vtx_vgrad_2_surfgrad


end module cellmesh3d_gradient_coupling_mod