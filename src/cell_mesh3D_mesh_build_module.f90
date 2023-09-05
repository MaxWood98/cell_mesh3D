!cell_mesh3d mesh building module
!Max Wood - mw16116@bristol.ac.uk
!Univeristy of Bristol - Department of Aerospace Engineering

!Version 4.0
!Updated 03-08-2023

!Module
module cellmesh3d_mesh_build_mod
use cellmesh3d_adtree_mod
use cell_mesh3d_octree_mod
use cellmesh3d_geometry_mod
contains 


!Volume mesh construction subroutine ('variable' surface format) ===========================
subroutine construct_mesh(volume_mesh,cell_keep,vtx_external,vtx_type1,ot_mesh,surface_mesh,surface_adtree,&
                          cm3dopt,cm3dfailure)
implicit none 

!Variables - Import
integer(in) :: cm3dfailure
integer(in), dimension(:) :: cell_keep,vtx_external,vtx_type1
type(cm3d_options) :: cm3dopt
type(octree_data) :: ot_mesh
type(surface_data) :: surface_mesh
type(vol_mesh_data) :: volume_mesh
type(tree_data) :: surface_adtree

!Variables - Local 
integer(in) :: type_tgt!,ii
type(vol_mesh_data) :: volume_mesh_full

!Set target edge/face type for mesh 
if (cm3dopt%meshinout == 'in') then !mesh internal 
    type_tgt = 3
elseif (cm3dopt%meshinout == 'out') then !mesh external 
    type_tgt = 2
end if 

!Display
if (cm3dopt%dispt == 1) then
    write(*,'(A)') '    {constructing base complete volume mesh}'
end if 

!Construct initial full mesh 
call build_full_mesh(volume_mesh_full,ot_mesh,cell_keep,cm3dopt)
call build_vmesh_edges(volume_mesh_full%nedge,volume_mesh_full%edges,volume_mesh_full,6_in)

!Debug export full mesh =======================
! open(11,file=cm3dopt%iopath//'mesh_test_vtx')
!     do ii=1,volume_mesh_full%nvtx
!         write(11,*) volume_mesh_full%vtx(ii,:)
!     end do 
! close(11)
! open(11,file=cm3dopt%iopath//'mesh_test_faces')
!     write(11,*) maxval(volume_mesh_full%faces(:)%nvtx)
!     do ii=1,volume_mesh_full%nface
!         write(11,*) volume_mesh_full%faces(ii)%nvtx, volume_mesh_full%faces(ii)%vertices(:)
!     end do 
! close(11)
! open(11,file=cm3dopt%iopath//'cell_lr')
!     do ii=1,volume_mesh_full%nface
!         write(11,*) volume_mesh_full%faces(ii)%cleft,volume_mesh_full%faces(ii)%cright
!     end do 
! close(11)
!Debug export full mesh =======================

!Intersect the full volume mesh with the surface geometry to construct the final mesh 
call clip_vmesh2surface(volume_mesh,volume_mesh_full,surface_adtree,surface_mesh,ot_mesh,cm3dopt,cell_keep,vtx_external,&
                        vtx_type1,type_tgt,cm3dfailure)

!Debug ==============================================
! !Debug export full mesh ---
! open(11,file=cm3dopt%iopath//'mesh_test_vtx')
!     do ii=1,volume_mesh_full%nvtx
!         write(11,*) volume_mesh_full%vtx(ii,:)
!     end do 
! close(11)
! open(11,file=cm3dopt%iopath//'mesh_test_faces')
!     write(11,*) maxval(volume_mesh_full%faces(:)%nvtx)
!     do ii=1,volume_mesh_full%nface
!         write(11,*) volume_mesh_full%faces(ii)%nvtx, volume_mesh_full%faces(ii)%vertices(:)
!     end do 
! close(11)

! open(11,file=cm3dopt%iopath//'cell_lr')
!     do ii=1,volume_mesh_full%nface
!         write(11,*) volume_mesh_full%faces(ii)%cleft,volume_mesh_full%faces(ii)%cright
!     end do 
! close(11)

! !Debug export final mesh ---
! open(11,file=cm3dopt%iopath//'mesh_test_vtx')
!     do ii=1,volume_mesh%nvtx
!         write(11,*) volume_mesh%vtx(ii,:)
!     end do 
! close(11)
! open(11,file=cm3dopt%iopath//'mesh_test_faces')
!     write(11,*) maxval(volume_mesh%faces(:)%nvtx)
!     do ii=1,volume_mesh%nface
!         write(11,*) volume_mesh%faces(ii)%nvtx, volume_mesh%faces(ii)%vertices(:)
!     end do 
! close(11)
! open(11,file=cm3dopt%iopath//'cell_lr')
!     do ii=1,volume_mesh%nface
!         write(11,*) volume_mesh%faces(ii)%cleft,volume_mesh%faces(ii)%cright
!     end do 
! close(11)
! open(11,file=cm3dopt%iopath//'mesh_FT_vtx')
!     do ii=1,volume_mesh%nvtx
!         write(11,*) volume_mesh%vtx(ii,:)
!     end do 
! close(11)
return 
end subroutine construct_mesh




!Clip volume mesh to geometry subroutine ===========================
subroutine clip_vmesh2surface(volume_mesh,volume_mesh_full,surface_adtree,surface_mesh,ot_mesh,cm3dopt,cell_keep,&
                              vtx_external,vtx_type1,type_tgt,cm3dfailure)
implicit none 

!Variables - Import
integer(in) :: type_tgt,cm3dfailure
integer(in), dimension(:) :: cell_keep,vtx_external,vtx_type1
type(cm3d_options) :: cm3dopt
type(octree_data) :: ot_mesh
type(surface_data) :: surface_mesh
type(vol_mesh_data) :: volume_mesh_full,volume_mesh
type(tree_data) :: surface_adtree

!Variables - Local 
character(len=2) :: vtx_tloc
integer(in) :: ii,ee,ff,vv,aa 
integer(in) :: vtx_idx_smesh,ttgt,v1,v2,etgt,vtgt,ftgt,exist,exist2
real(dp) :: enorm,ef
real(dp) :: ve1(3),ve2(3),vbc(3)
type(edgeint), dimension(:), allocatable :: edge_clip_vm,edge_clip_sm
type(triint), dimension(:), allocatable :: tri_clip_sm
type(vtx_intersect), dimension(:), allocatable :: vtx_clip_sm

!Build surface mesh edges
call build_smesh_edges(surface_mesh,surface_mesh%maxValence) 

!Build surface mesh V2E 
allocate(surface_mesh%V2E(surface_mesh%nvtx,surface_mesh%maxValence)) 
surface_mesh%V2E(:,:) = 0 
do ee=1,surface_mesh%nedge
    v1 = surface_mesh%edges(ee,1)
    v2 = surface_mesh%edges(ee,2)
    do ii=1,surface_mesh%maxValence
        if (surface_mesh%V2E(v1,ii) == 0) then 
            surface_mesh%V2E(v1,ii) = ee 
            exit 
        end if
    end do 
    do ii=1,surface_mesh%maxValence
        if (surface_mesh%V2E(v2,ii) == 0) then 
            surface_mesh%V2E(v2,ii) = ee 
            exit 
        end if
    end do 
end do 

!Build surface mesh E2F
allocate(surface_mesh%E2F(surface_mesh%nedge,2))
surface_mesh%E2F(:,:) = 0
do ff=1,surface_mesh%nfcs
    do ee=1,3
        etgt = surface_mesh%F2E(ff,ee)
        if ((surface_mesh%E2F(etgt,1) .NE. ff) .AND. (surface_mesh%E2F(etgt,2) .NE. ff)) then 
            if (surface_mesh%E2F(etgt,1) == 0) then 
                surface_mesh%E2F(etgt,1) = ff
            elseif (surface_mesh%E2F(etgt,2) == 0) then 
                surface_mesh%E2F(etgt,2) = ff
            end if 
        end if 
    end do 
end do 

!Build volume mesh V2E 
allocate(volume_mesh_full%V2E(volume_mesh_full%nvtx,6)) 
volume_mesh_full%V2E(:,:) = 0 
do ee=1,volume_mesh_full%nedge
    v1 = volume_mesh_full%edges(ee,1)
    v2 = volume_mesh_full%edges(ee,2)
    do ii=1,6
        if (volume_mesh_full%V2E(v1,ii) == 0) then 
            volume_mesh_full%V2E(v1,ii) = ee 
            exit 
        end if
    end do 
    do ii=1,6
        if (volume_mesh_full%V2E(v2,ii) == 0) then 
            volume_mesh_full%V2E(v2,ii) = ee 
            exit 
        end if
    end do 
end do 

!Build E2F for the volume mesh 
allocate(volume_mesh_full%E2F(volume_mesh_full%nedge,4))
volume_mesh_full%E2F(:,:) = 0 
do ff=1,volume_mesh_full%nface
    do ee=1,volume_mesh_full%faces(ff)%nvtx
        etgt = volume_mesh_full%faces(ff)%edges(ee)
        do ii=1,4
            if (volume_mesh_full%E2F(etgt,ii) == 0) then 
                volume_mesh_full%E2F(etgt,ii) = ff
                exit 
            end if
        end do 
    end do 
end do 

!Initialise clipped surface structure 
cm3dopt%max_int_size = cm3dopt%NintEmax 
allocate(vtx_clip_sm(surface_mesh%nvtx))
do vv=1,surface_mesh%nvtx
    vtx_clip_sm(vv)%nintersect = 0 
    vtx_clip_sm(vv)%face_int_idx(:) = 0 
end do 
allocate(edge_clip_sm(surface_mesh%nedge))
do ee=1,surface_mesh%nedge
    edge_clip_sm(ee)%nitem = cm3dopt%NintEmax
    edge_clip_sm(ee)%nint = 0 
    edge_clip_sm(ee)%type = 0 
    allocate(edge_clip_sm(ee)%intloc(cm3dopt%NintEmax,3))
    allocate(edge_clip_sm(ee)%vtx_idx(cm3dopt%NintEmax))
    allocate(edge_clip_sm(ee)%vmface(cm3dopt%NintEmax,8))
    allocate(edge_clip_sm(ee)%nfint(cm3dopt%NintEmax))
    edge_clip_sm(ee)%vtx_idx(:) = 0 
    edge_clip_sm(ee)%vmface(:,:) = 0 
    edge_clip_sm(ee)%nfint(:) = 0 
end do 
allocate(tri_clip_sm(surface_mesh%nfcs))
do ff=1,surface_mesh%nfcs
    tri_clip_sm(ff)%nitem = cm3dopt%NintEmax
    tri_clip_sm(ff)%nvtx = 0 
    tri_clip_sm(ff)%nedge = 0 
    allocate(tri_clip_sm(ff)%edges(cm3dopt%NintEmax,3))
    allocate(tri_clip_sm(ff)%vtx_idx(cm3dopt%NintEmax))
    allocate(tri_clip_sm(ff)%nfint(cm3dopt%NintEmax))
    allocate(tri_clip_sm(ff)%face_int_idx(cm3dopt%NintEmax,8))
    allocate(tri_clip_sm(ff)%edge_int_idx(cm3dopt%NintEmax))
    allocate(tri_clip_sm(ff)%intloc(cm3dopt%NintEmax,3))
    tri_clip_sm(ff)%vtx_idx(:) = 0 
    tri_clip_sm(ff)%nfint(:) = 0 
    tri_clip_sm(ff)%face_int_idx(:,:) = 0 
    tri_clip_sm(ff)%edge_int_idx(:) = 0 
    tri_clip_sm(ff)%edges(:,:) = 0 
end do 

!Set clipped surface mesh vertex index
vtx_idx_smesh = surface_mesh%nvtx

!Display
if (cm3dopt%dispt == 1) then
    write(*,'(A)') '    {intersecting volume mesh with surface geometry}'
end if 

!Intersect volume mesh edges with the surface mesh 
call clip_voledges2surffaces(edge_clip_vm,volume_mesh_full%nedge,volume_mesh_full%edges,vtx_external,vtx_type1,&
                             volume_mesh_full,surface_adtree,surface_mesh,cm3dopt,vtx_idx_smesh,cm3dfailure)

!Display
if (cm3dopt%dispt == 1) then
    write(*,'(A)') '    {intersecting surface mesh with volume mesh}'
end if 

!Clip surface mesh edges to volume mesh faces 
call clip_surfedges2volfaces(edge_clip_sm,cell_keep,volume_mesh_full,surface_mesh,surface_adtree,ot_mesh,&
                             cm3dopt,vtx_idx_smesh,vtx_clip_sm)

!Accumulate volume mesh surface intersect vertices to the clipped triangle surface mesh structure 
do ee=1,volume_mesh_full%nedge
    if (edge_clip_vm(ee)%nint .NE. 0) then 
        do ii=1,edge_clip_vm(ee)%nint

            !Target surface triangle 
            ttgt = edge_clip_vm(ee)%surfseg(ii)

            !Classify this intersections position in triangle ttgt using its barycentric intersection 
            vtx_tloc = vtx_bary_tri_loc(edge_clip_vm(ee)%intlocBC(ii,1),edge_clip_vm(ee)%intlocBC(ii,2),&
            edge_clip_vm(ee)%intlocBC(ii,3),cm3dopt%intcointol)

            !Add to correct structure
            if (vtx_tloc == 'in') then !within triangle

                !Set intersect
                if (tri_clip_sm(ttgt)%nvtx+1 .GT. tri_clip_sm(ttgt)%nitem) then !append
                    call append_tri_clip_entry(tri_clip_sm,tri_clip_sm(ttgt)%nitem*2,ttgt)
                    if (tri_clip_sm(ttgt)%nitem .GT. cm3dopt%max_int_size) then 
                        cm3dopt%max_int_size = tri_clip_sm(ttgt)%nitem
                    end if 
                end if
                tri_clip_sm(ttgt)%nvtx = tri_clip_sm(ttgt)%nvtx + 1
                tri_clip_sm(ttgt)%vtx_idx(tri_clip_sm(ttgt)%nvtx) = edge_clip_vm(ee)%vtx_idx(ii)
                tri_clip_sm(ttgt)%edge_int_idx(tri_clip_sm(ttgt)%nvtx) = ee 
                tri_clip_sm(ttgt)%intloc(tri_clip_sm(ttgt)%nvtx,:) = edge_clip_vm(ee)%intloc(ii,:)

                !Add all faces on this edge to the intersect
                do aa=1,4
                    ftgt = volume_mesh_full%E2F(ee,aa)
                    if (ftgt .GT. 0) then 
                        tri_clip_sm(ttgt)%nfint(tri_clip_sm(ttgt)%nvtx) = tri_clip_sm(ttgt)%nfint(tri_clip_sm(ttgt)%nvtx) + 1
                        tri_clip_sm(ttgt)%face_int_idx(tri_clip_sm(ttgt)%nvtx,tri_clip_sm(ttgt)%nfint(tri_clip_sm(ttgt)%nvtx)) &
                        = ftgt
                    end if 
                end do 
            elseif (vtx_tloc(1:1) == 'e') then !on an edge 

                !Identify edge
                if (vtx_tloc == 'e1') then !on edge 1
                    etgt = surface_mesh%F2E(ttgt,1)
                elseif (vtx_tloc == 'e2') then !on edge 2
                    etgt = surface_mesh%F2E(ttgt,2)
                elseif (vtx_tloc == 'e3') then !on edge 3
                    etgt = surface_mesh%F2E(ttgt,3)
                end if 

                !Edge ends
                ve1(:) = surface_mesh%vertices(surface_mesh%edges(etgt,1),:)
                ve2(:) = surface_mesh%vertices(surface_mesh%edges(etgt,2),:)

                !Edge norm 
                enorm = norm2(ve2(:) - ve1(:))

                !Check if this intersection already exists to within tollerance 
                exist = 0 
                do aa=1,edge_clip_sm(etgt)%nint
                    if (norm2(edge_clip_vm(ee)%intloc(ii,:) - edge_clip_sm(etgt)%intloc(aa,:)) &
                    .LE. set_cmp_prc(enorm*cm3dopt%intcointol,min_precision)) then 
                        exist = aa
                        exit 
                    end if
                end do 

                !Cases on existance 
                if (exist == 0) then !new intersection on this edge

                    !Set properties 
                    edge_clip_sm(etgt)%type = 4
                    edge_clip_sm(etgt)%nint = edge_clip_sm(etgt)%nint + 1
                    edge_clip_sm(etgt)%intloc(edge_clip_sm(etgt)%nint,:) = edge_clip_vm(ee)%intloc(ii,:)
                    edge_clip_sm(etgt)%vtx_idx(edge_clip_sm(etgt)%nint) = edge_clip_vm(ee)%vtx_idx(ii)

                    !Add intersection for each vmesh face on this volume mesh edge 
                    do aa=1,4
                        ftgt = volume_mesh_full%E2F(ee,aa)
                        if (ftgt .GT. 0) then 
                            edge_clip_sm(etgt)%nfint(edge_clip_sm(etgt)%nint) = edge_clip_sm(etgt)%nfint(edge_clip_sm(etgt)%nint) &
                            + 1
                            edge_clip_sm(etgt)%vmface(edge_clip_sm(etgt)%nint,edge_clip_sm(etgt)%nfint(edge_clip_sm(etgt)%nint)) &
                            = ftgt 
                        end if 
                    end do 
                else !existing intersection so add all new volume mesh faces from this edge to this intersection and link the volume mesh intersect index

                    !Set intersect index on the volume mesh edge as the surface vertex 
                    edge_clip_vm(ee)%vtx_idx(ii) = edge_clip_sm(etgt)%vtx_idx(exist) !set as positive as this points to the same vertex set as the surface mesh vertices 

                    !Add all new faces on this volume mesh edge to the intersection
                    do aa=1,4
                        ftgt = volume_mesh_full%E2F(ee,aa)
                        if (ftgt .GT. 0) then 
                            exist2 = 0 
                            do ff=1,edge_clip_sm(etgt)%nfint(exist)
                                if (edge_clip_sm(etgt)%vmface(exist,ff) == ftgt) then 
                                    exist2 = 1
                                    exit 
                                end if 
                            end do 
                            if (exist2 == 0) then 
                                edge_clip_sm(etgt)%nfint(exist) = edge_clip_sm(etgt)%nfint(exist) + 1
                                edge_clip_sm(etgt)%vmface(exist,edge_clip_sm(etgt)%nfint(exist)) = ftgt 
                            end if 
                        end if 
                    end do 
                end if 
            elseif (vtx_tloc(1:1) == 'v') then !on a vertex
                
                !Identify vertex
                if (vtx_tloc == 'v1') then !on vertex 1
                    vtgt = surface_mesh%connectivity(ttgt,1)
                elseif (vtx_tloc == 'v2') then !on vertex 2
                    vtgt = surface_mesh%connectivity(ttgt,2)
                elseif (vtx_tloc == 'v3') then !on vertex 3
                    vtgt = surface_mesh%connectivity(ttgt,3)
                end if 

                !Set intersect index on the volume mesh edge as the surface vertex
                edge_clip_vm(ee)%vtx_idx(ii) = vtgt !set as positive as this points to the same vertex set as the surface mesh vertices 

                !Add each face on this vertex
                do aa=1,4

                    !Target face
                    ftgt = volume_mesh_full%E2F(ee,aa)

                    !Check if this face already intersects this vertex
                    exist = 0 
                    do ff=1,vtx_clip_sm(vtgt)%nintersect
                        if (vtx_clip_sm(vtgt)%face_int_idx(ff) == ftgt) then 
                            exist = 1
                            exit 
                        end if 
                    end do 

                    !Add if new
                    if (exist == 0) then 
                        vtx_clip_sm(vtgt)%nintersect = vtx_clip_sm(vtgt)%nintersect + 1
                        vtx_clip_sm(vtgt)%face_int_idx(vtx_clip_sm(vtgt)%nintersect) = ftgt 
                    end if 
                end do 
            end if 
        end do 
    end if
end do 

!Construct edge mesh on any intersected surface mesh edges 
call build_surf_intersect_edge_mesh(edge_clip_sm,surface_mesh,cm3dopt)

!Debug write intersections =======================
! open(11,file='io/vmsint_edge.dat') 
! do ee=1,surface_mesh%nedge
!     if (edge_clip_sm(ee)%nint .GT. 0) then 
!         do ii=1,edge_clip_sm(ee)%nint
!             write(11,*) edge_clip_sm(ee)%intloc(ii,:),edge_clip_sm(ee)%vmface(ii,:)
!         end do 
!     end if 
! end do 
! close(11)
! open(11,file='io/vmsint_face.dat') 
! do ff=1,surface_mesh%nfcs
!     if (tri_clip_sm(ff)%nvtx .GT. 0) then 
!         do ii=1,tri_clip_sm(ff)%nvtx
!             write(11,*) tri_clip_sm(ff)%intloc(ii,:),volume_mesh_full%E2F(tri_clip_sm(ff)%edge_int_idx(ii),:)
!         end do 
!     end if 
! end do 
! close(11)   
! open(11,file='io/vmsint_vtx.dat') 
! do vv=1,surface_mesh%nvtx
!     if (vtx_clip_sm(vv)%nintersect .GT. 0) then 
!         write(11,*) surface_mesh%vertices(vv,:)
!     end if 
! end do 
! close(11)
!Debug write intersections =======================

!Build complete list of surface vertices including all clipping intersections 
surface_mesh%nvtxf = vtx_idx_smesh
allocate(surface_mesh%vertices_full(vtx_idx_smesh,3))
allocate(surface_mesh%vtx_rcurv_full(vtx_idx_smesh))
surface_mesh%vertices_full(1:surface_mesh%nvtx,:) = surface_mesh%vertices(:,:)
surface_mesh%vtx_rcurv_full(1:surface_mesh%nvtx) = surface_mesh%vtx_rcurv(:)
do ee=1,surface_mesh%nedge
    if (edge_clip_sm(ee)%nint .GT. 0) then 
        do ii=1,edge_clip_sm(ee)%nint

            !Assign intersection location
            surface_mesh%vertices_full(edge_clip_sm(ee)%vtx_idx(ii),:) = edge_clip_sm(ee)%intloc(ii,:)

            !Edge fraction
            ef = norm2(edge_clip_sm(ee)%intloc(ii,:) - surface_mesh%vertices(surface_mesh%edges(ee,1),:))/&
            norm2(surface_mesh%vertices(surface_mesh%edges(ee,2),:) - surface_mesh%vertices(surface_mesh%edges(ee,1),:))

            !Assign interpolated surface curvature 
            surface_mesh%vtx_rcurv_full(edge_clip_sm(ee)%vtx_idx(ii)) = &
            (1.0d0 - ef)*surface_mesh%vtx_rcurv(surface_mesh%edges(ee,1)) + &
            ef*surface_mesh%vtx_rcurv(surface_mesh%edges(ee,2))
        end do
    end if
end do 
do ff=1,surface_mesh%nfcs
    if (tri_clip_sm(ff)%nvtx .GT. 0) then 
        do ii=1,tri_clip_sm(ff)%nvtx

            !Assign intersection location
            surface_mesh%vertices_full(tri_clip_sm(ff)%vtx_idx(ii),:) = tri_clip_sm(ff)%intloc(ii,:)

            !Barycentric location 
            vbc = cart2baryc_onplane(tri_clip_sm(ff)%intloc(ii,:),surface_mesh%vertices(surface_mesh%connectivity(ff,1),:),&
            surface_mesh%vertices(surface_mesh%connectivity(ff,2),:),surface_mesh%vertices(surface_mesh%connectivity(ff,3),:))

            !Assign interpolated surface curvature
            surface_mesh%vtx_rcurv_full(tri_clip_sm(ff)%vtx_idx(ii)) = &
            surface_mesh%vtx_rcurv(surface_mesh%connectivity(ff,1))*vbc(1) +&
            surface_mesh%vtx_rcurv(surface_mesh%connectivity(ff,2))*vbc(2) +&
            surface_mesh%vtx_rcurv(surface_mesh%connectivity(ff,3))*vbc(3)
        end do 
    end if 
end do 

!Display
if (cm3dopt%dispt == 1) then
    write(*,'(A)') '    {constructing volume mesh clipped surface mesh edges}'
end if 

!Build surface triangle internal edges
call build_surftri_internal_edges(tri_clip_sm,edge_clip_sm,vtx_clip_sm,surface_mesh,volume_mesh_full,cm3dopt)

!Display
if (cm3dopt%dispt == 1) then
    write(*,'(A)') '    {constructing surface mesh clipped clipped volume mesh faces}'
end if 

!Build surface mesh clippled volume mesh faces 
call build_surface_clipped_vmesh_faces(volume_mesh_full,surface_mesh,surface_adtree,cm3dopt,tri_clip_sm,vtx_clip_sm,&
                                       edge_clip_sm,edge_clip_vm,cell_keep,vtx_idx_smesh,type_tgt,surface_mesh%vertices_full)

!Display
if (cm3dopt%dispt == 1) then
    write(*,'(A)') '    {constructing volume mesh clipped surface mesh faces}'
end if 

!Build volume mesh clipped surface triangles 
call build_volume_clipped_smesh_faces(surface_mesh,tri_clip_sm,edge_clip_sm,vtx_clip_sm,volume_mesh_full,ot_mesh,cm3dopt)

!Display
if (cm3dopt%dispt == 1) then
    write(*,'(A)') '    {constructing clipped volume mesh}'
end if 

!Construct complete final volume mesh 
call build_final_vmesh(volume_mesh,volume_mesh_full,surface_mesh,type_tgt)

!Set cell associations of all unclipped surface mesh faces
call set_surface_mesh_cells(volume_mesh)
return 
end subroutine clip_vmesh2surface




!Set surface mesh cells subroutine ===========================
subroutine set_surface_mesh_cells(volume_mesh)
implicit none 

!Variables - Import
type(vol_mesh_data) :: volume_mesh

!Variables - Local
integer(in) :: ff,ee,vv,fiter
integer(in) :: Nface_surf,e1,e2,v1,v2,vlceUB,ftgt,evalid,Nedge,eidx,nupdate,f1,f2
integer(in) :: face_surface_idx(volume_mesh%nface),valence(volume_mesh%nvtx)
integer(in), dimension(:), allocatable :: edge_invalid 
integer(in), dimension(:,:), allocatable :: vconnect,edge_idx,E2F

!Extract the surface mesh within the volume mesh 
Nface_surf = 0 
face_surface_idx(:) = 0 
do ff=1,volume_mesh%nface
    if ((volume_mesh%faces(ff)%cleft == -1) .OR. (volume_mesh%faces(ff)%cright == -1)) then 
        Nface_surf = Nface_surf + 1
        face_surface_idx(Nface_surf) = ff 
    end if 
end do 

!Find upper bound of valence
valence(:) = 0 
do ff=1,Nface_surf
    ftgt = face_surface_idx(ff)
    do ee=1,volume_mesh%faces(ftgt)%nvtx
        e1 = ee 
        e2 = mod(ee,volume_mesh%faces(ftgt)%nvtx) + 1
        v1 = volume_mesh%faces(ftgt)%vertices(e1)
        v2 = volume_mesh%faces(ftgt)%vertices(e2)
        valence(v1) = valence(v1) + 1
        valence(v2) = valence(v2) + 1
    end do 
end do 
vlceUB = 2*maxval(valence(:))

!Construct actual valence of each vertex and build edges
Nedge = 0 
allocate(vconnect(volume_mesh%nvtx,vlceUB))
allocate(edge_idx(volume_mesh%nvtx,vlceUB)) 
valence(:) = 0
vconnect(:,:) = 0 
edge_idx(:,:) = 0 
do ff=1,Nface_surf
    ftgt = face_surface_idx(ff)
    do ee=1,volume_mesh%faces(ftgt)%nvtx

        !Ends of this edge
        e1 = ee 
        e2 = mod(ee,volume_mesh%faces(ftgt)%nvtx) + 1
        v1 = volume_mesh%faces(ftgt)%vertices(e1)
        v2 = volume_mesh%faces(ftgt)%vertices(e2)

        !Check against vconnect 
        eidx = 0 
        evalid = 1
        do vv=1,vlceUB
            if (vconnect(v1,vv) == v2) then
                evalid = 0
                eidx = edge_idx(v1,vv)
                exit 
            end if 
            if (vconnect(v2,vv) == v1) then 
                evalid = 0
                eidx = edge_idx(v2,vv)
                exit 
            end if 
        end do 

        !Add valence if new edge
        if (evalid == 1) then 
            
            !Increment edge count 
            Nedge = Nedge + 1

            !Add edge 
            volume_mesh%faces(ftgt)%edges(ee) = Nedge

            !Increment valence on each vertex
            valence(v1) = valence(v1) + 1
            valence(v2) = valence(v2) + 1

            !Update vconnect
            do vv=1,vlceUB
                if (vconnect(v1,vv) == 0) then 
                    vconnect(v1,vv) = v2
                    edge_idx(v1,vv) = Nedge
                    exit
                end if 
            end do 
            do vv=1,vlceUB
                if (vconnect(v2,vv) == 0) then 
                    vconnect(v2,vv) = v1
                    edge_idx(v2,vv) = Nedge
                    exit 
                end if 
            end do 
        else
            volume_mesh%faces(ftgt)%edges(ee) = eidx
        end if 
    end do 
end do 

!Build E2F on the surface mesh 
allocate(E2F(Nedge,2))
E2F(:,:) = 0 
do ff=1,Nface_surf
    ftgt = face_surface_idx(ff)
    do ee=1,volume_mesh%faces(ftgt)%nvtx
        eidx = volume_mesh%faces(ftgt)%edges(ee)
        if ((E2F(eidx,1) .NE. ftgt) .AND. (E2F(eidx,2) .NE. ftgt)) then 
            if (E2F(eidx,1) == 0) then 
                E2F(eidx,1) = ftgt
            elseif (E2F(eidx,2) == 0) then 
                E2F(eidx,2) = ftgt
            end if 
        end if 
    end do 
end do 

!Tag all edges that are invalid to flood through as they form cell boundaries 
allocate(edge_invalid(Nedge))
edge_invalid(:) = 0 
do ff=1,volume_mesh%nface
    if ((volume_mesh%faces(ff)%cleft .NE. -1) .AND. (volume_mesh%faces(ff)%cright .NE. -1)) then 
        do ee=1,volume_mesh%faces(ff)%nvtx

            !Ends of this edge in this face
            e1 = ee 
            e2 = mod(ee,volume_mesh%faces(ff)%nvtx) + 1
            v1 = volume_mesh%faces(ff)%vertices(e1)
            v2 = volume_mesh%faces(ff)%vertices(e2)

            !Check if this vertex pair forms an edge in the surface mesh 
            eidx = 0 
            do vv=1,vlceUB
                if (vconnect(v1,vv) == v2) then
                    eidx = edge_idx(v1,vv)
                    exit 
                end if 
                if (vconnect(v2,vv) == v1) then 
                    eidx = edge_idx(v2,vv)
                    exit 
                end if 
            end do 

            !If they form an edge then tag this as invalid 
            if (eidx .NE. 0) then 
                edge_invalid(eidx) = 1
            end if 
        end do 
    end if 
end do 

!Flood cell associations from clipped faces 
do fiter=1,volume_mesh%nface
    nupdate = 0 
    do ee=1,Nedge
        if (edge_invalid(ee) == 0) then 
            f1 = E2F(ee,1)
            f2 = E2F(ee,2)
            if ((volume_mesh%faces(f1)%cright == 0) .AND. (volume_mesh%faces(f2)%cright .NE. 0)) then 
                volume_mesh%faces(f1)%cright = volume_mesh%faces(f2)%cright
                nupdate = nupdate + 1
            end if 
            if ((volume_mesh%faces(f2)%cright == 0) .AND. (volume_mesh%faces(f1)%cright .NE. 0)) then 
                volume_mesh%faces(f2)%cright = volume_mesh%faces(f1)%cright
                nupdate = nupdate + 1
            end if 
        end if 
    end do 
    if (nupdate == 0) then 
        exit 
    end if 
end do 

!Check for faces with no cell association 
do ff=1,volume_mesh%nface
    if ((volume_mesh%faces(ff)%cright == 0) .OR. (volume_mesh%faces(ff)%cleft == 0)) then 
        print *, '** surface face with no assigned cell ',ff 
    end if 
end do 
return 
end subroutine set_surface_mesh_cells




!Build final clipped volume mesh subroutine ===========================
subroutine build_final_vmesh(volume_mesh,volume_mesh_full,surface_mesh,type_tgt)
implicit none 

!Variables - Import
integer(in) :: type_tgt
type(surface_data) :: surface_mesh
type(vol_mesh_data) :: volume_mesh_full,volume_mesh

!Variables - Local
integer(in) :: vv,ff,aa 
integer(in) :: vins,face_idx,fcount,vtgt
integer(in), dimension(:), allocatable :: vtxmap

!Build complete vertex list 
volume_mesh%nvtx = volume_mesh_full%nvtx + surface_mesh%nvtxf
vins = volume_mesh_full%nvtx 
allocate(vtxmap(surface_mesh%nvtxf))
vtxmap(:) = 0 
do vv=1,surface_mesh%nvtxf
    vins = vins + 1
    vtxmap(vv) = vins 
end do 
allocate(volume_mesh%vtx(volume_mesh%nvtx,3)) 
allocate(volume_mesh%vtx_sRcurv(volume_mesh%nvtx))
allocate(volume_mesh%vmvtx_2_smvtx(volume_mesh%nvtx))
volume_mesh%vtx(1:volume_mesh_full%nvtx,:) = volume_mesh_full%vtx(:,:)
volume_mesh%vtx_sRcurv(:) = 0.0d0 
do vv=1,surface_mesh%nvtxf
    volume_mesh%vtx(vtxmap(vv),:) = surface_mesh%vertices_full(vv,:)
    volume_mesh%vtx_sRcurv(vtxmap(vv)) = surface_mesh%vtx_rcurv_full(vv)
end do
volume_mesh%vmvtx_2_smvtx(:) = 0 
do vv=1,surface_mesh%nvtx
    volume_mesh%vmvtx_2_smvtx(vtxmap(vv)) = vv
end do 

!Index complete face list
face_idx = 0 
do ff=1,volume_mesh_full%nface
    face_idx = face_idx + volume_mesh_full%faces(ff)%nfclipped
end do 
do ff=1,surface_mesh%nfcs
    face_idx = face_idx + surface_mesh%tri_clipped(ff)%nfclipped
end do 

!Allocate final mesh face array
volume_mesh%nface = face_idx
allocate(volume_mesh%faces(face_idx))

!Add clipped mesh faces
fcount = 0 
face_idx = 0 
do ff=1,volume_mesh_full%nface !volume mesh
    do aa=1,volume_mesh_full%faces(ff)%nfclipped

        !Initialise face
        fcount = fcount + 1
        face_idx = face_idx + 1
        volume_mesh%faces(face_idx)%nvtx = volume_mesh_full%faces(ff)%face_clipped(aa)%nvtx
        allocate(volume_mesh%faces(face_idx)%vertices(volume_mesh%faces(face_idx)%nvtx))
        allocate(volume_mesh%faces(face_idx)%edges(volume_mesh%faces(face_idx)%nvtx))
        volume_mesh%faces(face_idx)%edges(:) = 0 

        !Add vertices
        do vv=1,volume_mesh_full%faces(ff)%face_clipped(aa)%nvtx
            vtgt = volume_mesh_full%faces(ff)%face_clipped(aa)%vertices(vv)
            if (vtgt .LT. 0) then
                vtgt = vtxmap(abs(vtgt))
            end if  
            volume_mesh%faces(face_idx)%vertices(vv) = vtgt
        end do

        !Set adjacency
        volume_mesh%faces(face_idx)%cleft = volume_mesh_full%faces(ff)%cleft
        volume_mesh%faces(face_idx)%cright = volume_mesh_full%faces(ff)%cright
    end do 
end do 
do ff=1,surface_mesh%nfcs !surface mesh
    do aa=1,surface_mesh%tri_clipped(ff)%nfclipped

        !Initialise face
        fcount = fcount + 1
        face_idx = face_idx + 1
        volume_mesh%faces(face_idx)%nvtx = surface_mesh%tri_clipped(ff)%face_clipped(aa)%nvtx
        allocate(volume_mesh%faces(face_idx)%vertices(volume_mesh%faces(face_idx)%nvtx))
        allocate(volume_mesh%faces(face_idx)%edges(volume_mesh%faces(face_idx)%nvtx))
        volume_mesh%faces(face_idx)%edges(:) = 0 

        !Add vertices
        do vv=1,volume_mesh%faces(face_idx)%nvtx
            vtgt = surface_mesh%tri_clipped(ff)%face_clipped(aa)%vertices(vv)
            vtgt = vtxmap(vtgt) 
            volume_mesh%faces(face_idx)%vertices(vv) = vtgt
        end do

        !Set adjacency
        if (type_tgt == 2) then !mesh external 
            volume_mesh%faces(face_idx)%cleft = -1
            volume_mesh%faces(face_idx)%cright = surface_mesh%tri_clipped(ff)%face_clipped(aa)%cright
        elseif (type_tgt == 3) then !mesh internal 
            volume_mesh%faces(face_idx)%cleft = surface_mesh%tri_clipped(ff)%face_clipped(aa)%cright
            volume_mesh%faces(face_idx)%cright = -1
        end if 
        !print *, volume_mesh%faces(face_idx)%cright
        ! if (face_idx == 33978) then 
        !     print *, 'SMface = ',ff
        ! end if 
    end do 
end do 

!Set final cell count 
volume_mesh%ncell = volume_mesh_full%ncell

!Set cell levels 
allocate(volume_mesh%cell_level(volume_mesh%ncell))
volume_mesh%cell_level(:) = volume_mesh_full%cell_level(:)
return 
end subroutine build_final_vmesh




!Build surface mesh volume clipped faces ===========================
subroutine build_volume_clipped_smesh_faces(surface_mesh,tri_clip_sm,edge_clip_sm,vtx_clip_sm,volume_mesh_full,ot_mesh,cm3dopt)
implicit none 

!Variables - Import
type(surface_data) :: surface_mesh
type(vol_mesh_data) :: volume_mesh_full
type(triint), dimension(:) :: tri_clip_sm
type(edgeint), dimension(:) :: edge_clip_sm
type(vtx_intersect), dimension(:) :: vtx_clip_sm
type(cm3d_options) :: cm3dopt
type(octree_data) :: ot_mesh

!Variables - Local 
integer(in) :: ff,vv,ee,aa,ii,cc,ss 
integer(in) :: vintvalid1,vintvalid2,ebase,vbase,enew,vnew,all_te_subedge,diff_parent,nint_tri
integer(in) :: nc_ontri,cl,cr,ftgt,ctgt,etgt,vtgt,nedge_tface,ve1,ve2,coctree,v1,v2,face_nvtx,vc,vn,v0
integer(in) :: cell_on_tri(cm3dopt%max_int_size),edge_face_idx(2*cm3dopt%max_int_size),edges_on_tri(2*cm3dopt%max_int_size,4)
integer(in) :: cell_selected(volume_mesh_full%ncell),tedge_cell(3),tvtx_cell(3),eend_cell(3,2)
integer(in) :: vtx_face_pos(surface_mesh%nvtxf),fliparray(cm3dopt%max_int_size)
integer(in) :: v2e_flocal(surface_mesh%nvtxf,2)
real(dp) :: zxmin,zxmax,zymin,zymax,zzmin,zzmax
real(dp) :: emid(3),vtxref(3),Nft(3),Nf(3),vtxC(3),vtxN(3),fmid(3)

!Build sub clipped faces for each surface mesh face 
nedge_tface = 0
nc_ontri = 0 
fliparray(:) = 0 
edge_face_idx(:) = 0 
cell_on_tri(:) = 0 
vtx_face_pos(:) = 0 
cell_selected(:) = 0
v2e_flocal(:,:) = 0 
edges_on_tri(:,:) = 0 
allocate(surface_mesh%tri_clipped(surface_mesh%nfcs))
do ff=1,surface_mesh%nfcs

    !Initialise number of clipped faces
    surface_mesh%tri_clipped(ff)%nfclipped = 0 

    !Count all intersections in edges or the face of this triangle 
    nint_tri = tri_clip_sm(ff)%nedge
    do ee=1,3
        etgt = surface_mesh%F2E(ff,ee)
        nint_tri = nint_tri + edge_clip_sm(etgt)%nint
    end do 

    !Build sub-faces if this face has been clipped
    if (nint_tri .GT. 0) then !Build sub-faces as clipped by vm face

        !Normal vector to this triangle 
        Nf = newell_normal(3_in,surface_mesh%connectivity(ff,:),surface_mesh%vertices)

        !Accumulate all vm cells on this triangle 
        nc_ontri = 0 
        do vv=1,tri_clip_sm(ff)%nvtx !from intersections within this face 
            do aa=1,tri_clip_sm(ff)%nfint(vv)
                ftgt = tri_clip_sm(ff)%face_int_idx(vv,aa)
                cl = volume_mesh_full%faces(ftgt)%cleft
                cr = volume_mesh_full%faces(ftgt)%cright
                ctgt = cl 
                if (ctgt .GT. 0) then
                    if (cell_selected(ctgt) == 0) then 
                        nc_ontri = nc_ontri + 1
                        cell_on_tri(nc_ontri) = ctgt 
                        cell_selected(ctgt) = nc_ontri
                    end if 
                end if 
                ctgt = cr 
                if (ctgt .GT. 0) then
                    if (cell_selected(ctgt) == 0) then 
                        nc_ontri = nc_ontri + 1
                        cell_on_tri(nc_ontri) = ctgt 
                        cell_selected(ctgt) = nc_ontri
                    end if 
                end if 
            end do 
        end do 
        do ee=1,3 !from triangle edge intersections 
            etgt = surface_mesh%F2E(ff,ee)
            do ii=1,edge_clip_sm(etgt)%nint
                do aa=1,edge_clip_sm(etgt)%nfint(ii)
                    ftgt = edge_clip_sm(etgt)%vmface(ii,aa)
                    cl = volume_mesh_full%faces(ftgt)%cleft
                    cr = volume_mesh_full%faces(ftgt)%cright
                    ctgt = cl 
                    if (ctgt .GT. 0) then
                        if (cell_selected(ctgt) == 0) then 
                            nc_ontri = nc_ontri + 1
                            cell_on_tri(nc_ontri) = ctgt 
                            cell_selected(ctgt) = nc_ontri
                        end if 
                    end if 
                    ctgt = cr 
                    if (ctgt .GT. 0) then
                        if (cell_selected(ctgt) == 0) then 
                            nc_ontri = nc_ontri + 1
                            cell_on_tri(nc_ontri) = ctgt 
                            cell_selected(ctgt) = nc_ontri
                        end if 
                    end if 
                end do 
            end do 
        end do 
        do vv=1,3 !from vertex intersections 
            vtgt = surface_mesh%connectivity(ff,vv)
            do aa=1,vtx_clip_sm(vtgt)%nintersect
                ftgt = vtx_clip_sm(vtgt)%face_int_idx(aa)
                cl = volume_mesh_full%faces(ftgt)%cleft
                cr = volume_mesh_full%faces(ftgt)%cright
                ctgt = cl 
                if (ctgt .GT. 0) then
                    if (cell_selected(ctgt) == 0) then 
                        nc_ontri = nc_ontri + 1
                        cell_on_tri(nc_ontri) = ctgt 
                        cell_selected(ctgt) = nc_ontri
                    end if 
                end if 
                ctgt = cr 
                if (ctgt .GT. 0) then
                    if (cell_selected(ctgt) == 0) then 
                        nc_ontri = nc_ontri + 1
                        cell_on_tri(nc_ontri) = ctgt 
                        cell_selected(ctgt) = nc_ontri
                    end if 
                end if 
            end do 
        end do 
        cell_selected(cell_on_tri(1:nc_ontri)) = 0 

        !Assign a cell from cell_on_tri to any edge of this triangle that has no face intersections 
        tedge_cell(:) = 0 
        do ee=1,3
            etgt = surface_mesh%F2E(ff,ee)
            if (edge_clip_sm(etgt)%nint == 0) then !Pick cell that contains the edge midpoint 

                !Find edge midpoint 
                ve1 = surface_mesh%edges(etgt,1)
                ve2 = surface_mesh%edges(etgt,2)
                emid(1) = 0.5d0*(surface_mesh%vertices(ve1,1) + surface_mesh%vertices(ve2,1))
                emid(2) = 0.5d0*(surface_mesh%vertices(ve1,2) + surface_mesh%vertices(ve2,2))
                emid(3) = 0.5d0*(surface_mesh%vertices(ve1,3) + surface_mesh%vertices(ve2,3))

                !Test containement against each cell selected for this triangle 
                do cc=1,nc_ontri
                    ctgt = cell_on_tri(cc)
                    coctree = volume_mesh_full%cell_otidx(ctgt)
                    zxmin = ot_mesh%vtx(ot_mesh%cell_vcnr(coctree,1),1)
                    zxmax = ot_mesh%vtx(ot_mesh%cell_vcnr(coctree,7),1)
                    zymin = ot_mesh%vtx(ot_mesh%cell_vcnr(coctree,1),2)
                    zymax = ot_mesh%vtx(ot_mesh%cell_vcnr(coctree,7),2)
                    zzmin = ot_mesh%vtx(ot_mesh%cell_vcnr(coctree,1),3)
                    zzmax = ot_mesh%vtx(ot_mesh%cell_vcnr(coctree,7),3)
                    if ((emid(1) .GE. zxmin) .AND. (emid(1) .LE. zxmax)) then 
                        if ((emid(2) .GE. zymin) .AND. (emid(2) .LE. zymax)) then 
                            if ((emid(3) .GE. zzmin) .AND. (emid(3) .LE. zzmax)) then 
                                tedge_cell(ee) = ctgt 
                                exit 
                            end if 
                        end if 
                    end if 
                end do 
            end if
        end do 

        !Assign a cell from cell_on_tri to any vertex on this triangle that has no intersections 
        tvtx_cell(:) = 0 
        do vv=1,3
            vtxref(:) = surface_mesh%vertices(surface_mesh%connectivity(ff,vv),:)
            do cc=1,nc_ontri
                ctgt = cell_on_tri(cc)
                coctree = volume_mesh_full%cell_otidx(ctgt)
                zxmin = ot_mesh%vtx(ot_mesh%cell_vcnr(coctree,1),1)
                zxmax = ot_mesh%vtx(ot_mesh%cell_vcnr(coctree,7),1)
                zymin = ot_mesh%vtx(ot_mesh%cell_vcnr(coctree,1),2)
                zymax = ot_mesh%vtx(ot_mesh%cell_vcnr(coctree,7),2)
                zzmin = ot_mesh%vtx(ot_mesh%cell_vcnr(coctree,1),3)
                zzmax = ot_mesh%vtx(ot_mesh%cell_vcnr(coctree,7),3)
                if ((vtxref(1) .GE. zxmin) .AND. (vtxref(1) .LE. zxmax)) then 
                    if ((vtxref(2) .GE. zymin) .AND. (vtxref(2) .LE. zymax)) then 
                        if ((vtxref(3) .GE. zzmin) .AND. (vtxref(3) .LE. zzmax)) then 
                            tvtx_cell(vv) = ctgt 
                            exit 
                        end if 
                    end if 
                end if 
            end do 
        end do 

        !Assign edge end cell containements 
        eend_cell(:,:) = 0 
        do ee=1,3 !Triangle edges 
            etgt = surface_mesh%F2E(ff,ee)
            v1 = surface_mesh%edges(etgt,1)
            do vv=1,3
                if (v1 == surface_mesh%connectivity(ff,vv)) then 
                    eend_cell(ee,1) = tvtx_cell(vv)
                    exit 
                end if 
            end do 
            v2 = surface_mesh%edges(etgt,2)
            do vv=1,3
                if (v2 == surface_mesh%connectivity(ff,vv)) then 
                    eend_cell(ee,2) = tvtx_cell(vv)
                    exit 
                end if 
            end do 
        end do 

        !Allocate clipped face array to its maximum size
        allocate(surface_mesh%tri_clipped(ff)%face_clipped(nc_ontri))

        !Build sub-faces for each cell on this face 
        do cc=1,nc_ontri

            !Target cell 
            ctgt = cell_on_tri(cc)

            !Find all edges on this triangle that are associated with cell ctgt 
            nedge_tface = 0 
            do ee=1,tri_clip_sm(ff)%nedge !Triangle internal edges
                ftgt = tri_clip_sm(ff)%edges(ee,3)
                cl = volume_mesh_full%faces(ftgt)%cleft
                cr = volume_mesh_full%faces(ftgt)%cright
                if ((cl == ctgt) .OR. (cr == ctgt)) then 
                    nedge_tface = nedge_tface + 1 
                    edges_on_tri(nedge_tface,1:2) = tri_clip_sm(ff)%edges(ee,1:2)
                    edges_on_tri(nedge_tface,3) = 1 !tag as triangle internal edge 
                    edges_on_tri(nedge_tface,4) = ftgt !tag with parent face
                end if 
            end do 
            do ee=1,3 !Triangle edges 
                etgt = surface_mesh%F2E(ff,ee)
                if (edge_clip_sm(etgt)%nint == 0) then !Take full edge if assigned cell ctgt or if both vertices on the edge intersect the cell 
                    if (tedge_cell(ee) == ctgt) then !If edge is assigned to the current cell through containement 
                        nedge_tface = nedge_tface + 1 
                        edges_on_tri(nedge_tface,1) = surface_mesh%edges(etgt,1)
                        edges_on_tri(nedge_tface,2) = surface_mesh%edges(etgt,2)
                        edges_on_tri(nedge_tface,3) = 2 !tag as triangle full edge 
                        edges_on_tri(nedge_tface,4) = etgt !tag with parent edge
                    else !If no containement then check using end vertex associations 
                        v1 = surface_mesh%edges(etgt,1)
                        v2 = surface_mesh%edges(etgt,2)
                        vintvalid1 = 0 
                        do vv=1,vtx_clip_sm(v1)%nintersect
                            ftgt = vtx_clip_sm(v1)%face_int_idx(vv)
                            cl = volume_mesh_full%faces(ftgt)%cleft
                            cr = volume_mesh_full%faces(ftgt)%cright
                            if (cl == ctgt) then
                                vintvalid1 = 1
                                exit
                            end if 
                            if (cr == ctgt) then
                                vintvalid1 = 1
                                exit
                            end if 
                        end do 
                        vintvalid2 = 0 
                        do vv=1,vtx_clip_sm(v2)%nintersect
                            ftgt = vtx_clip_sm(v2)%face_int_idx(vv)
                            cl = volume_mesh_full%faces(ftgt)%cleft
                            cr = volume_mesh_full%faces(ftgt)%cright
                            if (cl == ctgt) then
                                vintvalid2 = 1
                                exit
                            end if 
                            if (cr == ctgt) then
                                vintvalid2 = 1
                                exit
                            end if 
                        end do 
                        if ((vintvalid1 == 1) .AND. (vintvalid2 == 1)) then !both ends of this segment intersect the face
                            nedge_tface = nedge_tface + 1
                            edges_on_tri(nedge_tface,1) = v1
                            edges_on_tri(nedge_tface,2) = v2
                            edges_on_tri(nedge_tface,3) = 2 !tag as triangle full edge
                            edges_on_tri(nedge_tface,4) = etgt !tag with parent edge
                        end if
                    end if 
                else !Take edge segments where both ends of the segment are associated with the cell ctgt 
                   
                    do ss=1,edge_clip_sm(etgt)%nint+1
                        v1 = edge_clip_sm(etgt)%edge_mesh(ss,1)
                        v2 = edge_clip_sm(etgt)%edge_mesh(ss,2)
                        if (ss == 1) then !end 1 to first intersect

                            vintvalid1 = 0 
                            do vv=1,vtx_clip_sm(v1)%nintersect
                                ftgt = vtx_clip_sm(v1)%face_int_idx(vv)
                                cl = volume_mesh_full%faces(ftgt)%cleft
                                cr = volume_mesh_full%faces(ftgt)%cright
                                if (cl == ctgt) then
                                    vintvalid1 = 1
                                    exit
                                end if 
                                if (cr == ctgt) then
                                    vintvalid1 = 1
                                    exit
                                end if 
                            end do 
                            if (eend_cell(ee,1) == ctgt) then 
                                vintvalid1 = 1
                            end if 

                            vintvalid2 = 0 
                            do aa=1,edge_clip_sm(etgt)%nfint(abs(v2))
                                ftgt = edge_clip_sm(etgt)%vmface(abs(v2),aa)
                                cl = volume_mesh_full%faces(ftgt)%cleft
                                cr = volume_mesh_full%faces(ftgt)%cright
                                if (cl == ctgt) then
                                    vintvalid2 = 1
                                    exit
                                end if 
                                if (cr == ctgt) then
                                    vintvalid2 = 1
                                    exit
                                end if 
                            end do 

                            if ((vintvalid1 == 1) .AND. (vintvalid2 == 1)) then !both ends of this segment intersect the face
                                nedge_tface = nedge_tface + 1
                                edges_on_tri(nedge_tface,1) = v1
                                edges_on_tri(nedge_tface,2) = edge_clip_sm(etgt)%vtx_idx(abs(v2))
                                edges_on_tri(nedge_tface,3) = 3 !tag as triangle sub-edge 
                                edges_on_tri(nedge_tface,4) = etgt !tag with parent edge
                                ! print *, 'v1 - int'
                            end if
                        elseif (ss == edge_clip_sm(etgt)%nint+1) then !last intersect to end 2

                            vintvalid1 = 0 
                            do aa=1,edge_clip_sm(etgt)%nfint(abs(v1))
                                ftgt = edge_clip_sm(etgt)%vmface(abs(v1),aa)
                                cl = volume_mesh_full%faces(ftgt)%cleft
                                cr = volume_mesh_full%faces(ftgt)%cright
                                if (cl == ctgt) then
                                    vintvalid1 = 1
                                    exit
                                end if 
                                if (cr == ctgt) then
                                    vintvalid1 = 1
                                    exit
                                end if 
                            end do 

                            vintvalid2 = 0 
                            do vv=1,vtx_clip_sm(v2)%nintersect
                                ftgt = vtx_clip_sm(v2)%face_int_idx(vv)
                                cl = volume_mesh_full%faces(ftgt)%cleft
                                cr = volume_mesh_full%faces(ftgt)%cright
                                if (cl == ctgt) then
                                    vintvalid2 = 1
                                    exit
                                end if 
                                if (cr == ctgt) then
                                    vintvalid2 = 1
                                    exit
                                end if 
                            end do 
                            if (eend_cell(ee,2) == ctgt) then 
                                vintvalid2 = 1
                            end if

                            if ((vintvalid1 == 1) .AND. (vintvalid2 == 1)) then !both ends of this segment intersect the face
                                nedge_tface = nedge_tface + 1
                                edges_on_tri(nedge_tface,1) = edge_clip_sm(etgt)%vtx_idx(abs(v1))
                                edges_on_tri(nedge_tface,2) = v2
                                edges_on_tri(nedge_tface,3) = 3 !tag as triangle sub-edge 
                                edges_on_tri(nedge_tface,4) = etgt !tag with parent edge
                                ! print *, 'int - v2'
                            end if
                        else !mid intersects 

                            vintvalid1 = 0 
                            do aa=1,edge_clip_sm(etgt)%nfint(abs(v1))
                                ftgt = edge_clip_sm(etgt)%vmface(abs(v1),aa)
                                cl = volume_mesh_full%faces(ftgt)%cleft
                                cr = volume_mesh_full%faces(ftgt)%cright
                                if (cl == ctgt) then
                                    vintvalid1 = 1
                                    exit
                                end if 
                                if (cr == ctgt) then
                                    vintvalid1 = 1
                                    exit
                                end if 
                            end do 

                            vintvalid2 = 0 
                            do aa=1,edge_clip_sm(etgt)%nfint(abs(v2))
                                ftgt = edge_clip_sm(etgt)%vmface(abs(v2),aa)
                                cl = volume_mesh_full%faces(ftgt)%cleft
                                cr = volume_mesh_full%faces(ftgt)%cright
                                if (cl == ctgt) then
                                    vintvalid2 = 1
                                    exit
                                end if 
                                if (cr == ctgt) then
                                    vintvalid2 = 1
                                    exit
                                end if 
                            end do 

                            if ((vintvalid1 == 1) .AND. (vintvalid2 == 1)) then !both ends of this segment intersect the face
                                nedge_tface = nedge_tface + 1
                                edges_on_tri(nedge_tface,1) = edge_clip_sm(etgt)%vtx_idx(abs(v1))
                                edges_on_tri(nedge_tface,2) = edge_clip_sm(etgt)%vtx_idx(abs(v2))
                                edges_on_tri(nedge_tface,3) = 3 !tag as triangle sub-edge 
                                edges_on_tri(nedge_tface,4) = etgt !tag with parent edge
                                ! print *, 'int - int'
                            end if
                        end if
                    end do 

                end if 
            end do 

            !Test if all edges_on_tri are sub edges from the same parent edge
            diff_parent = 1
            all_te_subedge = 1
            do ee=1,nedge_tface
                if (edges_on_tri(ee,3) .NE. 3) then 
                    all_te_subedge = 0 
                end if 
            end do 
            if (all_te_subedge == 1) then 
                diff_parent = 0 
                do ee=1,nedge_tface
                    if (edges_on_tri(ee,4) .NE. edges_on_tri(1,4)) then
                        diff_parent = 1
                        exit
                    end if 
                end do 
            end if 

            !Build clipped face if there are at least three clipped edges within this cell and all are not sub-edges from the same pareent surface edge
            !- if zero nedge_tface then this face clips the cell on a face vertex or cell edge so no clipping to this cell is required 
            !- if one nedge_tface then one edge of this face is parallel and coincident with an edge or face of the cell so no clipping to this cell is required 
            if ((nedge_tface .GE. 3) .AND. (diff_parent == 1)) then 

                !Build local v2e on these clipped edges 
                do ee=1,nedge_tface
                    v2e_flocal(edges_on_tri(ee,1),:) = 0
                    v2e_flocal(edges_on_tri(ee,2),:) = 0  
                end do 
                do ee=1,nedge_tface
                    if (v2e_flocal(edges_on_tri(ee,1),1) == 0) then 
                        v2e_flocal(edges_on_tri(ee,1),1) = ee 
                    else
                        v2e_flocal(edges_on_tri(ee,1),2) = ee 
                    end if 
                    if (v2e_flocal(edges_on_tri(ee,2),1) == 0) then 
                        v2e_flocal(edges_on_tri(ee,2),1) = ee 
                    else
                        v2e_flocal(edges_on_tri(ee,2),2) = ee 
                    end if 
                end do 

                !Debug print (if non contiguous face) =============================
                if ((minval(v2e_flocal(edges_on_tri(1:nedge_tface,1),1)) == 0) .OR. &
                    (minval(v2e_flocal(edges_on_tri(1:nedge_tface,1),2)) == 0) .OR. &
                    (minval(v2e_flocal(edges_on_tri(1:nedge_tface,2),1)) == 0) .OR. &
                    (minval(v2e_flocal(edges_on_tri(1:nedge_tface,2),2)) == 0)) then 
                ! if (ff == 25784) then 
                    print *, 'zero v2e at surface face = ', ff,' ====================='
                    print *, 'cell :',ctgt
                    print *, 'unclipped edge assignments :',tedge_cell(:)
                    print *, 'face edges'
                    do ee=1,nedge_tface
                        print *, edges_on_tri(ee,:)
                    end do 
                    open(11,file='io/face_cvtx.dat') 
                    do ee=1,nedge_tface
                        write(11,*) surface_mesh%vertices_full(edges_on_tri(ee,1),:),&
                                    surface_mesh%vertices_full(edges_on_tri(ee,2),:)
                    end do 
                    close(11)
                    cycle
                end if

                !Set base edge and vertex to mesh this face from 
                ebase = 1
                vbase = edges_on_tri(ebase,2) !set to end 2 for the correct direction 
                v0 = vbase
                edge_face_idx(ebase) = 1

                !Initialise this face 
                face_nvtx = 1
                vtx_face_pos(vbase) = 1
                
                !Accumulate vertices to this triangle sub-face
                edge_face_idx(1:nedge_tface) = 0  
                do ee=1,nedge_tface

                    !Find next edge
                    enew = 0 
                    ! vbase = face_ledge(ebase,1)
                    if (edge_face_idx(v2e_flocal(vbase,1)) == 0) then !take edge 1
                        enew = v2e_flocal(vbase,1)
                    elseif (edge_face_idx(v2e_flocal(vbase,2)) == 0) then !take edge 2
                        enew = v2e_flocal(vbase,2)
                    else !no unset edges so face is closed hence exit 
                        exit 
                    end if 

                    !If no new edge found then loop is complete to exit search 
                    if (enew == 0) then 
                        exit 
                    end if 

                    !Find next vertex
                    if (edges_on_tri(enew,1) == vbase) then 
                        vnew = edges_on_tri(enew,2)
                    else
                        vnew = edges_on_tri(enew,1)
                    end if 

                    !Exit if new vertex is the original base vertex as the loop is then complete
                    if (vnew == v0) then 
                        exit 
                    end if 

                    !Add vnew to the current face
                    face_nvtx = face_nvtx + 1
                    vtx_face_pos(vnew) = face_nvtx

                    !Add edge to this face index 
                    edge_face_idx(enew) = 1

                    !Update base 
                    ebase = enew 
                    vbase = vnew 
                end do

                !Store clipped face 
                surface_mesh%tri_clipped(ff)%nfclipped = surface_mesh%tri_clipped(ff)%nfclipped + 1
                surface_mesh%tri_clipped(ff)%face_clipped(surface_mesh%tri_clipped(ff)%nfclipped)%nvtx = face_nvtx
                allocate(surface_mesh%tri_clipped(ff)%face_clipped(surface_mesh%tri_clipped(ff)%nfclipped)%vertices(face_nvtx))

                !Set adjacency for this face 
                surface_mesh%tri_clipped(ff)%face_clipped(surface_mesh%tri_clipped(ff)%nfclipped)%cleft = -1 
                surface_mesh%tri_clipped(ff)%face_clipped(surface_mesh%tri_clipped(ff)%nfclipped)%cright = ctgt

                !Accumulate face 
                do ee=1,nedge_tface
                    vtgt = edges_on_tri(ee,1) 
                    if (vtgt .NE. 0) then 
                        surface_mesh%tri_clipped(ff)%face_clipped(surface_mesh%tri_clipped(ff)%nfclipped)%&
                        vertices(vtx_face_pos(vtgt)) = vtgt
                        edges_on_tri(ee,1) = 0 
                    end if 
                    vtgt = edges_on_tri(ee,2) 
                    if (vtgt .NE. 0) then 
                        surface_mesh%tri_clipped(ff)%face_clipped(surface_mesh%tri_clipped(ff)%nfclipped)%&
                        vertices(vtx_face_pos(vtgt)) = vtgt
                        edges_on_tri(ee,2) = 0 
                    end if 
                end do 

                !Check and set face orientation ------
                !Normal vector to the face 
                Nft(:) = 0.0d0 
                do vv=1,surface_mesh%tri_clipped(ff)%face_clipped(surface_mesh%tri_clipped(ff)%nfclipped)%nvtx
                    vc = vv 
                    vn = mod(vv,surface_mesh%tri_clipped(ff)%face_clipped(surface_mesh%tri_clipped(ff)%nfclipped)%nvtx) + 1
                    vc = surface_mesh%tri_clipped(ff)%face_clipped(surface_mesh%tri_clipped(ff)%nfclipped)%&
                    vertices(vc)
                    vn = surface_mesh%tri_clipped(ff)%face_clipped(surface_mesh%tri_clipped(ff)%nfclipped)%&
                    vertices(vn)
                    vtxC(:) = surface_mesh%vertices_full(vc,:)
                    vtxN(:) = surface_mesh%vertices_full(vn,:)
                    Nft(1) = Nft(1) - 0.5d0*(vtxN(3) + vtxC(3))*(vtxN(2) - vtxC(2))
                    Nft(2) = Nft(2) - 0.5d0*(vtxN(1) + vtxC(1))*(vtxN(3) - vtxC(3))
                    Nft(3) = Nft(3) - 0.5d0*(vtxN(2) + vtxC(2))*(vtxN(1) - vtxC(1))
                end do 

                !Flip if required 
                if (dot_product(Nf,Nft) .LT. 0.0d0) then 
                    fliparray(1:surface_mesh%tri_clipped(ff)%face_clipped(surface_mesh%tri_clipped(ff)%nfclipped)%nvtx) =&
                    surface_mesh%tri_clipped(ff)%face_clipped(surface_mesh%tri_clipped(ff)%nfclipped)%vertices(:)
                    do vv=1,surface_mesh%tri_clipped(ff)%face_clipped(surface_mesh%tri_clipped(ff)%nfclipped)%nvtx
                        vc = surface_mesh%tri_clipped(ff)%face_clipped(surface_mesh%tri_clipped(ff)%nfclipped)%nvtx - vv + 1
                        surface_mesh%tri_clipped(ff)%face_clipped(surface_mesh%tri_clipped(ff)%nfclipped)%vertices(vv) =&
                        fliparray(vc)
                    end do 
                    !print *, '*** flip face'
                end if 
            end if 

            !Debug 
            ! if (nedge_tface == 2) then 
            !     print *, nedge_tface, 'face : ',ff, '// cell :',ctgt
            !     do ee=1,nedge_tface
            !         print *, edges_on_tri(ee,:)
            !     end do 
            ! end if 
        end do
    else !Retain full triangle as not clipped if within the mesh domain 

        !Face midpoint
        fmid(1) = sum(surface_mesh%vertices(surface_mesh%connectivity(ff,:),1))/3.0d0
        fmid(2) = sum(surface_mesh%vertices(surface_mesh%connectivity(ff,:),2))/3.0d0
        fmid(3) = sum(surface_mesh%vertices(surface_mesh%connectivity(ff,:),3))/3.0d0

        !Keep if within domain 
        if ((fmid(1) .GE. cm3dopt%mesh_xmin) .AND. (fmid(1) .LE. cm3dopt%mesh_xmax)) then 
            if ((fmid(2) .GE. cm3dopt%mesh_ymin) .AND. (fmid(2) .LE. cm3dopt%mesh_ymax)) then 
                if ((fmid(3) .GE. cm3dopt%mesh_zmin) .AND. (fmid(3) .LE. cm3dopt%mesh_zmax)) then 
                    surface_mesh%tri_clipped(ff)%nfclipped = 1
                    allocate(surface_mesh%tri_clipped(ff)%face_clipped(1))
                    surface_mesh%tri_clipped(ff)%face_clipped(1)%nvtx = 3
                    allocate(surface_mesh%tri_clipped(ff)%face_clipped(1)%vertices(3))
                    surface_mesh%tri_clipped(ff)%face_clipped(1)%vertices(:) = surface_mesh%connectivity(ff,:)
                    surface_mesh%tri_clipped(ff)%face_clipped(surface_mesh%tri_clipped(ff)%nfclipped)%cleft = -1 
                    surface_mesh%tri_clipped(ff)%face_clipped(surface_mesh%tri_clipped(ff)%nfclipped)%cright = 0
                end if 
            end if 
        end if 
    end if
end do 
return 
end subroutine build_volume_clipped_smesh_faces




!Build volume mesh surface clipped faces ===========================
subroutine build_surface_clipped_vmesh_faces(volume_mesh_full,surface_mesh,surface_adtree,cm3dopt,tri_clip_sm,vtx_clip_sm,&
                                             edge_clip_sm,edge_clip_vm,cell_keep,vtx_idx_smesh,type_tgt,smesh_vtx_full)
implicit none 

!Variables - Import
integer(in) :: vtx_idx_smesh,type_tgt
integer(in), dimension(:) :: cell_keep
real(dp), dimension(:,:) :: smesh_vtx_full
type(triint), dimension(:) :: tri_clip_sm
type(edgeint), dimension(:) :: edge_clip_vm,edge_clip_sm
type(vtx_intersect), dimension(:) :: vtx_clip_sm
type(tree_data) :: surface_adtree
type(surface_data) :: surface_mesh
type(vol_mesh_data) :: volume_mesh_full
type(cm3d_options) :: cm3dopt

!Variables - Local
integer(in) :: ff,nn,kk,ii,ee,vv,ss,v1,v2,aa
integer(in) :: co_l,co_r,Tcl,Tcr,nselected,ttgt,nteselect,vs1,vs2,ve1,ve2,nedge_face,etgt,etype,nface_clip,e0
integer(in) :: ebase,vbase,enew,vnew,vtgt,vc,vn,edg_dup,nfaceb,vintvalid1,vintvalid2
integer(in) :: node_select(surface_adtree%nnode),edge_face_idx(surface_mesh%nfcs)
integer(in) :: face_nvtx(cm3dopt%max_int_size),fliparray(cm3dopt%max_int_size)
integer(in) :: te_select(surface_mesh%nfcs,2),face_ledge(surface_mesh%nfcs,2)
integer(in) :: v2e_flocal(-vtx_idx_smesh:volume_mesh_full%nvtx,2)
integer(in) :: vtx_face_idx(-vtx_idx_smesh:volume_mesh_full%nvtx)
integer(in) :: vtx_face_pos(-vtx_idx_smesh:volume_mesh_full%nvtx)
real(dp) :: cpadSZ,zxmin,zxmax,zymin,zymax,zzmin,zzmax
real(dp) :: Nf(3),Nft(3),vtxC(3),vtxN(3)

!Clip each face
nselected = 0 
nteselect = 0 
nedge_face = 0
face_nvtx(:) = 0  
fliparray(:) = 0 
node_select(:) = 0 
vtx_face_idx(:) = 0 
vtx_face_pos(:) = 0
edge_face_idx(:) = 0 
te_select(:,:) = 0 
face_ledge(:,:) = 0 
v2e_flocal(:,:) = 0 
do ff=1,volume_mesh_full%nface
    volume_mesh_full%faces(ff)%nfclipped = 0 
    co_l = volume_mesh_full%faces(ff)%cleft_ot
    co_r = volume_mesh_full%faces(ff)%cright_ot
    if (co_l .GT. 0) then 
        Tcl = cell_keep(co_l)
    else
        Tcl = 0 
    end if 
    if (co_r .GT. 0) then 
        Tcr = cell_keep(co_r)
    else
        Tcr = 0 
    end if 
    if ((Tcl == 1) .OR. (Tcr == 1)) then !If either octree cell on this face is geometry intersecting then search for intersections 

        !Build face normal vector 
        Nf = newell_normal(volume_mesh_full%faces(ff)%nvtx,volume_mesh_full%faces(ff)%vertices,volume_mesh_full%vtx)

        !Padding size 
        cpadSZ = norm2(Nf)*0.05d0

        !Intersection bounding box
        zxmin = minval(volume_mesh_full%vtx(volume_mesh_full%faces(ff)%vertices(:),1)) - cpadSZ !tgt bounding box -> xmin
        zxmax = maxval(volume_mesh_full%vtx(volume_mesh_full%faces(ff)%vertices(:),1)) + cpadSZ !tgt bounding box -> xmax
        zymin = minval(volume_mesh_full%vtx(volume_mesh_full%faces(ff)%vertices(:),2)) - cpadSZ !tgt bounding box -> ymin
        zymax = maxval(volume_mesh_full%vtx(volume_mesh_full%faces(ff)%vertices(:),2)) + cpadSZ !tgt bounding box -> ymax
        zzmin = minval(volume_mesh_full%vtx(volume_mesh_full%faces(ff)%vertices(:),3)) - cpadSZ !tgt bounding box -> zmin
        zzmax = maxval(volume_mesh_full%vtx(volume_mesh_full%faces(ff)%vertices(:),3)) + cpadSZ !tgt bounding box -> zmax

        !Identify any triangle bounding boxes that may overlap the face 
        call search_ADtree(nselected,node_select,surface_adtree,zxmin,zxmax,zymin,zymax,zzmin,zzmax)

        !If any overlaps 
        if (nselected .NE. 0) then 

            !Accumulate triangle internal edges on this face
            nteselect = 0 
            do nn=1,nselected
                do kk=1,surface_adtree%tree(node_select(nn))%nentry

                    !Target triangle
                    ttgt = surface_adtree%tree(node_select(nn))%entry(kk)

                    !Edges within the triangle
                    do ii=1,tri_clip_sm(ttgt)%nedge
                        if (tri_clip_sm(ttgt)%edges(ii,3) == ff) then 
                            nteselect = nteselect + 1
                            te_select(nteselect,:) = tri_clip_sm(ttgt)%edges(ii,1:2) 
                        end if 
                    end do 

                    !Edge segments on the edges of this triangle 
                    do ee=1,3

                        !Target edge 
                        etgt = surface_mesh%F2E(ttgt,ee)

                        !Check cases 
                        if (edge_clip_sm(etgt)%nint .NE. 0) then !Check edge mesh for segments where both ends intersect the face
                            do ss=1,edge_clip_sm(etgt)%nint+1
                                v1 = edge_clip_sm(etgt)%edge_mesh(ss,1)
                                v2 = edge_clip_sm(etgt)%edge_mesh(ss,2)
                                if (ss == 1) then !end 1 to first intersect

                                    vintvalid1 = 0 
                                    do vv=1,vtx_clip_sm(v1)%nintersect
                                        if (vtx_clip_sm(v1)%face_int_idx(vv) == ff) then 
                                            vintvalid1 = 1
                                            exit 
                                        end if
                                    end do 

                                    vintvalid2 = 0 
                                    do aa=1,edge_clip_sm(etgt)%nfint(abs(v2))
                                        if (edge_clip_sm(etgt)%vmface(abs(v2),aa) == ff) then 
                                            vintvalid2 = 1
                                            exit 
                                        end if 
                                    end do 

                                    if ((vintvalid1 == 1) .AND. (vintvalid2 == 1)) then !both ends of this segment intersect the face
                                        nteselect = nteselect + 1
                                        te_select(nteselect,1) = v1
                                        te_select(nteselect,2) = edge_clip_sm(etgt)%vtx_idx(abs(v2))
                                        ! print *, 'v1 - int'
                                    end if

                                elseif (ss == edge_clip_sm(etgt)%nint+1) then !last intersect to end 2

                                    vintvalid1 = 0 
                                    do aa=1,edge_clip_sm(etgt)%nfint(abs(v1))
                                        if (edge_clip_sm(etgt)%vmface(abs(v1),aa) == ff) then 
                                            vintvalid1 = 1
                                            exit 
                                        end if 
                                    end do 

                                    vintvalid2 = 0 
                                    do vv=1,vtx_clip_sm(v2)%nintersect
                                        if (vtx_clip_sm(v2)%face_int_idx(vv) == ff) then 
                                            vintvalid2 = 1
                                            exit 
                                        end if
                                    end do 

                                    if ((vintvalid1 == 1) .AND. (vintvalid2 == 1)) then !both ends of this segment intersect the face
                                        nteselect = nteselect + 1
                                        te_select(nteselect,1) = edge_clip_sm(etgt)%vtx_idx(abs(v1))
                                        te_select(nteselect,2) = v2
                                        ! print *, 'int - v2'
                                    end if

                                else !mid intersects 

                                    vintvalid1 = 0 
                                    do aa=1,edge_clip_sm(etgt)%nfint(abs(v1))
                                        if (edge_clip_sm(etgt)%vmface(abs(v1),aa) == ff) then 
                                            vintvalid1 = 1
                                            exit 
                                        end if 
                                    end do 

                                    vintvalid2 = 0 
                                    do aa=1,edge_clip_sm(etgt)%nfint(abs(v2))
                                        if (edge_clip_sm(etgt)%vmface(abs(v2),aa) == ff) then 
                                            vintvalid2 = 1
                                            exit 
                                        end if 
                                    end do

                                    if ((vintvalid1 == 1) .AND. (vintvalid2 == 1)) then !both ends of this segment intersect the face
                                        nteselect = nteselect + 1
                                        te_select(nteselect,1) = edge_clip_sm(etgt)%vtx_idx(abs(v1))
                                        te_select(nteselect,2) = edge_clip_sm(etgt)%vtx_idx(abs(v2))
                                        ! print *, 'int - v2'
                                    end if
                                   
                                end if
                            end do 
                        else !Check if both vertices on this edge intersect the face 
                            v1 = surface_mesh%edges(etgt,1)
                            v2 = surface_mesh%edges(etgt,2)
                            vintvalid1 = 0 
                            do vv=1,vtx_clip_sm(v1)%nintersect
                                if (vtx_clip_sm(v1)%face_int_idx(vv) == ff) then 
                                    vintvalid1 = 1
                                    exit 
                                end if
                            end do 
                            vintvalid2 = 0 
                            do vv=1,vtx_clip_sm(v2)%nintersect
                                if (vtx_clip_sm(v2)%face_int_idx(vv) == ff) then 
                                    vintvalid2 = 1
                                    exit 
                                end if
                            end do 
                            if ((vintvalid1 == 1) .AND. (vintvalid2 == 1)) then 
                                nteselect = nteselect + 1
                                te_select(nteselect,1) = v1
                                te_select(nteselect,2) = v2
                                ! print *, 'fedge vints'
                            end if 
                        end if
                    end do 
                end do 
            end do 
            
            !If any triangle edges then this face must be clipped to the surface
            if (nteselect .NE. 0) then 

                !Build local face plane mesh (combination of triangle internal edges 
                ! correct type: face edges and intersecting edge sub sections) ----------------------------
                !Build all edges of valid sections within the face
                nedge_face = 0 
                do ee=1,volume_mesh_full%faces(ff)%nvtx
                    etgt = volume_mesh_full%faces(ff)%edges(ee)
                    etype = edge_clip_vm(etgt)%type
                    if (etype == type_tgt) then !Add full edge 
                        if (volume_mesh_full%edges(etgt,1) == volume_mesh_full%faces(ff)%vertices(ee)) then 
                            ve1 = volume_mesh_full%edges(etgt,1)
                            ve2 = volume_mesh_full%edges(etgt,2)
                        else
                            ve1 = volume_mesh_full%edges(etgt,2)
                            ve2 = volume_mesh_full%edges(etgt,1)
                        end if 
                        nedge_face = nedge_face + 1
                        face_ledge(nedge_face,1) = ve1
                        face_ledge(nedge_face,2) = ve2
                    elseif (etype == 4) then !Add by intersection segment 

                        !Set direction 
                        if (volume_mesh_full%edges(etgt,1) == volume_mesh_full%faces(ff)%vertices(ee)) then !same as edge direction
                            vs1 = 1
                            vs2 = 2
                        else !opposite to edge direction 
                            vs1 = 2
                            vs2 = 1
                        end if

                        !Add segments 
                        do ii=1,edge_clip_vm(etgt)%nint+1 
                            if (edge_clip_vm(etgt)%int_seg_type(ii) == type_tgt) then 
                                nedge_face = nedge_face + 1
                                if (edge_clip_vm(etgt)%edge_mesh(ii,vs1) .GT. 0) then !mesh vertex
                                    face_ledge(nedge_face,1) = edge_clip_vm(etgt)%edge_mesh(ii,vs1)
                                else !edge intersection
                                    face_ledge(nedge_face,1) = -1*edge_clip_vm(etgt)%vtx_idx(abs(edge_clip_vm(etgt)%&
                                    edge_mesh(ii,vs1)))
                                    !print *, 'eintidx = ',edge_clip(etgt)%vtx_idx(abs(edge_clip(etgt)%edge_mesh(ii,vs1)))
                                end if 
                                if (edge_clip_vm(etgt)%edge_mesh(ii,vs2) .GT. 0) then !mesh vertex
                                    face_ledge(nedge_face,2) = edge_clip_vm(etgt)%edge_mesh(ii,vs2)
                                else !edge intersection
                                    face_ledge(nedge_face,2) = -1*edge_clip_vm(etgt)%vtx_idx(abs(edge_clip_vm(etgt)%&
                                    edge_mesh(ii,vs2)))
                                    !print *, 'eintidx = ',edge_clip(etgt)%vtx_idx(abs(edge_clip(etgt)%edge_mesh(ii,vs2)))
                                end if 
                            end if 
                        end do 
                    end if 
                end do 

                !Add triangle internal edge segments 
                nfaceb = nedge_face
                do ee=1,nteselect

                    !Check if dulicate 
                    edg_dup = 0 
                    do ii=nfaceb,nedge_face
                        if ((face_ledge(ii,1) == -1*te_select(ee,1)) .AND. (face_ledge(ii,2) == -1*te_select(ee,2))) then 
                            edg_dup = 1
                            exit 
                        elseif ((face_ledge(ii,2) == -1*te_select(ee,1)) .AND. (face_ledge(ii,1) == -1*te_select(ee,2))) then 
                            edg_dup = 1
                            exit 
                        end if 
                    end do 

                    !Add if new
                    if (edg_dup == 0) then 
                        nedge_face = nedge_face + 1 
                        face_ledge(nedge_face,1) = -1*te_select(ee,1)
                        face_ledge(nedge_face,2) = -1*te_select(ee,2)
                    else 
                        !print *, 'edupe ',ff
                    end if 
                end do 

                !Build v2e for the clipped face mesh 
                do ee=1,nedge_face
                    v2e_flocal(face_ledge(ee,1),:) = 0
                    v2e_flocal(face_ledge(ee,2),:) = 0  
                end do 
                do ee=1,nedge_face
                    if (v2e_flocal(face_ledge(ee,1),1) == 0) then 
                        v2e_flocal(face_ledge(ee,1),1) = ee 
                    else
                        v2e_flocal(face_ledge(ee,1),2) = ee 
                    end if 
                    if (v2e_flocal(face_ledge(ee,2),1) == 0) then 
                        v2e_flocal(face_ledge(ee,2),1) = ee 
                    else
                        v2e_flocal(face_ledge(ee,2),2) = ee 
                    end if 
                end do 

                !Debug print (if non contiguous face)=============================
                if ((minval(v2e_flocal(face_ledge(1:nedge_face,1),1)) == 0) .OR. &
                    (minval(v2e_flocal(face_ledge(1:nedge_face,1),2)) == 0) .OR. &
                    (minval(v2e_flocal(face_ledge(1:nedge_face,2),1)) == 0) .OR. &
                    (minval(v2e_flocal(face_ledge(1:nedge_face,2),2)) == 0)) then 
                ! if (ff == 25784) then 
                    print *, 'Face = ', ff,' ====================='
                    print *, 'face edges'
                    do ee=1,volume_mesh_full%faces(ff)%nvtx
                        etgt = volume_mesh_full%faces(ff)%edges(ee)
                        print *, '---- edge :',ee,'(',edge_clip_vm(etgt)%nint,')' 
                        print *, 'v12 = ',volume_mesh_full%edges(etgt,1),volume_mesh_full%edges(etgt,2)
                        if (edge_clip_vm(etgt)%nint .NE. 0) then
                            do nn=1,edge_clip_vm(etgt)%nint
                                print *, edge_clip_vm(etgt)%inttype(nn)
                                print *, edge_clip_vm(etgt)%vtx_idx(nn),edge_clip_vm(etgt)%intloc(nn,:)
                            end do 
                        end if 
                    end do 
                    print *, '----------- nedge_face = ',nedge_face 
                    do ee=1,nedge_face
                        print '(A,I0,A,I0,A,I0)', ' edge : ',ee,' ',face_ledge(ee,1),' ',face_ledge(ee,2)
                        print '(I0,A,I0,A,I0,A,I0)', v2e_flocal(face_ledge(ee,1),1),' ',v2e_flocal(face_ledge(ee,1),2),&
                        ' / ',v2e_flocal(face_ledge(ee,2),1),' ',v2e_flocal(face_ledge(ee,2),2)
                    end do 
                    open(11,file='io/face_cvtx.dat') 
                    do ee=1,nedge_face
                        if (face_ledge(ee,1) .GT. 0) then 
                            Nf(:) = volume_mesh_full%vtx(face_ledge(ee,1),:)
                        else
                            Nf(:) = smesh_vtx_full(abs(face_ledge(ee,1)),:) !volume_mesh_full%lsurface(ctgt)%vertices(abs(face_ledge(ee,1)),:)
                        end if 
                        if (face_ledge(ee,2) .GT. 0) then 
                            Nft(:) = volume_mesh_full%vtx(face_ledge(ee,2),:)
                        else
                            Nft(:) = smesh_vtx_full(abs(face_ledge(ee,2)),:)!volume_mesh_full%lsurface(ctgt)%vertices(abs(face_ledge(ee,2)),:)
                        end if 
                        write(11,*) Nf(:),Nft(:)
                    end do 
                    close(11)
                    ! print *, 'surf triangles: '
                    ! do nn=1,nselected
                    !     do kk=1,surface_adtree%tree(node_select(nn))%nentry
                    !         print *, surface_adtree%tree(node_select(nn))%entry(kk)
                    !     end do 
                    ! end do 
                    ! stop
                end if

                !Construct set of full faces from the clipped face mesh 
                nface_clip = 0 
                face_nvtx(:) = 0
                edge_face_idx(1:nedge_face) = 0 
                do ii=1,nedge_face

                    !Find starting edge 
                    e0 = 0 
                    do ee=1,nedge_face
                        if (edge_face_idx(ee) == 0) then 
                            e0 = ee 
                            exit 
                        end if
                    end do 

                    !Exit if no new edges located
                    if (e0 == 0) then 
                        exit 
                    end if 

                    !Increment face count 
                    nface_clip = nface_clip + 1

                    !Set initial parameters 
                    ebase = e0
                    vbase = face_ledge(ebase,2) !set to end 2 for the correct direction 
                    edge_face_idx(e0) = nface_clip 

                    !Initialise this curve 
                    vtx_face_pos(vbase) = 1
                    vtx_face_idx(vbase) = nface_clip
                    face_nvtx(nface_clip) = 1

                    !Mesh face 
                    do ee=1,nedge_face

                        !Find next edge
                        enew = 0 
                        ! vbase = face_ledge(ebase,1)
                        if (edge_face_idx(v2e_flocal(vbase,1)) == 0) then !take edge 1
                            enew = v2e_flocal(vbase,1)
                        elseif (edge_face_idx(v2e_flocal(vbase,2)) == 0) then !take edge 2
                            enew = v2e_flocal(vbase,2)
                        else !no unset edges so face is closed hence exit 
                            exit 
                        end if 

                        !If no new edge found then loop is complete to exit search 
                        if (enew == 0) then 
                            exit 
                        end if 

                        !Find next vertex
                        if (face_ledge(enew,1) == vbase) then 
                            vnew = face_ledge(enew,2)
                        else
                            vnew = face_ledge(enew,1)
                        end if 

                        !Add vnew to the current face
                        face_nvtx(nface_clip) = face_nvtx(nface_clip) + 1
                        vtx_face_pos(vnew) = face_nvtx(nface_clip)
                        vtx_face_idx(vnew) = nface_clip 

                        !Add edge to this face index 
                        edge_face_idx(enew) = nface_clip

                        !Update base 
                        ebase = enew 
                        vbase = vnew 
                    end do
                    ! print *, face_nvtx(nface_clip) 
                end do 

                !Debug =============================
                ! if (nface_clip .GT. 1) then 
                !     print *, 'fbase = ',ff,' - Nface = ',nface_clip
                !     print *, 'Face = ', ff
                !     do ee=1,4
                !         print *, volume_mesh_full%vtx(volume_mesh_full%faces(ff)%vertices(ee),:)
                !     end do 
    
                !     print *, '----------- nedge_face = ',nedge_face 
                !     do ee=1,nedge_face
                !         print '(A,I0,A,I0,A,I0)', ' edge : ',ee,' ',face_ledge(ee,1),' ',face_ledge(ee,2)
                !         print '(I0,A,I0,A,I0,A,I0)', v2e_flocal(face_ledge(ee,1),1),' ',v2e_flocal(face_ledge(ee,1),2),&
                !         ' / ',v2e_flocal(face_ledge(ee,2),1),' ',v2e_flocal(face_ledge(ee,2),2)
                !     end do 

                !     open(11,file='io/face_cvtx.dat') 
                !     do ee=1,nedge_face
                !         if (face_ledge(ee,1) .GT. 0) then 
                !             Nf(:) = volume_mesh_full%vtx(face_ledge(ee,1),:)
                !         else
                !             Nf(:) = smesh_vtx_full(abs(face_ledge(ee,1)),:) !volume_mesh_full%lsurface(ctgt)%vertices(abs(face_ledge(ee,1)),:)
                !         end if 
                !         if (face_ledge(ee,2) .GT. 0) then 
                !             Nft(:) = volume_mesh_full%vtx(face_ledge(ee,2),:)
                !         else
                !             Nft(:) = smesh_vtx_full(abs(face_ledge(ee,2)),:)!volume_mesh_full%lsurface(ctgt)%vertices(abs(face_ledge(ee,2)),:)
                !         end if 
                !         write(11,*) Nf(:),Nft(:)
                !     end do 
                !     close(11)

                !     print *,'Face : ',ff
                !     do vv=1,volume_mesh_full%faces(ff)%face_clipped(ii)%nvtx
                !         print *, volume_mesh_full%faces(ff)%face_clipped(ii)%vertices(vv),edge_face_idx(vv)
                !     end do 
                ! end if 

                !Build full clipped faces 
                volume_mesh_full%faces(ff)%nfclipped = nface_clip
                allocate(volume_mesh_full%faces(ff)%face_clipped(nface_clip))
                do ii=1,nface_clip

                    !Set properties
                    volume_mesh_full%faces(ff)%face_clipped(ii)%nvtx = face_nvtx(ii)
                    allocate(volume_mesh_full%faces(ff)%face_clipped(ii)%vertices(face_nvtx(ii)))
                    volume_mesh_full%faces(ff)%face_clipped(ii)%vertices(:) = 0 

                    !Accumulate face 
                    do ee=1,nedge_face
                        vtgt = face_ledge(ee,1) 
                        if (vtgt .NE. 0) then 
                            if (vtx_face_idx(vtgt) == ii) then 
                                volume_mesh_full%faces(ff)%face_clipped(ii)%vertices(vtx_face_pos(vtgt)) = vtgt
                                face_ledge(ee,1) = 0 
                            end if 
                        end if 
                        vtgt = face_ledge(ee,2) 
                        if (vtgt .NE. 0) then 
                            if (vtx_face_idx(vtgt) == ii) then 
                                volume_mesh_full%faces(ff)%face_clipped(ii)%vertices(vtx_face_pos(vtgt)) = vtgt
                                face_ledge(ee,2) = 0 
                            end if 
                        end if 
                    end do 

                    !Check and set face orientation ------
                    !Normal vector to the face 
                    Nft(:) = 0.0d0 
                    do vv=1,volume_mesh_full%faces(ff)%face_clipped(ii)%nvtx
                        vc = vv 
                        vn = mod(vv,volume_mesh_full%faces(ff)%face_clipped(ii)%nvtx) + 1
                        vc = volume_mesh_full%faces(ff)%face_clipped(ii)%vertices(vc)
                        vn = volume_mesh_full%faces(ff)%face_clipped(ii)%vertices(vn)
                        if (vc .GT. 0) then 
                            vtxC(:) = volume_mesh_full%vtx(vc,:)
                        else
                            vtxC(:) = smesh_vtx_full(abs(vc),:)!volume_mesh_full%lsurface(ctgt)%vertices(abs(vc),:)
                        end if 
                        if (vn .GT. 0) then 
                            vtxN(:) = volume_mesh_full%vtx(vn,:)
                        else
                            vtxN(:) = smesh_vtx_full(abs(vn),:)!volume_mesh_full%lsurface(ctgt)%vertices(abs(vn),:)
                        end if 
                        Nft(1) = Nft(1) - 0.5d0*(vtxN(3) + vtxC(3))*(vtxN(2) - vtxC(2))
                        Nft(2) = Nft(2) - 0.5d0*(vtxN(1) + vtxC(1))*(vtxN(3) - vtxC(3))
                        Nft(3) = Nft(3) - 0.5d0*(vtxN(2) + vtxC(2))*(vtxN(1) - vtxC(1))
                    end do 
                    ! print *,'---------'
                    ! print *,'Fdir = ', dot_product(Nf,Nft)
                    ! print *,Nf
                    ! print *,Nft

                    !Flip if required 
                    if (dot_product(Nf,Nft) .LT. 0.0d0) then 
                        fliparray(1:volume_mesh_full%faces(ff)%face_clipped(ii)%nvtx) = &
                        volume_mesh_full%faces(ff)%face_clipped(ii)%vertices(:)
                        do vv=1,volume_mesh_full%faces(ff)%face_clipped(ii)%nvtx
                            vc = volume_mesh_full%faces(ff)%face_clipped(ii)%nvtx - vv + 1
                            volume_mesh_full%faces(ff)%face_clipped(ii)%vertices(vv) = fliparray(vc)
                        end do 
                        ! print *, '*** flip face'
                    end if 
                end do 
            else !Not clipped by the surface so retain full face if all edges are the correct type

                !Check edge type on the face
                nedge_face = 0 
                do ee=1,volume_mesh_full%faces(ff)%nvtx
                    etgt = volume_mesh_full%faces(ff)%edges(ee)
                    etype = edge_clip_vm(etgt)%type
                    if (etype == type_tgt) then 
                        nedge_face = nedge_face + 1
                    end if
                end do 
                ! print *,volume_mesh_full%faces(ff)%nvtx,' -- ',nedge_face

                !If all edges are the correct type then add full face as is 
                if (nedge_face == volume_mesh_full%faces(ff)%nvtx) then 
                    nface_clip = 1
                    volume_mesh_full%faces(ff)%nfclipped = nface_clip
                    allocate(volume_mesh_full%faces(ff)%face_clipped(nface_clip))
                    volume_mesh_full%faces(ff)%face_clipped(1)%nvtx = volume_mesh_full%faces(ff)%nvtx
                    allocate(volume_mesh_full%faces(ff)%face_clipped(1)%vertices(volume_mesh_full%faces(ff)%nvtx))
                    volume_mesh_full%faces(ff)%face_clipped(1)%vertices(:) = volume_mesh_full%faces(ff)%vertices(:)
                else
                    volume_mesh_full%faces(ff)%nfclipped = 0
                end if 

            end if 
        else !Retain full face if each edge is the correct type type_tgt

            !Check edge type on the face
            nedge_face = 0 
            do ee=1,volume_mesh_full%faces(ff)%nvtx
                etgt = volume_mesh_full%faces(ff)%edges(ee)
                etype = edge_clip_vm(etgt)%type
                if (etype == type_tgt) then 
                    nedge_face = nedge_face + 1
                end if
            end do 
            ! print *,volume_mesh_full%faces(ff)%nvtx,' -- ',nedge_face

            !If all edges are the correct type then add full face as is 
            if (nedge_face == volume_mesh_full%faces(ff)%nvtx) then 
                nface_clip = 1
                volume_mesh_full%faces(ff)%nfclipped = nface_clip
                allocate(volume_mesh_full%faces(ff)%face_clipped(nface_clip))
                volume_mesh_full%faces(ff)%face_clipped(1)%nvtx = volume_mesh_full%faces(ff)%nvtx
                allocate(volume_mesh_full%faces(ff)%face_clipped(1)%vertices(volume_mesh_full%faces(ff)%nvtx))
                volume_mesh_full%faces(ff)%face_clipped(1)%vertices(:) = volume_mesh_full%faces(ff)%vertices(:)
            else
                volume_mesh_full%faces(ff)%nfclipped = 0
            end if 
        end if 
    else !No surface clipping so retain full face if each edge is the correct type type_tgt

        !Check edge type on the face
        nedge_face = 0 
        do ee=1,volume_mesh_full%faces(ff)%nvtx
            etgt = volume_mesh_full%faces(ff)%edges(ee)
            etype = edge_clip_vm(etgt)%type
            if (etype == type_tgt) then 
                nedge_face = nedge_face + 1
            end if
        end do 
        ! print *,volume_mesh_full%faces(ff)%nvtx,' -- ',nedge_face

        !If all edges are the correct type then add full face as is 
        if (nedge_face == volume_mesh_full%faces(ff)%nvtx) then 
            nface_clip = 1
            volume_mesh_full%faces(ff)%nfclipped = nface_clip
            allocate(volume_mesh_full%faces(ff)%face_clipped(nface_clip))
            volume_mesh_full%faces(ff)%face_clipped(1)%nvtx = volume_mesh_full%faces(ff)%nvtx
            allocate(volume_mesh_full%faces(ff)%face_clipped(1)%vertices(volume_mesh_full%faces(ff)%nvtx))
            volume_mesh_full%faces(ff)%face_clipped(1)%vertices(:) = volume_mesh_full%faces(ff)%vertices(:)
        else
            volume_mesh_full%faces(ff)%nfclipped = 0
        end if
    end if 
end do 
return 
end subroutine build_surface_clipped_vmesh_faces




!Subroutine to build clipped surface triangle internal edges ===========================
subroutine build_surftri_internal_edges(tri_clip_sm,edge_clip_sm,vtx_clip_sm,surface_mesh,volume_mesh_full,cm3dopt)
implicit none 

!Variables - Import
type(triint), dimension(:), allocatable :: tri_clip_sm
type(edgeint), dimension(:) :: edge_clip_sm
type(vtx_intersect), dimension(:) :: vtx_clip_sm
type(surface_data) :: surface_mesh
type(vol_mesh_data) :: volume_mesh_full
type(cm3d_options) :: cm3dopt

!Variables - Local
integer(in) :: ff,ee,ii,kk,nn,aa,vv
integer(in) :: nint_tri_tot,etgt,nftriint,ftgt,vtgt,nint_face,ie,iv,intnew,ev1,ev2,vvejoined
integer(in) :: faces_on_tri(surface_mesh%nfcs),ftselected(volume_mesh_full%nface)
integer(in) :: int_from_face(50,3)
real(dp) :: intloc(50,3)

!Build edges for each triangle 
nftriint = 0 
faces_on_tri(:) = 0 
ftselected(:) = 0 
int_from_face(:,:) = 0 
do ff=1,surface_mesh%nfcs

    !Count intersections of any kind on this triangle 
    nint_tri_tot = tri_clip_sm(ff)%nvtx
    do ee=1,3
        etgt = surface_mesh%F2E(ff,ee)
        nint_tri_tot = nint_tri_tot + edge_clip_sm(etgt)%nint
    end do 

    !If any intersections 
    if (nint_tri_tot .NE. 0) then 

        !Build list of volume mesh faces that have intersections on this triangle 
        nftriint = 0 
        do ii=1,tri_clip_sm(ff)%nvtx !volume mesh edge intersections if any
            etgt = tri_clip_sm(ff)%edge_int_idx(ii)
            do ee=1,4
                ftgt = volume_mesh_full%E2F(etgt,ee) 
                if (ftgt .GT. 0) then
                    if (ftselected(ftgt) == 0) then 
                        nftriint = nftriint + 1
                        faces_on_tri(nftriint) = ftgt 
                        ftselected(ftgt) = nftriint
                    end if 
                end if 
            end do 
        end do 
        do ee=1,3 !triangle edge intersections 
            etgt = surface_mesh%F2E(ff,ee)
            do ii=1,edge_clip_sm(etgt)%nint
                do aa=1,edge_clip_sm(etgt)%nfint(ii)
                    ftgt = edge_clip_sm(etgt)%vmface(ii,aa)
                    if (ftselected(ftgt) == 0) then 
                        nftriint = nftriint + 1
                        faces_on_tri(nftriint) = ftgt 
                        ftselected(ftgt) = nftriint
                    end if 
                end do 
            end do 
        end do 
        do vv=1,3 !vertex intersections 
            vtgt = surface_mesh%connectivity(ff,vv)
            do aa=1,vtx_clip_sm(vtgt)%nintersect
                ftgt = vtx_clip_sm(vtgt)%face_int_idx(aa)
                if (ftselected(ftgt) == 0) then 
                    nftriint = nftriint + 1
                    faces_on_tri(nftriint) = ftgt 
                    ftselected(ftgt) = nftriint
                end if 
            end do 
        end do 
        ftselected(faces_on_tri(1:nftriint)) = 0 

        !Build full list of intersection vertices on this triangle for each face with any on the triangle and add edges between them
        do kk=1,nftriint

            !Target face
            ftgt = faces_on_tri(kk)

            !Gather all intersections from this face 
            nint_face = 0 
            do ee=1,3 !from surface mesh edges
                etgt = surface_mesh%F2E(ff,ee)
                do ii=1,edge_clip_sm(etgt)%nint
                    do aa=1,edge_clip_sm(etgt)%nfint(ii)
                        if (edge_clip_sm(etgt)%vmface(ii,aa) == ftgt) then 
                            intnew = 1 
                            do nn=1,nint_face
                                if (int_from_face(nn,1) == edge_clip_sm(etgt)%vtx_idx(ii)) then 
                                    intnew = 0
                                    exit
                                end if 
                            end do 
                            if (intnew == 1) then 
                                nint_face = nint_face + 1
                                int_from_face(nint_face,1) = edge_clip_sm(etgt)%vtx_idx(ii)
                                int_from_face(nint_face,2) = 1 !tag as edge intersection
                                int_from_face(nint_face,3) = etgt
                                intloc(nint_face,:) = edge_clip_sm(etgt)%intloc(ii,:)
                            end if 
                        end if 
                    end do 

                end do 
            end do 
            do ii=1,tri_clip_sm(ff)%nvtx !from surface mesh face
                do aa=1,tri_clip_sm(ff)%nfint(ii)
                    if (tri_clip_sm(ff)%face_int_idx(ii,aa) == ftgt) then
                        intnew = 1 
                        do nn=1,nint_face
                            if (int_from_face(nn,1) == tri_clip_sm(ff)%vtx_idx(ii)) then 
                                intnew = 0
                                exit
                            end if 
                        end do
                        if (intnew == 1) then 
                            nint_face = nint_face + 1
                            int_from_face(nint_face,1) = tri_clip_sm(ff)%vtx_idx(ii)
                            int_from_face(nint_face,2) = 2 !tag as face intersection
                            int_from_face(nint_face,3) = tri_clip_sm(ff)%edge_int_idx(ii)
                            intloc(nint_face,:) = tri_clip_sm(ff)%intloc(ii,:)
                        end if 
                    end if 
                end do 
            end do 
            do vv=1,3 !from surface mesh face vertices 
                vtgt = surface_mesh%connectivity(ff,vv)
                do aa=1,vtx_clip_sm(vtgt)%nintersect
                    if (vtx_clip_sm(vtgt)%face_int_idx(aa) == ftgt) then 
                    
                        intnew = 1 
                        do nn=1,nint_face
                            if (int_from_face(nn,1) == vtgt) then 
                                intnew = 0
                                exit
                            end if 
                        end do

                        if (intnew == 1) then 
                            nint_face = nint_face + 1
                            int_from_face(nint_face,1) = vtgt
                            int_from_face(nint_face,2) = 3 !tag as vertex intersection
                            int_from_face(nint_face,3) = 0
                            intloc(nint_face,:) = surface_mesh%vertices(vtgt,:)
                        end if

                    end if 
                end do 
            end do 

            !Excess face-triangle intersection warning 
            if (nint_face .GT. 2) then 
                print *, ' ** warning - more than two mesh face - surface triangle intersections (vmf - smf nint) ',&
                ftgt,' - ',ff,nint_face
                do ii=1,nint_face
                    print *,int_from_face(ii,:),' // ',intloc(ii,:)
                end do 
                cycle
            end if 

            !Build edges on this triangle  
            if (nint_face == 2) then 
                if ((int_from_face(1,2) == 1) .AND. (int_from_face(2,2) == 1)) then !Both edge intersects
                    if (int_from_face(1,3) .NE. int_from_face(2,3)) then !Check if they both lie on the same surface edge as if so they will already be joined by a sub edge on there 
                        if (tri_clip_sm(ff)%nedge+1 .GT. tri_clip_sm(ff)%nitem) then !append
                            call append_tri_clip_entry(tri_clip_sm,tri_clip_sm(ff)%nitem*2,ff)
                            if (tri_clip_sm(ff)%nitem .GT. cm3dopt%max_int_size) then 
                                cm3dopt%max_int_size = tri_clip_sm(ff)%nitem
                            end if
                        end if
                        tri_clip_sm(ff)%nedge = tri_clip_sm(ff)%nedge + 1
                        tri_clip_sm(ff)%edges(tri_clip_sm(ff)%nedge,1) = int_from_face(1,1)
                        tri_clip_sm(ff)%edges(tri_clip_sm(ff)%nedge,2) = int_from_face(2,1)
                        tri_clip_sm(ff)%edges(tri_clip_sm(ff)%nedge,3) = ftgt
                    else !They both lie on the same mesh edge 
                        !Dont add as not a triangle internal edge
                    end if
                elseif ((int_from_face(1,2) == 2) .AND. (int_from_face(2,2) == 2)) then !Both face intersects
                    if (tri_clip_sm(ff)%nedge+1 .GT. tri_clip_sm(ff)%nitem) then !append
                        call append_tri_clip_entry(tri_clip_sm,tri_clip_sm(ff)%nitem*2,ff)
                        if (tri_clip_sm(ff)%nitem .GT. cm3dopt%max_int_size) then 
                            cm3dopt%max_int_size = tri_clip_sm(ff)%nitem
                        end if
                    end if
                    tri_clip_sm(ff)%nedge = tri_clip_sm(ff)%nedge + 1
                    tri_clip_sm(ff)%edges(tri_clip_sm(ff)%nedge,1) = int_from_face(1,1)
                    tri_clip_sm(ff)%edges(tri_clip_sm(ff)%nedge,2) = int_from_face(2,1)
                    tri_clip_sm(ff)%edges(tri_clip_sm(ff)%nedge,3) = ftgt
                elseif (((int_from_face(1,2) == 1) .AND. (int_from_face(2,2) == 3)) .OR. &
                        ((int_from_face(1,2) == 3) .AND. (int_from_face(2,2) == 1)) )then !One edge intersect and one vertex intersect -> add if they dont lie on the same mesh edge

                        !Find edge the edge vertex is on 
                        if (int_from_face(1,2) == 1) then 
                            ie = 1
                            iv = 2
                        else
                            ie = 2
                            iv = 1
                        end if
                        etgt = int_from_face(ie,3)

                        !If the vertex intersection is neither of the vertices on this edge then add 
                        if ((int_from_face(iv,1) .NE. surface_mesh%edges(etgt,1)) .AND.&
                            (int_from_face(iv,1) .NE. surface_mesh%edges(etgt,2))) then 
                                if (tri_clip_sm(ff)%nedge+1 .GT. tri_clip_sm(ff)%nitem) then !append
                                    call append_tri_clip_entry(tri_clip_sm,tri_clip_sm(ff)%nitem*2,ff)
                                    if (tri_clip_sm(ff)%nitem .GT. cm3dopt%max_int_size) then 
                                        cm3dopt%max_int_size = tri_clip_sm(ff)%nitem
                                    end if
                                end if
                                tri_clip_sm(ff)%nedge = tri_clip_sm(ff)%nedge + 1
                                tri_clip_sm(ff)%edges(tri_clip_sm(ff)%nedge,1) = int_from_face(1,1)
                                tri_clip_sm(ff)%edges(tri_clip_sm(ff)%nedge,2) = int_from_face(2,1)
                                tri_clip_sm(ff)%edges(tri_clip_sm(ff)%nedge,3) = ftgt
                        end if 
                elseif ((int_from_face(1,2) == 3) .AND. (int_from_face(2,2) == 3)) then !Both vertex intersects

                    !Check if both of these vertices are joined by an edge of the triangle
                    vvejoined = 0 
                    do ee=1,3
                        etgt = surface_mesh%F2E(ff,ee)
                        ev1 = surface_mesh%edges(etgt,1)
                        ev2 = surface_mesh%edges(etgt,2)
                        if ((ev1 == int_from_face(1,1)) .AND. (ev2 == int_from_face(2,1))) then 
                            vvejoined = 1
                            exit 
                        elseif ((ev1 == int_from_face(2,1)) .AND. (ev2 == int_from_face(1,1))) then 
                            vvejoined = 1
                            exit 
                        end if 
                    end do 

                    !If both are not joined by an edge of this triangle then add the internal edge 
                    if (vvejoined == 0) then 
                        if (tri_clip_sm(ff)%nedge+1 .GT. tri_clip_sm(ff)%nitem) then !append
                            call append_tri_clip_entry(tri_clip_sm,tri_clip_sm(ff)%nitem*2,ff)
                            if (tri_clip_sm(ff)%nitem .GT. cm3dopt%max_int_size) then 
                                cm3dopt%max_int_size = tri_clip_sm(ff)%nitem
                            end if
                        end if
                        tri_clip_sm(ff)%nedge = tri_clip_sm(ff)%nedge + 1
                        tri_clip_sm(ff)%edges(tri_clip_sm(ff)%nedge,1) = int_from_face(1,1)
                        tri_clip_sm(ff)%edges(tri_clip_sm(ff)%nedge,2) = int_from_face(2,1)
                        tri_clip_sm(ff)%edges(tri_clip_sm(ff)%nedge,3) = ftgt
                    end if 
                else !One of face and edge / face and vertex so add
                    if (tri_clip_sm(ff)%nedge+1 .GT. tri_clip_sm(ff)%nitem) then !append
                        call append_tri_clip_entry(tri_clip_sm,tri_clip_sm(ff)%nitem*2,ff)
                        if (tri_clip_sm(ff)%nitem .GT. cm3dopt%max_int_size) then 
                            cm3dopt%max_int_size = tri_clip_sm(ff)%nitem
                        end if
                    end if
                    tri_clip_sm(ff)%nedge = tri_clip_sm(ff)%nedge + 1
                    tri_clip_sm(ff)%edges(tri_clip_sm(ff)%nedge,1) = int_from_face(1,1)
                    tri_clip_sm(ff)%edges(tri_clip_sm(ff)%nedge,2) = int_from_face(2,1)
                    tri_clip_sm(ff)%edges(tri_clip_sm(ff)%nedge,3) = ftgt
                end if
            end if 
        end do 
        !print *, tri_clip_sm(ff)%nedge
    end if 
end do 
return 
end subroutine build_surftri_internal_edges




!Clip surface edges to volume faces subroutine ===========================
subroutine clip_surfedges2volfaces(edge_clip_sm,cell_keep,volume_mesh_full,surface_mesh,surface_adtree,ot_mesh,&
                                   cm3dopt,vtx_idx_smesh,vtx_clip_sm)
implicit none 

!Variables - Import
integer(in) :: vtx_idx_smesh
integer(in), dimension(:) :: cell_keep
type(edgeint), dimension(:), allocatable :: edge_clip_sm
type(vtx_intersect), dimension(:) :: vtx_clip_sm
type(tree_data) :: surface_adtree
type(surface_data) :: surface_mesh
type(vol_mesh_data) :: volume_mesh_full
type(octree_data) :: ot_mesh
type(cm3d_options) :: cm3dopt

!Variables - Local 
character(len=2) :: vtxi_loc_f1,vtxi_loc_f2
integer(in) :: ff,nn,kk,ee,tt,ii,aa   
integer(in) :: co_l,co_r,Tcl,Tcr,nselected,nedge_selected,ttgt,etgt,ftgt,coctree,foctree,f1,f2,vtgt1,vtgt2
integer(in) :: intstat,int_unique,exist,ev1,ev2,atend1,atend2,NintN
integer(in) :: edges2int(surface_mesh%nedge)
integer(in) :: node_select(surface_adtree%nnode),edge_select(surface_mesh%nedge)
integer(in) :: faces(6,4),vmf_tri(2,3)
integer(in), dimension(:), allocatable :: int_map,vtx_idx_temp,nfint_temp
integer(in), dimension(:,:), allocatable :: vmface_temp
real(dp) :: cpadSZ,zxmin,zxmax,zymin,zymax,zzmin,zzmax,idist,enorm
real(dp) :: Nf(3),ve1(3),ve2(3),vt1(3),vt2(3),vt3(3),vint(3),vmerge(3),vbc1(3),vbc2(3)
real(dp) :: vt1_s(3),vt2_s(3),vt3_s(3)
real(dp), dimension(:,:), allocatable :: intloc_temp

!Define octree cell faces 
faces(1,1) = 1
faces(1,2) = 5
faces(1,3) = 6
faces(1,4) = 2
faces(2,1) = 3
faces(2,2) = 7
faces(2,3) = 8
faces(2,4) = 4
faces(3,1) = 2
faces(3,2) = 6
faces(3,3) = 7
faces(3,4) = 3
faces(4,1) = 4
faces(4,2) = 8
faces(4,3) = 5
faces(4,4) = 1
faces(5,1) = 3
faces(5,2) = 4
faces(5,3) = 1
faces(5,4) = 2
faces(6,1) = 8
faces(6,2) = 7
faces(6,3) = 6
faces(6,4) = 5

!Clip each face
nselected = 0 
node_select(:) = 0 
edge_select(:) = 0
nedge_selected = 0 
edges2int(:) = 0 
vmf_tri(:,:) = 0
do ff=1,volume_mesh_full%nface
    co_l = volume_mesh_full%faces(ff)%cleft_ot
    co_r = volume_mesh_full%faces(ff)%cright_ot
    if (co_l .GT. 0) then 
        Tcl = cell_keep(co_l)
    else
        Tcl = 0 
    end if 
    if (co_r .GT. 0) then 
        Tcr = cell_keep(co_r)
    else
        Tcr = 0 
    end if 
    if ((Tcl == 1) .OR. (Tcr == 1)) then !If either octree cell on this face is geometry intersecting then search for intersections 

        !Build face normal vector 
        Nf = newell_normal(volume_mesh_full%faces(ff)%nvtx,volume_mesh_full%faces(ff)%vertices,volume_mesh_full%vtx)

        !Padding size 
        cpadSZ = norm2(Nf)*0.05d0

        !Intersection bounding box
        zxmin = minval(volume_mesh_full%vtx(volume_mesh_full%faces(ff)%vertices(:),1)) - cpadSZ !tgt bounding box -> xmin
        zxmax = maxval(volume_mesh_full%vtx(volume_mesh_full%faces(ff)%vertices(:),1)) + cpadSZ !tgt bounding box -> xmax
        zymin = minval(volume_mesh_full%vtx(volume_mesh_full%faces(ff)%vertices(:),2)) - cpadSZ !tgt bounding box -> ymin
        zymax = maxval(volume_mesh_full%vtx(volume_mesh_full%faces(ff)%vertices(:),2)) + cpadSZ !tgt bounding box -> ymax
        zzmin = minval(volume_mesh_full%vtx(volume_mesh_full%faces(ff)%vertices(:),3)) - cpadSZ !tgt bounding box -> zmin
        zzmax = maxval(volume_mesh_full%vtx(volume_mesh_full%faces(ff)%vertices(:),3)) + cpadSZ !tgt bounding box -> zmax

        !Identify any triangle bounding boxes that may overlap the face 
        call search_ADtree(nselected,node_select,surface_adtree,zxmin,zxmax,zymin,zymax,zzmin,zzmax)

        !If any overlaps then clip the edges of this triangle to the face exactly if they intersect
        if (nselected .NE. 0) then 

            !Octree cell and face of this cell
            coctree = volume_mesh_full%faces(ff)%cleft_ot
            foctree = volume_mesh_full%faces(ff)%fot

            !Build minimal triangulation of this volume mesh face from the octree structure 
            vmf_tri(1,1) = ot_mesh%cell_vcnr(coctree,faces(foctree,1))
            vmf_tri(1,2) = ot_mesh%cell_vcnr(coctree,faces(foctree,2))
            vmf_tri(1,3) = ot_mesh%cell_vcnr(coctree,faces(foctree,3))
            vmf_tri(2,1) = ot_mesh%cell_vcnr(coctree,faces(foctree,1))
            vmf_tri(2,2) = ot_mesh%cell_vcnr(coctree,faces(foctree,3))
            vmf_tri(2,3) = ot_mesh%cell_vcnr(coctree,faces(foctree,4))
            
            !Build list of edges on the targeted triangles 
            nedge_selected = 0 
            do nn=1,nselected
                do kk=1,surface_adtree%tree(node_select(nn))%nentry
                    ttgt = surface_adtree%tree(node_select(nn))%entry(kk)
                    do ee=1,3
                        etgt = surface_mesh%F2E(ttgt,ee)
                        if (edge_select(etgt) == 0) then 
                            nedge_selected = nedge_selected + 1
                            edges2int(nedge_selected) = etgt 
                            edge_select(etgt) = nedge_selected
                        end if 
                    end do 
                end do  
            end do 
            edge_select(edges2int(1:nedge_selected)) = 0 

            !Intersect selected edges with the volume mesh face 
            do ee=1,nedge_selected

                !Target edge 
                etgt = edges2int(ee)

                !Edge end vertices 
                ev1 = surface_mesh%edges(etgt,1)
                ev2 = surface_mesh%edges(etgt,2)
                ve1(:) = surface_mesh%vertices(surface_mesh%edges(etgt,1),:)
                ve2(:) = surface_mesh%vertices(surface_mesh%edges(etgt,2),:)
                enorm = norm2(ve2(:) - ve1(:))

                !Edge adjacent faces 
                f1 = surface_mesh%E2F(etgt,1)
                f2 = surface_mesh%E2F(etgt,2)

                !Test intersections
                do tt=1,2      

                    !Face triangulation vertices
                    vt1(:) = ot_mesh%vtx(vmf_tri(tt,1),:)
                    vt2(:) = ot_mesh%vtx(vmf_tri(tt,2),:)
                    vt3(:) = ot_mesh%vtx(vmf_tri(tt,3),:)

                    !Check 
                    intstat = line_tri_intersect_bool(ve1,ve2,vt1,vt2,vt3)

                    !Find intersection location and add if valid 
                    if (intstat == 1) then 

                        !Intersection location 
                        vint(:) = line_tri_intersect(ve1,ve2,vt1,vt2,vt3)
                        
                        !Initialise as unique 
                        int_unique = 0 

                        !Check if unique on this edge with other edge intersections
                        do ii=1,edge_clip_sm(etgt)%nint
                            idist = norm2(vint(:) - edge_clip_sm(etgt)%intloc(ii,:))
                            if (idist .LE. set_cmp_prc(enorm*cm3dopt%intcointol,min_precision)) then 
                                int_unique = -1
                            end if
                        end do 

                        !Check if either vertex at the end of this edge already intersects with this face
                        atend1 = 0 
                        do ii=1,vtx_clip_sm(ev1)%nintersect
                            if (vtx_clip_sm(ev1)%face_int_idx(ii) == ff) then 
                                atend1 = 1
                                exit 
                            end if 
                        end do 
                        atend2 = 0 
                        do ii=1,vtx_clip_sm(ev2)%nintersect
                            if (vtx_clip_sm(ev2)%face_int_idx(ii) == ff) then 
                                atend2 = 1
                                exit 
                            end if 
                        end do 
                        if ((atend1 == 1) .OR. (atend2 == 1)) then 
                            int_unique = -1
                        end if 

                        !If within tollerance of edge end then set the intersect vertex index as that end of the edge 
                        if (int_unique == 0) then !assign with barycentric coorinates on faces adjacent to this edge 

                            if (f1 .GT. 0) then 
                                vt1_s(:) = surface_mesh%vertices(surface_mesh%connectivity(f1,1),:)
                                vt2_s(:) = surface_mesh%vertices(surface_mesh%connectivity(f1,2),:)
                                vt3_s(:) = surface_mesh%vertices(surface_mesh%connectivity(f1,3),:)
                                vbc1 = cart2baryc_onplane(vint,vt1_s,vt2_s,vt3_s)
                                vtxi_loc_f1 = vtx_bary_tri_loc(vbc1(1),vbc1(2),vbc1(3),cm3dopt%intcointol)
                                if (vtxi_loc_f1(1:1) == 'v') then 
                                    if (vtxi_loc_f1 == 'v1') then !on vertex 1
                                        vtgt1 = surface_mesh%connectivity(f1,1)
                                    elseif (vtxi_loc_f1 == 'v2') then !on vertex 2
                                        vtgt1 = surface_mesh%connectivity(f1,2)
                                    elseif (vtxi_loc_f1 == 'v3') then !on vertex 3
                                        vtgt1 = surface_mesh%connectivity(f1,3)
                                    end if
                                    int_unique = vtgt1
                                    vmerge(:) = surface_mesh%vertices(vtgt1,:)
                                end if 
                            else
                                vt1_s(:) = surface_mesh%vertices(surface_mesh%connectivity(f2,1),:)
                                vt2_s(:) = surface_mesh%vertices(surface_mesh%connectivity(f2,2),:)
                                vt3_s(:) = surface_mesh%vertices(surface_mesh%connectivity(f2,3),:)
                                vbc2 = cart2baryc_onplane(vint,vt1_s,vt2_s,vt3_s)
                                vtxi_loc_f2 = vtx_bary_tri_loc(vbc2(1),vbc2(2),vbc2(3),cm3dopt%intcointol)
                                if (vtxi_loc_f2(1:1) == 'v') then 
                                    if (vtxi_loc_f2 == 'v1') then !on vertex 1
                                        vtgt2 = surface_mesh%connectivity(f2,1)
                                    elseif (vtxi_loc_f2 == 'v2') then !on vertex 2
                                        vtgt2 = surface_mesh%connectivity(f2,2)
                                    elseif (vtxi_loc_f2 == 'v3') then !on vertex 3
                                        vtgt2 = surface_mesh%connectivity(f2,3)
                                    end if 
                                    int_unique = vtgt2
                                    vmerge(:) = surface_mesh%vertices(vtgt2,:)
                                end if 
                            end if 
                        end if 

                        !Add to surface mesh edge 
                        if (int_unique == 0) then !Add new intersect
                            if (edge_clip_sm(etgt)%nint+1 .GT. edge_clip_sm(etgt)%nitem) then !append
                                call append_edge_clip_entry(edge_clip_sm,edge_clip_sm(etgt)%nitem*2,etgt)
                                if (edge_clip_sm(etgt)%nitem .GT. cm3dopt%max_int_size) then 
                                    cm3dopt%max_int_size = edge_clip_sm(etgt)%nitem
                                end if 
                            end if
                            edge_clip_sm(etgt)%type = 4
                            edge_clip_sm(etgt)%nint = edge_clip_sm(etgt)%nint + 1
                            edge_clip_sm(etgt)%intloc(edge_clip_sm(etgt)%nint,:) = vint(:)
                            vtx_idx_smesh = vtx_idx_smesh + 1
                            edge_clip_sm(etgt)%vtx_idx(edge_clip_sm(etgt)%nint) = vtx_idx_smesh
                            edge_clip_sm(etgt)%vmface(edge_clip_sm(etgt)%nint,1) = ff 
                            edge_clip_sm(etgt)%nfint(edge_clip_sm(etgt)%nint) = 1  
                        elseif (int_unique .GT. 0) then !Set as edge end vertex intersect 

                            !Check if this face already intersects this vertex
                            exist = 0 
                            do ii=1,vtx_clip_sm(int_unique)%nintersect
                                if (vtx_clip_sm(int_unique)%face_int_idx(ii) == ff) then 
                                    exist = 1
                                    exit 
                                end if 
                            end do 

                            !Add if new
                            if (exist == 0) then 
                                vtx_clip_sm(int_unique)%nintersect = vtx_clip_sm(int_unique)%nintersect + 1
                                vtx_clip_sm(int_unique)%face_int_idx(vtx_clip_sm(int_unique)%nintersect) = ff 
                                ! print *, 'add vint ',int_unique,'//',ff
                            end if 
                        end if 

                        !Exit search for this edge on this face to ensure there is only one intersection on this face 
                        exit 
                    end if 
                end do 
            end do 
        end if 
    end if 
end do 

!Allocate items 
allocate(nfint_temp(cm3dopt%max_int_size))
allocate(vtx_idx_temp(cm3dopt%max_int_size))
allocate(int_map(cm3dopt%max_int_size))
allocate(intloc_temp(cm3dopt%max_int_size,3))
allocate(vmface_temp(cm3dopt%max_int_size,8))

!Clean intersections -> remove intersections with each face ff on edges where one end has been tagged cooincident with face ff
do ee=1,surface_mesh%nedge
    if (edge_clip_sm(ee)%nint .NE. 0) then !if any intersections 

        !Edge ends
        ev1 = surface_mesh%edges(ee,1)
        ev2 = surface_mesh%edges(ee,2)

        !If either end has tagged vertex intersections 
        if ((vtx_clip_sm(ev1)%nintersect .GT. 0) .OR. (vtx_clip_sm(ev2)%nintersect .GT. 0)) then 
           
           !Identify any edge intersections added here that share a face with the vertex end intersections 
            NintN = 0 
            int_map(:) = 0 
            do ii=1,edge_clip_sm(ee)%nint
                ftgt = edge_clip_sm(ee)%vmface(ii,1)
                exist = 0 
                do aa=1,vtx_clip_sm(ev1)%nintersect
                    if (vtx_clip_sm(ev1)%face_int_idx(aa) == ftgt) then 
                        exist = 1
                        exit 
                    end if 
                end do 
                do aa=1,vtx_clip_sm(ev2)%nintersect
                    if (vtx_clip_sm(ev2)%face_int_idx(aa) == ftgt) then 
                        exist = 1
                        exit 
                    end if 
                end do 
                if (exist == 0) then 
                    NintN = NintN + 1
                    int_map(ii) = NintN
                end if 
            end do 

            !Remap if any have been removed
            if (NintN .LT. edge_clip_sm(ee)%nint) then 
                if (NintN == 0) then 
                    edge_clip_sm(ee)%type = 0
                    edge_clip_sm(ee)%nint = 0 
                else



                    nfint_temp(1:edge_clip_sm(ee)%nint) = edge_clip_sm(ee)%nfint(1:edge_clip_sm(ee)%nint) 
                    vmface_temp(1:edge_clip_sm(ee)%nint,:) = edge_clip_sm(ee)%vmface(1:edge_clip_sm(ee)%nint,:) 
                    vtx_idx_temp(1:edge_clip_sm(ee)%nint) = edge_clip_sm(ee)%vtx_idx(1:edge_clip_sm(ee)%nint) 
                    intloc_temp(1:edge_clip_sm(ee)%nint,:) = edge_clip_sm(ee)%intloc(1:edge_clip_sm(ee)%nint,:) 


                    do ii=1,edge_clip_sm(ee)%nint
                        if (int_map(ii) .NE. 0) then 
                            edge_clip_sm(ee)%nfint(int_map(ii)) = nfint_temp(ii)
                            edge_clip_sm(ee)%vmface(int_map(ii),:) = vmface_temp(ii,:)
                            edge_clip_sm(ee)%vtx_idx(int_map(ii)) = vtx_idx_temp(ii)
                            edge_clip_sm(ee)%intloc(int_map(ii),:) = intloc_temp(ii,:)
                        end if 
                    end do 
                    edge_clip_sm(ee)%nint = NintN
                end if 
            end if 
        end if 
    end if
end do 
return 
end subroutine clip_surfedges2volfaces




!Subroutine to build surface mesh intersecting edge mesh ===========================
subroutine build_surf_intersect_edge_mesh(edge_clip_sm,surface_mesh,cm3dopt)
implicit none 

!Variables - Import
type(edgeint), dimension(:) :: edge_clip_sm
type(surface_data) :: surface_mesh
type(cm3d_options) :: cm3dopt

!Variables - Local 
integer(in) :: ii,ee,kk
integer(in) :: ev1,ev2,intselect
integer(in) :: int_selected(cm3dopt%max_int_size)
integer(in) :: vtx_idx_temp(cm3dopt%max_int_size),vmface_temp(cm3dopt%max_int_size,8),nfint_temp(cm3dopt%max_int_size)
real(dp) :: ve1(3),ve2(3),edir(3)
real(dp) :: intdist(cm3dopt%max_int_size)
real(dp) :: intloc_temp(cm3dopt%max_int_size,3)

!Order intersections and build edge mesh where neccesary
intdist(:) = 0.0d0 
nfint_temp(:) = 0 
int_selected(:) = 0
vtx_idx_temp(:) = 0 
vmface_temp(:,:) = 0 
intloc_temp(:,:) = 0.0d0 
do ee=1,surface_mesh%nedge
    if (edge_clip_sm(ee)%nint .NE. 0) then !if any intersections 

        !Edge ends
        ev1 = surface_mesh%edges(ee,1)
        ev2 = surface_mesh%edges(ee,2)

        !Edge end vertices 
        ve1(:) = surface_mesh%vertices(ev1,:)
        ve2(:) = surface_mesh%vertices(ev2,:)

        !Edge direction 
        edir(:) = ve2(:) - ve1(:)

        !If more than one intersection then order intersections 
        if (edge_clip_sm(ee)%nint .GT. 1) then

            !Build intersection distances from end 1
            do ii=1,edge_clip_sm(ee)%nint
                intdist(ii) = norm2(edge_clip_sm(ee)%intloc(ii,:) - ve1(:))
            end do 

            !Order 
            int_selected(:) = 0
            nfint_temp(1:edge_clip_sm(ee)%nint) = edge_clip_sm(ee)%nfint(1:edge_clip_sm(ee)%nint) 
            vmface_temp(1:edge_clip_sm(ee)%nint,:) = edge_clip_sm(ee)%vmface(1:edge_clip_sm(ee)%nint,:) 
            vtx_idx_temp(1:edge_clip_sm(ee)%nint) = edge_clip_sm(ee)%vtx_idx(1:edge_clip_sm(ee)%nint) 
            intloc_temp(1:edge_clip_sm(ee)%nint,:) = edge_clip_sm(ee)%intloc(1:edge_clip_sm(ee)%nint,:) 
            do ii=1,edge_clip_sm(ee)%nint

                !Select intersection at maximum distance
                intselect = maxloc(intdist(1:edge_clip_sm(ee)%nint),1,&
                                   mask = int_selected(1:edge_clip_sm(ee)%nint) == 0)

                !Store
                edge_clip_sm(ee)%nfint(edge_clip_sm(ee)%nint-ii+1) = nfint_temp(intselect)
                edge_clip_sm(ee)%vmface(edge_clip_sm(ee)%nint-ii+1,:) = vmface_temp(intselect,:)
                edge_clip_sm(ee)%vtx_idx(edge_clip_sm(ee)%nint-ii+1) = vtx_idx_temp(intselect)
                edge_clip_sm(ee)%intloc(edge_clip_sm(ee)%nint-ii+1,:) = intloc_temp(intselect,:)

                !Set distance to zero and set to visited
                int_selected(intselect) = 1
                intdist(intselect) = 0.0d0 
            end do 
            ! print *, vtx_idx_temp(1:edge_clip_sm(ee)%nint)
            ! print *, edge_clip_sm(ee)%vtx_idx(1:edge_clip_sm(ee)%nint)
        end if 

        !Build edge sub-mesh ---- 
        !Allocate
        allocate(edge_clip_sm(ee)%edge_mesh(edge_clip_sm(ee)%nint+1,2))
            
        !Set first segment (nagative numbers index the intersections along the edge)
        edge_clip_sm(ee)%edge_mesh(1,1) = ev1
        edge_clip_sm(ee)%edge_mesh(1,2) = -1

        !Set middle segments 
        if (edge_clip_sm(ee)%nint .GE. 2) then 
            do kk=2,edge_clip_sm(ee)%nint
                edge_clip_sm(ee)%edge_mesh(kk,1) = -1*(kk - 1)
                edge_clip_sm(ee)%edge_mesh(kk,2) = -1*kk
            end do 
        end if 

        !Set final segment 
        edge_clip_sm(ee)%edge_mesh(edge_clip_sm(ee)%nint+1,1) = -1*edge_clip_sm(ee)%nint
        edge_clip_sm(ee)%edge_mesh(edge_clip_sm(ee)%nint+1,2) = ev2
    end if 
end do 
return 
end subroutine build_surf_intersect_edge_mesh




!Clip edges to surface subroutine ===========================
subroutine clip_voledges2surffaces(edge_clip,nedge,vmf_edges,vtx_external,vtx_type1,volume_mesh_full,surface_adtree,&
                                   surface_mesh,cm3dopt,vtx_idx_smesh,cm3dfailure)
implicit none 

!Variables - Import
integer(in) :: nedge,vtx_idx_smesh,cm3dfailure
integer(in), dimension(:) :: vtx_external,vtx_type1
integer(in), dimension(:,:) :: vmf_edges
type(vol_mesh_data) :: volume_mesh_full
type(edgeint), dimension(:), allocatable :: edge_clip
type(cm3d_options) :: cm3dopt
type(surface_data) :: surface_mesh
type(tree_data) :: surface_adtree

!Variables - Local 
integer(in) :: ee,nn,kk,ii
integer(in) :: etype_check,int_unique
integer(in) :: ev1,ev2,echeck,nselected,N_edge_intersect,intnew,trinew,int_type,Neintersect,intselect,sint_idx
integer(in) :: int_selected(cm3dopt%NintEmax),edge_intersect_type(cm3dopt%NintEmax)
integer(in) :: node_select(surface_adtree%nnode),edge_intersect_tri(cm3dopt%NintEmax)
integer(in) :: edge_check(nedge)
real(dp) :: cpadSZ,zxmin,zxmax,zymin,zymax,zzmin,zzmax,segnorm!,transpert
real(dp) :: ve1(3),ve2(3),ve1p(3),ve2p(3),segdir(3),vf1(3),vf2(3),vf3(3),vint(3),vintBC(3),tri_normal(3),vmerge(3)!,vs1(3),vs2(3)
real(dp) :: edge_intersect_norm(cm3dopt%NintEmax),edge_intersect(cm3dopt%NintEmax,3),edge_intersect_bcc(cm3dopt%NintEmax,3)

!Set perturbation for transition tests 
! transpert = 1.5d0*set_cmp_prc(cm3dopt%intcointol,min_precision) !cm3dopt%intcointol 

!Allocate edge_clip structure 
allocate(edge_clip(nedge))

!Initialise
N_edge_intersect = 0 
edge_intersect_tri(:) = 0 
edge_intersect_type(:) = 0 
edge_intersect(:,:) = 0.0d0 
edge_intersect_bcc(:,:) = 0.0d0 

!Find intersections and classify edges
sint_idx = 0 
Neintersect = 0 
nselected = 0 
node_select(:) = 0 
etype_check = 0 
do ee=1,nedge
    edge_clip(ee)%type = 0 
    edge_clip(ee)%nint = 0
end do 
do ee=1,nedge

    !Initialise clipping on this edge 
    edge_clip(ee)%type = 0 
    edge_clip(ee)%nint = 0

    !Edge vertices 
    ev1 = vmf_edges(ee,1)
    ev2 = vmf_edges(ee,2)

    !If either vertex belongs to a type 1 surface intersecting cell then check this edge 
    echeck = 0 
    if ((vtx_type1(ev1) == 1) .OR. (vtx_type1(ev2) == 1)) then 
        echeck = 1
    end if 

    !Clip edge if required 
    if (echeck == 0) then !Set as fully intetnal or external based on end vertex states 
        if ((vtx_external(ev1) == 1) .AND. (vtx_external(ev2) == 1)) then !set external
            edge_clip(ee)%type = 2
        elseif ((vtx_external(ev1) == 0) .AND. (vtx_external(ev2) == 0)) then !set internal 
            edge_clip(ee)%type = 3
        else !error edge type
            print *, '** unkonwn edge case - no clip check'
            cm3dfailure = 1
            return 
        end if 
    else !Clip edge

        !Edge end vertices
        ve1(:) = volume_mesh_full%vtx(ev1,:)
        ve2(:) = volume_mesh_full%vtx(ev2,:)

        !Pad segment length 
        segdir(:) = ve2(:) - ve1(:)
        segnorm = norm2(segdir(:))
        ve2p(:) = ve2(:) + segdir(:)*cm3dopt%elenpad
        ve1p(:) = ve1(:) - segdir(:)*cm3dopt%elenpad

        !Padding size 
        cpadSZ = segnorm

        !Intersection bounding box
        zxmin = min(ve1p(1),ve2p(1)) - cpadSZ*cm3dopt%ad_padding !tgt bounding box -> xmin
        zxmax = max(ve1p(1),ve2p(1)) + cpadSZ*cm3dopt%ad_padding !tgt bounding box -> xmax
        zymin = min(ve1p(2),ve2p(2)) - cpadSZ*cm3dopt%ad_padding !tgt bounding box -> ymin
        zymax = max(ve1p(2),ve2p(2)) + cpadSZ*cm3dopt%ad_padding !tgt bounding box -> ymax
        zzmin = min(ve1p(3),ve2p(3)) - cpadSZ*cm3dopt%ad_padding !tgt bounding box -> zmin
        zzmax = max(ve1p(3),ve2p(3)) + cpadSZ*cm3dopt%ad_padding !tgt bounding box -> zmax

        !Identify any triangle bounding boxes that may overlap the edge 
        call search_ADtree(nselected,node_select,surface_adtree,zxmin,zxmax,zymin,zymax,zzmin,zzmax)

        !If any nodes are identified then search these for intersections 
        N_edge_intersect = 0 
        if (nselected .NE. 0) then 
            edge_intersect_tri(:) = 0 
            edge_intersect_type(:) = 0 
            edge_intersect(:,:) = 0.0d0 
            do nn=1,nselected
                do kk=1,surface_adtree%tree(node_select(nn))%nentry

                    !Verticies of this face
                    vf1(:) = surface_mesh%vertices(surface_mesh%connectivity(&
                    surface_adtree%tree(node_select(nn))%entry(kk),1),:)
                    vf2(:) = surface_mesh%vertices(surface_mesh%connectivity(&
                    surface_adtree%tree(node_select(nn))%entry(kk),2),:)
                    vf3(:) = surface_mesh%vertices(surface_mesh%connectivity(&
                    surface_adtree%tree(node_select(nn))%entry(kk),3),:)

                    !Test intersection
                    int_type = line_tri_intersect_bool(ve1,ve2,vf1,vf2,vf3)
                    if (int_type == 1) then 

                        !Find intersection location (is within edge and triangle here)
                        vint = line_tri_intersect(ve1,ve2,vf1,vf2,vf3)
                        vintBC = baryc_line_tri_intersect(ve1,ve2,vf1,vf2,vf3)

                        !Test against existing intersections to ensure it is new and valid --- 
                        !Check if a new surface segment
                        trinew = 1 
                        do ii=1,N_edge_intersect
                            if (surface_adtree%tree(node_select(nn))%entry(kk) == edge_intersect_tri(N_edge_intersect)) then 
                                trinew = 0 
                            end if
                        end do  

                        !Check if co-incidint with another intersection on this edge 
                        intnew = 1
                        do ii=1,N_edge_intersect
                            if (norm2(vint(:)-edge_intersect(ii,:)) .LE. set_cmp_prc(segnorm*cm3dopt%intcointol,min_precision)) then
                                intnew = 0 
                            end if 
                        end do  

                        !If near an edge end then check if co-incidint with an intersection on another edge 
                        int_unique = 0 
                        ! efrac = norm2(vint(:) -  ve1(:))/segnorm
                        ! if ((efrac .GT. 0.99d0) .OR. (efrac .LT. 0.01d0)) then 
                        !     if (efrac .GT. 0.5d0) then 
                        !         vref = ev2
                        !     else
                        !         vref = ev1
                        !     end if 
                        !     do aa=1,6 
                        !         eadj = volume_mesh_full%V2E(vref,aa)
                        !         if ((eadj .GT. 0) .AND. (eadj .NE. ee)) then 
                        !             do ii=1,edge_clip(eadj)%nint
                        !                 idist = norm2(vint(:) - edge_clip(eadj)%intloc(ii,:))
                        !                 if (idist .LE. set_cmp_prc(segnorm*cm3dopt%intcointol,min_precision)) then 
                        !                     int_unique = edge_clip(eadj)%vtx_idx(ii)
                        !                     vmerge(:) = edge_clip(eadj)%intloc(ii,:)
                        !                     exit 
                        !                 end if
                        !             end do 
                        !             if (int_unique .GT. 0) then 
                        !                 exit 
                        !             end if 
                        !         end if 
                        !     end do 
                        ! end if 

                        !Add intersection if valid and new 
                        if ((trinew == 1) .AND. (intnew == 1)) then 
                            N_edge_intersect = N_edge_intersect + 1
                            if (N_edge_intersect .GT. cm3dopt%NintEmax) then 
                                cm3dfailure = 1
                                print *, '** maximum number of edge-geometry intersections exceeded on edge ',ee
                                return 
                            end if 
                            if (int_unique == 0) then !Unique or new so add
                                edge_intersect_type(N_edge_intersect) = -1
                                edge_intersect(N_edge_intersect,:) = vint(:)
                                edge_intersect_tri(N_edge_intersect) = surface_adtree%tree(node_select(nn))%entry(kk)
                                edge_intersect_bcc(N_edge_intersect,:) = vintBC(:)
                            else !Merge with vertex int_unique
                                edge_intersect_type(N_edge_intersect) = int_unique
                                edge_intersect(N_edge_intersect,:) = vmerge(:)
                                edge_intersect_tri(N_edge_intersect) = surface_adtree%tree(node_select(nn))%entry(kk)
                                edge_intersect_bcc(N_edge_intersect,:) = vintBC(:)
                            end if 
                        end if
                    end if 
                end do 
            end do 
        end if 

        !If any intersections 
        if (N_edge_intersect .GT. 0) then !Order intersections and identify internal and external sections of the edge 

            !Count intersecting edge 
            Neintersect = Neintersect + 1

            !Set edge type to intersecting 
            edge_clip(ee)%type = 4

            !Store count
            edge_clip(ee)%nint = N_edge_intersect

            !Store reference end 1 of edge 
            edge_clip(ee)%refend1 = ev1

            !Set intersection list to ordered from edge start refend1
            int_selected(:) = 0 
            edge_intersect_norm(:) = 0.0d0  
            do kk=1,N_edge_intersect
                edge_intersect_norm(kk) = norm2(edge_intersect(kk,:) - ve1(:))
            end do 
            allocate(edge_clip(ee)%intloc(N_edge_intersect,3))
            allocate(edge_clip(ee)%intlocBC(N_edge_intersect,3))
            allocate(edge_clip(ee)%vtx_idx(N_edge_intersect))
            allocate(edge_clip(ee)%inttype(N_edge_intersect))
            allocate(edge_clip(ee)%intidx(N_edge_intersect))  
            allocate(edge_clip(ee)%surfseg(N_edge_intersect)) 
            allocate(edge_clip(ee)%intfrac(N_edge_intersect)) 
            edge_clip(ee)%vtx_idx(:) = 0 !construct later when these vertices are added to the mesh
            do kk=1,N_edge_intersect

                !Select intersection at maximum distance
                intselect = maxloc(edge_intersect_norm(1:N_edge_intersect),1,&
                                   mask = int_selected(1:N_edge_intersect) == 0)

                !Store 
                edge_clip(ee)%intloc(N_edge_intersect-kk+1,:) = edge_intersect(intselect,:)
                edge_clip(ee)%intlocBC(N_edge_intersect-kk+1,:) = edge_intersect_bcc(intselect,:)
                edge_clip(ee)%inttype(N_edge_intersect-kk+1) = edge_intersect_type(intselect)
                sint_idx = sint_idx + 1
                edge_clip(ee)%intidx(N_edge_intersect-kk+1) = sint_idx
                edge_clip(ee)%surfseg(N_edge_intersect-kk+1) = edge_intersect_tri(intselect)
                edge_clip(ee)%intfrac(N_edge_intersect-kk+1) = real(edge_intersect_norm(intselect)/segnorm,dp) 

                !Increment vertex index if not a merge
                if (edge_intersect_type(intselect) .GT. 0) then !merge
                    edge_clip(ee)%vtx_idx(N_edge_intersect-kk+1) = edge_intersect_type(intselect)
                else !normal
                    vtx_idx_smesh = vtx_idx_smesh + 1 
                    edge_clip(ee)%vtx_idx(N_edge_intersect-kk+1) = vtx_idx_smesh
                end if 

                !Set distance to zero and set to visited
                int_selected(intselect) = 1
                edge_intersect_norm(intselect) = 0.0d0 
            end do 

            !Classify each intersection as entry or exit 
            allocate(edge_clip(ee)%int_inout(N_edge_intersect))
            edge_clip(ee)%int_inout(:) = 0 
            do kk=1,N_edge_intersect
                vf1(:) = surface_mesh%vertices(surface_mesh%connectivity(edge_clip(ee)%surfseg(kk),1),:)
                vf2(:) = surface_mesh%vertices(surface_mesh%connectivity(edge_clip(ee)%surfseg(kk),2),:)
                vf3(:) = surface_mesh%vertices(surface_mesh%connectivity(edge_clip(ee)%surfseg(kk),3),:)
                tri_normal = crossp3(vf2-vf1,vf3-vf1)
                if (dot_product(tri_normal,segdir) .GT. 0.0d0) then !Exit 
                    edge_clip(ee)%int_inout(kk) = -1
                elseif (dot_product(tri_normal,segdir) .LT. 0.0d0) then !Entry 
                    edge_clip(ee)%int_inout(kk) = 1
                else
                    cm3dfailure = 1
                    print *, '** unable to identify state of edge surface intersection'
                    return 
                end if 
            end do 
            ! do kk=1,N_edge_intersect

            !     !Perturbed vertices along segment direction from the intersection 
            !     vs1(:) = edge_clip(ee)%intloc(kk,:) - real(transpert*segdir(:),dp)
            !     vs2(:) = edge_clip(ee)%intloc(kk,:) + real(transpert*segdir(:),dp)

            !     !Check and set the transition state from these two vertices 
            !     surf_transition_state = check_edge_surf_transition(vs1,vs2,surface_mesh,surface_adtree,nselected,&
            !                                                        node_select,cm3dopt)
            !     edge_clip(ee)%int_inout(kk) = surf_transition_state
 
            !     !If identification failed then set on surface dot product 
            !     if (surf_transition_state == 0) then 
            !         vf1(:) = surface_mesh%vertices(surface_mesh%connectivity(edge_clip(ee)%surfseg(kk),1),:)
            !         vf2(:) = surface_mesh%vertices(surface_mesh%connectivity(edge_clip(ee)%surfseg(kk),2),:)
            !         vf3(:) = surface_mesh%vertices(surface_mesh%connectivity(edge_clip(ee)%surfseg(kk),3),:)
            !         tri_normal = crossp3(vf2-vf1,vf3-vf1)
            !         if (dot_product(tri_normal,segdir) .GT. 0.0d0) then !Exit 
            !             edge_clip(ee)%int_inout(kk) = -1
            !         elseif (dot_product(tri_normal,segdir) .LT. 0.0d0) then !Entry 
            !             edge_clip(ee)%int_inout(kk) = 1
            !         else
            !             cm3dfailure = 1
            !             print *, '** unable to identify state of edge surface intersection'
            !             return 
            !         end if 
            !     end if 
            ! end do 

            !Classify internal/external status of each segment of this edge ----------
            !Allocate
            allocate(edge_clip(ee)%int_seg_type(N_edge_intersect+1))
            edge_clip(ee)%int_seg_type(:) = 0 

            !Set first segment 
            if ((edge_clip(ee)%int_inout(1) == 1) .AND. (vtx_external(ev1) == 1)) then !entry so segment external 
                edge_clip(ee)%int_seg_type(1) = 2
            elseif ((edge_clip(ee)%int_inout(1) == -1) .AND. (vtx_external(ev1) == 0)) then !exit so segment internal 
                edge_clip(ee)%int_seg_type(1) = 3
            else 
                cm3dfailure = 1
                print *, '** edge intersection enter/exit disagreement with vtxinternal at edge: ',ee
                print *, edge_clip(ee)%nint
                print *, edge_clip(ee)%int_inout(1)
                print *, edge_clip(ee)%type
                print *, vtx_external(ev1)
                print *, volume_mesh_full%vtx(ev1,:)
                print *, volume_mesh_full%vtx(ev2,:)
                return
            end if 

            !Set remaining segments 
            do kk=2,N_edge_intersect+1
                if (edge_clip(ee)%int_seg_type(kk-1) == 2) then !currently external
                    if (edge_clip(ee)%int_inout(kk-1) == 1) then !entry -> current segment internal
                        edge_clip(ee)%int_seg_type(kk) = 3
                    elseif (edge_clip(ee)%int_inout(kk-1) == -1) then !exit from external -> warn (keep state)
                        edge_clip(ee)%int_seg_type(kk) = 2
                        print *, '** warning: edge exiting geometry from external state '
                        cm3dfailure = 1
                    elseif (edge_clip(ee)%int_inout(kk-1) == 0) then !no transition so keep state
                        edge_clip(ee)%int_seg_type(kk) = 2
                    end if 
                elseif (edge_clip(ee)%int_seg_type(kk-1) == 3) then !currently internal
                    if (edge_clip(ee)%int_inout(kk-1) == 1) then !entry -> warn (keep state)
                        edge_clip(ee)%int_seg_type(kk) = 3
                        print *, '** warning: edge entering geometry from internal state '
                        cm3dfailure = 1
                    elseif (edge_clip(ee)%int_inout(kk-1) == -1) then !exit from internal ->  current segment external
                        edge_clip(ee)%int_seg_type(kk) = 2
                    elseif (edge_clip(ee)%int_inout(kk-1) == 0) then !no transition so keep state
                        edge_clip(ee)%int_seg_type(kk) = 3
                    end if 
                end if 
            end do 

            !Build edge intersection mesh ----------
            !Allocate
            allocate(edge_clip(ee)%edge_mesh(N_edge_intersect+1,2))
            
            !Set first segment (nagative numbers index the intersections along the edge)
            edge_clip(ee)%edge_mesh(1,1) = ev1
            edge_clip(ee)%edge_mesh(1,2) = -1

            !Set middle segments 
            if (N_edge_intersect .GE. 2) then 
                do kk=2,N_edge_intersect
                    edge_clip(ee)%edge_mesh(kk,1) = -1*(kk - 1)
                    edge_clip(ee)%edge_mesh(kk,2) = -1*kk
                end do 
            end if 

            !Set final segment 
            edge_clip(ee)%edge_mesh(N_edge_intersect+1,1) = -N_edge_intersect
            edge_clip(ee)%edge_mesh(N_edge_intersect+1,2) = ev2
        else !Set as fully intetnal or external based on end vertex states 
            if ((vtx_external(ev1) == 1) .AND. (vtx_external(ev2) == 1)) then !set external
                edge_clip(ee)%type = 2
            elseif ((vtx_external(ev1) == 0) .AND. (vtx_external(ev2) == 0)) then !set internal 
                edge_clip(ee)%type = 3
            else !Set edge type to check subsequently      !error edge type
                edge_clip(ee)%type = 5
                etype_check = etype_check + 1
                edge_check(etype_check) = ee 
            end if 
        end if 
    end if 
end do 

!Check ambigous edges and set vtx_external on these edges if required 
if (etype_check .NE. 0) then 
    call set_amb_edge_state(edge_clip,vtx_external,etype_check,edge_check,vmf_edges)
end if 

!Debug ================================
! !Debug write error edges -------
! print *, 'Necheck = ',etype_check
! open(11,file='io/edgeerror.dat')
! do ee=1,etype_check
!     kk = edge_check(ee)
!     ev1 = vmf_edges(kk,1)
!     ev2 = vmf_edges(kk,2)
!     write(11,*) volume_mesh_full%vtx(ev1,:) , volume_mesh_full%vtx(ev2,:)
! end do 
! close(11)

! !Debug write reference external vertices ------- 
! open(11,file='io/vtxexternal.dat')
! do ee=1,volume_mesh_full%nvtx
!     if (vtx_external(ee) == 1) then 
!         write(11,*) volume_mesh_full%vtx(ee,:)
!     end if
! end do 
! close(11)
return 
end subroutine clip_voledges2surffaces




!Check and set state of ambugous edges subroutine ===========================
subroutine set_amb_edge_state(edge_clip,vtx_external,etype_check,edge_check,vmf_edges)
implicit none 

!Variables - Import
integer(in) :: etype_check
integer(in), dimension(:) :: edge_check,vtx_external
integer(in), dimension(:,:) :: vmf_edges
type(edgeint), dimension(:), allocatable :: edge_clip

!Variables - Local 
integer(in) :: ff,ee 
integer(in) :: nupdate,etgt,ev1,ev2

!Flood check ambugous edges 
do ff=1,etype_check !10*etype_check
    nupdate = 0 
    do ee=1,etype_check
        etgt = edge_check(ee)
        if (edge_clip(etgt)%type == 5) then 
            nupdate = nupdate + 1
            ev1 = vmf_edges(etgt,1)
            ev2 = vmf_edges(etgt,2)
            if (vtx_external(ev1) .NE. vtx_external(ev2)) then 
                if ((vtx_external(ev1) == 1) .OR. (vtx_external(ev2) == 1)) then 
                    edge_clip(etgt)%type = 2
                    vtx_external(ev1) = 1
                    vtx_external(ev2) = 1
                end if  
            elseif ((vtx_external(ev1) == 1) .AND. (vtx_external(ev2) == 1)) then 
                edge_clip(etgt)%type = 2
            end if 
        end if 
    end do 

    !Exit if no more edges to update 
    if (nupdate == 0) then 
        exit 
    end if 
end do 
return 
end subroutine set_amb_edge_state




!Build surface mesh edges subroutine ===========================
subroutine build_smesh_edges(surface_mesh,maxvalence)
implicit none 

!Variables - Import
integer(in) :: maxvalence
type(surface_data) :: surface_mesh

!Variables - Local 
integer(in) :: ff,ee,vv 
integer(in) :: ev1,ev2,evalid,edge_idx,nedge
integer(in), dimension(:,:), allocatable :: vconnect,edgeidx 

!Initialise
allocate(vconnect(surface_mesh%nvtx,maxvalence))
allocate(edgeidx(surface_mesh%nvtx,maxvalence))
vconnect(:,:) = 0
edgeidx(:,:) = 0

!Index edges 
nedge = 0 
do ff=1,surface_mesh%nfcs
    do ee=1,3

        !Edge end vertices
        ev1 = ee
        ev2 = mod(ee,3) + 1
        ev1 = surface_mesh%connectivity(ff,ev1)
        ev2 = surface_mesh%connectivity(ff,ev2)

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
            nedge = nedge + 1

            !Add edge to connection structure     
            do vv=1,maxvalence
                if (vconnect(ev1,vv) == 0) then 
                    vconnect(ev1,vv) = ev2 
                    edgeidx(ev1,vv) = nedge
                    exit 
                end if 
            end do 
            do vv=1,maxvalence
                if (vconnect(ev2,vv) == 0) then 
                    vconnect(ev2,vv) = ev1 
                    edgeidx(ev2,vv) = nedge
                    exit 
                end if 
            end do 
        end if 
    end do 
end do 

!Build edges 
surface_mesh%nedge = nedge
allocate(surface_mesh%edges(nedge,2))
allocate(surface_mesh%F2E(surface_mesh%nfcs,3))
surface_mesh%edges(:,:) = 0 
surface_mesh%F2E(:,:) = 0 
do ff=1,surface_mesh%nfcs
    do ee=1,3

        !Edge end vertices
        ev1 = ee
        ev2 = mod(ee,3) + 1
        ev1 = surface_mesh%connectivity(ff,ev1)
        ev2 = surface_mesh%connectivity(ff,ev2)

        !Check index of this edge 
        edge_idx = 0 
        do vv=1,maxvalence
            if (vconnect(ev1,vv) == ev2) then 
                edge_idx = edgeidx(ev1,vv)
                exit 
            end if 
        end do 

        !Construction cases
        if (edge_idx .GT. 0) then !If valid index then remove from edgeidx and construct 

            !Remove edge index 
            do vv=1,maxvalence
                if (edgeidx(ev1,vv) == edge_idx) then 
                    edgeidx(ev1,vv) = -1*edgeidx(ev1,vv)
                    exit 
                end if 
            end do 
            do vv=1,maxvalence
                if (edgeidx(ev2,vv) == edge_idx) then 
                    edgeidx(ev2,vv) = -1*edgeidx(ev2,vv)
                    exit 
                end if 
            end do 

            !Build edge 
            surface_mesh%edges(edge_idx,1) = ev1
            surface_mesh%edges(edge_idx,2) = ev2

            !Add edge to this face
            surface_mesh%F2E(ff,ee) = edge_idx
        elseif (edge_idx .NE. 0) then !Edge already constructed so add to edges on this face
            surface_mesh%F2E(ff,ee) = abs(edge_idx)
        end if 
    end do 
end do 
return 
end subroutine build_smesh_edges




!Build volume mesh edges subroutine ===========================
subroutine build_vmesh_edges(nedge,vmf_edges,volume_mesh,maxvalence)
implicit none 

!Variables - Import
integer(in) :: nedge,maxvalence
integer(in), dimension(:,:), allocatable :: vmf_edges
type(vol_mesh_data) :: volume_mesh

!Variables - Local 
integer(in) :: ff,ee,vv 
integer(in) :: ev1,ev2,evalid,edge_idx
integer(in), dimension(:,:), allocatable :: vconnect,edgeidx 

!Initialise
allocate(vconnect(volume_mesh%nvtx,maxvalence))
allocate(edgeidx(volume_mesh%nvtx,maxvalence))
vconnect(:,:) = 0
edgeidx(:,:) = 0

!Index edges 
nedge = 0 
do ff=1,volume_mesh%nface
    do ee=1,volume_mesh%faces(ff)%nvtx

        !Edge end vertices
        ev1 = ee
        ev2 = mod(ee,volume_mesh%faces(ff)%nvtx) + 1
        ev1 = volume_mesh%faces(ff)%vertices(ev1)
        ev2 = volume_mesh%faces(ff)%vertices(ev2)

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
            nedge = nedge + 1

            !Add edge to connection structure     
            do vv=1,maxvalence
                if (vconnect(ev1,vv) == 0) then 
                    vconnect(ev1,vv) = ev2 
                    edgeidx(ev1,vv) = nedge
                    exit 
                end if 
            end do 
            do vv=1,maxvalence
                if (vconnect(ev2,vv) == 0) then 
                    vconnect(ev2,vv) = ev1 
                    edgeidx(ev2,vv) = nedge
                    exit 
                end if 
            end do 
        end if 
    end do 
end do 

!Build edges 
allocate(vmf_edges(nedge,2))
vmf_edges(:,:) = 0 
do ff=1,volume_mesh%nface
    do ee=1,volume_mesh%faces(ff)%nvtx

        !Edge end vertices
        ev1 = ee
        ev2 = mod(ee,volume_mesh%faces(ff)%nvtx) + 1
        ev1 = volume_mesh%faces(ff)%vertices(ev1)
        ev2 = volume_mesh%faces(ff)%vertices(ev2)

        !Check index of this edge 
        edge_idx = 0 
        do vv=1,maxvalence
            if (vconnect(ev1,vv) == ev2) then 
                edge_idx = edgeidx(ev1,vv)
                exit 
            end if 
        end do 

        !Construction cases
        if (edge_idx .GT. 0) then !If valid index then remove from edgeidx and construct 

            !Remove edge index 
            do vv=1,maxvalence
                if (edgeidx(ev1,vv) == edge_idx) then 
                    edgeidx(ev1,vv) = -1*edgeidx(ev1,vv)
                    exit 
                end if 
            end do 
            do vv=1,maxvalence
                if (edgeidx(ev2,vv) == edge_idx) then 
                    edgeidx(ev2,vv) = -1*edgeidx(ev2,vv)
                    exit 
                end if 
            end do 

            !Build edge 
            vmf_edges(edge_idx,1) = ev1
            vmf_edges(edge_idx,2) = ev2

            !Add edge to this face
            volume_mesh%faces(ff)%edges(ee) = edge_idx
        elseif (edge_idx .NE. 0) then !Edge already constructed so add to edges on this face
            volume_mesh%faces(ff)%edges(ee) = abs(edge_idx)
        end if 
    end do 
end do 

!Debug check f2e ====================
! do ff=1,volume_mesh_full%nface
!     ev1 = 0 
!     do ee=1,volume_mesh_full%faces(ff)%nvtx
!         if (volume_mesh_full%faces(ff)%edges(ee) == 0) then 
!             ev1 = 1
!         end if 
!     end do 
!     if (ev1 == 1) then 
!         ! print *, '========================================='
!         print *, 'cl / cr = ', volume_mesh_full%faces(ff)%cleft_ot,volume_mesh_full%faces(ff)%cright_ot
!         print *, ff,'-- edg',volume_mesh_full%faces(ff)%edges(:)
!         print *, ff,'-- vtx',volume_mesh_full%faces(ff)%vertices(:)
!         ! do ee=1,volume_mesh_full%faces(ff)%nvtx
!         !     print *,volume_mesh_full%vtx((volume_mesh_full%faces(ff)%vertices(ee)),:)
!         ! end do 
!         print *, 'V2V - V1'
!         do ee=1,volume_mesh_full%faces(ff)%nvtx
!             print *, 'vtx = ',ee ,' --> ',volume_mesh_full%faces(ff)%vertices(ee)
!             ev1 = volume_mesh_full%faces(ff)%vertices(ee)
!             do vv=1,maxvalence
!                 print *, vconnect(ev1,vv) , volume_mesh_full%vtx(vconnect(ev1,vv),:)

!             end do 
!         end do 
!     end if 
! end do 
return 
end subroutine build_vmesh_edges




!Build base full mesh subroutine ===========================
subroutine build_full_mesh(volume_mesh_full,ot_mesh,cell_keep,cm3dopt)
implicit none 
    
!Variables - Import
integer(in), dimension(:) :: cell_keep
type(cm3d_options) :: cm3dopt
type(octree_data) :: ot_mesh
type(vol_mesh_data) :: volume_mesh_full

!Variables - Local 
integer(in) :: cc,ff,ee,vv
integer(in) :: nvtx_face,cl_face,cr_face,fidx,exist,edge_idx,pidx
integer(in) :: cadj,fop,nface,ncell,cadj_valid,vc1,vc2,cedge,Nsubcedge_vtx,etgt_adj,ctgt_adj,fadj,vselect,vins,etgt
integer(in) :: fopposite(6),ediagonal(12),edges(12,2),faces(6,4),face2edge(6,4),edge2face(12,2),edgeop_overface(12,2)
integer(in) :: cell_index(ot_mesh%cins-1),face_index(ot_mesh%cins-1,6),edge_nvtx(4)
integer(in) :: subcedge_vtx_visit(cm3dopt%NintEmax),subcedge_vtx(cm3dopt%NintEmax,3),edge_vtx_odr(cm3dopt%NintEmax,4) !SET SIZE BASED ON MAXIMUM ADJACENCY DELTA_R
real(dp) :: subcedge_vtx_dist(cm3dopt%NintEmax)

!Define opposite faces
fopposite(1) = 2
fopposite(2) = 1
fopposite(3) = 4
fopposite(4) = 3
fopposite(5) = 6
fopposite(6) = 5

!Define diagonal edges
ediagonal(1) = 11
ediagonal(2) = 12
ediagonal(3) = 9
ediagonal(4) = 10
ediagonal(5) = 7
ediagonal(6) = 8
ediagonal(7) = 5
ediagonal(8) = 6
ediagonal(9) = 3
ediagonal(10) = 4
ediagonal(11) = 1
ediagonal(12) = 2

!Define cell faces 
faces(1,1) = 1
faces(1,2) = 5
faces(1,3) = 6
faces(1,4) = 2
faces(2,1) = 3
faces(2,2) = 7
faces(2,3) = 8
faces(2,4) = 4
faces(3,1) = 2
faces(3,2) = 6
faces(3,3) = 7
faces(3,4) = 3
faces(4,1) = 4
faces(4,2) = 8
faces(4,3) = 5
faces(4,4) = 1
faces(5,1) = 3
faces(5,2) = 4
faces(5,3) = 1
faces(5,4) = 2
faces(6,1) = 8
faces(6,2) = 7
faces(6,3) = 6
faces(6,4) = 5

!Define cell edges
edges(1,1) = 1
edges(1,2) = 2
edges(2,1) = 2
edges(2,2) = 3
edges(3,1) = 3
edges(3,2) = 4
edges(4,1) = 4
edges(4,2) = 1
edges(5,1) = 1
edges(5,2) = 5
edges(6,1) = 2
edges(6,2) = 6
edges(7,1) = 3
edges(7,2) = 7
edges(8,1) = 4
edges(8,2) = 8
edges(9,1) = 5
edges(9,2) = 6
edges(10,1) = 6
edges(10,2) = 7
edges(11,1) = 7
edges(11,2) = 8
edges(12,1) = 8
edges(12,2) = 5

!Define face2edge
face2edge(1,1) = 5
face2edge(1,2) = 9
face2edge(1,3) = 6
face2edge(1,4) = 1
face2edge(2,1) = 7
face2edge(2,2) = 11
face2edge(2,3) = 8
face2edge(2,4) = 3
face2edge(3,1) = 6
face2edge(3,2) = 10
face2edge(3,3) = 7
face2edge(3,4) = 2
face2edge(4,1) = 8
face2edge(4,2) = 12
face2edge(4,3) = 5
face2edge(4,4) = 4
face2edge(5,1) = 3
face2edge(5,2) = 4
face2edge(5,3) = 1
face2edge(5,4) = 2
face2edge(6,1) = 11
face2edge(6,2) = 10
face2edge(6,3) = 9
face2edge(6,4) = 12

!Define edge2face
edge2face(1,1) = 1
edge2face(1,2) = 5
edge2face(2,1) = 3
edge2face(2,2) = 5
edge2face(3,1) = 2
edge2face(3,2) = 5
edge2face(4,1) = 4
edge2face(4,2) = 5
edge2face(5,1) = 1
edge2face(5,2) = 4
edge2face(6,1) = 1
edge2face(6,2) = 3
edge2face(7,1) = 2
edge2face(7,2) = 3
edge2face(8,1) = 2
edge2face(8,2) = 4
edge2face(9,1) = 1
edge2face(9,2) = 6
edge2face(10,1) = 3
edge2face(10,2) = 6
edge2face(11,1) = 2
edge2face(11,2) = 6
edge2face(12,1) = 4
edge2face(12,2) = 6

!Define opposite edge when crossing edge2face
edgeop_overface(1,1) = 3
edgeop_overface(1,2) = 9
edgeop_overface(2,1) = 4
edgeop_overface(2,2) = 10
edgeop_overface(3,1) = 1
edgeop_overface(3,2) = 11
edgeop_overface(4,1) = 2
edgeop_overface(4,2) = 12
edgeop_overface(5,1) = 8
edgeop_overface(5,2) = 6
edgeop_overface(6,1) = 7
edgeop_overface(6,2) = 5
edgeop_overface(7,1) = 6
edgeop_overface(7,2) = 8
edgeop_overface(8,1) = 5
edgeop_overface(8,2) = 7
edgeop_overface(9,1) = 11
edgeop_overface(9,2) = 1
edgeop_overface(10,1) = 12
edgeop_overface(10,2) = 2
edgeop_overface(11,1) = 9
edgeop_overface(11,2) = 3
edgeop_overface(12,1) = 10
edgeop_overface(12,2) = 4

!Set maximum possible valence of each vertex in the octree
allocate(ot_mesh%vtx_valence(ot_mesh%vins-1))
ot_mesh%vtx_valence(:) = 0 
do cc=1,ot_mesh%cins-1
    if (cell_keep(cc) .NE. 0) then 
        ot_mesh%vtx_valence(ot_mesh%cell_vcnr(cc,:)) = ot_mesh%vtx_valence(ot_mesh%cell_vcnr(cc,:)) + 3
    end if
end do 
ot_mesh%maxvalence = maxval(ot_mesh%vtx_valence)

!Build octree V2V and index edges 
ot_mesh%nedge = 0 
allocate(ot_mesh%edge_index(ot_mesh%vins-1,ot_mesh%maxvalence))
allocate(ot_mesh%cell2edge(ot_mesh%cins-1,12))
allocate(ot_mesh%V2V(ot_mesh%vins-1,ot_mesh%maxvalence))
ot_mesh%edge_index(:,:) = 0 
ot_mesh%cell2edge(:,:) = 0 
ot_mesh%V2V(:,:) = 0 
do cc=1,ot_mesh%cins-1
    if (cell_keep(cc) .NE. 0) then 
        do ee=1,12

            !Edge ends 
            vc1 = ot_mesh%cell_vcnr(cc,edges(ee,1))
            vc2 = ot_mesh%cell_vcnr(cc,edges(ee,2))

            !Check against V2V 
            exist = 0 
            do vv=1,ot_mesh%maxvalence
                if (ot_mesh%V2V(vc1,vv) == vc2) then 
                    exist = 1
                    exit 
                end if 
                if (ot_mesh%V2V(vc2,vv) == vc1) then 
                    exist = 1
                    exit 
                end if 
            end do 

            !Add if valid 
            if (exist == 0) then 

                !Increment edge count 
                ot_mesh%nedge = ot_mesh%nedge + 1

                !Add edge to connection structure     
                do vv=1,ot_mesh%maxvalence
                    if (ot_mesh%V2V(vc1,vv) == 0) then 
                        ot_mesh%V2V(vc1,vv) = vc2 
                        ot_mesh%edge_index(vc1,vv) = ot_mesh%nedge
                        exit 
                    end if 
                end do 
                do vv=1,ot_mesh%maxvalence
                    if (ot_mesh%V2V(vc2,vv) == 0) then 
                        ot_mesh%V2V(vc2,vv) = vc1 
                        ot_mesh%edge_index(vc2,vv) = ot_mesh%nedge
                        exit 
                    end if 
                end do 
            end if 
        end do 
    end if
end do 

!Build octree cell2edge
do cc=1,ot_mesh%cins-1
    if (cell_keep(cc) .NE. 0) then 
        do ee=1,12

            !Edge ends 
            vc1 = ot_mesh%cell_vcnr(cc,edges(ee,1))
            vc2 = ot_mesh%cell_vcnr(cc,edges(ee,2))

            !Check index of this edge 
            edge_idx = 0 
            do vv=1,ot_mesh%maxvalence
                if (ot_mesh%V2V(vc1,vv) == vc2) then 
                    edge_idx = ot_mesh%edge_index(vc1,vv)
                    exit 
                end if 
            end do 

            !Construction cases
            if (edge_idx .GT. 0) then !If valid index then remove from ot_mesh%edge_index and construct 

                !Remove edge index 
                do vv=1,ot_mesh%maxvalence
                    if (ot_mesh%edge_index(vc1,vv) == edge_idx) then 
                        ot_mesh%edge_index(vc1,vv) = -1*ot_mesh%edge_index(vc1,vv)
                        exit 
                    end if 
                end do 
                do vv=1,ot_mesh%maxvalence
                    if (ot_mesh%edge_index(vc2,vv) == edge_idx) then 
                        ot_mesh%edge_index(vc2,vv) = -1*ot_mesh%edge_index(vc2,vv)
                        exit 
                    end if 
                end do 

                !Add edge to this cell
                ot_mesh%cell2edge(cc,ee) = edge_idx
            elseif (edge_idx .NE. 0) then !Edge already constructed so add to edges on this cell
                ot_mesh%cell2edge(cc,ee) = abs(edge_idx)
            end if 
        end do 
    end if 
end do 

! !Debug ===========================
! do cc=1,ot_mesh%cins-1
!     if (cell_keep(cc) .NE. 0) then 
!     if (minval(ot_mesh%cell2edge(cc,:)) == 0) then 
!     print '(I0,A,I0,A,I0,A,I0,A,I0,A,I0,A,I0,A,I0,A,I0,A,I0,A,I0,A,I0)', ot_mesh%cell2edge(cc,1),' ',&
!     ot_mesh%cell2edge(cc,2),' ',&
!     ot_mesh%cell2edge(cc,3),' ',ot_mesh%cell2edge(cc,4),' ',ot_mesh%cell2edge(cc,5),' ',ot_mesh%cell2edge(cc,6),&
!     ' ',ot_mesh%cell2edge(cc,7),' ',ot_mesh%cell2edge(cc,8),' ',ot_mesh%cell2edge(cc,9),' ',&
!     ot_mesh%cell2edge(cc,10),' ',ot_mesh%cell2edge(cc,11),' ',ot_mesh%cell2edge(cc,12)
!     end if 
!     ! print *, ot_mesh%cell2edge(cc,:)
!     ! do ee=1,12
!     !     if (ot_mesh%cell2edge(cc,ee) == 0) then 
!     !         print *, ot_mesh%cell2edge(cc,:)
!     !         exit 
!     !     end if
!     ! end do 
!     end if
! end do 
! print *,'NOTedge = ', ot_mesh%nedge
! print *,'maxvalence = ', ot_mesh%maxvalence 

!Accumulate octree edge-mid vertices 
allocate(ot_mesh%edge_M_vtx(ot_mesh%nedge))
do ee=1,ot_mesh%nedge
    ot_mesh%edge_M_vtx(ee)%Nevtx = 0 
    allocate(ot_mesh%edge_M_vtx(ee)%e_vertices(cm3dopt%NintEmax))
    ot_mesh%edge_M_vtx(ee)%e_vertices(:) = 0 
end do 
do cc=1,ot_mesh%cins-1
    if (cell_keep(cc) .NE. 0) then 

        !Check edges of the parent of this cell if the edge is valid 
        pidx = ot_mesh%cell_parent(cc)
        if (pidx .GT. 0) then 
            do ee=1,12

                !Edge ends 
                vc1 = ot_mesh%cell_vcnr(pidx,edges(ee,1))
                vc2 = ot_mesh%cell_vcnr(pidx,edges(ee,2))

                !Find edge 
                etgt = 0 
                do vv=1,ot_mesh%maxvalence
                    if (ot_mesh%V2V(vc1,vv) == vc2) then 
                        etgt = ot_mesh%edge_index(vc1,vv)
                        exit 
                    end if 
                end do 
                etgt = abs(etgt)

                !Search if edge is valid 
                if (etgt .GT. 0) then
                    etgt_adj = ee 
                    call get_otmesh_edge_mid_vertices(ot_mesh%edge_M_vtx(etgt)%e_vertices,ot_mesh%edge_M_vtx(etgt)%Nevtx,&
                                                      pidx,etgt_adj,ot_mesh,edges)
                end if 
            end do 
        end if 

        !Look across each face of this cell 
        do ff=1,6

            !Adjacent cell 
            cadj = ot_mesh%cell_adjacent(cc,ff)

            !Each edge of this face
            do ee=1,4

                !Corresponding cell edge indecies  
                cedge = face2edge(ff,ee)
                etgt = ot_mesh%cell2edge(cc,cedge)

                !Cell across face 1 of edge etgt
                fadj = edge2face(cedge,1)
                ctgt_adj = ot_mesh%cell_adjacent(cc,fadj)
                etgt_adj = edgeop_overface(cedge,1)
                if (ctgt_adj .GT. 0) then 
                    if (ot_mesh%cell_level(ctgt_adj) == ot_mesh%cell_level(cc)) then 
                        call get_otmesh_edge_mid_vertices(ot_mesh%edge_M_vtx(etgt)%e_vertices,ot_mesh%edge_M_vtx(etgt)%Nevtx,&
                                                            ctgt_adj,etgt_adj,ot_mesh,edges)
                    end if 
                end if 

                !Cell across face 2 of edge 
                fadj = edge2face(cedge,2)
                ctgt_adj = ot_mesh%cell_adjacent(cc,fadj)
                etgt_adj = edgeop_overface(cedge,2)
                if (ctgt_adj .GT. 0) then 
                    if (ot_mesh%cell_level(ctgt_adj) == ot_mesh%cell_level(cc)) then 
                        call get_otmesh_edge_mid_vertices(ot_mesh%edge_M_vtx(etgt)%e_vertices,ot_mesh%edge_M_vtx(etgt)%Nevtx,&
                                                            ctgt_adj,etgt_adj,ot_mesh,edges)
                    end if 
                end if
            end do 
        end do 
    end if
end do 

!Index valid faces
nface = 0 
ncell = 0 
cell_index(:) = 0 
face_index(:,:) = 0 
do cc=1,ot_mesh%cins-1
    if (cell_keep(cc) .NE. 0) then 
        ncell = ncell + 1
        cell_index(cc) = ncell 
        do ff=1,6
            if (face_index(cc,ff) == 0) then 

                !Adjacent cell 
                cadj = ot_mesh%cell_adjacent(cc,ff)

                !Check valid cell and face 
                cadj_valid = 0 
                if (cadj .GT. 0) then 
                    if (ot_mesh%cell_level(cadj) .LE. ot_mesh%cell_level(cc)) then 
                        if (cell_keep(cadj) .NE. 0) then 
                            cadj_valid = 1
                        end if 
                    end if 
                elseif (cadj .LT. 0) then 
                    cadj_valid = 1
                end if

                !Index face if valid 
                if (cadj_valid == 1) then 

                    !Index
                    nface = nface + 1
                    face_index(cc,ff) = nface

                    !Push to adjacent cell if required 
                    if (cadj .GT. 0) then 
                        if (ot_mesh%cell_level(cadj) == ot_mesh%cell_level(cc)) then 
                            face_index(cadj,fopposite(ff)) = face_index(cc,ff)
                        end if 
                    end if 
                end if 
            end if 
        end do 
    end if 
end do 

!Allocate face array
allocate(volume_mesh_full%faces(nface))

!Build valid faces
do cc=1,ot_mesh%cins-1
    if (cell_keep(cc) .NE. 0) then 
        do ff=1,6
            if (face_index(cc,ff) .GT. 0) then 

                !Face index
                fidx = face_index(cc,ff)

                !Adjacent cell and opposite face
                cadj = ot_mesh%cell_adjacent(cc,ff)
                fop = fopposite(ff)

                !Initialise sub edge structure
                edge_nvtx(:) = 0 
                edge_vtx_odr(:,:) = 0 

                !Build each edge of the face (accumulate sub cell edge vertices)
                do ee=1,4

                    !Corresponding cell edge index 
                    cedge = face2edge(ff,ee)

                    !Octree edge index 
                    etgt = ot_mesh%cell2edge(cc,cedge)

                    !Cell vertices at the ends of this edge
                    vc1 = faces(ff,ee)
                    vc2 = faces(ff,mod(ee,4)+1)

                    !Extract sub-edge vertices
                    Nsubcedge_vtx = ot_mesh%edge_M_vtx(etgt)%Nevtx
                    subcedge_vtx(1:Nsubcedge_vtx,1) = ot_mesh%edge_M_vtx(etgt)%e_vertices(1:Nsubcedge_vtx)

                    !Order all sub edge midpoint vertices along the edge by distance for edge vertex 1 
                    if (Nsubcedge_vtx .GT. 1) then 

                        !Initialise
                        subcedge_vtx_visit(:) = 0 
                        subcedge_vtx_dist(:) = 0.0d0 

                        !Get distances
                        do vv=1,Nsubcedge_vtx
                            subcedge_vtx_dist(vv) = norm2(ot_mesh%vtx(subcedge_vtx(vv,1),:) - &
                            ot_mesh%vtx(ot_mesh%cell_vcnr(cc,vc1),:))
                        end do 

                        !Order
                        do vv=1,Nsubcedge_vtx

                            !Select vertex at maximum distance
                            vselect = maxloc(subcedge_vtx_dist(1:Nsubcedge_vtx),1, mask = subcedge_vtx_visit(1:Nsubcedge_vtx) == 0)

                            !Store 
                            edge_vtx_odr(Nsubcedge_vtx-vv+1,ee) = subcedge_vtx(vselect,1)
                            subcedge_vtx_visit(vselect) = 1
                            subcedge_vtx_dist(vselect) = 0.0d0 
                        end do 
                        edge_nvtx(ee) = Nsubcedge_vtx

                        !Debug ======
                        ! print *, '------------------'
                        ! print *, Nsubcedge_vtx,' || ',&
                        ! ot_mesh%vtx(ot_mesh%cell_vcnr(cc,vc1),:),' / ',ot_mesh%vtx(ot_mesh%cell_vcnr(cc,vc2),:)
                        ! print *, subcedge_vtx(1:Nsubcedge_vtx,1)
                        ! print *, edge_vtx_odr(1:Nsubcedge_vtx,ee)
                        ! do vv=1,Nsubcedge_vtx
                        !     print *, ot_mesh%vtx(edge_vtx_odr(vv,ee),:)
                        ! end do 
                    elseif (Nsubcedge_vtx == 1) then 
                        edge_nvtx(ee) = 1 
                        edge_vtx_odr(1,ee) = subcedge_vtx(1,1) 
                    end if 
                end do 

                !Set face properties 
                nvtx_face = 4 + sum(edge_nvtx(:))
                cl_face = cell_index(cc)
                if (cadj .GT. 0) then 
                    cr_face = cell_index(cadj)
                else
                    cr_face = cadj
                end if 
                volume_mesh_full%faces(fidx)%fot = ff 
                volume_mesh_full%faces(fidx)%cleft_ot = cc
                volume_mesh_full%faces(fidx)%cright_ot = cadj

                !Build mesh face
                vins = 0
                volume_mesh_full%faces(fidx)%nvtx = nvtx_face
                volume_mesh_full%faces(fidx)%cleft = cl_face
                volume_mesh_full%faces(fidx)%cright = cr_face
                allocate(volume_mesh_full%faces(fidx)%vertices(nvtx_face))
                allocate(volume_mesh_full%faces(fidx)%edges(nvtx_face))
                volume_mesh_full%faces(fidx)%edges(:) = 0 
                do ee=1,4

                    !Cell vertices at the ends of this edge
                    vc1 = ot_mesh%cell_vcnr(cc,faces(ff,ee))
                    ! vc2 = ot_mesh%cell_vcnr(cc,faces(ff,mod(ee,4)+1))

                    !Add edge to face
                    vins = vins + 1
                    volume_mesh_full%faces(fidx)%vertices(vins) = vc1
                    if (edge_nvtx(ee) .NE. 0) then 
                        do vv=1,edge_nvtx(ee)
                            vins = vins + 1
                            volume_mesh_full%faces(fidx)%vertices(vins) = edge_vtx_odr(vv,ee)
                        end do 
                    end if 
                    ! vins = vins + 1
                    ! volume_mesh_full%faces(fidx)%vertices(vins) = vc2
                end do 

                !Set face index to constructed state
                face_index(cc,ff) = -1*face_index(cc,ff)

                !Push constructed state to adjacent cell if required 
                if (cadj .GT. 0) then 
                    if (ot_mesh%cell_level(cadj) == ot_mesh%cell_level(cc)) then 
                        face_index(cadj,fopposite(ff)) = face_index(cc,ff)
                    end if 
                end if 
            end if
        end do 
    end if 
end do 

!Set full mesh properties 
volume_mesh_full%nvtx = ot_mesh%vins-1 
volume_mesh_full%ncell = ncell 
volume_mesh_full%nface = nface 
allocate(volume_mesh_full%vtx(ot_mesh%vins-1,3))
allocate(volume_mesh_full%vtx_surfseg(ot_mesh%vins-1))
allocate(volume_mesh_full%cell_level(ncell))
allocate(volume_mesh_full%cell_otidx(ncell))
volume_mesh_full%vtx(:,:) = ot_mesh%vtx(1:ot_mesh%vins-1,:)
volume_mesh_full%vtx_surfseg(:) = 0 
do cc=1,ot_mesh%cins-1
    if (cell_index(cc) .GT. 0) then 
        volume_mesh_full%cell_otidx(cell_index(cc)) = cc 
        volume_mesh_full%cell_level(cell_index(cc)) = ot_mesh%cell_level(cc)
    end if
end do 
return 
end subroutine build_full_mesh




!Get valence subroutine ===========================
subroutine get_vmesh_valence(valence,volume_mesh)
implicit none 

!Variables - Import 
integer(in), dimension(:) :: valence
type(vol_mesh_data) :: volume_mesh

!Variables - Local
integer(in) :: ff,ee,vv
integer(in) :: e1,e2,v1,v2,vlceUB,evalid
integer(in), dimension(:,:), allocatable :: vconnect

!Find upper bound of valence
valence(:) = 0 
do ff=1,volume_mesh%nface
    do ee=1,volume_mesh%faces(ff)%nvtx
        e1 = ee 
        e2 = mod(ee,volume_mesh%faces(ff)%nvtx) + 1
        v1 = volume_mesh%faces(ff)%vertices(e1)
        v2 = volume_mesh%faces(ff)%vertices(e2)
        valence(v1) = valence(v1) + 1
        valence(v2) = valence(v2) + 1
    end do 
end do 
vlceUB = 2*maxval(valence(:))

!Construct actual valence of each vertex
allocate(vconnect(volume_mesh%nvtx,vlceUB))
valence(:) = 0
vconnect(:,:) = 0 
do ff=1,volume_mesh%nface
    do ee=1,volume_mesh%faces(ff)%nvtx

        !Ends of this edge
        e1 = ee 
        e2 = mod(ee,volume_mesh%faces(ff)%nvtx) + 1
        v1 = volume_mesh%faces(ff)%vertices(e1)
        v2 = volume_mesh%faces(ff)%vertices(e2)

        !Check against vconnect 
        evalid = 1
        do vv=1,vlceUB
            if (vconnect(v1,vv) == v2) then
                evalid = 0
                exit 
            end if 
            if (vconnect(v2,vv) == v1) then 
                evalid = 0
                exit 
            end if 
        end do 

        !Add valence if new edge
        if (evalid == 1) then 
            
            !Increment valence on each vertex
            valence(v1) = valence(v1) + 1
            valence(v2) = valence(v2) + 1

            !Update vconnect
            do vv=1,vlceUB
                if (vconnect(v1,vv) == 0) then 
                    vconnect(v1,vv) = v2
                    exit
                end if 
            end do 
            do vv=1,vlceUB
                if (vconnect(v2,vv) == 0) then 
                    vconnect(v2,vv) = v1
                    exit 
                end if 
            end do 
        end if 
    end do 
end do 
return 
end subroutine get_vmesh_valence




!Append edge intersect subroutine ===========================
subroutine append_edge_clip_entry(edge_clip,lenN,etgt)
implicit none 

!Variables - Import
integer(in) :: lenN,etgt
type(edgeint), dimension(:), allocatable :: edge_clip

!Variables - Local
integer(in) :: vtx_idxT(edge_clip(etgt)%nitem)
integer(in) :: nfintT(edge_clip(etgt)%nitem)
integer(in) :: vmfaceT(edge_clip(etgt)%nitem,8)
real(dp) :: intlocT(edge_clip(etgt)%nitem,3)

!Store data for the edge etgt
intlocT(:,:) = edge_clip(etgt)%intloc(:,:)
vtx_idxT(:) = edge_clip(etgt)%vtx_idx(:)
vmfaceT(:,:) = edge_clip(etgt)%vmface(:,:)
nfintT(:) = edge_clip(etgt)%nfint(:)

!Allocate new arrays 
deallocate(edge_clip(etgt)%intloc)
deallocate(edge_clip(etgt)%vtx_idx)
deallocate(edge_clip(etgt)%vmface)
deallocate(edge_clip(etgt)%nfint)
allocate(edge_clip(etgt)%intloc(lenN,3))
allocate(edge_clip(etgt)%vtx_idx(lenN))
allocate(edge_clip(etgt)%vmface(lenN,8))
allocate(edge_clip(etgt)%nfint(lenN))

!Store data in the new structure 
edge_clip(etgt)%intloc(1:edge_clip(etgt)%nitem,:) = intlocT(:,:) 
edge_clip(etgt)%vtx_idx(1:edge_clip(etgt)%nitem) = vtx_idxT(:)  
edge_clip(etgt)%vmface(1:edge_clip(etgt)%nitem,:) = vmfaceT(:,:)  
edge_clip(etgt)%nfint(1:edge_clip(etgt)%nitem) = nfintT(:) 
edge_clip(etgt)%intloc(edge_clip(etgt)%nitem+1:lenN,:) = 0.0d0
edge_clip(etgt)%vtx_idx(edge_clip(etgt)%nitem+1:lenN) = 0
edge_clip(etgt)%vmface(edge_clip(etgt)%nitem+1:lenN,:) = 0
edge_clip(etgt)%nfint(edge_clip(etgt)%nitem+1:lenN) = 0

!Update size of this entry
edge_clip(etgt)%nitem = lenN
return 
end subroutine append_edge_clip_entry




!Append tri intersect subroutine ===========================
subroutine append_tri_clip_entry(tri_clip,lenN,ftgt)
implicit none 

!Variables - Import
integer(in) :: lenN,ftgt
type(triint), dimension(:), allocatable :: tri_clip

!Variables - Local
integer(in) :: vtx_idxT(tri_clip(ftgt)%nitem)
integer(in) :: nfintT(tri_clip(ftgt)%nitem)
integer(in) :: edge_int_idxT(tri_clip(ftgt)%nitem)
integer(in) :: face_int_idxT(tri_clip(ftgt)%nitem,8)
integer(in) :: edges(tri_clip(ftgt)%nitem,3)
real(dp) :: intlocT(tri_clip(ftgt)%nitem,3)

!Store data for the face ftgt
vtx_idxT(:) = tri_clip(ftgt)%vtx_idx(:) 
nfintT(:) = tri_clip(ftgt)%nfint(:)
edge_int_idxT(:) = tri_clip(ftgt)%edge_int_idx(:)
face_int_idxT(:,:) = tri_clip(ftgt)%face_int_idx(:,:)
edges(:,:) = tri_clip(ftgt)%edges(:,:)
intlocT(:,:) = tri_clip(ftgt)%intloc(:,:)

!Allocate new arrays 
deallocate(tri_clip(ftgt)%edges)
deallocate(tri_clip(ftgt)%vtx_idx)
deallocate(tri_clip(ftgt)%nfint)
deallocate(tri_clip(ftgt)%face_int_idx)
deallocate(tri_clip(ftgt)%edge_int_idx)
deallocate(tri_clip(ftgt)%intloc)
allocate(tri_clip(ftgt)%edges(lenN,3))
allocate(tri_clip(ftgt)%vtx_idx(lenN))
allocate(tri_clip(ftgt)%nfint(lenN))
allocate(tri_clip(ftgt)%face_int_idx(lenN,8))
allocate(tri_clip(ftgt)%edge_int_idx(lenN))
allocate(tri_clip(ftgt)%intloc(lenN,3))

!Store data in the new structure 
tri_clip(ftgt)%edges(1:tri_clip(ftgt)%nitem,:) = edges(:,:) 
tri_clip(ftgt)%vtx_idx(1:tri_clip(ftgt)%nitem) = vtx_idxT(:) 
tri_clip(ftgt)%nfint(1:tri_clip(ftgt)%nitem) = nfintT(:)
tri_clip(ftgt)%face_int_idx(1:tri_clip(ftgt)%nitem,:) = face_int_idxT(:,:)
tri_clip(ftgt)%edge_int_idx(1:tri_clip(ftgt)%nitem) = edge_int_idxT(:) 
tri_clip(ftgt)%intloc(1:tri_clip(ftgt)%nitem,:) = intlocT(:,:)
tri_clip(ftgt)%edges(tri_clip(ftgt)%nitem+1:lenN,:) = 0
tri_clip(ftgt)%vtx_idx(tri_clip(ftgt)%nitem+1:lenN) = 0
tri_clip(ftgt)%nfint(tri_clip(ftgt)%nitem+1:lenN) = 0
tri_clip(ftgt)%face_int_idx(tri_clip(ftgt)%nitem+1:lenN,:) = 0
tri_clip(ftgt)%edge_int_idx(tri_clip(ftgt)%nitem+1:lenN) = 0
tri_clip(ftgt)%intloc(tri_clip(ftgt)%nitem+1:lenN,:) = 0.0d0 

!Update size of this entry
tri_clip(ftgt)%nitem = lenN
return 
end subroutine append_tri_clip_entry


end module cellmesh3d_mesh_build_mod