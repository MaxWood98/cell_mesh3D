!cell_mesh3d mesh building module
!Max Wood - mw16116@bristol.ac.uk
!Univeristy of Bristol - Department of Aerospace Engineering

!Version 6.2
!Updated 02-12-2023

!Module
module cellmesh3d_mesh_build_mod
use cellmesh3d_adtree_mod
use cell_mesh3d_octree_mod
use cellmesh3d_geometry_mod
contains 


!Volume mesh construction subroutine ('variable' surface format) ===========================
subroutine construct_mesh(volume_mesh,ot_mesh,surface_mesh,surface_adtree,cm3dopt,cm3dfailure)
implicit none 

!Variables - Import
integer(in) :: cm3dfailure
type(cm3d_options) :: cm3dopt
type(octree_data) :: ot_mesh
type(surface_data) :: surface_mesh
type(vol_mesh_data) :: volume_mesh
type(tree_data) :: surface_adtree

!Variables - Local 
integer(in) :: type_tgt
type(vol_mesh_data) :: volume_mesh_full

!Set target edge/face type for mesh 
if (cm3dopt%meshinout == 'in') then !mesh internal 
    type_tgt = 3
elseif (cm3dopt%meshinout == 'out') then !mesh external 
    type_tgt = 2
end if 

!Display
if (cm3dopt%dispt == 1) then
    write(*,'(A)') '--> constructing base complete volume mesh '
end if

!Construct initial full mesh 
call build_full_mesh(volume_mesh_full,ot_mesh,cm3dopt)
call build_vmesh_edges(volume_mesh_full%nedge,volume_mesh_full%edges,volume_mesh_full,6_in)

!Intersect the full volume mesh with the surface geometry to construct the final volume mesh 
call clip_vmesh2surface(volume_mesh,volume_mesh_full,surface_adtree,surface_mesh,cm3dopt,type_tgt,cm3dfailure)
return 
end subroutine construct_mesh




!Clip volume mesh to geometry subroutine ===========================
subroutine clip_vmesh2surface(volume_mesh,volume_mesh_full,surface_adtree,surface_mesh,cm3dopt,type_tgt,cm3dfailure)
implicit none 

!Variables - Import
integer(in) :: type_tgt,cm3dfailure
type(cm3d_options) :: cm3dopt
type(surface_data) :: surface_mesh
type(vol_mesh_data) :: volume_mesh_full,volume_mesh
type(tree_data) :: surface_adtree

!Variables - Local 
integer(in) :: ii,ee,ff,vv 
integer(in) :: vtx_idx_smesh,etgt,ev1,ev2
real(dp) :: ef
real(dp) :: vbc(3)
integer(in), dimension(:), allocatable :: vtx_external
type(edgeint), dimension(:), allocatable :: edge_clip_vm,edge_clip_sm
type(triint), dimension(:), allocatable :: tri_clip_sm
type(vtx_intersect), dimension(:), allocatable :: vtx_clip_sm

!Build volume mesh V2E 
allocate(volume_mesh_full%V2E(volume_mesh_full%nvtx,6))
allocate(volume_mesh_full%valence(volume_mesh_full%nvtx))
volume_mesh_full%V2E(:,:) = 0 
volume_mesh_full%valence(:) = 0 
do ee=1,volume_mesh_full%nedge
    ev1 = volume_mesh_full%edges(ee,1)
    ev2 = volume_mesh_full%edges(ee,2)
    volume_mesh_full%valence(ev1) = volume_mesh_full%valence(ev1) + 1
    volume_mesh_full%valence(ev2) = volume_mesh_full%valence(ev2) + 1
    if ((volume_mesh_full%valence(ev1) .GT. 6) .OR. (volume_mesh_full%valence(ev2) .GT. 6)) then 
        print *, '** octree valence 6 exceeded on edge : ',ev1,ev2
        if (volume_mesh_full%valence(ev1) .GT. 6) then 
            print *, '** vertex excessive valence : ',ev1
            do ii=1,6
                etgt = volume_mesh_full%V2E(ev1,ii)
                print *, volume_mesh_full%edges(etgt,1),volume_mesh_full%edges(etgt,2)
            end do 
        end if 
        if (volume_mesh_full%valence(ev2) .GT. 6) then 
            print *, '** vertex excessive valence : ',ev2
            do ii=1,6
                etgt = volume_mesh_full%V2E(ev2,ii)
                print *, volume_mesh_full%edges(etgt,1),volume_mesh_full%edges(etgt,2)
            end do 
        end if 
        cm3dfailure = 1
        return 
    end if 
    volume_mesh_full%V2E(ev1,volume_mesh_full%valence(ev1)) = ee
    volume_mesh_full%V2E(ev2,volume_mesh_full%valence(ev2)) = ee 
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

!Initialise clipped volume structure 
allocate(edge_clip_vm(volume_mesh_full%nedge))
do ee=1,volume_mesh_full%nedge
    edge_clip_vm(ee)%nitem = 0 
    edge_clip_vm(ee)%nint = 0 
    edge_clip_vm(ee)%type = 0 
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
    allocate(edge_clip_sm(ee)%int_tri_loc(cm3dopt%NintEmax))
    allocate(edge_clip_sm(ee)%surfseg(cm3dopt%NintEmax))
    edge_clip_sm(ee)%vtx_idx(:) = 0 
    edge_clip_sm(ee)%vmface(:,:) = 0 
    edge_clip_sm(ee)%nfint(:) = 0 
    edge_clip_sm(ee)%surfseg(:) = 0 
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
    write(*,'(A)') '--> perturbing surface and mesh co-incident vertices '
end if

!Perturb volume mesh vertices within tollerance of surfaces 
call perturb_onsurf_vmesh_vertices(volume_mesh_full,surface_adtree,surface_mesh,cm3dopt)

!Perturb surface mesh vertices within tollerance of volume mesh faces 
call perturb_onmesh_smesh_vertices(volume_mesh_full,surface_adtree,surface_mesh,cm3dopt)

!Display
if (cm3dopt%dispt == 1) then
    write(*,'(A)') '--> intersecting surface mesh with the base volume mesh '
end if

!Intersect the surface mesh edges with the volume mesh faces
if (cm3dopt%dispt == 1) then
    write(*,'(A)') '    {intersecting surface edges with the volume mesh}'
end if
call clip_surfedges2volfaces(edge_clip_sm,vtx_clip_sm,edge_clip_vm,volume_mesh_full,surface_mesh,surface_adtree,&
                             cm3dopt,vtx_idx_smesh)

!Intersect volume mesh edges with the surface mesh faces
if (cm3dopt%dispt == 1) then
    write(*,'(A)') '    {intersecting volume edges with the surface mesh}'
end if                             
call clip_voledges2surffaces(edge_clip_vm,edge_clip_sm,vtx_clip_sm,tri_clip_sm,volume_mesh_full,surface_adtree,surface_mesh,&
                             cm3dopt,vtx_idx_smesh)

!Construct edge meshes 
if (cm3dopt%dispt == 1) then
    write(*,'(A)') '    {constructing intersecting edge internal meshes}'
end if    
call build_vol_intersect_edge_mesh(edge_clip_vm,volume_mesh_full,cm3dopt)
call build_surf_intersect_edge_mesh(edge_clip_sm,surface_mesh,cm3dopt)

!Construct vtx_external and classify edge types 
if (cm3dopt%dispt == 1) then
    write(*,'(A)') '--> evaluating edge and vertex containment '
end if
call classify_edge_vertex_geom_containment(vtx_external,edge_clip_vm,volume_mesh_full,surface_mesh,cm3dopt,cm3dfailure)

!Assign the volume mesh cells of each surface mesh vertex
if (cm3dopt%dispt == 1) then
    write(*,'(A)') '    {assigning surface mesh vertex containing volume mesh cells}'
end if
call assign_smvtx_vmcells(edge_clip_sm,vtx_clip_sm,surface_mesh,volume_mesh_full)

!Display maximum number of intersections 
if (cm3dopt%dispt == 1) then
    write(*,'(A,I0,A)') '    {maximum intersections on a volume edge = ',maxval(edge_clip_vm(:)%nint),'}'
    write(*,'(A,I0,A)') '    {maximum intersections on a surface vertex = ',maxval(vtx_clip_sm(:)%nintersect),'}'
    write(*,'(A,I0,A)') '    {maximum intersections on a surface edge = ',maxval(edge_clip_sm(:)%nint),'}'
    write(*,'(A,I0,A)') '    {maximum intersections on a surface triangle = ',maxval(tri_clip_sm(:)%nvtx),'}'
end if

!Set maximum intersection count on a surface triangle 
cm3dopt%stnintmax = 3*(maxval(vtx_clip_sm(:)%nintersect) + maxval(edge_clip_sm(:)%nint)) + maxval(tri_clip_sm(:)%nvtx)

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

!Debug export external vertices =======================
! open(11,file='io/vtxexternal.dat')
! do ee=1,volume_mesh_full%nvtx
!     if (vtx_external(ee) == 1) then 
!         write(11,*) volume_mesh_full%vtx(ee,:)
!     end if
! end do 
! close(11)
!Debug export external vertices =======================

!Debug export intersecting vm edges =======================
! open(11,file='io/vmeint.dat')
! do ee=1,volume_mesh_full%nedge
!     if ((edge_clip_vm(ee)%nint .GT. 0) .AND. (edge_clip_vm(ee)%type .NE. 4)) then 
!         print *, 'ret'
!     end if 

!     if (edge_clip_vm(ee)%type == 4) then 
!         write(11,*) volume_mesh_full%vtx(volume_mesh_full%edges(ee,1),:),volume_mesh_full%vtx(volume_mesh_full%edges(ee,2),:)
!     end if
! end do 
! close(11)
!Debug export intersecting vm edges =======================

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
            vbc = get_barycentric_coordinates(tri_clip_sm(ff)%intloc(ii,:),&
            surface_mesh%vertices(surface_mesh%connectivity(ff,1),:),&
            surface_mesh%vertices(surface_mesh%connectivity(ff,2),:),&
            surface_mesh%vertices(surface_mesh%connectivity(ff,3),:))

            !Assign interpolated surface curvature
            surface_mesh%vtx_rcurv_full(tri_clip_sm(ff)%vtx_idx(ii)) = &
            surface_mesh%vtx_rcurv(surface_mesh%connectivity(ff,1))*vbc(1) +&
            surface_mesh%vtx_rcurv(surface_mesh%connectivity(ff,2))*vbc(2) +&
            surface_mesh%vtx_rcurv(surface_mesh%connectivity(ff,3))*vbc(3)
        end do 
    end if 
end do 

!Build surface triangle internal edges
if (cm3dopt%dispt == 1) then
    write(*,'(A)') '--> constructing volume mesh clipped surface mesh edges '
end if
call build_surftri_internal_edges(tri_clip_sm,edge_clip_sm,vtx_clip_sm,surface_mesh,volume_mesh_full,cm3dopt)

!Build surface mesh clippled volume mesh faces 
if (cm3dopt%dispt == 1) then
    write(*,'(A)') '--> constructing surface mesh clipped volume mesh faces '
end if
call build_surface_clipped_vmesh_faces(volume_mesh_full,surface_mesh,surface_adtree,cm3dopt,tri_clip_sm,vtx_clip_sm,&
                                       edge_clip_sm,edge_clip_vm,vtx_idx_smesh,type_tgt,surface_mesh%vertices_full)

!Display
if (cm3dopt%dispt == 1) then
    write(*,'(A)') '--> constructing volume mesh clipped surface mesh faces '
end if 

!Build volume mesh clipped surface triangles 
call build_volume_clipped_smesh_faces(surface_mesh,tri_clip_sm,edge_clip_sm,vtx_clip_sm,volume_mesh_full,cm3dopt)

!Display
if (cm3dopt%dispt == 1) then
    write(*,'(A)') '--> constructing final clipped volume mesh '
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
integer(in) :: vins,face_idx,fcount,vtgt,fidx,flparent
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
    do aa=1,volume_mesh_full%faces(ff)%nfclipped !construct faces

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
        volume_mesh%faces(face_idx)%cleft = volume_mesh_full%faces(ff)%cright
        volume_mesh%faces(face_idx)%cright = volume_mesh_full%faces(ff)%cleft

        !Store final index of this face
        volume_mesh_full%faces(ff)%face_clipped(aa)%fidx_final = face_idx
    end do 
    do aa=1,volume_mesh_full%faces(ff)%nfclipped !set face properties 

        !Index of this face 
        fidx = volume_mesh_full%faces(ff)%face_clipped(aa)%fidx_final

        !Set properties 
        if (volume_mesh_full%faces(ff)%face_clipped(aa)%ftype == 'int') then !set as sub face to tag during cell bisection 

            !Local clipping parent index 
            flparent = volume_mesh_full%faces(ff)%face_clipped(aa)%fparent

            !Set properties 
            volume_mesh%faces(fidx)%ftype = 'int'
            volume_mesh%faces(fidx)%fparent = volume_mesh_full%faces(ff)%face_clipped(flparent)%fidx_final 
        else !set as normal face
            volume_mesh%faces(fidx)%ftype = 'nrm'
            volume_mesh%faces(fidx)%fparent = 0 
        end if 
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
subroutine build_volume_clipped_smesh_faces(surface_mesh,tri_clip_sm,edge_clip_sm,vtx_clip_sm,volume_mesh_full,cm3dopt)
implicit none 

!Variables - Import
type(surface_data) :: surface_mesh
type(vol_mesh_data) :: volume_mesh_full
type(triint), dimension(:) :: tri_clip_sm
type(edgeint), dimension(:) :: edge_clip_sm
type(vtx_intersect), dimension(:) :: vtx_clip_sm
type(cm3d_options) :: cm3dopt

!Variables - Local 
integer(in) :: ff,vv,ee,aa,ii,cc,ss 
integer(in) :: vintvalid1,vintvalid2,ebase,vbase,enew,vnew,all_te_subedge,diff_parent,nint_tri
integer(in) :: nc_ontri,cl,cr,ftgt,ctgt,etgt,vtgt,nedge_tface,v1,v2,face_nvtx,vc,vn,v0
integer(in) :: cell_on_tri(cm3dopt%max_int_size),edge_face_idx(2*cm3dopt%max_int_size),edges_on_tri(2*cm3dopt%max_int_size,4)
integer(in) :: cell_selected(volume_mesh_full%ncell),tedge_cell(3),eend_cell(3,2)
integer(in) :: vtx_face_pos(surface_mesh%nvtxf),fliparray(cm3dopt%max_int_size)
integer(in) :: v2e_flocal(surface_mesh%nvtxf,2)
real(dp) :: Nft(3),Nf(3),vtxC(3),vtxN(3),fmid(3)

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
        vtgt = surface_mesh%connectivity(ff,ee)
        nint_tri = nint_tri + vtx_clip_sm(vtgt)%nintersect
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

        !Assign a cell containement for and edge of this triangle that has no face intersections 
        tedge_cell(:) = 0 
        do ee=1,3
            etgt = surface_mesh%F2E(ff,ee)
            if (edge_clip_sm(etgt)%nint == 0) then !Pick cell that contains the edge midpoint 

                !Edge ends 
                v1 = surface_mesh%edges(etgt,1)
                v2 = surface_mesh%edges(etgt,2)

                !If they are both within the same cell by vertex assignment take that cell for this edge 
                if (surface_mesh%vtx_vmesh_cell(v1) == surface_mesh%vtx_vmesh_cell(v2)) then !both clearly within the cell 
                    tedge_cell(ee) = surface_mesh%vtx_vmesh_cell(v1)
                else !one or both ends lies on the cell boundary 
                    if (surface_mesh%vtx_vmesh_cell(v1) .NE. 0) then     
                        tedge_cell(ee) = surface_mesh%vtx_vmesh_cell(v1)
                    elseif (surface_mesh%vtx_vmesh_cell(v2) .NE. 0) then     
                        tedge_cell(ee) = surface_mesh%vtx_vmesh_cell(v2)
                    else
                        print *, '** ambiguous cell internal surface mesh edge assignment at surface face : ',ff
                        print *, surface_mesh%vtx_vmesh_cell(v1),surface_mesh%vtx_vmesh_cell(v2)
                    end if 
                end if 
            end if 
        end do 

        !Assign triangle edge end cell containements 
        eend_cell(:,:) = 0 
        do ee=1,3 !Triangle edges 
            etgt = surface_mesh%F2E(ff,ee)
            v1 = surface_mesh%edges(etgt,1)
            v2 = surface_mesh%edges(etgt,2)
            eend_cell(ee,1) = surface_mesh%vtx_vmesh_cell(v1)
            eend_cell(ee,2) = surface_mesh%vtx_vmesh_cell(v2)
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
                    if (ff == 709) then 
                        stop 
                    end if 
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
                do ee=1,nedge_tface*2

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
                                             edge_clip_sm,edge_clip_vm,vtx_idx_smesh,type_tgt,smesh_vtx_full)
implicit none 

!Variables - Import
integer(in) :: vtx_idx_smesh,type_tgt
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
integer(in) :: nselected,ttgt,ftgt,nteselect,vs1,vs2,ve1,ve2,nedge_face,etgt,etype,nface_clip,e0,is_clipped,vlnceflag
integer(in) :: ebase,vbase,v0face,enew,vnew,vtgt,vc,vn,edg_dup,nfaceb,vintvalid1,vintvalid2,vtx_in_base_face,fsubparent
integer(in) :: node_select(surface_adtree%nnode),edge_face_idx(surface_mesh%nfcs)
integer(in) :: face_nvtx(cm3dopt%max_int_size),vmface_is_intersected(volume_mesh_full%nface)
integer(in) :: te_select(surface_mesh%nfcs,2),face_ledge(surface_mesh%nfcs,3)
integer(in) :: v2e_flocal(-vtx_idx_smesh:volume_mesh_full%nvtx,4),valence_flocal(-vtx_idx_smesh:volume_mesh_full%nvtx)
integer(in) :: vtx_face_idx(-vtx_idx_smesh:volume_mesh_full%nvtx,4)
integer(in) :: vtx_face_pos(-vtx_idx_smesh:volume_mesh_full%nvtx,4)
real(dp) :: cpadSZ,zxmin,zxmax,zymin,zymax,zzmin,zzmax
real(dp) :: Nf(3),Nft(3),vtxC(3),vtxN(3)

!Tag all volume mesh faces that have any form of intersection
vmface_is_intersected(:) = 0 
do ff=1,volume_mesh_full%nface
    is_clipped = 0 
    do ee=1,volume_mesh_full%faces(ff)%nvtx
        etgt = volume_mesh_full%faces(ff)%edges(ee)
        if (edge_clip_vm(etgt)%nint .GT. 0) then 
            is_clipped = 1
            exit 
        end if 
    end do 
    if (is_clipped == 1) then 
        vmface_is_intersected(ff) = 1
    end if 
end do 
do ee=1,surface_mesh%nedge
    if (edge_clip_sm(ee)%nint .GT. 0) then 
        do ii=1,edge_clip_sm(ee)%nint
            do aa=1,edge_clip_sm(ee)%nfint(ii)
                ftgt = edge_clip_sm(ee)%vmface(ii,aa)
                vmface_is_intersected(ftgt) = 1
            end do 
        end do 
    end if 
end do 
do vv=1,surface_mesh%nvtx
    do ii=1,vtx_clip_sm(vv)%nintersect
        ftgt = vtx_clip_sm(vv)%face_int_idx(ii)
        vmface_is_intersected(ftgt) = 1
    end do
end do 

!Clip each face
nselected = 0 
nteselect = 0 
nedge_face = 0
face_nvtx(:) = 0  
node_select(:) = 0 
vtx_face_idx(:,:) = 0 
vtx_face_pos(:,:) = 0
edge_face_idx(:) = 0 
te_select(:,:) = 0 
face_ledge(:,:) = 0 
v2e_flocal(:,:) = 0 
valence_flocal(:) = 0 
do ff=1,volume_mesh_full%nface

    !If this face has any edge intersections with the surface mesh then clip 
    is_clipped = 0 
    do ee=1,volume_mesh_full%faces(ff)%nvtx
        etgt = volume_mesh_full%faces(ff)%edges(ee)
        if (edge_clip_vm(etgt)%nint .GT. 0) then 
            is_clipped = 1
            exit 
        end if 
    end do 

    !If this face has no intersections on its edges but surface mesh faces intersect this face then set as internally clipped 
    !(this must construct the face(s) formed by the clipped patch with the opposite orientation to the parent face to subtract from it)
    if ((is_clipped == 0) .AND. (vmface_is_intersected(ff) == 1)) then 
        if (cm3dopt%dispt == 1) then
            write(*,'(A,I0,A)') '    {face internally clipped : ',ff,'}'
        end if
        is_clipped = 2
    end if 

    !Process face if clipped
    if (is_clipped .NE. 0) then !If any edge on the face is clipped by the surface geometry

        !Build face normal vector 
        Nf = newell_normal(volume_mesh_full%faces(ff)%nvtx,volume_mesh_full%faces(ff)%vertices,volume_mesh_full%vtx)

        !Padding size 
        cpadSZ = norm2(Nf)*0.005d0

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
                        face_ledge(nedge_face,3) = etgt
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
                                face_ledge(nedge_face,3) = etgt
                            end if 
                        end do 
                    end if 
                end do 

                !Add triangle internal edge segments 
                nfaceb = nedge_face
                do ee=1,nteselect

                    !Check if dulicate if any selected base volume mesh edges or edge segments are selected 
                    edg_dup = 0 
                    if (nfaceb .GE. 1) then 
                        do ii=nfaceb,nedge_face
                            if ((face_ledge(ii,1) == -1*te_select(ee,1)) .AND. (face_ledge(ii,2) == -1*te_select(ee,2))) then 
                                edg_dup = 1
                                exit 
                            elseif ((face_ledge(ii,2) == -1*te_select(ee,1)) .AND. (face_ledge(ii,1) == -1*te_select(ee,2))) then 
                                edg_dup = 1
                                exit 
                            end if 
                        end do 
                    end if 

                    !Add if new
                    if (edg_dup == 0) then 
                        nedge_face = nedge_face + 1 
                        face_ledge(nedge_face,1) = -1*te_select(ee,1)
                        face_ledge(nedge_face,2) = -1*te_select(ee,2)
                        face_ledge(nedge_face,3) = 0
                    else 
                        !print *, 'edupe ',ff
                    end if 
                end do 

                !Initialise arrays for this face 
                do ee=1,nedge_face
                    valence_flocal(face_ledge(ee,1)) = 0
                    valence_flocal(face_ledge(ee,2)) = 0  
                    v2e_flocal(face_ledge(ee,1),:) = 0
                    v2e_flocal(face_ledge(ee,2),:) = 0  
                    vtx_face_pos(face_ledge(ee,1),:) = 0
                    vtx_face_pos(face_ledge(ee,2),:) = 0  
                end do 

                !Build v2e for the clipped face mesh and tag high valence if present 
                vlnceflag = 0 
                do ee=1,nedge_face
                    valence_flocal(face_ledge(ee,1)) = valence_flocal(face_ledge(ee,1)) + 1
                    v2e_flocal(face_ledge(ee,1),valence_flocal(face_ledge(ee,1))) = ee 
                    valence_flocal(face_ledge(ee,2)) = valence_flocal(face_ledge(ee,2)) + 1
                    v2e_flocal(face_ledge(ee,2),valence_flocal(face_ledge(ee,2))) = ee 
                end do 
                if ((maxval(valence_flocal(face_ledge(1:nedge_face,1))) .GT. 2) .OR. &
                    (maxval(valence_flocal(face_ledge(1:nedge_face,2))) .GT. 2)) then !if maxvalence is greater than two tag as a high valence case
                    vlnceflag = 1
                    ! print *, '=== maxvalence e1 / e2 = ',maxval(valence_flocal(face_ledge(1:nedge_face,1))),&
                    !                                  maxval(valence_flocal(face_ledge(1:nedge_face,2))) 
                    ! vlnceflag = face_ledge(maxloc(valence_flocal(face_ledge(1:nedge_face,1)),1),1)
                    ! print *,'vh 1 = ', vlnceflag
                    ! if (vlnceflag .GT. 0) then
                    !     print *, volume_mesh_full%vtx(vlnceflag,:)
                    ! else
                    !     print *, smesh_vtx_full(abs(vlnceflag),:)
                    ! end if 
                    ! vlnceflag = face_ledge(maxloc(valence_flocal(face_ledge(1:nedge_face,2)),1),2)
                    ! print *,'vh 2 = ', vlnceflag
                    ! if (vlnceflag .GT. 0) then
                    !     print *, volume_mesh_full%vtx(vlnceflag,:)
                    ! else
                    !     print *, smesh_vtx_full(abs(vlnceflag),:)
                    ! end if 
                    ! vlnceflag = 1
                end if 

                !Debug print (if non contiguous face)=============================
                if ((minval(v2e_flocal(face_ledge(1:nedge_face,1),1)) == 0) .OR. &
                    (minval(v2e_flocal(face_ledge(1:nedge_face,1),2)) == 0) .OR. &
                    (minval(v2e_flocal(face_ledge(1:nedge_face,2),1)) == 0) .OR. &
                    (minval(v2e_flocal(face_ledge(1:nedge_face,2),2)) == 0)) then 
                ! if (ff == 1156) then 
                    print *, 'Zero V2E Face = ', ff,' ====================='
                    print *, 'V2E zero at: '
                    do ee=1,nedge_face
                        if ((v2e_flocal(face_ledge(ee,1),1) == 0) .OR. (v2e_flocal(face_ledge(ee,1),2) == 0 )) then 
                            print *, 'vtx = ',face_ledge(ee,1),' edge = ',face_ledge(ee,3)
                            if (face_ledge(ee,1) .GT. 0) then 
                                Nf(:) = volume_mesh_full%vtx(face_ledge(ee,1),:)
                            else
                                Nf(:) = smesh_vtx_full(abs(face_ledge(ee,1)),:) !volume_mesh_full%lsurface(ctgt)%vertices(abs(face_ledge(ee,1)),:)
                            end if 
                            print *, Nf(:)
                        end if 
                        if ((v2e_flocal(face_ledge(ee,2),1) == 0) .OR. (v2e_flocal(face_ledge(ee,2),2) == 0 )) then 
                            print *, 'vtx = ',face_ledge(ee,2),' edge = ',face_ledge(ee,3)
                            if (face_ledge(ee,2) .GT. 0) then 
                                Nf(:) = volume_mesh_full%vtx(face_ledge(ee,2),:)
                            else
                                Nf(:) = smesh_vtx_full(abs(face_ledge(ee,2)),:) !volume_mesh_full%lsurface(ctgt)%vertices(abs(face_ledge(ee,1)),:)
                            end if 
                            print *, Nf(:)
                        end if 
                    end do 
                    print *, 'face edges'
                    do ee=1,volume_mesh_full%faces(ff)%nvtx
                        etgt = volume_mesh_full%faces(ff)%edges(ee)
                        print *, '---- edge :',ee,' / ',etgt,'(',edge_clip_vm(etgt)%nint,') type = ', edge_clip_vm(etgt)%type
                        print *, 'v12 = ',volume_mesh_full%edges(etgt,1),volume_mesh_full%edges(etgt,2)
                        if (edge_clip_vm(etgt)%nint .NE. 0) then
                            do nn=1,edge_clip_vm(etgt)%nint
                                ! print *, edge_clip_vm(etgt)%inttype(nn)
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
                    stop
                end if

                !Construct set of full faces from the clipped face mesh 
                nface_clip = 0 
                face_nvtx(:) = 0
                edge_face_idx(1:nedge_face) = 0 
                do ii=1,nedge_face

                    !Find starting edge 
                    e0 = 0 
                    if (vlnceflag == 0) then !all vertices in the clipping set have valence = 2 so any edge is a valid start point 
                        do ee=1,nedge_face
                            if (edge_face_idx(ee) == 0) then 
                                e0 = ee 
                                exit 
                            end if
                        end do 
                    else !some vertices have valence > 2 so pick an edge containting one of these to start from if any are available 

                        !Find and select a high valence ended edge vlnceflag
                        do ee=1,nedge_face
                            if (edge_face_idx(ee) == 0) then 
                                if ((valence_flocal(face_ledge(ee,1)) .GT. 2) .OR. (valence_flocal(face_ledge(ee,2)) .GT. 2)) then 
                                    e0 = ee 
                                    exit 
                                end if 
                            end if
                        end do 

                        !If no high valence ended edges are avalable then pick a normal edge and/or set the valence flag to zero 
                        if (e0 == 0) then 

                            !Set flag to zero 
                            vlnceflag = 0 

                            !Select edge if any 
                            do ee=1,nedge_face
                                if (edge_face_idx(ee) == 0) then 
                                    e0 = ee 
                                    exit 
                                end if
                            end do 
                        end if 
                    end if 

                    !Exit if no new edges located
                    if (e0 == 0) then 
                        exit 
                    end if 

                    !Increment face count 
                    nface_clip = nface_clip + 1

                    !Set initial parameters 
                    ebase = e0
                    if (vlnceflag == 0) then !set to end 2 for the correct direction 
                        vbase = face_ledge(ebase,2) 
                        v0face = face_ledge(ebase,1)
                    else !pick an end without a high valence vertex
                        if (valence_flocal(face_ledge(ebase,2)) .GT. 2) then 
                            vbase = face_ledge(ebase,1) 
                            v0face = face_ledge(ebase,2)
                        else
                            vbase = face_ledge(ebase,2) 
                            v0face = face_ledge(ebase,1)
                        end if 
                    end if 
                    edge_face_idx(ebase) = nface_clip 

                    !Initialise this curve 
                    do aa=1,4 
                        if (vtx_face_pos(vbase,aa) == 0) then 
                            vtx_face_pos(vbase,aa) = 1
                            vtx_face_idx(vbase,aa) = nface_clip
                            exit
                        end if 
                    end do 
                    face_nvtx(nface_clip) = 1

                    !Mesh face 
                    do ee=1,nedge_face

                        !Find next edge
                        enew = 0 
                        do aa=1,valence_flocal(vbase)
                            if (edge_face_idx(v2e_flocal(vbase,aa)) == 0) then 
                                enew = v2e_flocal(vbase,aa)
                                exit 
                            end if  
                        end do 
                        
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
                        do aa=1,4 
                            if (vtx_face_pos(vnew,aa) == 0) then 
                                vtx_face_pos(vnew,aa) = face_nvtx(nface_clip)
                                vtx_face_idx(vnew,aa) = nface_clip
                                exit
                            end if 
                        end do 

                        !Add edge to this face index 
                        edge_face_idx(enew) = nface_clip

                        !Update base location
                        ebase = enew 
                        vbase = vnew 

                        !If the new base vertex is the starting vertex then exit as the face is closed 
                        if (vbase == v0face) then 
                            exit 
                        end if 

                        !If the base vertex is high valence and the tag is still active then exit 
                        if (vlnceflag == 1) then 
                            if (valence_flocal(vbase) .GT. 2) then 
                                exit 
                            end if 
                        end if 
                    end do
                    !print *, face_nvtx(nface_clip) 
                end do 
                !print *, 'nface_clip = ',nface_clip

                !Debug =============================
                ! if (nface_clip .GT. 1) then 
                !     print *, 'fbase = ',ff,' - Nface = ',nface_clip
                !     print *, 'Face = ', ff
                !     ! do ee=1,4
                !     !     print *, volume_mesh_full%vtx(volume_mesh_full%faces(ff)%vertices(ee),:)
                !     ! end do 
    
                !     print *, '----------- nedge_face = ',nedge_face 
                !     do ee=1,nedge_face
                !         print '(A,I0,A,I0,A,I0)', ' edge : ',ee,' ',face_ledge(ee,1),' ',face_ledge(ee,2)
                !         print '(I0,A,I0,A,I0,A,I0)', v2e_flocal(face_ledge(ee,1),1),' ',v2e_flocal(face_ledge(ee,1),2),&
                !         ' / ',v2e_flocal(face_ledge(ee,2),1),' ',v2e_flocal(face_ledge(ee,2),2)
                !     end do 

                !     open(11,file='io/face_cvtx.dat') 
                !     do ee=1,nedge_face
                !         if (edge_face_idx(ee) == 3) then 
                !             if (face_ledge(ee,1) .GT. 0) then 
                !                 Nf(:) = volume_mesh_full%vtx(face_ledge(ee,1),:)
                !             else
                !                 Nf(:) = smesh_vtx_full(abs(face_ledge(ee,1)),:) !volume_mesh_full%lsurface(ctgt)%vertices(abs(face_ledge(ee,1)),:)
                !             end if 
                !             if (face_ledge(ee,2) .GT. 0) then 
                !                 Nft(:) = volume_mesh_full%vtx(face_ledge(ee,2),:)
                !             else
                !                 Nft(:) = smesh_vtx_full(abs(face_ledge(ee,2)),:)!volume_mesh_full%lsurface(ctgt)%vertices(abs(face_ledge(ee,2)),:)
                !             end if 
                !             write(11,*) Nf(:),Nft(:)
                !         end if 
                !     end do 
                !     close(11)

                !     ! print *,'Face : ',ff
                !     ! do vv=1,volume_mesh_full%faces(ff)%face_clipped(ii)%nvtx
                !     !     print *, volume_mesh_full%faces(ff)%face_clipped(ii)%vertices(vv),edge_face_idx(vv)
                !     ! end do 
                ! end if 

                !Build full clipped faces 
                volume_mesh_full%faces(ff)%nfclipped = nface_clip
                allocate(volume_mesh_full%faces(ff)%face_clipped(nface_clip))
                do ii=1,nface_clip

                    !Set properties
                    volume_mesh_full%faces(ff)%face_clipped(ii)%ftype = 'clp'
                    volume_mesh_full%faces(ff)%face_clipped(ii)%fparent = ff
                    volume_mesh_full%faces(ff)%face_clipped(ii)%nvtx = face_nvtx(ii)
                    allocate(volume_mesh_full%faces(ff)%face_clipped(ii)%vertices(face_nvtx(ii)))
                    volume_mesh_full%faces(ff)%face_clipped(ii)%vertices(:) = 0 

                    !Accumulate face 
                    do ee=1,nedge_face
                        vtgt = face_ledge(ee,1) 
                        do aa=1,4   
                            if (vtx_face_idx(vtgt,aa) == ii) then 
                                volume_mesh_full%faces(ff)%face_clipped(ii)%vertices(vtx_face_pos(vtgt,aa)) = vtgt
                                exit 
                            end if 
                        end do 
                        vtgt = face_ledge(ee,2) 
                        do aa=1,4   
                            if (vtx_face_idx(vtgt,aa) == ii) then 
                                volume_mesh_full%faces(ff)%face_clipped(ii)%vertices(vtx_face_pos(vtgt,aa)) = vtgt
                                exit 
                            end if 
                        end do 
                    end do 

                    !Check and set clipped face orientation ------
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

                    !Flip if required 
                    if (dot_product(Nf,Nft) .GT. 0.0d0) then 
                        volume_mesh_full%faces(ff)%face_clipped(ii)%vertices(:) = &
                        flip_face(volume_mesh_full%faces(ff)%face_clipped(ii)%nvtx,&
                        volume_mesh_full%faces(ff)%face_clipped(ii)%vertices(:))
                    end if 
                end do 

                !If the face is internally clipped and a face is formed from the base volume mesh edges of this face (so at least two faces have been built)
                !flip the direction of all internally clipped faces to subtract from the main face 
                !and set the required properties to avoid bisecting this cell by not tagging the subtracted face 
                if ((is_clipped == 2) .AND. (nface_clip .GE. 2)) then 

                    !Find the face in the clip structure with all vertices from the base face 
                    fsubparent = 0 
                    do ii=1,nface_clip
                        vtx_in_base_face = 0
                        do vv=1,volume_mesh_full%faces(ff)%face_clipped(ii)%nvtx
                            vtgt = volume_mesh_full%faces(ff)%face_clipped(ii)%vertices(vv)
                            do aa=1,volume_mesh_full%faces(ff)%nvtx
                                if (volume_mesh_full%faces(ff)%vertices(aa) == vtgt) then 
                                    vtx_in_base_face = vtx_in_base_face + 1
                                end if 
                            end do 
                        end do
                        if ((vtx_in_base_face == volume_mesh_full%faces(ff)%nvtx) .AND. &
                        (volume_mesh_full%faces(ff)%face_clipped(ii)%nvtx == volume_mesh_full%faces(ff)%nvtx)) then 
                            fsubparent = ii 
                            exit 
                        end if 
                    end do 

                    !Set the properties of all the internal faces
                    do ii=1,nface_clip

                        !Check if face contains any vertices from the base volume mesh face 
                        vtx_in_base_face = 0
                        do vv=1,volume_mesh_full%faces(ff)%face_clipped(ii)%nvtx
                            vtgt = volume_mesh_full%faces(ff)%face_clipped(ii)%vertices(vv)
                            do aa=1,volume_mesh_full%faces(ff)%nvtx
                                if (volume_mesh_full%faces(ff)%vertices(aa) == vtgt) then 
                                    vtx_in_base_face = 1
                                    exit 
                                end if 
                            end do 
                            if (vtx_in_base_face == 1) then 
                                exit 
                            end if 
                        end do

                        !If this face has no vertices from the base face
                        if (vtx_in_base_face == 0) then 

                            !Flip face
                            volume_mesh_full%faces(ff)%face_clipped(ii)%vertices(:) = &
                            flip_face(volume_mesh_full%faces(ff)%face_clipped(ii)%nvtx,&
                            volume_mesh_full%faces(ff)%face_clipped(ii)%vertices(:))

                            !Tag as internally subtracted face and tag parent as the correct clipped face in this structure 
                            volume_mesh_full%faces(ff)%face_clipped(ii)%ftype = 'int'
                            volume_mesh_full%faces(ff)%face_clipped(ii)%fparent = fsubparent
                        end if 
                    end do 
                end if 
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
                    ! volume_mesh_full%faces(ff)%face_clipped(1)%vertices(:) = volume_mesh_full%faces(ff)%vertices(:)
                    volume_mesh_full%faces(ff)%face_clipped(1)%vertices(:) = &
                    flip_face(volume_mesh_full%faces(ff)%nvtx,volume_mesh_full%faces(ff)%vertices(:))
                    volume_mesh_full%faces(ff)%face_clipped(1)%ftype = 'ret'
                    volume_mesh_full%faces(ff)%face_clipped(1)%fparent = ff
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
                ! volume_mesh_full%faces(ff)%face_clipped(1)%vertices(:) = volume_mesh_full%faces(ff)%vertices(:)
                volume_mesh_full%faces(ff)%face_clipped(1)%vertices(:) = &
                flip_face(volume_mesh_full%faces(ff)%nvtx,volume_mesh_full%faces(ff)%vertices(:))
                volume_mesh_full%faces(ff)%face_clipped(1)%ftype = 'ret'
                volume_mesh_full%faces(ff)%face_clipped(1)%fparent = ff
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
            ! volume_mesh_full%faces(ff)%face_clipped(1)%vertices(:) = volume_mesh_full%faces(ff)%vertices(:)
            volume_mesh_full%faces(ff)%face_clipped(1)%vertices(:) = &
            flip_face(volume_mesh_full%faces(ff)%nvtx,volume_mesh_full%faces(ff)%vertices(:))
            volume_mesh_full%faces(ff)%face_clipped(1)%ftype = 'ret'
            volume_mesh_full%faces(ff)%face_clipped(1)%fparent = ff
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
integer(in) :: faces_on_tri(8*cm3dopt%max_int_size),ftselected(volume_mesh_full%nface)
integer(in) :: int_from_face(cm3dopt%stnintmax+1,3)
real(dp) :: intloc(cm3dopt%stnintmax+1,3)

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
        vtgt = surface_mesh%connectivity(ff,ee)
        nint_tri_tot = nint_tri_tot + vtx_clip_sm(vtgt)%nintersect
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
                print *,'surface face vertices = '
                do ii=1,3
                    print *,surface_mesh%vertices(surface_mesh%connectivity(ff,ii),:)
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




!Subroutine to assign volume mesh cells to each surface mesh vertex ===========================
subroutine assign_smvtx_vmcells(edge_clip_sm,vtx_clip_sm,surface_mesh,volume_mesh_full)
implicit none 

!Variables - Import
type(surface_data) :: surface_mesh
type(edgeint), dimension(:) :: edge_clip_sm
type(vtx_intersect), dimension(:) :: vtx_clip_sm
type(vol_mesh_data) :: volume_mesh_full

!Variables - Local
integer(in) :: ee,ff
integer(in) :: ftgt,ev1,ev2,cleft,cright,nupdate,nfi,cell_selected
integer(in) :: cell_from_face(12),cell_select(volume_mesh_full%ncell)
real(dp) :: dotval
real(dp) :: Nf(3),edir(3)

!Assign a cell to each vertex on an intersecting edge 
cell_select(:) = 0 
cell_from_face(:) = 0 
allocate(surface_mesh%vtx_vmesh_cell(surface_mesh%nvtx))
surface_mesh%vtx_vmesh_cell(:) = 0 
do ee=1,surface_mesh%nedge
    if (edge_clip_sm(ee)%nint .GE. 1) then 

        !Vertices 
        ev1 = surface_mesh%edges(ee,1)
        ev2 = surface_mesh%edges(ee,2)

        !Edge direction 
        edir(:) = surface_mesh%vertices(ev2,:) - surface_mesh%vertices(ev1,:)

        !Assign start vertex
        if (surface_mesh%vtx_vmesh_cell(ev1) == 0) then 
            if (vtx_clip_sm(ev1)%nintersect == 0) then !Not an intersecting vertex 
                if (edge_clip_sm(ee)%nfint(1) == 1) then !If a proper intersection with a single face 

                    !Target face 
                    ftgt = edge_clip_sm(ee)%vmface(1,1)

                    !Cells 
                    cleft = volume_mesh_full%faces(ftgt)%cleft
                    cright = volume_mesh_full%faces(ftgt)%cright

                    !Face normal vector 
                    Nf = newell_normal(volume_mesh_full%faces(ftgt)%nvtx,volume_mesh_full%faces(ftgt)%vertices,volume_mesh_full%vtx)

                    !Assign 
                    dotval = dot_product(edir,Nf)
                    if (dotval .GT. 0.0d0) then 
                        surface_mesh%vtx_vmesh_cell(ev1) = cleft
                    elseif (dotval .LT. 0.0d0) then 
                        surface_mesh%vtx_vmesh_cell(ev1) = cright
                    else
                        !print *, '** (fi) zero dot intersection sm edge -> vm face at edge ',ee
                    end if 
                elseif (edge_clip_sm(ee)%nfint(1) .GT. 1) then !If a proper intersection with a multiple faces

                    !Classify cell on all faces 
                    nfi = edge_clip_sm(ee)%nfint(1)
                    do ff=1,edge_clip_sm(ee)%nfint(1)

                        !Target face 
                        ftgt = edge_clip_sm(ee)%vmface(1,ff)

                        !Cells 
                        cleft = volume_mesh_full%faces(ftgt)%cleft
                        cright = volume_mesh_full%faces(ftgt)%cright

                        !Face normal vector 
                        Nf = newell_normal(volume_mesh_full%faces(ftgt)%nvtx,volume_mesh_full%faces(ftgt)%vertices,&
                        volume_mesh_full%vtx)

                        !Assign 
                        dotval = dot_product(edir,Nf)
                        if (dotval .GT. 0.0d0) then 
                            cell_from_face(ff) = cleft
                        elseif (dotval .LT. 0.0d0) then 
                            cell_from_face(ff) = cright
                        else
                            !print *, '** (ei) zero dot intersection sm edge -> vm face at edge ',ee
                            cell_from_face(ff) = 0
                        end if 
                    end do

                    !Select cell with two agreeing faces 
                    cell_selected = 0 
                    do ff=1,nfi 
                        if (cell_from_face(ff) .GT. 0) then 
                            cell_select(cell_from_face(ff)) = cell_select(cell_from_face(ff)) + 1
                        end if 
                    end do 
                    do ff=1,nfi 
                        if (cell_from_face(ff) .GT. 0) then 
                            if (cell_select(cell_from_face(ff)) == 2) then 
                                cell_selected = cell_from_face(ff)
                                exit 
                            end if 
                        end if 
                    end do 

                    !Store 
                    if (cell_selected .GT. 0) then 
                        surface_mesh%vtx_vmesh_cell(ev1) = cell_selected
                    else    
                        print *, '** no cell selected at edge ',ee
                    end if 

                    !Reset selected cell counts 
                    do ff=1,nfi 
                        if (cell_from_face(ff) .GT. 0) then 
                            cell_select(cell_from_face(ff)) = 0
                        end if 
                    end do 
          
                    !Debug print 
                    ! print *, '============ mfint 1',ev1
                    ! print *, 'edge = ',ev1,ev2,' // ',ee
                    ! print *, 'nint = ', edge_clip_sm(ee)%nfint(1)
                    ! print *, cell_from_face(1:edge_clip_sm(ee)%nfint(1))
                    ! print *, edge_clip_sm(ee)%vmface(1,1:edge_clip_sm(ee)%nfint(1))
                end if 
            end if 
        end if 

        !Assign end vertex
        if (surface_mesh%vtx_vmesh_cell(ev2) == 0) then 
            if (vtx_clip_sm(ev2)%nintersect == 0) then !Not an intersecting vertex 
                if (edge_clip_sm(ee)%nfint(edge_clip_sm(ee)%nint) == 1) then !If a proper intersection with a face 

                    !Target face 
                    ftgt = edge_clip_sm(ee)%vmface(edge_clip_sm(ee)%nint,1)

                    !Cells 
                    cleft = volume_mesh_full%faces(ftgt)%cleft
                    cright = volume_mesh_full%faces(ftgt)%cright

                    !Face normal vector 
                    Nf = newell_normal(volume_mesh_full%faces(ftgt)%nvtx,volume_mesh_full%faces(ftgt)%vertices,volume_mesh_full%vtx)

                    !Assign 
                    dotval = dot_product(edir,Nf)
                    if (dotval .LT. 0.0d0) then 
                        surface_mesh%vtx_vmesh_cell(ev2) = cleft
                    elseif (dotval .GT. 0.0d0) then 
                        surface_mesh%vtx_vmesh_cell(ev2) = cright
                    else
                        !print *, '** (fi) zero dot intersection sm edge -> vm face at edge ',ee
                    end if 
                elseif (edge_clip_sm(ee)%nfint(edge_clip_sm(ee)%nint) .GT. 1) then !If a proper intersection with a multiple faces

                    !Classify cell on all faces 
                    nfi = edge_clip_sm(ee)%nfint(edge_clip_sm(ee)%nint)
                    do ff=1,edge_clip_sm(ee)%nfint(edge_clip_sm(ee)%nint)

                        !Target face 
                        ftgt = edge_clip_sm(ee)%vmface(edge_clip_sm(ee)%nint,ff)

                        !Cells 
                        cleft = volume_mesh_full%faces(ftgt)%cleft
                        cright = volume_mesh_full%faces(ftgt)%cright

                        !Face normal vector 
                        Nf = newell_normal(volume_mesh_full%faces(ftgt)%nvtx,volume_mesh_full%faces(ftgt)%vertices,&
                        volume_mesh_full%vtx)

                        !Assign 
                        dotval = dot_product(edir,Nf)
                        if (dotval .LT. 0.0d0) then 
                            cell_from_face(ff) = cleft
                        elseif (dotval .GT. 0.0d0) then 
                            cell_from_face(ff) = cright
                        else
                            !print *, '** (ei) zero dot intersection sm edge -> vm face at edge ',ee
                            cell_from_face(ff) = 0
                        end if 
                    end do 

                    !Select cell with two agreeing faces 
                    cell_selected = 0 
                    do ff=1,nfi 
                        if (cell_from_face(ff) .GT. 0) then 
                            cell_select(cell_from_face(ff)) = cell_select(cell_from_face(ff)) + 1
                        end if 
                    end do 
                    do ff=1,nfi 
                        if (cell_from_face(ff) .GT. 0) then 
                            if (cell_select(cell_from_face(ff)) == 2) then 
                                cell_selected = cell_from_face(ff)
                                exit 
                            end if 
                        end if 
                    end do 

                    !Store 
                    if (cell_selected .GT. 0) then 
                        surface_mesh%vtx_vmesh_cell(ev2) = cell_selected
                    else    
                        print *, '** no cell selected at edge ',ee
                    end if 

                    !Reset selected cell counts 
                    do ff=1,nfi 
                        if (cell_from_face(ff) .GT. 0) then 
                            cell_select(cell_from_face(ff)) = 0
                        end if 
                    end do 

                    !Debug print 
                    ! print *, '============ mfint 2', ev2
                    ! print *, 'edge = ',ev1,ev2,' // ',ee
                    ! print *, 'nint = ', edge_clip_sm(ee)%nfint(edge_clip_sm(ee)%nint)
                    ! print *, cell_from_face(1:edge_clip_sm(ee)%nfint(edge_clip_sm(ee)%nint))
                    ! print *, edge_clip_sm(ee)%vmface(edge_clip_sm(ee)%nint,1:edge_clip_sm(ee)%nfint(edge_clip_sm(ee)%nint))
                end if 
            end if 
        end if 
    end if 
end do 

!Flood assignments though non clipped edges 
do ff=1,surface_mesh%nedge
    nupdate = 0 
    do ee=1,surface_mesh%nedge
        if (edge_clip_sm(ee)%nint == 0) then 

            !Vertices 
            ev1 = surface_mesh%edges(ee,1)
            ev2 = surface_mesh%edges(ee,2)

            !Flood
            if ((surface_mesh%vtx_vmesh_cell(ev1) == 0) .AND. (surface_mesh%vtx_vmesh_cell(ev2) .NE. 0)) then 
                if (vtx_clip_sm(ev1)%nintersect == 0) then 
                    surface_mesh%vtx_vmesh_cell(ev1) = surface_mesh%vtx_vmesh_cell(ev2) 
                    nupdate = nupdate + 1
                end if 
            elseif ((surface_mesh%vtx_vmesh_cell(ev1) .NE. 0) .AND. (surface_mesh%vtx_vmesh_cell(ev2) == 0)) then 
                if (vtx_clip_sm(ev2)%nintersect == 0) then 
                    surface_mesh%vtx_vmesh_cell(ev2) = surface_mesh%vtx_vmesh_cell(ev1) 
                    nupdate = nupdate + 1
                end if 
            end if 
        end if 
    end do 
    !print *, nupdate
    if (nupdate == 0) then 
        exit 
    end if 
end do 
return 
end subroutine assign_smvtx_vmcells




!Clip surface edges to volume faces subroutine ===========================
subroutine clip_surfedges2volfaces(edge_clip_sm,vtx_clip_sm,edge_clip_vm,volume_mesh_full,surface_mesh,surface_adtree,&
                                   cm3dopt,vtx_idx_smesh)
implicit none 

!Variables - Import
integer(in) :: vtx_idx_smesh
type(edgeint), dimension(:), allocatable :: edge_clip_sm,edge_clip_vm
type(vtx_intersect), dimension(:) :: vtx_clip_sm
type(tree_data) :: surface_adtree
type(surface_data) :: surface_mesh
type(vol_mesh_data) :: volume_mesh_full
type(cm3d_options) :: cm3dopt

!Variables - Local 
character(len=2) :: vtxi_loc_sme,vtxi_loc_vmf,edgeloc,vtx_loc
integer(in) :: ff,ff2,nn,kk,ee,ii,aa
integer(in) :: nselected,nedge_selected,vm_intvalid,fadj,fadj_vm,ftgt
integer(in) :: ttgt,etgt,etgt_vm,vtgt,f1,f2,ev1,ev2,vft2,vft3,int_type,exist
integer(in) :: edges2int(surface_mesh%nedge),node_select(surface_adtree%nnode),edge_select(surface_mesh%nedge)
real(dp) :: cpadSZ,zxmin,zxmax,zymin,zymax,zzmin,zzmax,zitol,zitol_bc
real(dp) :: Nf(3),ve1(3),ve2(3),vt1(3),vt2(3),vt3(3),vint(3)
real(dp) :: vt1_s(3),vt2_s(3),vt3_s(3),vmf_midp(3),vint_bcv(3),vint_bcs(3),esurfnormal(3)
type(smvmintlist) :: smvm_pint(surface_mesh%nedge)

!Set comparison tollerance 
zitol = cm3dopt%intcointol !Set as intersection tollerance
zitol_bc = cm3dopt%baryloctol !Set as barycentric tollerance 

!Initialise selected nodes 
nselected = 0 
node_select(:) = 0 

!Initialise selected face type 
do ee=1,surface_mesh%nedge
    smvm_pint(ee)%nface = 0 
    smvm_pint(ee)%nentry = cm3dopt%NintEmax
    allocate(smvm_pint(ee)%faces(cm3dopt%NintEmax))
    allocate(smvm_pint(ee)%fint_sm_loc(cm3dopt%NintEmax))
    allocate(smvm_pint(ee)%fint_vm_loc(cm3dopt%NintEmax))
    allocate(smvm_pint(ee)%fint_vm_edge(cm3dopt%NintEmax))
    allocate(smvm_pint(ee)%fint_intersect_invalid(cm3dopt%NintEmax))
    allocate(smvm_pint(ee)%fint_int_loc(cm3dopt%NintEmax,3))
    smvm_pint(ee)%faces(:) = 0 
    smvm_pint(ee)%fint_intersect_invalid(:) = 0 
    smvm_pint(ee)%fint_int_loc(:,:) = 0.0d0 
end do 

!Initialise 
vtgt = 0 
nedge_selected = 0 
edge_select(:) = 0
edges2int(:) = 0 

!Accumulate volume mesh faces to clip against to each surface mesh edge 
do ff=1,volume_mesh_full%nface

    !Build face normal vector 
    Nf = newell_normal(volume_mesh_full%faces(ff)%nvtx,volume_mesh_full%faces(ff)%vertices,volume_mesh_full%vtx)

    !Padding size 
    cpadSZ = norm2(Nf)*0.005d0

    !Intersection bounding box
    zxmin = minval(volume_mesh_full%vtx(volume_mesh_full%faces(ff)%vertices(:),1)) - cpadSZ !tgt bounding box -> xmin
    zxmax = maxval(volume_mesh_full%vtx(volume_mesh_full%faces(ff)%vertices(:),1)) + cpadSZ !tgt bounding box -> xmax
    zymin = minval(volume_mesh_full%vtx(volume_mesh_full%faces(ff)%vertices(:),2)) - cpadSZ !tgt bounding box -> ymin
    zymax = maxval(volume_mesh_full%vtx(volume_mesh_full%faces(ff)%vertices(:),2)) + cpadSZ !tgt bounding box -> ymax
    zzmin = minval(volume_mesh_full%vtx(volume_mesh_full%faces(ff)%vertices(:),3)) - cpadSZ !tgt bounding box -> zmin
    zzmax = maxval(volume_mesh_full%vtx(volume_mesh_full%faces(ff)%vertices(:),3)) + cpadSZ !tgt bounding box -> zmax
   
    !Identify any triangle bounding boxes that may overlap the face 
    call search_ADtree(nselected,node_select,surface_adtree,zxmin,zxmax,zymin,zymax,zzmin,zzmax)

    !Clip if any overlapping surface faces have been found 
    if (nselected .NE. 0) then 

        !Build list of surface mesh edges on the targeted triangles 
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

        !Volume mesh face midpoint 
        vmf_midp(1) = sum(volume_mesh_full%vtx(volume_mesh_full%faces(ff)%vertices(:),1))/real(volume_mesh_full%faces(ff)%nvtx,dp)
        vmf_midp(2) = sum(volume_mesh_full%vtx(volume_mesh_full%faces(ff)%vertices(:),2))/real(volume_mesh_full%faces(ff)%nvtx,dp)
        vmf_midp(3) = sum(volume_mesh_full%vtx(volume_mesh_full%faces(ff)%vertices(:),3))/real(volume_mesh_full%faces(ff)%nvtx,dp)

        !Set as base triangle vertex
        vt1(:) = vmf_midp(:)

        !Check each edge for intersection on this volume mesh face 
        do ee=1,nedge_selected

            !Target edge 
            etgt = edges2int(ee)

            !Edge vertices 
            ev1 = surface_mesh%edges(etgt,1)
            ev2 = surface_mesh%edges(etgt,2)

            !Edge end vertices
            ve1(:) = surface_mesh%vertices(ev1,:)
            ve2(:) = surface_mesh%vertices(ev2,:)

            !Adjacent surface faces to this edge 
            f1 = surface_mesh%E2F(etgt,1)
            f2 = surface_mesh%E2F(etgt,2)

            !Approximate surface normal at this edge 
            if ((f1 .GT. 0) .AND. (f2 .GT. 0)) then 
                fadj = f1 
                esurfnormal(:) = 0.5d0*(surface_mesh%face_normal(f1,:) +  surface_mesh%face_normal(f2,:))
            elseif (f1 .GT. 0) then
                fadj = f1
                esurfnormal(:) = surface_mesh%face_normal(f1,:)
            elseif (f2 .GT. 0) then
                fadj = f2
                esurfnormal(:) = surface_mesh%face_normal(f2,:)
            end if 

            !Check triangulation of volume mesh face
            do aa=1,volume_mesh_full%faces(ff)%nvtx

                !Face edge end vertices 
                vft2 = aa 
                vft3 = mod(aa,volume_mesh_full%faces(ff)%nvtx) + 1
                vft2 = volume_mesh_full%faces(ff)%vertices(vft2)
                vft3 = volume_mesh_full%faces(ff)%vertices(vft3)
                vt2(:) = volume_mesh_full%vtx(vft2,:)
                vt3(:) = volume_mesh_full%vtx(vft3,:)

                !Test intersection 
                int_type = line_tri_intersect_bool(ve1,ve2,vt1,vt2,vt3)
                if (int_type == 1) then !add this vm face to this sm edge 

                    !Find intersection location 
                    vint = line_tri_intersect(ve1,ve2,vt1,vt2,vt3)
                    vint_bcv = baryc_line_tri_intersect(ve1,ve2,vt1,vt2,vt3)

                    !Classify the location on the volume mesh face 
                    vtxi_loc_vmf = vtx_bary_tri_loc(vint_bcv(1),vint_bcv(2),vint_bcv(3),zitol_bc)

                    !Set internal to volume mesh face if not on the outer edge (e2) or its vertices (v2-v3)
                    if ((vtxi_loc_vmf .NE. 'e2') .AND. (vtxi_loc_vmf .NE. 'v2') .AND. &
                    (vtxi_loc_vmf .NE. 'v3') .AND. (vtxi_loc_vmf .NE. 'ot')) then
                        vtxi_loc_vmf = 'in'
                        ! print *, vtxi_loc_vmf
                    end if 
                    if (vtxi_loc_vmf == 'ot') then 
                        print *, vint_bcv,' edge out = ',etgt
                    end if 

                    !Build virtual face with this edge and the surface normal
                    vt1_s(:) = ve1(:)
                    vt2_s(:) = ve2(:)
                    vt3_s(:) = 0.5d0*(ve1(:) + ve2(:)) + esurfnormal(:)*norm2(ve2(:) - ve1(:))

                    !Classify the location on this surface mesh edge by constructing a virtual face
                    vint_bcs = project_to_barycentric(vint,vt1_s,vt2_s,vt3_s)
                    vtxi_loc_sme = vtx_bary_tri_loc(vint_bcs(1),vint_bcs(2),vint_bcs(3),zitol_bc)

                    !Add to intersection structure 
                    if (smvm_pint(etgt)%nface+1 .GT. smvm_pint(etgt)%nentry) then !append if required
                        call append_smvm_pint_entry(smvm_pint,etgt,2*smvm_pint(etgt)%nentry)
                    end if
                    smvm_pint(etgt)%nface = smvm_pint(etgt)%nface + 1
                    smvm_pint(etgt)%faces(smvm_pint(etgt)%nface) = ff
                    smvm_pint(etgt)%fint_vm_edge(smvm_pint(etgt)%nface) = aa 
                    smvm_pint(etgt)%fint_vm_loc(smvm_pint(etgt)%nface) = vtxi_loc_vmf
                    smvm_pint(etgt)%fint_sm_loc(smvm_pint(etgt)%nface) = vtxi_loc_sme
                    smvm_pint(etgt)%fint_int_loc(smvm_pint(etgt)%nface,:) = vint(:)

                    !Exit search for this edge 
                    exit
                end if
            end do 
        end do 
    end if 
end do 

!Set invalid intersections on faces adjacent to a volume mesh edge that an intersection will be snapped to 
do ee=1,surface_mesh%nedge
    if (smvm_pint(ee)%nface .GT. 0) then 

        !Initialise validity 
        smvm_pint(ee)%fint_intersect_invalid(:) = 0 

        !Collapse any sets of intersects where at least one is within tollerance of a vm edge into a single intersect on that vm edge 
        do ff=1,smvm_pint(ee)%nface 
            if ((smvm_pint(ee)%fint_sm_loc(ff)(1:1) == 'v') .AND. (smvm_pint(ee)%fint_intersect_invalid(ff) == 0)) then !search only vertex intersects on the sm edge first to ensure these are prioritised
                if ((smvm_pint(ee)%fint_vm_loc(ff) == 'e2') .OR. (smvm_pint(ee)%fint_vm_loc(ff) == 'v2') &
                .OR. (smvm_pint(ee)%fint_vm_loc(ff) == 'v3')) then !If within tollerance of a volume mesh edge 

                    !Volume mesh edge 
                    etgt_vm = volume_mesh_full%faces(smvm_pint(ee)%faces(ff))%edges(smvm_pint(ee)%fint_vm_edge(ff))

                    !Check if any other faces on E2f for this vm edge are in smvm_pint(ee)%faces(ff) for this surface mesh edge 
                    do kk=1,4   

                        !Adjacent face 
                        fadj_vm = volume_mesh_full%E2F(etgt_vm,kk)

                        !If this is a valid face and not the current face 
                        if ((fadj_vm .GT. 0) .AND. (fadj_vm .NE. smvm_pint(ee)%faces(ff))) then 
                            do ff2=1,smvm_pint(ee)%nface !If fadj_vm is within smvm_pint(ee)%faces(:) then tag the intersect on fadj_vm (position ff2) as invalid 
                                if ((ff2 .NE. ff) .AND. (smvm_pint(ee)%faces(ff2) == fadj_vm)) then 
                                    smvm_pint(ee)%fint_intersect_invalid(ff2) = 1
                                end if 
                            end do 
                        end if 
                    end do 
                end if 
            end if 
        end do 
        do ff=1,smvm_pint(ee)%nface 
            if ((smvm_pint(ee)%fint_sm_loc(ff)(1:1) == 'e') .AND. (smvm_pint(ee)%fint_intersect_invalid(ff) == 0)) then !now search only edge intersects on the sm edge 
                if ((smvm_pint(ee)%fint_vm_loc(ff) == 'e2') .OR. (smvm_pint(ee)%fint_vm_loc(ff) == 'v2') &
                .OR. (smvm_pint(ee)%fint_vm_loc(ff) == 'v3')) then !If within tollerance of a volume mesh edge 

                    !Volume mesh edge 
                    etgt_vm = volume_mesh_full%faces(smvm_pint(ee)%faces(ff))%edges(smvm_pint(ee)%fint_vm_edge(ff))

                    !Check if any other faces on E2f for this vm edge are in smvm_pint(ee)%faces(ff) for this surface mesh edge 
                    do kk=1,4   

                        !Adjacent face 
                        fadj_vm = volume_mesh_full%E2F(etgt_vm,kk)

                        !If this is a valid face and not the current face 
                        if ((fadj_vm .GT. 0) .AND. (fadj_vm .NE. smvm_pint(ee)%faces(ff))) then 
                            do ff2=1,smvm_pint(ee)%nface !If fadj_vm is within smvm_pint(ee)%faces(:) then tag the intersect on fadj_vm (position ff2) as invalid 
                                if ((ff2 .NE. ff) .AND. (smvm_pint(ee)%faces(ff2) == fadj_vm)) then 
                                    smvm_pint(ee)%fint_intersect_invalid(ff2) = 1
                                end if 
                            end do 
                        end if 
                    end do 
                end if 
            end if 
        end do 
    end if 
end do 

!Add intersections on surface mesh edge vertices 
do ee=1,surface_mesh%nedge
    if (smvm_pint(ee)%nface .GT. 0) then 

        !Edge vertices 
        ev1 = surface_mesh%edges(ee,1)
        ev2 = surface_mesh%edges(ee,2)

        !Adjacent surface faces to this edge 
        f1 = surface_mesh%E2F(ee,1)
        f2 = surface_mesh%E2F(ee,2)

        !Approximate surface normal at this edge 
        if ((f1 .GT. 0) .AND. (f2 .GT. 0)) then 
            fadj = f1 
        elseif (f1 .GT. 0) then
            fadj = f1
        elseif (f2 .GT. 0) then
            fadj = f2
        end if 

        !Find the position of edge ee in fadj 
        do aa=1,3
            if (surface_mesh%F2E(fadj,aa) == ee) then 
                if (aa == 1) then 
                    edgeloc = 'e1'
                elseif (aa == 2) then 
                    edgeloc = 'e2'
                elseif (aa == 3) then 
                    edgeloc = 'e3'
                end if 
                exit
            end if 
        end do 

        !Accumulate valid intersects onto the surface mesh edge and push to the volume mesh edges where required 
        do ff=1,smvm_pint(ee)%nface 
            if (smvm_pint(ee)%fint_intersect_invalid(ff) == 0) then !if valid intersection
                if (smvm_pint(ee)%fint_sm_loc(ff)(1:1) == 'v') then !intersect lies on a vertex at the end of the edge (add to vtx_clip_sm)

                    !Set this intersection to invalid so its not added again on the next edge loop 
                    smvm_pint(ee)%fint_intersect_invalid(ff) = 1

                    !Select vertex to add to 
                    if (smvm_pint(ee)%fint_sm_loc(ff) == 'v1') then !on vertex 1
                        vtgt = ev1
                    elseif (smvm_pint(ee)%fint_sm_loc(ff) == 'v2') then !on vertex 2
                        vtgt = ev2
                    elseif (smvm_pint(ee)%fint_sm_loc(ff) == 'v3') then !on vertex 3
                        vtgt = 0
                        print *, '** edge-vertex selection error - picked vertex 3 on virtual face'
                    end if 

                    !Find this vertices position on fadj 
                    do aa=1,3
                        if (surface_mesh%connectivity(fadj,aa) == vtgt) then 
                            if (aa == 1) then 
                                vtx_loc = 'v1'
                            elseif (aa == 2) then 
                                vtx_loc = 'v2'
                            elseif (aa == 3) then 
                                vtx_loc = 'v3'
                            end if 
                            exit 
                        end if 
                    end do 

                    !Check if this face is new on this vertex 
                    exist = 0
                    do ii=1,vtx_clip_sm(vtgt)%nintersect
                        if (vtx_clip_sm(vtgt)%face_int_idx(ii) == smvm_pint(ee)%faces(ff)) then 
                            exist = 1
                            exit 
                        end if 
                    end do 

                    !Add if new
                    if (exist == 0) then 
                        vtx_clip_sm(vtgt)%nintersect = vtx_clip_sm(vtgt)%nintersect + 1
                        vtx_clip_sm(vtgt)%face_int_idx(vtx_clip_sm(vtgt)%nintersect) = smvm_pint(ee)%faces(ff) 
                    end if 

                    !Push to the volume mesh edge if required 
                    if ((smvm_pint(ee)%fint_vm_loc(ff) == 'e2') .OR. (smvm_pint(ee)%fint_vm_loc(ff) == 'v2') &
                    .OR. (smvm_pint(ee)%fint_vm_loc(ff) == 'v3')) then !If within tollerance of a volume mesh edge 

                        !Target volume mesh edge 
                        etgt_vm = volume_mesh_full%faces(smvm_pint(ee)%faces(ff))%edges(smvm_pint(ee)%fint_vm_edge(ff))

                        !Validate 
                        vm_intvalid = validate_vme_intersect(edge_clip_vm,surface_mesh,vtx_loc,fadj,etgt_vm)

                        !Add if valid 
                        if (vm_intvalid == 1) then 
                            if (edge_clip_vm(etgt_vm)%nint+1 .GT. edge_clip_vm(etgt_vm)%nitem) then !append
                                if (edge_clip_vm(etgt_vm)%nitem == 0) then 
                                    call append_edge_clip_entry(edge_clip_vm,cm3dopt%NintEmax,etgt_vm) !allocate
                                else
                                    call append_edge_clip_entry(edge_clip_vm,edge_clip_vm(etgt_vm)%nitem*2,etgt_vm) !extend
                                end if 
                                if (edge_clip_vm(etgt_vm)%nitem .GT. cm3dopt%max_int_size) then 
                                    cm3dopt%max_int_size = edge_clip_vm(etgt_vm)%nitem
                                end if 
                            end if

                            !Set edge type 
                            edge_clip_vm(etgt_vm)%type = 4 

                            !Increment count 
                            edge_clip_vm(etgt_vm)%nint = edge_clip_vm(etgt_vm)%nint + 1

                            !Store 
                            edge_clip_vm(etgt_vm)%intloc(edge_clip_vm(etgt_vm)%nint,:) = smvm_pint(ee)%fint_int_loc(ff,:)
                            edge_clip_vm(etgt_vm)%surfseg(edge_clip_vm(etgt_vm)%nint) = fadj
                            edge_clip_vm(etgt_vm)%int_tri_loc(edge_clip_vm(etgt_vm)%nint) = vtx_loc
                            edge_clip_vm(etgt_vm)%vtx_idx(edge_clip_vm(etgt_vm)%nint) = vtgt
                            edge_clip_vm(etgt_vm)%intfrac(edge_clip_vm(etgt_vm)%nint) = 0.0d0 

                            !Push this intersection to all volume mesh faces on this volume mesh edge in the surface mesh structure
                            do kk=1,4   

                                !Adjacent face 
                                fadj_vm = volume_mesh_full%E2F(etgt_vm,kk)

                                !If non zero face 
                                if (fadj_vm .GT. 0) then 

                                    !Check if this exists on this surface mesh edge
                                    exist = 0 
                                    do nn=1,vtx_clip_sm(vtgt)%nintersect 
                                        if (vtx_clip_sm(vtgt)%face_int_idx(nn) == fadj_vm) then 
                                            exist = 1
                                            exit 
                                        end if 
                                    end do 

                                    !Add if new 
                                    if (exist == 0) then 
                                        vtx_clip_sm(vtgt)%nintersect = vtx_clip_sm(vtgt)%nintersect  + 1
                                        vtx_clip_sm(vtgt)%face_int_idx(vtx_clip_sm(vtgt)%nintersect) = fadj_vm
                                    end if 
                                end if 
                            end do 
                        end if 
                    end if 
                end if 
            end if 
        end do 
    end if 
end do 

!Add intersections internal to surface mesh edges
do ee=1,surface_mesh%nedge
    if (smvm_pint(ee)%nface .GT. 0) then 

        !Edge vertices 
        ev1 = surface_mesh%edges(ee,1)
        ev2 = surface_mesh%edges(ee,2)

        !Adjacent surface faces to this edge 
        f1 = surface_mesh%E2F(ee,1)
        f2 = surface_mesh%E2F(ee,2)

        !Approximate surface normal at this edge 
        if ((f1 .GT. 0) .AND. (f2 .GT. 0)) then 
            fadj = f1 
        elseif (f1 .GT. 0) then
            fadj = f1
        elseif (f2 .GT. 0) then
            fadj = f2
        end if 

        !Find the position of edge ee in fadj 
        do aa=1,3
            if (surface_mesh%F2E(fadj,aa) == ee) then 
                if (aa == 1) then 
                    edgeloc = 'e1'
                elseif (aa == 2) then 
                    edgeloc = 'e2'
                elseif (aa == 3) then 
                    edgeloc = 'e3'
                end if 
                exit
            end if 
        end do 

        !Set invalid any intersections within the edge with faces where there is already an intersection at a vertex of this edge with this face
        !and this intersect will not be snapped to a volume mesh face edge
        do ff=1,smvm_pint(ee)%nface 
            if ((smvm_pint(ee)%fint_sm_loc(ff)(1:1) == 'e') .AND. (smvm_pint(ee)%fint_intersect_invalid(ff) == 0) &
            .AND. (smvm_pint(ee)%fint_vm_loc(ff) == 'in')) then

                !Target volume mesh face
                ftgt = smvm_pint(ee)%faces(ff)

                !Check vertex intersection structures on each end of this edge 
                exist = 0 
                do nn=1,vtx_clip_sm(ev1)%nintersect 
                    if (vtx_clip_sm(ev1)%face_int_idx(nn) == ftgt) then 
                        exist = 1
                        exit 
                    end if 
                end do 
                do nn=1,vtx_clip_sm(ev2)%nintersect 
                    if (vtx_clip_sm(ev2)%face_int_idx(nn) == ftgt) then 
                        exist = 1
                        exit 
                    end if 
                end do 

                !Set invalid if there is an intersection with this face at either end of the edge 
                if (exist == 1) then 
                    smvm_pint(ee)%fint_intersect_invalid(ff) = 1
                end if 
            end if 
        end do 

        !Accumulate valid intersects onto the surface mesh edge and push to the volume mesh edges where required 
        do ff=1,smvm_pint(ee)%nface 
            if (smvm_pint(ee)%fint_intersect_invalid(ff) == 0) then !if valid intersection
                if (smvm_pint(ee)%fint_sm_loc(ff)(1:1) == 'e') then !intersect lies on the edge (add to edge_clip_sm)

                    !Set this intersection to invalid so its not added again on the next edge loop 
                    smvm_pint(ee)%fint_intersect_invalid(ff) = 1

                    !Add this intersection to the edge structure 
                    if (edge_clip_sm(ee)%nint+1 .GT. edge_clip_sm(ee)%nitem) then !append
                        call append_edge_clip_entry(edge_clip_sm,edge_clip_sm(ee)%nitem*2,ee)
                        if (edge_clip_sm(ee)%nitem .GT. cm3dopt%max_int_size) then 
                            cm3dopt%max_int_size = edge_clip_sm(ee)%nitem
                        end if 
                    end if
                    edge_clip_sm(ee)%type = 4
                    edge_clip_sm(ee)%nint = edge_clip_sm(ee)%nint + 1
                    edge_clip_sm(ee)%intloc(edge_clip_sm(ee)%nint,:) = smvm_pint(ee)%fint_int_loc(ff,:)
                    vtx_idx_smesh = vtx_idx_smesh + 1
                    edge_clip_sm(ee)%vtx_idx(edge_clip_sm(ee)%nint) = vtx_idx_smesh
                    edge_clip_sm(ee)%nfint(edge_clip_sm(ee)%nint) = &
                    edge_clip_sm(ee)%nfint(edge_clip_sm(ee)%nint) + 1 
                    edge_clip_sm(ee)%vmface(edge_clip_sm(ee)%nint,&
                    edge_clip_sm(ee)%nfint(edge_clip_sm(ee)%nint)) = smvm_pint(ee)%faces(ff) 
        
                    !Push to the volume mesh edge if required 
                    if ((smvm_pint(ee)%fint_vm_loc(ff) == 'e2') .OR. (smvm_pint(ee)%fint_vm_loc(ff) == 'v2') &
                    .OR. (smvm_pint(ee)%fint_vm_loc(ff) == 'v3')) then !If within tollerance of a volume mesh edge 

                        !Target volume mesh edge 
                        etgt_vm = volume_mesh_full%faces(smvm_pint(ee)%faces(ff))%edges(smvm_pint(ee)%fint_vm_edge(ff))

                        !Validate 
                        vm_intvalid = validate_vme_intersect(edge_clip_vm,surface_mesh,edgeloc,fadj,etgt_vm)

                        !Add if valid 
                        if (vm_intvalid == 1) then 
                            if (edge_clip_vm(etgt_vm)%nint+1 .GT. edge_clip_vm(etgt_vm)%nitem) then !append
                                if (edge_clip_vm(etgt_vm)%nitem == 0) then 
                                    call append_edge_clip_entry(edge_clip_vm,cm3dopt%NintEmax,etgt_vm) !allocate
                                else
                                    call append_edge_clip_entry(edge_clip_vm,edge_clip_vm(etgt_vm)%nitem*2,etgt_vm) !extend
                                end if 
                                if (edge_clip_vm(etgt_vm)%nitem .GT. cm3dopt%max_int_size) then 
                                    cm3dopt%max_int_size = edge_clip_vm(etgt_vm)%nitem
                                end if 
                            end if

                            !Set edge type 
                            edge_clip_vm(etgt_vm)%type = 4 

                            !Increment count 
                            edge_clip_vm(etgt_vm)%nint = edge_clip_vm(etgt_vm)%nint + 1

                            !Store 
                            edge_clip_vm(etgt_vm)%intloc(edge_clip_vm(etgt_vm)%nint,:) = smvm_pint(ee)%fint_int_loc(ff,:)
                            edge_clip_vm(etgt_vm)%surfseg(edge_clip_vm(etgt_vm)%nint) = fadj 
                            edge_clip_vm(etgt_vm)%int_tri_loc(edge_clip_vm(etgt_vm)%nint) = edgeloc 
                            edge_clip_vm(etgt_vm)%vtx_idx(edge_clip_vm(etgt_vm)%nint) = vtx_idx_smesh
                            edge_clip_vm(etgt_vm)%intfrac(edge_clip_vm(etgt_vm)%nint) = 0.0d0 

                            !Push this intersection to all volume mesh faces on this volume mesh edge in the surface mesh structure
                            do kk=1,4   

                                !Adjacent face 
                                fadj_vm = volume_mesh_full%E2F(etgt_vm,kk)

                                !If non zero face 
                                if (fadj_vm .GT. 0) then 

                                    !Check if this exists on this surface mesh edge
                                    exist = 0 
                                    do nn=1,edge_clip_sm(ee)%nfint(edge_clip_sm(ee)%nint)
                                        if (edge_clip_sm(ee)%vmface(edge_clip_sm(ee)%nint,nn) == fadj_vm) then 
                                            exist = 1
                                            exit 
                                        end if 
                                    end do 

                                    !Add if new 
                                    if (exist == 0) then 
                                        edge_clip_sm(ee)%nfint(edge_clip_sm(ee)%nint) = &
                                        edge_clip_sm(ee)%nfint(edge_clip_sm(ee)%nint) + 1
                                        edge_clip_sm(ee)%vmface(edge_clip_sm(ee)%nint,&
                                        edge_clip_sm(ee)%nfint(edge_clip_sm(ee)%nint)) = fadj_vm 
                                    end if 
                                end if 
                            end do 
                        end if 
                    end if 
                end if 
            end if 
        end do 
    end if 
end do 

!Debug check for valid unused edge intersections ====
! do ee=1,surface_mesh%nedge
!     if (smvm_pint(ee)%nface .GT. 0) then 
!         exist = 0 
!         do ff=1,smvm_pint(ee)%nface
!             if (smvm_pint(ee)%fint_intersect_invalid(ff) == 0) then 
!                 exist = 1
!             end if 
!         end do 
!         if (exist == 1) then 
!             print *, '=== eint missed -> edg = ',ee
!             do ff=1,smvm_pint(ee)%nface
!                 if (smvm_pint(ee)%fint_intersect_invalid(ff) == 0) then 
!                     print *, smvm_pint(ee)%faces(ff),smvm_pint(etgt)%fint_sm_loc(ff)
!                 end if 
!             end do 
!         end if 
!     end if 
! end do 

!Debug write intersections on volume mesh ====
! open(11,file='io/vmedge_seint.dat')
! do ee=1,volume_mesh_full%nedge
!     if (edge_clip_vm(ee)%nint .GE. 1) then 
!         do ii=1,edge_clip_vm(ee)%nint
!             write(11,*) edge_clip_vm(ee)%intloc(ii,:)
!         end do 
!     end if 
! end do 
! close(11)
return 
end subroutine clip_surfedges2volfaces




!Function to validate volume mesh edge intersection ===========================
function validate_vme_intersect(edge_clip_vm,surface_mesh,int_location_on_tri,tritgt,vmedge) result(intvalid)
implicit none 

!Variables - Import
integer(in) :: vmedge,tritgt,intvalid
character(len=2) :: int_location_on_tri
type(surface_data) :: surface_mesh
type(edgeint), dimension(:), allocatable :: edge_clip_vm

!Variables - Local 
integer(in) :: nn,aa
integer(in) :: etgt,vtgt1,vtgt2,fadj,eadj,vadj,f_connected
character(len=2) :: int_location_on_tri_adj

!If no intersection so far then set this intersection to valid and return 
if (edge_clip_vm(vmedge)%nint == 0) then 
    intvalid = 1
    return 
end if 

!Initialise 
intvalid = 1

!Test 
if (int_location_on_tri .NE. 'in') then !if on the triangle boundary 
    if (int_location_on_tri(1:1) == 'e') then !on triangle edge -> check face across edge and those on both vertices at the edge ends 

        !Target edge and vertices of the new intersect 
        if (int_location_on_tri(2:2) == '1') then 
            etgt = 1 
            vtgt1 = 1
            vtgt2 = 2
        elseif (int_location_on_tri(2:2) == '2') then 
            etgt = 2 
            vtgt1 = 2
            vtgt2 = 3
        elseif (int_location_on_tri(2:2) == '3') then 
            etgt = 3 
            vtgt1 = 3
            vtgt2 = 1
        end if 
        etgt = surface_mesh%F2E(tritgt,etgt)
        vtgt1 = surface_mesh%connectivity(tritgt,vtgt1)
        vtgt2 = surface_mesh%connectivity(tritgt,vtgt2)

        !Invalidate if there are intersections on adjacent faces co-incident with this one 
        do nn=1,edge_clip_vm(vmedge)%nint

            !Adjacent face 
            fadj = edge_clip_vm(vmedge)%surfseg(nn)

            !Test if this face contains one of etgt/vtgt1/vtgt2
            f_connected = 0 
            do aa=1,3
                if ((surface_mesh%connectivity(fadj,aa) == vtgt1) .OR. &
                (surface_mesh%connectivity(fadj,aa) == vtgt2)) then 
                    f_connected = 1
                    exit 
                end if 
                if (surface_mesh%F2E(fadj,aa) == etgt) then 
                    f_connected = 1
                    exit 
                end if 
            end do  

            !If this face is connected check the intersection location of the other intersection and remove if its tagged as the same
            !(either if its on the same edge or on a vertex at the end of the edge as the two must be co-incident)
            if (f_connected == 1) then 
                int_location_on_tri_adj = edge_clip_vm(vmedge)%int_tri_loc(nn)
                if (int_location_on_tri_adj(1:1) == 'e') then !on edge
                    if (int_location_on_tri_adj(2:2) == '1') then 
                        eadj = 1 
                    elseif (int_location_on_tri_adj(2:2) == '2') then 
                        eadj = 2 
                    elseif (int_location_on_tri_adj(2:2) == '3') then 
                        eadj = 3 
                    end if
                    eadj = surface_mesh%F2E(fadj,eadj)
                    if (eadj == etgt) then !on the same edge so remove 
                        intvalid = 0 
                        return 
                    end if 
                elseif (int_location_on_tri_adj(1:1) == 'v') then !on vertex
                    if (int_location_on_tri_adj(2:2) == '1') then 
                        vadj = 1
                    elseif (int_location_on_tri_adj(2:2) == '2') then 
                        vadj = 2
                    elseif (int_location_on_tri_adj(2:2) == '3') then 
                        vadj = 3
                    end if
                    vadj = surface_mesh%connectivity(fadj,vadj)
                    if ((vadj == vtgt1) .OR. (vadj == vtgt2)) then !on a vertex on this edge in the other triangle so remove 
                        intvalid = 0 
                        return 
                    end if 
                end if 
            end if 
        end do 
    elseif (int_location_on_tri(1:1) == 'v') then !on triangle vertex -> check triangles on this vertex 

        !Target vertex
        if (int_location_on_tri(2:2) == '1') then 
            vtgt1 = 1
        elseif (int_location_on_tri(2:2) == '2') then 
            vtgt1 = 2
        elseif (int_location_on_tri(2:2) == '3') then 
            vtgt1 = 3
        end if
        vtgt1 = surface_mesh%connectivity(tritgt,vtgt1)

        !Invalidate intersections on adjacent faces that are co-incident with this one 
        do nn=1,edge_clip_vm(vmedge)%nint

            !Adjacent face 
            fadj = edge_clip_vm(vmedge)%surfseg(nn)

            !Test if this face contains vtgt1
            f_connected = 0 
            do aa=1,3
                if (surface_mesh%connectivity(fadj,aa) == vtgt1) then 
                    f_connected = 1
                    exit 
                end if 
            end do  

            !If this face is connected check the intersection location of the other intersection and remove if its tagged as the same
            !(if its on the same vertex)
            if (f_connected == 1) then 
                int_location_on_tri_adj = edge_clip_vm(vmedge)%int_tri_loc(nn)
                if (int_location_on_tri_adj(2:2) == '1') then 
                    vadj = 1
                elseif (int_location_on_tri_adj(2:2) == '2') then 
                    vadj = 2
                elseif (int_location_on_tri_adj(2:2) == '3') then 
                    vadj = 3
                end if
                vadj = surface_mesh%connectivity(fadj,vadj)
                if (vadj == vtgt1) then !on the same vertex so remove 
                    intvalid = 0 
                    return  
                end if 
            end if 
        end do 
    end if 
end if 
return 
end function validate_vme_intersect




!Clip volume mesh edges to surface subroutine ===========================
subroutine clip_voledges2surffaces(edge_clip_vm,edge_clip_sm,vtx_clip_sm,tri_clip_sm,volume_mesh_full,surface_adtree,surface_mesh,&
                                   cm3dopt,vtx_idx_smesh)
implicit none 

!Variables - Import
integer(in) :: vtx_idx_smesh
type(vol_mesh_data) :: volume_mesh_full
type(edgeint), dimension(:), allocatable :: edge_clip_sm,edge_clip_vm
type(vtx_intersect), dimension(:) :: vtx_clip_sm
type(triint), dimension(:), allocatable :: tri_clip_sm
type(cm3d_options) :: cm3dopt
type(surface_data) :: surface_mesh
type(tree_data) :: surface_adtree

!Variables - Local 
integer(in) :: ii,jj,ee,nn,kk,aa,vv
integer(in) :: ftgt,etgt,ev1,ev2,nselected,vadj,eadj,fadj_vm,fadj_sm,vtgt
integer(in) :: int_type,Neintersect,sint_idx,vmint_on_this_tri_exist
integer(in) :: smt_nint_on_vmface(4),node_select(surface_adtree%nnode)
real(dp) :: cpadSZ,zxmin,zxmax,zymin,zymax,zzmin,zzmax,segnorm,zitol,zitol_bc
real(dp) :: ve1(3),ve2(3),segdir(3),vt1(3),vt2(3),vt3(3),vint(3),vint_bc_sm(3)
character(len=2) :: vtxi_loc_smf

!Set comparison tollerance
zitol = cm3dopt%intcointol !Set as intersection tollerance
zitol_bc = cm3dopt%baryloctol !Set as barycentric tollerance 

!Initialise selected nodes 
nselected = 0 
node_select(:) = 0 

!Find intersections and classify edges
sint_idx = 0 
Neintersect = 0 
do ee=1,volume_mesh_full%nedge

    !Edge vertices 
    ev1 = volume_mesh_full%edges(ee,1)
    ev2 = volume_mesh_full%edges(ee,2)

    !Edge end vertices
    ve1(:) = volume_mesh_full%vtx(ev1,:)
    ve2(:) = volume_mesh_full%vtx(ev2,:)

    !Edge length and direction 
    segdir(:) = ve2(:) - ve1(:)
    segnorm = norm2(segdir(:))

    !Padding size 
    cpadSZ = segnorm*0.005d0 !10.0d0*zitol

    !Intersection bounding box
    zxmin = min(ve1(1),ve2(1)) - cpadSZ !tgt bounding box -> xmin
    zxmax = max(ve1(1),ve2(1)) + cpadSZ !tgt bounding box -> xmax
    zymin = min(ve1(2),ve2(2)) - cpadSZ !tgt bounding box -> ymin
    zymax = max(ve1(2),ve2(2)) + cpadSZ !tgt bounding box -> ymax
    zzmin = min(ve1(3),ve2(3)) - cpadSZ !tgt bounding box -> zmin
    zzmax = max(ve1(3),ve2(3)) + cpadSZ !tgt bounding box -> zmax

    !Identify any triangle bounding boxes that may overlap the edge 
    call search_ADtree(nselected,node_select,surface_adtree,zxmin,zxmax,zymin,zymax,zzmin,zzmax)

    !Clip if there are any potential surface overlaps 
    if (nselected .NE. 0) then 

        !Search all targeted faces 
        do nn=1,nselected
            do kk=1,surface_adtree%tree(node_select(nn))%nentry

                !Target face 
                ftgt = surface_adtree%tree(node_select(nn))%entry(kk)

                !Verticies of this surface mesh face
                vt1(:) = surface_mesh%vertices(surface_mesh%connectivity(ftgt,1),:)
                vt2(:) = surface_mesh%vertices(surface_mesh%connectivity(ftgt,2),:)
                vt3(:) = surface_mesh%vertices(surface_mesh%connectivity(ftgt,3),:)

                !Test intersection
                int_type = line_tri_intersect_bool(ve1,ve2,vt1,vt2,vt3)
                if (int_type == 1) then 

                    

                    !Count intersections of this surface face with each volume face on this volume mesh edge
                    smt_nint_on_vmface(:) = 0 
                    do aa=1,4

                        !Vm face
                        fadj_vm = volume_mesh_full%E2F(ee,aa)

                        !Search surface mesh triangle vertices 
                        do vv=1,3
                            vtgt = surface_mesh%connectivity(ftgt,vv)
                            do ii=1,vtx_clip_sm(vtgt)%nintersect
                                if (vtx_clip_sm(vtgt)%face_int_idx(ii) == fadj_vm) then 
                                    smt_nint_on_vmface(aa) = smt_nint_on_vmface(aa) + 1
                                    exit 
                                end if 
                            end do
                        end do 

                        !Serch surface mesh triangle edges 
                        do vv=1,3
                            etgt = surface_mesh%F2E(ftgt,vv)
                            do ii=1,edge_clip_sm(etgt)%nint
                                do jj=1,edge_clip_sm(etgt)%nfint(ii)
                                    if (edge_clip_sm(etgt)%vmface(ii,jj) ==  fadj_vm) then 
                                        smt_nint_on_vmface(aa) = smt_nint_on_vmface(aa) + 1
                                        exit 
                                    end if 
                                end do 
                            end do 
                        end do 

                        !Search surface mesh triangle internal? -> none will exist yet for this edge so not required
                    end do 

                    !Identify if any volume mesh edge intersections exist on this triangles edges and vertices 
                    vmint_on_this_tri_exist = 0 
                    do ii=1,edge_clip_vm(ee)%nint
                        fadj_sm = edge_clip_vm(ee)%surfseg(ii)
                        if (edge_clip_vm(ee)%int_tri_loc(ii)(1:1) == 'e') then 
                            if (edge_clip_vm(ee)%int_tri_loc(ii)(2:2) == '1') then 
                                eadj = 1 
                            elseif (edge_clip_vm(ee)%int_tri_loc(ii)(2:2) == '2') then 
                                eadj = 2 
                            elseif (edge_clip_vm(ee)%int_tri_loc(ii)(2:2) == '3') then 
                                eadj = 3 
                            end if
                            eadj = surface_mesh%F2E(fadj_sm,eadj)
                            if ((eadj == surface_mesh%F2E(ftgt,1)) .OR. &
                            (eadj == surface_mesh%F2E(ftgt,2)) .OR. &
                            (eadj == surface_mesh%F2E(ftgt,3))) then 
                                vmint_on_this_tri_exist = 1
                                exit 
                            end if 
                        elseif (edge_clip_vm(ee)%int_tri_loc(ii)(1:1) == 'v') then 
                            if (edge_clip_vm(ee)%int_tri_loc(ii)(2:2) == '1') then 
                                vadj = 1
                            elseif (edge_clip_vm(ee)%int_tri_loc(ii)(2:2) == '2') then 
                                vadj = 2
                            elseif (edge_clip_vm(ee)%int_tri_loc(ii)(2:2) == '3') then 
                                vadj = 3
                            end if
                            vadj = surface_mesh%connectivity(fadj_sm,vadj)
                            if ((vadj == surface_mesh%connectivity(ftgt,1)) .OR. &
                            (vadj == surface_mesh%connectivity(ftgt,2)) .OR. &
                            (vadj == surface_mesh%connectivity(ftgt,3))) then 
                                vmint_on_this_tri_exist = 1
                                exit 
                            end if 
                        end if 
                    end do 

                    !If there is requirement for an internal intersection then construct one 
                    !(if there is at most one intersection with vm faces on this vm edge on this triangle and there is no intersection on this vm edge with any edge or vertex of this triangle)
                    if ((maxval(smt_nint_on_vmface) .LE. 1) .AND. (vmint_on_this_tri_exist == 0)) then 

                        !Find intersection location 
                        vint = line_tri_intersect(ve1,ve2,vt1,vt2,vt3)
                        vint_bc_sm = baryc_line_tri_intersect(ve1,ve2,vt1,vt2,vt3)

                        !Classify the location on the surface mesh face 
                        vtxi_loc_smf = vtx_bary_tri_loc(vint_bc_sm(1),vint_bc_sm(2),vint_bc_sm(3),zitol_bc)
                        
                        !Add to volume mesh edge ------
                        !Append or allocate edge clip entry if required 
                        if (edge_clip_vm(ee)%nint+1 .GT. edge_clip_vm(ee)%nitem) then !append
                            if (edge_clip_vm(ee)%nitem == 0) then 
                                call append_edge_clip_entry(edge_clip_vm,cm3dopt%NintEmax,ee) !allocate
                            else
                                call append_edge_clip_entry(edge_clip_vm,edge_clip_vm(ee)%nitem*2,ee) !extend
                            end if 
                            if (edge_clip_vm(ee)%nitem .GT. cm3dopt%max_int_size) then 
                                cm3dopt%max_int_size = edge_clip_vm(ee)%nitem
                            end if 
                        end if

                        !Set vm edge type 
                        edge_clip_vm(ee)%type = 4 

                        !Increment count 
                        edge_clip_vm(ee)%nint = edge_clip_vm(ee)%nint + 1

                        !Store 
                        edge_clip_vm(ee)%intloc(edge_clip_vm(ee)%nint,:) = vint(:)
                        edge_clip_vm(ee)%surfseg(edge_clip_vm(ee)%nint) = ftgt
                        edge_clip_vm(ee)%int_tri_loc(edge_clip_vm(ee)%nint) = 'in' !vtxi_loc_smf
                        edge_clip_vm(ee)%intfrac(edge_clip_vm(ee)%nint) = 0.0d0 
                        vtx_idx_smesh = vtx_idx_smesh + 1
                        edge_clip_vm(ee)%vtx_idx(edge_clip_vm(ee)%nint) = vtx_idx_smesh
                        !Add to volume mesh edge ------

                        !Add to surface mesh triangle ------
                        !Append tri clip entry if required 
                        if (tri_clip_sm(ftgt)%nvtx+1 .GT. tri_clip_sm(ftgt)%nitem) then !append
                            call append_tri_clip_entry(tri_clip_sm,tri_clip_sm(ftgt)%nitem*2,ftgt)
                            if (tri_clip_sm(ftgt)%nitem .GT. cm3dopt%max_int_size) then 
                                cm3dopt%max_int_size = tri_clip_sm(ftgt)%nitem
                            end if 
                        end if
                        
                        !Store 
                        tri_clip_sm(ftgt)%nvtx = tri_clip_sm(ftgt)%nvtx + 1
                        tri_clip_sm(ftgt)%vtx_idx(tri_clip_sm(ftgt)%nvtx) = edge_clip_vm(ee)%vtx_idx(edge_clip_vm(ee)%nint)
                        tri_clip_sm(ftgt)%edge_int_idx(tri_clip_sm(ftgt)%nvtx) = ee 
                        tri_clip_sm(ftgt)%intloc(tri_clip_sm(ftgt)%nvtx,:) = edge_clip_vm(ee)%intloc(ii,:)

                        !Add all faces on this volume mesh edge to the intersect
                        do aa=1,4
                            fadj_vm = volume_mesh_full%E2F(ee,aa)
                            if (fadj_vm .GT. 0) then 
                                tri_clip_sm(ftgt)%nfint(tri_clip_sm(ftgt)%nvtx) = &
                                tri_clip_sm(ftgt)%nfint(tri_clip_sm(ftgt)%nvtx) + 1
                                tri_clip_sm(ftgt)%face_int_idx(tri_clip_sm(ftgt)%nvtx,&
                                tri_clip_sm(ftgt)%nfint(tri_clip_sm(ftgt)%nvtx)) = fadj_vm
                            end if 
                        end do 
                        !Add to surface mesh triangle ------
                    end if 
                end if 
            end do 
        end do 
    end if 
end do 

!Debug write all vm-surface intersections ------- 
! open(11,file='io/vm_surf_intersect_in.dat')
! open(12,file='io/vm_surf_intersect_out.dat')
! do ee=1,volume_mesh_full%nedge
!     if (edge_clip(ee)%type == 4) then 
!         ! print *, 'et4'
!         do aa=1,edge_clip(ee)%nint
!             if (edge_clip(ee)%int_inout(aa)(1:1) == 'i') then 
!                 write(11,*) edge_clip(ee)%intloc(aa,:)
!             elseif (edge_clip(ee)%int_inout(aa)(1:1) == 'o') then 
!                 write(12,*) edge_clip(ee)%intloc(aa,:)
!             end if 
!         end do 
!     end if
! end do 
! close(11)
! close(12)
return 
end subroutine clip_voledges2surffaces




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
real(dp) :: ve1(3),ve2(3)
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




!Subroutine to build volume mesh intersecting edge mesh ===========================
subroutine build_vol_intersect_edge_mesh(edge_clip_vm,volume_mesh_full,cm3dopt)
implicit none 

!Variables - Import
type(edgeint), dimension(:) :: edge_clip_vm
type(vol_mesh_data) :: volume_mesh_full
type(cm3d_options) :: cm3dopt

!Variables - Local 
integer(in) :: ii,ee,kk
integer(in) :: ev1,ev2,intselect
integer(in) :: int_selected(cm3dopt%max_int_size)
integer(in) :: vtx_idx_temp(cm3dopt%max_int_size)
integer(in) :: surfseg_temp(cm3dopt%max_int_size)
real(dp) :: ve1(3),ve2(3)
real(dp) :: intdist(cm3dopt%max_int_size)
real(dp) :: intloc_temp(cm3dopt%max_int_size,3)
character(len=2) :: int_tri_loc_temp(cm3dopt%max_int_size)

!Order intersections and build edge mesh where neccesary
intdist(:) = 0.0d0 
int_selected(:) = 0
surfseg_temp(:) = 0 
vtx_idx_temp(:) = 0 
intloc_temp(:,:) = 0.0d0 
do ee=1,volume_mesh_full%nedge
    if (edge_clip_vm(ee)%nint .NE. 0) then !if any intersections 

        !Edge ends
        ev1 = volume_mesh_full%edges(ee,1)
        ev2 = volume_mesh_full%edges(ee,2)

        !Edge end vertices 
        ve1(:) = volume_mesh_full%vtx(ev1,:)
        ve2(:) = volume_mesh_full%vtx(ev2,:)

        !If more than one intersection then order intersections 
        if (edge_clip_vm(ee)%nint .GT. 1) then

            !Build intersection distances from end 1
            do ii=1,edge_clip_vm(ee)%nint
                intdist(ii) = norm2(edge_clip_vm(ee)%intloc(ii,:) - ve1(:))
            end do 

            !Order 
            int_selected(:) = 0
            intloc_temp(1:edge_clip_vm(ee)%nint,:) = edge_clip_vm(ee)%intloc(1:edge_clip_vm(ee)%nint,:) 
            surfseg_temp(1:edge_clip_vm(ee)%nint) = edge_clip_vm(ee)%surfseg(1:edge_clip_vm(ee)%nint) 
            int_tri_loc_temp(1:edge_clip_vm(ee)%nint) = edge_clip_vm(ee)%int_tri_loc(1:edge_clip_vm(ee)%nint) 
            vtx_idx_temp(1:edge_clip_vm(ee)%nint) = edge_clip_vm(ee)%vtx_idx(1:edge_clip_vm(ee)%nint) 
            do ii=1,edge_clip_vm(ee)%nint

                !Select intersection at maximum distance
                intselect = maxloc(intdist(1:edge_clip_vm(ee)%nint),1,&
                                   mask = int_selected(1:edge_clip_vm(ee)%nint) == 0)

                !Store
                edge_clip_vm(ee)%intloc(edge_clip_vm(ee)%nint-ii+1,:) = intloc_temp(intselect,:)
                edge_clip_vm(ee)%surfseg(edge_clip_vm(ee)%nint-ii+1) = surfseg_temp(intselect)
                edge_clip_vm(ee)%int_tri_loc(edge_clip_vm(ee)%nint-ii+1) = int_tri_loc_temp(intselect)
                edge_clip_vm(ee)%vtx_idx(edge_clip_vm(ee)%nint-ii+1) = vtx_idx_temp(intselect)
                edge_clip_vm(ee)%intfrac(edge_clip_vm(ee)%nint-ii+1) = intdist(intselect)

                !Set distance to zero and set to visited
                int_selected(intselect) = 1
                intdist(intselect) = 0.0d0 
            end do 
            ! print *, vtx_idx_temp(1:edge_clip_sm(ee)%nint)
            ! print *, edge_clip_vm(ee)%vtx_idx(1:edge_clip_sm(ee)%nint)
        end if 

        !Build edge sub-mesh ---- 
        !Allocate
        allocate(edge_clip_vm(ee)%edge_mesh(edge_clip_vm(ee)%nint+1,2))
            
        !Set first segment (nagative numbers index the intersections along the edge)
        edge_clip_vm(ee)%edge_mesh(1,1) = ev1
        edge_clip_vm(ee)%edge_mesh(1,2) = -1

        !Set middle segments 
        if (edge_clip_vm(ee)%nint .GE. 2) then 
            do kk=2,edge_clip_vm(ee)%nint
                edge_clip_vm(ee)%edge_mesh(kk,1) = -1*(kk - 1)
                edge_clip_vm(ee)%edge_mesh(kk,2) = -1*kk
            end do 
        end if 

        !Set final segment 
        edge_clip_vm(ee)%edge_mesh(edge_clip_vm(ee)%nint+1,1) = -1*edge_clip_vm(ee)%nint
        edge_clip_vm(ee)%edge_mesh(edge_clip_vm(ee)%nint+1,2) = ev2
    end if 
end do 
return 
end subroutine build_vol_intersect_edge_mesh




!Subroutine to perturb volume mesh vertices that lie exactly on geometry surface to off the surface to avoid degenerate intersections ===========================
subroutine perturb_onsurf_vmesh_vertices(volume_mesh_full,surface_adtree,surface_mesh,cm3dopt)
implicit none 

!Variables - Import
type(vol_mesh_data) :: volume_mesh_full
type(cm3d_options) :: cm3dopt
type(surface_data) :: surface_mesh
type(tree_data) :: surface_adtree

!Variables - Local 
integer(in) :: vv,aa,nn,kk
integer(in) :: ftgt,etgt,nselected,surf_coin,tri_tgt,Npert
integer(in) :: node_select(surface_adtree%nnode)
real(dp) :: cpadSZ,zxmin,zxmax,zymin,zymax,zzmin,zzmax,enorm,enorm_max,zitol
real(dp) :: vf1(3),vf2(3),vf3(3),vid(4)

!Set comparison tollerance 
zitol = 2.0d0*cm3dopt%intcointol !twice the intersection tollerance 

!Initialise selected nodes 
nselected = 0 
node_select(:) = 0 

!Identify and perturb vertices 
Npert = 0
do vv=1,volume_mesh_full%nvtx

    !Padding size 
    enorm_max = 0.0d0 
    do aa=1,volume_mesh_full%valence(vv)
        etgt = volume_mesh_full%v2e(vv,aa)
        enorm = norm2(volume_mesh_full%vtx(volume_mesh_full%edges(etgt,2),:) - &
        volume_mesh_full%vtx(volume_mesh_full%edges(etgt,1),:))
        if (enorm .GT. enorm_max) then 
            enorm_max = enorm 
        end if 
    end do 
    cpadSZ = enorm_max*0.005d0 !10.0d0!*zitol*enorm_max

    !Intersection bounding box
    zxmin = volume_mesh_full%vtx(vv,1) - cpadSZ !tgt bounding box -> xmin
    zxmax = volume_mesh_full%vtx(vv,1) + cpadSZ !tgt bounding box -> xmax
    zymin = volume_mesh_full%vtx(vv,2) - cpadSZ !tgt bounding box -> ymin
    zymax = volume_mesh_full%vtx(vv,2) + cpadSZ !tgt bounding box -> ymax
    zzmin = volume_mesh_full%vtx(vv,3) - cpadSZ !tgt bounding box -> zmin
    zzmax = volume_mesh_full%vtx(vv,3) + cpadSZ !tgt bounding box -> zmax

    !Identify any triangle bounding boxes that may overlap the vertex 
    call search_ADtree(nselected,node_select,surface_adtree,zxmin,zxmax,zymin,zymax,zzmin,zzmax)

    !If any triangles 
    if (nselected .NE. 0) then 

        !Check if triangle is surface co-incident 
        surf_coin = 0 
        do nn=1,nselected
            do kk=1,surface_adtree%tree(node_select(nn))%nentry

                !Target face 
                ftgt = surface_adtree%tree(node_select(nn))%entry(kk)

                !Verticies of this face
                vf1(:) = surface_mesh%vertices(surface_mesh%connectivity(ftgt,1),:)
                vf2(:) = surface_mesh%vertices(surface_mesh%connectivity(ftgt,2),:)
                vf3(:) = surface_mesh%vertices(surface_mesh%connectivity(ftgt,3),:)

                !Distance to this triangle 
                vid = min_dist_point_to_tri(volume_mesh_full%vtx(vv,:),vf1,vf2,vf3)

                !If within tollerance of the triangle then set this as surface co-incident 
                if (vid(4) .LE. zitol) then 
                    tri_tgt = ftgt
                    surf_coin = 1
                    exit 
                end if 
            end do
            if (surf_coin == 1) then 
                exit 
            end if 
        end do 

        !If surface co-incident then perturb away from the surface 
        if (surf_coin == 1) then 
            volume_mesh_full%vtx(vv,:) = volume_mesh_full%vtx(vv,:) + 0.01d0*enorm_max*surface_mesh%face_normal(tri_tgt,:)
            Npert = Npert + 1
            ! print *, vv
        end if 
    end if 
end do 

!Display
if (cm3dopt%dispt == 1) then
    write(*,'(A,I0,A)') '    {perturbed ',Npert,' volume vertices}'
end if
return 
end subroutine perturb_onsurf_vmesh_vertices




!Subroutine to perturb surface mesh vertices so they dont lie exactly on volume mesh faces ===========================
subroutine perturb_onmesh_smesh_vertices(volume_mesh_full,surface_adtree,surface_mesh,cm3dopt)
implicit none 

!Variables - Import
type(vol_mesh_data) :: volume_mesh_full
type(cm3d_options) :: cm3dopt
type(surface_data) :: surface_mesh
type(tree_data) :: surface_adtree

!Variables - Local 
integer(in) :: ff,vv,aa,nn,kk,ee
integer(in) :: etgt,nselected,surf_coin,Npert,nvtx_selected,ttgt,vtgt,vft1,vft2,vft3
integer(in) :: node_select(surface_adtree%nnode),vtx2int(surface_mesh%nvtx),vtx_select(surface_mesh%nvtx)
real(dp) :: cpadSZ,zxmin,zxmax,zymin,zymax,zzmin,zzmax,zitol
real(dp) :: vid(4),Nf(3),vt1(3),vt2(3),vt3(3),vp(3)

!Set comparison tollerance 
zitol = 2.0d0*cm3dopt%intcointol !twice the intersection tollerance 

!Initialise selected nodes 
nselected = 0 
node_select(:) = 0 

!Identify and perturb vertices 
Npert = 0
vtx2int(:) = 0
vtx_select(:) = 0 
do ff=1,volume_mesh_full%nface

    !Build face normal vector 
    Nf = newell_normal(volume_mesh_full%faces(ff)%nvtx,volume_mesh_full%faces(ff)%vertices,volume_mesh_full%vtx)

    !Padding size 
    cpadSZ = norm2(Nf)*0.005d0 !*10.0d0*zitol

    !Intersection bounding box
    zxmin = minval(volume_mesh_full%vtx(volume_mesh_full%faces(ff)%vertices(:),1)) - cpadSZ !tgt bounding box -> xmin
    zxmax = maxval(volume_mesh_full%vtx(volume_mesh_full%faces(ff)%vertices(:),1)) + cpadSZ !tgt bounding box -> xmax
    zymin = minval(volume_mesh_full%vtx(volume_mesh_full%faces(ff)%vertices(:),2)) - cpadSZ !tgt bounding box -> ymin
    zymax = maxval(volume_mesh_full%vtx(volume_mesh_full%faces(ff)%vertices(:),2)) + cpadSZ !tgt bounding box -> ymax
    zzmin = minval(volume_mesh_full%vtx(volume_mesh_full%faces(ff)%vertices(:),3)) - cpadSZ !tgt bounding box -> zmin
    zzmax = maxval(volume_mesh_full%vtx(volume_mesh_full%faces(ff)%vertices(:),3)) + cpadSZ !tgt bounding box -> zmax

    !Identify any triangle bounding boxes that may overlap the face 
    call search_ADtree(nselected,node_select,surface_adtree,zxmin,zxmax,zymin,zymax,zzmin,zzmax)

    !If any triangles 
    if (nselected .NE. 0) then 

        !Build list of vertices on these triangles 
        nvtx_selected = 0 
        do nn=1,nselected
            do kk=1,surface_adtree%tree(node_select(nn))%nentry
                ttgt = surface_adtree%tree(node_select(nn))%entry(kk)
                do ee=1,3
                    etgt = surface_mesh%F2E(ttgt,ee)
                    vtgt = surface_mesh%edges(etgt,1)
                    if (vtx_select(vtgt) == 0) then 
                        nvtx_selected = nvtx_selected + 1
                        vtx2int(nvtx_selected) = vtgt 
                        vtx_select(vtgt) = nvtx_selected
                    end if
                    vtgt = surface_mesh%edges(etgt,2)
                    if (vtx_select(vtgt) == 0) then 
                        nvtx_selected = nvtx_selected + 1
                        vtx2int(nvtx_selected) = vtgt 
                        vtx_select(vtgt) = nvtx_selected
                    end if
                end do 
            end do  
        end do 
        vtx_select(vtx2int(1:nvtx_selected)) = 0 
        
        !Set volume mesh face triangulation base vertex
        vft1 = volume_mesh_full%faces(ff)%vertices(1)
        vt1(:) = volume_mesh_full%vtx(vft1,:)

        !Perturb vertices co-incicdent with this volume mesh face 
        do vv=1,nvtx_selected

            !Target vertex
            vtgt = vtx2int(vv)
            vp(:) = surface_mesh%vertices(vtgt,:)

            !Search for proximity with this volume mesh face
            surf_coin = 0 
            do aa=2,volume_mesh_full%faces(ff)%nvtx-1

                !Triangle on this face 
                vft2 = aa 
                vft3 = aa + 1
                vft2 = volume_mesh_full%faces(ff)%vertices(vft2)
                vft3 = volume_mesh_full%faces(ff)%vertices(vft3)

                !Vertices 
                vt2(:) = volume_mesh_full%vtx(vft2,:)
                vt3(:) = volume_mesh_full%vtx(vft3,:)

                !Check for proximity 
                vid = min_dist_point_to_tri(vp,vt1,vt2,vt3)
                if (vid(4) .LE. zitol) then 
                    surf_coin = 1
                    exit 
                end if 
            end do 

            !If surface co-incident then perturb away from the surface 
            if (surf_coin == 1) then 
                surface_mesh%vertices(vtgt,:) = surface_mesh%vertices(vtgt,:) + 2.0d0*cpadSZ*surface_mesh%vtx_normal(vtgt,:)
                Npert = Npert + 1
            end if 
        end do 
    end if 
end do 

!Display
if (cm3dopt%dispt == 1) then
    write(*,'(A,I0,A)') '    {perturbed ',Npert,' surface vertices}'
end if
return 
end subroutine perturb_onmesh_smesh_vertices




!Classify edge-vertex gometry containment status subroutine ===========================
subroutine classify_edge_vertex_geom_containment(vtx_external,edge_clip_vm,volume_mesh_full,surface_mesh,cm3dopt,cm3dfailure)
implicit none 

!Variables - Import
integer(in) :: cm3dfailure
integer(in), dimension(:), allocatable :: vtx_external
type(edgeint), dimension(:), allocatable :: edge_clip_vm
type(vol_mesh_data) :: volume_mesh_full
type(surface_data) :: surface_mesh
type(cm3d_options) :: cm3dopt

!Variables - Local 
integer(in) :: ee,ff,vv,kk
integer(in) :: ev1,ev2,nupdate,nflooditer,Nexternal,Ninternal,ftgt,stri,vf1,vf2,vf3
real(dp) :: ve1(3),ve2(3),vt1(3),vt2(3),vt3(3)
character(len=2) int_type

!Classify each volume mesh surface intersection as entry/exit 
do ee=1,volume_mesh_full%nedge
    if (edge_clip_vm(ee)%type == 4) then 

        !Edge ends 
        ev1 = volume_mesh_full%edges(ee,1)
        ev2 = volume_mesh_full%edges(ee,2)

        !Edge vertices 
        ve1(:) = volume_mesh_full%vtx(ev1,:)
        ve2(:) = volume_mesh_full%vtx(ev2,:)

        !Classify each intersection
        allocate(edge_clip_vm(ee)%int_inout(edge_clip_vm(ee)%nint))
        do kk=1,edge_clip_vm(ee)%nint

            !Triangle vertices 
            ftgt = edge_clip_vm(ee)%surfseg(kk)
            vt1(:) = surface_mesh%vertices(surface_mesh%connectivity(ftgt,1),:)
            vt2(:) = surface_mesh%vertices(surface_mesh%connectivity(ftgt,2),:)
            vt3(:) = surface_mesh%vertices(surface_mesh%connectivity(ftgt,3),:)

            !Set intersect state 
            edge_clip_vm(ee)%int_inout(kk) = line_tri_intersect_type(ve1,ve2,vt1,vt2,vt3)
            ! print *, edge_clip(ee)%int_inout(kk)
        end do 
    end if 
end do 

!Initialise vtx_external
allocate(vtx_external(volume_mesh_full%nvtx))
vtx_external(:) = 0 

!Classify the vertices on ends of surface intersecting edges 
do ee=1,volume_mesh_full%nedge
    if (edge_clip_vm(ee)%type == 4) then 

        !Edge ends 
        ev1 = volume_mesh_full%edges(ee,1)
        ev2 = volume_mesh_full%edges(ee,2)

        !Classify end 1 
        if (edge_clip_vm(ee)%int_inout(1)(1:1) == 'i') then !intersection is in -> vertex 1 is external 
            vtx_external(ev1) = 1
        elseif (edge_clip_vm(ee)%int_inout(1)(1:1) == 'o') then !intersection is out -> vertex 1 is internal 
            vtx_external(ev1) = -1
        end if 

        !Classify end 2 
        if (edge_clip_vm(ee)%int_inout(edge_clip_vm(ee)%nint)(1:1) == 'i') then !intersection is in -> vertex 2 is internal 
            vtx_external(ev2) = -1
        elseif (edge_clip_vm(ee)%int_inout(edge_clip_vm(ee)%nint)(1:1) == 'o') then !intersection is out -> vertex 2 is external 
            vtx_external(ev2) = 1
        end if 
    end if 
end do 

!Flood vtx_external status through the mesh 
nflooditer = 0 
do ff=1,volume_mesh_full%nedge
    nupdate = 0 
    do ee=1,volume_mesh_full%nedge
        if (edge_clip_vm(ee)%type == 0) then  
        
            !Edge ends 
            ev1 = volume_mesh_full%edges(ee,1)
            ev2 = volume_mesh_full%edges(ee,2)

            !If mixed state
            if ((vtx_external(ev1) == 1) .AND. (vtx_external(ev2) == 0)) then 
                vtx_external(ev2) = 1
                nupdate = nupdate + 1
            elseif ((vtx_external(ev1) == 0) .AND. (vtx_external(ev2) == 1)) then 
                vtx_external(ev1) = 1
                nupdate = nupdate + 1
            elseif ((vtx_external(ev1) == -1) .AND. (vtx_external(ev2) == 0)) then 
                vtx_external(ev2) = -1
                nupdate = nupdate + 1
            elseif ((vtx_external(ev1) == 0) .AND. (vtx_external(ev2) == -1)) then 
                vtx_external(ev1) = -1
                nupdate = nupdate + 1
            end if 
        end if 
    end do 
    nflooditer = nflooditer + 1
    if (nupdate == 0) then 
        exit 
    end if 
end do 
Nexternal = 0 
Ninternal = 0 
do vv=1,volume_mesh_full%nvtx
    if (vtx_external(vv) == 1) then 
        Nexternal = Nexternal + 1
    elseif (vtx_external(vv) == -1) then 
        Ninternal = Ninternal + 1
    end if
end do 
if (cm3dopt%dispt == 1) then
    write(*,'(A,I0,A)') '    {vertex state flood completed after ',nflooditer,' iterations}'
    write(*,'(A,I0,A)') '    {identified ',Nexternal,' external vertices}'
    write(*,'(A,I0,A)') '    {identified ',Ninternal,' internal vertices}'
end if

!Classify edge and edge section states 
do ee=1,volume_mesh_full%nedge

    !Edge ends 
    ev1 = volume_mesh_full%edges(ee,1)
    ev2 = volume_mesh_full%edges(ee,2)

    !Classification state
    if (edge_clip_vm(ee)%type == 0) then !unclassified non-intersecting edge -> classify full edge
        
        !Classify full edge
        if ((vtx_external(ev1) == 1) .AND. (vtx_external(ev2) == 1)) then !external edge
            edge_clip_vm(ee)%type = 2
        elseif ((vtx_external(ev1) == -1) .AND. (vtx_external(ev2) == -1)) then !internal edge
            edge_clip_vm(ee)%type = 3
        else
            ! print *, '** type 0 intersecting edge identified: ',vtx_external(ev1),' -> ',vtx_external(ev2)
            ! print *, volume_mesh_full%vtx(ev1,:)
            ! print *, volume_mesh_full%vtx(ev2,:)
        end if 
    elseif (edge_clip_vm(ee)%type == 4) then !intersecting edge -> classify edge segments between intersections 

        !Allocate type structure 
        allocate(edge_clip_vm(ee)%int_seg_type(edge_clip_vm(ee)%nint+1))
        edge_clip_vm(ee)%int_seg_type(:) = 0 

        !Set first segment 
        if ((edge_clip_vm(ee)%int_inout(1)(1:1) == 'i') .AND. (vtx_external(ev1) == 1)) then !entry so segment external 
            edge_clip_vm(ee)%int_seg_type(1) = 2
        elseif ((edge_clip_vm(ee)%int_inout(1)(1:1) == 'o') .AND. (vtx_external(ev1) == -1)) then !exit so segment internal 
            edge_clip_vm(ee)%int_seg_type(1) = 3
        else 
            cm3dfailure = 1
            print *, '** edge intersection enter/exit disagreement with vtxinternal at edge: ',ee
            print *, edge_clip_vm(ee)%nint
            print *, edge_clip_vm(ee)%int_inout(1)
            print *, edge_clip_vm(ee)%type
            print *, vtx_external(ev1)
            print *, volume_mesh_full%vtx(ev1,:)
            print *, volume_mesh_full%vtx(ev2,:)
            ! return
        end if 

        !Set remaining segments 
        do kk=2,edge_clip_vm(ee)%nint+1
            if (edge_clip_vm(ee)%int_seg_type(kk-1) == 2) then !currently external
                if (edge_clip_vm(ee)%int_inout(kk-1)(1:1) == 'i') then !entry -> current segment internal
                    edge_clip_vm(ee)%int_seg_type(kk) = 3
                elseif (edge_clip_vm(ee)%int_inout(kk-1)(1:1) == 'o') then !exit from external -> warn (keep state)
                    edge_clip_vm(ee)%int_seg_type(kk) = 2
                    print *, '** warning: edge exiting geometry from external state '
                    print '(A,I0)', '========== edge = ',ee
                    print '(A,I0)', 'nint = ',edge_clip_vm(ee)%nint
                    do ff=1,edge_clip_vm(ee)%nint
                        print *,'intloc = ',edge_clip_vm(ee)%intloc(ff,:)
                        stri = edge_clip_vm(ee)%surfseg(ff)
                        print *, stri
                        vf1 = surface_mesh%connectivity(stri,1)
                        vf2 = surface_mesh%connectivity(stri,2)
                        vf3 = surface_mesh%connectivity(stri,3)
                        int_type = line_tri_intersect_type(volume_mesh_full%vtx(ev1,:),&
                        volume_mesh_full%vtx(ev2,:),&
                        surface_mesh%vertices(vf1,:),&
                        surface_mesh%vertices(vf2,:),&
                        surface_mesh%vertices(vf3,:))
                        print *,'int_type = ',int_type
                        print *,'loc on tri = ',edge_clip_vm(ee)%int_tri_loc(ff)
                    end do 
                    cm3dfailure = 1
                end if 
            elseif (edge_clip_vm(ee)%int_seg_type(kk-1) == 3) then !currently internal
                if (edge_clip_vm(ee)%int_inout(kk-1)(1:1) == 'i') then !entry -> warn (keep state)
                    edge_clip_vm(ee)%int_seg_type(kk) = 3
                    print *, '** warning: edge entering geometry from internal state '
                    print '(A,I0)', '========== edge = ',ee
                    print '(A,I0)', 'nint = ',edge_clip_vm(ee)%nint
                    do ff=1,edge_clip_vm(ee)%nint
                        print *,'intloc = ',edge_clip_vm(ee)%intloc(ff,:)
                        stri = edge_clip_vm(ee)%surfseg(ff)
                        print *, stri
                        vf1 = surface_mesh%connectivity(stri,1)
                        vf2 = surface_mesh%connectivity(stri,2)
                        vf3 = surface_mesh%connectivity(stri,3)
                        int_type = line_tri_intersect_type(volume_mesh_full%vtx(ev1,:),&
                        volume_mesh_full%vtx(ev2,:),&
                        surface_mesh%vertices(vf1,:),&
                        surface_mesh%vertices(vf2,:),&
                        surface_mesh%vertices(vf3,:))
                        print *,'int_type = ',int_type
                        print *,'loc on tri = ',edge_clip_vm(ee)%int_tri_loc(ff)
                    end do 
                    cm3dfailure = 1
                elseif (edge_clip_vm(ee)%int_inout(kk-1)(1:1) == 'o') then !exit from internal ->  current segment external
                    edge_clip_vm(ee)%int_seg_type(kk) = 2
                ! elseif (edge_clip(ee)%int_inout(kk-1) == 0) then !no transition so keep state
                !     edge_clip(ee)%int_seg_type(kk) = 3
                end if 
            end if 
        end do 
    end if
end do 

!Debug write reference external vertices ------- 
! open(11,file='io/vtxexternal.dat')
! do ee=1,volume_mesh_full%nvtx
!     if (vtx_external(ee) == 1) then 
!         write(11,*) volume_mesh_full%vtx(ee,:)
!     end if
! end do 
! close(11)
return 
end subroutine classify_edge_vertex_geom_containment




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
if (allocated(vmf_edges)) then 
    deallocate(vmf_edges)
end if 
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

! !Build v2e
! allocate(volume_mesh%v2e(volume_mesh%nvtx,6))
! allocate(volume_mesh%valence(volume_mesh%nvtx))
! volume_mesh%v2e(:,:) = 0 
! volume_mesh%valence(:) = 0 
! do ee=1,nedge
!     ev1 = vmf_edges(ee,1)
!     ev2 = vmf_edges(ee,2)
!     volume_mesh%valence(ev1) = volume_mesh%valence(ev1) + 1
!     volume_mesh%valence(ev2) = volume_mesh%valence(ev2) + 1
!     volume_mesh%v2e(ev1,volume_mesh%valence(ev1)) = ee
!     volume_mesh%v2e(ev2,volume_mesh%valence(ev2)) = ee 
! end do 

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
subroutine build_full_mesh(volume_mesh_full,ot_mesh,cm3dopt)
implicit none 
    
!Variables - Import
type(cm3d_options) :: cm3dopt
type(octree_data) :: ot_mesh
type(vol_mesh_data) :: volume_mesh_full

!Variables - Local 
integer(in) :: cc,ff,ee,vv
integer(in) :: nvtx_face,cl_face,cr_face,fidx,exist,edge_idx,pidx,maxDR,testDR,max_nevtx
integer(in) :: cadj,fop,nface,ncell,cadj_valid,vc1,vc2,cedge,Nsubcedge_vtx,etgt_adj,ctgt_adj,fadj,vselect,vins,etgt
integer(in) :: fopposite(6),ediagonal(12),edges(12,2),faces(6,4),face2edge(6,4),edge2face(12,2),edgeop_overface(12,2)
integer(in) :: cell_index(ot_mesh%cins-1),face_index(ot_mesh%cins-1,6),edge_nvtx(4)
integer(in), dimension(:), allocatable :: subcedge_vtx_visit
integer(in), dimension(:,:), allocatable :: edge_vtx_odr,subcedge_vtx
real(dp), dimension(:), allocatable :: subcedge_vtx_dist 

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

!Find maximum refinement level delta in the octree 
maxDR = 0 
do cc=1,ot_mesh%cins-1
    if (ot_mesh%cell_keep(cc) .NE. 0) then 
        testDR = max_adjacent_dR(ot_mesh,cc)
        if (testDR .GT. maxDR) then 
            maxDR = testDR
        end if 
    end if 
end do 
ot_mesh%maxDR = maxDR
max_nevtx = (2**maxDR - 1)*4

!Display
if (cm3dopt%dispt == 1) then
    write(*,'(A,I0,A)') '    {maximum adjacent cell refinement delta = ',maxDR,'}'
end if 

!Allocate adjacency arrays based on maxDR
allocate(subcedge_vtx_visit(max_nevtx))
allocate(subcedge_vtx(max_nevtx,3))
allocate(edge_vtx_odr(max_nevtx,4))
allocate(subcedge_vtx_dist(max_nevtx))

!Set maximum possible valence of each vertex in the octree
allocate(ot_mesh%vtx_valence(ot_mesh%vins-1))
ot_mesh%vtx_valence(:) = 0 
do cc=1,ot_mesh%cins-1
    if (ot_mesh%cell_keep(cc) .NE. 0) then 
        ot_mesh%vtx_valence(ot_mesh%cell_vcnr(cc,:)) = ot_mesh%vtx_valence(ot_mesh%cell_vcnr(cc,:)) + 3
    end if
end do 
ot_mesh%maxvalence = maxval(ot_mesh%vtx_valence)

!Build octree V2V and index edges 
ot_mesh%nedge = 0 
allocate(ot_mesh%edge_index(ot_mesh%vins-1,ot_mesh%maxvalence))
allocate(ot_mesh%cell2edge(ot_mesh%cins-1,12))
allocate(ot_mesh%V2V(ot_mesh%vins-1,ot_mesh%maxvalence))
allocate(ot_mesh%valence(ot_mesh%vins-1))
ot_mesh%edge_index(:,:) = 0 
ot_mesh%cell2edge(:,:) = 0 
ot_mesh%V2V(:,:) = 0 
ot_mesh%valence(:) = 0 
do cc=1,ot_mesh%cins-1
    if (ot_mesh%cell_keep(cc) .NE. 0) then 
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
                        ot_mesh%valence(vc1) = ot_mesh%valence(vc1) + 1
                        exit 
                    end if 
                end do 
                do vv=1,ot_mesh%maxvalence
                    if (ot_mesh%V2V(vc2,vv) == 0) then 
                        ot_mesh%V2V(vc2,vv) = vc1 
                        ot_mesh%edge_index(vc2,vv) = ot_mesh%nedge
                        ot_mesh%valence(vc2) = ot_mesh%valence(vc2) + 1
                        exit 
                    end if 
                end do 
            end if 
        end do 
    end if
end do 
!print *, 'max ot_mesh%valence => ',maxval(ot_mesh%valence(:),1)

!Build octree cell2edge
do cc=1,ot_mesh%cins-1
    if (ot_mesh%cell_keep(cc) .NE. 0) then 
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
            elseif (edge_idx .LT. 0) then !Edge already constructed so add to edges on this cell
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
    allocate(ot_mesh%edge_M_vtx(ee)%e_vertices(max_nevtx))
    ot_mesh%edge_M_vtx(ee)%e_vertices(:) = 0 
end do 
do cc=1,ot_mesh%cins-1
    if (ot_mesh%cell_keep(cc) .NE. 0) then 

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
    if (ot_mesh%cell_keep(cc) .NE. 0) then 
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
                        if (ot_mesh%cell_keep(cadj) .NE. 0) then 
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
    if (ot_mesh%cell_keep(cc) .NE. 0) then 
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




!Append edge intersect subroutine ===========================
subroutine append_edge_clip_entry(edge_clip,lenN,etgt)
implicit none 

!Variables - Import
integer(in) :: lenN,etgt
type(edgeint), dimension(:), allocatable :: edge_clip

!Variables - Local
integer(in) :: vtx_idxT(edge_clip(etgt)%nitem)
integer(in) :: nfintT(edge_clip(etgt)%nitem),surfsegT(edge_clip(etgt)%nitem)
integer(in) :: vmfaceT(edge_clip(etgt)%nitem,8)
real(dp) :: intlocT(edge_clip(etgt)%nitem,3)
character(len=2) :: int_tri_loc(edge_clip(etgt)%nitem)


!Store data for the edge etgt
if (edge_clip(etgt)%nitem .GT. 0) then 
    intlocT(:,:) = edge_clip(etgt)%intloc(:,:)
    vtx_idxT(:) = edge_clip(etgt)%vtx_idx(:)
    vmfaceT(:,:) = edge_clip(etgt)%vmface(:,:)
    nfintT(:) = edge_clip(etgt)%nfint(:)
    int_tri_loc(:) = edge_clip(etgt)%int_tri_loc(:)
    surfsegT(:) = edge_clip(etgt)%surfseg(:)
end if 

!Allocate new arrays 
if (allocated(edge_clip(etgt)%intloc)) then 
    deallocate(edge_clip(etgt)%intloc)
end if 
if (allocated(edge_clip(etgt)%vtx_idx)) then 
    deallocate(edge_clip(etgt)%vtx_idx)
end if 
if (allocated(edge_clip(etgt)%vmface)) then 
    deallocate(edge_clip(etgt)%vmface)
end if 
if (allocated(edge_clip(etgt)%nfint)) then 
    deallocate(edge_clip(etgt)%nfint)
end if 
if (allocated(edge_clip(etgt)%int_tri_loc)) then 
    deallocate(edge_clip(etgt)%int_tri_loc)
end if 
if (allocated(edge_clip(etgt)%intfrac)) then 
    deallocate(edge_clip(etgt)%intfrac)
end if
if (allocated(edge_clip(etgt)%surfseg)) then 
    deallocate(edge_clip(etgt)%surfseg)
end if
allocate(edge_clip(etgt)%intloc(lenN,3))
allocate(edge_clip(etgt)%vtx_idx(lenN))
allocate(edge_clip(etgt)%vmface(lenN,8))
allocate(edge_clip(etgt)%nfint(lenN))
allocate(edge_clip(etgt)%int_tri_loc(lenN))
allocate(edge_clip(etgt)%intfrac(lenN))
allocate(edge_clip(etgt)%surfseg(lenN))

!Store data in the new structure 
if (edge_clip(etgt)%nitem .GT. 0) then 
    edge_clip(etgt)%intloc(1:edge_clip(etgt)%nitem,:) = intlocT(:,:) 
    edge_clip(etgt)%vtx_idx(1:edge_clip(etgt)%nitem) = vtx_idxT(:)  
    edge_clip(etgt)%vmface(1:edge_clip(etgt)%nitem,:) = vmfaceT(:,:)  
    edge_clip(etgt)%nfint(1:edge_clip(etgt)%nitem) = nfintT(:) 
    edge_clip(etgt)%int_tri_loc(1:edge_clip(etgt)%nitem) = int_tri_loc(:)
    edge_clip(etgt)%surfseg(1:edge_clip(etgt)%nitem) = surfsegT(:)
end if 
edge_clip(etgt)%intloc(edge_clip(etgt)%nitem+1:lenN,:) = 0.0d0
edge_clip(etgt)%vtx_idx(edge_clip(etgt)%nitem+1:lenN) = 0
edge_clip(etgt)%vmface(edge_clip(etgt)%nitem+1:lenN,:) = 0
edge_clip(etgt)%nfint(edge_clip(etgt)%nitem+1:lenN) = 0
edge_clip(etgt)%surfseg(edge_clip(etgt)%nitem+1:lenN) = 0

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




!Append smvm_pint subroutine ===========================
subroutine append_smvm_pint_entry(smvm_pint,etgt,lenN)
implicit none

!Variables - Import
integer(in) :: etgt,lenN
type(smvmintlist), dimension(:) :: smvm_pint

!Variables - local 
integer(in) :: facesT(smvm_pint(etgt)%nentry)
integer(in) :: fint_vm_edgeT(smvm_pint(etgt)%nentry)
integer(in) :: fint_intersect_invalidT(smvm_pint(etgt)%nentry)
real(dp) :: fint_int_locT(smvm_pint(etgt)%nentry,3)
character(len=2) :: fint_vm_locT(smvm_pint(etgt)%nentry),fint_sm_locT(smvm_pint(etgt)%nentry)

!Store data
facesT(:) = smvm_pint(etgt)%faces(:)
fint_vm_edgeT(:) = smvm_pint(etgt)%fint_vm_edge(:)
fint_vm_locT(:) = smvm_pint(etgt)%fint_vm_loc(:)
fint_sm_locT(:) = smvm_pint(etgt)%fint_sm_loc(:)
fint_intersect_invalidT(:) = smvm_pint(etgt)%fint_intersect_invalid(:)
fint_int_locT(:,:) = smvm_pint(etgt)%fint_int_loc(:,:)

!Reallocate and resize
deallocate(smvm_pint(etgt)%faces)
deallocate(smvm_pint(etgt)%fint_vm_edge)
deallocate(smvm_pint(etgt)%fint_vm_loc)
deallocate(smvm_pint(etgt)%fint_sm_loc)
deallocate(smvm_pint(etgt)%fint_intersect_invalid)
deallocate(smvm_pint(etgt)%fint_int_loc)
allocate(smvm_pint(etgt)%faces(lenN))
allocate(smvm_pint(etgt)%fint_vm_edge(lenN))
allocate(smvm_pint(etgt)%fint_vm_loc(lenN))
allocate(smvm_pint(etgt)%fint_sm_loc(lenN))
allocate(smvm_pint(etgt)%fint_intersect_invalid(lenN))
allocate(smvm_pint(etgt)%fint_int_loc(lenN,3))

!Update data
smvm_pint(etgt)%faces(1:smvm_pint(etgt)%nentry) = facesT(:)
smvm_pint(etgt)%faces(smvm_pint(etgt)%nentry+1:lenN) = 0 
smvm_pint(etgt)%fint_vm_edge(1:smvm_pint(etgt)%nentry) = fint_vm_edgeT(:)
smvm_pint(etgt)%fint_vm_edge(smvm_pint(etgt)%nentry+1:lenN) = 0 
smvm_pint(etgt)%fint_vm_loc(1:smvm_pint(etgt)%nentry) = fint_vm_locT(:)
smvm_pint(etgt)%fint_sm_loc(1:smvm_pint(etgt)%nentry) = fint_sm_locT(:)
smvm_pint(etgt)%fint_intersect_invalid(1:smvm_pint(etgt)%nentry) = fint_intersect_invalidT(:)
smvm_pint(etgt)%fint_intersect_invalid(smvm_pint(etgt)%nentry+1:lenN) = 0 
smvm_pint(etgt)%fint_int_loc(1:smvm_pint(etgt)%nentry,:) = fint_int_locT(:,:)
smvm_pint(etgt)%fint_int_loc(smvm_pint(etgt)%nentry+1:lenN,:) = 0.0d0 
smvm_pint(etgt)%nentry = lenN
return 
end subroutine append_smvm_pint_entry




!Remove intersection from surface mesh edge subroutine ===========================
subroutine remove_surfedge_intersection(edge_clip,etgt,int2rem)
implicit none 

!Variables - Import
integer(in) :: etgt,int2rem
type(edgeint), dimension(:) :: edge_clip

!Variables - Local 
integer(in) :: ii 
integer(in) :: NintB,int_ins
integer(in) :: int_keep(edge_clip(etgt)%nitem)
type(edgeint) :: edge_clip_temp

!Safe return with no intersections
if (edge_clip(etgt)%nint == 0) then 
    return 
end if 

!Store clipped structure in temporary array 
edge_clip_temp%type = edge_clip(etgt)%type
edge_clip_temp%nint = edge_clip(etgt)%nint
allocate(edge_clip_temp%intloc(edge_clip(etgt)%nint,3))
edge_clip_temp%intloc(:,:) = edge_clip(etgt)%intloc(1:edge_clip(etgt)%nint,:)
allocate(edge_clip_temp%vtx_idx(edge_clip(etgt)%nint))
edge_clip_temp%vtx_idx(:) = edge_clip(etgt)%vtx_idx(1:edge_clip(etgt)%nint) 
allocate(edge_clip_temp%nfint(edge_clip(etgt)%nint))
edge_clip_temp%nfint(:) = edge_clip(etgt)%nfint(1:edge_clip(etgt)%nint) 
allocate(edge_clip_temp%vmface(edge_clip(etgt)%nint,8))
edge_clip_temp%vmface(:,:) = edge_clip(etgt)%vmface(1:edge_clip(etgt)%nint,:)

!Tag intersection to remove 
int_keep(:) = 1 
int_keep(int2rem) = 0

!Update intersection count 
NintB = edge_clip(etgt)%nint
edge_clip(etgt)%nint = edge_clip(etgt)%nint - 1

!Update cases 
if (edge_clip(etgt)%nint .LE. 0) then !no intersections remain 
    edge_clip(etgt)%nint = 0 
    edge_clip(etgt)%type = 0 
    edge_clip(etgt)%intloc(:,:) = 0.0d0 
    edge_clip(etgt)%vtx_idx(:) = 0 
    edge_clip(etgt)%nfint(:) = 0 
    edge_clip(etgt)%vmface(:,:) = 0 
else !intersections remain 
    
    !Reset arrays 
    edge_clip(etgt)%intloc(:,:) = 0.0d0 
    edge_clip(etgt)%vtx_idx(:) = 0 
    edge_clip(etgt)%nfint(:) = 0 
    edge_clip(etgt)%vmface(:,:) = 0 

    !Repopulate from temporary array 
    int_ins = 0 
    do ii=1,NintB
        if (int_keep(ii) == 1) then 
            int_ins = int_ins + 1
            edge_clip(etgt)%intloc(int_ins,:) = edge_clip_temp%intloc(ii,:)
            edge_clip(etgt)%vtx_idx(int_ins) = edge_clip_temp%vtx_idx(ii) 
            edge_clip(etgt)%nfint(int_ins) = edge_clip_temp%nfint(ii) 
            edge_clip(etgt)%vmface(int_ins,:) = edge_clip_temp%vmface(ii,:)
        end if 
    end do 
end if 
return 
end subroutine remove_surfedge_intersection




!Remove face reference from an intersection on a surface mesh edge subroutine ===========================
subroutine remove_surfedge_intersection_fref(edge_clip,etgt,ftgt,inttgt)
implicit none 

!Variables - Import
integer(in) :: etgt,inttgt,ftgt
type(edgeint), dimension(:) :: edge_clip

!Variables - Local 
integer(in) :: ii,fexist 
integer(in) :: NintB,int_ins
integer(in) :: fint_keep(8)
type(edgeint) :: edge_clip_temp

!Store data in temporary array 
allocate(edge_clip_temp%nfint(edge_clip(etgt)%nint))
edge_clip_temp%nfint(:) = edge_clip(etgt)%nfint(1:edge_clip(etgt)%nint) 
allocate(edge_clip_temp%vmface(edge_clip(etgt)%nint,8))
edge_clip_temp%vmface(:,:) = edge_clip(etgt)%vmface(1:edge_clip(etgt)%nint,:)

!Remove reference to face ftgt from intersection inttgt
fexist = 0 
fint_keep(:) = 1
do ii=1,edge_clip(etgt)%nfint(inttgt)
    if (edge_clip(etgt)%vmface(inttgt,ii) == ftgt) then 
        fint_keep(ii) = 0 
        fexist = 1
    end if 
end do 
if (fexist == 1) then 

    !If face exists on this edge then remove 
    NintB = edge_clip(etgt)%nfint(inttgt)
    edge_clip(etgt)%nfint(inttgt) = edge_clip(etgt)%nfint(inttgt) - 1
    edge_clip(etgt)%vmface(inttgt,:) = 0 
    int_ins = 0 
    do ii=1,NintB
        if (fint_keep(ii) == 1) then 
            int_ins = int_ins + 1
            edge_clip(etgt)%vmface(inttgt,int_ins) = edge_clip_temp%vmface(inttgt,ii)
        end if 
    end do 
    if (edge_clip(etgt)%nfint(inttgt) .LE. 0) then 
        edge_clip(etgt)%nfint(inttgt) = 0 
    end if 

    ! !If this results in no faces linked to intersection inttgt then remove the intersection inttgt
    ! if (edge_clip(etgt)%nfint(inttgt) .LE. 0) then 
    !     call remove_surfedge_intersection(edge_clip,etgt,inttgt)
    ! end if 
end if 
return 
end subroutine remove_surfedge_intersection_fref




!Flip mesh face function ===========================
function flip_face(Nvtx,face_initial) result(face_flipped)
implicit none 

!Variables - Import 
integer(in) :: Nvtx
integer(in) :: face_flipped(Nvtx)
integer(in) :: face_initial(Nvtx)

!Variables - Local 
integer(in) :: vv,vc

!Flip face
do vv=1,Nvtx
    vc = Nvtx - vv + 1
    face_flipped(vv) = face_initial(vc)
end do 
return 
end function flip_face


end module cellmesh3d_mesh_build_mod