!Cell Mesh 3D Cutcell Mesh Generation Module 
!Max Wood - mw16116@bristol.ac.uk
!Univeristy of Bristol - Department of Aerospace Engineering

!Version 1.1
!Updated 16-10-2023

!Module
module cellmesh3d_mesh_generation_mod
use cellmesh3d_io_mod
use cellmesh3d_adtree_mod
use cellmesh3d_surface_mod
use cell_mesh3d_octree_mod
use cellmesh3d_mesh_build_mod
use cellmesh3d_postprocess_mod
contains

!cell_mesh3d meshing subroutine ================================================
subroutine cell_mesh3d_mesh(volume_mesh,surface_mesh,cm3dopt,cm3dfailure)
implicit none 

!Variables - Import ----
integer(in) :: cm3dfailure
type(cm3d_options) :: cm3dopt
type(surface_data) :: surface_mesh
type(vol_mesh_data) :: volume_mesh

!Variables - Local ----
!Counters 
integer(in) :: ii

!Object data
real(dp) :: obj_max_x,obj_max_y,obj_max_z,obj_min_x,obj_min_y,obj_min_z,obj_cx,obj_cy,obj_cz
real(dp), dimension(:,:), allocatable :: tvtx

!Mesh data
integer(in) :: MaxValence,Nmerge,Nmerge_fail,nmiter,nifail
real(dp), dimension(:), allocatable :: Cvol

!AD_Tree data
integer(in) :: Ndim,node_minDIVsize
real(dp) :: global_target_pad
type(tree_data) :: surface_adtree

!Octree 
type(octree_data) :: ot_mesh

!Initialisation -------------------------------------------------------
if (cm3dopt%dispt == 1) then
    write(*,'(A)') ' '
    write(*,'(A)')'+--------------------------------------------+'
    write(*,'(A)')'|                Cell Mesh 3D                |'
    write(*,'(A)')'|         3D Cut-Cell Mesh Generator         |'
    write(*,'(A)')'|        Version 0.0.2 || 16/10/2023         |'
    write(*,'(A)')'|                 Max Wood                   |'
    write(*,'(A)')'|           University of Bristol            |'
    write(*,'(A)')'|    Department of Aerospace Engineering     |'
    write(*,'(A)')'+--------------------------------------------+'
    write(*,'(A)') ' '
end if

!Set hardcoded parameters -------------------------------------------------------
!Set MaxValence parameter 
MaxValence = 6

!Set global object bounding box padding for adtree node containement
global_target_pad = 0.0d0

!Adtree number of dimensions (6)
Ndim = 6

!Set minimum divisible node size within the AD tree 
node_minDIVsize = 10

!Initialise failure tag
cm3dfailure = 0

!Set initial count state
volume_mesh%nvtx = 0 
volume_mesh%nedge = 0 
volume_mesh%nface = 0 
volume_mesh%ncell = 0 

!Property calculation -------------------------------------------------------
if (cm3dopt%dispt == 1) then
    write(*,'(A)') '--> constructing surface geometry parameters'
end if

!Orient surface mesh for positive object volume (this ensures normal vector convention is correct)
call orient_surface(surface_mesh,cm3dopt)

!Build surface mesh connectivity 
call get_valence(surface_mesh%valence,surface_mesh%maxValence,surface_mesh%connectivity,surface_mesh%nfcs,surface_mesh%nvtx)
call construct_edges(surface_mesh%nedge,surface_mesh%edges,surface_mesh%valence,surface_mesh%nvtx,surface_mesh%nfcs,&
surface_mesh%connectivity)
call get_connectivity(surface_mesh%V2E,surface_mesh%V2F,surface_mesh%F2E,surface_mesh%E2F,surface_mesh%valence,&
surface_mesh%maxvalence,surface_mesh%nedge,surface_mesh%nfcs,surface_mesh%nvtx,surface_mesh%edges,surface_mesh%connectivity)
call construct_surface_normals(surface_mesh)

!Evaluate surface mesh face curvature 
call evaluate_surf_rcurv(surface_mesh,cm3dopt)

!Global object bounds
obj_max_x = maxval(surface_mesh%vertices(:,1))
obj_max_y = maxval(surface_mesh%vertices(:,2))
obj_max_z = maxval(surface_mesh%vertices(:,3))
obj_min_x = minval(surface_mesh%vertices(:,1))
obj_min_y = minval(surface_mesh%vertices(:,2))
obj_min_z = minval(surface_mesh%vertices(:,3))

!Global object centre
obj_cx = 0.5d0*(obj_max_x + obj_min_x)
obj_cy = 0.5d0*(obj_max_y + obj_min_y)
obj_cz = 0.5d0*(obj_max_z + obj_min_z)

!Displays
if (cm3dopt%dispt == 1) then
    write(*,'(A)') ' '
    write(*,'(A)') '== Properties =========================='
    write(*,"(A,I0)") '   Imported surface vertices = ', surface_mesh%nvtx
    write(*,"(A,I0)") '   Imported surface faces = ', surface_mesh%nfcs
    write(*,"(A,F12.6,A,F12.6,A,F12.6)") '   Geometry dimensions (x/y/z) = ', obj_max_x - obj_min_x,' ', &
    obj_max_y - obj_min_y,' ', obj_max_z - obj_min_z
    write(*,'(A)') '== Properties =========================='
    write(*,'(A)') ' '
end if

!Construct AD_Tree on the target mesh ------------------------------------------
if (cm3dopt%dispt == 1) then
    write(*,'(A)') '--> constructing AD-tree on target geometry surface mesh'
end if

!Construct 6D bounding box coordinates for each face in the surface mesh
allocate(tvtx(surface_mesh%nfcs,6))
do ii=1,surface_mesh%nfcs
    tvtx(ii,1) = minval(surface_mesh%vertices(surface_mesh%connectivity(ii,:),1)) !xmin
    tvtx(ii,2) = minval(surface_mesh%vertices(surface_mesh%connectivity(ii,:),2)) !ymin
    tvtx(ii,3) = minval(surface_mesh%vertices(surface_mesh%connectivity(ii,:),3)) !zmin
    tvtx(ii,4) = maxval(surface_mesh%vertices(surface_mesh%connectivity(ii,:),1)) !xmax
    tvtx(ii,5) = maxval(surface_mesh%vertices(surface_mesh%connectivity(ii,:),2)) !ymax
    tvtx(ii,6) = maxval(surface_mesh%vertices(surface_mesh%connectivity(ii,:),3)) !zmax
end do

!Construct ad_tree
call build_ADtree(surface_adtree,ndim,cm3dopt%max_depth,node_minDIVsize,tvtx,global_target_pad,cm3dopt%dispt)

!Octree construction -------------------------------------------------------
if (cm3dopt%dispt == 1) then
    write(*,'(A)') '--> constructing mesh octree'
end if

!Octree refinement to geometry
call octree_mesh_refine(ot_mesh,surface_mesh,surface_adtree,cm3dopt,obj_cx,obj_cy,obj_cz,cm3dfailure)
if (cm3dfailure == 1) then
    return 
end if 

!Construct volume mesh -------------------------------------------------------
!Build mesh
! if (cm3dopt%dispt == 1) then
!     write(*,'(A)') '--> constructing full volume mesh: '
! end if
call construct_mesh(volume_mesh,ot_mesh,surface_mesh,surface_adtree,cm3dopt,cm3dfailure)
if (cm3dfailure == 1) then 
    return 
end if

!Ensure contiguous cell indecies 
call remap_cell_indecies(volume_mesh)

!Remap vertices 




!Postprocess complete mesh ---------------------------------------------------
! ! !Apply custom boundary conditions conditions in target regions 
! ! if (cm2dopt%set_customBCs == 1) then 
! !     call set_custom_bcs(volume_mesh,cm2dopt) 
! ! end if 

! ! !Remove mesh zones with any far field boundary conditions (for internal duct flows -> remove all regions outside the duct)
! ! if (cm2dopt%remFFzones == 1) then 
! !     call remove_farfield_adjacent(volume_mesh) !Remove duct external regions of the mesh
! !     if (volume_mesh%nedge == 0) then !Check for complete mesh removal 
! !         cm2dfailure = 1
! !         print *, '** mesh partitioning error -> no retained edges [remove far field zones]'
! !         return 
! !     end if
! ! end if

! ! !Remove mesh zones that are only connected to wall boundary conditions (isolated regions of the mesh with no inflow/outflow possible)
! ! if (cm2dopt%remISzones == 1) then 
! !     call remove_isolated_regions(volume_mesh)
! !     if (volume_mesh%nedge == 0) then !Check for complete mesh removal 
! !         cm2dfailure = 1
! !         print *, '** mesh partitioning error -> no retained edges [remove isolated zones]'
! !         return 
! !     end if
! ! end if 

!Display
if (cm3dopt%dispt == 1) then
    write(*,'(A)') '--> cleaning volume mesh'
end if

!Check for and correct bisected cells 
call correct_bisected_cells(volume_mesh,cm3dopt)

!Remove mesh internal valence two vertices 
call clean_vlnc2_vertices(volume_mesh,cm3dopt,0_in) 

!Simplify the mesh surface within each cell if simplified surface requested
if (cm3dopt%surface_type == 0) then     
    if (cm3dopt%dispt == 1) then
        write(*,'(A)') '--> constructing simplified surface geometry'
    end if
    call simplify_surface(volume_mesh,cm3dopt,cm3dfailure)
    call clean_vlnc2_vertices(volume_mesh,cm3dopt,1_in) !surface vertices
    if (cm3dfailure == 1) then 
        return 
    end if 
end if 

!Clean mesh by removing sliver cells 
nmiter = 10 !volume_mesh%ncell !set maximum merging iterations 
nifail = 0 
do ii=1,nmiter
    call clean_mesh_sliverC(volume_mesh,Cvol,cm3dopt,Nmerge,Nmerge_fail)
    if (Nmerge_fail .NE. 0) then 
        nifail = nifail + 1
    end if 
    if (Nmerge .NE. 0) then 
        nifail = nifail - 1
    end if 
    if ((Nmerge == 0) .AND. (Nmerge_fail .NE. 0) .AND. (nifail == 1)) then 
        exit 
    end if 
    if ((Nmerge == 0) .AND. (Nmerge_fail == 0)) then 
        exit 
    end if
end do 

!Remap cell indecies 
call remap_cell_indecies(volume_mesh)

! ! !Build surface vertex mappings 
! ! call build_surface_links(volume_mesh,surface_mesh)

!Flip surface and boundary normals if requested (both are constructed in state 'in')
if (cm3dopt%surface_dir == 'out') then 
    call flip_set_faces(volume_mesh,-1_in)
end if 
if (cm3dopt%boundary_dir == 'out') then 
    call flip_set_faces(volume_mesh,-2_in)
    call flip_set_faces(volume_mesh,-3_in)
    call flip_set_faces(volume_mesh,-4_in)
    call flip_set_faces(volume_mesh,-5_in)
    call flip_set_faces(volume_mesh,-6_in)
end if 

!Restructure mesh for SU2 output formats 
if (cm3dopt%meshfrmat == 'su2_dual') then 
    call construct_dual_mesh(volume_mesh)
end if 

!Remove any zero area faces 
call remove_zero_area_faces(volume_mesh,cm3dopt) !========================================

!Evaluate cell volumes to check for invalid cells 
call get_cell_volumes(Cvol,volume_mesh)

!Set mesh and surface face counts
volume_mesh%nface_surface = 0 
do ii=1,volume_mesh%nface
    if ((volume_mesh%faces(ii)%cleft .LT. 0) .OR. (volume_mesh%faces(ii)%cright .LT. 0)) then
        volume_mesh%nface_surface = volume_mesh%nface_surface + 1
    end if 
end do 
volume_mesh%nface_mesh = volume_mesh%nface - volume_mesh%nface_surface

!Completion of mesh construction display ----------------------
if (cm3dopt%dispt == 1) then
    write(*,'(A)') '--> mesh construction completed'
    write(*,'(A,I0,A)') '    {cells: ',volume_mesh%ncell,'}'
    write(*,'(A,I0,A)') '    {faces: ',volume_mesh%nface,'}'
    write(*,'(A,I0,A)') '    {mesh faces: ',volume_mesh%nface_mesh,'}'
    write(*,'(A,I0,A)') '    {surface faces: ',volume_mesh%nface_surface,'}'
    write(*,'(A,I0,A)') '    {vertices: ',volume_mesh%nvtx,'}'
    ! write(*,'(A,I0,A)') '    {surface vertices: ',volume_mesh%nvtx_surf,'}'
    write(*,'(A,E11.5,A,E11.5,A,E11.5,A)') '    {cell volume (max/min/total) : ',maxval(Cvol),' / ',minval(Cvol),' / ',sum(Cvol),'}'
end if
! print *, minloc(Cvol)
if (minval(Cvol) .LE. 0.0d0) then 
    cm3dfailure = 1
    do ii=1,volume_mesh%ncell
        if (Cvol(ii) .LT. 0.0d0) then 
            print '(A,I0)', '** negative volume cell identified: ',ii
        elseif (Cvol(ii) == 0.0d0) then  
            print '(A,I0)', '** zero volume cell identified: ',ii
        end if 
    end do 
end if  
return 
end subroutine cell_mesh3d_mesh


end module cellmesh3d_mesh_generation_mod