!3D Octree Cutcell Meshing Program
!Max Wood - mw16116@bristol.ac.uk
!Univeristy of Bristol - Department of Aerospace Engineering

!Version 2.1
!Updated 13-02-2023

!Main
program cell_mesh3d
use cellmesh3d_mesh_generation_mod  
implicit none

!Variables
logical :: is_selfintersecting
integer(in) :: cm3dfailure
type(cm3d_options) :: cm3dopt
type(surface_data) :: surface_mesh
type(vol_mesh_data) :: volume_mesh
real(dp), dimension(:,:), allocatable :: gradient_vol,gradient_surf

!Set option defaults 
call set_default_options(cm3dopt)

!Get command arguments 
call get_process_arguments(cm3dopt)

!Import options
call cm3d_import_options(cm3dopt)

!Initialisation display
if (cm3dopt%dispt == 1) then
    write(*,'(A)') ' '
    write(*,'(A)')'+--------------------------------------------+'
    write(*,'(A)')'|                Cell Mesh 3D                |'
    write(*,'(A)')'|         3D Cut-Cell Mesh Generator         |'
    write(*,'(A)')'|        Version 0.4.1 || 13/02/2024         |'
    write(*,'(A)')'|                 Max Wood                   |'
    write(*,'(A)')'|           University of Bristol            |'
    write(*,'(A)')'|    Department of Aerospace Engineering     |'
    write(*,'(A)')'+--------------------------------------------+'
    write(*,'(A)') ' '
end if

!Process requested mode 
if (cm3dopt%mode == 'check') then !geometry check mode 

    !Display
    if (cm3dopt%dispt == 1) then
        write(*,'(A)') '         == geometry checking mode =='
    end if

    !Load object surface data
    if (cm3dopt%dispt == 1) then
        write(*,'(A)') '--> importing geometry data'
    end if
    call import_surface_geometry(surface_mesh,cm3dopt)

    !Preprocess surface geometry 
    if (cm3dopt%dispt == 1) then
        write(*,'(A)') '--> preprocessing geometry data'
    end if
    call preprocess_surface_mesh(surface_mesh,cm3dopt)

    !Check for self intersections 
    if (cm3dopt%dispt == 1) then
        write(*,'(A)') '--> testing for self intersections'
    end if
    is_selfintersecting = is_self_intersecting(surface_mesh,cm3dopt)
    if (cm3dopt%dispt == 1) then
        write(*,'(A,L,A)') '    {is self intersecting: ',is_selfintersecting,'}'
    end if
    
    !Export check results
    call export_geometry_check(is_selfintersecting,cm3dopt)
elseif (cm3dopt%mode == 'mesh') then !mesh generation mode

    !Display
    if (cm3dopt%dispt == 1) then
        write(*,'(A)') '        == mesh contruction mode =='
    end if

    !Import custom boundary condition zones
    if (cm3dopt%set_customBCs == 1) then 
        call cm3d_import_customBC_zones(cm3dopt)
    end if
    
    !Load object surface data
    if (cm3dopt%dispt == 1) then
        write(*,'(A)') '--> importing geometry data'
    end if
    call import_surface_geometry(surface_mesh,cm3dopt)
    
    !Construct mesh
    call cell_mesh3d_mesh(volume_mesh,surface_mesh,cm3dopt,cm3dfailure)
    
    !Export items 
    call export_status(cm3dopt,cm3dfailure)
    call export_volume_mesh(volume_mesh,cm3dopt)
    call export_volume_mesh_flux(volume_mesh,cm3dopt)
    if (cm3dopt%glink_con == 1) then
        call export_vs2s_interpstruc(volume_mesh,surface_mesh,cm3dopt)
    end if 
elseif (cm3dopt%mode == 'project') then !volume to surface gradient projection mode 

    !Display
    if (cm3dopt%dispt == 1) then
        write(*,'(A)') '       == gradient projection mode =='
    end if

    !Load object surface data
    if (cm3dopt%dispt == 1) then
        write(*,'(A)') '--> importing geometry data'
    end if
    call import_surface_geometry(surface_mesh,cm3dopt)

    !Load volume mesh 
    if (cm3dopt%dispt == 1) then
        write(*,'(A)') '--> importing volume mesh'
    end if
    call import_volume_mesh_flux(volume_mesh,cm3dopt)

    !Load flow gradient file 
    if (cm3dopt%dispt == 1) then
        write(*,'(A)') '--> importing volume flow gradients'
    end if
    call import_flow_gradients(gradient_vol,volume_mesh,cm3dopt)

    !Preprocess surface mesh 
    if (cm3dopt%dispt == 1) then
        write(*,'(A)') '--> pre-processing the surface geometety'
    end if
    call preprocess_surface_mesh(surface_mesh,cm3dopt)

    !Project volume gradients to the surface 
    if (cm3dopt%dispt == 1) then
        write(*,'(A)') '--> projecting volume flow gradients to the surface geometry'
    end if
    if (cm3dopt%glink_type == 'rbf') then 
        call project_gradients_RBF(gradient_surf,gradient_vol,volume_mesh,surface_mesh,6_in,&
        cm3dopt%ADTminNodedivsize,0.0d0,cm3dopt)
    elseif (cm3dopt%glink_type == 'int') then 

    end if 

    !Export surface projected gradients 
    if (cm3dopt%dispt == 1) then
        write(*,'(A)') '--> exporting surface gradients'
    end if
    call export_surface_gradients(gradient_surf,surface_mesh,cm3dopt)
else
    write(*,'(A,A)') '** unknown mode option requested : ',cm3dopt%mode
    stop
end if 

!End
if (cm3dopt%dispt == 1) then
    stop
end if
end program cell_mesh3d