!Quadtree Cutcell Surface Snapping Cell Containment Meshing Program
!Max Wood - mw16116@bristol.ac.uk
!Univeristy of Bristol - Department of Aerospace Engineering

!Version 1.0
!Updated 05-01-2023

!Main
program cell_mesh3d
use cellmesh3d_mesh_generation_mod  
implicit none

!Variables
integer(in) :: cm3dfailure
type(cm3d_options) :: cm3dopt
type(surface_data) :: surface_mesh
type(vol_mesh_data) :: volume_mesh

!Get command arguments 
call get_process_arguments(cm3dopt)

!Import options
call cm3d_import_options(cm3dopt)

!Import custom boundary condition zones
if (cm3dopt%set_customBCs == 1) then 
    call cm3d_import_customBC_zones(cm3dopt)
end if

!Load object surface data
if (cm3dopt%dispt == 1) then
    write(*,*) '--> importing geometry data'
end if
call import_surface_geometry(surface_mesh,cm3dopt)

!Construct mesh
call cell_mesh3d_mesh(volume_mesh,surface_mesh,cm3dopt,cm3dfailure)

!Export items 
call export_status(cm3dopt,cm3dfailure)
call export_volume_mesh(volume_mesh,cm3dopt)
call export_volume_mesh_flux(volume_mesh,cm3dopt)

!End
if (cm3dopt%dispt == 1) then
    stop
end if
end program cell_mesh3d