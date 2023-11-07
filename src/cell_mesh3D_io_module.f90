!cell_mesh2d io module
!Max Wood - mw16116@bristol.ac.uk
!Univeristy of Bristol - Department of Aerospace Engineering

!Version 6.2
!Updated 07-11-2023

!Module
module cellmesh3d_io_mod
use cellmesh3d_data_mod
contains

!Read command arguments subroutine ===========================
subroutine get_process_arguments(cm3dopt)
implicit none

!Variables - Import
type(cm3d_options) :: cm3dopt
    
!Variables - Local 
integer(in) :: nargs
integer(in32) :: arglen,argstat

!Check and process supplied command arguments 
nargs = command_argument_count()
if (nargs == 0) then 
    write(*,'(A)') '** at least one argument must be supplied [mode / surface name / io path / options path]'
    stop
elseif (nargs == 1) then !Use mode

    !Read mode 
    call get_command_argument(number=1, length=arglen)
    allocate(character(len=arglen) :: cm3dopt%mode)
    call get_command_argument(number=1, value=cm3dopt%mode, status=argstat)

    !Set default paths and surface name 
    allocate(character(len=3) :: cm3dopt%iopath)
    cm3dopt%iopath = 'io/'
    allocate(character(len=3) :: cm3dopt%optpath)
    cm3dopt%optpath = 'io/'
    allocate(character(len=23) :: cm3dopt%surfacename)
    cm3dopt%surfacename = 'cell_mesh3d_surface.dat'
elseif (nargs == 2) then !Use mode and surface name

    !Read mode 
    call get_command_argument(number=1, length=arglen)
    allocate(character(len=arglen) :: cm3dopt%mode)
    call get_command_argument(number=1, value=cm3dopt%mode, status=argstat)

    !Read surface name 
    call get_command_argument(number=2, length=arglen)
    allocate(character(len=arglen) :: cm3dopt%surfacename)
    call get_command_argument(number=2, value=cm3dopt%surfacename, status=argstat)

    !Set default paths
    allocate(character(len=3) :: cm3dopt%iopath)
    cm3dopt%iopath = 'io/'
    allocate(character(len=3) :: cm3dopt%optpath)
    cm3dopt%optpath = 'io/'
elseif (nargs == 3) then !Use mode, surface name and io path

    !Read mode 
    call get_command_argument(number=1, length=arglen)
    allocate(character(len=arglen) :: cm3dopt%mode)
    call get_command_argument(number=1, value=cm3dopt%mode, status=argstat)

    !Read surface name 
    call get_command_argument(number=2, length=arglen)
    allocate(character(len=arglen) :: cm3dopt%surfacename)
    call get_command_argument(number=2, value=cm3dopt%surfacename, status=argstat)

    !Read io path
    call get_command_argument(number=3, length=arglen)
    allocate(character(len=arglen) :: cm3dopt%iopath)
    call get_command_argument(number=3, value=cm3dopt%iopath, status=argstat)

    !Set default options path
    allocate(character(len=3) :: cm3dopt%optpath)
    cm3dopt%optpath = 'io/'
elseif (nargs == 4) then !Use mode, surface name, io path and options path

    !Read mode 
    call get_command_argument(number=1, length=arglen)
    allocate(character(len=arglen) :: cm3dopt%mode)
    call get_command_argument(number=1, value=cm3dopt%mode, status=argstat)

    !Read surface name 
    call get_command_argument(number=2, length=arglen)
    allocate(character(len=arglen) :: cm3dopt%surfacename)
    call get_command_argument(number=2, value=cm3dopt%surfacename, status=argstat)

    !Read io path
    call get_command_argument(number=3, length=arglen)
    allocate(character(len=arglen) :: cm3dopt%iopath)
    call get_command_argument(number=3, value=cm3dopt%iopath, status=argstat)

    !Read options path
    call get_command_argument(number=4, length=arglen)
    allocate(character(len=arglen) :: cm3dopt%optpath)
    call get_command_argument(number=4, value=cm3dopt%optpath, status=argstat)
else
    write(*,'(A)') '** too many command arguments supplied'
    stop
end if 
return 
end subroutine get_process_arguments




!Options import subroutine ===========================
subroutine cm3d_import_options(cm3dopt)
implicit none

!Variables - Import
type(cm3d_options) :: cm3dopt

!Variables - Local 
character(len=100) :: rtemp 

!Open file
open(11,file=cm3dopt%optpath//'cell_mesh3d_options.dat')

!Import options  
read(11,*) !skip
read(11,*) !skip
read(11,*) !skip
read(11,*) !skip
read(11,*) cm3dopt%dispt
read(11,*) !skip
read(11,*) !skip
read(11,*) cm3dopt%normDconv
read(11,*) !skip
read(11,*) !skip
read(11,*) !skip

read(11,*) cm3dopt%meshtype
read(11,*) !skip
read(11,*) !skip
read(11,*) rtemp
allocate(character(len=len_trim(rtemp)) :: cm3dopt%meshinout)
cm3dopt%meshinout = rtemp(1:len_trim(rtemp))
read(11,*) !skip
read(11,*) !skip
read(11,*) rtemp
allocate(character(len=len_trim(rtemp)) :: cm3dopt%surface_dir)
cm3dopt%surface_dir = rtemp(1:len_trim(rtemp))
read(11,*) !skip
read(11,*) !skip
read(11,*) rtemp
allocate(character(len=len_trim(rtemp)) :: cm3dopt%boundary_dir)
cm3dopt%boundary_dir = rtemp(1:len_trim(rtemp))
read(11,*) !skip
read(11,*) !skip
read(11,*) rtemp
allocate(character(len=len_trim(rtemp)) :: cm3dopt%meshfrmat)
cm3dopt%meshfrmat = rtemp(1:len_trim(rtemp))
read(11,*) !skip
read(11,*) !skip
read(11,*) !skip

read(11,*) cm3dopt%Nrefine
read(11,*) !skip
read(11,*) !skip
read(11,*) cm3dopt%NrefineB
read(11,*) !skip
read(11,*) !skip
read(11,*) cm3dopt%Ncell_max
read(11,*) !skip
read(11,*) !skip
read(11,*) cm3dopt%Nrefine_flood_i
read(11,*) !skip
read(11,*) !skip
read(11,*) cm3dopt%Nrefine_flood_f
read(11,*) !skip
read(11,*) !skip
read(11,*) cm3dopt%Nrefine_flood_B
read(11,*) !skip
read(11,*) !skip
read(11,*) cm3dopt%far_field_bound
read(11,*) !skip
read(11,*) !skip
read(11,*) cm3dopt%om_offset_x
read(11,*) cm3dopt%om_offset_y
read(11,*) cm3dopt%om_offset_z
read(11,*) !skip
read(11,*) !skip
read(11,*) !skip

read(11,*) cm3dopt%FminArea
read(11,*) !skip
read(11,*) !skip 
read(11,*) cm3dopt%CminVol
read(11,*) !skip
read(11,*) !skip
read(11,*) !skip

read(11,*) cm3dopt%NintEmax
read(11,*) !skip
read(11,*) !skip
read(11,*) cm3dopt%intcointol
read(11,*) !skip
read(11,*) !skip
read(11,*) cm3dopt%baryloctol
read(11,*) !skip
read(11,*) !skip
read(11,*) !skip

read(11,*) cm3dopt%surface_type
read(11,*) !skip
read(11,*) !skip 
read(11,*) cm3dopt%surf_force_simplify
read(11,*) !skip
read(11,*) !skip 
read(11,*) cm3dopt%surfRcurvM
read(11,*) !skip
read(11,*) !skip
read(11,*) !skip

read(11,*) cm3dopt%Nsstype
read(11,*) !skip
read(11,*) !skip
read(11,*) cm3dopt%nlpflood
read(11,*) !skip
read(11,*) !skip
read(11,*) cm3dopt%nlpsmooth
read(11,*) !skip
read(11,*) !skip
read(11,*) !skip

read(11,*) cm3dopt%ad_padding
read(11,*) !skip
read(11,*) !skip
read(11,*) cm3dopt%max_depth
read(11,*) !skip
read(11,*) !skip
read(11,*) !skip

read(11,*) cm3dopt%glink_con
read(11,*) !skip
read(11,*) !skip
read(11,*) cm3dopt%glink_nnn
read(11,*) !skip
read(11,*) !skip
read(11,*) cm3dopt%glink_nsmooth
read(11,*) !skip
read(11,*) !skip
read(11,*) cm3dopt%RBF_relax
read(11,*) !skip
read(11,*) !skip
read(11,*) !skip

read(11,*) cm3dopt%set_customBCs
read(11,*) !skip
read(11,*) !skip
read(11,*) cm3dopt%remFFzones
read(11,*) !skip
read(11,*) !skip
read(11,*) cm3dopt%remISzones
read(11,*) !skip
read(11,*) !skip
read(11,*) cm3dopt%bc_xmin
read(11,*) cm3dopt%bc_xmax
read(11,*) cm3dopt%bc_ymin
read(11,*) cm3dopt%bc_ymax
read(11,*) cm3dopt%bc_zmin
read(11,*) cm3dopt%bc_zmax

!Close file 
close(11)
return 
end subroutine cm3d_import_options




!Surface data import subroutine =========================== 
subroutine import_surface_geometry(surface_mesh,cm3dopt)
implicit none 

!Variables - Import
type(surface_data) :: surface_mesh
type(cm3d_options) :: cm3dopt

!Variables - Local
integer(in) :: ii

!Open file 
open(11,file=cm3dopt%iopath//cm3dopt%surfacename)

!Import item quantities
read(11,*) surface_mesh%nvtx,surface_mesh%nfcs

!Allocate
allocate(surface_mesh%vertices(surface_mesh%nvtx,3))
allocate(surface_mesh%connectivity(surface_mesh%nfcs,3))

!Import surface verticies
do ii=1,surface_mesh%nvtx
    read(11,*) surface_mesh%vertices(ii,:) !x || y || z
end do

!Import surface connectivity
do ii=1,surface_mesh%nfcs
    read(11,*) surface_mesh%connectivity(ii,:) !triangle v1 -> v2 -> v3
end do

!Close file
close(11)

!Display
if (cm3dopt%dispt == 1) then
    write(*,'(A)') '    {complete}'
end if
return 
end subroutine import_surface_geometry




!Custom boundary condition zone import subroutine ===========================
subroutine cm3d_import_customBC_zones(cm3dopt)
implicit none

!Variables - Import
type(cm3d_options) :: cm3dopt

!Variables - Local
integer(in) :: ii

!Read file
open(11,file=cm3dopt%optpath//'cell_mesh3d_bcond_zones.dat')
    read(11,*) cm3dopt%Nzone_cBC
    allocate(cm3dopt%BC_zone_bc(cm3dopt%Nzone_cBC))
    allocate(cm3dopt%BC_zone_coords(cm3dopt%Nzone_cBC,6))
    do ii=1,cm3dopt%Nzone_cBC
        read(11,*) cm3dopt%BC_zone_bc(ii),cm3dopt%BC_zone_coords(ii,:)
    end do 
close(11)
return 
end subroutine cm3d_import_customBC_zones




!Export volume mesh subroutines ===========================
subroutine export_volume_mesh(volume_mesh,cm3dopt)
implicit none 

!Variables - Import
type(cm3d_options) :: cm3dopt
type(vol_mesh_data) :: volume_mesh

!Variables - Local
integer(in) :: ii,vv

!Display
if (cm3dopt%dispt == 1) then
    write(*,'(A)') '--> writing mesh to file'
end if 

!Mesh 
open(11,file=cm3dopt%iopath//'grid') !mesh file write(11,'(*(I0,A))')
    write(11,'(I0,A,I0,A,I0,A,I0,A,I0)') volume_mesh%ncell,' ',volume_mesh%nface,' ',volume_mesh%nvtx,&
                                         ' ',volume_mesh%nface_mesh,' ',volume_mesh%nface_surface !Properties
    do ii=1,volume_mesh%nface
        write(11,'(I0)') volume_mesh%faces(ii)%nvtx
    end do 
    do ii=1,volume_mesh%nface
        do vv=1,volume_mesh%faces(ii)%nvtx
            write(11,'(I0,A)',advance="no") volume_mesh%faces(ii)%vertices(vv),' '
        end do 
        write(11,'(A)')
    end do 
    do ii=1,volume_mesh%nface
        write(11,'(I0,A,I0)') volume_mesh%faces(ii)%cleft,' ',volume_mesh%faces(ii)%cright
        ! write(11,'(I0,A,I0)') volume_mesh%faces(ii)%cright,' ',volume_mesh%faces(ii)%cleft
    end do 
    do ii=1,volume_mesh%nvtx
        write(11,'(F18.12,A,F18.12,A,F18.12)') volume_mesh%vtx(ii,1),' ',volume_mesh%vtx(ii,2),' ',volume_mesh%vtx(ii,3) !Vertices
    end do
close(11)

!Surface link data 
! open(11,file=cm3dopt%iopath//'grid_surf_link') !mesh surface to geometry surface links  
!     do ii=1,volume_mesh%nvtx_surf
!         write(11,'(I10,A,I10,A,F18.12)') volume_mesh%surf_vtx(ii),' ',volume_mesh%surf_vtx_seg(ii),' ',&
!                                          volume_mesh%surf_vtx_segfrac(ii) !vertex | link | fraction
!     end do
! close(11)

!Display
if (cm3dopt%dispt == 1) then
    write(*,'(A)') '    {complete}'
end if
return 
end subroutine export_volume_mesh

subroutine export_volume_mesh_flux(volume_mesh,cm3dopt)
implicit none 

!Variables - Import
type(cm3d_options) :: cm3dopt
type(vol_mesh_data) :: volume_mesh

!Variables - Local
integer(in) :: ii,vv

!Display
if (cm3dopt%dispt == 1) then
    write(*,'(A)') '--> writing mesh to file'
end if 

!Mesh 
open(11,file=cm3dopt%iopath//'grid_flux') !mesh file write(11,'(*(I0,A))')
    write(11,'(I0,A,I0,A,I0)') volume_mesh%ncell,' ',volume_mesh%nface,' ',volume_mesh%nvtx !Properties
    do ii=1,volume_mesh%nface
        write(11,'(I0,A)',advance="no") volume_mesh%faces(ii)%nvtx,' '
    end do 
    write(11,'(A)')
    do ii=1,volume_mesh%nface
        do vv=1,volume_mesh%faces(ii)%nvtx
            write(11,'(I0,A)',advance="no") volume_mesh%faces(ii)%vertices(vv),' '
        end do 
        write(11,'(A)')
    end do 
    do ii=1,volume_mesh%nface
        write(11,'(I0,A)',advance="no") volume_mesh%faces(ii)%cleft,' '
    end do
    write(11,'(A)')
    do ii=1,volume_mesh%nface
        write(11,'(I0,A)',advance="no") volume_mesh%faces(ii)%cright,' '
    end do
    write(11,'(A)')
    do ii=1,volume_mesh%nvtx
        write(11,'(F18.12,A,F18.12,A,F18.12)') volume_mesh%vtx(ii,1),' ',volume_mesh%vtx(ii,2),' ',volume_mesh%vtx(ii,3) !Vertices
    end do
close(11)

!Surface link data 
! open(11,file=cm3dopt%iopath//'grid_surf_link') !mesh surface to geometry surface links  
!     do ii=1,volume_mesh%nvtx_surf
!         write(11,'(I10,A,I10,A,F18.12)') volume_mesh%surf_vtx(ii),' ',volume_mesh%surf_vtx_seg(ii),' ',&
!                                          volume_mesh%surf_vtx_segfrac(ii) !vertex | link | fraction
!     end do
! close(11)

!Display
if (cm3dopt%dispt == 1) then
    write(*,'(A)') '    {complete}'
end if
return 
end subroutine export_volume_mesh_flux




!Import volume mesh subroutine ===========================
subroutine import_volume_mesh_flux(volume_mesh,cm3dopt)
implicit none 

!Variables - Import
type(cm3d_options) :: cm3dopt
type(vol_mesh_data) :: volume_mesh

!Variables - Local
integer(in) :: ii

!Open file 
open(11,file=cm3dopt%iopath//'grid_flux')

!Read item quantities 
read(11,*) volume_mesh%ncell,volume_mesh%nface,volume_mesh%nvtx

!Read mesh faces 
allocate(volume_mesh%faces(volume_mesh%nface))
read(11,*) volume_mesh%faces(:)%nvtx
do ii=1,volume_mesh%nface
    allocate(volume_mesh%faces(ii)%vertices(volume_mesh%faces(ii)%nvtx))
    read(11,*) volume_mesh%faces(ii)%vertices(:)
end do 
read(11,*) volume_mesh%faces(:)%cleft
read(11,*) volume_mesh%faces(:)%cright

!Read mesh vertices 
allocate(volume_mesh%vtx(volume_mesh%nvtx,3))
do ii=1,volume_mesh%nvtx
    read(11,*) volume_mesh%vtx(ii,:)
end do 

!Close file 
close(11)

!Display
if (cm3dopt%dispt == 1) then
    write(*,'(A)') '    {complete}'
end if
return 
end subroutine import_volume_mesh_flux



!Export PPoU volume-surface to surface interpolation structure ===========================
subroutine export_vs2s_interpstruc(volume_mesh,surface_mesh,cm3dopt)
implicit none 

!Variables - Import
type(cm3d_options) :: cm3dopt
type(vol_mesh_data) :: volume_mesh
type(surface_data) :: surface_mesh

!Variables - Local 
integer(in) :: vv,ii,jj

!Display
if (cm3dopt%dispt == 1) then
    write(*,'(A)') '--> writing volume-surface to surface interpolation structure to file'
end if 

!Open 
open(11,file=cm3dopt%iopath//'vs2s_interp') 

!Write data
write(11,'(I0,A,I0,A,I0)') surface_mesh%nvtx,' ',cm3dopt%glink_nnn,' ',cm3dopt%glink_nsmooth+1 !number of surface vertices | number of interpolation points per surface vertex | number of smoothing points per surface vertex
write(11,'(A)') ' '

!Write each interpolation structure for each surface vertex
do vv=1,surface_mesh%nvtx
    write(11,'(I0)') vv !surface point index 
    write(11,'(I0,A,I0)') volume_mesh%vsinterp(vv)%npnts_vi,' ',volume_mesh%vsinterp(vv)%npnts_ss !number of unique vertices in the interpolation and smoothing structures
    do ii=1,cm3dopt%glink_nsmooth+1 !write smoothing point list for this surface point 
        write(11,'(I0,A)',advance='no') volume_mesh%vsinterp(vv)%surf_smooth_pnts(ii),' '
    end do 
    write(11,'(A)') !skip to next line 
    do ii=1,cm3dopt%glink_nsmooth+1 !write RBF dependance for this surface point on each surface smoothing point
        write(11,'(A,A)',advance='no') real2F0_Xstring(volume_mesh%vsinterp(vv)%surf_smoothRBF(ii),12_in),' '
    end do 
    write(11,'(A)') !skip to next line 
    do ii=1,cm3dopt%glink_nnn !write volume point list for this surface point 
        write(11,'(I0,A)',advance='no') volume_mesh%vsinterp(vv)%vol_pnts(ii),' '
    end do 
    write(11,'(A)') !skip to next line 
    do ii=1,cm3dopt%glink_nnn !write RBF dependance for this surface point on each volume point 
        write(11,'(A,A)',advance='no') real2F0_Xstring(volume_mesh%vsinterp(vv)%surf2volRBF(ii),12_in),' '
    end do 
    write(11,'(A)') !skip to next line 
    do ii=1,cm3dopt%glink_nnn !write interpolation matrix 
        do jj=1,cm3dopt%glink_nnn 
            write(11,'(A,A)',advance='no') real2F0_Xstring(volume_mesh%vsinterp(vv)%Ri(ii,jj),12_in),' '
        end do 
        write(11,'(A)') !skip to next line
    end do 
    write(11,'(A)') ' '
end do 

!Close
close(11)

!Display
if (cm3dopt%dispt == 1) then
    write(*,'(A)') '    {complete}'
end if
return 
end subroutine export_vs2s_interpstruc




!Import flow gradients subroutine ===========================
subroutine import_flow_gradients(gradient_vol,volume_mesh,cm3dopt)
implicit none 

!Variables - Import
real(dp), dimension(:,:), allocatable :: gradient_vol
type(cm3d_options) :: cm3dopt
type(vol_mesh_data) :: volume_mesh

!Variables - Local
integer(in) :: ii

!Allocate gradient array 
allocate(gradient_vol(volume_mesh%nvtx,3))

!Open file 
open(11,file=cm3dopt%iopath//'gradient.dat')

!Read gradients 
do ii=1,volume_mesh%nvtx
    read(11,*) gradient_vol(ii,:)
end do 

!Close file 
close(11)

!Display
if (cm3dopt%dispt == 1) then
    write(*,'(A)') '    {complete}'
end if
return 
end subroutine import_flow_gradients




!Export surface gradients subroutine ===========================
subroutine export_surface_gradients(gradient_surf,surface_mesh,cm3dopt)
implicit none 

!Variables - Import
real(dp), dimension(:,:) :: gradient_surf
type(cm3d_options) :: cm3dopt
type(surface_data) :: surface_mesh

!Variables - Local
integer(in) :: ii

!Open file 
open(11,file=cm3dopt%iopath//'gradient_surf.dat')

!Write gradients 
do ii=1,surface_mesh%nvtx
    write(11,'(A,A,A,A,A)') real2F0_Xstring(gradient_surf(ii,1),12_in),' ',&
    real2F0_Xstring(gradient_surf(ii,2),12_in),' ',real2F0_Xstring(gradient_surf(ii,3),12_in)
end do 

!Close file 
close(11)

!Display
if (cm3dopt%dispt == 1) then
    write(*,'(A)') '    {complete}'
end if
return 
end subroutine export_surface_gradients




!Write status subroutine ===========================
subroutine export_status(cm3dopt,cm3dfailure)
implicit none 

!Variables - Import
integer(in) :: cm3dfailure
type(cm3d_options) :: cm3dopt

!Write
open(11,file=cm3dopt%iopath//'cm3d_status') 
    write(11,'(i2)') cm3dfailure
close(11)
return 
end subroutine export_status




!Export TECPLOT mesh file with cell data subroutine =========================  
subroutine write_cell_dataPLT(filename,volume_mesh,cell_data)
use ieee_arithmetic
implicit none 


!Variables - Import
real(dp), dimension(:) :: cell_data
character(*), intent(in) :: filename
type(vol_mesh_data) :: volume_mesh

!Variables - Local 
integer(in) :: fh,i,nperline

!Open file
fh = 11
open(fh,file=filename//'.plt',status='unknown')

!TECPLOT formatting
write(fh,'(A)',advance="no") 'VARIABLES="X" "Y" "data"'
write(fh,*) 'ZONE T="CellData"'
write(fh,'(A)',advance="no") 'VARLOCATION=([1,2]=NODAL,[3]=CELLCENTERED)'
write(fh,*) 'ZONETYPE=FEPOLYGON'
write(fh,'(A,I8)') ' Nodes=',volume_mesh%nvtx
write(fh,'(A,I8)') ' Elements=',volume_mesh%ncell
write(fh,'(A,I8)') ' Faces=',volume_mesh%nedge
write(fh,*) 'NumConnectedBoundaryFaces=0 '
write(fh,*) 'TotalNumBoundaryConnections=0 '

! These loops are because tecplot has a maximum number of characters per line
nperline = 100
write(fh,*) ( volume_mesh%vtx(i:min(i+nperline-1,volume_mesh%nvtx),1),NEW_LINE('A') , i=1,volume_mesh%nvtx,nperline )
write(fh,*) ( volume_mesh%vtx(i:min(i+nperline-1,volume_mesh%nvtx),2),NEW_LINE('A') , i=1,volume_mesh%nvtx,nperline )
write(fh,*) ( cell_data(i:min(i+nperline-1,volume_mesh%ncell)),NEW_LINE('A') , i=1,volume_mesh%ncell,nperline )
do i=1,volume_mesh%nedge
    write(fh,*) volume_mesh%edge(i,1),volume_mesh%edge(i,2)  
end do
write(fh,*) ( max(0,volume_mesh%edge(i:min(i+nperline-1,volume_mesh%nedge),3)),NEW_LINE('A') , i=1,volume_mesh%nedge,nperline )
write(fh,*) ( max(0,volume_mesh%edge(i:min(i+nperline-1,volume_mesh%nedge),4)),NEW_LINE('A') , i=1,volume_mesh%nedge,nperline )

!Close file
close(fh)
return 
end subroutine write_cell_dataPLT




!F0.X format with leading zero function =========================
function real2F0_Xstring(val,X) result(str)

!Result 
character(len=:), allocatable :: str

!Variables - Import 
character(len=10) :: frmtI
character(len=:), allocatable :: frmt,str_I
integer(in) :: X,len_frmt,len_str
real(dp) :: val

!Set format descriptor
write(frmtI,'(I0)') X
len_frmt = len_trim(frmtI)
allocate(character(len=len_frmt) :: frmt)
frmt = frmtI(1:len_frmt)
frmt = frmt//')'
frmt = '(F0.'//frmt

!Allocate initial character
allocate(character(len=4*X) :: str_I)

!Write data to return charachter
write(str_I,frmt) val

!Allocate return character
len_str = len_trim(str_I)
allocate(character(len=len_str) :: str)
str = str_I(1:len_str)

!Assign leading zero if required
if (str(1:1) == '.') then 
    str = '0'//str
elseif (str(1:2) == '-.') then 
    str = '-0.'//str(3:len_str)
end if 
end function real2F0_Xstring


end module cellmesh3d_io_mod