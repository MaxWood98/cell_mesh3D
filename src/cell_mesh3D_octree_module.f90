!cell_mesh2d quadtree module
!Max Wood - mw16116@bristol.ac.uk
!Univeristy of Bristol - Department of Aerospace Engineering

!Version 3.0
!Updated 18-10-2023

!Module
module cell_mesh3d_octree_mod
use cellmesh3d_adtree_mod
use cellmesh3d_geometry_mod
contains 


!Octree mesh refinement subroutine ===========================
subroutine octree_mesh_refine(ot_mesh,surface_mesh,surface_adtree,cm3dopt,obj_cx,obj_cy,obj_cz,cm3dfailure)
implicit none 

!Variables - Import
integer(in) :: cm3dfailure
real(dp) :: obj_cx,obj_cy,obj_cz
type(cm3d_options) :: cm3dopt
type(octree_data) :: ot_mesh
type(surface_data) :: surface_mesh
type(tree_data) :: surface_adtree

!Variables - Local 
integer(in) :: aa,ff,cc,ccf,rr
integer(in) :: cins,vins,nselected,nrefineN,nrefine,cidx,cadj,Nflood,ncexit,reflevel
integer(in) :: Nflood_refiter(cm3dopt%Nrefine),node_select(surface_adtree%nnode)
integer(in) :: Qrefine(cm3dopt%Ncell_max),QfloodR(cm3dopt%Ncell_max),QfloodRN(cm3dopt%Ncell_max)
integer(in) :: fopposite(6),faces(6,4),edges(12,2),face2edge(6,4),face2opedge(6,4),ccadj(8,6)
integer(in) :: ediagonal(12),ccadj_diag(8,12),ccadj_diag_dir(8,12)
real(dp) :: Nflood_crr,Nflood_dcr,cpadSZ
real(dp) :: zxmin,zxmax,zymin,zymax,zzmin,zzmax

!Initialise failure tag
cm3dfailure = 0 

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

!Define face2opedge
face2opedge(1,1) = 8
face2opedge(1,2) = 11
face2opedge(1,3) = 7
face2opedge(1,4) = 3
face2opedge(2,1) = 6
face2opedge(2,2) = 9
face2opedge(2,3) = 5
face2opedge(2,4) = 1
face2opedge(3,1) = 5
face2opedge(3,2) = 12
face2opedge(3,3) = 8
face2opedge(3,4) = 4
face2opedge(4,1) = 7
face2opedge(4,2) = 10
face2opedge(4,3) = 6
face2opedge(4,4) = 2
face2opedge(5,1) = 11
face2opedge(5,2) = 12
face2opedge(5,3) = 9
face2opedge(5,4) = 10
face2opedge(6,1) = 3
face2opedge(6,2) = 2
face2opedge(6,3) = 1
face2opedge(6,4) = 4

!Define cell child-child adjacency
ccadj(1,1) = -4
ccadj(1,2) = 4
ccadj(1,3) = 2
ccadj(1,4) = -2
ccadj(1,5) = -5
ccadj(1,6) = 5
ccadj(2,1) = -3
ccadj(2,2) = 3
ccadj(2,3) = -1
ccadj(2,4) = 1
ccadj(2,5) = -6
ccadj(2,6) = 6
ccadj(3,1) = 2
ccadj(3,2) = -2
ccadj(3,3) = -4
ccadj(3,4) = 4
ccadj(3,5) = -7
ccadj(3,6) = 7
ccadj(4,1) = 1
ccadj(4,2) = -1
ccadj(4,3) = 3
ccadj(4,4) = -3
ccadj(4,5) = -8
ccadj(4,6) = 8
ccadj(5,1) = -8
ccadj(5,2) = 8
ccadj(5,3) = 6
ccadj(5,4) = -6
ccadj(5,5) = 1
ccadj(5,6) = -1
ccadj(6,1) = -7
ccadj(6,2) = 7
ccadj(6,3) = -5
ccadj(6,4) = 5
ccadj(6,5) = 2
ccadj(6,6) = -2
ccadj(7,1) = 6
ccadj(7,2) = -6
ccadj(7,3) = -8
ccadj(7,4) = 8
ccadj(7,5) = 3
ccadj(7,6) = -3
ccadj(8,1) = 5
ccadj(8,2) = -5
ccadj(8,3) = 7
ccadj(8,4) = -7
ccadj(8,5) = 4
ccadj(8,6) = -4

!Define cell child-child diagonal adjacency
ccadj_diag(1,:) = (/8,6,8,6,3,3,3,3,8,6,8,6/)
ccadj_diag(2,:) = (/7,5,7,5,4,4,4,4,7,5,7,5/)
ccadj_diag(3,:) = (/6,8,6,8,1,1,1,1,6,8,6,8/)
ccadj_diag(4,:) = (/5,7,5,7,2,2,2,2,5,7,5,7/)
ccadj_diag(5,:) = (/4,2,4,2,7,7,7,7,4,2,4,2/)
ccadj_diag(6,:) = (/3,1,3,1,8,8,8,8,3,1,3,1/)
ccadj_diag(7,:) = (/2,4,2,4,5,5,5,5,2,4,2,4/)
ccadj_diag(8,:) = (/1,3,1,3,6,6,6,6,1,3,1,3/)

!Define cell child-child diagonal adjacency search location (0 = internal | -1 = diagonal on same edge | positive number x = in cell_adjacent through face x of parent)
ccadj_diag_dir(1,:) = (/-1,5,5,-1,-1,1,0,4,1,0,0,4/)
ccadj_diag_dir(2,:) = (/-1,-1,5,5,1,-1,3,0,1,3,0,0/)
ccadj_diag_dir(3,:) = (/5,-1,-1,5,0,3,-1,2,0,3,2,0/)
ccadj_diag_dir(4,:) = (/5,5,-1,-1,4,0,2,-1,0,0,2,4/)
ccadj_diag_dir(5,:) = (/1,0,0,4,-1,1,0,4,-1,6,6,-1/)
ccadj_diag_dir(6,:) = (/1,3,0,0,1,-1,3,0,-1,-1,6,6/)
ccadj_diag_dir(7,:) = (/0,3,2,0,0,3,-1,2,6,-1,-1,6/)
ccadj_diag_dir(8,:) = (/0,0,2,4,4,0,2,-1,6,6,-1,-1/)

!Initialise refinement arrays
Qrefine(:) = 0
QfloodR(:) = 0
QfloodRN(:) = 0
nrefine = 0 
nrefineN = 0

!Allocate octree mesh arrays
allocate(ot_mesh%cell_level(cm3dopt%Ncell_max))
allocate(ot_mesh%cell_parent(cm3dopt%Ncell_max))
allocate(ot_mesh%cell_adjacent(cm3dopt%Ncell_max,6)) !-2 = far field || -1 = object boundary || else cell index
allocate(ot_mesh%cell_diagonal(cm3dopt%Ncell_max,12)) !-2 = far field || -1 = object boundary || else cell index
allocate(ot_mesh%cell_vcnr(cm3dopt%Ncell_max,8))
allocate(ot_mesh%cell_vcmid(cm3dopt%Ncell_max))
allocate(ot_mesh%cell_vfmid(cm3dopt%Ncell_max,6))
allocate(ot_mesh%cell_vemid(cm3dopt%Ncell_max,12))
allocate(ot_mesh%cell_child(cm3dopt%Ncell_max,8))
allocate(ot_mesh%vtx(2*cm3dopt%Ncell_max,3))
ot_mesh%cell_level(:) = 0
ot_mesh%cell_parent(:) = 0 
ot_mesh%cell_adjacent(:,:) = 0
ot_mesh%cell_diagonal(:,:) = 0
ot_mesh%cell_vcnr(:,:) = 0
ot_mesh%cell_vcmid(:) = 0
ot_mesh%cell_vfmid(:,:) = 0
ot_mesh%cell_vemid(:,:) = 0
ot_mesh%cell_child(:,:) = 0
ot_mesh%vtx(:,:) = 0.0d0

!Set initial item counter values
cins = 2 !cell
vins = 9 !vertex 

!Construct first cell to define global mesh domain
ot_mesh%vtx(1,1) = obj_cx - cm3dopt%far_field_bound + cm3dopt%om_offset_x
ot_mesh%vtx(1,2) = obj_cy - cm3dopt%far_field_bound + cm3dopt%om_offset_y
ot_mesh%vtx(1,3) = obj_cz - cm3dopt%far_field_bound + cm3dopt%om_offset_z
ot_mesh%vtx(2,1) = obj_cx - cm3dopt%far_field_bound + cm3dopt%om_offset_x
ot_mesh%vtx(2,2) = obj_cy + cm3dopt%far_field_bound + cm3dopt%om_offset_y
ot_mesh%vtx(2,3) = obj_cz - cm3dopt%far_field_bound + cm3dopt%om_offset_z
ot_mesh%vtx(3,1) = obj_cx + cm3dopt%far_field_bound + cm3dopt%om_offset_x
ot_mesh%vtx(3,2) = obj_cy + cm3dopt%far_field_bound + cm3dopt%om_offset_y
ot_mesh%vtx(3,3) = obj_cz - cm3dopt%far_field_bound + cm3dopt%om_offset_z
ot_mesh%vtx(4,1) = obj_cx + cm3dopt%far_field_bound + cm3dopt%om_offset_x
ot_mesh%vtx(4,2) = obj_cy - cm3dopt%far_field_bound + cm3dopt%om_offset_y
ot_mesh%vtx(4,3) = obj_cz - cm3dopt%far_field_bound + cm3dopt%om_offset_z
ot_mesh%vtx(5,1) = obj_cx - cm3dopt%far_field_bound + cm3dopt%om_offset_x
ot_mesh%vtx(5,2) = obj_cy - cm3dopt%far_field_bound + cm3dopt%om_offset_y
ot_mesh%vtx(5,3) = obj_cz + cm3dopt%far_field_bound + cm3dopt%om_offset_z
ot_mesh%vtx(6,1) = obj_cx - cm3dopt%far_field_bound + cm3dopt%om_offset_x
ot_mesh%vtx(6,2) = obj_cy + cm3dopt%far_field_bound + cm3dopt%om_offset_y
ot_mesh%vtx(6,3) = obj_cz + cm3dopt%far_field_bound + cm3dopt%om_offset_z
ot_mesh%vtx(7,1) = obj_cx + cm3dopt%far_field_bound + cm3dopt%om_offset_x
ot_mesh%vtx(7,2) = obj_cy + cm3dopt%far_field_bound + cm3dopt%om_offset_y
ot_mesh%vtx(7,3) = obj_cz + cm3dopt%far_field_bound + cm3dopt%om_offset_z
ot_mesh%vtx(8,1) = obj_cx + cm3dopt%far_field_bound + cm3dopt%om_offset_x
ot_mesh%vtx(8,2) = obj_cy - cm3dopt%far_field_bound + cm3dopt%om_offset_y
ot_mesh%vtx(8,3) = obj_cz + cm3dopt%far_field_bound + cm3dopt%om_offset_z
ot_mesh%cell_vcnr(1,1) = 1
ot_mesh%cell_vcnr(1,2) = 2
ot_mesh%cell_vcnr(1,3) = 3
ot_mesh%cell_vcnr(1,4) = 4
ot_mesh%cell_vcnr(1,5) = 5
ot_mesh%cell_vcnr(1,6) = 6
ot_mesh%cell_vcnr(1,7) = 7
ot_mesh%cell_vcnr(1,8) = 8
ot_mesh%cell_adjacent(1,1) = -2
ot_mesh%cell_adjacent(1,2) = -2
ot_mesh%cell_adjacent(1,3) = -2
ot_mesh%cell_adjacent(1,4) = -2
ot_mesh%cell_adjacent(1,5) = -2
ot_mesh%cell_adjacent(1,6) = -2
ot_mesh%cell_diagonal(1,1) = -2
ot_mesh%cell_diagonal(1,2) = -2
ot_mesh%cell_diagonal(1,3) = -2
ot_mesh%cell_diagonal(1,4) = -2
ot_mesh%cell_diagonal(1,5) = -2
ot_mesh%cell_diagonal(1,6) = -2
ot_mesh%cell_diagonal(1,7) = -2
ot_mesh%cell_diagonal(1,8) = -2
ot_mesh%cell_diagonal(1,9) = -2
ot_mesh%cell_diagonal(1,10) = -2
ot_mesh%cell_diagonal(1,11) = -2
ot_mesh%cell_diagonal(1,12) = -2
ot_mesh%cell_level(1) = 1
ot_mesh%cell_parent(1) = 0 

!Store the global mesh domain 
cm3dopt%mesh_xmin = minval(ot_mesh%vtx(1:8,1))
cm3dopt%mesh_xmax = maxval(ot_mesh%vtx(1:8,1))
cm3dopt%mesh_ymin = minval(ot_mesh%vtx(1:8,2))
cm3dopt%mesh_ymax = maxval(ot_mesh%vtx(1:8,2))
cm3dopt%mesh_zmin = minval(ot_mesh%vtx(1:8,3))
cm3dopt%mesh_zmax = maxval(ot_mesh%vtx(1:8,3))

!Set custom boundary conditions on base cell faces
if (cm3dopt%set_customBCs == 1) then 
    ot_mesh%cell_adjacent(1,1) = cm3dopt%bc_xmin
    ot_mesh%cell_adjacent(1,2) = cm3dopt%bc_xmax
    ot_mesh%cell_adjacent(1,3) = cm3dopt%bc_ymax
    ot_mesh%cell_adjacent(1,4) = cm3dopt%bc_ymin
    ot_mesh%cell_adjacent(1,5) = cm3dopt%bc_zmin
    ot_mesh%cell_adjacent(1,6) = cm3dopt%bc_zmax
end if 

!Construct refinement adjacency flooding value for each iteration 
Nflood_dcr = (real(cm3dopt%Nrefine_flood_i,dp) - real(cm3dopt%Nrefine_flood_f,dp))/real(cm3dopt%Nrefine-1,dp)
Nflood_refiter(:) = 0
Nflood_refiter(1) = cm3dopt%Nrefine_flood_i
Nflood_crr = real(cm3dopt%Nrefine_flood_i,dp)
do rr=2,cm3dopt%Nrefine
    Nflood_crr = Nflood_crr - Nflood_dcr
    Nflood_refiter(rr) = nint(Nflood_crr,in)
end do
Nflood_refiter(cm3dopt%Nrefine) = cm3dopt%Nrefine_flood_f

!Construct octree mesh -------------------------------------------------------
ncexit = 0 
if (cm3dopt%dispt == 1) then
    write(*,'(A)') '--> refining cells: '
    write(*,'(A)') '  level  |  cells  |  Nflood'
    write(*,'(A)') ' ----------------------------'
    write(*,'(A,I2,A,I8,A,I4)') '    ',0,'   ',1,'     ',Nflood_refiter(1)
end if
reflevel = 0 
do rr=1,cm3dopt%Nrefine + cm3dopt%NrefineB

    !Increment refinement level 
    reflevel = reflevel + 1

    !Set padding range base size -> size of cell at this level 
    cpadSZ = 2.0d0*cm3dopt%far_field_bound/(2.0d0**(rr - 1))

    !Identify cells to refine 
    Qrefine(1:cins-1) = 0
    QfloodR(1:nrefine) = 0 
    nrefine = 0
    do cc=1,cins-1
        if (ot_mesh%cell_level(cc) == rr) then !At current maximum level 
            if (ot_mesh%cell_child(cc,1) == 0) then 

                !Intersection bounding box
                zxmin = ot_mesh%vtx(ot_mesh%cell_vcnr(cc,1),1) - cpadSZ*cm3dopt%ADTpadding !tgt bounding box -> xmin
                zxmax = ot_mesh%vtx(ot_mesh%cell_vcnr(cc,7),1) + cpadSZ*cm3dopt%ADTpadding !tgt bounding box -> xmax
                zymin = ot_mesh%vtx(ot_mesh%cell_vcnr(cc,1),2) - cpadSZ*cm3dopt%ADTpadding !tgt bounding box -> ymin
                zymax = ot_mesh%vtx(ot_mesh%cell_vcnr(cc,7),2) + cpadSZ*cm3dopt%ADTpadding !tgt bounding box -> ymax
                zzmin = ot_mesh%vtx(ot_mesh%cell_vcnr(cc,1),3) - cpadSZ*cm3dopt%ADTpadding !tgt bounding box -> zmin
                zzmax = ot_mesh%vtx(ot_mesh%cell_vcnr(cc,7),3) + cpadSZ*cm3dopt%ADTpadding !tgt bounding box -> zmax

                !Identify any triangle bounding boxes that may overlap the cell 
                call search_ADtree(nselected,node_select,surface_adtree,zxmin,zxmax,zymin,zymax,zzmin,zzmax)

                !Tag cell to refine 
                if (nselected .NE. 0) then 
                    if (rr .LE. cm3dopt%Nrefine) then !Check for and tag geometry surface overlap (normal refinement)
                        Qrefine(cc) = face_cell_ovlp_bool(ot_mesh,surface_mesh,surface_adtree,node_select,nselected,cc)
                    else !Check for and tag geometry surface overlap with smaller than cell radius of curvature (boosted refinement)
                        Qrefine(cc) = face_cell_ovlp_wrcurv_bool(ot_mesh,surface_mesh,surface_adtree,node_select,nselected,cc)
                    end if 
                    if (Qrefine(cc) == 1) then 
                        nrefine = nrefine + 1
                        QfloodR(nrefine) = cc
                    end if
                end if 
            end if 
        end if 
    end do 

    !Exit with no refinement requested
    if (nrefine == 0) then 
        if (cm3dopt%dispt == 1) then
            write(*,'(A)')  '    {refinement completed as no cells tagged}'
        end if
        ncexit = 1 
        exit 
    end if 

    !Flood refinement status 
    if (rr .LE. cm3dopt%Nrefine) then 
        Nflood = Nflood_refiter(rr)
    else
        Nflood = cm3dopt%Nrefine_flood_b 
    end if 
    do ff=1,Nflood
        nrefineN = 0 
        do ccf=1,nrefine
            cidx = QfloodR(ccf)
            do aa=1,6 
                cadj = ot_mesh%cell_adjacent(cidx,aa)
                if (cadj .GT. 0) then 
                    if (ot_mesh%cell_child(cadj,1) == 0) then 
                        if (Qrefine(cadj) == 0) then 
                            nrefineN = nrefineN + 1
                            QfloodRN(nrefineN) = cadj
                            Qrefine(cadj) = 1
                        end if 
                    end if 
                end if 
            end do 
        end do 
        QfloodR(1:nrefineN) = QfloodRN(1:nrefineN)
        QfloodRN(1:nrefineN) = 0
        nrefine = nrefineN
    end do 

    !Refine tagged cells 
    do cc=1,cins-1
        if (Qrefine(cc) == 1) then
            call octree_cell_refine(ot_mesh,cm3dopt,cm3dfailure,fopposite,faces,edges,face2edge,face2opedge,ediagonal,&
                                    ccadj_diag,ccadj_diag_dir,ccadj,vins,cins,cc)
            if (cm3dfailure == 1) then 
                return 
            end if 
        end if 
    end do 

    !Ensure all vertex adjacencies are updated
    ! do cc=1,cins-1

    !     !Pull adjacent face midpoint vertices 
    !     do aa=1,6 
    !         cadj = ot_mesh%cell_adjacent(cc,aa)
    !         if (ot_mesh%cell_vfmid(cc,aa) == 0) then 
    !             if (cadj .GT. 0) then 
    !                 if (ot_mesh%cell_level(cadj) == ot_mesh%cell_level(cc)) then 
    !                     ot_mesh%cell_vfmid(cc,aa) = ot_mesh%cell_vfmid(cadj,fopposite(aa))
    !                 end if 
    !             end if
    !         end if 
    !     end do 

    !     !Pull adjacent edge midpoint vertices 
    !     do ff=1,6
    !         cadj = ot_mesh%cell_adjacent(cc,ff)
    !         if (cadj .GT. 0) then 
    !             if (ot_mesh%cell_level(cadj) == ot_mesh%cell_level(cc)) then 
    !                 do aa=1,4
    !                     if ((ot_mesh%cell_vemid(cc,face2edge(ff,aa)) == 0) .AND. &
    !                     (ot_mesh%cell_vemid(cadj,face2opedge(ff,aa)) .NE. 0)) then 
    !                         ot_mesh%cell_vemid(cc,face2edge(ff,aa)) = ot_mesh%cell_vemid(cadj,face2opedge(ff,aa))
    !                     end if 
    !                 end do 
    !             end if 
    !         end if 
    !     end do 

    !     !Push face midpoint vertices
    !     do aa=1,6 
    !         cadj = ot_mesh%cell_adjacent(cc,aa)
    !         if (ot_mesh%cell_vfmid(cc,aa) == 0) then 
    !             if (cadj .GT. 0) then 
    !                 if (ot_mesh%cell_level(cadj) == ot_mesh%cell_level(cc)) then 
    !                     if ((ot_mesh%cell_vfmid(cadj,fopposite(aa)) == 0) .AND. (ot_mesh%cell_vfmid(cc,aa) .NE. 0)) then 
    !                         ot_mesh%cell_vfmid(cadj,fopposite(aa)) = ot_mesh%cell_vfmid(cc,aa)
    !                     end if 
    !                 end if 
    !             end if
    !         end if 
    !     end do 

    !     !Push edge midpoint vertices
    !     do ff=1,6
    !         cadj = ot_mesh%cell_adjacent(cc,ff)
    !         if (cadj .GT. 0) then 
    !             if (ot_mesh%cell_level(cadj) == ot_mesh%cell_level(cc)) then 
    !                 do aa=1,4
    !                     if ((ot_mesh%cell_vemid(cc,face2edge(ff,aa)) .NE. 0) .AND. &
    !                     (ot_mesh%cell_vemid(cadj,face2opedge(ff,aa)) == 0)) then 
    !                         ot_mesh%cell_vemid(cadj,face2opedge(ff,aa)) = ot_mesh%cell_vemid(cc,face2edge(ff,aa))
    !                     end if 
    !                 end do 
    !             end if 
    !         end if 
    !     end do 
    ! end do 

    !Display
    if (cm3dopt%dispt == 1) then
        write(*,'(A,I2,A,I8,A,I4)') '    ',rr,'   ',cins-1,'     ',Nflood
    end if

    !Display completion of base refinement 
    if (rr == cm3dopt%Nrefine) then 
        if (cm3dopt%dispt == 1) then
            write(*,'(A)')  '    {base refinement completed}'
        end if
    end if 
end do 
if ((rr == cm3dopt%Nrefine + cm3dopt%NrefineB + 1) .AND. (ncexit == 0)) then 
    if (cm3dopt%dispt == 1) then
        write(*,'(A)')  '    {refinement completed at maximum level}'
    end if
end if 

!Assign refined item quantities
ot_mesh%cins = cins 
ot_mesh%vins = vins 

!Initialise cell_keep to tag retained cells 
allocate(ot_mesh%cell_keep(ot_mesh%cins-1))
ot_mesh%cell_keep(:) = 0 
do cc=1,ot_mesh%cins-1
    if (ot_mesh%cell_child(cc,1) == 0) then !Has no children
        ot_mesh%cell_keep(cc) = 1
    end if 
end do

!Debug =============
! !write vertices
! open(11,file='io/vertices.dat')
! do cc=1,ot_mesh%vins-1
!     write(11,*) ot_mesh%vtx(cc,:)
! end do 
! close(11)
return 
end subroutine octree_mesh_refine




!Octree refinement subroutine ===========================
subroutine octree_cell_refine(ot_mesh,cm3dopt,cm3dfailure,fopposite,faces,edges,face2edge,face2opedge,ediagonal,&
                              ccadj_diag,ccadj_diag_dir,ccadj,vins,cins,ctgt)
implicit none 

!Variables - Import
integer(in) :: cm3dfailure,vins,cins,ctgt 
integer(in) :: fopposite(6),ediagonal(12),faces(6,4),edges(12,2),face2edge(6,4),face2opedge(6,4),ccadj(8,6)
integer(in) :: ccadj_diag(8,12),ccadj_diag_dir(8,12)
type(octree_data) :: ot_mesh
type(cm3d_options) :: cm3dopt

!Variables - Local 
integer(in) :: aa,ff,ccf,ee
integer(in) :: vtxCC,cadj,child_idx

!Build cell centre vertex 
vtxCC = vins 
ot_mesh%cell_vcmid(ctgt) = vins
ot_mesh%vtx(vins,1) = 0.125d0*sum(ot_mesh%vtx(ot_mesh%cell_vcnr(ctgt,:),1))
ot_mesh%vtx(vins,2) = 0.125d0*sum(ot_mesh%vtx(ot_mesh%cell_vcnr(ctgt,:),2))
ot_mesh%vtx(vins,3) = 0.125d0*sum(ot_mesh%vtx(ot_mesh%cell_vcnr(ctgt,:),3))
vins = vins + 1
if (vins .GE. cm3dopt%Ncell_max) then
    cm3dfailure = 1
    write(*,'(A)') '** maximum vertex count exceeded'
    return 
end if

!Pull adjacent face midpoint vertices 
do aa=1,6 
    cadj = ot_mesh%cell_adjacent(ctgt,aa)
    if (ot_mesh%cell_vfmid(ctgt,aa) == 0) then 
        if (cadj .GT. 0) then 
            if (ot_mesh%cell_level(cadj) == ot_mesh%cell_level(ctgt)) then 
                ot_mesh%cell_vfmid(ctgt,aa) = ot_mesh%cell_vfmid(cadj,fopposite(aa))
            end if 
        end if
    end if 
end do 

!Pull adjacent edge midpoint vertices 
do ff=1,6 !though faces
    cadj = ot_mesh%cell_adjacent(ctgt,ff)
    if (cadj .GT. 0) then 
        if (ot_mesh%cell_level(cadj) == ot_mesh%cell_level(ctgt)) then 
            do aa=1,4
                if ((ot_mesh%cell_vemid(ctgt,face2edge(ff,aa)) == 0) .AND. &
                (ot_mesh%cell_vemid(cadj,face2opedge(ff,aa)) .NE. 0)) then 
                    ot_mesh%cell_vemid(ctgt,face2edge(ff,aa)) = ot_mesh%cell_vemid(cadj,face2opedge(ff,aa))
                end if 
            end do 
        end if 
    end if 
end do 
do ee=1,12 !though edges
    cadj = ot_mesh%cell_diagonal(ctgt,ee)
    if (cadj .GT. 0) then 
        if (ot_mesh%cell_level(cadj) == ot_mesh%cell_level(ctgt)) then 
            if ((ot_mesh%cell_vemid(ctgt,ee) == 0) .AND. (ot_mesh%cell_vemid(cadj,ediagonal(ee)) .NE. 0)) then 
                ot_mesh%cell_vemid(ctgt,ee) = ot_mesh%cell_vemid(cadj,ediagonal(ee))
            end if 
        end if 
    end if 
end do 

!Build required face midpoint vertices
do aa=1,6 
    if (ot_mesh%cell_vfmid(ctgt,aa) == 0) then 
        ot_mesh%cell_vfmid(ctgt,aa) = vins 
        ot_mesh%vtx(vins,1) = 0.25d0*sum(ot_mesh%vtx(ot_mesh%cell_vcnr(ctgt,faces(aa,:)),1))
        ot_mesh%vtx(vins,2) = 0.25d0*sum(ot_mesh%vtx(ot_mesh%cell_vcnr(ctgt,faces(aa,:)),2))
        ot_mesh%vtx(vins,3) = 0.25d0*sum(ot_mesh%vtx(ot_mesh%cell_vcnr(ctgt,faces(aa,:)),3))
        vins = vins + 1
        if (vins .GE. cm3dopt%Ncell_max) then
            cm3dfailure = 1
            write(*,'(A)') '** maximum vertex count exceeded'
            return 
        end if
    end if 
end do 

!Build required edge midpoint vertices
do aa=1,12
    if (ot_mesh%cell_vemid(ctgt,aa) == 0) then 
        ot_mesh%cell_vemid(ctgt,aa) = vins 
        ot_mesh%vtx(vins,1) = 0.5d0*sum(ot_mesh%vtx(ot_mesh%cell_vcnr(ctgt,edges(aa,:)),1))
        ot_mesh%vtx(vins,2) = 0.5d0*sum(ot_mesh%vtx(ot_mesh%cell_vcnr(ctgt,edges(aa,:)),2))
        ot_mesh%vtx(vins,3) = 0.5d0*sum(ot_mesh%vtx(ot_mesh%cell_vcnr(ctgt,edges(aa,:)),3)) 
        vins = vins + 1
        if (vins .GE. cm3dopt%Ncell_max) then
            cm3dfailure = 1
            write(*,'(A)') '** maximum vertex count exceeded'
            return 
        end if
    end if
end do 

!Push face midpoint vertices
do aa=1,6 
    cadj = ot_mesh%cell_adjacent(ctgt,aa)
    if (ot_mesh%cell_vfmid(ctgt,aa) == 0) then 
        if (cadj .GT. 0) then 
            if (ot_mesh%cell_level(cadj) == ot_mesh%cell_level(ctgt)) then 
                if ((ot_mesh%cell_vfmid(cadj,fopposite(aa)) == 0) .AND. (ot_mesh%cell_vfmid(ctgt,aa) .NE. 0)) then 
                    ot_mesh%cell_vfmid(cadj,fopposite(aa)) = ot_mesh%cell_vfmid(ctgt,aa)
                end if 
            end if 
        end if
    end if 
end do 

!Push edge midpoint vertices
do ff=1,6 !though faces
    cadj = ot_mesh%cell_adjacent(ctgt,ff)
    if (cadj .GT. 0) then 
        if (ot_mesh%cell_level(cadj) == ot_mesh%cell_level(ctgt)) then 
            do aa=1,4
                if ((ot_mesh%cell_vemid(ctgt,face2edge(ff,aa)) .NE. 0) .AND. &
                (ot_mesh%cell_vemid(cadj,face2opedge(ff,aa)) == 0)) then 
                    ot_mesh%cell_vemid(cadj,face2opedge(ff,aa)) = ot_mesh%cell_vemid(ctgt,face2edge(ff,aa))
                end if 
            end do 
        end if 
    end if 
end do 
do ee=1,12 !though edges
    cadj = ot_mesh%cell_diagonal(ctgt,ee)
    if (cadj .GT. 0) then 
        if (ot_mesh%cell_level(cadj) == ot_mesh%cell_level(ctgt)) then 
            if ((ot_mesh%cell_vemid(ctgt,ee) .NE. 0) .AND. (ot_mesh%cell_vemid(cadj,ediagonal(ee)) == 0)) then 
                ot_mesh%cell_vemid(cadj,ediagonal(ee)) = ot_mesh%cell_vemid(ctgt,ee)
            end if 
        end if 
    end if 
end do 

!Construct new child cells ----
!child 1
ot_mesh%cell_vcnr(cins,1) = ot_mesh%cell_vcnr(ctgt,1)
ot_mesh%cell_vcnr(cins,2) = ot_mesh%cell_vemid(ctgt,1)
ot_mesh%cell_vcnr(cins,3) = ot_mesh%cell_vfmid(ctgt,5)
ot_mesh%cell_vcnr(cins,4) = ot_mesh%cell_vemid(ctgt,4)
ot_mesh%cell_vcnr(cins,5) = ot_mesh%cell_vemid(ctgt,5)
ot_mesh%cell_vcnr(cins,6) = ot_mesh%cell_vfmid(ctgt,1)
ot_mesh%cell_vcnr(cins,7) = vtxCC
ot_mesh%cell_vcnr(cins,8) = ot_mesh%cell_vfmid(ctgt,4)
ot_mesh%cell_child(ctgt,1) = cins
ot_mesh%cell_level(cins) = ot_mesh%cell_level(ctgt) + 1
ot_mesh%cell_parent(cins) = ctgt 
cins = cins + 1
if (cins .GE. cm3dopt%Ncell_max) then
    cm3dfailure = 1
    write(*,'(A)') 'maximum cell count exceeded'
    return 
end if

!child 2
ot_mesh%cell_vcnr(cins,1) = ot_mesh%cell_vemid(ctgt,1)
ot_mesh%cell_vcnr(cins,2) = ot_mesh%cell_vcnr(ctgt,2)
ot_mesh%cell_vcnr(cins,3) = ot_mesh%cell_vemid(ctgt,2)
ot_mesh%cell_vcnr(cins,4) = ot_mesh%cell_vfmid(ctgt,5)
ot_mesh%cell_vcnr(cins,5) = ot_mesh%cell_vfmid(ctgt,1)
ot_mesh%cell_vcnr(cins,6) = ot_mesh%cell_vemid(ctgt,6)
ot_mesh%cell_vcnr(cins,7) = ot_mesh%cell_vfmid(ctgt,3)
ot_mesh%cell_vcnr(cins,8) = vtxCC
ot_mesh%cell_child(ctgt,2) = cins
ot_mesh%cell_level(cins) = ot_mesh%cell_level(ctgt) + 1
ot_mesh%cell_parent(cins) = ctgt 
cins = cins + 1
if (cins .GE. cm3dopt%Ncell_max) then
    cm3dfailure = 1
    write(*,'(A)') 'maximum cell count exceeded'
    return 
end if

!child 3
ot_mesh%cell_vcnr(cins,1) = ot_mesh%cell_vfmid(ctgt,5)
ot_mesh%cell_vcnr(cins,2) = ot_mesh%cell_vemid(ctgt,2)
ot_mesh%cell_vcnr(cins,3) = ot_mesh%cell_vcnr(ctgt,3)
ot_mesh%cell_vcnr(cins,4) = ot_mesh%cell_vemid(ctgt,3)
ot_mesh%cell_vcnr(cins,5) = vtxCC
ot_mesh%cell_vcnr(cins,6) = ot_mesh%cell_vfmid(ctgt,3)
ot_mesh%cell_vcnr(cins,7) = ot_mesh%cell_vemid(ctgt,7)
ot_mesh%cell_vcnr(cins,8) = ot_mesh%cell_vfmid(ctgt,2)
ot_mesh%cell_child(ctgt,3) = cins
ot_mesh%cell_level(cins) = ot_mesh%cell_level(ctgt) + 1
ot_mesh%cell_parent(cins) = ctgt 
cins = cins + 1
if (cins .GE. cm3dopt%Ncell_max) then
    cm3dfailure = 1
    write(*,'(A)') 'maximum cell count exceeded'
    return 
end if

!child 4
ot_mesh%cell_vcnr(cins,1) = ot_mesh%cell_vemid(ctgt,4)
ot_mesh%cell_vcnr(cins,2) = ot_mesh%cell_vfmid(ctgt,5)
ot_mesh%cell_vcnr(cins,3) = ot_mesh%cell_vemid(ctgt,3)
ot_mesh%cell_vcnr(cins,4) = ot_mesh%cell_vcnr(ctgt,4)
ot_mesh%cell_vcnr(cins,5) = ot_mesh%cell_vfmid(ctgt,4)
ot_mesh%cell_vcnr(cins,6) = vtxCC
ot_mesh%cell_vcnr(cins,7) = ot_mesh%cell_vfmid(ctgt,2)
ot_mesh%cell_vcnr(cins,8) = ot_mesh%cell_vemid(ctgt,8)
ot_mesh%cell_child(ctgt,4) = cins
ot_mesh%cell_level(cins) = ot_mesh%cell_level(ctgt) + 1
ot_mesh%cell_parent(cins) = ctgt 
cins = cins + 1
if (cins .GE. cm3dopt%Ncell_max) then
    cm3dfailure = 1
    write(*,'(A)') 'maximum cell count exceeded'
    return 
end if

!child 5
ot_mesh%cell_vcnr(cins,1) = ot_mesh%cell_vemid(ctgt,5)
ot_mesh%cell_vcnr(cins,2) = ot_mesh%cell_vfmid(ctgt,1)
ot_mesh%cell_vcnr(cins,3) = vtxCC
ot_mesh%cell_vcnr(cins,4) = ot_mesh%cell_vfmid(ctgt,4)
ot_mesh%cell_vcnr(cins,5) = ot_mesh%cell_vcnr(ctgt,5)
ot_mesh%cell_vcnr(cins,6) = ot_mesh%cell_vemid(ctgt,9)
ot_mesh%cell_vcnr(cins,7) = ot_mesh%cell_vfmid(ctgt,6)
ot_mesh%cell_vcnr(cins,8) = ot_mesh%cell_vemid(ctgt,12)
ot_mesh%cell_child(ctgt,5) = cins
ot_mesh%cell_level(cins) = ot_mesh%cell_level(ctgt) + 1
ot_mesh%cell_parent(cins) = ctgt 
cins = cins + 1
if (cins .GE. cm3dopt%Ncell_max) then
    cm3dfailure = 1
    write(*,'(A)') 'maximum cell count exceeded'
    return 
end if

!child 6
ot_mesh%cell_vcnr(cins,1) = ot_mesh%cell_vfmid(ctgt,1)
ot_mesh%cell_vcnr(cins,2) = ot_mesh%cell_vemid(ctgt,6)
ot_mesh%cell_vcnr(cins,3) = ot_mesh%cell_vfmid(ctgt,3)
ot_mesh%cell_vcnr(cins,4) = vtxCC
ot_mesh%cell_vcnr(cins,5) = ot_mesh%cell_vemid(ctgt,9)
ot_mesh%cell_vcnr(cins,6) = ot_mesh%cell_vcnr(ctgt,6)
ot_mesh%cell_vcnr(cins,7) = ot_mesh%cell_vemid(ctgt,10)
ot_mesh%cell_vcnr(cins,8) = ot_mesh%cell_vfmid(ctgt,6)
ot_mesh%cell_child(ctgt,6) = cins
ot_mesh%cell_level(cins) = ot_mesh%cell_level(ctgt) + 1
ot_mesh%cell_parent(cins) = ctgt 
cins = cins + 1
if (cins .GE. cm3dopt%Ncell_max) then
    cm3dfailure = 1
    write(*,'(A)') 'maximum cell count exceeded'
    return 
end if

!child 7
ot_mesh%cell_vcnr(cins,1) = vtxCC
ot_mesh%cell_vcnr(cins,2) = ot_mesh%cell_vfmid(ctgt,3)
ot_mesh%cell_vcnr(cins,3) = ot_mesh%cell_vemid(ctgt,7)
ot_mesh%cell_vcnr(cins,4) = ot_mesh%cell_vfmid(ctgt,2)
ot_mesh%cell_vcnr(cins,5) = ot_mesh%cell_vfmid(ctgt,6)
ot_mesh%cell_vcnr(cins,6) = ot_mesh%cell_vemid(ctgt,10)
ot_mesh%cell_vcnr(cins,7) = ot_mesh%cell_vcnr(ctgt,7)
ot_mesh%cell_vcnr(cins,8) = ot_mesh%cell_vemid(ctgt,11)
ot_mesh%cell_child(ctgt,7) = cins
ot_mesh%cell_level(cins) = ot_mesh%cell_level(ctgt) + 1
ot_mesh%cell_parent(cins) = ctgt 
cins = cins + 1
if (cins .GE. cm3dopt%Ncell_max) then
    cm3dfailure = 1
    write(*,'(A)') 'maximum cell count exceeded'
    return 
end if

!child 8
ot_mesh%cell_vcnr(cins,1) = ot_mesh%cell_vfmid(ctgt,4)
ot_mesh%cell_vcnr(cins,2) = vtxCC
ot_mesh%cell_vcnr(cins,3) = ot_mesh%cell_vfmid(ctgt,2)
ot_mesh%cell_vcnr(cins,4) = ot_mesh%cell_vemid(ctgt,8)
ot_mesh%cell_vcnr(cins,5) = ot_mesh%cell_vemid(ctgt,12)
ot_mesh%cell_vcnr(cins,6) = ot_mesh%cell_vfmid(ctgt,6)
ot_mesh%cell_vcnr(cins,7) = ot_mesh%cell_vemid(ctgt,11)
ot_mesh%cell_vcnr(cins,8) = ot_mesh%cell_vcnr(ctgt,8)
ot_mesh%cell_child(ctgt,8) = cins
ot_mesh%cell_level(cins) = ot_mesh%cell_level(ctgt) + 1
ot_mesh%cell_parent(cins) = ctgt 
cins = cins + 1
if (cins .GE. cm3dopt%Ncell_max) then
    cm3dfailure = 1
    write(*,'(A)') 'maximum cell count exceeded'
    return 
end if

!Set and update cell adjacency in this region 
call set_local_adjacency(ot_mesh,ccadj,fopposite,ctgt) !through faces
call set_local_adjacency_diag(ot_mesh,ccadj_diag,ccadj_diag_dir,ediagonal,ctgt) !through edges

!Flood any edge or face midpoint vertices through the new child cells
do ccf=1,8

    !Child index
    child_idx = ot_mesh%cell_child(ctgt,ccf)

    !Pull adjacent face midpoint vertices 
    do aa=1,6 
        cadj = ot_mesh%cell_adjacent(child_idx,aa)
        if (ot_mesh%cell_vfmid(child_idx,aa) == 0) then 
            if (cadj .GT. 0) then 
                if (ot_mesh%cell_level(cadj) == ot_mesh%cell_level(child_idx)) then 
                    ot_mesh%cell_vfmid(child_idx,aa) = ot_mesh%cell_vfmid(cadj,fopposite(aa))
                end if 
            end if
        end if 
    end do 

    !Pull adjacent edge midpoint vertices 
    do ff=1,6
        cadj = ot_mesh%cell_adjacent(child_idx,ff)
        if (cadj .GT. 0) then 
            if (ot_mesh%cell_level(cadj) == ot_mesh%cell_level(child_idx)) then 
                do aa=1,4
                    if ((ot_mesh%cell_vemid(child_idx,face2edge(ff,aa)) == 0) .AND. &
                    (ot_mesh%cell_vemid(cadj,face2opedge(ff,aa)) .NE. 0)) then 
                        ot_mesh%cell_vemid(child_idx,face2edge(ff,aa)) = ot_mesh%cell_vemid(cadj,face2opedge(ff,aa))
                    end if 
                end do 
            end if 
        end if 
    end do 

    !Push face midpoint vertices
    do aa=1,6 
        cadj = ot_mesh%cell_adjacent(child_idx,aa)
        if (ot_mesh%cell_vfmid(child_idx,aa) == 0) then 
            if (cadj .GT. 0) then 
                if (ot_mesh%cell_level(cadj) == ot_mesh%cell_level(child_idx)) then 
                    if (ot_mesh%cell_vfmid(cadj,fopposite(aa)) == 0) then 
                        ot_mesh%cell_vfmid(cadj,fopposite(aa)) = ot_mesh%cell_vfmid(child_idx,aa)
                    end if 
                end if 
            end if
        end if 
    end do 

    !Push edge midpoint vertices
    do ff=1,6
        cadj = ot_mesh%cell_adjacent(child_idx,ff)
        if (cadj .GT. 0) then 
            if (ot_mesh%cell_level(cadj) == ot_mesh%cell_level(child_idx)) then 
                do aa=1,4
                    if ((ot_mesh%cell_vemid(child_idx,face2edge(ff,aa)) .NE. 0) .AND. &
                    (ot_mesh%cell_vemid(cadj,face2opedge(ff,aa)) == 0)) then 
                        ot_mesh%cell_vemid(cadj,face2opedge(ff,aa)) = ot_mesh%cell_vemid(child_idx,face2edge(ff,aa))
                    end if 
                end do 
            end if 
        end if 
    end do 
end do 
return 
end subroutine octree_cell_refine




!Cell child adjacency subroutine (edges) ===========================
subroutine set_local_adjacency_diag(ot_mesh,ccadj_diag,ccadj_diag_dir,ediagonal,cidx)
implicit none

!Variables - Import
integer(in) :: cidx
integer(in) :: ediagonal(12),ccadj_diag(8,12),ccadj_diag_dir(8,12)
type(octree_data) :: ot_mesh

!Variables - Local
integer(in) :: cc,aa,child_idx,cadj,cadjidx,etgt_op 

!Set child diagonal adjacency for each child cell
do cc=1,8
    child_idx = ot_mesh%cell_child(cidx,cc)
    do aa=1,12
        if (ccadj_diag_dir(cc,aa) == 0) then !Parent internal cell
            ot_mesh%cell_diagonal(child_idx,aa) = ot_mesh%cell_child(cidx,ccadj_diag(cc,aa))
        else
            if (ccadj_diag_dir(cc,aa) .GT. 0) then !Cell through parent face
                cadj = ot_mesh%cell_adjacent(cidx,ccadj_diag_dir(cc,aa)) 
            else
                cadj = ot_mesh%cell_diagonal(cidx,aa) 
            end if 
            if (cadj .GT. 0) then !Cell
                if (ot_mesh%cell_child(cadj,1) == 0) then !Adjacent cell has no children
                    ot_mesh%cell_diagonal(child_idx,aa) = cadj
                else

                    !Find index of adjacent child cell adjacent to this new child cell
                    cadjidx = ot_mesh%cell_child(cadj,ccadj_diag(cc,aa))

                    !Set adjacency of this parents child to the adjacent cells child
                    ot_mesh%cell_diagonal(child_idx,aa) = cadjidx

                    !Diagonaly opposite edge to the targeted child edge
                    etgt_op = ediagonal(aa)

                    !Set adjacency of adjacent child cadjidx to this parents child on the opposite edge
                    ot_mesh%cell_diagonal(cadjidx,etgt_op) = child_idx

                    !Cascade the adjacency of all children of this adjacent cell cadjidx in the direction etgt_op to child_idx if they are tagged with adjacency cidx on edge etgt_op currently
                    call cascade_adjacency_diag(ot_mesh,cadjidx,etgt_op,child_idx,cidx)
                end if 
            else !External boundary
                ot_mesh%cell_diagonal(child_idx,aa) = -2 
            end if 
        end if 
    end do 
end do 
end subroutine set_local_adjacency_diag




!Cell child adjacency cascade subroutine (edges) ===========================
recursive subroutine cascade_adjacency_diag(ot_mesh,cell2up,etgt,newadj,oldadj)
implicit none 

!Variables - Import
integer(in) :: cell2up,etgt,newadj,oldadj
type(octree_data), intent(inout) :: ot_mesh

!Variables - Local
integer(in) :: ch,cell2upN 

!Update current cell cell2up
if (ot_mesh%cell_diagonal(cell2up,etgt) == oldadj) then 
    ot_mesh%cell_diagonal(cell2up,etgt) = newadj
end if 

!Check for and update children of cell cell2up
do ch=1,8
    if (ot_mesh%cell_child(cell2up,ch) .NE. 0) then 
        cell2upN = ot_mesh%cell_child(cell2up,ch)
        call cascade_adjacency_diag(ot_mesh,cell2upN,etgt,newadj,oldadj)
    end if
end do 
return 
end subroutine cascade_adjacency_diag




!Cell child adjacency subroutine (faces) ===========================
subroutine set_local_adjacency(ot_mesh,ccadj,fopposite,cidx)
implicit none

!Variables - Import
integer(in) :: cidx
integer(in) :: fopposite(6),ccadj(8,6)
type(octree_data) :: ot_mesh

!Variables - Local
integer(in) :: cc,aa,child_idx,cadj,cadjidx,ftgt_op 

!Set child adjacency for each child cell
do cc=1,8
    child_idx = ot_mesh%cell_child(cidx,cc)
    do aa=1,6
        if (ccadj(cc,aa) .GT. 0) then !Parent internal cell
            ot_mesh%cell_adjacent(child_idx,aa) = ot_mesh%cell_child(cidx,ccadj(cc,aa))
        else !Parent external cell
            cadj = ot_mesh%cell_adjacent(cidx,aa) 
            if (cadj .GT. 0) then !Cell
                if (ot_mesh%cell_child(cadj,1) == 0) then !Adjacent cell has no children
                    ot_mesh%cell_adjacent(child_idx,aa) = cadj
                else !Adjacent cell has children

                    !Find index of adjacent child cell adjacent to this new child cell
                    cadjidx = ot_mesh%cell_child(cadj,abs(ccadj(cc,aa)))

                    !Set adjacency of this parents child to the adjacent cells child
                    ot_mesh%cell_adjacent(child_idx,aa) = cadjidx

                    !Opposite face to the targeted child face
                    ftgt_op = fopposite(aa)

                    !Set adjacency of adjacent child cadjidx to this parents child on the opposite face
                    ot_mesh%cell_adjacent(cadjidx,ftgt_op) = child_idx

                    !Cascade the adjacency of all children of this adjacent cell cadjidx in the direction ftgt to child_idx if they are tagged with adjacency cidx on face ftgt currently
                    call cascade_adjacency(ot_mesh,cadjidx,ftgt_op,child_idx,cidx)
                end if 
            else !External boundary
                ot_mesh%cell_adjacent(child_idx,aa) = cadj
            end if 
        end if 
    end do 
end do 
end subroutine set_local_adjacency




!Cell child adjacency cascade subroutine (faces) ===========================
recursive subroutine cascade_adjacency(ot_mesh,cell2up,ftgt,newadj,oldadj)
implicit none 

!Variables - Import
integer(in) :: cell2up,ftgt,newadj,oldadj
type(octree_data), intent(inout) :: ot_mesh

!Variables - Local
integer(in) :: ch,cell2upN 

!Update current cell cell2up
if (ot_mesh%cell_adjacent(cell2up,ftgt) == oldadj) then 
    ot_mesh%cell_adjacent(cell2up,ftgt) = newadj
end if 

!Check for and update children of cell cell2up
do ch=1,8
    if (ot_mesh%cell_child(cell2up,ch) .NE. 0) then 
        cell2upN = ot_mesh%cell_child(cell2up,ch)
        call cascade_adjacency(ot_mesh,cell2upN,ftgt,newadj,oldadj)
    end if
end do 
return 
end subroutine cascade_adjacency




!Face - cell overlap testing function ===========================
function face_cell_ovlp_bool(ot_mesh,surface_mesh,adtree,node_select,nselected,cidx) result(olstat)
implicit none 

!Variables - Import
integer(in) :: olstat,nselected,cidx
integer(in), dimension(:) :: node_select
type(tree_data) :: adtree
type(surface_data) :: surface_mesh
type(octree_data) :: ot_mesh

!Variables - Local
integer(in) :: nn,kk,exist 
real(dp) :: vf1(3),vf2(3),vf3(3),v1(3),v2(3),v3(3),v4(3),v5(3),v6(3),v7(3),v8(3)

!Test
exist = 0 
olstat = 0 
v1(:) = ot_mesh%vtx(ot_mesh%cell_vcnr(cidx,1),:)
v2(:) = ot_mesh%vtx(ot_mesh%cell_vcnr(cidx,2),:)
v3(:) = ot_mesh%vtx(ot_mesh%cell_vcnr(cidx,3),:)
v4(:) = ot_mesh%vtx(ot_mesh%cell_vcnr(cidx,4),:)
v5(:) = ot_mesh%vtx(ot_mesh%cell_vcnr(cidx,5),:)
v6(:) = ot_mesh%vtx(ot_mesh%cell_vcnr(cidx,6),:)
v7(:) = ot_mesh%vtx(ot_mesh%cell_vcnr(cidx,7),:)
v8(:) = ot_mesh%vtx(ot_mesh%cell_vcnr(cidx,8),:)
do nn=1,nselected
    do kk=1,adtree%tree(node_select(nn))%nentry

        !Verticies of this face
        vf1(:) = surface_mesh%vertices(surface_mesh%connectivity(adtree%tree(node_select(nn))%entry(kk),1),:)
        vf2(:) = surface_mesh%vertices(surface_mesh%connectivity(adtree%tree(node_select(nn))%entry(kk),2),:)
        vf3(:) = surface_mesh%vertices(surface_mesh%connectivity(adtree%tree(node_select(nn))%entry(kk),3),:)

        !Check for overlap of surface triangle face with cell
        exist = tri_cube_intersect_bool(v1,v2,v3,v4,v5,v6,v7,v8,vf1,vf2,vf3)

        !Exit if overlap found
        if (exist == 1) then 
            olstat = 1
            exit 
        end if
    end do
    if (exist == 1) then
        exit 
    end if
end do 
return 
end function face_cell_ovlp_bool




!Face - cell edge intersect testing function ===========================
function face_celledge_intersect_bool(ot_mesh,surface_mesh,adtree,node_select,nselected,cidx) result(olstat)
implicit none 

!Variables - Import
integer(in) :: olstat,nselected,cidx
integer(in), dimension(:) :: node_select
type(tree_data) :: adtree
type(surface_data) :: surface_mesh
type(octree_data) :: ot_mesh

!Variables - Local
integer(in) :: nn,kk,exist 
real(dp) :: vf1(3),vf2(3),vf3(3),v1(3),v2(3),v3(3),v4(3),v5(3),v6(3),v7(3),v8(3)

!Test
exist = 0 
olstat = 0 
v1(:) = ot_mesh%vtx(ot_mesh%cell_vcnr(cidx,1),:)
v2(:) = ot_mesh%vtx(ot_mesh%cell_vcnr(cidx,2),:)
v3(:) = ot_mesh%vtx(ot_mesh%cell_vcnr(cidx,3),:)
v4(:) = ot_mesh%vtx(ot_mesh%cell_vcnr(cidx,4),:)
v5(:) = ot_mesh%vtx(ot_mesh%cell_vcnr(cidx,5),:)
v6(:) = ot_mesh%vtx(ot_mesh%cell_vcnr(cidx,6),:)
v7(:) = ot_mesh%vtx(ot_mesh%cell_vcnr(cidx,7),:)
v8(:) = ot_mesh%vtx(ot_mesh%cell_vcnr(cidx,8),:)
do nn=1,nselected
    do kk=1,adtree%tree(node_select(nn))%nentry

        !Verticies of this face
        vf1(:) = surface_mesh%vertices(surface_mesh%connectivity(adtree%tree(node_select(nn))%entry(kk),1),:)
        vf2(:) = surface_mesh%vertices(surface_mesh%connectivity(adtree%tree(node_select(nn))%entry(kk),2),:)
        vf3(:) = surface_mesh%vertices(surface_mesh%connectivity(adtree%tree(node_select(nn))%entry(kk),3),:)

        !Check for overlap of surface triangle face with cell
        exist = tri_cubeedge_intersect_bool(v1,v2,v3,v4,v5,v6,v7,v8,vf1,vf2,vf3)

        !Exit if overlap found
        if (exist == 1) then 
            olstat = 1
            exit 
        end if
    end do
    if (exist == 1) then
        exit 
    end if
end do 
return 
end function face_celledge_intersect_bool




!Segment - cell overlap with surface radius of curvature check testing function ===========================
function face_cell_ovlp_wrcurv_bool(ot_mesh,surface_mesh,adtree,node_select,nselected,cidx) result(olstat)
implicit none 

!Variables - Import
integer(in) :: olstat,nselected,cidx
integer(in), dimension(:) :: node_select
type(tree_data) :: adtree
type(surface_data) :: surface_mesh
type(octree_data) :: ot_mesh

!Variables - Local
integer(in) :: nn,kk,vv,exist,curv_valid 
real(dp) :: sRcurv,cdx,cdy,cdz,reflen
real(dp) :: vf1(3),vf2(3),vf3(3),v1(3),v2(3),v3(3),v4(3),v5(3),v6(3),v7(3),v8(3)

!Test
exist = 0 
olstat = 0 
v1(:) = ot_mesh%vtx(ot_mesh%cell_vcnr(cidx,1),:)
v2(:) = ot_mesh%vtx(ot_mesh%cell_vcnr(cidx,2),:)
v3(:) = ot_mesh%vtx(ot_mesh%cell_vcnr(cidx,3),:)
v4(:) = ot_mesh%vtx(ot_mesh%cell_vcnr(cidx,4),:)
v5(:) = ot_mesh%vtx(ot_mesh%cell_vcnr(cidx,5),:)
v6(:) = ot_mesh%vtx(ot_mesh%cell_vcnr(cidx,6),:)
v7(:) = ot_mesh%vtx(ot_mesh%cell_vcnr(cidx,7),:)
v8(:) = ot_mesh%vtx(ot_mesh%cell_vcnr(cidx,8),:)
cdx = abs(max(v1(1),v2(1),v3(1),v4(1),v5(1),v6(1),v7(1),v8(1)) - min(v1(1),v2(1),v3(1),v4(1),v5(1),v6(1),v7(1),v8(1)))
cdy = abs(max(v1(2),v2(2),v3(2),v4(2),v5(2),v6(2),v7(2),v8(2)) - min(v1(2),v2(2),v3(2),v4(2),v5(2),v6(2),v7(2),v8(2)))
cdz = abs(max(v1(3),v2(3),v3(3),v4(3),v5(3),v6(3),v7(3),v8(3)) - min(v1(3),v2(3),v3(3),v4(3),v5(3),v6(3),v7(3),v8(3)))
reflen = (cdx*cdy*cdz)**(1.0d0/3.0d0)
do nn=1,nselected
    do kk=1,adtree%tree(node_select(nn))%nentry

        !Verticies of this face
        vf1(:) = surface_mesh%vertices(surface_mesh%connectivity(adtree%tree(node_select(nn))%entry(kk),1),:)
        vf2(:) = surface_mesh%vertices(surface_mesh%connectivity(adtree%tree(node_select(nn))%entry(kk),2),:)
        vf3(:) = surface_mesh%vertices(surface_mesh%connectivity(adtree%tree(node_select(nn))%entry(kk),3),:)

        !Check for overlap of surface triangle face with cell
        exist = tri_cube_intersect_bool(v1,v2,v3,v4,v5,v6,v7,v8,vf1,vf2,vf3)

        !Check if surface radius of curvature here is smaller than the cell size curv_valid
        curv_valid = 0 
        sRcurv = surface_mesh%face_rcurv(adtree%tree(node_select(nn))%entry(kk))
        if (sRcurv .LE. reflen) then 
            curv_valid = 1 
        end if 
        sRcurv = surface_mesh%face_maxcurv(adtree%tree(node_select(nn))%entry(kk))
        if (sRcurv .LE. reflen) then 
            curv_valid = 1 
        end if
        do vv=1,3
            sRcurv = surface_mesh%vtx_rcurv(surface_mesh%connectivity(adtree%tree(node_select(nn))%entry(kk),vv))
            if (sRcurv .LE. reflen) then 
                curv_valid = 1 
            end if 
        end do 

        !Exit if overlap found
        if ((exist == 1) .AND. (curv_valid == 1)) then 
            olstat = 1
            exit 
        end if
    end do
    if ((exist == 1) .AND. (curv_valid == 1)) then 
        exit 
    end if
end do 
return 
end function face_cell_ovlp_wrcurv_bool




!Adjacency difference calculation function ===========================
function max_adjacent_dR(ot_mesh,cidx) result(dR)
implicit none 

!Variables - Import
integer(in) :: dR,cidx
type(octree_data) :: ot_mesh

!Variables - Local
integer(in) :: kk,lc 
integer(in) :: deltaR(6)

!Check adjacent cells 
lc = ot_mesh%cell_level(cidx)
do kk=1,6
    if (ot_mesh%cell_adjacent(cidx,kk) .GT. 0) then
        deltaR(kk) = ot_mesh%cell_level(ot_mesh%cell_adjacent(cidx,kk)) - lc
    else
        deltaR(kk) = 0 
    end if
end do 
dR = maxval(abs(deltaR(:)))
return 
end function max_adjacent_dR




!Identify adjacent child cells cascade subroutine ===========================
recursive subroutine find_adjacent_children(cell_keep,frontN,NfrontN,Nexternal,cell_tgt,cell2check,ftgt_op,ot_mesh,cell_nc)
implicit none 

!System data
type(octree_data), intent(in) :: ot_mesh

!Variables - Import
integer(in), intent(inout) :: NfrontN,Nexternal
integer(in) :: cell_tgt,cell2check,ftgt_op
integer(in) :: cell_nc(ot_mesh%cins-1)
integer(in), intent(inout) :: cell_keep(ot_mesh%cins-1),frontN(ot_mesh%cins-1)

!Variables - Local 
integer(in) :: ch,child_idx

!Check children of cell cell2check for those that are adjacent to cell_tgt that have no children 
do ch=1,8
    if (ot_mesh%cell_child(cell2check,ch) .NE. 0) then 
        child_idx = ot_mesh%cell_child(cell2check,ch)
        if (ot_mesh%cell_adjacent(child_idx,ftgt_op) == cell_tgt) then !If adjacent to current target cell 
            if (cell_nc(child_idx) == 1) then !If cell child_idx has no children 
                if (cell_keep(child_idx) == 0) then !Not yet visited -> hence add to new front and set as external
                    NfrontN = NfrontN + 1
                    Nexternal = Nexternal + 1
                    cell_keep(child_idx) = 2
                    frontN(NfrontN) = child_idx
                end if
            else !Search this cells children for valid cells if any -> recursive call 
                call find_adjacent_children(cell_keep,frontN,NfrontN,Nexternal,cell_tgt,child_idx,ftgt_op,ot_mesh,cell_nc)
            end if
        end if
    end if 
end do 
return 
end subroutine find_adjacent_children




!Get diagonal cell subroutine ===========================
subroutine get_diagonal_cell(cdiag,cbase,ebase,eopposite_diag,f1,f2,edges,ot_mesh)
implicit none 

!Variables - Import 
integer(in) :: cdiag,cbase,ebase,eopposite_diag,f1,f2
integer(in) :: edges(12,2)
type(octree_data) :: ot_mesh

!Variables - Local 
integer(in) :: aa,cadj1,cadj2,ctgt_search,cadj_test
integer(in) :: v1_eadj,v2_eadj,v1_ebase,v2_ebase,ematched

!Initialise 
cdiag = 0 

!Adjacent cells across each face
cadj1 = ot_mesh%cell_adjacent(cbase,f1)
cadj2 = ot_mesh%cell_adjacent(cbase,f2)

!Pick cell to check against
ctgt_search = 0 
if (cadj1 .GT. 0) then 
    if (ot_mesh%cell_level(cadj1) == ot_mesh%cell_level(cbase)) then 
        ctgt_search = cadj1
    end if 
end if 
if (ctgt_search == 0) then 
    if (cadj2 .GT. 0) then 
        if (ot_mesh%cell_level(cadj2) == ot_mesh%cell_level(cbase)) then 
            ctgt_search = cadj2
        end if 
    end if 
end if 

!If no valid cells -> return no diagonal cell 
if (ctgt_search == 0) then 
    cdiag = 0 
    return 
end if 

!Base edge vertices
v1_ebase = ot_mesh%cell_vcnr(cbase,edges(ebase,1))
v2_ebase = ot_mesh%cell_vcnr(cbase,edges(ebase,2))

!Search cell ctgt_search for an adjacent cell that is not cbase that shares the correct edge 
ematched = 0 
do aa=1,6
    cadj_test = ot_mesh%cell_adjacent(ctgt_search,aa)
    if ((cadj_test .NE. cbase) .AND. (cadj_test .GT. 0)) then 
        if (ot_mesh%cell_level(cadj_test) == ot_mesh%cell_level(cbase)) then 
            v1_eadj = ot_mesh%cell_vcnr(cadj_test,edges(eopposite_diag,1))
            v2_eadj = ot_mesh%cell_vcnr(cadj_test,edges(eopposite_diag,2))
            if ((v1_eadj == v1_ebase) .AND. (v2_eadj == v2_ebase)) then
                ematched = 1
            elseif ((v2_eadj == v1_ebase) .AND. (v1_eadj == v2_ebase)) then 
                ematched = 1
            end if 
            if (ematched == 1) then 
                cdiag = cadj_test
                exit 
            end if 
        end if 
    end if 
end do 
return 
end subroutine get_diagonal_cell




!Subedge vertex identification subroutine ===========================
recursive subroutine get_otmesh_edge_mid_vertices(subcedge_vtx,Nsubcedge_vtx,ctgt,etgt,ot_mesh,edges)
implicit none 

!Variables - Import 
integer(in), intent(inout) :: ctgt,etgt,Nsubcedge_vtx
integer(in) :: edges(12,2)
integer(in), dimension(:) :: subcedge_vtx
type(octree_data) :: ot_mesh

!Variables - Local 
integer(in) :: cc,vv,vtgt,exist,c_child,etgt_v1,etgt_v2,echld_v1,echld_v2  

!Check if target edge contains a midpoint vertex
if (ot_mesh%cell_vemid(ctgt,etgt) .NE. 0) then 

    !Target midpoint vertex
    vtgt = ot_mesh%cell_vemid(ctgt,etgt)

    !Check if vertex is new 
    exist = 0 
    do vv=1,Nsubcedge_vtx
        if (subcedge_vtx(vv) == vtgt) then 
            exist = 1 
            exit 
        end if 
    end do 

    !Add this edge midpoint vertex to the list if new
    if (exist == 0) then 
        Nsubcedge_vtx = Nsubcedge_vtx + 1
        subcedge_vtx(Nsubcedge_vtx) = vtgt
    end if 

    !If this cell has children then check the target edge in the children that contain either vertex of etgt
    if (ot_mesh%cell_child(ctgt,1) .NE. 0) then 

        !Vertices in the current cell on the target edge 
        etgt_v1 = ot_mesh%cell_vcnr(ctgt,edges(etgt,1))
        etgt_v2 = ot_mesh%cell_vcnr(ctgt,edges(etgt,2)) 
        
        !Check children
        do cc=1,8

            !Child cell 
            c_child = ot_mesh%cell_child(ctgt,cc)

            !Target edge ends in child cell
            echld_v1 = ot_mesh%cell_vcnr(c_child,edges(etgt,1))
            echld_v2 = ot_mesh%cell_vcnr(c_child,edges(etgt,2))

            !Check this cell and edge if either vertex is shared with the parent on the targeted edge
            if ((echld_v1 == etgt_v1) .OR. (echld_v2 == etgt_v2) .OR. (echld_v1 == etgt_v2) .OR. (echld_v2 == etgt_v1)) then 
                call get_otmesh_edge_mid_vertices(subcedge_vtx,Nsubcedge_vtx,c_child,etgt,ot_mesh,edges)
            end if 
        end do 
    end if 
end if 
return 
end subroutine get_otmesh_edge_mid_vertices


end module cell_mesh3d_octree_mod