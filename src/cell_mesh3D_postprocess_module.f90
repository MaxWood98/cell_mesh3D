!cell_mesh3d mesh postprocessing module
!Max Wood - mw16116@bristol.ac.uk
!Univeristy of Bristol - Department of Aerospace Engineering

!Version 2.2
!Updated 03-08-2023

!Module
module cellmesh3d_postprocess_mod
! use cellmesh3d_data_mod
! use cellmesh3d_geometry_mod
use cellmesh3d_mesh_build_mod
use ieee_arithmetic, only: ieee_value,IEEE_POSITIVE_INF
contains 


!Cell index remapping subroutine ===========================
subroutine remap_cell_indecies(volume_mesh)
implicit none 

!Variables - Import
type(vol_mesh_data) :: volume_mesh

!Variables - Local 
integer(in) :: cc,ff,cl,cr,NcellN
integer(in) :: cell_indexN(volume_mesh%ncell)
integer(in), dimension(:), allocatable :: cell_levelN

!Remap cells
NcellN = 0 
cell_indexN(:) = 0 
do ff=1,volume_mesh%nface
    cl = volume_mesh%faces(ff)%cleft
    cr = volume_mesh%faces(ff)%cright
    if (cl .GT. 0) then 
        if (cell_indexN(cl) == 0) then
            NcellN = NcellN + 1
            cell_indexN(cl) = NcellN
        end if 
    end if
    if (cr .GT. 0) then 
        if (cell_indexN(cr) == 0) then
            NcellN = NcellN + 1
            cell_indexN(cr) = NcellN
        end if 
    end if
end do 

!Reasign indecies on edges
do ff=1,volume_mesh%nface
    cl = volume_mesh%faces(ff)%cleft
    cr = volume_mesh%faces(ff)%cright
    if (cl .GT. 0) then 
        if (cell_indexN(cl) .GT. 0) then
            volume_mesh%faces(ff)%cleft = cell_indexN(cl)
        end if 
    end if
    if (cr .GT. 0) then 
        if (cell_indexN(cr) .GT. 0) then
            volume_mesh%faces(ff)%cright = cell_indexN(cr)
        end if 
    end if
end do 

!Reassign cell levels 
allocate(cell_levelN(NcellN))
do cc=1,volume_mesh%ncell
    if (cell_indexN(cc) .GT. 0) then 
        cell_levelN(cell_indexN(cc)) = volume_mesh%cell_level(cc)
    end if 
end do 
deallocate(volume_mesh%cell_level)
allocate(volume_mesh%cell_level(NcellN))
volume_mesh%cell_level(:) = cell_levelN(:)

!Assign new cell count 
volume_mesh%ncell = NcellN
return 
end subroutine remap_cell_indecies




!Remove zero area faces subroutine ===========================
subroutine remove_zero_area_faces(volume_mesh,cm3dopt)
implicit none 

!Variables - Import
type(cm3d_options) :: cm3dopt
type(vol_mesh_data) :: volume_mesh

!Variables - Local 
integer(in) :: ff,NFvalid,NFrem
integer(in) :: findex(volume_mesh%nface)
real(dp) :: Aface
real(dp) :: Nf(3)
type(vol_mesh_data) :: volume_mesh_temp

!Index valid faces
NFrem = 0 
NFvalid = 0 
findex(:) = 0 
do ff=1,volume_mesh%nface

    !Normal vector of this face
    Nf = newell_normal(volume_mesh%faces(ff)%nvtx,volume_mesh%faces(ff)%vertices,volume_mesh%vtx)

    !Area of this face
    Aface = norm2(Nf)

    !Tag to remove if zero area
    if (Aface == 0.0d0) then 
        NFrem = NFrem + 1
    else
        NFvalid = NFvalid + 1
        findex(ff) = NFvalid
    end if 
end do 

!Build new face list 
allocate(volume_mesh_temp%faces(NFvalid))
do ff=1,volume_mesh%nface
    if (findex(ff) .NE. 0) then 
        volume_mesh_temp%faces(findex(ff))%nvtx = volume_mesh%faces(ff)%nvtx
        allocate(volume_mesh_temp%faces(findex(ff))%vertices(volume_mesh%faces(ff)%nvtx))
        volume_mesh_temp%faces(findex(ff))%vertices(:) = volume_mesh%faces(ff)%vertices(:)
        volume_mesh_temp%faces(findex(ff))%cleft = volume_mesh%faces(ff)%cleft
        volume_mesh_temp%faces(findex(ff))%cright = volume_mesh%faces(ff)%cright
    end if 
end do 
deallocate(volume_mesh%faces)
allocate(volume_mesh%faces(NFvalid))
volume_mesh%nface = NFvalid
do ff=1,NFvalid
    volume_mesh%faces(ff)%nvtx = volume_mesh_temp%faces(ff)%nvtx
    allocate(volume_mesh%faces(ff)%vertices(volume_mesh%faces(ff)%nvtx))
    allocate(volume_mesh%faces(ff)%edges(volume_mesh%faces(ff)%nvtx))
    volume_mesh%faces(ff)%edges(:) = 0 
    volume_mesh%faces(ff)%vertices(:) = volume_mesh_temp%faces(ff)%vertices(:) 
    volume_mesh%faces(ff)%cleft = volume_mesh_temp%faces(ff)%cleft
    volume_mesh%faces(ff)%cright = volume_mesh_temp%faces(ff)%cright
end do 

!Display
if (NFrem .NE. 0) then 
    if (cm3dopt%dispt == 1) then
        write(*,'(A,I0,A)') '    {removed ',NFrem,' zero area faces}'
    end if
end if 
return 
end subroutine remove_zero_area_faces




!Correct bisected cells subroutine ===========================
subroutine correct_bisected_cells(volume_mesh,cm3dopt)
implicit none 

!Variables - Import
type(cm3d_options) :: cm3dopt
type(vol_mesh_data) :: volume_mesh

!Variables - Local 
integer(in) :: ff,ee,cc,ii,jj
integer(in) :: nedge,cl,cr,maxcface,maxeface,etgt,fbase,ftgt,nftag,nbisect,ncellN,fsurf,hassubfaces
integer(in) :: cellnface(volume_mesh%ncell),cell_tag(volume_mesh%nface)
integer(in) :: cell_levelN(2*volume_mesh%ncell)
integer(in), dimension(:), allocatable :: edgenface
integer(in), dimension(:,:), allocatable :: vmf_edges,cell2face,edge2face

!Build complete volume mesh edges 
call build_vmesh_edges(nedge,vmf_edges,volume_mesh,20_in)
allocate(edgenface(nedge))

!Build cell2face
cellnface(:) = 0 
do ff=1,volume_mesh%nface
    cl = volume_mesh%faces(ff)%cleft
    cr = volume_mesh%faces(ff)%cright
    if (cl .GT. 0) then 
        cellnface(cl) = cellnface(cl) + 1
    end if 
    if (cr .GT. 0) then 
        cellnface(cr) = cellnface(cr) + 1
    end if 
end do 
maxcface = maxval(cellnface)
allocate(cell2face(volume_mesh%ncell,maxcface))
cellnface(:) = 0 
cell2face(:,:) = 0 
do ff=1,volume_mesh%nface
    cl = volume_mesh%faces(ff)%cleft
    cr = volume_mesh%faces(ff)%cright
    if (cl .GT. 0) then 
        cellnface(cl) = cellnface(cl) + 1
        cell2face(cl,cellnface(cl)) = ff 
    end if 
    if (cr .GT. 0) then 
        cellnface(cr) = cellnface(cr) + 1
        cell2face(cr,cellnface(cr)) = ff 
    end if 
end do 

!Build edge2face 
edgenface(:) = 0 
do ff=1,volume_mesh%nface
    do ee=1,volume_mesh%faces(ff)%nvtx
        etgt = volume_mesh%faces(ff)%edges(ee)
        edgenface(etgt) = edgenface(etgt) + 1
    end do 
end do 
maxeface = maxval(edgenface)
allocate(edge2face(nedge,maxeface))
edgenface(:) = 0 
edge2face(:,:) = 0 
do ff=1,volume_mesh%nface
    do ee=1,volume_mesh%faces(ff)%nvtx
        etgt = volume_mesh%faces(ff)%edges(ee)
        edgenface(etgt) = edgenface(etgt) + 1
        edge2face(etgt,edgenface(etgt)) = ff 
    end do 
end do 

!Check for and correct bisected cells 
nbisect = 0 
cell_tag(:) = 0 
ncellN = volume_mesh%ncell
cell_levelN(1:volume_mesh%ncell) = volume_mesh%cell_level(:)
do cc=1,volume_mesh%ncell
    
    !If cell intersects the surface 
    fsurf = 0 
    do ii=1,cellnface(cc)
        ftgt = cell2face(cc,ii)
        if ((volume_mesh%faces(ftgt)%cleft == -1) .OR. (volume_mesh%faces(ftgt)%cright == -1)) then 
            fsurf = 1
            exit 
        end if 
    end do 
    if (fsurf == 1) then 

        !Flood all faces in this cell from the first 
        nftag = 1
        cell_tag(cell2face(cc,1)) = 1
        do ii=1,cellnface(cc) !Flood iteration
            do jj=1,cellnface(cc) !Loop all faces
                if (cell_tag(cell2face(cc,jj)) == 1) then !Flood tag to adjacent vertices in this cell along edges that border this cell 
                    
                    !Face to flood from 
                    fbase = cell2face(cc,jj)

                    !Flood across all edges of this cell 
                    do ee=1,volume_mesh%faces(fbase)%nvtx

                        !Edge
                        etgt = volume_mesh%faces(fbase)%edges(ee)

                        !Each face on this edge 
                        do ff=1,edgenface(etgt)
                            ftgt = edge2face(etgt,ff)
                            if (cell_tag(ftgt) == 0) then 
                                cl = volume_mesh%faces(ftgt)%cleft
                                cr = volume_mesh%faces(ftgt)%cright
                                if ((cl == cc) .OR. (cr == cc)) then 
                                    cell_tag(ftgt) = 1
                                    nftag = nftag + 1
                                end if  
                            end if
                        end do 
                    end do 
                end if
            end do 
        end do 

        !Flood from any tagged sub faces if any are present 
        hassubfaces = 0 
        do ii=1,cellnface(cc)
            ftgt = cell2face(cc,ii)
            if (volume_mesh%faces(ftgt)%ftype == 'int') then 
                hassubfaces = 1
                exit 
            end if 
        end do 
        if (hassubfaces == 1) then 

            !Find and tag any sub-faces in this cell if their parent face has been tagged  
            do ii=1,cellnface(cc)
                ftgt = cell2face(cc,ii)
                if (volume_mesh%faces(ftgt)%ftype == 'int') then 
                    cell_tag(ftgt) = cell_tag(volume_mesh%faces(ftgt)%fparent)
                    ! print *, 'intface'
                    ! print *, cell_tag(volume_mesh%faces(ftgt)%fparent)
                end if 
            end do 

            !Flood state though from the tagged sub faces
            do ii=1,cellnface(cc) !Flood iteration
                do jj=1,cellnface(cc) !Loop all faces
                    if (cell_tag(cell2face(cc,jj)) == 1) then !Flood tag to adjacent vertices in this cell along edges that border this cell 
                        
                        !Face to flood from 
                        fbase = cell2face(cc,jj)
    
                        !Flood across all edges of this cell 
                        do ee=1,volume_mesh%faces(fbase)%nvtx
    
                            !Edge
                            etgt = volume_mesh%faces(fbase)%edges(ee)
    
                            !Each face on this edge 
                            do ff=1,edgenface(etgt)
                                ftgt = edge2face(etgt,ff)
                                if (cell_tag(ftgt) == 0) then 
                                    cl = volume_mesh%faces(ftgt)%cleft
                                    cr = volume_mesh%faces(ftgt)%cright
                                    if ((cl == cc) .OR. (cr == cc)) then 
                                        cell_tag(ftgt) = 1
                                        nftag = nftag + 1
                                    end if  
                                end if
                            end do 
                        end do 
                    end if
                end do 
            end do 
        end if 

        !Check if all faces in the cell have been tagged - if not then cell is bisected
        if (nftag .NE. cellnface(cc)) then 

            !Increment counts 
            ncellN = ncellN + 1
            nbisect = nbisect + 1

            !Store level of this new cell 
            cell_levelN(ncellN) = cell_levelN(cc)

            !Tag adjacency for this cell of all untagged faces in this cell with the new cell index 
            do ff=1,cellnface(cc)
                ftgt = cell2face(cc,ff)
                if (cell_tag(ftgt) == 0) then 
                    if (volume_mesh%faces(ftgt)%cleft == cc) then 
                        volume_mesh%faces(ftgt)%cleft = ncellN
                    end if 
                    if (volume_mesh%faces(ftgt)%cright == cc) then 
                        volume_mesh%faces(ftgt)%cright = ncellN
                    end if 
                end if 
            end do 
            ! print *, 'cell bisected ',cc ,' -- ',nftag,' / ',cellnface(cc)
        end if

        !Reset cell tags
        cell_tag(cell2face(cc,1:cellnface(cc))) = 0
    end if 
end do 

!Assign new cell count and cell levels 
volume_mesh%ncell = ncellN
deallocate(volume_mesh%cell_level)
allocate(volume_mesh%cell_level(ncellN))
volume_mesh%cell_level(:) = cell_levelN(1:ncellN)

!Tag faces that map to the same cell on both sides to be removed
do ff=1,volume_mesh%nface
    if (volume_mesh%faces(ff)%cleft == volume_mesh%faces(ff)%cright) then 
        print *, '** double face : ', ff 
    end if
end do 

!Display
if (nbisect .GT. 0) then 
    if (cm3dopt%dispt == 1) then
        write(*,'(A,I0,A)') '    {identified and split ',nbisect,' bisected cells}'
    end if 
end if 
return 
end subroutine correct_bisected_cells




!Mesh sliver cell cleaning subroutine ===========================
subroutine clean_mesh_sliverC(volume_mesh,Cvol,cm3dopt,Nmerge,Nmerge_fail)
implicit none 

!Variables - Import
integer(in) :: Nmerge,Nmerge_fail
real(dp), dimension(:), allocatable :: Cvol
type(cm3d_options) :: cm3dopt
type(vol_mesh_data) :: volume_mesh
    
!Variables - Local 
integer(in) :: cc,ff,aa,aa2,aa3 
integer(in) :: Ncremove,Ninvalid,maxcface,cl,cr,ctgt,ftgt,cselect,fselect,merge_invalid,NcellN,NfaceN
integer(in) :: ctgta,ftgta,ctgta2,ftgta2,invalid_trigger
integer(in) :: cell_surfadj(volume_mesh%ncell),cell_remove(volume_mesh%ncell),cellnface(volume_mesh%ncell)
integer(in) :: face_remove(volume_mesh%nface),cell_map(volume_mesh%ncell),face_map(volume_mesh%nface)
integer(in), dimension(:,:), allocatable :: cell2face
real(dp) :: vref,vmax
type(vol_mesh_data) :: volume_mesh_temp

!Initialise 
Nmerge = 0
Ncremove = 0 
Ninvalid = 0 
Nmerge_fail = 0 

!Tag surface adjacent cells 
cell_surfadj(:) = 0 
do ff=1,volume_mesh%nface
    if ((volume_mesh%faces(ff)%cleft == -1)) then 
        cell_surfadj(volume_mesh%faces(ff)%cright) = 1
    end if 
    if ((volume_mesh%faces(ff)%cright == -1)) then 
        cell_surfadj(volume_mesh%faces(ff)%cleft) = 1
    end if 
end do 

!Evaluate cell volumes 
call get_cell_volumes(Cvol,volume_mesh)

!Tag cells to remove 
cell_remove(:) = 0 
do cc=1,volume_mesh%ncell 
    if (cell_surfadj(cc) == 1) then
        vref = (2.0d0*cm3dopt%far_field_bound/(2.0d0**(volume_mesh%cell_level(cc) - 1)))**3
        if (Cvol(cc) .LE. cm3dopt%CminVol*vref) then 
            Ncremove = Ncremove + 1
            cell_remove(cc) = 1
        end if
    end if 
end do 

!Remove cells if required
if (Ncremove .NE. 0) then 

    !Build cell2face
    cellnface(:) = 0 
    do ff=1,volume_mesh%nface
        cl = volume_mesh%faces(ff)%cleft
        cr = volume_mesh%faces(ff)%cright
        if (cl .GT. 0) then 
            cellnface(cl) = cellnface(cl) + 1
        end if 
        if (cr .GT. 0) then 
            cellnface(cr) = cellnface(cr) + 1
        end if 
    end do 
    maxcface = maxval(cellnface)
    allocate(cell2face(volume_mesh%ncell,maxcface))
    cellnface(:) = 0 
    cell2face(:,:) = 0 
    do ff=1,volume_mesh%nface
        cl = volume_mesh%faces(ff)%cleft
        cr = volume_mesh%faces(ff)%cright
        if (cl .GT. 0) then 
            cellnface(cl) = cellnface(cl) + 1
            cell2face(cl,cellnface(cl)) = ff 
        end if 
        if (cr .GT. 0) then 
            cellnface(cr) = cellnface(cr) + 1
            cell2face(cr,cellnface(cr)) = ff 
        end if 
    end do 

    !Build cell mapping for retained cells
    NcellN = 0 
    cell_map(:) = 0 
    do cc=1,volume_mesh%ncell 
        if (cell_remove(cc) .NE. 1) then 
            NcellN = NcellN + 1
            cell_map(cc) = NcellN
        end if
    end do

    !Tag faces to remove to eliminate sliver cells and map removed cells to their replacements
    face_remove(:) = 0 
    do cc=1,volume_mesh%ncell 
        if (cell_remove(cc) == 1) then 
            if (Cvol(cc) == 0.0d0) then !Force merge with any adjacent cell 

                !Find any valid adjacent cell to merge with 
                cselect = 0 
                fselect = 0 
                vmax = 0.0d0 
                invalid_trigger = 0 
                do aa=1,cellnface(cc)
                    ftgt = cell2face(cc,aa)
                    if ((volume_mesh%faces(ftgt)%cleft == cc)) then 
                        ctgt = volume_mesh%faces(ftgt)%cright
                    else
                        ctgt = volume_mesh%faces(ftgt)%cleft
                    end if 
                    if (ctgt .GT. 0) then 
                        if (Cvol(ctgt) .GT. vmax) then 
                            if (cell_remove(ctgt) == 0) then 
                                fselect = ftgt   
                                cselect = ctgt
                                exit 
                            end if
                        end if
                    end if
                end do 

                !Tag if valid merge identified
                if (cselect .GT. 0) then !Valid merge found -> Map cell for replacement and tag faces joining ii and ctgt for removal
                    Nmerge = Nmerge + 1
                    cell_map(cc) = cell_map(cselect) 
                    face_remove(fselect) = 1
                    do aa=1,cellnface(cc)
                        ftgt = cell2face(cc,aa)
                        if ((volume_mesh%faces(ftgt)%cleft == cc)) then 
                            ctgt = volume_mesh%faces(ftgt)%cright
                        else
                            ctgt = volume_mesh%faces(ftgt)%cleft
                        end if 
                        if (ctgt == cselect) then 
                            face_remove(ftgt) = 1
                        end if 
                    end do
                else !No valid merge case identified
                    cell_remove(cc) = 2
                    Nmerge_fail = Nmerge_fail + 1
                end if 
            else !Merge as normal to largest adjacent cell 

                !Find adjacent cell with largest volume to merge with 
                cselect = 0 
                fselect = 0 
                vmax = 0.0d0 
                invalid_trigger = 0 
                do aa=1,cellnface(cc)
                    ftgt = cell2face(cc,aa)
                    if ((volume_mesh%faces(ftgt)%cleft == cc)) then 
                        ctgt = volume_mesh%faces(ftgt)%cright
                    else
                        ctgt = volume_mesh%faces(ftgt)%cleft
                    end if 
                    if (ctgt .GT. 0) then 
                        if (Cvol(ctgt) .GT. vmax) then 
                            if (cell_remove(ctgt) == 0) then 
                                if (volume_mesh%cell_level(ctgt) == volume_mesh%cell_level(cc)) then 

                                    !Check if there is a cell adjacent to the target merge cell ctgt that also borders the original cell cc
                                    merge_invalid = 0 
                                    do aa2=1,cellnface(ctgt) !cells adjacent to the target merge cell 
                                        ftgta = cell2face(cc,aa)
                                        if ((volume_mesh%faces(ftgta)%cleft == ctgt)) then 
                                            ctgta = volume_mesh%faces(ftgta)%cright
                                        else
                                            ctgta = volume_mesh%faces(ftgta)%cleft
                                        end if 
                                        if ((ctgta .NE. cc) .AND. (ctgta .GT. 0)) then 
                                            do aa3=1,cellnface(cc) !cells adjacent to the base cell 
                                                ftgta2 = cell2face(cc,aa3)
                                                if ((volume_mesh%faces(ftgta2)%cleft == cc)) then 
                                                    ctgta2 = volume_mesh%faces(ftgta2)%cright
                                                else
                                                    ctgta2 = volume_mesh%faces(ftgta2)%cleft
                                                end if
                                                if ((ctgta2 .NE. cc) .AND. (ctgta2 .GT. 0)) then 
                                                    if (ctgta2 == ctgta) then 
                                                        merge_invalid = 1
                                                        invalid_trigger = 1
                                                        exit 
                                                    end if 
                                                end if 
                                            end do 
                                        end if 
                                    end do 

                                    !If valid merge then select cell ctgt
                                    if (merge_invalid == 0) then 
                                        vmax = Cvol(ctgt)
                                        fselect = ftgt   
                                        cselect = ctgt
                                    end if 
                                end if 
                            end if
                        end if 
                    end if 
                end do 

                !Count invalid cells 
                if ((cselect == 0) .AND. (invalid_trigger == 1)) then 
                    Ninvalid = Ninvalid + 1
                end if

                !Tag if valid merge identified
                if (cselect .GT. 0) then !Valid merge found -> Map cell for replacement and tag faces joining ii and ctgt for removal
                    Nmerge = Nmerge + 1
                    cell_map(cc) = cell_map(cselect) 
                    face_remove(fselect) = 1
                    do aa=1,cellnface(cc)
                        ftgt = cell2face(cc,aa)
                        if ((volume_mesh%faces(ftgt)%cleft == cc)) then 
                            ctgt = volume_mesh%faces(ftgt)%cright
                        else
                            ctgt = volume_mesh%faces(ftgt)%cleft
                        end if 
                        if (ctgt == cselect) then 
                            face_remove(ftgt) = 1
                        end if 
                    end do
                else !No valid merge case identified
                    cell_remove(cc) = 2
                    Nmerge_fail = Nmerge_fail + 1
                end if 
            end if 
        end if 
    end do 

    !Update mapping for retained cells 
    do cc=1,volume_mesh%ncell 
        if (cell_remove(cc) == 2) then 
            NcellN = NcellN + 1
            cell_map(cc) = NcellN 
        end if
    end do
    
    !Tag faces that map to the same cell on both sides to be removed
    do ff=1,volume_mesh%nface
        if ((volume_mesh%faces(ff)%cleft .GT. 0) .AND. (volume_mesh%faces(ff)%cright .GT. 0)) then 
            if (cell_map(volume_mesh%faces(ff)%cleft) == cell_map(volume_mesh%faces(ff)%cright)) then 
                face_remove(ff) = 1
            end if
        end if
    end do 

    !Build face mapping for retained faces
    NfaceN = 0 
    face_map(:) = 0 
    do ff=1,volume_mesh%nface
        if (face_remove(ff) .NE. 1) then 
            NfaceN = NfaceN + 1
            face_map(ff) = NfaceN
        end if
    end do 

    !Rebuild mesh
    allocate(volume_mesh_temp%faces(NfaceN))
    allocate(volume_mesh_temp%cell_level(NcellN))
    do cc=1,volume_mesh%ncell 
        if (cell_map(cc) .NE. 0) then 
            volume_mesh_temp%cell_level(cell_map(cc)) = volume_mesh%cell_level(cc)
        end if 
    end do 
    do ff=1,volume_mesh%nface
        if (face_remove(ff) .NE. 1) then 
            volume_mesh_temp%faces(face_map(ff))%nvtx = volume_mesh%faces(ff)%nvtx
            allocate(volume_mesh_temp%faces(face_map(ff))%vertices(volume_mesh%faces(ff)%nvtx))
            volume_mesh_temp%faces(face_map(ff))%vertices(:) = volume_mesh%faces(ff)%vertices(:) 
            if (volume_mesh%faces(ff)%cleft .GT. 0) then 
                volume_mesh_temp%faces(face_map(ff))%cleft = cell_map(volume_mesh%faces(ff)%cleft)
            else
                volume_mesh_temp%faces(face_map(ff))%cleft = volume_mesh%faces(ff)%cleft
            end if 
            if (volume_mesh%faces(ff)%cright .GT. 0) then 
                volume_mesh_temp%faces(face_map(ff))%cright = cell_map(volume_mesh%faces(ff)%cright)
            else
                volume_mesh_temp%faces(face_map(ff))%cright = volume_mesh%faces(ff)%cright
            end if
        end if 
    end do 
    volume_mesh%ncell = NcellN
    volume_mesh%nface = NfaceN
    deallocate(volume_mesh%faces)
    deallocate(volume_mesh%cell_level)
    allocate(volume_mesh%faces(NfaceN))
    allocate(volume_mesh%cell_level(NcellN))
    do ff=1,NfaceN
        volume_mesh%faces(ff)%nvtx = volume_mesh_temp%faces(ff)%nvtx
        allocate(volume_mesh%faces(ff)%vertices(volume_mesh_temp%faces(ff)%nvtx))
        allocate(volume_mesh%faces(ff)%edges(volume_mesh_temp%faces(ff)%nvtx))
        volume_mesh%faces(ff)%vertices(:) = volume_mesh_temp%faces(ff)%vertices(:)
        volume_mesh%faces(ff)%edges(:) = 0 
        volume_mesh%faces(ff)%cleft = volume_mesh_temp%faces(ff)%cleft
        volume_mesh%faces(ff)%cright = volume_mesh_temp%faces(ff)%cright
    end do 
    volume_mesh%cell_level(:) = volume_mesh_temp%cell_level(:)

    !Display
    if (cm3dopt%dispt == 1) then
        write(*,'(A,I0,A,I0,A)') '    {identified ',Ncremove,' sliver cells and merged ',Nmerge,' cells}'
        if (Ninvalid .NE. 0) then 
            write(*,'(A,I0,A)') '    {',Ninvalid,' merges prevented due to double adjacency}'
        end if 
    end if 
end if
return 
end subroutine clean_mesh_sliverC




!Subroutine to flip faces with a target boundary condition ===========================
subroutine flip_set_faces(volume_mesh,fcondition)
implicit none 

!Variables - Import
integer(in) :: fcondition
type(vol_mesh_data) :: volume_mesh

!Variables - Local 
integer(in) :: ff,vv
integer(in) :: fl,fr
integer(in) :: face_vtx_temp(volume_mesh%nface)

!Flip edges with the target condition fcondition as cl or cr 
face_vtx_temp(:) = 0 
do ff=1,volume_mesh%nface
    fl = volume_mesh%faces(ff)%cleft
    fr = volume_mesh%faces(ff)%cright
    if ((fl == fcondition) .OR. (fr == fcondition)) then 
        face_vtx_temp(1:volume_mesh%faces(ff)%nvtx) = volume_mesh%faces(ff)%vertices(:)
        do vv=1,volume_mesh%faces(ff)%nvtx
            volume_mesh%faces(ff)%vertices(vv) = face_vtx_temp(volume_mesh%faces(ff)%nvtx-vv+1)
        end do 
        volume_mesh%faces(ff)%cleft = fr
        volume_mesh%faces(ff)%cright = fl
    end if 
end do 
return 
end subroutine flip_set_faces

    


!Subroutine to simplify the volume mesh surface ===========================
subroutine simplify_surface(volume_mesh,cm3dopt,cm3dfailure)
implicit none 

!Variables - Import
integer(in) :: cm3dfailure
type(cm3d_options) :: cm3dopt
type(vol_mesh_data) :: volume_mesh

!Variables - Local 
integer(in) :: cc,ff,ee,aa,vv,kk,fn,fiter
integer(in) :: nedge_vm,maxvlnce,vtgt,etgt,ftgt,fadded,maxNsf,is_on_surf,is_in_vol,Nsurfcell,csidx,NfaceN_surf,face_ins
integer(in) :: NfaceM,fbase,fadj,nupdate,NfaceN,NedgeinF,evalid,Nvfnew,v_current,v_new,e_current,e_new,e_new1,e_new2,NvtxN
integer(in) :: v1,v2,v1t,v2t,ev1,ev2,fdir
integer(in) :: valence(volume_mesh%nvtx),surf_cell_idx(volume_mesh%ncell),sface_cidx(volume_mesh%nface)
integer(in) :: face_new_acc(volume_mesh%nvtx),vtxmap(volume_mesh%nvtx)
integer(in) :: V2E_fmerge(volume_mesh%nvtx,2)
integer(in), dimension(:), allocatable :: edge_svol,cellNsface,fmerge,edge_in_face,edge_in_face_visited,edge_fidx
integer(in), dimension(:,:), allocatable :: cell2face,vmf_edges,E2F
real(dp) :: mincurv,curvtest,reflen
real(dp) :: vertices_temp(volume_mesh%nvtx,3)
type(face_data) :: faces_surf_new(volume_mesh%nface),faces_temp(volume_mesh%nface)

!Build vmesh edges 
call get_vmesh_valence(valence,volume_mesh)
maxvlnce = maxval(valence(:))
call build_vmesh_edges(nedge_vm,vmf_edges,volume_mesh,maxvlnce)

!Build E2F 
allocate(E2F(nedge_vm,6))
E2F(:,:) = 0 
do ff=1,volume_mesh%nface 
    do ee=1,volume_mesh%faces(ff)%nvtx
        etgt = volume_mesh%faces(ff)%edges(ee)
        fadded = 0 
        do aa=1,6
            if (E2F(etgt,aa) == 0) then 
                E2F(etgt,aa) = ff 
                fadded = 1
                exit 
            end if   
        end do 
        if (fadded == 0) then 
            print *, '** E2F valence exceeded '
            cm3dfailure = 1
            return 
        end if 
    end do 
end do 

!Tag all surface edges that are also in a volume face
allocate(edge_svol(nedge_vm))
edge_svol(:) = 0 
do ee=1,nedge_vm
    is_in_vol = 0 
    is_on_surf = 0 
    do aa=1,6
        ftgt = E2F(ee,aa)
        if (ftgt .GT. 0) then
            if ((volume_mesh%faces(ftgt)%cleft == -1) .OR. (volume_mesh%faces(ftgt)%cright == -1)) then !is surface face
                is_on_surf = 1
            end if 
            if ((volume_mesh%faces(ftgt)%cleft .GT. 0) .AND. (volume_mesh%faces(ftgt)%cright .GT. 0)) then !not surface face
                is_in_vol = 1
            end if
            if ((volume_mesh%faces(ftgt)%cleft .LT. -1) .AND. (volume_mesh%faces(ftgt)%cright .LT. -1)) then !not surface face
                is_in_vol = 1
            end if
        end if 
    end do 
    if ((is_on_surf == 1) .AND. (is_in_vol == 1)) then 
        edge_svol(ee) = 1
    end if 
end do 

!Identify surface adjacent mesh cells
Nsurfcell = 0 
surf_cell_idx(:) = 0 
do ff=1,volume_mesh%nface 
    if ((volume_mesh%faces(ff)%cleft == -1) .OR. (volume_mesh%faces(ff)%cright == -1)) then !is surface face
        if (volume_mesh%faces(ff)%cleft .GT. 0) then 
            if (surf_cell_idx(volume_mesh%faces(ff)%cleft) == 0) then 
                Nsurfcell = Nsurfcell + 1
                surf_cell_idx(volume_mesh%faces(ff)%cleft) = Nsurfcell
            end if 
        end if
        if (volume_mesh%faces(ff)%cright .GT. 0) then 
            if (surf_cell_idx(volume_mesh%faces(ff)%cright) == 0) then 
                Nsurfcell = Nsurfcell + 1
                surf_cell_idx(volume_mesh%faces(ff)%cright) = Nsurfcell
            end if 
        end if
    end if
end do 

!Assign surface mesh faces to each cell
allocate(cellNsface(Nsurfcell))
cellNsface(:) = 0 
do ff=1,volume_mesh%nface 
    if ((volume_mesh%faces(ff)%cleft == -1) .OR. (volume_mesh%faces(ff)%cright == -1)) then !is surface face
        if (volume_mesh%faces(ff)%cleft .GT. 0) then 
            cellNsface(surf_cell_idx(volume_mesh%faces(ff)%cleft)) = cellNsface(surf_cell_idx(volume_mesh%faces(ff)%cleft)) + 1
        end if 
        if (volume_mesh%faces(ff)%cright .GT. 0) then 
            cellNsface(surf_cell_idx(volume_mesh%faces(ff)%cright)) = cellNsface(surf_cell_idx(volume_mesh%faces(ff)%cright)) + 1
        end if 
    end if 
end do 
maxNsf = maxval(cellNsface(:))
allocate(cell2face(Nsurfcell,maxNsf))
sface_cidx(:) = 0 
cellNsface(:) = 0 
cell2face(:,:) = 0 
do ff=1,volume_mesh%nface 
    if ((volume_mesh%faces(ff)%cleft == -1) .OR. (volume_mesh%faces(ff)%cright == -1)) then !is surface face
        if (volume_mesh%faces(ff)%cleft .GT. 0) then 
            cellNsface(surf_cell_idx(volume_mesh%faces(ff)%cleft)) = cellNsface(surf_cell_idx(volume_mesh%faces(ff)%cleft)) + 1
            cell2face(surf_cell_idx(volume_mesh%faces(ff)%cleft),cellNsface(surf_cell_idx(volume_mesh%faces(ff)%cleft))) = ff 
            sface_cidx(ff) = cellNsface(surf_cell_idx(volume_mesh%faces(ff)%cleft))
        end if 
        if (volume_mesh%faces(ff)%cright .GT. 0) then 
            cellNsface(surf_cell_idx(volume_mesh%faces(ff)%cright)) = cellNsface(surf_cell_idx(volume_mesh%faces(ff)%cright)) + 1
            cell2face(surf_cell_idx(volume_mesh%faces(ff)%cright),cellNsface(surf_cell_idx(volume_mesh%faces(ff)%cright))) = ff 
            sface_cidx(ff) = cellNsface(surf_cell_idx(volume_mesh%faces(ff)%cright))
        end if 
    end if
end do 

!Filter surface cells for simplification when there is a small feature in the cell 
if (cm3dopt%surf_force_simplify == 0) then 
    do cc=1,volume_mesh%ncell 
        if (surf_cell_idx(cc) .GT. 0) then 

            !Surface cell index
            csidx = surf_cell_idx(cc)

            !Check curvature
            mincurv = ieee_value(1.0d0,IEEE_POSITIVE_INF)
            do ff=1,cellNsface(csidx)
                ftgt = cell2face(csidx,ff)
                do vv=1,volume_mesh%faces(ftgt)%nvtx
                    vtgt = volume_mesh%faces(ftgt)%vertices(vv)
                    curvtest = volume_mesh%vtx_sRcurv(vtgt)
                    if (curvtest .LT. mincurv) then 
                        mincurv = curvtest
                    end if 
                end do 
            end do

            !If minimum curvature is less than the reference length then dont simplify this cell -> set negative surface index so it is not simplified 
            reflen = 2.0d0*cm3dopt%far_field_bound/(2.0d0**(volume_mesh%cell_level(cc) - 1))
            if (mincurv .LE. reflen) then 
                surf_cell_idx(cc) = -surf_cell_idx(cc) 
            end if 
        end if
    end do 
end if 

!Group surface faces in each cell into new faces divided only by volume face attached edges 
NfaceN = 0 
face_new_acc(:) = 0 
V2E_fmerge(:,:) = 0 
allocate(edge_fidx(nedge_vm))
allocate(edge_in_face(nedge_vm)) 
allocate(edge_in_face_visited(nedge_vm)) 
edge_fidx(:) = 0
edge_in_face(:) = 0
edge_in_face_visited(:) = 0
do cc=1,volume_mesh%ncell 
    if (surf_cell_idx(cc) .GT. 0) then !simplify surface

        !Surface cell index
        csidx = surf_cell_idx(cc)

        !Allocate face merging array 
        allocate(fmerge(cellNsface(csidx)))
        fmerge(:) = 0 

        !Flood face merges
        NfaceM = 0 
        do kk=1,cellNsface(csidx)

            !Find untagged surface face
            fbase = 0 
            do ff=1,cellNsface(csidx)
                ftgt = cell2face(csidx,ff)
                ftgt = sface_cidx(ftgt)
                if (fmerge(ftgt) == 0) then 
                    fbase = ftgt 
                    exit 
                end if 
            end do 

            !If no faces found then exit as all are tagged 
            if (fbase == 0) then 
                exit 
            end if 

            !Increment contained face count
            NfaceM = NfaceM + 1 

            !Tag base face with the face index
            fmerge(fbase) = NfaceM

            !Flood face tag across all non volume attached edges in this cell 
            do fiter=1,4*cellNsface(csidx)

                !Flood on valid edges on all faces in this cell
                nupdate = 0 
                do ff=1,cellNsface(csidx)
                    ftgt = cell2face(csidx,ff)
                    if (fmerge(sface_cidx(ftgt)) .NE. 0) then !if tagged with a merged face
                        do ee=1,volume_mesh%faces(ftgt)%nvtx !edges on this face
                            etgt = volume_mesh%faces(ftgt)%edges(ee) 
                            if (edge_svol(etgt) == 0) then !if not attached to a volume face
                                do aa=1,6 !flood to faces on this edge  
                                    if (E2F(etgt,aa) .NE. 0) then 
                                        fadj = E2F(etgt,aa)
                                        if ((fadj .NE. ftgt) .AND. (fadj .GT. 0))then 
                                            if (sface_cidx(fadj) .GT. 0) then 
                                                if (fmerge(sface_cidx(fadj)) == 0) then 
                                                    fmerge(sface_cidx(fadj)) = NfaceM
                                                    nupdate = nupdate + 1 
                                                end if 
                                            end if 
                                        end if 
                                    end if 
                                end do 
                            end if 
                        end do 
                    end if 
                end do 

                !Exit at no updates
                if (nupdate == 0) then 
                    exit 
                end if 
            end do 
        end do 

        !Gather edges of each merged face to build new faces 
        do fn=1,NfaceM

            !Initialise 
            NedgeinF = 0 

            !Gather edges to this face
            do ff=1,cellNsface(csidx)

                !Face
                ftgt = cell2face(csidx,ff)

                !If within this merged face
                if (fmerge(sface_cidx(ftgt)) == fn) then 

                    !Accumulate each edge to this merged face that is adjacent to a different merged face
                    do ee=1,volume_mesh%faces(ftgt)%nvtx
                        
                        !Edge
                        etgt = volume_mesh%faces(ftgt)%edges(ee) 

                        !Test if valid 
                        evalid = 0 
                        if (edge_svol(etgt) == 1) then !is on a volume face
                            evalid = 1
                        else
                            do aa=1,6 
                                if (E2F(etgt,aa) .NE. 0) then 
                                    fadj = E2F(etgt,aa)
                                    if (fadj .GT. 0) then 
                                        if (fadj .NE. ftgt)then 
                                            if (sface_cidx(fadj) .GT. 0) then !if a surface face in this cell 
                                                if (fmerge(sface_cidx(fadj)) .NE. fmerge(sface_cidx(ftgt))) then !is edge to different merged face
                                                    evalid = 1 
                                                    exit 
                                                end if 
                                            elseif ((volume_mesh%faces(fadj)%cleft .LT. -1) .OR. &
                                                    (volume_mesh%faces(fadj)%cright .LT. -1)) then !is a global mesh boundary condition
                                                evalid = 1 
                                                exit 
                                            end if 
                                        end if 
                                    end if 
                                end if 
                            end do 
                        end if 

                        !Add edge if valid
                        if (evalid == 1) then 
                            if (edge_fidx(etgt) == 0) then 
                                NedgeinF = NedgeinF + 1
                                edge_fidx(etgt) = NedgeinF
                                edge_in_face(NedgeinF) = etgt 
                            end if 
                        end if 
                    end do 
                end if 
            end do 

            !Debug print link data ====================
            ! print *,'cell = ', cc
            ! open(11,file='io/fedges')
            !     do ee=1,NedgeinF    
            !         v1 = vmf_edges(edge_in_face(ee),1)
            !         v2 = vmf_edges(edge_in_face(ee),2)
            !         write(11,*) v1,v2
            !     end do 
            ! close(11)
            ! open(11,file='io/f_faces')
            !     do ff=1,cellNsface(csidx) 
            !         ftgt = cell2face(csidx,ff)
            !         write(11,*) ftgt
            !     end do 
            ! close(11)
            !Debug print link data ====================

            !Build v2e for these edges
            do ee=1,NedgeinF
                v1 = vmf_edges(edge_in_face(ee),1)
                v2 = vmf_edges(edge_in_face(ee),2)
                if (V2E_fmerge(v1,1) == 0) then 
                    V2E_fmerge(v1,1) = ee 
                else
                    V2E_fmerge(v1,2) = ee 
                end if 
                if (V2E_fmerge(v2,1) == 0) then 
                    V2E_fmerge(v2,1) = ee 
                else
                    V2E_fmerge(v2,2) = ee 
                end if 
            end do 

            !Accumulate face 
            v_current = vmf_edges(edge_in_face(1),1)
            e_current = 1 
            edge_in_face_visited(e_current) = 1
            Nvfnew = 1 
            face_new_acc(1) = v_current
            do ee=1,NedgeinF

                !Find new edge 
                e_new1 = V2E_fmerge(v_current,1)
                e_new2 = V2E_fmerge(v_current,2)
                if (edge_in_face_visited(e_new1) == 0) then !e_new1 new
                    e_new = e_new1
                elseif (edge_in_face_visited(e_new2) == 0) then !e_new2 new
                    e_new = e_new2
                else !Both visited so the face is closed so exit
                    exit 
                end if 

                !Tag edge as visited
                edge_in_face_visited(e_new) = 1

                !Find new vertex
                if (vmf_edges(edge_in_face(e_new),1) == v_current) then 
                    v_new = vmf_edges(edge_in_face(e_new),2) 
                else
                    v_new = vmf_edges(edge_in_face(e_new),1) 
                end if 

                !Add new vertex to face
                Nvfnew = Nvfnew + 1
                face_new_acc(Nvfnew) = v_new

                !Update current state
                e_current = e_new 
                v_current = v_new
            end do 

            !Reset linking arrays
            edge_in_face_visited(1:NedgeinF) = 0 
            do ee=1,NedgeinF
                v1 = vmf_edges(edge_in_face(ee),1)
                v2 = vmf_edges(edge_in_face(ee),2)
                V2E_fmerge(v1,:) = 0 
                V2E_fmerge(v2,:) = 0 
            end do 

            !Check face is aligned correctly by ensuring one sub-edge is in the same direction as the merged face winding order
            v1 = face_new_acc(1)
            v2 = face_new_acc(2)
            etgt = 0 
            ftgt = 0 
            fdir = 0 
            do ee=1,NedgeinF
                if ((vmf_edges(edge_in_face(ee),1) == v1) .AND. (vmf_edges(edge_in_face(ee),2) == v2)) then 
                    etgt = edge_in_face(ee)
                    exit 
                elseif ((vmf_edges(edge_in_face(ee),1) == v2) .AND. (vmf_edges(edge_in_face(ee),2) == v1)) then 
                    etgt = edge_in_face(ee)
                    exit 
                end if
            end do 
            do aa=1,6 
                if (E2F(etgt,aa) .NE. 0) then 
                    fadj = E2F(etgt,aa)
                    if (fadj .GT. 0)then 
                        if ((volume_mesh%faces(fadj)%cleft == -1) .OR. (volume_mesh%faces(fadj)%cright == -1)) then !is surface face
                            if ((volume_mesh%faces(fadj)%cleft == cc) .OR. (volume_mesh%faces(fadj)%cright == cc)) then !is in this cell cc
                                ftgt = fadj 
                            end if 
                        end if 
                    end if 
                end if 
            end do 
            do ee=1,volume_mesh%faces(ftgt)%nvtx
                ev1 = ee
                ev2 = mod(ee,volume_mesh%faces(ftgt)%nvtx) + 1
                v1t = volume_mesh%faces(ftgt)%vertices(ev1)
                v2t = volume_mesh%faces(ftgt)%vertices(ev2)
                if ((v1 == v1t) .AND. (v2 == v2t)) then 
                    fdir = 1
                    exit 
                elseif ((v1 == v2t) .AND. (v2 == v1t)) then 
                    fdir = -1 
                    exit 
                end if 
            end do 
            
            !Store new face
            NfaceN = NfaceN + 1
            faces_surf_new(NfaceN)%nvtx = Nvfnew
            allocate(faces_surf_new(NfaceN)%vertices(Nvfnew))
            if (fdir == 1) then 
                faces_surf_new(NfaceN)%vertices(:) = face_new_acc(1:Nvfnew)
            elseif (fdir == -1) then 
                do ee=1,Nvfnew
                    faces_surf_new(NfaceN)%vertices(ee) = face_new_acc(Nvfnew-ee+1)
                end do 
            else    
                print *, '** ambiguous simplified face direction'
            end if 
            if (cm3dopt%meshinout == 'out') then !mesh external 
                faces_surf_new(NfaceN)%cleft = -1 
                faces_surf_new(NfaceN)%cright = cc 
            elseif (cm3dopt%meshinout == 'in') then !mesh internal 
                faces_surf_new(NfaceN)%cleft = cc 
                faces_surf_new(NfaceN)%cright = -1
            end if 

            !Reset edge tags 
            edge_fidx(edge_in_face(1:NedgeinF)) = 0 
            edge_in_face(1:NedgeinF) = 0 
        end do 

        !Deallocate face merging array 
        deallocate(fmerge)
    elseif (surf_cell_idx(cc) .LT. 0) then !retain triangulated surface 

        !Surface cell index
        csidx = abs(surf_cell_idx(cc))

        !Copy all faces in this cell
        do ff=1,cellNsface(csidx)

            !Target base face
            ftgt = cell2face(csidx,ff)

            !Store new face
            NfaceN = NfaceN + 1
            faces_surf_new(NfaceN)%nvtx = volume_mesh%faces(ftgt)%nvtx
            allocate(faces_surf_new(NfaceN)%vertices(faces_surf_new(NfaceN)%nvtx))
            faces_surf_new(NfaceN)%vertices(:) = volume_mesh%faces(ftgt)%vertices(:)
            faces_surf_new(NfaceN)%cleft = volume_mesh%faces(ftgt)%cleft
            faces_surf_new(NfaceN)%cright = volume_mesh%faces(ftgt)%cright
        end do 
    end if 
end do 
NfaceN_surf = NfaceN

!Count final number of new faces 
do ff=1,volume_mesh%nface 
    if ((volume_mesh%faces(ff)%cleft .NE. -1) .AND. (volume_mesh%faces(ff)%cright .NE. -1)) then !not surface face
        if (volume_mesh%faces(ff)%nvtx .NE. 0) then 
            NfaceN = NfaceN + 1
        end if 
    end if
end do 

!Update volume mesh faces 
face_ins = 0 
do ff=1,volume_mesh%nface 
    faces_temp(ff)%nvtx = volume_mesh%faces(ff)%nvtx
    faces_temp(ff)%cleft = volume_mesh%faces(ff)%cleft
    faces_temp(ff)%cright = volume_mesh%faces(ff)%cright
    allocate(faces_temp(ff)%vertices(faces_temp(ff)%nvtx))
    faces_temp(ff)%vertices(:) = volume_mesh%faces(ff)%vertices(:)
    deallocate(volume_mesh%faces(ff)%vertices)
end do 
deallocate(volume_mesh%faces)
allocate(volume_mesh%faces(NfaceN))
do ff=1,volume_mesh%nface 
    if ((faces_temp(ff)%cleft .NE. -1) .AND. (faces_temp(ff)%cright .NE. -1)) then !not surface face
        if (faces_temp(ff)%nvtx .NE. 0) then 
            face_ins = face_ins + 1
            volume_mesh%faces(face_ins)%nvtx = faces_temp(ff)%nvtx 
            volume_mesh%faces(face_ins)%cleft = faces_temp(ff)%cleft
            volume_mesh%faces(face_ins)%cright = faces_temp(ff)%cright
            allocate(volume_mesh%faces(face_ins)%vertices(volume_mesh%faces(face_ins)%nvtx))
            allocate(volume_mesh%faces(face_ins)%edges(volume_mesh%faces(face_ins)%nvtx))
            volume_mesh%faces(face_ins)%vertices(:) = faces_temp(ff)%vertices(:)
            volume_mesh%faces(face_ins)%edges(:) = 0 
        end if 
    end if 
end do 
do ff=1,NfaceN_surf
    face_ins = face_ins + 1
    volume_mesh%faces(face_ins)%nvtx = faces_surf_new(ff)%nvtx 
    volume_mesh%faces(face_ins)%cleft = faces_surf_new(ff)%cleft
    volume_mesh%faces(face_ins)%cright = faces_surf_new(ff)%cright
    allocate(volume_mesh%faces(face_ins)%vertices(volume_mesh%faces(face_ins)%nvtx))
    allocate(volume_mesh%faces(face_ins)%edges(volume_mesh%faces(face_ins)%nvtx))
    volume_mesh%faces(face_ins)%vertices(:) = faces_surf_new(ff)%vertices(:)
    volume_mesh%faces(face_ins)%edges(:) = 0 
end do 
volume_mesh%nface = face_ins

!Update mesh vertices 
NvtxN = 0
vtxmap(:) = 0 
do ff=1,volume_mesh%nface 
    do vv=1,volume_mesh%faces(ff)%nvtx
        if (vtxmap(volume_mesh%faces(ff)%vertices(vv)) == 0) then 
            NvtxN = NvtxN + 1 
            vtxmap(volume_mesh%faces(ff)%vertices(vv)) = NvtxN
        end if 
    end do 
end do 
vertices_temp(:,:) = volume_mesh%vtx(:,:)
deallocate(volume_mesh%vtx)
allocate(volume_mesh%vtx(NvtxN,3))
do vv=1,volume_mesh%nvtx
    if (vtxmap(vv) .NE. 0) then 
        volume_mesh%vtx(vtxmap(vv),:) = vertices_temp(vv,:)
    end if 
end do 
volume_mesh%nvtx = NvtxN

!Remap faces
do ff=1,volume_mesh%nface 
    do vv=1,volume_mesh%faces(ff)%nvtx
        volume_mesh%faces(ff)%vertices(vv) = vtxmap(volume_mesh%faces(ff)%vertices(vv))
    end do  
end do 
return 
end subroutine simplify_surface




!Subroutine to remove mesh surface or internal valence two vertices by merging adjacent edges ===========================
subroutine clean_vlnc2_vertices(volume_mesh,cm3dopt,surf_int)
implicit none 

!Variables - Import
integer(in) :: surf_int
type(cm3d_options) :: cm3dopt
type(vol_mesh_data) :: volume_mesh

!Variables - Local 
integer(in) :: ff,vv,ee
integer(in) :: maxvalence,face_nvtxn,nvtx_rem,nvtxN
integer(in) :: vtx_bndry(volume_mesh%nvtx),valence(volume_mesh%nvtx),vtx_rem(volume_mesh%nvtx)
integer(in) :: fvtx_temp(volume_mesh%nvtx),vtx_map(volume_mesh%nvtx)
real(dp) :: vertices_temp(volume_mesh%nvtx,3)

!Approximate maximum valence 
valence(:) = 0
do ff=1,volume_mesh%nface 
    do vv=1,volume_mesh%faces(ff)%nvtx
        valence(volume_mesh%faces(ff)%vertices(vv)) = valence(volume_mesh%faces(ff)%vertices(vv)) + 1
    end do   
end do 
maxvalence = maxval(valence)

!Build mesh edges 
if (allocated(volume_mesh%edges)) then 
    volume_mesh%nedge = 0 
    deallocate(volume_mesh%edges)
end if 
call build_vmesh_edges(volume_mesh%nedge,volume_mesh%edges,volume_mesh,maxvalence)

!Evaluate true valence 
valence(:) = 0
do ee=1,volume_mesh%nedge
    valence(volume_mesh%edges(ee,:)) = valence(volume_mesh%edges(ee,:)) + 1
end do 

!Tag surface vertices 
vtx_bndry(:) = 0
do ff=1,volume_mesh%nface 
    if ((volume_mesh%faces(ff)%cleft .LT. 0) .OR. (volume_mesh%faces(ff)%cright .LT. 0)) then
        do vv=1,volume_mesh%faces(ff)%nvtx
            vtx_bndry(volume_mesh%faces(ff)%vertices(vv)) = 1
        end do 
    end if
end do 

!Tag vertices that are valid for removal 
nvtx_rem = 0 
vtx_rem(:) = 0 
do vv=1,volume_mesh%nvtx
    if (surf_int == 1) then !remove surface vertices
        if ((vtx_bndry(vv) == 1) .AND. (valence(vv) == 2)) then 
            vtx_rem(vv) = 1
            nvtx_rem = nvtx_rem + 1
        end if 
    elseif (surf_int == 0) then !remove volume vertices
        if ((vtx_bndry(vv) == 0) .AND. (valence(vv) == 2)) then 
            vtx_rem(vv) = 1
            nvtx_rem = nvtx_rem + 1
        end if 
    end if 
end do 

!Update each face by removing tagged vertices from the faces
fvtx_temp(:) = 0 
do ff=1,volume_mesh%nface 
    face_nvtxn = 0 
    do vv=1,volume_mesh%faces(ff)%nvtx
        if (vtx_rem(volume_mesh%faces(ff)%vertices(vv)) == 0) then 
            face_nvtxn = face_nvtxn + 1
            fvtx_temp(face_nvtxn) = volume_mesh%faces(ff)%vertices(vv)
        end if 
    end do 
    deallocate(volume_mesh%faces(ff)%vertices)
    deallocate(volume_mesh%faces(ff)%edges)
    volume_mesh%faces(ff)%nvtx = face_nvtxn
    allocate(volume_mesh%faces(ff)%vertices(face_nvtxn))
    allocate(volume_mesh%faces(ff)%edges(face_nvtxn))
    volume_mesh%faces(ff)%vertices(:) = fvtx_temp(1:face_nvtxn) 
    volume_mesh%faces(ff)%edges(:) = 0 
end do 

!Build vertex map 
nvtxN = 0 
vtx_map(:) = 0 
do vv=1,volume_mesh%nvtx
    if (vtx_rem(vv) == 0) then 
        nvtxN = nvtxN + 1
        vtx_map(vv) = nvtxN
    end if 
end do 

!Update vertex list eliminating all removed vertices 
vertices_temp(:,:) = volume_mesh%vtx(:,:)
deallocate(volume_mesh%vtx)
allocate(volume_mesh%vtx(nvtxN,3))
do vv=1,volume_mesh%nvtx
    if (vtx_map(vv) .NE. 0) then 
        volume_mesh%vtx(vtx_map(vv),:) = vertices_temp(vv,:)
    end if 
end do 
volume_mesh%nvtx = nvtxN
do ff=1,volume_mesh%nface 
    do vv=1,volume_mesh%faces(ff)%nvtx
        volume_mesh%faces(ff)%vertices(vv) = vtx_map(volume_mesh%faces(ff)%vertices(vv))
    end do 
end do 

!Display 
if (cm3dopt%dispt == 1) then
    if (surf_int == 1) then !remove surface vertices
        write(*,'(A,I0,A)') '    {collapsed ',nvtx_rem,' surface valence two vertices}'
    elseif (surf_int == 0) then !remove volume vertices    
        write(*,'(A,I0,A)') '    {collapsed ',nvtx_rem,' internal valence two vertices}'
    end if 
end if 
return 
end subroutine clean_vlnc2_vertices




!Subroutine to construct dual mesh ===========================
subroutine construct_dual_mesh(volume_mesh)
implicit none 

!Variables - Import
type(vol_mesh_data) :: volume_mesh

!Variables - Local 
integer(in) :: ff,vv,cc,ee
integer(in) :: maxfvalence,maxvalence,maxebcf,Ncnew,Nvnew,BCbase
integer(in) :: vtx_bndry(volume_mesh%nvtx),fvalence(volume_mesh%nvtx),valence(volume_mesh%nvtx)
integer(in) :: cell_vtxidx(volume_mesh%ncell),face_vtxidx(volume_mesh%nface),vtx_cidx(volume_mesh%nvtx)
integer(in), dimension(:), allocatable :: edge_vtxidx,nBCfedge
integer(in), dimension(:,:), allocatable :: V2F,E2F

!Tag surface vertices 
vtx_bndry(:) = 0
do ff=1,volume_mesh%nface 
    if ((volume_mesh%faces(ff)%cleft .LT. 0) .OR. (volume_mesh%faces(ff)%cright .LT. 0)) then
        do vv=1,volume_mesh%faces(ff)%nvtx
            vtx_bndry(volume_mesh%faces(ff)%vertices(vv)) = 1
        end do 
    end if
end do 

!Build V2F for boundary condition vertices and faces
fvalence(:) = 0 
do ff=1,volume_mesh%nface 
    if ((volume_mesh%faces(ff)%cleft .LT. 0) .OR. (volume_mesh%faces(ff)%cright .LT. 0)) then
        do vv=1,volume_mesh%faces(ff)%nvtx
            fvalence(volume_mesh%faces(ff)%vertices(vv)) = fvalence(volume_mesh%faces(ff)%vertices(vv)) + 1
        end do 
    end if 
end do 
maxfvalence = maxval(fvalence(:))
allocate(V2F(volume_mesh%nvtx,maxfvalence))
V2F(:,:) = 0 
fvalence(:) = 0 
do ff=1,volume_mesh%nface 
    if ((volume_mesh%faces(ff)%cleft .LT. 0) .OR. (volume_mesh%faces(ff)%cright .LT. 0)) then
        do vv=1,volume_mesh%faces(ff)%nvtx
            fvalence(volume_mesh%faces(ff)%vertices(vv)) = fvalence(volume_mesh%faces(ff)%vertices(vv)) + 1
            V2F(volume_mesh%faces(ff)%vertices(vv),fvalence(volume_mesh%faces(ff)%vertices(vv))) = ff
        end do 
    end if 
end do 

!Estimate mesh valence 
valence(:) = 0
do ff=1,volume_mesh%nface 
    do vv=1,volume_mesh%faces(ff)%nvtx
        valence(volume_mesh%faces(ff)%vertices(vv)) = valence(volume_mesh%faces(ff)%vertices(vv)) + 1
    end do   
end do 
maxvalence = maxval(valence(:))

!Build mesh edges 
if (allocated(volume_mesh%edges)) then 
    volume_mesh%nedge = 0 
    deallocate(volume_mesh%edges)
end if 
call build_vmesh_edges(volume_mesh%nedge,volume_mesh%edges,volume_mesh,maxvalence)

!Build E2F for boundary condition edges and faces
allocate(nBCfedge(volume_mesh%nedge))
nBCfedge(:) = 0 
do ff=1,volume_mesh%nface 
    if ((volume_mesh%faces(ff)%cleft .LT. 0) .OR. (volume_mesh%faces(ff)%cright .LT. 0)) then
        do vv=1,volume_mesh%faces(ff)%nvtx
            nBCfedge(volume_mesh%faces(ff)%edges(vv)) = nBCfedge(volume_mesh%faces(ff)%edges(vv)) + 1 
        end do 
    end if
end do 
maxebcf = maxval(nBCfedge)
allocate(E2F(volume_mesh%nedge,maxebcf))
E2F(:,:) = 0 
nBCfedge(:) = 0 
do ff=1,volume_mesh%nface 
    if ((volume_mesh%faces(ff)%cleft .LT. 0) .OR. (volume_mesh%faces(ff)%cright .LT. 0)) then
        do vv=1,volume_mesh%faces(ff)%nvtx
            nBCfedge(volume_mesh%faces(ff)%edges(vv)) = nBCfedge(volume_mesh%faces(ff)%edges(vv)) + 1 
            E2F(volume_mesh%faces(ff)%edges(vv),nBCfedge(volume_mesh%faces(ff)%edges(vv))) = ff
        end do 
    end if 
end do

!Index new cells at mesh vertex locations 
Ncnew = 0 
vtx_cidx(:) = 0 
do vv=1,volume_mesh%nvtx
    Ncnew = Ncnew + 1
    vtx_cidx(vv) = Ncnew
end do 

!Index new vertices at edge, face and cell locations 
Nvnew = 0 
cell_vtxidx(:) = 0 
face_vtxidx(:) = 0 
allocate(edge_vtxidx(volume_mesh%nedge))
edge_vtxidx(:) = 0 
do cc=1,volume_mesh%ncell
    Nvnew = Nvnew + 1
    cell_vtxidx(cc) = Nvnew
end do 
do ff=1,volume_mesh%nface
    if ((volume_mesh%faces(ff)%cleft .LT. 0) .OR. (volume_mesh%faces(ff)%cright .LT. 0)) then
        Nvnew = Nvnew + 1
        face_vtxidx(ff) = Nvnew
    end if
end do 
do ee=1,volume_mesh%nedge
    if (nBCfedge(ee) .GE. 2) then 
        BCbase = 0
        do ff=2,nBCfedge(ee)
            if (volume_mesh%faces(E2F(ee,ff))%cleft .LT. 0) then 
                BCbase = volume_mesh%faces(E2F(ee,ff))%cleft
                exit 
            elseif (volume_mesh%faces(E2F(ee,ff))%cright .LT. 0) then 
                BCbase = volume_mesh%faces(E2F(ee,ff))%cright
                exit 
            end if 
        end do  


        


        do ff=2,nBCfedge(ee)
        

        end do 
    end if 
end do 






return 
end subroutine construct_dual_mesh


end module cellmesh3d_postprocess_mod