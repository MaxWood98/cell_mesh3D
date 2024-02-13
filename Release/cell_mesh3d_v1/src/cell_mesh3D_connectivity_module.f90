!Cell Mesh 3D Mesh Connectivity Module
!Max Wood - mw16116@bristol.ac.uk
!Univeristy of Bristol - Department of Aerospace Engineering

!Version 1.1
!Updated 19-10-2023

!Module
module cellmesh3d_connectivity_mod
use cellmesh3d_data_mod
contains 

!Subroutine to evaluate valence in a triangular mesh ===========================
subroutine get_tri_valence(valence,maxValence,faces,nface,nvtx)
implicit none 

!Variables - Import
integer(in) :: maxValence,nface,nvtx
integer(in), dimension(:), allocatable :: valence
integer(in), dimension(:,:) :: faces

!Variables - Local
integer(in) :: ff,ee,vv
integer(in) :: ev1,ev2,ubValence,evalid
integer(in), dimension(:,:), allocatable :: vconnect

!Initialse valence array
if (allocated(valence)) then 
    deallocate(valence)
end if 
allocate(valence(nvtx))
valence(:) = 0 

!Upper bound of maximum valence 
valence(:) = 0 
do ff=1,nface
    do ee=1,3

        !Edge end vertices
        ev1 = ee
        ev2 = mod(ee,3) + 1
        ev1 = faces(ff,ev1) 
        ev2 = faces(ff,ev2) 

        !Accumulate valence 
        valence(ev1) = valence(ev1) + 1
        valence(ev2) = valence(ev2) + 1
    end do 
end do 
ubValence = 2*maxval(valence)

!Construct actual valence of each vertex
allocate(vconnect(nvtx,ubValence))
vconnect(:,:) = 0 
valence(:) = 0
do ff=1,nface
    do ee=1,3

        !Edge end vertices
        ev1 = ee
        ev2 = mod(ee,3) + 1
        ev1 = faces(ff,ev1) 
        ev2 = faces(ff,ev2) 

        !Check against vconnect 
        evalid = 1
        do vv=1,ubValence
            if (vconnect(ev1,vv) == ev2) then 
                evalid = 0
                exit
            end if 
            if (vconnect(ev2,vv) == ev1) then 
                evalid = 0
                exit 
            end if 
        end do 

        !Add valence if new edge
        if (evalid == 1) then 
            
            !Increment valence on each vertex
            valence(ev1) = valence(ev1) + 1
            valence(ev2) = valence(ev2) + 1

            !Update vconnect
            do vv=1,ubValence
                if (vconnect(ev1,vv) == 0) then 
                    vconnect(ev1,vv) = ev2 
                    exit 
                end if 
            end do 
            do vv=1,ubValence
                if (vconnect(ev2,vv) == 0) then 
                    vconnect(ev2,vv) = ev1 
                    exit 
                end if 
            end do 
        end if 
    end do 
end do 

!Set maximum valence
maxValence = maxval(valence(:))
return 
end subroutine get_tri_valence




!Subroutine to evaluate valence in a volume mesh structure ===========================
subroutine get_vm_valence(valence,maxValence,volume_mesh)
implicit none 

!Variables - Import
integer(in) :: maxValence
integer(in), dimension(:), allocatable :: valence
type(vol_mesh_data) :: volume_mesh

!Variables - Local
integer(in) :: ff,ee,vv
integer(in) :: ev1,ev2,ubValence,evalid
integer(in), dimension(:,:), allocatable :: vconnect

!Initialse valence array
if (allocated(valence)) then 
    deallocate(valence)
end if 
allocate(valence(volume_mesh%nvtx))
valence(:) = 0 

!Upper bound of maximum valence 
valence(:) = 0 
do ff=1,volume_mesh%nface
    do ee=1,volume_mesh%faces(ff)%nvtx

        !Edge end vertices
        ev1 = ee
        ev2 = mod(ee,volume_mesh%faces(ff)%nvtx) + 1
        ev1 = volume_mesh%faces(ff)%vertices(ev1) 
        ev2 = volume_mesh%faces(ff)%vertices(ev2) 

        !Accumulate valence 
        valence(ev1) = valence(ev1) + 1
        valence(ev2) = valence(ev2) + 1
    end do 
end do 
ubValence = 2*maxval(valence)

!Construct actual valence of each vertex
allocate(vconnect(volume_mesh%nvtx,ubValence))
vconnect(:,:) = 0 
valence(:) = 0
do ff=1,volume_mesh%nface
    do ee=1,volume_mesh%faces(ff)%nvtx

        !Edge end vertices
        ev1 = ee
        ev2 = mod(ee,volume_mesh%faces(ff)%nvtx) + 1
        ev1 = volume_mesh%faces(ff)%vertices(ev1) 
        ev2 = volume_mesh%faces(ff)%vertices(ev2)  

        !Check against vconnect 
        evalid = 1
        do vv=1,ubValence
            if (vconnect(ev1,vv) == ev2) then 
                evalid = 0
                exit
            end if 
            if (vconnect(ev2,vv) == ev1) then 
                evalid = 0
                exit 
            end if 
        end do 

        !Add valence if new edge
        if (evalid == 1) then 
            
            !Increment valence on each vertex
            valence(ev1) = valence(ev1) + 1
            valence(ev2) = valence(ev2) + 1

            !Update vconnect
            do vv=1,ubValence
                if (vconnect(ev1,vv) == 0) then 
                    vconnect(ev1,vv) = ev2 
                    exit 
                end if 
            end do 
            do vv=1,ubValence
                if (vconnect(ev2,vv) == 0) then 
                    vconnect(ev2,vv) = ev1 
                    exit 
                end if 
            end do 
        end if 
    end do 
end do 

!Set maximum valence
maxValence = maxval(valence(:))
return 
end subroutine get_vm_valence




!Subroutine to evaluate valence in an arbitrary mesh ===========================
subroutine get_arb_valence(valence,maxValence,faces,nface,nvtx)
implicit none 

!Variables - Import
integer(in) :: maxValence,nface,nvtx
integer(in), dimension(:), allocatable :: valence
type(face_data), dimension(:) :: faces

!Variables - Local
integer(in) :: ff,ee,vv
integer(in) :: ev1,ev2,ubValence,evalid
integer(in), dimension(:,:), allocatable :: vconnect

!Initialse valence array
if (allocated(valence)) then 
    deallocate(valence)
end if 
allocate(valence(nvtx))
valence(:) = 0 

!Upper bound of maximum valence 
valence(:) = 0 
do ff=1,nface
    do ee=1,faces(ff)%nvtx

        !Edge end vertices
        ev1 = ee
        ev2 = mod(ee,faces(ff)%nvtx) + 1
        ev1 = faces(ff)%vertices(ev1) 
        ev2 = faces(ff)%vertices(ev2) 

        !Accumulate valence 
        valence(ev1) = valence(ev1) + 1
        valence(ev2) = valence(ev2) + 1
    end do 
end do 
ubValence = 2*maxval(valence)

!Construct actual valence of each vertex
allocate(vconnect(nvtx,ubValence))
vconnect(:,:) = 0 
valence(:) = 0
do ff=1,nface
    do ee=1,faces(ff)%nvtx

        !Edge end vertices
        ev1 = ee
        ev2 = mod(ee,faces(ff)%nvtx) + 1
        ev1 = faces(ff)%vertices(ev1) 
        ev2 = faces(ff)%vertices(ev2)  

        !Check against vconnect 
        evalid = 1
        do vv=1,ubValence
            if (vconnect(ev1,vv) == ev2) then 
                evalid = 0
                exit
            end if 
            if (vconnect(ev2,vv) == ev1) then 
                evalid = 0
                exit 
            end if 
        end do 

        !Add valence if new edge
        if (evalid == 1) then 
            
            !Increment valence on each vertex
            valence(ev1) = valence(ev1) + 1
            valence(ev2) = valence(ev2) + 1

            !Update vconnect
            do vv=1,ubValence
                if (vconnect(ev1,vv) == 0) then 
                    vconnect(ev1,vv) = ev2 
                    exit 
                end if 
            end do 
            do vv=1,ubValence
                if (vconnect(ev2,vv) == 0) then 
                    vconnect(ev2,vv) = ev1 
                    exit 
                end if 
            end do 
        end if 
    end do 
end do 

!Set maximum valence
maxValence = maxval(valence(:))
return 
end subroutine get_arb_valence




!Subroutine to construct edges from triangular faces in a mesh ===========================
subroutine construct_edges(Nedge,edges,valence,Nvtx,Nface,faces)
implicit none 

!Variables - Import
integer(in) :: Nedge,Nvtx,Nface
integer(in), dimension(:) :: valence
integer(in), dimension(:,:) :: faces
integer(in), dimension(:,:),allocatable :: edges

!Variables - Local 
integer(in) :: ff,ee,vv 
integer(in) :: maxvalence,ev1,ev2,evalid
integer(in), dimension(:,:), allocatable :: vconnect,edge_index

!Set maxvalence
maxvalence = maxval(valence)

!Allocate vconnect and edge_index
allocate(vconnect(Nvtx,maxvalence))
allocate(edge_index(Nvtx,maxvalence))
vconnect(:,:) = 0 
edge_index(:,:) = 0 

!Build vconnect 
Nedge = 0 
do ff=1,Nface
    do ee=1,3

        !Edge end vertices
        ev1 = ee
        ev2 = mod(ee,3) + 1
        ev1 = faces(ff,ev1) 
        ev2 = faces(ff,ev2) 

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
            Nedge = Nedge + 1

            !Add edge to connection structure     
            do vv=1,maxvalence
                if (vconnect(ev1,vv) == 0) then 
                    vconnect(ev1,vv) = ev2 
                    edge_index(ev1,vv) = Nedge
                    exit 
                end if 
            end do 
            do vv=1,maxvalence
                if (vconnect(ev2,vv) == 0) then 
                    vconnect(ev2,vv) = ev1 
                    edge_index(ev2,vv) = Nedge
                    exit 
                end if 
            end do 
        end if 
    end do 
end do 

!Allocate edge array 
allocate(edges(Nedge,2))
edges(:,:) = 0 
do ff=1,Nface
    do ee=1,3

        !Edge end vertices
        ev1 = ee
        ev2 = mod(ee,3) + 1
        ev1 = faces(ff,ev1) 
        ev2 = faces(ff,ev2) 

        !Find connection and build edge 
        do vv=1,maxvalence
            if (vconnect(ev1,vv) == ev2) then 
                edges(edge_index(ev1,vv),1) = ev1 
                edges(edge_index(ev1,vv),2) = ev2 
            end if
        end do 
    end do 
end do
return 
end subroutine construct_edges




!Subroutine to construct connectivity ===========================
subroutine get_connectivity(V2E,V2F,F2E,E2F,valence,maxvalence,Nedge,Nface,Nvtx,edges,faces)
implicit none 

!Variables - Import
integer(in) :: Nedge,Nface,Nvtx
integer(in), dimension(:,:) :: edges,faces
integer(in), dimension(:) :: valence
integer(in), dimension(:,:), allocatable :: V2E,V2F,F2E,E2F

!Variables - Local 
integer(in) :: ff,ee,vv
integer(in) :: ev1,ev2,edgc,vtxc,maxvalence

!Build V2F
allocate(V2F(Nvtx,maxvalence))
V2F(:,:) = 0 
do ff=1,Nface
    do ee=1,3
        vtxc = faces(ff,ee)
        do vv=1,maxvalence
            if (V2F(vtxc,vv) == 0) then 
                V2F(vtxc,vv) = ff
                exit
            end if 
        end do 
    end do 
end do 

!Build V2E
allocate(V2E(Nvtx,maxvalence))
V2E(:,:) = 0 
do ee=1,Nedge
    
    !Edge end vertices
    ev1 = edges(ee,1) 
    ev2 = edges(ee,2) 

    !Add to V2E
    do vv=1,maxvalence
        if (V2E(ev1,vv) == 0) then 
            V2E(ev1,vv) = ee
            exit
        end if 
    end do 
    do vv=1,maxvalence
        if (V2E(ev2,vv) == 0) then 
            V2E(ev2,vv) = ee
            exit
        end if 
    end do
end do 

!Build F2E (ordered)
allocate(F2E(Nface,3))
F2E(:,:) = 0 
do ff=1,Nface
    do ee=1,3

        !Edge end vertices
        ev1 = ee
        ev2 = mod(ee,3) + 1
        ev1 = faces(ff,ev1) 
        ev2 = faces(ff,ev2) 

        !Find edge joining ev1-ev2 in this face from V2E of ev1
        do vv=1,valence(ev1)
            edgc = V2E(ev1,vv)
            if ((edges(edgc,1) == ev2) .OR. (edges(edgc,2) == ev2)) then 
                F2E(ff,ee) = edgc
                exit
            end if 
        end do 
    end do 
end do 

!Build E2F
allocate(E2F(Nedge,2))
E2F(:,:) = 0 
do ff=1,Nface
    do ee=1,3
        edgc = F2E(ff,ee)
        if ((E2F(edgc,1) .NE. ff) .AND. (E2F(edgc,2) .NE. ff)) then 
            if (E2F(edgc,1) == 0) then 
                E2F(edgc,1) = ff
            elseif (E2F(edgc,2) == 0) then 
                E2F(edgc,2) = ff
            end if 
        end if 
    end do 
end do 
return 
end subroutine get_connectivity




!Build V2F subroutine ===========================
subroutine build_V2F(V2F,maxValence,faces,nface,nvtx)
implicit none 

!Variables - Import
integer(in) :: maxValence,nface,nvtx
integer(in), dimension(:,:), allocatable :: V2F
type(face_data), dimension(:) :: faces

!Variables - Local 
integer(in) :: ff,ee,vv 
integer(in) :: vtxc

!Initialise 
if (allocated(V2F)) then 
    deallocate(V2F)
end if 

!Build V2F
allocate(V2F(Nvtx,maxvalence))
V2F(:,:) = 0 
do ff=1,Nface
    do ee=1,faces(ff)%nvtx
        vtxc = faces(ff)%vertices(ee)
        do vv=1,maxvalence
            if (V2F(vtxc,vv) == 0) then 
                V2F(vtxc,vv) = ff
                exit
            end if 
        end do 
    end do 
end do 
return 
end subroutine build_V2F


end module cellmesh3d_connectivity_mod