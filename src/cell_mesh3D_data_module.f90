!cell_mesh3d data types module - (derived from cell_mesh2d data types module v4.0)
!Max Wood - mw16116@bristol.ac.uk
!Univeristy of Bristol - Department of Aerospace Engineering

!Version 2.3
!Updated 12-02-2024

!Module
module cellmesh3d_data_mod

!Integer data type 
use ISO_FORTRAN_ENV, only: in32=>int32
use ISO_FORTRAN_ENV, only: in=>int64

!Real data types 
use ISO_FORTRAN_ENV, only: dp=>real64     
use ISO_FORTRAN_ENV, only: qp=>real128     

!Set round off error bound 
real(dp), parameter :: min_precision = 2.0d0*epsilon(1.0d0)

!Surface from volume interpolation type
type surfFvolinterp
    integer(in) :: npnts_vi,npnts_ss
    integer(in), dimension(:), allocatable :: vol_pnts,surf_smooth_pnts
    real(dp), dimension(:), allocatable :: surf2volRBF,surf_smoothRBF
    real(dp), dimension(:,:), allocatable :: Ri
end type surfFvolinterp

!Shell curve data type
type shellcurvetype
    integer(in) :: nvtx
    integer(in), dimension(:), allocatable :: vertices,vtag
end type shellcurvetype

!Octree edge vertex data type
type octree_evtx
    integer(in) :: Nevtx 
    integer(in), dimension(:), allocatable :: e_vertices
end type octree_evtx

!Mesh data data type 
type meshdata 
    integer(in), dimension(:), allocatable :: vtx_vlnc,cell_nedge,bc_active
    integer(in), dimension(:,:), allocatable :: vtx_2_cell,vtx_v2v,vtx_v2e,cell2edge
    real(dp), dimension(:), allocatable :: cell_area,edgedx,edgedy,edgedxd,edgedyd
    real(dp), dimension(:,:), allocatable :: cell_midp,edge_midp,cell_vtx_W
end type meshdata

!Edge intersection data type
type edgeint 
    integer(in) :: type
    integer(in) :: nitem,nint,refend1
    integer(in), dimension(:), allocatable :: vtx_idx,inttype,intidx,surfseg,int_seg_type
    integer(in), dimension(:), allocatable :: vncell,nfint
    integer(in), dimension(:,:), allocatable :: edge_mesh,vcell,vmface
    real(dp), dimension(:), allocatable :: intfrac
    real(dp), dimension(:,:), allocatable :: intloc,intlocBC 
    character(len=2), dimension(:), allocatable :: int_tri_loc,int_inout
end type edgeint 

!Edge intersection data type
type triint 
    integer(in) :: nitem,nedge,nvtx 
    integer(in), dimension(:), allocatable :: vtx_idx,edge_int_idx,nfint
    integer(in), dimension(:,:), allocatable :: edges,face_int_idx
    real(dp), dimension(:,:), allocatable :: intloc
end type triint 

!Surface intersection type
type surfint 
    integer(in) :: nint,refend1
    integer(in), dimension(:), allocatable :: vtx_idx,edg_idx
    real(dp), dimension(:), allocatable :: int_frac
end type surfint 

!Vertex intersect type
type vtx_intersect
    integer(in) :: nintersect
    integer(in) :: face_int_idx(8)
end type vtx_intersect

!Surface mesh edge potential volume mesh intersection list type 
type smvmintlist
    integer(in) :: nface,nentry 
    integer(in), dimension(:), allocatable :: faces,fint_vm_edge,fint_intersect_invalid
    real(dp), dimension(:,:), allocatable :: fint_int_loc
    character(len=2), dimension(:), allocatable :: fint_sm_loc,fint_vm_loc
end type smvmintlist

!Clipped face data type
type clipped_face_data
    character(len=3) :: ftype
    integer(in) :: nvtx,cleft,cright,fparent,fidx_final
    integer(in), dimension(:), allocatable :: vertices
end type clipped_face_data

!Face data type
type face_data
    character(len=3) :: ftype
    integer(in) :: nvtx,cleft,cright,nfclipped,cleft_ot,cright_ot,fot,fparent
    integer(in), dimension(:), allocatable :: vertices,edges
    type(clipped_face_data), dimension(:), allocatable :: face_clipped
end type face_data

!Local surface data type
type lsurface
    integer(in) :: nvtx,nface,nedge,neshell,nshellcurve,nface_simp,nvtx_simp 
    integer(in), dimension(:), allocatable :: eshell,valence,vtx_map,vtx_map2simp
    integer(in), dimension(:,:), allocatable :: faces,edges,V2E,F2E,E2F
    real(dp), dimension(:,:), allocatable :: vertices,vertices_simp 
    type(shellcurvetype), dimension(:), allocatable :: shellcurve
    type(face_data), dimension(:), allocatable :: faces_simp
end type lsurface

!Options data type
type cm3d_options
    character(len=:), allocatable :: iopath,optpath,surfacename,surface_dir,boundary_dir,meshinout,meshfrmat,mode,glink_type
    integer(in) :: Nrefine,NrefineB,Ncell_max,Nrefine_flood_i,Nrefine_flood_f,Nrefine_flood_b,meshtype,surfRcurvNpts
    integer(in) :: normDconv,NintEmax,glink_con,glink_nnn,glink_nsmooth,max_int_size,surf_force_simplify,ADTminNodedivsize
    integer(in) :: ADTmax_depth,dispt,nlpflood,nlpsmooth,surface_type,set_customBCs,remFFzones,zbndsiter,set_mbounds
    integer(in) :: Nlevel,Nsmooth_Norm,Nsmooth_front,Ls_smooth_front,Nsstype,Nzone_cBC,remISzones,remNCzones,stnintmax
    integer(in) :: bc_xmin,bc_xmax,bc_ymin,bc_ymax,bc_zmin,bc_zmax
    integer(in), dimension(:), allocatable :: BC_zone_bc
    real(dp) :: ADTpadding,far_field_bound,FminArea,CminVol,surfRcurvM,RBF_relaxP,RBF_relaxD,RBF_rsup
    real(dp) :: dsolve_kd,dsolve_cfl,dsolve_resconv,zangbnd,zlenmax,zlenmin,elenpad,otr_cellpad
    real(dp) :: cell1h,cellh_gr,om_offset_x,om_offset_y,om_offset_z,intcointol,baryloctol
    real(dp), dimension(:,:), allocatable :: BC_zone_coords
    real(dp) :: mesh_xmin,mesh_xmax,mesh_ymin,mesh_ymax,mesh_zmin,mesh_zmax
end type cm3d_options

!Mesh data type 
type vol_mesh_data
    integer(in) :: nvtx,nedge,ncell,nface,nvtx_surf,nface_mesh,nface_surface,maxValence
    integer(in), dimension(:), allocatable :: cell_level,vtx_surfseg,cell_otidx,vtx_type 
    integer(in), dimension(:), allocatable :: surf_vtx,surf_vtx_seg,vmvtx_2_smvtx,valence
    integer(in), dimension(:,:), allocatable :: edge,edges,E2F,V2E 
    real(dp), dimension(:), allocatable :: surf_vtx_segfrac,vtx_sRcurv
    real(dp), dimension(:,:), allocatable :: vtx,edge_normal
    type(face_data), dimension(:), allocatable :: faces
    type(lsurface), dimension(:), allocatable :: lsurface
    type(surfFvolinterp), dimension(:), allocatable :: vsinterp
end type vol_mesh_data

!Surface data type 
type surface_data
    integer(in) :: nvtx,nfcs,nobj,nedge,nvtxf,maxValence
    integer(in), dimension(:), allocatable :: vertex_obj,vncell,valence,vtx_vmesh_cell
    integer(in), dimension(:,:), allocatable :: connectivity,connectivityM,fesharp,edges,V2F,F2E,V2E,E2F,vcell
    real(dp), dimension(:), allocatable :: face_rcurv,face_maxcurv,vtx_rcurv_full,face_area
    real(dp), dimension(:), allocatable :: vtx_rcurv,vtx_meancurv,vtx_gausscurv,vtx_k1,vtx_k2
    real(dp), dimension(:,:), allocatable :: vertices,face_ecurv,vertices_full,face_normal,vtx_normal
    type(face_data), dimension(:), allocatable :: tri_clipped
end type surface_data

!Octree data type 
type octree_data
    integer(in) :: cins,vins,maxvalence,nedge,maxDR 
    integer(in), dimension(:), allocatable :: cell_level,cell_vcmid,vtx_valence,cell_parent,cell_keep,valence
    integer(in), dimension(:,:), allocatable :: V2V,edge_index,cell2edge
    integer(in), dimension(:,:), allocatable :: cell_vcnr,cell_vemid,cell_vfmid,cell_child,cell_adjacent,cell_diagonal
    real(dp), dimension(:,:), allocatable :: vtx
    type(octree_evtx), dimension(:), allocatable :: edge_M_vtx
end type octree_data

end module cellmesh3d_data_mod