#build arguments  #-Wall -fbounds-check -ffpe-trap=underflow,zero -fopt-info-optimized=$@_opt.dat
# buildargs = -O0 -Wall -fbounds-check -ffpe-trap=underflow,zero,invalid 
# buildargs = -O0 -Wall -fbounds-check
buildargs = -O2 -fopt-info-optimized=$@_opt.dat

#build settings
buildsettings = -J obj 

#object file directory
OBJDIR = obj/

#list object file names -> 1 for each .f90 in the correct compilation order
OBJS = $(addprefix $(OBJDIR), \
		cell_mesh3D_data_module.o\
		cell_mesh3D_io_module.o\
		cell_mesh3D_AD_tree_module.o\
		cell_mesh3D_geometry_module.o\
		cell_mesh3D_surface_mesh_module.o\
		cell_mesh3D_octree_module.o\
		cell_mesh3D_mesh_build_module.o\
		cell_mesh3D_postprocess_module.o\
		cell_mesh3D_mesh_generation_module.o\
		)

#object patturn rule -> for every file in $(OBJDIR) of the form var.o make it from src/var.f90
$(OBJDIR)%.o : src/%.f90
	gfortran $(buildsettings) $(buildargs) -c $< -o $@

#main build procedure
build: cm2d_link

#linking procedure
cm2d_link: $(OBJS) $(addprefix $(OBJDIR), cell_mesh3D_main.o)
	gfortran -o cell_mesh3d $^ $(buildsettings) -I obj $(buildargs) 

#clean procedure 
clean: 
	rm obj/*.mod
	rm obj/*.o 
	rm obj/*.dat