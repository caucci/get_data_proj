################### PROJECT SETUP ###################

# Creates a new project and remove any project data that already exists
open_project -reset "proj_contr_grid"

# Adds design source files to the current project and specify some flags
add_files contr_grid.cpp -cflags "-std=c++11 -O2 -Wall -Wno-misleading-indentation -Wno-parentheses -Wstrict-aliasing -Wno-unknown-pragmas -Wno-unused-variable -Wno-unused-label"

# Adds a test bench program and files
add_files -tb contr_grid_test.cpp
add_files -tb ../data/camera0_79x79_1.5mm_tc99m_mean
add_files -tb ../data/camera0_thresh.dat
add_files -tb ../data/camera0_79x79_1.5mm_tc99m_gains
add_files -tb ../data/ResPhantom022516-0mm_00.dat

# Sets the top-level function
set_top "contr_grid_kernel"

################### SOLUTION SETUP ###################

# Creates or, if already existing, opens a new solution
open_solution -reset "solution_U250"

# Sets the target FPGA device
set_part {xcu250-figd2104-2L-e}

# Sets the clock
create_clock -period 10 -name default

config_sdx -target xocc
config_export -vivado_optimization_level 3
config_schedule -effort high

# Compiles and runs pre-synthesis C simulation using the provided C test bench
csim_design -O -clean

# Synthesizes the Vivado HLS database for the active solution
csynth_design

# Exports and packages the synthesized design in RTL
export_design -rtl verilog -format ip_catalog -xo contr_grid.xo

# Closes the current project
close_project
