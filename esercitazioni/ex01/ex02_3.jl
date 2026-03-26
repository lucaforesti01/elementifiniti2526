# Author: Ivan Bioli (https://github.com/IvanBioli)

import Pkg
Pkg.activate("elementifinitiunipv_pkg")
using Revise
includet("../../modules/Meshing.jl")
using Plots

function add_triangle(p1, p2, p3)
    x = [p1[1], p2[1], p3[1], p1[1]]  # x-coordinates (closed loop)
    y = [p1[2], p2[2], p3[2], p1[2]]  # y-coordinates (closed loop)
    plot!(x, y, label="", linewidth=1, color=:black)
end

# Function to plot the mesh
function plot_mesh(T, p)
    fig = plot() # Initialize an empty plot
    
    ###################################
    ########## PLOT THE MESH ##########
    ####################################

    return fig 
end

# Create the mesh
h = 0.05
out_file = mesh_square(h)
T, p = get_nodes_connectivity(out_file)

# Plot the mesh
fig = plot_mesh(T, p)
plot!(fig, aspect_ratio=1)
savefig("figures/ex02_1.pdf")
display(fig)

# Compare with Meshes.jl
import Meshes
mesh = to_Meshes(T, p)
Meshes.viz(mesh, showsegments = true)
