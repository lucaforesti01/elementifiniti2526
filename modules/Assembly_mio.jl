# Author: Ivan Bioli (https://github.com/IvanBioli)
# Inspired by code written by Jochen Hinz (https://github.com/JochenHinz) for MATH-451 @ EPFL

using Memoize
using SparseArrays

"""
    initialize_assembly!(mesh::Mesh)

Initialize the assembly process for the given mesh by computing the necessary geometric quantities.

# Arguments
- `mesh::Mesh`: The mesh object for which the assembly is initialized.
"""
function initialize_assembly!(mesh::Mesh)
    get_Bk!(mesh)
    get_detBk!(mesh)
    get_invBk!(mesh)
end

########################### GLOBAL ASSEMBLER ########################### 
"""
    assemble_global(mesh::Mesh, local_assembler!)

Assemble the global stiffness matrix and force vector for the given mesh using the provided local assembler function.

# Arguments
- `mesh::Mesh`: The mesh object.
- `local_assembler!`: A function that assembles the local stiffness matrix and force vector.

# Returns
- `K::SparseMatrixCSC`: The global stiffness matrix.
- `f::Vector`: The global force vector.
"""
function assemble_global(mesh::Mesh, local_assembler!)
    ###########################################################################
    ############################ ADD CODE HERE ################################
    ########################################################################### 
end

########################################################################
########################### LOCAL ASSEMBLERS ###########################
########################################################################

########################### POISSON PROBLEM ###########################
"""
    shapef_2DLFE(quadrule::TriQuad)

Compute the shape functions for the Poisson problem.

# Arguments
- `quadrule::TriQuad`: The quadrature rule.

# Returns
- `shapef`: The shape functions evaluated at the quadrature points.
"""
function shapef_2DLFE(quadrule::TriQuad)
    ###########################################################################
    # ....
    ########################################################################### 
    PQ = quadrule.points;
    shapef = zeros(3, size(PQ, 2));
    shapef[1, :] = 1 .- PQ[1, :] .- PQ[2, :];
    shapef[2, :] = PQ[1, :];
    shapef[3, :] = PQ[2, :];

    return shapef 
end

"""
    ∇shapef_2DLFE(quadrule::TriQuad)

Compute the gradients of the shape functions for the Poisson problem.

# Arguments
- `quadrule::TriQuad`: The quadrature rule.

# Returns
- `∇shapef`: The gradients of the shape functions evaluated at the quadrature points.
"""
function ∇shapef_2DLFE(quadrule::TriQuad)
    ###########################################################################
    # ... dipende solo dal numero dei punti di quadratura del triangolo di riferimento 
    # restituisce una matrice tridimensionale in cui ogni faccia è una copia dei gradienti 
    # delle funzioni di base (che sono costanti quindi ripetuti)
    ###########################################################################
    PQ = quadrule.points;
    face = [-1 1 0; -1 0 1];
    ∇shapef = repeat(face, 1, 1, size(PQ, 2));
    return ∇shapef
end

"""
    poisson_assemble_local!(Ke::Matrix, fe::Vector, mesh::Mesh, cell_index::Integer, f)

Assemble the local stiffness matrix and force vector for the Poisson problem.

# Arguments
- `Ke::Matrix`: The local stiffness matrix to be assembled.
- `fe::Vector`: The local force vector to be assembled.
- `mesh::Mesh`: The mesh object.
- `cell_index::Integer`: The index of the current cell.
- `f`: The source term function.

# Returns
- `Ke`: The assembled local stiffness matrix.
- `fe`: The assembled local force vector.
"""
function poisson_assemble_local!(Ke::Matrix, fe::Vector, mesh::Mesh, cell_index::Integer, f)
    ###########################################################################
    ############################ ADD CODE HERE ################################
    ########################################################################### 
end