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
    # restituisce i valori delle funzioni di base valutate sui punti di quadratura 
    # di un metodo scelto 
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
    # restituisce i valori dei gradienti delle funzioni di base valutate sui punti 
    # di quadratura di un metodo scelto.
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
    ###########################################################################
    # richiama i punti della mesh nella matrice T
    T = mesh.T;
    p = mesh.p;
    Bk, ak = get_Bk!(mesh);
    Bk_inv = get_invBk!(mesh);
    detBk = get_detBk!(mesh);

    # TRIANGOLO LOCALE
    a = ak[cell_index];
    B = Bk[cell_index];
    B_inv = Bk_inv[cell_index]; # inversa della matrice del triangolo dato in input
    detB = detBk[cell_index]; # determinante della matrice del triangolo dato in input

    # Seleziona le etichette dei vertici del triangolo dato in input in un vettore N
    # N = T[:, cell_index];

    #inizializza gli output
    Ke = zeros(3,3);
    fe = zeros(3);

    # richiama i gradienti delle funzioni di base sul baricentro del triangolo di riferimento,
    # per usarle nella quadratura per Ke
    grad_Q0 = ∇shapef_2DLFE(Q0_ref);

    # richiama i valori delle funzioni di base sui punti di quadratura di Q2_ref per calcolare fe
    base_Q2 = shapef_2DLFE(Q2_ref);
    for i =1:3
        
        Q2_W = Q2_ref.weights;
        Q2_P = Q2_ref.points;
        # fe[i] = (Q2_W)'*( f.( B.* Q2_P .+ a) .* base_Q2[i,:] )



        for j = 1:3
            # seleziona il gradiente dell'i-esima e j-esima funzione di base nel baricentro del triangolo di riferimento
            ∇Q0_i = grad_Q0[:, i , :];
            ∇Q0_j = grad_Q0[:, j , :];

            # formula di quadratura punto medio (non dipende dal punto medio perché 
            # i gradienti delle funzioni di base sul triangolo di riferimento sono costanti)
            # 1/2 = area triangolo di riferimento 
            Ke[i,j] = 1/2 * abs(detB) * (    (B_inv' * ∇Q0_i)'*(B_inv' * ∇Q0_j)   );
        end
    end



end