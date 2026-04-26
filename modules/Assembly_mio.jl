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
    ########################################################################### 
    T = mesh.T;
    p = mesh.p;
    Ntri = size(T,2);
    Npoints = size(p,2);

    rows = []
    cols = []
    data = Float64[]
    
    row_fe = []
    data_fe = Float64[]

    A_loc = zeros(3,3);
    fe_loc = zeros(3);
    for k in 1:Ntri
        local_assembler!(A_loc, fe_loc, mesh, k);
        indices = T[:, k];
        for i in 1:3
            i_glob = indices[i];
            for j in 1:3
                j_glob = indices[j];
                push!(rows, i_glob); # indici di riga
                push!(cols, j_glob); # indici di colonna
                push!(data, A_loc[i,j]); # dati corrispondenti sommati direttamente
            end

            push!(row_fe, i_glob);
            push!(data_fe, fe_loc[i]);
        end
    end
    A_glob = sparse(rows, cols, data, Npoints, Npoints);
    fe_glob = Matrix(sparse(row_fe, ones(size(row_fe)), data_fe));

    return A_glob, fe_glob 

end

"""
    impose_dirichlet(A, b, g, mesh)

Impose Dirichlet boundary conditions on the system.

# Arguments
- `A`: The global stiffness matrix.
- `b`: The global force vector.
- `g`: The Dirichlet boundary condition function.
- `mesh::Mesh`: The mesh object.

# Returns
- `A_cond`: The modified stiffness matrix with Dirichlet conditions imposed.
- `b_cond`: The modified force vector with Dirichlet conditions imposed.
- `uh`: The solution vector with Dirichlet conditions applied.
"""
function impose_dirichlet(A, b, g, mesh)
    ###########################################################################
    F = mesh.freedofs;
    D = mesh.dirichletdofs;
    T = mesh.T
    p = mesh.p
    pd = p[:,T[D]]; # punti di bordo (coordinate)

    A_cond = A[F,F];  # sottomatrice di stiffness corrispondente ai punti interni alla mesh
    b_cond = b[F] - A[F,D] * g.(eachcol(pd)); # load vector corrispondente ai punti interni alla mesh

    u_f = A_cond\b_cond;

    # Definisce il vettore uh inserendo negli indici F il vettore u_f e negli indici D la funzione g calcolata sui punti del bordo di Dirichlet
    uh = zeros(size(T,2));
    uh[F] = u_f
    uh[D]= g.(eachcol(pd))

    return A_cond, b_cond, uh
    

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
    shapef[1, :] = 1 .- PQ[1, :] .- PQ[2, :]; # φ1 = 1 - x - y
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
    a = ak[:, cell_index];
    B = Bk[:, :, cell_index];
    B_inv = Bk_inv[:, :, cell_index]; # inversa della matrice del triangolo dato in input
    detB = detBk[cell_index]; # determinante della matrice del triangolo dato in input

    # Seleziona le etichette dei vertici del triangolo dato in input in un vettore N
    # N = T[:, cell_index];

    #inizializza gli output
    fill!(Ke, 0)
    fill!(fe, 0)


    # richiama i gradienti delle funzioni di base sul baricentro del triangolo di riferimento,
    # per usarle nella quadratura per Ke
    grad_Q0 = ∇shapef_2DLFE(Q0_ref);

    # richiama i valori delle funzioni di base sui punti di quadratura di Q2_ref per calcolare fe
    base_Q2 = shapef_2DLFE(Q2_ref);
    W2 = Q2_ref.weights # pesi i quadratura
    P2 = Q2_ref.points # punti di quadratura

    for i =1:3
        # calcola la Fe[i] sommando i valori moltiplicati per i pesi di
        for k in eachindex(W2)
            fe[i]+= W2[k] * f(B * P2[:, k] + a) * base_Q2[i, k] * abs(detB);
        end 

        for j = 1:3
            # seleziona il gradiente dell'i-esima e j-esima funzione di base nel baricentro del triangolo di riferimento
            ∇Q0_i = grad_Q0[:, i , :];
            ∇Q0_j = grad_Q0[:, j , :];

            # formula di quadratura punto medio (non dipende dal punto medio perché 
            # i gradienti delle funzioni di base sul triangolo di riferimento sono costanti)
            # 1/2 = area triangolo di riferimento 
            Ke[i,j] = 1/2 * abs(detB) * dot(    (B_inv' * ∇Q0_i), (B_inv' * ∇Q0_j)   );
        end
    end



end


########################### TRANSPORT PROBLEM ###########################
"""
    transport_assemble_local!(Ke::Matrix, fe::Vector, mesh::Mesh, cell_index::Integer, f, k, β; stab = nothing, δ = 0.5)

Assemble the local stiffness matrix and force vector for the transport problem.

# Arguments
- `Ke::Matrix`: The local stiffness matrix to be assembled.
- `fe::Vector`: The local force vector to be assembled.
- `mesh::Mesh`: The mesh object.
- `cell_index::Integer`: The index of the current cell.
- `f`: The source term function.
- `k`: The diffusion coefficient function.
- `β`: The advection velocity function.
- `stab`: The stabilization method (optional).
- `δ`: The stabilization parameter (optional).

# Returns
- `Ke`: The assembled local stiffness matrix.
- `fe`: The assembled local force vector.
"""
function transport_assemble_local!(Ke::Matrix, fe::Vector, mesh::Mesh, cell_index::Integer, f, k, β; stab = nothing, δ = 0.5)
    ###########################################################################
    ############################ ADD CODE HERE ################################
    ########################################################################### 
end

########################### DARCY PROBLEM ###########################
"""
    darcy_assemble_local!(Ke::Matrix, fe::Vector, mesh::Mesh, cell_index::Integer, f, k)

Assemble the local stiffness matrix and force vector for the Darcy problem.

# Arguments
- `Ke::Matrix`: The local stiffness matrix to be assembled.
- `fe::Vector`: The local force vector to be assembled.
- `mesh::Mesh`: The mesh object.
- `cell_index::Integer`: The index of the current cell.
- `f`: The source term function.
- `k`: The permeability coefficient function.

# Returns
- `Ke`: The assembled local stiffness matrix.
- `fe`: The assembled local force vector.
"""
function darcy_assemble_local!(Ke::Matrix, fe::Vector, mesh::Mesh, cell_index::Integer, f, k)
    ###########################################################################
    ############################ ADD CODE HERE ################################
    ########################################################################### 
end
