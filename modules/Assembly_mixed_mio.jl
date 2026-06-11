# Author: Ivan Bioli (https://github.com/IvanBioli)
# Inspired by code written by Jochen Hinz (https://github.com/JochenHinz) for MATH-451 @ EPFL

using Memoize
using SparseArrays

"""
    initialize_assembly_mixed!(mesh::Mesh)

Initializes the assembly process for a mixed finite element method on the given `mesh`.

This function performs the following steps:
- Computes the element transformation matrices by calling `get_Bk!`.
- Calculates the determinants of the transformation matrices by calling `get_detBk!`.

# Arguments
- `mesh::Mesh`: The mesh data structure on which the assembly is to be initialized.

# Side Effects
Modifies the internal state of `mesh` by updating its transformation matrices and their determinants.

# See Also
- [`get_Bk!`](@ref)
- [`get_detBk!`](@ref)
"""
function initialize_assembly_mixed!(mesh::Mesh)
    get_Bk!(mesh)
    get_detBk!(mesh)
end

########################### GLOBAL ASSEMBLER ########################### 
"""
    assemble_global_mixed(mesh::Mesh, local_assembler!)

Assembles the global mixed finite element system matrices and right-hand side vector for a given mesh.

# Arguments
- `mesh::Mesh`: The mesh data structure containing information about elements, edges, and connectivity.
- `local_assembler!`: A function that computes the local element matrices and force vector for a given cell. It should have the signature `local_assembler!(Ae, Be, fe, mesh, cell_index)`.

# Returns
- `K`: The global system matrix assembled as a block matrix, where the upper-left block corresponds to the velocity matrix, the upper-right and lower-left blocks correspond to the coupling (pressure) matrices, and the lower-right block is a zero matrix.
- `b`: The global right-hand side vector, concatenating zeros for velocity DOFs and the assembled force vector for pressure DOFs.
"""
function assemble_global_mixed(mesh::Mesh, local_assembler!)
    ###########################################################################
    # matrice [ A  tB ]
    #         [ B   0 ]  è una matrice (#(funz base totali) + N_tri) x (#(funz base totali) + N_tri)
    # f = vettori degli integrali di f sui triangoli 
    ########################################################################### 
    T = mesh.T;
    p = mesh.p;
    Ntri = size(T,2);
    Npoints = size(p,2);

    A_loc = zeros(3,3);B_loc = zeros(1,3); fe_loc = [0.0];


    A_glob = zeros(Npoints, Npoints);
    B_glob = zeros(Ntri, Npoints);
    O = zeros(Ntri, Ntri);
    b = zeros(Npoints + Ntri);

    for k = 1:Ntri
        local_assembler!(A_loc, B_loc, fe_loc, mesh, k);
        indici = T[:, k] # indici dei vertici di k 
       
        b[Npoints + k] = fe_loc[1]
        A_glob[indici, indici] .+= A_loc
        B_glob[k,indici] = B_loc


    end

    K = [A_glob transpose(B_glob); B_glob O]



    # rows_A = []
    # cols_A = []
    # data_A = Float64[]
    
    # rows_B = []
    # cols_B = []
    # data_B = Float64[]

    # row_fe = []
    # data_fe = Float64[]

    # # inizializza le matrici e vettori locali 
    # A_loc = zeros(3,3);
    # fe_loc = [0.0];
    # B_loc = zeros(1,3);

    # for k in 1:Ntri
    #     local_assembler!(A_loc, B_loc, fe_loc, mesh, k);
    #     indices = T[:, k];
    #     append!(rows_B, fill(k, 3)) # indice di colonna di B (numero di triangoli)

    #     push!(row_fe, k); # sistema gli indici per fe
    #     push!(data_fe, fe_loc[1]); # sistema i dati di fe

    #     for i in 1:3
    #         i_glob = indices[i];
    #         push!(cols_B, i_glob); # sistema gli indici di B 
    #         push!(data_B, B_loc[1, i]); # sistema i dati di B 


    #         for j in 1:3
    #             j_glob = indices[j];
    #             push!(rows_A, i_glob); # indici di riga
    #             push!(cols_A, j_glob); # indici di colonna
    #             push!(data_A, A_loc[i,j]); # dati corrispondenti sommati direttamente
    #         end

            
    #     end
    # end



    # A_glob = sparse(rows_A, cols_A, data_A, Npoints, Npoints);
    # B_glob = sparse(rows_B, cols_B, data_B, Ntri, Npoints)
    # fe_glob = Matrix(sparse(row_fe, ones(size(row_fe)), data_fe));

    # # assembla la matrice a blocchi
    # C = spzeros(Ntri, Ntri)
    # K = sparse( [A_glob transpose(B_glob); B_glob    C])
    # b = [zeros(Npoints); fe_glob]
    
    return K, b
    











end


"""
    shapef_2D_RT0FE(quadrule::TriQuad)

Compute the Raviart-Thomas RT0 vector-valued shape functions on the reference triangle at the given quadrature points.

# Arguments
- `quadrule::TriQuad`: A quadrature rule object for triangles, containing the quadrature points as a 2×n matrix.

# Returns
- `shapef::Array{Float64,3}`: A 3D array of size (2, 3, n), where:
    - The first dimension (2) corresponds to the vector components (x and y).
    - The second dimension (3) corresponds to the three RT0 basis functions.
    - The third dimension (n) corresponds to the number of quadrature points.

# Details
The RT0 basis functions on the reference triangle are:
- `f1(x, y) = [x; y - 1]`
- `f2(x, y) = [x; y]`
- `f3(x, y) = [x - 1; y]`

The function evaluates these basis functions at each quadrature point and returns them in a single array for efficient use in finite element assembly.

# Memoization
The function is memoized to cache results for repeated calls with the same quadrature rule.
"""
# FIXME: PUT MEMOIZE BACK AFTER IMPLEMENTATION
# @memoize function shapef_2D_RT0FE(quadrule::TriQuad)
function shapef_2D_RT0FE(quadrule::TriQuad)
    ###########################################################################
    # valuta le funzioni di base dei Raviart - Thomas sui punti di quadratura 
    # del triangolo di riferimento
    ########################################################################### 
    Quad_points = quadrule.points
    N_q = size(Quad_points, 2) # numero di punti di quadratura 

    # inizializza la matrice delle valutazioni
    shapef = zeros(2,3,N_q)
    for t = 1:N_q
        q = Quad_points[:, t] # t-esimo punto di quadratura  
        v1 = [q[1], q[1], q[1] - 1] # primo vettore [x, x, x-1]
        v2 = [q[2]-1, q[2], q[2]] # secondo vettore [y-1, y, y]
        # assembla la faccia t-esima della matrice con le valutazioni sul punto
        shapef[1, :, t] = v1
        shapef[2, :, t] = v2
    end 

    return shapef

end


"""
    divshapef_2D_RT0FE(quadrule::TriQuad)

Compute the divergence of the Raviart-Thomas RT0 vector-valued basis functions on a 2D triangle at the given quadrature points.

# Arguments
- `quadrule::TriQuad`: A quadrature rule object containing the quadrature points for integration over a triangle.

# Returns
- `divshapef::Array{Int,3}`: A `(1, 3, n)` array, where `n` is the number of quadrature points. Each entry contains the divergence of the three RT0 basis functions, which is constant and equal to 2 for each function.

# Notes
- The RT0 basis functions on the reference triangle are:
    - `f₁(x, y) = [x; y - 1]`
    - `f₂(x, y) = [x; y]`
    - `f₃(x, y) = [x - 1; y]`
- The divergence of each basis function is constant and equal to 2.
- The result is repeated for each quadrature point.
"""
# FIXME: PUT MEMOIZE BACK AFTER IMPLEMENTATION
# @memoize function divshapef_2D_RT0FE(quadrule::TriQuad)
function divshapef_2D_RT0FE(quadrule::TriQuad)
    ###########################################################################
    # valutazione dei gradienti delle funzioni di base RT0 sui punti di quadratura 
    # del triangolo di riferimento. oss: div(φ1) = 2; div(φ2) = 2; div(φ3) = 2
    ########################################################################### 
    N_q = size(quadrule.points,2)
    v = [2, 2, 2]
    divshapef = repeat(v', 1, 1, N_q) # il trasposto su v serve perchè se no avrei una matrice 3x1xN_q (v colonna)
    return divshapef
end

"""
    darcy_assemble_local_mixed!(Ae::Matrix, Be::Matrix, fe::Vector, mesh::Mesh, cell_index::Integer, f, μ)

Assembles the local matrices and vector for the mixed finite element formulation of the Darcy problem on a single cell.

# Arguments
- `Ae::Matrix`: Local stiffness matrix (to be filled in-place).
- `Be::Matrix`: Local divergence matrix (to be filled in-place).
- `fe::Vector`: Local right-hand side vector (to be filled in-place).
- `mesh::Mesh`: Mesh data structure containing geometry and connectivity.
- `cell_index::Integer`: Index of the current cell in the mesh.
- `f`: Function representing the source term, evaluated as `f(x)` at a point `x`.
- `μ`: Function representing the permeability coefficient, evaluated as `μ(x)` at a point `x`.

# Description
This function computes the local contributions to the global system for the mixed finite element discretization of the Darcy problem. It uses RT0 (Raviart-Thomas of lowest order) basis functions and Piola transformation for mapping reference shape functions to the physical element. The function performs numerical integration using a quadrature rule, and assembles the local matrices and vector by looping over quadrature points and basis functions.

# Notes
- The function assumes that the shape functions are oriented consistently using `elems2orientation`.
- The local matrices and vector are reset to zero at the beginning of the function.
- The function modifies `Ae`, `Be`, and `fe` in-place and returns them for convenience.

# Returns
- `(Ae, Be, fe)`: The assembled local stiffness matrix, divergence matrix, and right-hand side vector.
"""
########################### DARCY PROBLEM ###########################
function darcy_assemble_local_mixed!(Ae::Matrix, Be::Matrix, fe::Vector, mesh::Mesh, cell_index::Integer, f, μ)
    ###########################################################################
    # Assemblaggio locale delle matrici A e B e del load vector f
    # A_loc: matrice 3x3 (3 = numero di funzioni base sul triangolo di riferimento)
    # B_loc: matrice 1 x 3 (per ogni triangolo le tre divergenze delle funzioni di base che ho sul triangolo) 
    # nella matrice B avrò su ogni triangolo 3 entrate non nulle (per il triangolo fissato) e tutte le altre nulle 
    # alla fine ho 
    # matrice [ A  tB ]
    #         [ B   0 ]  è una matrice (#(funz base totali)=#(vertici totali mesh) + N_tri) x (#(funz base totali) + N_tri)
    # fe_loc: integrale sul triangolo fissato della funzione f (per me è uno scalare)
    ########################################################################### 
    # inizializza a zero gli output, stiffness locale e load vector
    fill!(Ae, 0.0)
    fill!(Be, 0.0)
    fill!(fe, 0.0)

    # richiama la formula di quadratura da utilizzare
    quadrule = Q2_ref

    # richiama le trasformazioni dal triangolo di riferimento al triangolo della mesh
    Bk, ak = get_Bk!(mesh);
    detBk = get_detBk!(mesh);

    # sul triangolo della mesh
    B = Bk[:, :, cell_index];
    a = ak[:, cell_index];
    detB = detBk[cell_index];


    # trasforma i punti di quadratura sul triangolo selezionato
    pe = B * quadrule.points .+ a;

    # valutazioni delle funzioni di base RT e divergenze sui punti di quadratura di riferimento
    shapef = shapef_2D_RT0FE(quadrule) # 2 x 3 x N_quadpoints
    divshapef = divshapef_2D_RT0FE(quadrule) # 1 x 2 x N_quadpoints

    # ciclo sui punti di quadratura 
    for p in eachindex(axes(pe,2))
        # peso di quadratura riscalato moltiplicando per il valore assoluto del determinante della matrice di trasformazione
        dΩ = quadrule.weights[p] * abs(detB);

        # aggiorna il load vector: è un vettore di N_tri componenti che sulla componente k = cell_index 
        # restituisce l'integrale della funzione f sul triangolo k. Quindi in teoria qui è solo una componente 
        fe .-= (f(pe[:, p]) * dΩ);

        # ciclo sulle funzioni di base dgli RT0
        for i = 1:3
            # calcola il segno della normale sull'i-esimo lato del triangolo cell_index
            sign_i = mesh.elems2orientation[i, cell_index] 
            # restituisce le valutazioni delle funzioni di base sui punti di quadratura trasformati con Piola
            v = (sign_i/detB) * B * shapef[:, i, p]
            # restituisce le divergenze valutate sui punti di quadratura trasformati con Piola
            div_v = (sign_i/detB) * divshapef[1, i, p]

            # assembla la riga k-esima della matrice B: la matrice B è tale che sulla riga k l'entrata della colonna i contiene l'integrale
            # sul triangolo k (dovuto all'indicatrice) della divergenza della funzione di base φi quindi qui, siccome è fissato k = cell_index
            # costruisco al variare di i = 1...#(funzioni di base)=3 gli integrali in riga. 
            # Insomma assemblo la k-esima riga di quella che sarà la matrice B, dove k = cell_index. 
            Be[1, i] -= div_v * dΩ  

            # NB: i varia sulle colonne (così ho le tre funzioni di base sulle colonne)
            # il ciclo va su tutti i punti di quadratura p per calcolare l'integrale  

            for j = 1:3
                # calcola il segno della normale sul j-esimo lato del triangolo cell_index
                sign_j = mesh.elems2orientation[j, cell_index] 
                # restituisce la valutazione della funzione di base sui punti di quadratura trasformati con Piola
                u = (sign_j/detB) * B * shapef[:, j, p]
                Ae[i,j] += μ(pe[:, p]) * (v ⋅ u) *dΩ
            end 



        end 



    end 

    return Ae, Be, fe

end


########################## DEFINE FUNCTIONS TO COMPUTE ERROR ##########################
"""
    L2error_mixed_p(p::Function, ph::Vector, mesh::Mesh, ref_quad::TriQuad) -> Float64

Compute the L² error between an exact function `p` and its discrete approximation `ph` over a given mesh using quadrature.

# Arguments
- `p::Function`: The exact function to be evaluated.
- `ph::Vector`: Vector of discrete approximations (one per element).
- `mesh::Mesh`: The mesh structure containing element connectivity and geometry.
- `ref_quad::TriQuad`: Quadrature rule on the reference triangle, providing points and weights.

# Returns
- `Float64`: The computed L² error over the mesh.
"""
function L2error_mixed_p(p::Function, ph::Vector, mesh::Mesh, ref_quad::TriQuad)
    ###########################################################################
    ############################ ADD CODE HERE ################################
    ########################################################################### 
end

"""
    H1diverror_mixed_u(u::Function, divu::Function, uh::Vector, mesh::Mesh, ref_quad::TriQuad) -> Float64

Compute the H(div) error norm between the exact solution `(u, divu)` and the finite element solution `uh` 
for a mixed finite element method using Raviart-Thomas (RT0) elements on a triangular mesh.

# Arguments
- `u::Function`: The exact vector-valued solution function, accepting a point and returning a 2D vector.
- `divu::Function`: The exact divergence of the solution, accepting a point and returning a scalar.
- `uh::Vector`: The vector of degrees of freedom for the finite element solution (RT0 coefficients).
- `mesh::Mesh`: The mesh data structure, containing element connectivity, geometry, and orientation.
- `ref_quad::TriQuad`: Quadrature rule on the reference triangle, providing points and weights.

# Returns
- `Float64`: The combined H(div) error norm, i.e., `sqrt(∫|u - uh|^2) + sqrt(∫|divu - div(uh)|^2)` over the domain.
"""
function H1diverror_mixed_u(u::Function, divu::Function, uh::Vector, mesh::Mesh, ref_quad::TriQuad)
    ###########################################################################
    ############################ ADD CODE HERE ################################
    ########################################################################### 
end