# Author: Ivan Bioli (https://github.com/IvanBioli)
# Inspired by code written by Jochen Hinz (https://github.com/JochenHinz) for MATH-451 @ EPFL

"""
    struct TriQuad

A structure representing a triangular quadrature rule.

# Fields
- `name::String`: The name of the quadrature rule.
- `order::Integer`: The order of the quadrature rule.
- `points::Matrix`: The quadrature points.
- `weights::Array`: The quadrature weights.
"""
struct TriQuad
    name::String
    order::Integer
    points::Matrix
    weights::Array
end

###########################################################################
# Definisco gli elementi della struttura "TriQuad" contenente le formule di quadratura:
# name: è il nome della regola di quadratura (ad esempio, "Q0", "Q1", "Q2");
# order: è il grado di precisione della quadratura;
# points: è una matrice 2×q che contiene le coordinate dei punti di quadratura sull’elemento di riferimento;
# weights: è un vettore di lunghezza q che contiene i pesi associati ai punti di quadratura
###########################################################################

###########################################################################
# Definisco le matrici dei punti di quadratura per le formule Q0, Q1, Q2.
###########################################################################
M0 = zeros(2,1);
M0[1,1]=1/3;
M0[2,1]=1/3;

M1 = [0 1 0; 0 0 1];

M2 = [0.5 0 0.5; 0.5 0.5 0];

###########################################################################
# Definisco i vettori dei pesi di quadratura per le formule Q0, Q1, Q2.
###########################################################################

W0 = [0.5];
W1 = (1/6)*ones(3);
W2 = (1/6)*ones(3);

Q0_ref = TriQuad("Q0", 2, M0, W0);
Q1_ref = TriQuad("Q1", 2, M1, W1);
Q2_ref = TriQuad("Q2", 3, M2, W2);

"""
    Quadrature(u, mesh::Mesh, ref_quad::TriQuad)

Perform numerical integration of a function over a mesh using a given quadrature rule.

# Arguments
- `u`: The function to be integrated.
- `mesh::Mesh`: The mesh object.
- `ref_quad::TriQuad`: The reference quadrature rule.

# Returns
- `I_approx::Float64`: The approximate integral of the function over the mesh.
"""









function Quadrature(u, mesh::Mesh, ref_quad::TriQuad)
    ###########################################################################
    # Calcola l'integrale di una funzione sul dominio triangolato dalla mesh
    ###########################################################################
    # Inizializza l'integrale = 0
    I = 0;

    # richiama i punti e pesi di quadratura del triangolo di riferimento (da TriQuad)
    PQ = ref_quad.points;
    WQ = ref_quad.weights;
    q = length(WQ);

    # richiama l'array di matrici di trasformazioni Bk e ak (da mesh) e i determinanti detBk
    Bk, ak = get_Bk!(mesh);
    detBk = get_detBk!(mesh);
    Tri = mesh.T;

    for i in eachindex(axes(Tri,2))
        T=0; # integrale sul singolo triangolo
        B = Bk[i];
        a = ak[i];
        d = abs(detBk[i]);
        for j in 1:q
            p = PQ[:,j];
            w = WQ[j];
            # Calcola la trasformazione del punto di quadratura tramite mappa affine Bx + a
            pT = B*p + a;
            T += w * u(pT);
            
        end
        I += d*T
    end

    return I


end

# Evaluation of a function
"""
    eval_u(u::Function, points_elem::Matrix, mesh::Mesh, tri_idx::Integer, quadrule::TriQuad)

Evaluate a function at given points within an element.

# Arguments
- `u::Function`: The function to be evaluated.
- `points_elem::Matrix`: The points at which to evaluate the function.
- `mesh::Mesh`: The mesh object (ignored).
- `tri_idx::Integer`: The index of the current element (ignored).
- `quadrule::TriQuad`: The quadrature rule (ignored).

# Returns
- `u_evals::Matrix`: The evaluated function values at the given points.
"""
function eval_u(u::Function, points_elem::Matrix, mesh::Mesh, tri_idx::Integer, quadrule::TriQuad)
    ###########################################################################
    # Prende una matrice di punti e restituisce una matrice che contiene la 
    # valutazione della funzione u su ognuno dei punti
    ###########################################################################
    l = size(points_elem, 2);
    u_evals = zeros(1,l);
    for i in 1:l
        p = points_elem[:, i];
        u_evals[1, i] = u(p);
    end

    return u_evals

end

"""
    eval_u(uh::Vector, points_elem::Matrix, mesh::Mesh, tri_idx::Integer, quadrule::TriQuad)

Evaluate a linear finite element solution at given quadrature points within an element.

# Arguments
- `uh::Vector`: The finite element solution vector.
- `points_elem::Matrix`: The points at which to evaluate the solution (ignored).
- `mesh::Mesh`: The mesh object.
- `tri_idx::Integer`: The index of the current element.
- `quadrule::TriQuad`: The quadrature rule.

# Returns
- `uh_evals::Matrix`: The evaluated solution values at the given points.
"""
function eval_u(uh::Vector, points_elem::Matrix, mesh::Mesh, tri_idx::Integer, quadrule::TriQuad)
    ###########################################################################
    # Prende il vettore uh che corrisponde all'approssimazione lineare della soluzione 
    # su un elemento, valuta questa funzione sui nodi di quadratura del dato elemento
    ###########################################################################


    # Punti di quadratura del metodo per il triangolo di riferimento:
    PQ = quadrule.points

    # Punti di quadratura del metodo per il triangolo generico trasformato: prende 
    # ogni punto della matrice PQ e lo trasforma secondo la trasformazione dell'elemento 
    # scelto: dipende dall'indice del triangolo tri_idx
    T = mesh.T;
    PT = PQ;
    Bk, ak = get_Bk!(mesh);
    for i in 1:size(PQ,2)
        PT[:, i]= Bk[tri_idx]*PQ[:, i] + ak[tri_idx];
    end

    # Di tutti i valori di uh estraggo quelli che mi interessano, cioè dell'elemento (loc)
    # corrispondente a tri_idx: 
    uh_loc = [uh[i] for i in T[:, tri_idx]];

    # Inizializza 
    uh_evals = zeros(1,size(PT, 2));



end

