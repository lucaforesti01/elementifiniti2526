


struct TriQuad
    name::String
    order::Integer
    points::Matrix
    weights::Array
end


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


∇shapef_2DLFE(Q0_ref)