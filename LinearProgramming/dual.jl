include("LP_matrix.jl")

function dual(LP::LP_matrix)
  A_transpose = - (LP.A)';
  b           = - (LP.c);
  c           = - (LP.b);
end
