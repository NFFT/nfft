function Z = Z(tau)
Z = kron([1,0;0,0;0,1;0,0],speye(2^tau));