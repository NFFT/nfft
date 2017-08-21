%Returns the index of f_hat from the Wigner D-function D^l_{m,k}
function my_ind = nfsoft_index (l, m, k)
my_ind = NFSOFT_F_HAT_SIZE(l-1) + (k+l).*(2*l+1) + (m+l) +1;