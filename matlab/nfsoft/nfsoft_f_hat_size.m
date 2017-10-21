%Returns size of vector f_hat for degree N
function f = nfsoft_f_hat_size(N)
f = ((N+1).*(4*(N+1).*(N+1)-1)/3);
return;