%Get Fourier coefficients from plan
function f_hat = nfsoft_get_f_hat (plan)
f_hat = nfsoftmex('get_f_hat',plan);