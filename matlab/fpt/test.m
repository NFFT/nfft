t = 2;
N = 2^t;
set = fptmex('init',t,0);

%%
alpha = ((2*(0:N-1)+1)./((0:N-1)+1))';
beta = zeros(N,1);
gamma = -(0:N-1)';
fptmex('precompute',set,alpha,beta,gamma,0)

%%
a = rand(4,1);
b = fptmex('trafo',set,a,length(a),0);

%%
fptmex('finalize',set)