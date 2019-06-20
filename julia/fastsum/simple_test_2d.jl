push!(LOAD_PATH, pwd())
using fastsum
using LinearAlgebra

println("2d fastsum test")

# set the parameters:
N = 100;
M = 100;
kernel = "multiquadric";
c = (1/sqrt(N),3.5);
m = 4;
p = 7;
n = 100*128;
eps_I = p/n;
eps_B =  1/32;

# create a Plan-Object in Julia
pt = fastsumplan(N,M,n,m,p,kernel,c,eps_I,eps_B);

# generate source nodes in circle of radius 0.25-eps_B/2
r = sqrt.(rand(N))*(0.25-eps_B/2);
phi = rand(N)*2*pi;
xhat = zeros(2,N);
xhat[1,:]=r.*cos.(phi);
xhat[2,:]=r.*sin.(phi);
pt.x = xhat;

# generate coefficients alpha_k
alpha_temp = [1/(1+k) for k in 0:N-1] +im*[1/(1+k^2) for k in 0:N-1];

pt.alpha = alpha_temp;

# generate target nodes in circle of radius 0.25-eps_B/2
r = sqrt.(rand(M))*(0.25-eps_B/2);
phi = rand(M)*2*pi;
y_temp=zeros(2,M);
y_temp[1,:] = r.*cos.(phi)
y_temp[2,:] = r.*sin.(phi);
pt.y = y_temp;
# Start the Transformation
fastsum.trafo(pt);

f1 = pt.f;
f_alg = zeros(N)+im*zeros(N);
for i=1:M
    f_alg[i] = f1[i];
end

fastsum.trafo_exact(pt);
f2 = pt.f_exact;
error_vector = f_alg-f2;

E_2 = norm(error_vector)/norm(f1)
E_infty = norm(error_vector, Inf)/norm(pt.alpha,1)
println( "E_2 error:" )
println( E_2 )
println( "E_infty error:" )
println( E_infty );
