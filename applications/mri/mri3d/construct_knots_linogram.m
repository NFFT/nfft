function [M] = construct_knots_linogram ( N,Z )

l=((-N/2):(N/2-1));
k=((-N):(N-1))';

x=[reshape(k*l/N^2,2*N^2,1) kron(ones(N,1),k/2/N);
   kron(ones(N,1),k/2/N) reshape(-k*l/N^2,2*N^2,1)];

file=[];
for z=0:Z-1, 
  file=[file; x ones(size(x,1),1)*(z/Z-0.5)];
end

M=size(file,1);

% feel free to plot the knots by uncommenting
% plot(file(1:M/Z,1),file(1:M/Z,2),'.');


save knots.dat -ascii file

