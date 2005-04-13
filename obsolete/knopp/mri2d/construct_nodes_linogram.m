function [M] = construct_nodes_linogram ( N )

l=((-N/2):(N/2-1));
k=((-N):(N-1))';

x=[reshape(k*l/N^2,2*N^2,1) kron(ones(N,1),k/2/N);
   kron(ones(N,1),k/2/N) reshape(-k*l/N^2,2*N^2,1)];
   
% feel free to plot the nodes by uncommenting
% plot(x(:,1),x(:,2),'r.')

save nodes.dat -ascii x

M=size(x,1);