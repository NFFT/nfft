function [M] = construct_nodes_radial( N )

A=N*2;
P=410/512*A;
M=A*P;
file=zeros(P*A,2);
for i=1:P
  for j=1:A,
   if mod(i,2) == 0,
     r=(j-1)/A - 1/2;
   else
     r=-(j-1)/A + 1/2; 
   end
   file((i-1)*A+j,1)=r*sin((i-1)*pi/P);
   file((i-1)*A+j,2)=r*cos((i-1)*pi/P);
  end
end

% feel free to plot the nodes by uncommenting
% plot(file(:,1),file(:,2),'.-');

save nodes.dat -ascii file
