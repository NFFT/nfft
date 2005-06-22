function [M] = construct_nodes_spiral ( N,Z )
B=1;
M=N^2*Z;
file = zeros(M,3);

A=0.5;

w=N/64*50;

for z=0:Z-1,
for b=0:B-1,
for i=1:M/B/Z,
    t=(i/(M/B/Z))^(1/2);
    file(z*M/Z+b*M/B+i,1) = A*t*cos(2*pi*w*t+b);
    file(z*M/Z+b*M/B+i,2) = A*t*sin(2*pi*w*t+b);
    file(z*M/Z+b*M/B+i,3) = z/Z-0.5;
end
end
end

% feel free to plot the nodes by uncommenting
% plot(file(1:M/Z,1),file(1:M/Z,2),'.-');

save nodes.dat -ascii file
