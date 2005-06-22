function [M] = construct_knots_spiral ( N )
B=1;
M=N^2;
file = zeros(M,2);

A=0.5;

w=N/64*50;

for b=0:B-1,
for i=1:M/B,
    t=(i/(M/B))^(1/2);
    file(b*M/B+i,1) = A*t*cos(2*pi*w*t+b);

    file(b*M/B+i,2) = A*t*sin(2*pi*w*t+b);
end
end

% feel free to plot the knots by uncommenting
% plot(file(:,1),file(:,2),'.-');

save knots.dat -ascii file
