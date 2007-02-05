function [M] = construct_knots_rose( N )

N=ceil(1.5*N); 
B=1;
M=N^2;
file = zeros(M,2);

A=0.5;

w=N/64*50;

for b=0:B-1,
for i=1:M/B,
    t=i/(M/B);
    file(b*M/B+i,1) = A*cos(2*pi*w*t)*cos(2*pi*t+b*2*pi/B);

    file(b*M/B+i,2) = A*cos(2*pi*w*t)*sin(2*pi*t+b*2*pi/B);
end
end

% feel free to plot the knots by uncommenting
% plot(file(:,1),file(:,2),'.-');

save knots.dat -ascii file
