function [M] = construct_knots_spiral ( N,arms )
M=N^2;
file = zeros(M,2);

A=0.5;

w=N/64*50;

hold on;

for b=0:arms-1,
  for i=1:M/arms,
    t=((i-1)/(M/arms))^(1/2);
    file(b*M/arms+i,1) = A*t*cos(2*pi*(w*t+b/arms));

    file(b*M/arms+i,2) = A*t*sin(2*pi*(w*t+b/arms));
  end
  % plot(file(b*M/arms+1:(b+1)*M/arms,1),file(b*M/arms+1:(b+1)*M/arms,2),'-');
end
hold off;

% feel free to plot the knots by uncommenting
% plot(file(:,1),file(:,2),'.');

save knots.dat -ascii file
