function [M] = construct_knots_radial( N,Z )

A=N*2;
P=410/512*A;
M=A*P*Z;
file=zeros(M,3);

for z=0:Z-1,
for i=1:P,
  for j=1:A,
   if mod(i,2) == 0,
     r=(j-1)/A - 1/2;
   else
     r=-(j-1)/A + 1/2; 
   end
   file(z*A*P+(i-1)*A+j,1)=r*sin((i-1)*pi/P);
   file(z*A*P+(i-1)*A+j,2)=r*cos((i-1)*pi/P);
   file(z*A*P+(i-1)*A+j,3)=z/Z-0.5;
  end
end
end

% feel free to plot the knots by uncommenting
plot(file(1:M/Z,1),file(1:M/Z,2),'.-');


save knots.dat -ascii file
