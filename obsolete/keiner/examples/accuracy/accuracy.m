function accuracy(M)
% ACCURACY

% Generate exponents.
T_MAX = ceil(log2(M));
% generate bandwidths
M = (2*ones(1,T_MAX+1)).^(0:T_MAX);
% Genrate Gauss-Legendre nodes and weights.
generateNodes(M);

% Execute test program
!./accuracy > temp.out

% Open output file.
file = fopen('temp.out','r');

% Read error data
err = fscanf(file,'%E\n')

% Create figure
figure
% Plot data.
loglog(M,err)
% Annotate plot
xlabel('M')
ylabel(texlabel('error'))

% Safe figure to file
print('-depsc','-r1200','accuracy.eps');
