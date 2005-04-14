function [y] = PermT(x)

% Funktion PermT
%
% Aufruf: [y] = PermT(x)  -------  Vergessen Sie nicht, den richtigen Pfad einzustellen !!!!!
%
% ==============================================================
% Multiplikation mit der (n/2)-Permutationsmatrix P_n^T: 
% ==============================================================
% 
% Eingabeparameter:   x     zu permutierender Vektor der Länge n (= Zweierpotenz)
%
% Ausgangsparameter:  y     P_n^T * x   
%

n  = length(x);
n1 = n/2;
i=1;
u = ones(n,1);
for i = 1:n1, u(2*i-1) = x(i); u(2*i) = x(i+n1);
end

y = u;

