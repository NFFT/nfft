function [y] = Perm(x)

% Funktion Perm
%
% Aufruf: [y] = Perm(x)  -------  Vergessen Sie nicht, den richtigen Pfad einzustellen !!!!!
%
% ==============================================================
% Multiplikation mit der Gerade-Ungerade-Permutationsmatrix P_n: 
% ==============================================================
% 
% Eingabeparameter:   x     zu permutierender Vektor der Länge n (= Zweierpotenz)
%
% Ausgangsparameter:  y     P_n * x   
%

n  = length(x);
n1 = n/2;

u1 = x(1:2:n-1);
u2 = x(2:2:n);

y = [u1' u2']';