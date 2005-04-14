function [y] = Butterfly(x)

% Funktion ButterflyT
%
% Aufruf: [y] = ButterflyT(x)    ... x bitte als Spaltenvektor
%
% ================================
% inverse Butterfly-Operation
% ================================
%
% Eingabeparameter: x ... zu transformierender Vektor der Länge n (= Zweierpotenz)
%
% Ausgabeparameter: y = B_n^T * x
%

% =========
% Funktion:
% =========

n  = length(x);
n1 = n/2;

u1 = x(1:n1); u2 = x(n1+1:n);
v1 = u1 + u2;
v2 = u1 - u2;
w2 = v2(n1:-1:1);

y = [v1' w2']';