function [y] = Butterfly(x)

% Funktion Butterfly
%
% Aufruf: [y] = Butterfly(x)    ... x bitte als Spaltenvektor
%
% ================================
% Butterfly-Operation
% ================================
%
% Eingabeparameter: x ... zu transformierender Vektor der Länge n (= Zweierpotenz)
%
% Ausgabeparameter: y = B_n * x
%

% =========
% Funktion:
% =========

n  = length(x);
n1 = n/2;

u1 = x(1:n1); u2 = x(n:-1:n1+1);
v1 = u1 + u2;
v2 = u1 - u2;

y = [v1' v2']';