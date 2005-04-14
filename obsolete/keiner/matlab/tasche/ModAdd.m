function [y] = ModAdd(x)

% Funktion ModAdd
%
% Aufruf: [y] = ModAdd(x)  -------  Vergessen Sie nicht, den richtigen Pfad einzustellen !!!!!
%
% ================================================
% Multiplikation mit der Matrix A_n (1) : 
% ================================================
% 
% Eingabeparameter:   x     zu transformierender Vektor der Länge n (= Zweierpotenz)
%
% Ausgangsparameter:  y     A_n (1) * x   
%

n  = length(x);
n1 = n/2;

u1 = x(1:n1); u2 = x(n:-1:n1+1);   % (I_n1 \oplus J_n1)

u2(2:2:n1) = -u2(2:2:n1);          % (I_n1 \oplus \Sigma_n1)

v1 = u1(2:n1); v2 = u2(1:n1-1);
w1 = v1 + v2;
w2 = v1 - v2;

z1 =  sqrt(2)*u1(1);
z2 = -sqrt(2)*u2(n1);

y = [z1 w1' w2' z2]'/sqrt(2);

