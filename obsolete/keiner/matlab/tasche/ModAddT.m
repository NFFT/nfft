function [y] = ModAddT(x)

% Funktion ModAddT
%
% Aufruf: [y] = ModAddT(x)  -------  Vergessen Sie nicht, den richtigen Pfad einzustellen !!!!!
%
% ================================================
% Multiplikation mit der Matrix A_n^T (1) : 
% ================================================
% 
% Eingabeparameter:   x     zu transformierender Vektor der Länge n (= Zweierpotenz)
%
% Ausgangsparameter:  y     A_n (1)^T * x   
%

n  = length(x);
n1 = n/2;

h1 =  sqrt(2)*x(1);
h2 =  x(2:n1);
h3 =  x(n1+1:n-1);
h4 = -sqrt(2)*x(n);

u1 = h2 + h3;
u2 = h2 - h3;

v1 = [h1 u1']';
v2 = [u2' h4]';

v2(2:2:n1) = -v2(2:2:n1);    % (I_{n1} \oplus \Sigma_{n1})

w2 = v2(n1:-1:1);            % (I_{n1} \oplus J_{n1})

y = [v1' w2']'/sqrt(2);

