function [y] = Twiddle(x)

% Funktion Twiddle
%
% Aufruf: [y] = Twiddle(x)  -------  Vergessen Sie nicht, den richtigen Pfad einzustellen !!!!!
%
% ================================================
% Multiplikation mit der Matrix T_n (1): 
% ================================================
% 
% Eingabeparameter:   x     zu transformierender Vektor der Länge n (= Zweierpotenz)
%
% Ausgangsparameter:  y     T_n (1) * x   
%

n  = length(x);
n1 = n/2;
t=ones(n1,1);
cn1=ones(n1,1);
sn1=ones(n1,1);
k=1;
for k=1:n1,
   t(k) = (2*k-1)*pi/(8*n1);
end;

cn1 = sqrt(2)*cos(t);
sn1 = sqrt(2)*sin(t);

u1 = x(1:n1); u2 = x(n:-1:n1+1);
v1 =  cn1.*u1 + sn1.*u2;

v2 = -sn1.*u1 + cn1.*u2;
w2 = v2(n1:-1:1);
w2(2:2:n1) = - w2(2:2:n1);           % \Sigma_{n1}

y = [v1' w2']';


