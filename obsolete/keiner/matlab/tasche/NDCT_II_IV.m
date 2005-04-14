function [y] = NDCT_II_IV(x,k)

% Funktion NDCT_II_IV
%
% Aufruf: [y] = NDCT_II_IV(x,k)  -------  Vergessen Sie nicht, den richtigen Pfad einzustellen !!!!!
%
% ================================================
% Rekursiver DCT-II(n) bzw. DCT-IV(n) Algorithmus : 
% ================================================
% 
% Eingabeparameter:   x     zu transformierender Vektor
%                     k     Typ der Transformation (2 oder 4)
%
% Ausgangsparameter:  y     sqrt{n} * C_n^{II} * x     im Fall k=0
%                           sqrt(n) * C_n^{IV} * x     im Fall k=1
%


% ========================
% Definition von C_2^{IV}:
% ========================

qw2 = sqrt(2);
el1 = qw2 * cos(pi/8);
el2 = qw2 * sin(pi/8);
C4 = [el1 el2; el2 -el1];


% ==================
% Hauptprogramm :
% ==================

n  = length(x);
n1 = n/2;

switch k, 
  case 0, switch n1,
     case 1, [y] = Butterfly(x);
     otherwise, u = Butterfly(x); u1 = u(1:n1); u2 = u(n1+1:n); 
        [v1] = NDCT_II_IV(u1,0);
        [v2] = NDCT_II_IV(u2,1);
          v  = [v1',v2']';
          y  = PermT(v);
     end;      
  case 1, switch n1,
     case 1, y = C4 * x;
     otherwise, u = Twiddle(x); u1 = u(1:n1); u2 = u(n1+1:n);
        [v1] = NDCT_II_IV(u1,0);
        [v2] = NDCT_II_IV(u2,0);
          z  = [v1',v2']';
          w  = ModAdd(z);
          y  = PermT(w);
     end;      
  otherwise, y=0; disp('Es ist ein Fehler bei der Eingabe von k aufgetreten.');
  end;
  
