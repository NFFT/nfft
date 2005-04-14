function [y] = NDCT_III_IV(x,k)

% Funktion NDCT_III_IV
%
% Aufruf: [y] = NDCT_III_IV(x,k)  -------  Vergessen Sie nicht, 
%                                           den richtigen Pfad einzustellen !!!!!
%
% ================================================
% Rekursiver DCT-III(n) bzw. DCT-IV(n) Algorithmus : 
% ================================================
% 
% Eingabeparameter:   x     zu transformierender Vektor der Länge n (= Zweierpotenz)
%                     k     Typ der Transformation (3 oder 4)
%
% Ausgangsparameter:  y     sqrt{n} * C_n^{III} * x     im Fall k=0
%                           sqrt(n) * C_n^{IV} * x      im Fall k=1
%

% ===========================================================
% An dieser Stelle müssen die benötigten Werte bzw. Matrizen 
% zur Verfügung gestellt werden :
% ===========================================================


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
     case 1, y = ButterflyT(x);
     otherwise, t=2; u = Perm(x); u1 = u(1:n1); u2 = u(n1+1:n); 
        [v1] = NDCT_III_IV(u1,0);
        [v2] = NDCT_III_IV(u2,1);
          v  = [v1',v2']';
          y  = ButterflyT(v);
     end;      
  case 1, switch n1,
     case 1, t=1; y = C4 * x;
     otherwise, w = Perm(x); u = ModAddT(w); u1 = u(1:n1); u2 = u(n1+1:n);
        [v1] = NDCT_III_IV(u1,0);
        [v2] = NDCT_III_IV(u2,0);
          v  = [v1',v2']';
          y  = TwiddleT(v);
     end;      
  otherwise, y=zeros(n,1); disp('Es ist ein Fehler bei der Eingabe von k aufgetreten.');
  end
  
