function s = char(p)
% Conversion to string for f_hat class
%
% Copyright 2007 Jens Keiner

if (p.N == -1)
  s = 'empty f_hat object';
else
  s = ['eegree ',num2str(p.N),' f_hat object'];
end
