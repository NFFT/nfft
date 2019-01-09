classdef nfsoftTestcaseDelegate
  %NFSFTTESTCASEDELEGATE Summary of this class goes here
  properties(Hidden=true,SetAccess='protected',GetAccess='public')
    N = [];
    M = [];
    x = [];
    f_hat = [];
    f = [];
  end
  
  methods
    function h = setup(h)
    end
    
    function h = destroy(h)
    end
  end
  
end
