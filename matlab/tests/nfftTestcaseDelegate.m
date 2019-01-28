classdef nfftTestcaseDelegate
  %NFFTTESTCASEDELEGATE Summary of this class goes here
  %   Detailed explanation goes here
  %
  %   properties(Hidden=true,SetAccess='private',GetAccess='private');
  %     setup_function = [];
  %     destroy_function = [];
  %   end
  %
  %   methods
  %     function h = nfft(setup, destroy)
  %       if isa(setup, 'function_handle')
  %         h.setup_function = setup;
  %       end
  %       if isa(destroy, 'function_handle')
  %         h.destroy_function = destroy;
  %       end
  %     end
  %
  %     function setup(h,d, N, NN, M, x, f_hat, f)
  %       if isa(h.setup_function, 'function_handle')
  %         h.setup_function(d, N, NN, M, x, f_hat, f);
  %       end
  %     end
  %
  %     function destroy(h,x, f_hat, f)
  %       if isa(h.destroy_function, 'function_handle')
  %         h.destroy(x, f_hat, f);
  %       end
  %     end
  %   end
  properties(Hidden=true,SetAccess='protected',GetAccess='public')
    d = [];
    N = [];
    NN = [];
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