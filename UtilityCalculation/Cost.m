% ================================= %
%            Cost.m                 %
%           Wenbin Liu              %
%       wl2581@columbia.edu         %
% ================================= %

classdef (Abstract) Cost
    properties(Access = private)
        % Python code for metaclass: __metaclass__ = ABCMeta
        % costMeta is meta.class object in MATLAB
        costMeta = ?Cost;
    end
    
    methods (Abstract = true)
        r1 = cost(obj);
        r2 = price(obj);
    end
    
    methods
        %constructor
        function obj = Cost()
        end
    end
 end