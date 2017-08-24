% ================================= %
%        BigStorageTree.m           %
%           Wenbin Liu              %
%       wl2581@columbia.edu         %
% ================================= %

% Assume following given:
%   subinterval_len: int
%   decision_times: array  

classdef BigStorageTree < BaseStorageTree   
    properties
         first_period_intervals
         subinterval_len
    end
    
    methods
         function obj = BigStorageTree(subinterval_len, array) 
            obj@BaseStorageTree(array);
            obj.subinterval_len = subinterval_len;
            obj.periods = 0:obj.subinterval_len:(obj.decision_times(end));
            obj.init_tree()
        end
        function r = get.first_period_intervals(obj)
            r = (obj.decision_times(2) - obj.decision_times(1))/obj.subinterval_len;
        end
        
        function r = get_next_period_array(obj, period)
            if period + obj.subinterval_len <= obj.decision_times(end)
                r = obj.tree(period + obj.subinterval_len);
            else
                Error('Period is not a valid period or too large!');
            end      
        end
        
        % Note: as MATLAB index starts with 1, so the return of this
        % function will be: 1 + returnOfPythonCode
        function r = between_decision_times(obj, period)
            if period == 0
                r = 0;
            end
            for i = 1:(length(obj.information_times))
                if obj.decision_times(i) <= period && period < obj.decision_times(i+1)
                    r = i;
                end                
            end            
        end
        
        % Note: as MATLAB index starts with 1, so the return of this
        % function will be: 1 + returnOfPythonCode
        function r = decision_interval(obj, period)
            if period == 0
                r = 0;
            end
            
            for i = 2:length(obj.decision_times) 
                if obj.decision_times(i-1) < period && period <= obj.decision_times(i)
                    r = i;
                end               
            end
        end
    end
    
end

