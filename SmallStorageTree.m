% ================================= %
%        SmallStorageTree.m         %
%           Wenbin Liu              %
%       wl2581@columbia.edu         %
% ================================= %

% Assume following given:
%   decision_times: array

classdef SmallStorageTree < BaseStorageTree
    properties
    end
  
    methods 
        function obj = SmallStorageTree(array)
            obj@BaseStorageTree(array);
            obj.periods = obj.decision_times;
            obj.init_tree();
        end         
        function r = get_next_period_array(obj, period)
        	if obj.is_real_decision_period(period)   
                for i=1:length(obj.decision_times)
                    if period == obj.decision_times(i)
                        j = i;                
                    end
                end
                index = obj.decision_times(j+1);
                r = obj.tree(index);
            else
                Error('Given period is not in real decision times!');
            end
        end
        
        function r = index_below(obj, period)
            if ismember(period, obj.decision_times(2:end)) == 1
                for i = 1:(length(obj.decision_times)-1)
                     if period == obj.decision_times(i)
                        j = i;                
                    end
                end
                r = obj.decision_times(j-1);
            else
                Error('Period not in decision times or first period!');
            end            
        end
    end
end