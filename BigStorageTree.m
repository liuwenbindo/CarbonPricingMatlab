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
         last = []
         nodes
         tree  = containers.Map('KeyType','int32', 'ValueType','any')
    end   
    
    
    properties(Access = private)
        tree_p  = containers.Map('KeyType','int32', 'ValueType','any')
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
        
        
        function r = get.tree(obj)
              r = obj.tree_p;
        end
        
        
        function obj = set.tree_p(obj, value)
            obj.tree_p = value;            
        end
        

        function r = get_next_period_array(obj, period)
            if period + obj.subinterval_len <= obj.decision_times(end)
                r = obj.tree(period + obj.subinterval_len);
            else
                Error('Period is not a valid period or too large!');
            end      
        end
               

        function obj = set_value(obj, period, values)
            if ismember(period, obj.periods) == 0
                error('Not a Valid Period!');              
            end
            if size(values) ~= size(obj.tree_p(period))
                 error('Shape %s and Shape %s are not aligned.', size(values), size(obj.tree_p(period)));
            end
            % Need refinement, but it works now          
            obj.tree_p(period) = values;
        end

                
        function r = between_decision_times(obj, period)
            
            if period == 0
                r = 0;
            end
            
            for i = 1:(length(obj.decision_times)-1)
                if (obj.decision_times(i) <= period) && (period < obj.decision_times(i+1))
                    r = i-1;
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
        
        
        % From BaseStorageTree
        
        % set property 3
        function r =  get.nodes(obj)
            temp_tree = obj.tree;
            n = 0;
            key = keys(temp_tree);
            value = values(temp_tree);
            for i = 1:length(key)
                n = n + length(value{i});
            end
            r= n;
        end
        
        %set property 1
        function r = get.last(obj)
            r =  obj.tree(obj.decision_times(end));
        end

        
%         % method 5
%         % USE "write_column_csv" function in tools
%         function write_column(obj, filename, header, start_year, delimeter)
%             start_year = 2015;
%             delimeter = ';';
%             
%             if exists(filename,'file')
%                 obj.write_column_existing(filename, header);
%             else
%                real_times = obj.decision_times(1:end-1);
%                
%                years = zeros(length(real_times),1);
%                
%                temp_nodes = [];
%                
%                tmp_treevalue = values(obj.tree);
%                
%                 %dims = [];
%                 %for t=1:length(real_times)
%                 %    dims(t,1) = 2^(t-1);
%                 %end
%                 %output_list = cell(dims);
%                output_list = [];
%                
%                k = 0;
%                for t = 1:length(real_times)
%                    years(t) = start_year + real_times(t);
%                    curnt_value = tmp_treevalue{t};
%                    %temp_value = values(obj.tree, t);
%                    
%                    for n = 1:length(obj.tree(t))                      
%                        temp_nodes = [temp_nodes, k]; 
%                        output_list = [output_list, curnt_value(n)];
%                        k = k+1;                
%                    end                   
%                end
%                final_header = ['Year','Node',header];
%                fopen(filename,'wt');
%                dlmwrite(filename, final_header, delimeter);
%                dlmwrite(filename, output_list, delimeter,'-append'); % index = [year, node]
%             end
%         end
        
%         % method 6
%         function write_column_existing(obj, filename, header, delimeter)
%             if nargin < 4 || isempty(delimeter)
%                 delimeter = ';';
%             end
%             
%             output_list = [];
%             tmp_treevalue = values(obj.tree);
%             
%              for t = 1:length(real_times)
%                  curnt_value = tmp_treevalue{t};               
%                  output_list = [output_list, curnt_value];
%              end
%              
%              fopen(filename,'wt');
%              dlmwrite(filename, header, delimeter);
%              dlmwrite(filename, output_list, delimeter,'-append');
%         end
        
    end
    
%     methods (Access = private)
%         %private method 1
%         function r = len(obj)
%             r =  size(obj.tree, 1); %number of items in tree
%         end
%         %private method 2
%         function r = get_item(obj, key)
%             if rem(key,1) == 0 || isfloat(key)
%                 r = obj.tree(key);
%             else
%                 error('Error: Index must be int, not %d.',class(key));
%             end    
%         end
%     end
  
    methods (Access = private)      
        function init_tree(obj)
            i = 0;
            for key = 1:length(obj.periods)
                obj.tree_p(obj.periods(key)) = zeros(2^i, 1); %every value is one row vector
                %obj.tree(obj.periods(key)) = zeros(2^i, 1);
                if ismember(obj.periods(key), obj.information_times) == 1
                    i = i + 1;
                end
            end          
        end          
    end
    
end

