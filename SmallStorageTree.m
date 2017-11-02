% ================================= %
%        SmallStorageTree.m         %
%           Wenbin Liu              %
%       wl2581@columbia.edu         %
% ================================= %

% Assume following given:
%   decision_times: array

classdef SmallStorageTree < BaseStorageTree
    
    properties        
        last =[]
        nodes
        tree  = containers.Map('KeyType','int32', 'ValueType','any')
    end
    
    properties(Access = private)
        tree_p  = containers.Map('KeyType','int32', 'ValueType','any')
    end
  
    methods 
        
        function obj = SmallStorageTree(array)
            obj@BaseStorageTree(array);
            obj.periods = obj.decision_times;
            obj.init_tree();
        end 
        
        % GET for tree property
        function r = get.tree(obj)
            r = obj.tree_p;
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
                for i = 2:length(obj.decision_times)
                    if period == obj.decision_times(i)
                        j = i;                
                    end
                end
                r = obj.decision_times(j-1);
            else
                Error('Period not in decision times or first period!');
            end            
        end
        
        
%         % From BaseStorageTree
%         
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
       
        
        %method 1
        function obj = set_value(obj, period, values)
            if ismember(period, obj.periods) == 0
                error('Not a Valid Period!');              
            end
            if size(values) ~= size(obj.tree(period))
                 error('Shape %s and Shape %s are not aligned.', size(values), size(obj.tree(period)));
            end
            % Need refinement, but it works now
            obj.tree_p(period) = values;
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
%         
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
  
    methods (Access = protected)
        
        function init_tree(obj)
            i = 0;
            for key = 1:length(obj.periods)
                obj.tree_p(obj.periods(key)) = zeros(2^i, 1); %every value is one row vector
                if ismember(obj.periods(key), obj.information_times) == 1
                    i = i + 1;
                end
            end          
        end          
    end
    
end