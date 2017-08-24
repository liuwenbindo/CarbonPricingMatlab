% ================================= %
%       BaseStorageTree.m           %
%           Wenbin Liu              %
%       wl2581@columbia.edu         %
% ================================= %

% Assume following given:
%   decision_times: array

classdef (Abstract) BaseStorageTree
    properties(Access = private)
        %Python code for metaclass: __metaclass__ = ABCMeta
        %baseStorageMeta is meta.class object in MATLAB
        baseStorageMeta = ?BaseStorageTree;
    end
    
    properties
        decision_times = []
        information_times = []
        periods = []
        tree = containers.Map('KeyType','int32', 'ValueType','any')
        last = []
        last_period
        nodes
    end

    methods (Access = private)
        %private method 1
        function r = len(obj)
            r =  size(obj.tree, 1); %number of items in tree
        end
        %private method 2
        function r = get_item(obj, key)
            if rem(key,1) == 0 || isfloat(key)
                r = obj.tree(key);
            else
                error('Error: Index must be int, not %d.',class(key));
            end    
        end
    end

    methods (Access = protected)
        
        function init_tree(obj)
            i = 0;
            for key = 1:length(obj.periods)
                obj.tree(obj.periods(key)) = zeros(2^i, 1); %every value is one row vector
                if ismember(obj.periods(key), obj.information_times) == 1
                    i = i + 1;
                end
            end          
        end          
    end

    methods (Abstract = true)
        r = get_next_period_array(obj)       
    end

    methods
        %constructor
        function obj = BaseStorageTree(array)
            obj.decision_times =  array;
            obj.information_times = obj.decision_times(1:end-2);
        end
        %set property 1
        function r = get.last(obj)
            r =  obj.tree(obj.decision_times(end));
        end
        %set property 2
        function r = get.last_period(obj)
            r = obj.decision_times(end);
        end
        %set property 3
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
        %method 1
        function set_value(obj, period, values)
            if ismember(period, obj.periods) == 0
                error('Not a Valid Period!');              
            end
            if size(values) ~= size(obj.tree(period))
                 error('Shape %s and Shape %s are not aligned.', size(values), size(obj.tree(period)));
            end
            % Need refinement, but it works now
            obj.tree(period) = values;
        end
        
        %method 2
        function r = is_decision_period(obj, time_period)
            r = ismember(time_period, obj.decision_times);
        end
        
        %method 3
        function r = is_real_decision_period(obj, time_period)
            r =  ismember(time_period, obj.decision_times(1:end-1));
        end
        
        % method 4
        % USE "find_path" function in tools
        function write_tree(obj, filename, header, delimeter)
            real_times = obj.decision_times(1:end-1);
            size = length(obj.tree(real_times(end)));
            output_list = [];
            prev_k = size;
            
            for t = 1:length(real_times)
                 temp_list = zeros(size*2, 1);
                 current_time = real_times(t);
                 k = round( size / length(obj.tree(current_time)) );
                 temp_list(k+1 : prev_k : length(temp_list)) = obj.tree(current_time);
                 % Append temp_list after output_list
                 output_list(:,t) = temp_list;
                 prev_k = k;
            end
         
            % ??? Question mark here. Correct or not?
            fopen(filename,'wt');
            dlmwrite(filename, header, delimeter);
            dlmwrite(filename, output_list, delimeter,'-append');
        end
        
        % method 5
        % USE "write_column_csv" function in tools
        function write_column(obj, filename, header, start_year, delimeter)
            start_year = 2015;
            delimeter = ';';
            
            if exists(filename,'file')
                obj.write_column_existing(filename, header);
            else
               real_times = obj.decision_times(1:end-1);
               
               years = zeros(length(real_times),1);
               
               temp_nodes = [];
               
               tmp_treevalue = values(obj.tree);
               
                %dims = [];
                %for t=1:length(real_times)
                %    dims(t,1) = 2^(t-1);
                %end
                %output_list = cell(dims);
               output_list = [];
               
               k = 0;
               for t = 1:length(real_times)
                   years(t) = start_year + real_times(t);
                   curnt_value = tmp_treevalue{t};
                   %temp_value = values(obj.tree, t);
                   
                   for n = 1:length(obj.tree(t))                      
                       temp_nodes = [temp_nodes, k]; 
                       output_list = [output_list, curnt_value(n)];
                       k = k+1;                
                   end                   
               end
               final_header = ['Year','Node',header];
               fopen(filename,'wt');
               dlmwrite(filename, final_header, delimeter);
               dlmwrite(filename, output_list, delimeter,'-append'); % index = [year, node]
            end
        end
        
        % method 6
        function write_column_existing(obj, filename, header, delimeter)
            if nargin < 4 || isempty(delimeter)
                delimeter = ';';
            end
            
            output_list = [];
            tmp_treevalue = values(obj.tree);
            
             for t = 1:length(real_times)
                 curnt_value = tmp_treevalue{t};               
                 output_list = [output_list, curnt_value];
             end
             
             fopen(filename,'wt');
             dlmwrite(filename, header, delimeter);
             dlmwrite(filename, output_list, delimeter,'-append');
        end
    end  
end