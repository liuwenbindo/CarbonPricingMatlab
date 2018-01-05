% ================================= %
%          TreeModel.m              %
%           Wenbin Liu              %
%       wl2581@columbia.edu         %
% ================================= %

% Assume following given:
%   decision_times: array
%   prob_scale:     int=1.0

classdef TreeModel
    properties
        decision_times = []
        final_state_probs = []
        all_node_probs = []
    end

    properties(SetAccess = private)
        num_periods
        num_decision_nodes
        num_final_states
    end

    properties(Constant)
      prob_scale = 1.0
    end

    methods(Access = private)
        function obj = set_otherprops(obj)
            % Dependent property 1
            obj.num_periods = length(obj.decision_times)-1;

            % Dependent property 2
            obj.num_decision_nodes = (2^obj.num_periods)-1;

            % Dependent property 3
            obj.num_final_states = 2^(obj.num_periods-1);

            % Dependent property 4
            temp_final_state = zeros(1, obj.num_final_states);
            temp_final_state(1) = 1.0;
            next_prob = 1.0;
            n = obj.num_final_states;

            for i = 1:n
                next_prob = next_prob * obj.prob_scale^(1/i);
                temp_final_state(i) = next_prob;
            end
            temp_final_state = temp_final_state / sum(temp_final_state);
            obj.final_state_probs = temp_final_state;

            % Dependent property 5
            temp_node_probs = zeros(1, obj.num_decision_nodes);
            temp_node_probs(obj.num_final_states:end) = obj.final_state_probs;

            for period = obj.num_periods-2:-1:0
                for state  =  0:(2^period-1) %state begins with 0
                    node_pos =  obj.get_node(period, state); %node_pos begins with 0
                    temp_node_probs(node_pos + 1) = temp_node_probs(2*node_pos + 2) + temp_node_probs(2*node_pos + 3);
                end
            end
            obj.all_node_probs = temp_node_probs;
        end
    end


    methods
        %constructor
        function obj = TreeModel() %array)
            %obj.decision_times = array;
        end

        % set method for decision_times
        function obj = set.decision_times(obj, arr)
            obj.decision_times = arr;
            obj = obj.set_otherprops;
        end

        function obj = set.final_state_probs(obj, v)
            obj.final_state_probs = v;
        end

        function obj = set.all_node_probs(obj, v)
            obj.all_node_probs = v;
        end

%         % set method for all_node_probs
%         function obj = set.all_node_probs(obj, value)
%             obj.all_node_probs = value;
%         end

%         %dependent property 1: num_periods
%         function r = get.num_periods(obj)
%             r = length(obj.decision_times)-1;
%         end

%         %dependent property 2: num_decision_nodes
%         function r = get.num_decision_nodes(obj)
%             r = (2^obj.num_periods)-1;
%         end

%         %dependent property 3: num_final_states
%         function r = get.num_final_states(obj)
%             r = 2^(obj.num_periods-1);
%         end

%         %method 1.1: _create_prob: final states
%         function r = get.final_state_probs(obj)
%             temp_final_state = zeros(1, obj.num_final_states);
%             temp_final_state(1) = 1.0;
%             next_prob = 1.0;
%             n = obj.num_final_states;
%
%             for i = 1:n
%                 next_prob = next_prob * obj.prob_scale^(1/i);
%                 temp_final_state(i) = next_prob;
%             end
%             temp_final_state = temp_final_state / sum(temp_final_state);
%             r = temp_final_state;
%         end

%         %method 1.2: _create_prob: all states
%         function r = get.all_node_probs(obj)
%             temp_node_probs = zeros(1, obj.num_decision_nodes);
%             temp_node_probs(obj.num_final_states:end) = obj.final_state_probs;
%
%             for period = obj.num_periods-2:-1:0
%                 for state  =  0:(2^period-1) %state begins with 0
%                     node_pos =  obj.get_node(period, state); %node_pos begins with 0
%                     temp_node_probs(node_pos + 1) = temp_node_probs(2*node_pos + 2) + temp_node_probs(2*node_pos + 3);
%                 end
%             end
%             r = temp_node_probs;
%         end

        %method 2: get_num_nodes_period
        function r = get_num_nodes_period(obj, period)
            if period >= obj.num_periods
                r = 2^(obj.num_periods - 1);
            else
                r = 2^(period);
            end
        end

        %method 3: get_nodes_in_period
        function r = get_nodes_in_period(obj, period)
            temp_node_range = [0,0];
            if period >= obj.num_periods
                period = obj.num_periods - 1;
            end
            num_nodes = obj.get_num_nodes_period(period);
            temp_node_range(1) = obj.get_node(period, 0);
            temp_node_range(2) = temp_node_range(1) + num_nodes - 1;
            r = temp_node_range;
        end

        %method 4: get_node
        function r = get_node(obj, period, state)
            if period > obj.num_periods
                error('Error: given period is larger than number of periods!');
            end
            if state >= 2^period
                error('Error: no such state in period %d.',period);
            end
            r = 2^period + state - 1;
        end

        %method 5: get_state
        %note: number of node starts from 0, number of state start from 0
        function r = get_state(obj, node)
            if node > obj.num_decision_nodes
                r =  node - obj.num_decision_nodes;
            end
            %Q: There's no need of input "period"
            period = obj.get_period(node);
            r = node - (2^period - 1);
        end

        %method 6: get_period
        function r = get_period(obj, node)
            if node > obj.num_decision_nodes
                r =  obj.num_periods;
            end
            for period = 0:obj.num_periods
                if node >= 2^period - 1 && node <= 2^(period+1) - 2
                    r = period;
                    break;
                end
            end
        end

        %method 7: get_parentnode
        function r = get_parentnode(obj, child)
            if child == 0
                r = 0;
            end
%             if child > obj.num_decision_nodes
%                 r = child - obj.num_final_states; %Q: Why?
%             end
            if rem(child, 2) == 0 && child < obj.num_decision_nodes
                r = (child - 2) / 2;
            elseif child < obj.num_decision_nodes
                r = (child - 1)/2;
            elseif child >= obj.num_decision_nodes
                r = child - obj.num_final_states;
            end
        end

        %method 8: get_path
        function r = get_path(obj, node)
            period = obj.get_period(node);
            temp_path = zeros(1, period+1);
            temp_path(1) = node;
            for i = 1:period
                temp_path(i+1) = obj.get_parentnode(temp_path(i));
            end
            r = fliplr(temp_path);
        end

        %method 9: get_probs_in_period
        function r = get_probs_in_period(obj, period)
            temp_node_range = obj.get_nodes_in_period(period);
            %r =  temp_node_range;
            r = obj.all_node_probs( (temp_node_range(1)+1) : (temp_node_range(2)+1) );
        end

        %method 10: reachable_end_state
        function r = reachable_end_state(obj,node)
            period =  obj.get_period(node);
            state = obj.get_state(node);
            temp = zeros(1,2);
            if period >= obj.num_periods
                temp(1) = node - obj.num_decision_nodes;
                temp(2) = temp(1);
            else
                k = (obj.num_final_states)/(2^period);
                temp(1) = k*(state);
                temp(2) = k*(state + 1) - 1;
            end
            r = temp;
        end
    end
end
