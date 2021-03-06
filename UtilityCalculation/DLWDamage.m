classdef DLWDamage < Damage

    properties
        ghg_levels
        cons_growth
        dnum
        subinterval_len
        cum_forcings
        d
        d_rcomb
        emit_pct
        damage_coefs
        d_import_flag
    end

    methods

        % Constructor
        function obj = DLWDamage(tree, bau, cons_growth, ghg_levels, subinterval_len)
            obj@Damage(tree, bau);
            obj.ghg_levels = ghg_levels;
            obj.cons_growth = cons_growth;
            obj.dnum = length(ghg_levels);
            obj.subinterval_len = subinterval_len;
            obj.cum_forcings = NaN;
            obj.d = NaN;
            obj.d_rcomb = NaN;
            obj.emit_pct = NaN;
            obj.damage_coefs = NaN;
            obj.d_import_flag = 0;
        end


        % Get/Set Methods for Properties
%         % Get method of property d_rcomb
%         function r = get.d_rcomb(obj)
%             r = obj.d;
%         end
        % Set method of property d_rcomb
        function obj = set.d_rcomb(obj,value)
            obj.d_rcomb = value;
        end
        % Set method of property emit_pct
        function obj = set.emit_pct(obj,value)
            obj.emit_pct = value;
        end
        % Set method of property d
        function obj = set.d(obj,value)
            obj.d = value;
        end
        % Set method of property cum_forcings
        function obj = set.cum_forcings(obj,value)
            obj.cum_forcings = value;
        end
        % Set method of property damage_coefs
        function obj = set.damage_coefs(obj,value)
            obj.damage_coefs = value;
        end
        % Set method of property d_import_flag
        function obj = set.d_import_flag(obj,value)
            obj.d_import_flag = value;
        end


        function obj = recombine_nodes(obj)

%       Creating damage coefficients for recombining tree. The state reached by an up-down move is
% 		separate from a down-up move because in general the two paths will lead to different degrees of
% 		mitigation and therefore of GHG level. A 'recombining' tree is one in which the movement from
% 		one state to the next through time is nonetheless such that an up move followed by a down move
% 		leads to the same fragility.

            nperiods = obj.tree.num_periods;
            sum_class = zeros(1, nperiods);
            new_state = {};

            % Shallow copy in original Python script, should I keep it?
            temp_prob = obj.tree.final_state_probs;
            obj.d_rcomb = obj.d;

            for old_state = 1:obj.tree.num_final_states
                temp = old_state;
                n = nperiods-2; % last period before recombining
                d_class = 0;

                while n >= 0
                    if temp - 1 >= 2^n % modify the lower half of all final states
                        temp = temp - 2^n;
                        d_class = d_class + 1;
                    end
                    n = n - 1;
                end
                sum_class(d_class + 1) = sum_class(d_class + 1) + 1;

                new_state{d_class + 1}(sum_class(d_class + 1)) = old_state-1;
            end

%           new_state:
%           0
%           1  2  4  8  16
%           3  5  6  9  10 12 17 18 20 24
%           7 11  13  14 19 21 22 25 26 28
%           15 23 27 29 30
%           31
%           sum_class: [1 5 10 10 5 1]
%           sum_nodes: [ 0, 1, 6, 16, 26, 31, 32]
%           prob_sum: sums up the probabilities of final states partitioned by sum_nodes

            sum_nodes = [];
            sum_nodes(1) = 0;
            sum_nodes(2:length(sum_class)+1) = cumsum(sum_class);

            prob_sum = zeros(1, length(sum_nodes)-1);
            for i = 1:length(sum_nodes)-1
                prob_sum(i) = sum( (obj.tree.final_state_probs(sum_nodes(i)+1 : sum_nodes(i+1)))');
            end

%           Please uncomment below for testing
%              sum_class
%              sum_nodes
%              new_state
%              new_state{2}
%              new_state{3}
%              new_state{4}
%              new_state{5}
%              prob_sum

             for period = 1:nperiods   %look into each period
                 for k = 1:obj.dnum   %look into each simulated scenario
                     d_sum = zeros(1, nperiods);    %store the probability-weighted sum of the simulated damage for each category
                     old_state = 1;    %direct to the first node in the the category set, signal


                     for d_class = 1 : nperiods
                          d_sum(d_class) = sum( (obj.tree.final_state_probs( old_state : old_state+sum_class(d_class)-1 ) ...
						 			                              .* obj.d_rcomb(old_state : old_state+sum_class(d_class)-1, period, k)')' );
                          old_state = old_state + sum_class(d_class);

                          % SET.final_state_probs (TreeModel.m)
                          final_state_idx  = new_state{d_class}(1:sum_class(d_class)) + 1;
                          obj.tree.final_state_probs( final_state_idx ) = temp_prob(1);
                     end

                     for d_class = 1:nperiods
                         % SET.d_rcomb
                         obj.d_rcomb(new_state{d_class}(1:sum_class(d_class))+1, period, k) = d_sum(d_class) / prob_sum(d_class);
                         % Alternative way:
                         % tmp_d_rcomb = obj.d_rcomb;
                         % tmp_d_rcomb(k, new_state(d_class, 1:sum_class(class)), period) = d_sum(d_class) / prob_sum(d_class);
                         % obj.d_rcomb = tmp_d_rcomb;
                     end
                 end
             end

             % After final_states_prob being updated
             last_num_elements = length(obj.tree.final_state_probs);

             % find the probability-weighted average damage
             % SET.tree; SET.all_node_probs
             tmp_tree = obj.tree;
             tmp_tree.all_node_probs(end-(last_num_elements-1):end) = obj.tree.final_state_probs;
             obj.tree = tmp_tree;

             for p =2:nperiods-1
                 nodes = obj.tree.get_nodes_in_period(p-1);
                 % look into all the nodes for the period
                 for node = nodes(1):nodes(2)
                    end_states = obj.tree.reachable_end_state(node); % determine reachable final states
                    worst_end_state = end_states(1);
                    best_end_state = end_states(2);

                    % SET.tree
                    tmp_tree = obj.tree;
                    tmp_tree.all_node_probs(node+1) = sum( obj.tree.final_state_probs(worst_end_state+1:best_end_state+1)' );
                    obj.tree = tmp_tree;
                 end
             end
        end


        function obj = damage_interpolation(obj)

        % Modified based on Python version.
        % Create the interpolation coefficients used in `damage_function`.

            if obj.d_import_flag == 0
                %fprintf('Importing stored damage simulation...\n');
                obj = obj.import_damages();
                %fprintf('Finished importing.\n')
            end

            obj = obj.recombine_nodes();

            if isnan(obj.emit_pct)
                bau_emission = obj.bau.ghg_end - obj.bau.ghg_start;
                tmp_emit_pct = 1.0 - (obj.ghg_levels - obj.bau.ghg_start) / bau_emission;
                obj.emit_pct = tmp_emit_pct;
            end

            tmp_damage_coefs = zeros(obj.dnum-1, obj.dnum, obj.tree.num_periods, obj.tree.num_final_states);
            amat = ones(obj.dnum, obj.dnum, obj.tree.num_periods);
            bmat = ones(obj.tree.num_periods, obj.dnum);

            tmp_damage_coefs(obj.dnum-1, obj.dnum, :, :) = obj.d_rcomb(:, :, obj.dnum)';
            tmp_damage_coefs(obj.dnum-1, obj.dnum-1, :, :) = (obj.d_rcomb(:, :, obj.dnum-1)' - obj.d_rcomb(:, :, obj.dnum)') / obj.emit_pct(1, end-1);
            % for p = 1:obj.tree.num_periods
            %     for s = 1:obj.tree.num_final_states
            %         tmp_damage_coefs(obj.dnum-1, obj.dnum, p, s) = obj.d_rcomb(s, p, obj.dnum);
            %         tmp_damage_coefs(obj.dnum-1, obj.dnum-1, p, s) = (obj.d_rcomb(s, p, obj.dnum-1) - obj.d_rcomb(s, p, obj.dnum)) / obj.emit_pct(1, end-1);
            %     end
            % end
            obj.damage_coefs =  tmp_damage_coefs;

            amat(1, 1, :) = 2.0 * obj.emit_pct(end-1);
            amat(1, end, :) = 0.0;

            for i = 1:obj.tree.num_periods
                amat(2:end, 1, i) = (obj.emit_pct(1:end-1).^2)';
                amat(2:end, 2, i) = (obj.emit_pct(1:end-1))';
            end

            for state = 1:obj.tree.num_final_states
              bmat(:, 1) = obj.damage_coefs(end, end-1, :, state) * obj.emit_pct(end-1);
              bmat(:, 2:end) = obj.d_rcomb(state, :, 1:end-1);
              for j = 1:obj.tree.num_periods
                 obj.damage_coefs(1, :, j, state) = linsolve(amat(:,:,j), bmat(j,:)');
              end
            end
        end


        function obj = import_damages(obj, file_name)

%       Import saved simulated damages. File must be saved in 'data' directory
% 		inside current working directory. Save imported values in `d`.
%
% 		Parameters
% 		----------
% 		file_name : str, optional
% 			name of file of saved simulated damages
%
% 		Raises
% 		------
% 		IOError
% 			If file does not exist.

            if nargin == 1 || isempty(file_name)
                file_name = 'simulated_damages.csv';
            end

            %fileID = fopen(file_name);

            try
                d_fromcsv = dlmread('simulated_damages.csv',';');
            catch
                message = ferror(fileID);
                error(message);
            end

            %fclose(file_name);

            n = obj.tree.num_final_states;
            tmp_d = zeros(n,6,obj.dnum);
            for i = 1:obj.dnum
                tmp_d(:,:,i) = d_fromcsv(n*(i-1)+1 : n*(i),:);
            end
            obj.d = tmp_d;
            obj.d_import_flag = 1;

        end


        function r = damage_simulation(	obj, draws, peak_temp, disaster_tail, tip_on, ...
                                        temp_map, temp_dist_params, maxh, save_simulation)

%       Initializion and simulation of damages, given by :mod:`ez_climate.DamageSimulation`.
%
% 		Parameters
% 		----------
% 		draws : int
% 			number of Monte Carlo draws
% 		peak_temp : float, optional
% 			tipping point parameter
% 		disaster_tail : float, optional
% 			curvature of tipping point
% 		tip_on : bool, optional
% 			flag that turns tipping points on or off
% 		temp_map : int, optional
% 			mapping from GHG to temperature
% 		        * 0: implies Pindyck displace gamma
% 		        * 1: implies Wagner-Weitzman normal
% 		        * 2: implies Roe-Baker
% 		        * 3: implies user-defined normal
% 		        * 4: implies user-defined gamma
% 		temp_dist_params : ndarray or list, optional
% 			if temp_map is either 3 or 4, user needs to define the distribution parameters
% 		maxh : float, optional
% 			time parameter from Pindyck which indicates the time it takes for temp to get half
% 		        way to its max value for a given level of ghg
% 		cons_growth : float, optional
% 			yearly growth in consumption
% 		save_simulation : bool, optional
% 			True if simulated values should be save, False otherwise
%
% 		Returns
% 		-------
% 		ndarray
% 			simulated damages

             if nargin == 2
                 peak_temp=9.0;
                 disaster_tail=12.0;
                 tip_on=true;
                 temp_map=1;
                 temp_dist_params=[];
                 maxh=100.0;
                 save_simulation=true;
             end

             % get simulated damage from damage_simulation class
             ds = DamageSimulation(obj.tree, obj.ghg_levels, peak_temp, disaster_tail, tip_on, temp_map,...
                                    temp_dist_params, maxh, obj.cons_growth);

             fprintf('Starting damage simulation...');
             % self.d, self.parameter_list= ds.simulate(draws,self.change, write_to_file = False)
             tmp_d = ds.simulate(draws, save_simulation);
             obj.d = tmp_d;
             fprintf('Done.');
             obj.damage_interpolation();
             r = obj.d;

        end


        function obj = forcing_init(obj)
        % Initialize `Forcing` object and cum_forcings used in calculating the force mitigation up to a node.

            if isnan(obj.emit_pct)
                bau_emission = obj.bau.ghg_end - obj.bau.ghg_start;
                tmp_emit_pct = 1.0 - (obj.ghg_levels - obj.bau.ghg_start) / bau_emission;
                obj.emit_pct = tmp_emit_pct;
            end

            % array like [x,y] contains number of the periods and the path
            tmp_cum_forcings = zeros(obj.tree.num_periods, obj.dnum);

            tmp_miti = ones(obj.dnum, obj.tree.num_decision_nodes);
            for i = 1:obj.dnum
                tmp_miti(i, :) = tmp_miti(i, :) * obj.emit_pct(i);
            end
            mitigation = tmp_miti;
            % mitigaion is indexed with number of path (row) and number of nodes(column), and is the emit_pct

            % look into each simulation
            for i  = 1:obj.dnum
              for n = 2:obj.tree.num_periods+1
                  node = obj.tree.get_node(n-1, 0);
                  % return the forcing on the node
                  forcingobj = Forcing();
                  tmp_cum_forcings(n-1, i) = forcingobj.forcing_at_node(mitigation(i,:), node, obj.tree, obj.bau, obj.subinterval_len);
              end
            end
            obj.cum_forcings = tmp_cum_forcings;
        end


        function [obj, r] = forcing_based_mitigation(obj, forcing, period)

%       Calculation of mitigation based on forcing up to period. Interpolating between the forcing associated
% 		with the constant degree of mitigation consistent with the damage simulation scenarios.

        % this whole function is based on a new theory
            if period == 0
                p = period;
            else
                p = period - 1;
            end

            %obj = obj.forcing_init();

            if forcing > obj.cum_forcings(p+1, 2)
                weight_on_sim2 = (obj.cum_forcings(p+1, 3) - forcing) / (obj.cum_forcings(p+1, 3) - obj.cum_forcings(p+1, 2));
                weight_on_sim3 = 0;
            elseif forcing > obj.cum_forcings(p+1, 1)
                weight_on_sim2 = (forcing - obj.cum_forcings(p+1,1)) / (obj.cum_forcings(p+1,2) - obj.cum_forcings(p+1,1));
                weight_on_sim3 = (obj.cum_forcings(p+1, 2) - forcing) / (obj.cum_forcings(p+1, 2) - obj.cum_forcings(p+1, 1));
            else
                weight_on_sim2 = 0;
                weight_on_sim3 = 1.0 + (obj.cum_forcings(p+1, 1) - forcing) / obj.cum_forcings(p+1,1);
            end

            r = weight_on_sim2 * obj.emit_pct(2) + weight_on_sim3 * obj.emit_pct(1);
        end


        function r = average_mitigation_node(obj, m, node, period)

%       Calculate the average mitigation up to a given node.
%
% 		Parameters
% 		----------
% 		m : ndarray or list
% 			array of mitigation
% 		node : int
% 			node for which average mitigation is to be calculated for
% 		period : int, optional
% 			the period the node is in
%
% 		Returns
% 		-------
% 		float
% 			average mitigation

            if nargin < 4 || isempty(period)
                period = NaN;
            end

            if period == 0
                r = 0;
            elseif isempty(period)
                period = obj.tree.get_period(node);
            else
                state = obj.tree.get_state(node);
                path = obj.tree.get_path(node);
                new_m = m(path(1:end-1)+1); % mitigation on the path until this node
                period_len = obj.tree.decision_times(2:period+1) - obj.tree.decision_times(1:period);
                bau_emissions = obj.bau.emission_by_decisions(1:period); % emission levels at each decision point
                total_emission = sum( (bau_emissions .* period_len)' ); % total emission: sum of emissions during each period
                tmp = bau_emissions .* period_len;

                ave_mitigation = sum( (new_m .* tmp)' );
                r = ave_mitigation / total_emission;
            end
        end


        function r = average_mitigation(obj, m, period)

%       Calculate the average mitigation for all nodes in a period.
%
% 		m : ndarray or list
% 			array of mitigation
% 		period : int
% 			period to calculate average mitigation for
%
% 		Returns
% 		-------
% 		ndarray
% 			average mitigations

            nodes = obj.tree.get_num_nodes_period(period); % number of nodes for a given period
            ave_mitigation = cons(zeros(1, nodes), m);
            for i = 1:nodes
                node = obj.tree.get_node(period, i-1);
                ave_mitigation(i) = obj.average_mitigation_node(m, node, period);
            end
            r = ave_mitigation;
        end


        function r = ghg_level_node(obj, m, node)

            forcingobj = Forcing();
            r = forcingobj.ghg_level_at_node(m, node, obj.tree, obj.bau, obj.subinterval_len);

        end


        function r = ghg_level_period(obj, m, period, nodes)

%       Calculate the GHG levels corresponding to the given mitigation.
% 		Need to provide either `period` or `nodes`.
%
% 		Parameters
% 		----------
% 		m : ndarray or list
% 			array of mitigation
% 		period : int, optional
% 			what period to calculate GHG levels for
% 		nodes : ndarray or list, optional
% 			the nodes to calculate GHG levels for
%
% 		Returns
% 		-------
% 		ndarray
% 			GHG levels

            if nargin == 2
                period = NaN;
                nodes = NaN;
            end

            if isempty(nodes) && ~isempty(period)
                [start_node, end_node] = obj.tree.get_nodes_in_period(period);
                if period >= obj.tree.num_periods
                    add = end_node - start_node + 1;
                    start_node = start_node + add;
                    end_node = end_node + add;
                end
                nodes = start_node:end_node;
            end

            if isempty(period) && isempty(nodes)
                error('Need to give function either nodes or the period.');
            end

            ghg_level = zeros(length(nodes));
            for i = 1:length(nodes)
                ghg_level(i) = obj.ghg_level_node(m, nodes(i));
            end
            r = ghg_level;
        end


        function r = ghg_level(obj, m, periods)

%       Calculate the GHG levels for more than one period.
%
% 		Parameters
% 		----------
% 		m : ndarray or list
% 			array of mitigation
% 		periods : int, optional
% 			number of periods to calculate GHG levels for
%
% 		Returns
% 		-------
% 		ndarray
% 			GHG levels

            if nargin == 2
                periods = NaN;
            end

            if isempty(periods)
                periods = obj.tree.num_periods-1;
            end

            if periods >= obj.tree.num_periods
                ghg_level = zeros(1, obj.tree.num_decision_nodes + obj.tree.num_final_states);
            else
                ghg_level = zeros(1, obj.tree.num_decision_nodes);
            end

            % for each period, find the right start and end node within this period and add ghg_level for this nodes
            for period = 1:periods+1
                [start_node, end_node] = obj.tree.get_nodes_in_period(period);
                if period >= obj.tree.num_periods
                    add = end_node - start_node + 1;
                    start_node = start_node + add;
                    end_node = end_node + add;
                end
                nodes = start_node : end_node;
                ghg_level(nodes+1) = obj.ghg_level_period(m, NaN, nodes);
            end

            r = ghg_level;
        end


        function [obj, r] = damage_function_node(obj, m, node)

%           Calculate the damage at any given node, based on mitigation actions in `m`.
            if isnan(obj.damage_coefs)
                obj = obj.damage_interpolation();
            end

            if isempty(obj.cum_forcings)
                obj = obj.forcing_init();
            end

            if node == 0
                r = 0.0;
            else

            period = obj.tree.get_period(node);

            if period == 0
                       p = 1;
            else
                       p = period;
            end

            forcingobj = Forcing();
            rarray = forcingobj.forcing_and_ghg_at_node(m, node, obj.tree, obj.bau, obj.subinterval_len, 'both');
            forcing = rarray(1);
            ghg_level = rarray(2);

            [obj, force_mitigation] = obj.forcing_based_mitigation(forcing, period);
            ghg_extension = 1.0 / (1 + exp(0.05 * (ghg_level-200)));

            rarray = obj.tree.reachable_end_state(node);
            worst_end_state = rarray(1);
            best_end_state = rarray(2);

            probs = obj.tree.final_state_probs(worst_end_state+1 : best_end_state+1);

            % new damage formulas
            if force_mitigation < obj.emit_pct(2)

                vec_to_sum = cons(zeros(1, best_end_state - worst_end_state + 1),m);
                for it = 1 : (best_end_state - worst_end_state + 1)
                    vec_to_sum(it) = probs(it) .* ( obj.damage_coefs(2, 2, period, worst_end_state + it) .* force_mitigation + ...
                        + obj.damage_coefs(2, 3, period, worst_end_state + it ) );
                end

                damage = sum(vec_to_sum');

%                 damage = sum( probs .* ( obj.damage_coefs(2, 2, period, worst_end_state+1:best_end_state+1) .* force_mitigation...
%                     + obj.damage_coefs(2, 3, period, worst_end_state+1:best_end_state+1) ) );

            elseif force_mitigation < obj.emit_pct(1)

                vec_to_sum = cons(zeros(1, best_end_state - worst_end_state + 1), m);
                for it = 1 : (best_end_state - worst_end_state + 1)
                    vec_to_sum(it) = probs(it) .* ( obj.damage_coefs(1, 1, period, worst_end_state + it) .* force_mitigation.^2 ...
                            + obj.damage_coefs(1, 2, period, worst_end_state + it) .* force_mitigation ...
                            + obj.damage_coefs(1, 3, p, worst_end_state + it));
                end

                damage = sum(vec_to_sum');
%                 damage = sum( prob * (obj.damage_coefs(1, 1, period, worst_end_state+1 : best_end_state+1) * force_mitigation^2 ...
%                             + obj.damage_coefs(1, 2, period, worst_end_state+1 : best_end_state+1) * force_mitigation ...
%                             + obj.damage_coefs(1, 3, p, worst_end_state+1 : best_end_state+1)));
            else
                damage = 0.0;
                i = 0;
                 for state = worst_end_state:best_end_state
                     if obj.d_rcomb(state+1, p, 1) > 1e-5
                         deriv = 2.0 * obj.damage_coefs(1, 1, p, state+1) * obj.emit_pct(1) ...
 							+ obj.damage_coefs(1, 2, p, state+1);

                         decay_scale = deriv / (obj.d_rcomb(state+1, p, 1) *log(0.5));

                         dist = force_mitigation - obj.emit_pct(1) + log(obj.d_rcomb(state+1, p, 1)) ...
                            / (log(0.5) * decay_scale);

                         damage = damage + probs(i+1) * (0.5^(decay_scale * dist) * exp(-(force_mitigation - obj.emit_pct(1)).^2/60.0));
                     end
                     i = i+1;
                end
            end
            %damage
            r = (damage ./ sum(probs')) + ghg_extension;
            end
        end


        function [obj, r] = damage_function(obj, m, period)

%       Calculate the damage for every node in a period, based on mitigation actions `m`.
%
% 		Parameters
% 		----------
% 		m : ndarray or list
% 			array of mitigation
% 		period : int
% 			period to calculate damages for
%
% 		Returns
% 		-------
% 		ndarray
% 			damages

            nodes = obj.tree.get_num_nodes_period(period);
            damages = cons(zeros(1, nodes),m);
            for i = 1:nodes
                node = obj.tree.get_node(period, i-1);
                [obj, damages(i)] = obj.damage_function_node(m, node);
            end
            r = damages;
        end

    end

end
