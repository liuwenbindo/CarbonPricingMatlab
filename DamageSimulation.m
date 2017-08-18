% ================================= %
%       DamageSimulation.m          %
%           Wenbin Liu              %
%       wl2581@columbia.edu         %
% ================================= %


% Simulation of damages for the EZ-Climate model.
% 
% The damage function simulation is a key input into the pricing engine. Damages are
% represented in arrays of dimension n x p, where n = num states and p = num periods.
% The arrays are created by Monte Carlo simulation. Each array specifies for each state
% and time period a damage coefficient.
% 
% Up to a point, the Monte Carlo follows Pindyck (2012) 'Uncertain Outcomes and Climate Change
% Policy':
% 
%     * There is a gamma distribution for temperature
%     * There is a gamma distribution for economic impact (conditional on temperature)
% 
% However, in addition, this program adds a probability of a tipping point (conditional on temperature).
% This probability is a decreasing function of the parameter `peak_temp`, conditional on a tipping
% point. Damage itself is a decreasing function of the parameter `disaster_tail`.
% 
% Parameters
% ----------
% tree : `TreeModel` object
%     tree structure used
% ghg_levels : ndarray or list
%     end GHG level for each path
% peak_temp : float
%     tipping point parameter
% disaster_tail : float
%     curvature of tipping point
% tip_on : bool
%     flag that turns tipping points on or off
% temp_map : int
%     mapping from GHG to temperature
% 
%         * 0: implies Pindyck displace gamma
%         * 1: implies Wagner-Weitzman normal
%         * 2: implies Roe-Baker
%         * 3: implies user-defined normal
%         * 4: implies user-defined gamma
% 
% temp_dist_params : ndarray or list
%     if temp_map is either 3 or 4, user needs to define the distribution parameters
% maxh : float
%     time paramter from Pindyck which indicates the time it takes for temp to get half
%     way to its max value for a given level of ghg
% cons_growth : float
%     yearly growth in consumption
% 
% Attributes
% ----------
% tree : `TreeModel` object
%     tree structure used
% ghg_levels : ndarray or list
%     end GHG level for each path
% peak_temp : float
%     tipping point parameter
% disaster_tail : float
%     curvature of tipping point
% tip_on : bool
%     flag that turns tipping points on or off
% temp_map : int
%     mapping from GHG to temperature
% temp_dist_params : ndarray or list
%     if temp_map is either 3 or 4, user needs to define the distribution parameters
% maxh : float
%     time paramter from Pindyck which indicates the time it takes for temp to get half
%     way to its max value for a given level of ghg
% cons_growth : float
%     yearly growth in consumption
% d : ndarray
%     simulated damages


classdef DamageSimulation
    properties
        tree
        ghg_levels
        peak_temp
        disaster_tail
        tip_on
        temp_map
        temp_dist_params
        maxh
        cons_growth
        d
        draws
    end
    
    methods
        %Constructor
        function obj = DamageSimulation (tree, ghg_levels, peak_temp, disaster_tail, tip_on...
                , temp_map, temp_dist_params, maxh, cons_growth)
            obj.tree = tree;
            obj.ghg_levels = ghg_levels;
            obj.peak_temp = peak_temp;
            obj.disaster_tail = disaster_tail;
            obj.tip_on = tip_on;
            obj.temp_map = temp_map;
            obj.temp_dist_params = temp_dist_params;
            obj.maxh = maxh;
            obj.cons_growth = cons_growth;
            obj.d = [];
        end
        
       
        function r = simulate(obj, draws, write_to_file)
%         Create damage function values in 'p-period' version of the Summers - Zeckhauser model.
% 
%         Parameters
%         ----------
%         draws : int
%             number of samples drawn in Monte Carlo simulation.
%         write_to_file : bool, optional
%             wheter to save simulated values
% 
%         Returns
%         -------
%         ndarray
%             3D-array of simulated damages # it should be 2D : self.tree.num_final_states, self.tree.num_periods
% 
%         Raises
%         ------
%         Error
%             If temp_map is not in the interval 0-4.

            if nargin < 3 || isempty(write_to_file)
                write_to_file = true;
            end

            dnum = length(obj.ghg_levels);
            obj.draws = draws;
            obj.peak_cons = exp(obj.cons_growth * obj.tree.decision_times(2:end));

            if obj.temp_map == 0
                temperature = obj.pindyck_simulation();
            elseif obj.temp_map == 1
                temperature = obj.ww_simulation();
            elseif obj.temp_map == 2
                temperature = obj.rb_simulation();
            elseif obj.temp_map == 3
                temperature = obj.normal_simulation();
            elseif obj.temp_map == 4
                temperature = obj.gamma_simulation();
            else
                Error('temp_map not in interval 0-4');
            end
            
            % Note: 
            % When testing, MATLAB keeps shooting the error:
            % 'Elements of cell array input A must be strings (row vectors
            % of class char).'
            % I found the solution of removing the
            % \MATLAB\R2015a\toolbox\local\pathdef.m 
            % from the paths, but I don't have the permission of doing so.
            % Can you help me test this code?
            
            job = batch(@obj.run_path, 1, {temperature}, 'Pool', dnum);
            wait(job);
            load(job);
            r = fetchOutputs(job); % Get results into a cell array
            obj.d = r{1};          % Save the result
            
            if write_to_file
                obj.write_to_file();
            end
            
            r = obj.d;
        end
    end

    
    methods (Access = protected)
        function write_to_file(obj) 
            filename = 'simulated_damages';
            delimeter = ',';
            %first = obj.d(1,:); % d[0] in python code stands for the first row
            %rest = obj.d(2:end,:);
            d_write = obj.d;
            
            fopen(filename,'wt');
            %dlmwrite(filename, first',delimeter);
            %dlmwrite(filename, rest', delimeter);
            dlmwrite(filename, d_write', delimeter);
        end
        
        %dimension here is an array, like [1,5]
        function r = gamma_array(shape, rate, dimension)
            r = gamrnd(shape, 1/rate, dimension);
        end
        
        function r = normal_array(mean, stdev, dimension)
            r = normrnd(mean, stdev, dimension);
        end
        
        function r = uniform_array(dimension)
            r = unifrnd(0,1,dimension);
        end
        
        function r = sort_array(obj, array)
            col = obj.tree.num_periods;    
            [value, order] = sort(array(:,col));
            r = array(order,:);
        end
        
        function r = normal_distribution(obj)
           if (length(obj.temp_dist_params) == 2)
               ave = obj.temp_dist_params(:,1);
               std = obj.temp_dist_params(:,2);
               n = length(ave);
               for i = 1:n
                   temperature(i,:) = obj.normal_array(ave(i), std(i), obj.draws); %obj.draws is array like [1, 50]                  
               end           
               r = exp(temperature);
           end
        end
         
        function r = gamma_simulation(obj)
            if (length(obj.temp_dist_params) == 3)
               k = obj.temp_dist_params(:,1);
               theta = obj.temp_dist_params(:,2);
               displace = obj.temp_dist_params(:,3);
               n = length(k);
               for i = 1:n
                   gamma_temp(i,:) = obj.gamma_array(k(i), theta(i), obj.draws); %obj.draws is array like [1, 50]
                   r(1,:) = gamma_temp(i,:) + displace(i);
               end         
           end
        end
        
        function r = generate_parameter(mean, stdev, times)
            para_list = [];
            for i = 1:times
                temp = normrnd(mean, stdev);
                para_list = [para_list, temp];
            end
            r = para_list;            
        end
        
        function r = pindyck_impact_simulation(obj)
            pindyck_impact_k = 4.5;
            pindyck_impact_theta = 21341.0;
            pindyck_impact_displace = -0.0000746;
            r = obj.gamma_array(pindyck_impact_k, pindyck_impact_theta, obj.draws) + ...
                 pindyck_impact_displace;           
        end
        
        function r = ww_simulation(obj)
            ww_temp_ave = [0.573, 1.148, 1.563];
            ww_temp_stddev = [0.462, 0.441, 0.432];
            for i = 1:3
                temperature(i,:) = obj.normal_array(ww_temp_ave(i), ww_temp_stddev(i), obj.draws);
            end      
            r = exp(temperature);
        end
        
        function r = rb_simulation(obj)
            rb_fbar = [0.75233, 0.844652, 0.858332];
            rb_sigf = [0.049921, 0.033055, 0.042408];
            rb_theta = [2.304627, 3.333599, 2.356967];
    
            for i = 1:3
                temperature(i,:) = obj.normal_array(rb_fbar(i), rb_sigf(i), obj.draws);
                item_to_compare(i,:) = (1 / ( 1 - temperature(i,:))) - rb_theta(i);
            end                       
            r = max(0, item_to_compare);
        end
        
        function r = disaster_simulation(obj)
            dim_unif = [obj.draws, obj.tree.num_periods];
            r = obj.uniform_array(dim_unif);
        end
        
        function r = disaster_cons_simulation(obj)
            r = obj.gamma_array(1.0, obj.disaster_tail, obj.draws); %obj.disaster_tail
        end
         
        function r = interpolation_of_temp(obj, temperature)
            for i = 1:size(temperature,1)
                r(i,:) = temperature(i,:) * 2.0 * (1 -  0.5^(obj.tree.decision_times(i+1)/obj.maxh));
            end
        end
        
        function r = economic_impact_of_temp(obj, temperature)
            impact = obj.pindyck_impact_simulation();
            for i = 1:length(impact)
                term1(i,:) = -2.0 * impact(i,:) .* temperature(i,:) * obj.maxh / log(0.5); 
                term2(i,:) = obj.cons_growth - 2.0 .* impact(i,:) .* temperature(i,:) * obj.tree.decision_times(i+1);
                term3(i,:) = 2.0 * impact(i,:) .* temperature(i,:) * 0.5^(obj.tree.decision_times(i+1)/obj.maxh)...
                             * obj.maxh / log(0.5);
            end
            r = exp(term1 + term2 + term3);
        end
        
        function r = tipping_point_update(obj, tmp, consump, peak_temp_interval)
             if nargin < 4 || isempty(peak_temp_interval)
                peak_temp_interval = 30.0;
             end
            
             draws = size(tmp, 1);
             disaster = obj.disaster_simulation();
             disaster_cons = obj.disaster_cons_simulation();
             period_lengths = obj.tree.decision_times(2:end) - obj.tree.decision_times(1:end-1);
             tmp_scale = max(obj.peak_temp,tmp);
             ave_prob_of_survival = 1 - (tmp/tmp_scale).^2;
             prob_of_survival = (ave_prob_of_survival).^(period_lengths / peak_temp_interval);
             
             res = prob_of_survival < disaster;
             [rows, cols] = find(res); %Find the rows and cols of nonzero elements
             row = unique(rows);
             counts = histc(rows, row); %Count of unique values
             cum_counts = cumsum(counts);
             insert_idx = zeros(1, length(cum_counts)); %insert_idx has same length with cum_counts
             insert_idx(1) = 1;
             for i = 1:length(cum_counts)-1 
                 insert_idx(i+1) = cum_counts(i) + 1;
             end
             
             cols_select = cols(insert_idx); %get the #col of first non-zero element in every row
             
             first_occurance = {};
             for i = 1:length(cols_select)
                 first_occurance{i} = [row(i), cols_select(i)];
             end
             for i = 1:length(first_occurance)
                 tmp_array = first_occurance{i};
                 for j = tmp_array(1):length(tmp_array)
                 consump(tmp_array(0),j) = consump(tmp_array(0),j) * exp(-disaster_cons(tmp_array(0)));
                 end
             end
             r = consump;    
        end
        
        function r = run_path(obj, temperature)
%         Calculate the distribution of damage for specific GHG-path. Implementation of
%         the temperature and economic impacts from Pindyck [2012] page 6.
%         
%         Remark
%         -------------
%         final states given periods can give us a specific state in that period since a child only have one parent

        % d_temp is d in python code
        d_temp = zeros(obj.tree.num_final_states, obj.tree.num_periods);
        tmp = obj.interpolation_of_temp(temperature);
        consump = obj.economic_impact_of_temp(temperature);
        peak_cons = exp(obj.cons_growth * obj.tree.decision_times(2:end));
        
        % adding tipping points
        if obj.tip_on
            consump = obj.tipping_point_update(tmp, consump);
        end
        
        % sort based on outcome of simulation
        consump = obj.sort_array(consump);
        damage = 1.0 - (consump / peak_cons);
        weights = obj.tree.final_states_prob*(obj.draws);
        weights = floor(cumsum(weights));
        
        d_temp(1,:) = mean(damage(1:weights(1),:)); %calculating mean value by column
        for n = 2:obj.tree.num_final_states
            d_temp(n,:) = max(0.0,  mean(damage(weights(n-1):weights(n),:)));
        end
        r = d_temp;
        end       
    end
end