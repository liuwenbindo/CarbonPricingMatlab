classdef EZUtility
    
   properties
       tree
       damage
       cost
       period_len
       decision_times
       cons_growth
       growth_term
       r
       a
       b
       potential_cons
       add_penalty_cost
       max_penalty
       penalty_scale
   end
   
   
   methods
       
       % Constructor
       function obj = EZUtility (tree, damage, cost, period_len, eis, ra, time_pref, add_penalty_cost,...
               max_penalty, penalty_scale)
           % Set default values of input arguments
           if nargin == 4
                eis = 0.9;
                ra = 7.0;
                time_pref = 0.005;
                add_penalty_cost = false;
                max_penalty = 0.0;
                penalty_scale = 1.0;
           end
           
           obj.tree = tree;
           obj.damage = damage;
           obj.cost = cost;
           obj.period_len = period_len;
           obj.decision_times = tree.decision_times;
           obj.cons_growth = damage.cons_growth;
           obj.growth_term = 1.0 + obj.cons_growth;
           obj.r = 1.0 - 1.0/eis;
           obj.a = 1.0 - ra;
           obj.b = (1.0 - time_pref)^period_len;
           obj.potential_cons =  (ones(size(obj.tree.decision_times)) + obj.cons_growth) .^ (obj.decision_times);
           obj.add_penalty_cost = add_penalty_cost;
           obj.max_penalty = max_penalty;
           obj.penalty_scale = penalty_scale;          
       end
       
       
       function [utility_tree, cons_tree, cost_tree] = end_period_utility(obj, m, utility_tree, cons_tree, cost_tree)
           % Calculate the terminal utility
           period_ave_mitigation = obj.damage.average_mitigation(m, obj.tree.num_periods);
           % Average_mitigation of all the possibilities in the final period (an array containing all the possible ave_mitigation)
           [obj.damage, period_damage] = obj.damage.damage_function(m, obj.tree.num_periods);
           damage_nodes = obj.tree.get_nodes_in_period(obj.tree.num_periods);
           
           period_mitigation = m(damage_nodes(1)+1:damage_nodes(2)+1);
           period_cost = obj.cost.cost(obj.tree.num_periods, period_mitigation, period_ave_mitigation);
           continuation = (1.0 / (1.0 - obj.b * (obj.growth_term^obj.r)))^(1.0 / obj.r);
           
           % Set value for cost_tree.last_period
           cost_tree = cost_tree.set_value(cost_tree.last_period, period_cost');
           
           % Set value for cons_tree.last_period
           period_consumption = obj.potential_cons(end) * (1.0 - period_damage);
           period_consumption(period_consumption <= 0) = 1e-18;           
           cons_tree = cons_tree.set_value(cons_tree.last_period, period_consumption');            
           
           % Set value for utility_tree.last_period
           utilivalue = (1.0 - obj.b)^(1.0 / obj.r) * cons_tree.last * continuation;          
           utility_tree = utility_tree.set_value(utility_tree.last_period, utilivalue);
       end
       
       
       function [ce_tree, mu_tree_0, mu_tree_1] = end_period_marginal_utility(obj, mu_tree_0, mu_tree_1, ce_tree, utility_tree, cons_tree)
           
           % Calculate the terminal marginal utility
           
           % the utility of the final state can be seperated into two parts: the part that is generated by consumption on T and the part is generated in the future.
           % the U_T can be wrote in to U_T = [(1-beta)c_T^rau+(1-beta)/(1-beta*(1+g)^rau*c_T*rau]^(1/rau)
           % the margin is future - now which is calculated by following codes.
           
           ce_term = utility_tree.last ^ obj.r - (1.0 - obj.b) * cons_tree.last^obj.r;
           ce_tree = ce_tree.set_value(ce_tree.last_period, ce_term);
           
           mu_0_last = (1.0 - obj.b) * (utility_tree(utility_tree.last_period - obj.period_len) / cons_tree.last)^(1.0 - obj.r);
           mu_tree_0 = mu_tree_0.set_value(mu_tree_0.last_period, mu_0_last);
           mu_0 = obj.mu_0(cons_tree(cons_tree.last_period - obj.period_len), ce_tree(ce_tree.last_period - obj.period_len));
           mu_tree_0 = mu_tree_0.set_value(mu_tree_0.last_period - obj.period_len, mu_0);
           
           next_term = obj.b * (1.0 - obj.b) / (1.0 - obj.b * obj.growth_term^obj.r);
           mu_1 = utility_tree(utility_tree.last_period - obj.period_len)^(1 - obj.r) * next_term...
                  * cons_tree.last^(obj.r - 1.0);
           mu_tree_1 = mu_tree_1.set_value(mu_tree_1.last_period - obj.period_len, mu_1);
       end
       
       
       function r = certain_equivalence(obj, period, damage_period, utility_tree)
           
%       Caclulate certainty equivalence utility. If we are between decision nodes, i.e. no branching,
% 		then certainty equivalent utility at time period depends only on the utility next period
% 		given information known today. Otherwise the certainty equivalent utility is the ability
% 		weighted sum of next period utility over the partition reachable
% 		from the state.

           if utility_tree.is_information_period(period)
               damage_nodes = obj.tree.get_nodes_in_period(damage_period + 1);
               probs = obj.tree.all_node_probs(damage_nodes(1)+1:damage_nodes(2)+1);
               even_probs = probs(1:2:end); % Odd # elements in MATLAB, same as even # of elements in Python.
               odd_probs = probs(2:2:end); % Same as above.
               
               utility_tree_getnexttmp = utility_tree.get_next_period_array(period);
               even_util = (utility_tree_getnexttmp(1:2:end).^(obj.a))' .* even_probs;
               
               odd_util = (utility_tree_getnexttmp(2:2:end).^(obj.a))' .* odd_probs;
               
               ave_util = (even_util + odd_util) ./ (even_probs + odd_probs);            
               cert_equiv = ave_util.^(1.0/obj.a);               
           else
               % no branching implies certainty equivalent utility at time period depends only on
               % the utility next period given information known today
               cert_equiv = (utility_tree.get_next_period_array(period))';
           end
           
           if size(cert_equiv,1) == 1
               r = cert_equiv;
           else
               r = cert_equiv';
           end         
       end
       
       
       function [utility_tree, r1_u, r2_period] = utility_generator(obj, m, utility_tree, cons_tree, cost_tree, ce_tree, cons_adj)
           
           % Generator for calculating utility for each utility period
           % besides the terminal utility.
           
           % there are two kinds of periods: make decision/not make decision.
           
           if nargin < 7 || isempty(cons_adj)
               cons_adj = 0.0;
           end
           
           r1_u = {};
           r2_period = {};
           k = 1; %iterator
           
           periods = fliplr(utility_tree.periods);                     
           
           fprintf('Setting Utility Value...');
           for i= 2:length(periods)              
                period = periods(i);
                damage_period = utility_tree.between_decision_times(period);
                cert_equiv = obj.certain_equivalence(period, damage_period, utility_tree);              
                
                if utility_tree.is_decision_period( period + obj.period_len )
                    
                    damage_nodes = obj.tree.get_nodes_in_period(damage_period);
                    period_mitigation = m(damage_nodes(1)+1:damage_nodes(2)+1);                  
                    period_ave_mitigation = obj.damage.average_mitigation(m, damage_period);                    
                    period_cost = obj.cost.cost(damage_period, period_mitigation, period_ave_mitigation);
                    [obj.damage, period_damage] = obj.damage.damage_function(m, damage_period);
                    
                    indexbelow = int32(cost_tree.index_below(period + obj.period_len));
                    cost_tree = cost_tree.set_value(indexbelow, period_cost');  
                    
                end
                
                period_consumption = obj.potential_cons(damage_period + 1) ... %NOTICE: damage_period + 1
                    .* (1.0 - period_damage) .* (1.0 - period_cost);               
                period_consumption(period_consumption <= 0) = 1e-18;
                
                % if not a decision time
                if ~utility_tree.is_decision_period( period )
                    
                    next_consumption = cons_tree.get_next_period_array( period );
                    segment = period - utility_tree.decision_times(damage_period + 1); %NOTICE: damage_period + 1
                    interval = segment + utility_tree.subinterval_len;                  
                    
                    % if the next period is a decision period
                    if utility_tree.is_decision_period(period + obj.period_len)
                        if period < utility_tree.decision_times(end - 1)                             
                            next_cost = cost_tree.tree( period + obj.period_len);                   
                            next_consumption = next_consumption' .* (1.0 - repelem(period_cost, 2))./(1.0 - next_cost');
                            next_consumption(next_consumption <= 0) = 1e-18;
                        end
                    end
                    
                    if size(next_consumption,1) == 1
                    else
                        next_consumption = next_consumption';                       
                    end
                    
                    % if the information is not gained i55  n next period, the consumption next period is the same to now
                    if period <  utility_tree.decision_times(end - 1)
                        
                        temp_consumption = next_consumption ./ repelem(period_consumption, 2);
                        period_consumption = sign(temp_consumption) .* (abs(temp_consumption).^(segment./double(interval)))...
                            .* repelem(period_consumption, 2);
                    else                        
                        temp_consumption = next_consumption./(period_consumption); % for the final stage, we know all the info, and thus the period_consumption ???
                        period_consumption = sign(temp_consumption) .* (abs(temp_consumption).^(segment./double(interval)))...
                            .* period_consumption;
                    end
                end
                
                if period == 0
                    period_consumption = period_consumption + cons_adj;
                end
                
                if size(period_consumption, 1) == 1                 
                else
                    period_consumption = period_consumption';
                end
                               
                ce_term = obj.b * cert_equiv.^(obj.r);               
                
                ce_tree.set_value(period, ce_term');                
                cons_tree.set_value(period, period_consumption');
                              
                u = ((1.0 - obj.b) * (period_consumption).^(obj.r) + ce_term).^(1.0/obj.r);
                % yield u, period
                r1_u{k} = u;
                r2_period{k} = period;
                utility_tree = utility_tree.set_value(period, u');
                k = k+1;
                
            end
       end
       
       
       function [r0, r1_utility, r2_cons, r3_cost, r4_ce] = utility(obj, m, return_trees)
           
%       Calculating utility for the specific mitigation decisions `m`.
% 
% 		Parameters
% 		----------
% 		m : ndarray or list
% 			array of mitigations
% 		return_trees : bool
% 			True if methid should return trees calculated in producing the utility
% 
% 		Returns
% 		-------
% 		ndarray or tuple
% 			tuple of `BaseStorageTree` if return_trees else ndarray with utility at period 0
% 
% 		Examples:
% 		---------
% 		Assuming we have declared a EZUtility object as 'ezu' and have a mitigation array 'm'
% 
% 		>>> ezu.utility(m)
% 		array([ 9.83391921])
% 		>>> utility_tree, cons_tree, cost_tree, ce_tree = ezu.utility(m, return_trees=True)

           
           if nargin < 3 || imempty(return_trees)
               return_trees = false;
           end
           
           utility_tree = BigStorageTree(obj.period_len, obj.decision_times);
           cons_tree = BigStorageTree(obj.period_len, obj.decision_times);
           ce_tree = BigStorageTree(obj.period_len, obj.decision_times);
           cost_tree = SmallStorageTree(obj.decision_times);
           
           [utility_tree, cons_tree, cost_tree] = obj.end_period_utility(m, utility_tree, cons_tree, cost_tree);
           
           [utility_tree, ~, ~] = obj.utility_generator(m, utility_tree, cons_tree, cost_tree, ce_tree);
             
%            for k = 1:length(u)
%                utility_tree = utility_tree.set_value(period{k}, (u{k})');
%            end         
           
           if return_trees
               r1_utility = utilty_tree;
               r2_cons = cons_tree;
               r3_cost = cost_tree;
               r4_ce = ce_tree;
               r0 = NaN;
           else
               r0 = utility_tree.tree(0);
               r1_utility = NaN;
               r2_cons = NaN;
               r3_cost = NaN;
               r4_ce = NaN;
           end                              
       end
       
       
       function [r0, r1_utility, r2_cons, r3_cost, r4_ce] = adjusted_utility(obj, m, period_cons_eps, node_cons_eps, final_cons_eps, first_period_consadj, return_trees)
           
%       Calculating adjusted utility for sensitivity analysis. Used e.g. to find zero-coupon bond price.
% 		Values in parameters are used to adjusted the utility in different ways.
% 
% 		Parameters
% 		----------
% 		m : ndarray
% 			array of mitigations
% 		period_cons_eps : ndarray, optional
% 			array of increases in consumption per period
% 		node_cons_eps : `SmallStorageTree`, optional
% 			increases in consumption per node
% 		final_cons_eps : float, optional
% 			value to increase the final utilities by
% 		first_period_consadj : float, optional
% 			value to increase consumption at period 0 by
% 		return_trees : bool, optional
% 			True if method should return trees calculculated in producing the utility
% 
% 		Returns
% 		-------
% 		ndarray or tuple
% 			tuple of `BaseStorageTree` if return_trees else ndarray with utility at period 0
% 
% 		Examples
% 		---------
% 		Assuming we have declared a EZUtility object as 'ezu' and have a mitigation array 'm'
% 
% 		>>> ezu.adjusted_utility(m, final_cons_eps=0.1)
% 		array([ 9.83424045])
% 		>>> utility_tree, cons_tree, cost_tree, ce_tree = ezu.adjusted_utility(m, final_cons_eps=0.1, return_trees=True)
% 
% 		>>> arr = np.zeros(int(ezu.decision_times[-1]/ezu.period_len) + 1)
% 		>>> arr[-1] = 0.1
% 		>>> ezu.adjusted_utility(m, period_cons_eps=arr)
% 		array([ 9.83424045])
% 
% 		>>> bst = BigStorageTree(5.0, [0, 15, 45, 85, 185, 285, 385])
% 		>>> bst.set_value(bst.last_period, np.repeat(0.01, len(bst.last)))
% 		>>> ezu.adjusted_utility(m, node_cons_eps=bst)
% 		array([ 9.83391921])
% 
% 		The last example differs from the rest in that the last values of the `node_cons_eps` will never be
% 		used. Hence if you want to update the last period consumption, use one of these two methods.
% 
% 		>>> ezu.adjusted_utility(m, first_period_consadj=0.01)
% 		array([ 9.84518772])

           if nargin < 7
               if isempty(period_cons_eps)
                 period_cons_eps = NaN;
               end
               if isempty(node_cons_eps)
                 node_cons_eps = NaN;
               end
               if isempty(final_cons_eps)
                 final_cons_eps = 0.0;
               end
               if isempty(first_period_consadj)
                 first_period_consadj = 0.0;
               end
               if isempty(return_trees)
                 return_trees = false;
               end
           end
           
           utility_tree = BigStorageTree(obj.period_len, obj.decision_times);
           cons_tree = BigStorageTree(obj.period_len, obj.decision_times);
           ce_tree = BigStorageTree(obj.period_len, obj.decision_times);
           cost_tree = SmallStorageTree(obj.decision_times);
           
           periods = fliplr(utility_tree.periods);
           
           if isnan(period_cons_eps)
               period_cons_eps = zeros(1, length(periods));               
           end
           if isnan(node_cons_eps)
               node_cons_eps = BigStorageTree(obj.period_len, obj.decision_times);
           end
           
           [utility_tree, cons_tree, cost_tree] = obj.end_period_utility(m, utility_tree, cons_tree, cost_tree);          
           %[utility_tree, u, period] = obj.utility_generator(m, utility_tree, cons_tree, cost_tree, ce_tree, first_period_consadj);
           
           
           
           %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    
           periods = fliplr(utility_tree.periods);                     
           
           fprintf('Setting Utility Value...');
           
           j = length(utility_tree.tree) - 2;
           
           for i = 2:length(periods) 
               
                period = periods(i);
                damage_period = utility_tree.between_decision_times(period);
                cert_equiv = obj.certain_equivalence(period, damage_period, utility_tree);              
                
                if utility_tree.is_decision_period( period + obj.period_len )
                    
                    damage_nodes = obj.tree.get_nodes_in_period(damage_period);
                    period_mitigation = m(damage_nodes(1)+1:damage_nodes(2)+1);                  
                    period_ave_mitigation = obj.damage.average_mitigation(m, damage_period);                    
                    period_cost = obj.cost.cost(damage_period, period_mitigation, period_ave_mitigation);
                    [obj.damage, period_damage] = obj.damage.damage_function(m, damage_period);
                    
                    indexbelow = int32(cost_tree.index_below(period + obj.period_len));
                    cost_tree = cost_tree.set_value(indexbelow, period_cost');  
                    
                end
                
                period_consumption = obj.potential_cons(damage_period + 1) ... %NOTICE: damage_period + 1
                    .* (1.0 - period_damage) .* (1.0 - period_cost);               
                period_consumption(period_consumption <= 0) = 1e-18;
                
                % if not a decision time
                if ~utility_tree.is_decision_period( period )
                    
                    next_consumption = cons_tree.get_next_period_array( period );
                    segment = period - utility_tree.decision_times(damage_period + 1); %NOTICE: damage_period + 1
                    interval = segment + utility_tree.subinterval_len;                  
                    
                    % if the next period is a decision period
                    if utility_tree.is_decision_period(period + obj.period_len)
                        if period < utility_tree.decision_times(end - 1)                             
                            next_cost = cost_tree.tree( period + obj.period_len);                   
                            next_consumption = next_consumption' .* (1.0 - repelem(period_cost, 2))./(1.0 - next_cost');
                            next_consumption(next_consumption <= 0) = 1e-18;
                        end
                    end
                    
                    if size(next_consumption,1) == 1
                    else
                        next_consumption = next_consumption';                       
                    end
                    
                    % if the information is not gained i55  n next period, the consumption next period is the same to now
                    if period <  utility_tree.decision_times(end - 1)
                        
                        temp_consumption = next_consumption ./ repelem(period_consumption, 2);
                        period_consumption = sign(temp_consumption) .* (abs(temp_consumption).^(segment./double(interval)))...
                            .* repelem(period_consumption, 2);
                    else                        
                        temp_consumption = next_consumption./(period_consumption); % for the final stage, we know all the info, and thus the period_consumption ???
                        period_consumption = sign(temp_consumption) .* (abs(temp_consumption).^(segment./double(interval)))...
                            .* period_consumption;
                    end
                end
                
                if period == 0
                    period_consumption = period_consumption + first_period_consadj;
                end
                
                if size(period_consumption, 1) == 1                 
                else
                    period_consumption = period_consumption';
                end
                               
                ce_term = obj.b * cert_equiv.^(obj.r);               
                
                ce_tree.set_value(period, ce_term');                
                cons_tree.set_value(period, period_consumption');
                              
                u = ((1.0 - obj.b) * (period_consumption).^(obj.r) + ce_term).^(1.0/obj.r);
                % yield u, period
                % utility_tree = utility_tree.set_value(period, u');
                
               if period == periods(2)                  
                   mu_0 = (1.0 - obj.b) * (u./cons_tree.tree(period)').^(1.0 - obj.r);                  
                   next_term = obj.b * (1.0 - obj.b) / (1.0 - obj.b * obj.growth_term^(obj.r));
                   mu_1 = (u.^(1.0 - obj.r)) * next_term * (cons_tree.last.^(obj.r - 1.0));                  
                   u = u + ((final_cons_eps + period_cons_eps(end) + node_cons_eps.last) * mu_1)';                 
                   u = u + (period_cons_eps(j+1) + node_cons_eps.tree(period))' .* mu_0;                   
                   utility_tree = utility_tree.set_value(period, u');
               else
                   [mu_0, mu_1, mu_2] = obj.period_marginal_utility(mu_0, mu_1, m, period, utility_tree, cons_tree, ce_tree);                                   
                   u = u + ((period_cons_eps(j+1) + node_cons_eps.tree(period)) .* mu_0)';              
                   utility_tree  = utility_tree.set_value(period, u');            
               end
               j = j - 1;   
               
            end
           
           
           %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
           
          
          if return_trees
               r1_utility = utilty_tree;
               r2_cons = cons_tree;
               r3_cost = cost_tree;
               r4_ce = ce_tree;
               r0 = NaN;
           else
               r0 = utility_tree.tree(0);
               r1_utility = NaN;
               r2_cons = NaN;
               r3_cost = NaN;
               r4_ce = NaN;
           end                
       end
       
       
       function r = mu_0(obj, cons, ce_term)
		% Marginal utility with respect to consumption function.
		t1 = (1.0 - obj.b)*cons.^(obj.r-1.0);
		t2 = (ce_term - (obj.b-1.0) * cons.^obj.r).^((1.0/obj.r)-1.0);
		r =  t1 .* t2;
       end
       
       
       function rvalue = mu_1(obj, cons, prob, cons_1, cons_2, ce_1, ce_2, do_print)
		% Marginal utility with respect to consumption next period.
        if nargin < 8 || isempty(do_print)
            do_print = false;
        end
        
		t1 = (1.0-obj.b) * obj.b * prob .* (cons_1.^(obj.r-1.0))';
        
		t2 = (ce_1 - (obj.b-1.0) * cons_1.^obj.r ).^((obj.a/obj.r)-1);
        % turn into row vector
        t2 = t2';
        
		t3 = (prob' .* (ce_1 - (obj.b .* (cons_1.^obj.r)) + cons_1.^obj.r).^(obj.a/obj.r) ...
			 + (1.0-prob)' .* (ce_2 - (obj.b-1.0) * cons_2.^obj.r).^(obj.a/obj.r)).^((obj.r/obj.a)-1.0);
        t3 = t3';
         
		t4 = prob' .* (ce_1 - obj.b*(cons_1.^obj.r) + cons_1.^obj.r).^(obj.a/obj.r) ...
			 + (1.0-prob)' .* (ce_2 - obj.b * (cons_2.^obj.r) + cons_2.^obj.r).^(obj.a/obj.r);
        % turn into row vector
        t4 = t4';
         
		t5 = (obj.b * t4.^(obj.r/obj.a) - ((obj.b-1.0) * cons.^obj.r)' ).^((1.0/obj.r)-1.0);

		rvalue = t1 .* t2 .* t3 .* t5;
       end
       
       
       function rvalue = mu_2(obj, cons, prev_cons, ce_term)
		% Marginal utility with respect to last period consumption.
		t1 = (1.0-obj.b) * obj.b * prev_cons.^(obj.r-1.0);
		t2 = ((1.0 - obj.b) * cons.^obj.r - (obj.b - 1.0) * obj.b ...
		     * prev_cons.^obj.r + obj.b * ce_term).^((1.0/obj.r)-1.0);
		rvalue = t1 .* t2;
       end

       
       function [r1, r2, r3] = period_marginal_utility(obj, prev_mu_0, prev_mu_1, m, period, utility_tree, cons_tree, ce_tree)
           % Marginal utility for each node in a period.
           damage_period = utility_tree.between_decision_times(period);
           mu_0 = obj.mu_0(cons_tree.tree(period), ce_tree.tree(period)); 

           prev_ce = ce_tree.get_next_period_array(period);
           prev_cons = cons_tree.get_next_period_array(period);
           if utility_tree.is_information_period(period)
               
                probs = obj.tree.get_probs_in_period(damage_period+1);
                
                %iterator for up_prob
                up_it = 1;
                for i = 1:2:length(probs)                    
                    up_prob(1, up_it) = probs(i) / (probs(i) + probs(i+1));
                    up_it = up_it + 1;
                end          
                down_prob = 1.0 - up_prob;
                
%               up_cons = prev_cons[::2]
%               down_cons = prev_cons[1::2]
%               up_ce = prev_ce[::2]
%               down_ce = prev_ce[1::2]
                up_cons = prev_cons(1:2:end);
                down_cons = prev_cons(2:2:end);
                up_ce = prev_ce(1:2:end);
                down_ce = prev_ce(2:2:end);
                
                % cons_tree(period + 1): it's a dictionary (map)
                mu_1 = obj.mu_1(cons_tree.tree(period), up_prob, up_cons, down_cons, up_ce, down_ce);
                mu_2 = obj.mu_1(cons_tree.tree(period), down_prob, down_cons, up_cons, down_ce, up_ce);
                
                r1 = mu_0;
                r2 = mu_1;
                r3 = mu_2;
           else
               mu_1 = obj.mu_2(cons_tree.tree(period), prev_cons, prev_ce);
               
               r1 = mu_0;
               r2 = mu_1;
               r3 = NaN;
           end  
       end
       
    
       function [r1, r2, r3] = marginal_utlity(obj, m, utility_tree, cons_tree, cost_tree, ce_tree)
           
%      	Calculating marginal utility for sensitivity analysis, e.g. in the SSC decomposition.
% 
% 		Parameters
% 		----------
% 		m : ndarray
% 			array of mitigations
% 		utility_tree : `BigStorageTree` object
% 			utility values from using mitigation `m`
% 		cons_tree : `BigStorageTree` object
% 			consumption values from using mitigation `m`
% 		cost_tree : `SmallStorageTree` object
% 			cost values from using mitigation `m`
% 		ce_tree : `BigStorageTree` object
% 			certain equivalence values from using mitigation `m`
% 
% 		Returns
% 		-------
% 		tuple
% 			marginal utility tree
% 
% 		Examples
% 		--------
% 		Assuming we have declared a EZUtility object as 'ezu' and have a mitigation array 'm'.
% 		>>>
% 		>>> utility_tree, cons_tree, cost_tree, ce_tree = ezu.utility(m, return_trees=True)
% 		>>> mu_0_tree, mu_1_tree, mu_2_tree = ezu.marginal_utility(m, utility_tree, cons_tree, cost_tree, ce_tree)
% 		>>> mu_0_tree[0] # value at period 0
% 		array([ 0.33001256])
% 		>>> mu_1_tree[0] # value at period 0
% 		array([ 0.15691619])
% 		>>> mu_2_tree[0] # value at period 0
% 		array([ 0.13948175])
                     
        mu_tree_0 = BigStorageTree(obj.period_len, obj.decision_times);
        mu_tree_1 = BigStorageTree(obj.period_len, obj.decision_times);
        mu_tree_2 = SmallStorageTree(obj.decision_times);

        obj.end_period_marginal_utility(mu_tree_0, mu_tree_1, ce_tree, utility_tree, cons_tree);
        periods = fliplr(utility_tree.periods);

            for period = periods(3):periods(end)
                [mu_0, mu_1, mu_2] = obj.period_marginal_utility(mu_tree_0.get_next_period_array(period), ...
                    mu_tree_1.get_next_period_array(period), m, period, utility_tree, cons_tree, ce_tree);
                mu_tree_0.set_value(period, mu_0);
                mu_tree_1.set_value(period, mu_1);

                if ~isnan(mu_2)
                    mu_tree_2.set_value(period, mu_2);
                end

                r1 = mu_tree_0;
                r2 = mu_tree_1;
                r3 = mu_tree_2;
            end
       end
       
        
       function r = partial_grad(obj, m, i, delta)
           
%       Calculate the ith element of the gradient vector.
% 
% 		Parameters
% 		----------
% 		m : ndarray
% 			array of mitigations
% 		i : int
% 			node to calculate partial grad for
% 
% 		Returns
% 		-------
% 		float
% 			gradient element
           
           if nargin < 4 || isempty(delta)
               delta = 1e-8;
           end
           
           % m_copy = m.copy(); 
           % Why do we need to shallow copy instead of using m directly?
           
           m_copy = m;
           m_copy(i) = m_copy(i) - delta;
           minus_utility = obj.utility(m_copy);
           m_copy(i) = m_copy(i) + 2*delta;
           plus_utility = obj.utility(m_copy);
           grad = (plus_utility-minus_utility) / (2*delta);
           r = grad;
       end
       
    end
       
end