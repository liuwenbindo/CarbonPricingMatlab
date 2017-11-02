% ================================= %
%            DLWCost.m              %
%           Wenbin Liu              %
%       wl2581@columbia.edu         %
% ================================= %

classdef DLWCost < Cost
    
%   Class to evaluate the cost curve for the EZ-Climate model.
% 
% 	Parameters
% 	----------
% 	tree : `TreeModel` object
% 		tree structure used
% 	emit_at_0 : float
% 		initial GHG emission level
% 	g : float --> const of k = gx^a
% 		initial scale of the cost function
% 	a : float --> alpha
% 		curvature of the cost function
% 	join_price : float -->tau_*
% 		price at which the cost curve is extended
% 	max_price : float --> tau_tilda
% 		price at which carbon dioxide can be removed from atmosphere in unlimited scale
% 	tech_const : float  --> alpha_0
% 		determines the degree of exogenous technological improvement over time. A number
% 			of 1.0 implies 1 percent per year lower cost
% 	tech_scale : float  --> alpha_1
% 		determines the sensitivity of technological change to previous mitigation
% 	cons_at_0 : float   --> c_bar
% 		initial consumption. Default $30460bn based on US 2010 values.
% 
% 	Attributes cbs: cost as a fraction of baseline consumption
% 	----------
% 	tree : `TreeModel` object
% 		tree structure used
% 	g : float
% 		initial scale of the cost function
% 	a : float
% 		curvature of the cost function
% 	max_price : float
% 		price at which carbon dioxide can be removed from atmosphere in unlimited scale
% 	tech_const : float
% 		determines the degree of exogenous technological improvement over time. A number
% 			of 1.0 implies 1 percent per year lower cost
% 	tech_scale : float
% 		determines the sensitivity of technological change to previous mitigation
% 	cons_at_0 : float
% 		initial consumption. Default $30460bn based on US 2010 values.
% 	cbs_level : float
% 		constant
% 	cbs_deriv : float
% 		constant
% 	cbs_b : float
% 		constant
% 	cbs_k : float
% 		constant
% 	cons_per_ton : float
% 		constant
    
    properties
        tree
        g
        a
        max_price
        tech_const
        tech_scale
        cbs_level
        cbs_deriv
        cbs_b
        cbs_k
        cons_per_ton
    end
    
    methods
        %constructor
        function obj = DLWCost(tree, emit_at_0, g, a, join_price, max_price, tech_const, tech_scale, cons_at_0)
            obj@Cost();
            obj.tree = tree; % tree is a TreeModel object
            obj.g = g; %92.08
            obj.a = a; %3.413
            obj.max_price = max_price;
            obj.tech_const = tech_const;
            obj.tech_scale = tech_scale;
            obj.cbs_level = (join_price / (g*a))^(1/(a-1));
            obj.cbs_deriv = obj.cbs_level / (join_price * (a - 1.0));
            obj.cbs_b = obj.cbs_deriv * (max_price - join_price) / obj.cbs_level;
            obj.cbs_k = obj.cbs_level * (max_price - join_price)^obj.cbs_b;
            obj.cons_per_ton = cons_at_0 / emit_at_0;
        end
        
        
        function r = cost(obj, period, mitigation, ave_mitigation)
            
%       Calculates the mitigation cost for the period. For details about the cost function
% 		see DLW-paper.
% 
% 		Parameters
% 		----------
% 		period : int
% 			period in tree for which mitigation cost is calculated
% 		mitigation : ndarray
% 			current mitigation values for period
% 		ave_mitigation : ndarray
% 			average mitigation up to this period for all nodes in the period
% 
% 		Returns
% 		-------
% 		ndarray
% 			cost
            
            % Period as input in Python code.
            years = obj.tree.decision_times(period + 1);
            tech_term = (1.0 - ((obj.tech_const + obj.tech_scale * ave_mitigation)/100)).^years;
            % cbs is a power function of mitigation
            cbs = obj.g * (mitigation.^obj.a);
            % check if backstop technology is implemented
            bool_arr = int32(mitigation < obj.cbs_level);
            
            % cost of traditional mitigation
            if all(bool_arr)
                r = (cbs .* tech_term) / obj.cons_per_ton;
            % cost with backstop technology
            else
                %cost of normal mitigation
                base_cbs = obj.g * obj.cbs_level.^(obj.a);
                bool_arr2 = int32(mitigation > obj.cbs_level);
                extension = ((mitigation - obj.cbs_level) .* obj.max_price...
                             - obj.cbs_b .* mitigation .* (obj.cbs_k./mitigation).^(1/obj.cbs_b) ./ (obj.cbs_b-1.0)...
                             + obj.cbs_b*obj.cbs_level .* (obj.cbs_k./obj.cbs_level).^(1/obj.cbs_b) ./ (obj.cbs_b-1.0));
                r = (cbs .* double(bool_arr) + (base_cbs + extension) .* double(bool_arr2)) .* tech_term ./ obj.cons_per_ton;
            end
        end


        function r = price(obj, years, mitigation, ave_mitigation)
            
%       Inverse of the cost function. Gives emissions price for any given
% 		degree of mitigation, average_mitigation, and horizon.
% 
% 		Parameters
% 		----------
% 		years : int y
% 			years of technological change so far
% 		mitigation : float
% 			mitigation value in node
% 		ave_mitigation : float
% 			average mitigation up to this period
% 
% 		Returns
% 		-------
% 		float
% 			the price.
            
            tech_term = (1 - ((obj.tech_const + obj.tech_scale*ave_mitigation)/100))^years;
            if mitigation < obj.cbs_level
                r = obj.g * obj.a * (mitigation^(obj.a - 1)) * tech_term;
            else
                r = (obj.max_price - (obj.cbs_k/mitigation)^(1/obj.cbs_b)) * tech_term;
            end
        end
  
    end
end