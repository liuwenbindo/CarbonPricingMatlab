classdef DLWCost < Cost
    
    %addded 07/28/17
    properties
        tree
        g
        a
        max_price
        tech_const
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
            obj.tree = tree; %tree is TreeModel object
            obj.g = g; %92.08
            obj.a = a; %3.413
            obj.max_price = max_price;
            obj.tech_const = tech_const;
            obj.tech_scale = tech_scale;
            obj.cbs_level = join_price / ((g*a)^(1/(a-1)));
            obj.cbs_deriv = obj.cbs_level / (join_price * (a - 1));
            obj.cbs_b = obj.cbs_deriv * (max_price - join_price) / obj.cbs_level;
            obj.cbs_k = obj.cbs_level * (max_price - join_price)^obj.cbs_b;
            obj.cons_per_ton = cons_at_0 / emit_at_0;
        end
        %method 1
        function r = cost(obj, period, mitigation, ave_mitigation)
            years = obj.tree.decision_times(period);
            tech_term = (1 - (obj.tech_const + obj.tech_scale*ave_mitigation)/100)^years;
            cbs = obj.g * (mitigation^obj.a);
            bool_arr = int32(mitigation > obj.cbs_level);
            if all(bool_arr)
                r = (cbs * tech_term) / obj.cons_per_ton;
            else
                base_cbs = obj.g * obj.cbs_level^(obj.a);
                bool_arr2 = int32(mitigation > obj.cbs_level);
                extension = ((mitigation - obj.cbs_level) * obj.max_price...
                             - obj.cbs_b * mitigation * (self.cbs_k/mitigation)^(1/self.cbs_b) / (self.cbs_b-1.0)...
                             + obj.cbs_b*obj.cbs_level * (obj.cbs_k/obj.cbs_level)^(1/obj.cbs_b) / (obj.cbs_b-1.0));
                r = (cbs * bool_arr + (base_cbs + extension)*bool_arr2) * tech_term / obj.cons_per_ton;
            end
        end
        %method 2
        function r = price(obj, years, mitigation, ave_mitigation)
            tech_term = (1 - ((obj.tech_const + obj.tech_scale*ave_mitigation)/100))^years;
            if mitigation < obj.cbs_level
                r = obj.g * obj.a * (mitigation^(obj.a - 1)) * tech_term;
            else
                r = (obj.max_price - (obj.cbs_k/mitigation)^(1/obj.cbs_b)) * tech_term;
            end
        end
    end
end