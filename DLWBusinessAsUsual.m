classdef DLWBusinessAsUsual < BusinessAsUsual
    
%     Business-as-usual scenario of emissions. Emissions growth is assumed to slow down 
%     exogenously - these assumptions represent an attempt to model emissions growth in a 
%     business-as-usual scenario that is in the absence of incentives.
% 
%     Parameters
%     ----------
%     ghg_start : float
%         today's GHG-level
%     ghg_end : float
%         GHG-level in the last period
%     emit_time : ndarray or list
%         time, in years, from now when emissions occurs
%     emit_level : ndarray or list
%         emission levels in future times `emit_time`
% 
%     Attributes
%     ----------
%     ghg_start : float
%         today's GHG-level
%     ghg_end : float
%         GHG-level in the last period
%     emission_by_decisions : ndarray
%         emissions at decision time periods??
%     emission_per_period : ndarray
%         total emission at decision time period??
%     emission_to_ghg : ndarray
%         GHG levels in decision time period??
%     emission_to_bau : float
%         constant for converting GHG to emission??
%     emit_time : ndarray or list
%         time, in years, from now when emissions occurs
%     emit_level : ndarray or list
%         emission levels in future times `emit_time`
   
    properties
        emit_time
        emit_level
        tree
    end
    
    properties(Dependent)
        emission_by_decisions
        emission_per_period
        emission_to_ghg
        emission_to_bau
        bau_path
    end
    
    methods        
        % Constructor
        function obj = DLWBusinessAsUsual(tree, ghg_start, ghg_end, emit_time, emit_level)
            if nargin == 1
                ghg_start = 400.0;
                ghg_end = 1000.0;
                emit_time = [0,30,60];
                emit_level = [52.0, 70.0, 81.4];
            end
            
            obj@BusinessAsUsual(ghg_start, ghg_end); 
            obj.tree = tree;
            obj.emit_time = emit_time;
            obj.emit_level = emit_level;            
        end
    
        
        function r = emission_by_time(obj, time)
            
%         Returns the BAU emissions at any time
% 
%         Parameters
%         ----------
%         time : int 
%             future time period in years
% 
%         Returns
%         -------
%         float
%             emission
            
            if time < obj.emit_time(2)
                emissions = obj.emit_level(1) + double(time) / (obj.emit_time(2) - obj.emit_time(1)) * (obj.emit_level(2) - ...
                    obj.emit_level(1));
            elseif time < obj.emit_time(3)
                emissions = obj.emit_level(2) + double(time - obj.emit_time(2)) / (obj.emit_time(3)... 
                        - obj.emit_time(2)) * (obj.emit_level(3) - obj.emit_level(2));
            else
                emissions = obj.emit_level(3);
            end
            r = emissions;
        end
        
        
        % Set the emission_by_decisions property
        function r = get.emission_by_decisions(obj)
            num_periods = obj.tree.num_periods;
            r = zeros(1, num_periods);
            r(1) = obj.emission_by_time(obj.tree.decision_times(1));
           
            for n = 2:num_periods
                r(n) = obj.emission_by_time(obj.tree.decision_times(n));                
            end  
        end
        
        % Set the emission_per_period property
        function r = get.emission_per_period(obj)
            num_periods = obj.tree.num_periods;
            r = zeros(1, num_periods);            
            period_len = obj.tree.decision_times(2:end) - obj.tree.decision_times(1:end-1);
            
            for n = 2:num_periods
                r(n) = period_len(n) * (mean(obj.emission_by_decisions(n-1)));                
            end  
        end
        
        % Set the emission_to_ghg property
        function r = get.emission_to_ghg(obj)
            r = (obj.ghg_end - obj.ghg_start) * obj.emission_per_period / sum(obj.emission_per_period);
        end
        
        % Set the emission_to_bau property
        function r = get.emission_to_bau(obj)
            r = obj.emission_to_ghg(end) / obj.emission_per_period(end);
        end
        
        % Set the bau_path property
        function r = get.bau_path(obj)
            num_periods = obj.tree.num_periods;
            r = zeros(1, num_periods);
            r(1) = obj.ghg_start;
            
            for n = 2:num_periods
                r(n) = r(n-1) + obj.emission_per_period(n) * obj.emission_to_bau;
            end
        end
        
        
%         function bau_emissions_setup(obj, tree)
%             
% %         Create default business as usual emissions path. The emission rate in each period is 
% %         assumed to be the average of the emissions at the beginning and at the end of the period.
% % 
% %         Parameters
% %         ----------
% %         tree : `TreeModel` object
% %             provides the tree structure used
%                        
%             num_periods = tree.num_periods; 
%             obj.emission_by_decisions = zeros(1, num_periods);            
%             obj.emission_per_period = zeros(1, num_periods);
%             obj.bau_path = zeros(1, num_periods);
%             obj.bau_path(1) = obj.ghg_start;
%             obj.emission_by_decisions(1) = obj.emission_by_time(tree.decision_times(1));
%             period_len = tree.decision_times(2:end) - tree.decision_times(1:end-1);
%             
%             for n = 2:num_periods
%                 obj.emission_by_decisions(n) = obj.emission_by_time(tree.decision_times(n));
%                 obj.emission_per_period(n) = period_len(n) * (mean(obj.emission_by_decisions(n-1:n)));
%             end
%             
%             % The total increase in ghg level of 600 (from 400 to 1000) in the bau path is allocated over time
%             obj.emission_to_ghg = (obj.ghg_end - obj.ghg_start) * obj.emission_per_period / sum(obj.emission_per_period);
%             obj.emission_to_bau = obj.emission_to_ghg(end) / obj.emission_per_period(end);
%             for n = 2:num_periods
%                 obj.bau_path(n) = obj.bau_path(n-1) + obj.emission_per_period(n) * obj.emission_to_bau;
%             end
%         end
%         
        
    end
end