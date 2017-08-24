classdef (Abstract) BusinessAsUsual
    
%     Abstract BAU class for the EZ-Climate model.
% 
%     Parameters
%     ----------
%     ghg_start : float
%         today's GHG-level
%     ghg_end : float
%         GHG-level in the last period
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
       
%     properties(Access = private)
%         
%     end
 
    properties (Access = private)
        % Python code for metaclass: __metaclass__ = ABCMeta
        % bauMeta is meta.class object in MATLAB
        bauMeta = ?BusinessAsUsual;
        
        emission_by_decisions
        emission_per_period
        emission_to_ghg
        emission_to_bau
        bau_path
    end
    
    properties 
        ghg_start
        ghg_end
    end
    
    methods (Abstract = true)
        r = emission_by_time(obj);
    end
          
    methods
        % Constructor
        function obj = BusinessAsUsual(ghg_start, ghg_end)            
            obj.ghg_start = ghg_start;
            obj.ghg_end = ghg_end;                      
        end
       
    end
       
end