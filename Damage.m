classdef (Abstract) Damage
    
%   Abstract damage class for the EZ-Climate model.
% 
% 	Parameters
% 	----------
% 	tree : `TreeModel` object
% 		provides the tree structure used
% 	bau : `BusinessAsUsual` object
% 		business-as-usual scenario of emissions
% 
% 	Attributes
% 	----------
% 	tree : `TreeModel` object
% 		provides the tree structure used
% 	bau : `BusinessAsUsual` object
% 		business-as-usual scenario of emissions
	
     properties(Access = private)
         
        % Python code for metaclass: __metaclass__ = ABCMeta
        % damageMeta is meta.class object in MATLAB
        damageMeta = ?Damage;        
     end
     
     properties
         tree
         bau
     end
     
     methods (Abstract = true)
%        The average_mitigation function should return a 1D array of the
%        average mitigation for every node in the period.
         r1 = average_mitigation(obj);
         
%       The damage_function should return a 1D array of the damages for
% 		every node in the period.         
         r2 = damage_function(obj);
     end
     
     methods
         % Constructor
         function obj = Damage(tree, bau)
             obj.tree = tree;
             obj.bau = bau;
         end
         
        % Set property of tree.final_states_prob
        function obj = set.tree(obj,value)
            obj.tree = value;
        end      
        
     end  
end