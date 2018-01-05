classdef Storage

    properties
      tree 
    end

    methods        
  
      function obj = Storage()
            obj.tree = containers.Map('KeyType','int32', 'ValueType','any');
      end
      
      function obj = set_tree(obj, period, value)
          obj.tree(period) = value;
      end
      
    end

end
