function u_value = utilitycalc( m, const )
     
    if nargin<2 || isempty(const)
        const = 1;
    end
    
    t = TreeModel([0, 15, 45, 85, 185, 285, 385]);
    bau_default_model = DLWBusinessAsUsual(t);
    c = DLWCost(t, bau_default_model.emit_level(1), 92.08, 3.413, 2000.0, 2500.0,...
					1.5, 0.0, 30460.0);
    df = DLWDamage(t, bau_default_model, 0.015, [450, 650, 1000], 5);
    [df, ~] = df.damage_function(m,4);
    u = EZUtility(t, df, c, 5.0);
    u_value = u.utility(m);

end

