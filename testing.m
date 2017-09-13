t = TreeModel([0, 15, 45, 85, 185, 285, 385]);

bau_default_model = DLWBusinessAsUsual(t);

c = DLWCost(t, bau_default_model.emit_level(1), 92.08, 3.413, 2000.0, 2500.0,...
					1.5, 0.0, 30460.0);

df = DLWDamage(t, bau_default_model, 0.015, [450, 650, 1000], 5);

ds =  DamageSimulation(t, [450, 650, 1000], 6.0, 18.0, true, 0, NaN, 100.0, 0.015);
%simtmp = ds.pindyck_simulation();

utilobj = EZUtility (t, df, c, 5);
mitistrat = ones(1, t.num_decision_nodes);

%df.average_mitigation(mitistrat,5)
df.damage_function(mitistrat,4)
                