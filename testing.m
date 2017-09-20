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
[df, damage] = df.damage_function(mitistrat,4);

u = EZUtility(t, df, c, 5.0);

utility_tree = BigStorageTree(u.period_len, u.decision_times);
cons_tree = BigStorageTree(u.period_len, u.decision_times);
ce_tree = BigStorageTree(u.period_len, u.decision_times);
cost_tree = SmallStorageTree(u.decision_times);

% Verify utility: Incorrect
% x = u.utility(mitistrat);

% Verify end_period_utility: Correct
% u.end_period_utility(mitistrat, utility_tree, cons_tree, cost_tree);

% Verify certain_equivalence: Correct
% period = 45;
% damage_period = utility_tree.between_decision_times(period);
% cetest = u.certain_equivalence(period, damage_period, utility_tree);
% cetest

% Verify adjusted_utility: Incorrect
u.adjusted_utility(mitistrat,NaN,NaN,0.1,0.0,false)


