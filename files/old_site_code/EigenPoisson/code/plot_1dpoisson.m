computed_sol = dlmread('out_1dpoisson.out', '', 2, 0);

N = numel(computed_sol);
xx = linspace(0,1,N)';
true_sol = sin(pi*xx)/pi^2;
printf('Error: %f\n', max(abs(true_sol - computed_sol)));


figure();
title('Comparison of solutions')
hold on;
plot(xx, true_sol, 'DisplayName', 'True Solution');
plot(xx, computed_sol, 'DisplayName', 'Computed Solution');
hold off;
legend();
