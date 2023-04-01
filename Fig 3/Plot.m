figure;

load 'Qkappa-6.mat'
plot(kappa, Quantumness, '--o', 'LineWidth', 2.0);
hold on
load 'Qkappa-mix.mat'
plot(kappa, Quantumness, '--o', 'LineWidth', 2.0);
load 'Qkappa-cat.mat'
plot(kappa, Quantumness, '--o', 'LineWidth', 2.0);
load 'Qkappa-coh.mat'
plot(kappa, Quantumness, '--o', 'LineWidth', 2.0);
legend('ket 6', 'mix', 'cat', 'coherent');
xlabel('$\kappa$', extraInputs{:});
ylabel('Quantumness', extraInputs{:});