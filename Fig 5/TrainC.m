%Grand Plot

figure();

extraInputs = {'interpreter','latex','fontsize',15};

load('TE1.mat');
plot(log(Training_size), mean(Er), '--o', 'LineWidth',2.0);
hold on

load('TE-Cat.mat');
plot(log(Training_size), mean(Er), '--*', 'LineWidth',1.0);
hold on

load('TE-mix.mat');
plot(log(Training_size), mean(Er), '--s', 'LineWidth',1.0);
hold on

load('TE-coherent.mat');
plot(log(Training_size), mean(Er), '--d', 'LineWidth',1.0);
hold on

load('TE-classic.mat');
plot(log(Training_size), mean(Er), '-.o', 'LineWidth',1.0);

legend('ket 6', 'cat', 'mix', 'coherent state', 'classical model');
xlabel('$\log($Training Size$)$', extraInputs{:});
ylabel('Test Error', extraInputs{:});