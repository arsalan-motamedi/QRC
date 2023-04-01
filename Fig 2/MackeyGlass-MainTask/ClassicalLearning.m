% Learning With the Classical Model

clear 
clc
close all

extraInputs = {'interpreter','latex','fontsize',14};
dt = 0.1;  %integration time step 
noise_lvl = 0; %set non-zero (0.1 say) for introducing error to input

%% Oscillator parameters
K = 0.05; kappa = 0.1; alpha = 1.2;
a_0 = 1+1j; %initial state

%% Simple supervisor, basic oscillator
load('LargeTdata.mat');
X = X - mean(X);
sup_mg = X(1:floor(.5*length(X))); %Using half of data for training (you could use as much as you wish!)
nt_mg = length(sup_mg);
%% Reservoir's dynamics

chunk_length = 200; %input size at each step
pred_length = 100; %how many samples to predict
N_chunk = floor(nt_mg/chunk_length)-1;
output = zeros(chunk_length, N_chunk);

Y = zeros(N_chunk, pred_length);
for i = 1:N_chunk
    Y(i,:) = sup_mg(i*chunk_length + 1: i*chunk_length + pred_length);
end
Y = Y';
Tspan = (0:chunk_length-1)*dt;


param = [K, kappa, alpha];
for i = 1:N_chunk
    input = sup_mg(1+(i-1)*chunk_length: i*chunk_length); %picking N_chunk number of samples from the signal
    [T, a] = ode45(@(t,a) Classical_osc(t,a,Tspan,input,param), Tspan, a_0);
    output(:,i) = real(a);
end
%% Training

gamma = 0.01;

W = Y * output' /( output * output' + gamma * eye(chunk_length) );
%% Data Reproduction
Train = floor(0.85 * length(X));
test_MG = X(Train+1:length(X));
test = zeros(1,length(test_MG));
test(1:chunk_length) = test_MG(1:chunk_length);

n_chunk = floor((length(test) - chunk_length)/pred_length); % #chunks in test data

for i = 1:n_chunk
   input = test(1+(i-1)*pred_length: (i-1) * pred_length + chunk_length);
   [T, a] = ode45(@(t,a) Classical_osc(t,a,Tspan,input,param), Tspan, a_0);
   test((i-1)*pred_length+chunk_length+1: i*pred_length+chunk_length) = W * real(a);
end
num_of_datapoints = n_chunk * pred_length + chunk_length;
%%
figure();
time = dt*(1:num_of_datapoints);
plot(time(1:1200),test(1:1200));
hold on
plot(time(1:1200),test_MG(1:1200));
hold on
legend('Produced by reservoir','Actual test data');
line_x = linspace(200*dt,200*dt,10);
line_y = linspace(-.5,.5,10);
plot(line_x,line_y,'--');
xlabel('time (s)', extraInputs{:});
ylabel('prediction', extraInputs{:});
hold off
%% Phase Diagrams
delay = 17; %Choose as you wish
[x,y] = phase_space(test_MG(1:1000),17);
[z,w] = phase_space(test(1:1000),17);
figure();
subplot(1,2,1);
plot(x,y);
xlabel('x(t-1.7)', extraInputs{:}); ylabel('x(t)', extraInputs{:});
subplot(1,2,2);
plot(z,w);
xlabel('x(t-1.7)', extraInputs{:}); ylabel('x(t)', extraInputs{:});
