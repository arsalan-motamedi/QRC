%Arbitrary Functions
clear 
clc
close all

extraInputs = {'interpreter','latex','fontsize',15};

d = 20;
a = diag(sqrt(1:d-1),1);
rho_0 = state_prep(d, 3, 2, 1);

T = 1000;  %total simulation time 
dt = 0.05;  %integration time step 
nt = round(T/dt); 
time = (1:nt)*dt;
noise_lvls = 0.3;
Er = zeros(1,length(noise_lvls));
%Oscillator parameters
K = 0.05; kappa = 0.05; alpha = 0.1;% .05, .1, .1
param = [K, kappa, alpha];

X_hat = a+a';
P_hat =1j*(a-a');
[XV,XD] = eig(X_hat);
XDD = diag(XD);
XD = diag(tanh(XDD));

[PV,PD] = eig(P_hat);
PDD = diag(PD);
PD = diag(tanh(PDD));
theta = pi/2; %rand * pi;
Observable = cos(theta) *XV * XD * XV' + sin(theta) * PV * PD * PV';

%% Simple supervisor, basic oscillator
freq = 0.4;
X = sin(2*pi*freq*time);
%load('Data.mat');
%X = Data(1:400000);
sup_mg = X(1:floor(.85*length(X)));
nt_mg = length(sup_mg);
%% Reservoir's dynamics
for count = 1:length(noise_lvls)

noise_lvl = noise_lvls(count);
chunk_length = 300; %input size at each step
pred_length = 150; %how many samples to predict
N_chunk = floor(nt_mg/chunk_length)-1;
output = zeros(chunk_length, N_chunk);

Y = zeros(N_chunk, pred_length);
for i = 1:N_chunk
    Y(i,:) = sup_mg(i*chunk_length + 1: i*chunk_length + pred_length);
end
Y = Y';
Tspan = (0:chunk_length-1)*dt;


for i = 1:N_chunk
    input = sup_mg(1+(i-1)*chunk_length: i*chunk_length) + noise_lvl * randn(1, chunk_length); %picking N_chunk number of samples from the signal
    [T,rho] = ode45(@(t,rho) Quantum_osc(t,rho,Tspan,input,a, param), Tspan, reshape(rho_0,[d^2,1]));
    for j = 1:chunk_length
        D = reshape(rho(j,:), [d,d]);
        output(j,i) = real(trace(D*(Observable)));
    end
end

%% Training

gamma = 0.01;

W = Y * output' /( output * output' + gamma * eye(chunk_length) );
%% Evaluating the accuracy

e = 0;
for i = 1:N_chunk
    e =  e + norm(Y(:,i) - W * output(:,i),2);
end

e = e/(N_chunk*sqrt(chunk_length));

fprintf('Average training error = %d\n', e)

% Visulalization

y = reshape(Y,[pred_length*N_chunk,1]);
y_pred = reshape(W*output, [pred_length*N_chunk,1]);

%figure()
%plot(y);
%hold on
%plot(y_pred,'--')
%h = xlabel('t');
%set(h, 'FontSize', 15);
%h = ylabel('x(t), x_e(t)');
%set(h,'FontSize',15);
%title('Training accuracy for MG data');
%% Data reproduction

Train = floor(0.8 * length(X));
test_MG = X(Train+1:length(X));
test = zeros(1,length(test_MG));
test(1:chunk_length) = test_MG(1:chunk_length);
recorded_input = zeros(1, length(test_MG));
n_chunk = floor((length(test) - chunk_length)/pred_length); % #chunks in test data

for i = 1:n_chunk
    output2 = zeros(chunk_length,1);
    input = test(1+(i-1)*pred_length: (i-1) * pred_length + chunk_length) + noise_lvl * randn(1,chunk_length);
    [T,rho] = ode45(@(t,rho) Quantum_osc(t,rho,Tspan,input,a,param), Tspan, reshape(rho_0,[d^2,1]));
        for j = 1:chunk_length
            D = reshape(rho(j,:),[d,d]);
            output2(j) = real(trace(D * (Observable)));
        end
    test((i-1)*pred_length+chunk_length+1: i*pred_length+chunk_length) = W * output2;
end

%%
t = (1:length(test_MG))*dt;
%
num_of_datapoints = n_chunk * pred_length + chunk_length;
g = figure;
subplot(2,1,1)
plot(t(1:num_of_datapoints), test_MG(1:num_of_datapoints)+noise_lvl*randn(1,num_of_datapoints),'LineWidth', 2.0);
xlim([60, 100])
title('Noisy input to the reservoir', extraInputs{:})
xlabel('$t(s)$', extraInputs{:});
ylim([-2,2]);

subplot(2,1,2)
plot(t(1:num_of_datapoints), test(1:num_of_datapoints), 'LineWidth', 2.0)
xlim([60, 100])
title('Reservoir''s output', extraInputs{:})
xlabel('$t(s)$', extraInputs{:});
ylim([-2,2]);

%%
%noisy_F = figure;
%plot(t(1:num_of_datapoints), test_MG(1:num_of_datapoints),'-red');
%hold on
%plot(t(1:num_of_datapoints), test(1:num_of_datapoints),'blue');
%hold on
%plot(t(1:num_of_datapoints), test_MG(1:num_of_datapoints) - test(1:num_of_datapoints),'black');
%hold off
%extraInputs = {'interpreter','latex','fontsize',14};
%xlabel('$t$', extraInputs{:});
%xlim([60, 100])
%ylabel('Function Values', extraInputs{:});
%title('Robustness with respect to error', extraInputs{:});

%% Evaluation of Test Error
test_MG = X(Train+1:length(X));
slide = 150;

windows = floor((length(test_MG) - pred_length - chunk_length)/slide);
for i = 1:windows
    output2 = zeros(chunk_length,1);
    input = test_MG(1+(i-1)*slide: (i-1) * slide + chunk_length) + noise_lvl * randn(1,chunk_length);
    [T,rho] =  ode45(@(t,rho) Quantum_osc(t,rho,Tspan,input,a,param), Tspan, reshape(rho_0,[d^2,1]));
        for j = 1:chunk_length
            D = reshape(rho(j,:),[d,d]);
            output2(j) = real(trace(D * (Observable)));
        end
    temp = W * output2;
    Er(count) = Er(count) + norm( (temp(1:150))' - test_MG((i-1) * slide + chunk_length+ 1 :(i-1) * slide + chunk_length + 150 ) ,2) ;
end

Er(count) = Er(count) / (windows * sqrt(150));

end

%%
extraInputs = {'interpreter','latex','fontsize',14};
NoisePlt = figure();
plot(noise_lvls, Er, '--o');
xlabel('Variance of input noise', extraInputs{:});
ylabel('Test Error', extraInputs{:});