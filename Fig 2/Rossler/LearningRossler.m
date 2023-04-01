%(QRC - time series prediction)
clear 
clc
close all

extraInputs = {'interpreter','latex','fontsize',20};

Q = 1; %Set it equal to 0 to go classic!
dt = 0.1;  %integration time step 
noise_lvl = 0;
%Oscillator parameters
K = 0.01; kappa = 0.2; alpha = 1; %lambda = 0;% .05, .1, 1
param = [K, kappa, alpha];
%% Simple supervisor, basic oscillator
load('RosslerData.mat');
%load('Input.mat');
%X = Y;
%clear Y;
%X = Data(1:400000);
X = X1 - mean(X1);
X = X/10;
sup_mg = X(1:floor(.85*length(X)));
nt_mg = length(sup_mg);
%% Reservoir's dynamics

chunk_length = 220; %input size at each step
pred_length = 110; %how many samples to predict
N_chunk = floor(nt_mg/chunk_length)-1;
output = zeros(chunk_length, N_chunk);

Y = zeros(N_chunk, pred_length);
for i = 1:N_chunk
    Y(i,:) = sup_mg(i*chunk_length + 1: i*chunk_length + pred_length);
end
Y = Y';
Tspan = (0:chunk_length-1)*dt;

if(Q == 1)
% Osc. parameters
d = 20;
a = diag(sqrt(1:d-1),1);
rho_0 = zeros(d,d);%squeezed_state(1+1j, 0, pi, d);
rho_0(7,7) = 1;
%rho_0 =  expm((1+1j) * a' - (1-1j) * a) * rho_0 * (expm((1+1j) * a' - (1-1j) * a))';
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

for i = 1:N_chunk
    input = sup_mg(1+(i-1)*chunk_length: i*chunk_length); %picking N_chunk number of samples from the signal
    [T,rho] = ode45(@(t,rho) Quantum_osc(t,rho,Tspan,input,a, param), Tspan, reshape(rho_0,[d^2,1]));
    for j = 1:chunk_length
        D = reshape(rho(j,:), [d,d]);
        output(j,i) = real(trace(D*(Observable)));
    end
end

else
a_0 = 1+1j;    
for i = 1:N_chunk
    input = sup_mg(1+(i-1)*chunk_length: i*chunk_length); %picking N_chunk number of samples from the signal
    [T, a] = ode45(@(t,a) Classical_osc(t,a,Tspan,input,param), Tspan, a_0);
    output(:,i) = tanh(imag(a));
end
end
%% Training

gamma = 0.01;

W = Y * output' /( output * output' + gamma * eye(chunk_length) );
%% Evaluating the accuracy

e = 0;
for i = 1:N_chunk
    e =  e + norm(Y(:,i) - W * output(:,i),1);
end

e = e/(N_chunk*chunk_length);

fprintf('Average error = %d', e)

% Visulalization

y = reshape(Y,[pred_length*N_chunk,1]);
y_pred = reshape(W*output, [pred_length*N_chunk,1]);

figure()
plot(y);
hold on
plot(y_pred,'--')
xlabel('$t$', extraInputs{:});
ylabel('x(t), x_e(t)', extraInputs{:});
title('Training accuracy for MG data', extraInputs{:});
%% Data reproduction

Train = floor(0.85 * length(X) - 200);
test_MG = X(Train+1:Train + 1000);
test = zeros(1,length(test_MG));
test(1:chunk_length) = test_MG(1:chunk_length);

n_chunk = floor((length(test) - chunk_length)/pred_length); % #chunks in test data

if (Q == 1)

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

else
    
   for i = 1:n_chunk
       input = test(1+(i-1)*pred_length: (i-1) * pred_length + chunk_length);
       [T, a] = ode45(@(t,a) Classical_osc(t,a,Tspan,input,param), Tspan, a_0);
       test((i-1)*pred_length+chunk_length+1: i*pred_length+chunk_length) = W * tanh(imag(a));
   end
end
%%
t = (1:length(test_MG))*dt;
num_of_datapoints = n_chunk * pred_length + chunk_length;

Line_y = linspace(min(X), max(X), 10);
Line_x = linspace(chunk_length * dt, chunk_length * dt, 10);

g = figure;
subplot(2,1,1)
plot(t(1:num_of_datapoints), test_MG(1:num_of_datapoints))
xlim([1, t(num_of_datapoints)])
title('Actual series; test data', extraInputs{:})
xlabel('$t$', extraInputs{:});
ylabel('$x(t)$', extraInputs{:});
hold on
plot(Line_x, Line_y, '--red');
hold off

subplot(2,1,2)
plot(t(1:num_of_datapoints), test(1:num_of_datapoints))
xlim([1, t(num_of_datapoints)])
title('Predicted series', extraInputs{:})
xlabel('$t$', extraInputs{:});
ylabel('$\hat{x}(t)$', extraInputs{:});
hold on
plot(Line_x, Line_y, '--red');
hold off
%%
figure();
plot(t(1:num_of_datapoints), test_MG(1:num_of_datapoints) - (test(1:num_of_datapoints)));
%%
g2 = figure;
plot(t(1:num_of_datapoints), test_MG(1:num_of_datapoints), 'blue', 'LineWidth', 2.0)
xlim([1, t(num_of_datapoints)])
%title('Comparison on Test Data', extraInputs{:})
xlabel('$t$', extraInputs{:});
ylabel('$x(t)$', extraInputs{:});
hold on
plot(t(1:num_of_datapoints), test(1:num_of_datapoints),'red', 'LineWidth', 2.0)
xlim([1, t(num_of_datapoints)])
xlabel('$t$', extraInputs{:});
ylabel('$\hat{x}(t)$', extraInputs{:});
ylim([min(X), max(X)+.3]);
legend('Produced by reservoir','Actual test data');
plot(Line_x, Line_y, '--k', 'LineWidth', 1.5);
hold off
%%
Z1 = test(1:num_of_datapoints); Xp1 = test_MG(1:num_of_datapoints);
save('Fir.mat','Z1','Xp1');