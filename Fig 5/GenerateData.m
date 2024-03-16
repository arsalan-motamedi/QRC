
%Codes for Fig. 6 of QRC paper.
%Autor: Arsalan Motamedi
tic
clear 
clc
close all
extraInputs = {'interpreter','latex','fontsize',14};
rng(1);
Q = 1; %Set it equal to 0 to go classic!
T = 3800;  %total simulation time 
dt = 0.1;  %integration time step 
nt = round(T/dt); 
time = (1:nt)*dt;
noise_lvl = 0;

%Oscillator parameters
K = 0.12; kappa = 0; alpha = 1.2;% .05, .1, 1

%ORIGINAL VLAS: K = 0.05; kappa = .1; alpha = 1.2;
param = [K, kappa, alpha];

chunk_length = 200; %input size at each step
pred_length = 100; %how many samples to predict

load('LargeTdata.mat');
X = X - mean(X);
test_MG = X(floor(0.85*length(X)):length(X)); %Test data

sup = X(1:(ceil(0.5*length(X))));
nt = length(sup);
N_chunk = floor(nt/chunk_length)-1;
output = zeros(chunk_length, N_chunk);
%%
Y = zeros(N_chunk, pred_length);
for i = 1:N_chunk
Y(i,:) = sup(i*chunk_length + 1: i*chunk_length + pred_length);
end

Y = Y';

d = 25;
a = diag(sqrt(1:d-1),1);
rho_0 = zeros(d,d);
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
Cmax = 5; % number of points per dimension
dims = 4:10;
Quant = zeros(length(dims),Cmax); %To store quantumness values
mrho = zeros(length(dims),Cmax, d^2);

Er = zeros(length(dims),Cmax); %TODO: add this line to main GitHub repo

Tspan = (0:chunk_length-1)*dt;
for dim = 1:length(dims)

for count = 1:Cmax
A = rand(d,d);
psi = randn(dims(dim),1) + 1j* randn(dims(dim),1);
psi = psi/norm(psi);
rho_0(1:dims(dim),1:dims(dim)) = psi*psi';
mrho(dim,count,:) = reshape(rho_0, [1,d^2]);
clear I;
I = zeros(1,N_chunk);
for i = 1:N_chunk
%picking N_chunk number of samples from the signal
input = sup(1+(i-1)*chunk_length: i*chunk_length) + noise_lvl * ...
    randn(chunk_length,1);
%integrating the PDE describing the dynamics
[T,rho] = ode45(@(t,rho) Quantum_osc(t,rho,Tspan,input,a, param), Tspan,...
    reshape(rho_0,[d^2,1]));
for j = 1:chunk_length
D = reshape(rho(j,:), [d,d]);
output(j,i) = real(trace(D*(Observable)));
end
LJ = LeeJeong(rho);
I(i) = sum(sum(LJ(LJ > 0))); %Computing the area above I=0 in the evolution
clear LJ;
end
Quant(dim,count) = mean(I);
%Linear regression to find W
gamma = 0.01; %TODO: change for K=.12, k=0.02
W = Y * output' /( output * output' + gamma * ...
    eye(chunk_length) );

%TODO: training error

%Test Error evaluation
slide = 100; windows = floor((length(test_MG) - pred_length - ...
    chunk_length)/slide);

for i = 1:windows
output2 = zeros(chunk_length,1);
input = test_MG(1+(i-1)*slide: (i-1) * slide + chunk_length);
[T,rho] =  ode45(@(t,rho) Quantum_osc(t,rho,Tspan,input,a,param),...
    Tspan, reshape(rho_0,[d^2,1]));
for j = 1:chunk_length
D = reshape(rho(j,:),[d,d]);
output2(j) = real(trace(D * (Observable)));
end
temp = W * output2;
Er(dim, count) = Er(dim, count) +...
norm( (temp(1:100)) - test_MG((i-1) * slide + chunk_length+ 1 :...
(i-1) * slide + chunk_length + 100 ) ,2);
end

end
end
save('CorrectedI.mat', 'Quant');
%%
Er = Er / (windows * sqrt(100));
figure();
plot(Quant,Er, 'ob');
xlabel('Quantumness', extraInputs{:});
ylabel('Test Error', extraInputs{:});
title('Random State Taining', extraInputs{:});
save('QuantumnessData.mat', 'Quant','Er', 'mrho');
R = corrcoef(Er,Quant);
display(R(1,2))% Shows the r-value for qurntumness-error correlation.
toc