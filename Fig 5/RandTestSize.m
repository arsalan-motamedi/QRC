% Trining Size

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
%Oscillator parameters
K = 0.05; kappa = 0.1; alpha = 1.2;% .05, .1, 1
param = [K,kappa,alpha];
fractions = [0.01, 0.02, 0.03, 0.04, 0.1, 0.2, 0.5];


load('LargeTdata.mat');
X = X - mean(X);
Training_size = floor(fractions * length(X));

chunk_length = 200; %input size at each step
pred_length = 100; %how many samples to predict
N_max = floor(floor(.85 * length(X))/chunk_length)-1;

% Osc. parameters
d = 20;
a = diag(sqrt(1:d-1),1);
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

%---- Coherent State
Er = zeros(5,length(Training_size));
rho_0 = coherent_state(1+1j, d);
for count = 1:length(Training_size)
    sup_mg = X(1:Training_size(count));
    nt_mg = Training_size(count);
    N_chunk = floor(nt_mg/chunk_length)-1;
    output = zeros(chunk_length, N_chunk);
    for ran=1:5
        display('Hurray');
        Ran = randperm(N_max);
        index = Ran(1:N_chunk);    
        Y = zeros(N_chunk, pred_length);
        for i = 1:N_chunk
            Y(i,:) = X(index(i)*chunk_length + 1: index(i)*chunk_length + pred_length);
        end
        Y = Y';
        Tspan = (0:chunk_length-1)*dt;


        for i = 1:N_chunk
            input = X(1+(index(i)-1)*chunk_length: index(i)*chunk_length); %picking N_chunk number of samples from the signal
            [T,rho] = ode45(@(t,rho) Quantum_osc(t,rho,Tspan,input,a, param), Tspan, reshape(rho_0,[d^2,1]));
            for j = 1:chunk_length
                D = reshape(rho(j,:), [d,d]);
                output(j,i) = real(trace(D*(Observable)));
            end
        end
        gamma = 0.01; W = Y * output' /( output * output' + gamma * eye(chunk_length));
        
        %Test Error
        Train = floor(0.85 * length(X) - 200);
        test_MG = X(Train+1:length(X));
        slide = 100;
        windows = floor((length(test_MG) - pred_length - chunk_length)/slide);
        for i = 1:windows
            output2 = zeros(chunk_length,1);
            input = test_MG(1+(i-1)*slide: (i-1) * slide + chunk_length);
            [T,rho] =  ode45(@(t,rho) Quantum_osc(t,rho,Tspan,input,a,param), Tspan, reshape(rho_0,[d^2,1]));
                for j = 1:chunk_length
                    D = reshape(rho(j,:),[d,d]);
                    output2(j) = real(trace(D * (Observable)));
                end
            temp = W * output2;
            Er(ran,count) = Er(ran,count) + norm( (temp) - test_MG((i-1) * slide + chunk_length+ 1 :(i-1) * slide + chunk_length + 100 ) ,2) ;
        end
    end
end
Er = Er / (windows * sqrt(100));
mEr_coh = sum(Er)/5;
%% ---- Cat State
Er = zeros(5,length(Training_size));
rho_0 = state_prep(d,1.2,2,1);
for count = 1:length(Training_size)
    sup_mg = X(1:Training_size(count));
    nt_mg = Training_size(count);
    N_chunk = floor(nt_mg/chunk_length)-1;
    output = zeros(chunk_length, N_chunk);
    for ran=1:5
        display('Hurray');
        Ran = randperm(N_max);
        index = Ran(1:N_chunk);    
        Y = zeros(N_chunk, pred_length);
        for i = 1:N_chunk
            Y(i,:) = X(index(i)*chunk_length + 1: index(i)*chunk_length + pred_length);
        end
        Y = Y';
        Tspan = (0:chunk_length-1)*dt;


        for i = 1:N_chunk
            input = X(1+(index(i)-1)*chunk_length: index(i)*chunk_length); %picking N_chunk number of samples from the signal
            [T,rho] = ode45(@(t,rho) Quantum_osc(t,rho,Tspan,input,a, param), Tspan, reshape(rho_0,[d^2,1]));
            for j = 1:chunk_length
                D = reshape(rho(j,:), [d,d]);
                output(j,i) = real(trace(D*(Observable)));
            end
        end
        gamma = 0.01; W = Y * output' /( output * output' + gamma * eye(chunk_length));
        
        %Test Error
        Train = floor(0.85 * length(X) - 200);
        test_MG = X(Train+1:length(X));
        slide = 100;
        windows = floor((length(test_MG) - pred_length - chunk_length)/slide);
        for i = 1:windows
            output2 = zeros(chunk_length,1);
            input = test_MG(1+(i-1)*slide: (i-1) * slide + chunk_length);
            [T,rho] =  ode45(@(t,rho) Quantum_osc(t,rho,Tspan,input,a,param), Tspan, reshape(rho_0,[d^2,1]));
                for j = 1:chunk_length
                    D = reshape(rho(j,:),[d,d]);
                    output2(j) = real(trace(D * (Observable)));
                end
            temp = W * output2;
            Er(ran,count) = Er(ran,count) + norm( (temp) - test_MG((i-1) * slide + chunk_length+ 1 :(i-1) * slide + chunk_length + 100 ) ,2) ;
        end
    end
end
Er = Er / (windows * sqrt(100));
mEr_cat = sum(Er)/5;
%% ---- Mixed state
Er = zeros(5,length(Training_size));
rho_0 = state_prep(d,1.2,2,0);
%rho_0(7,7) = 1;
for count = 1:length(Training_size)
    sup_mg = X(1:Training_size(count));
    nt_mg = Training_size(count);
    N_chunk = floor(nt_mg/chunk_length)-1;
    output = zeros(chunk_length, N_chunk);
    for ran=1:5
        display('Hurray');
        Ran = randperm(N_max);
        index = Ran(1:N_chunk);    
        Y = zeros(N_chunk, pred_length);
        for i = 1:N_chunk
            Y(i,:) = X(index(i)*chunk_length + 1: index(i)*chunk_length + pred_length);
        end
        Y = Y';
        Tspan = (0:chunk_length-1)*dt;


        for i = 1:N_chunk
            input = X(1+(index(i)-1)*chunk_length: index(i)*chunk_length); %picking N_chunk number of samples from the signal
            [T,rho] = ode45(@(t,rho) Quantum_osc(t,rho,Tspan,input,a, param), Tspan, reshape(rho_0,[d^2,1]));
            for j = 1:chunk_length
                D = reshape(rho(j,:), [d,d]);
                if(real(min(eig(D))) < -0.01)
                    display('Shit!');
                end
                output(j,i) = real(trace(D*(Observable)));
            end
        end
        gamma = 0.01; W = Y * output' /( output * output' + gamma * eye(chunk_length));
        
        %Test Error
        Train = floor(0.85 * length(X) - 200);
        test_MG = X(Train+1:length(X));
        slide = 100;
        windows = floor((length(test_MG) - pred_length - chunk_length)/slide);
        for i = 1:windows
            output2 = zeros(chunk_length,1);
            input = test_MG(1+(i-1)*slide: (i-1) * slide + chunk_length);
            [T,rho] =  ode45(@(t,rho) Quantum_osc(t,rho,Tspan,input,a,param), Tspan, reshape(rho_0,[d^2,1]));
                for j = 1:chunk_length
                    D = reshape(rho(j,:),[d,d]);
                    output2(j) = real(trace(D * (Observable)));
                end
            temp = W * output2;
            Er(ran,count) = Er(ran,count) + norm( (temp) - test_MG((i-1) * slide + chunk_length+ 1 :(i-1) * slide + chunk_length + 100 ) ,2) ;
        end
    end
end
Er = Er / (windows * sqrt(100));
mEr_mix = sum(Er)/5;
%% ---- Ket 6
Er = zeros(5,length(Training_size));
rho_0 = zeros(d,d);
rho_0(7,7) = 1;
for count = 1:length(Training_size)
    sup_mg = X(1:Training_size(count));
    nt_mg = Training_size(count);
    N_chunk = floor(nt_mg/chunk_length)-1;
    output = zeros(chunk_length, N_chunk);
    for ran=1:5
        display('Hurray');
        Ran = randperm(N_max);
        index = Ran(1:N_chunk);    
        Y = zeros(N_chunk, pred_length);
        for i = 1:N_chunk
            Y(i,:) = X(index(i)*chunk_length + 1: index(i)*chunk_length + pred_length);
        end
        Y = Y';
        Tspan = (0:chunk_length-1)*dt;


        for i = 1:N_chunk
            input = X(1+(index(i)-1)*chunk_length: index(i)*chunk_length); %picking N_chunk number of samples from the signal
            [T,rho] = ode45(@(t,rho) Quantum_osc(t,rho,Tspan,input,a, param), Tspan, reshape(rho_0,[d^2,1]));
            for j = 1:chunk_length
                D = reshape(rho(j,:), [d,d]);
                if(real(min(eig(D))) < -0.01)
                    display('Shit!');
                end
                output(j,i) = real(trace(D*(Observable)));
            end
        end
        gamma = 0.01; W = Y * output' /( output * output' + gamma * eye(chunk_length));
        
        %Test Error
        Train = floor(0.85 * length(X) - 200);
        test_MG = X(Train+1:length(X));
        slide = 100;
        windows = floor((length(test_MG) - pred_length - chunk_length)/slide);
        for i = 1:windows
            output2 = zeros(chunk_length,1);
            input = test_MG(1+(i-1)*slide: (i-1) * slide + chunk_length);
            [T,rho] =  ode45(@(t,rho) Quantum_osc(t,rho,Tspan,input,a,param), Tspan, reshape(rho_0,[d^2,1]));
                for j = 1:chunk_length
                    D = reshape(rho(j,:),[d,d]);
                    output2(j) = real(trace(D * (Observable)));
                end
            temp = W * output2;
            Er(ran,count) = Er(ran,count) + norm( (temp) - test_MG((i-1) * slide + chunk_length+ 1 :(i-1) * slide + chunk_length + 100 ) ,2) ;
        end
    end
end
Er = Er / (windows * sqrt(100));
mEr_6 = sum(Er)/5;
%% Classical Model
a_0 = 1+1j;
Er = zeros(5,length(Training_size));
for count = 1:length(Training_size)
    sup_mg = X(1:Training_size(count));
    nt_mg = Training_size(count);
    N_chunk = floor(nt_mg/chunk_length)-1;
    output = zeros(chunk_length, N_chunk);
    for ran=1:5
        display('Hurray');
        Ran = randperm(N_max);
        index = Ran(1:N_chunk);    
        Y = zeros(N_chunk, pred_length);
        for i = 1:N_chunk
            Y(i,:) = X(index(i)*chunk_length + 1: index(i)*chunk_length + pred_length);
        end
        Y = Y';
        Tspan = (0:chunk_length-1)*dt;


        for i = 1:N_chunk
            input = X(1+(index(i)-1)*chunk_length: index(i)*chunk_length); %picking N_chunk number of samples from the signal
            [T, a] = ode45(@(t,a) Classical_osc(t,a,Tspan,input,param), Tspan, a_0);
            output(:,i) = real(a);
        end
        gamma = 0.01; W = Y * output' /( output * output' + gamma * eye(chunk_length));
        
        %Test Error
        Train = floor(0.85 * length(X) - 200);
        test_MG = X(Train+1:length(X));
        slide = 100;
        windows = floor((length(test_MG) - pred_length - chunk_length)/slide);
        for i = 1:windows
            output2 = zeros(chunk_length,1);
            input = test_MG(1+(i-1)*slide: (i-1) * slide + chunk_length);
            [T, a] = ode45(@(t,a) Classical_osc(t,a,Tspan,input,param), Tspan, a_0);
            temp = W * real(a);
            Er(ran,count) = Er(ran,count) + norm( (temp) - test_MG((i-1) * ...
                slide + chunk_length+ 1 :(i-1) * slide + chunk_length + 100 ) ,2) ;
        end
    end
end
Er = Er / (windows * sqrt(100));
mEr_classic = sum(Er)/5;
%%
figure();
plot(log(Training_size), mEr_coh, '--o', 'LineWidth',2.0);
hold on
plot(log(Training_size), mEr_cat, '--o', 'LineWidth',2.0);
plot(log(Training_size), mEr_mix, '--o', 'LineWidth',2.0);
plot(log(Training_size), mEr_6, '--o', 'LineWidth',2.0);
plot(log(Training_size), mEr_classic, '--o', 'LineWidth',2.0);
xlabel('$\log($Training Size$)$', extraInputs{:});
legend('Coherent', 'Cat', 'Mixed', 'ket 6', 'Classical')
ylabel('Test Error', extraInputs{:});