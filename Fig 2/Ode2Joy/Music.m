% Music -- Quantumly!


clear
close all
clc 

T = 5000;  %Network parameters and total simulation time 
dt = 0.04; %Time step 
nt = round(T/dt);  %number of time steps

K = 0.01; kappa = 0.01; alpha = 0.03; %Oscillator's Parameters
param = [K, kappa, alpha];
%% Convert the sequence of notes and half notes in the ode2joyshort.mat into a teaching signal. 
% the file ode2joy short.mat contains 2 matrices, J and HN.  The matrix
% length corresponds to the number of notes while the matrix width
% corresponds to the note type and is the dimension of the teaching signal.
% J indicates the presence of a note while HN indicates the presence of a
% half note.  
%(This section is just copy-pasted from IZFORCESONG.m)
freq = 40; 
load ode2joyshort.mat;
nnotes = length(J);
nchord = min(size(J));
ds = (1000/freq)*nnotes; n1 = round(ds/dt);
ZS = abs(sin(pi*(1:1:n1)*dt*nnotes/(ds)));
ZS = repmat(ZS,nchord,1);
song = J'; 
nn = size(song);  

j = 1 ;
for i = 1:1:n1 
    if mod(i,round(1000/(freq*dt)))==0;
    j = j + 1;
    if j > nn(2); break; end 
    end
    ZS(:,i) = ZS(:,i).*song(:,j);
end 

for i = 1:1:15;
    if HN(i,1) > 0 
        q = length((i-1)*(1000/(freq*dt)):(i+1)*(1000/(freq*dt)));
        w = find(J(i,:)>0);
       ZS(w,(i-1)*(1000/(freq*dt)):(i+1)*(1000/(freq*dt)))= sin(pi*(1:1:q)/q);
    end 
end
zx = repmat(ZS,1,ceil(T/ds));

%% Oscillator's ado

K1 = 3; %Number of notes the oscillator sees
K2 = 3; %Number of notes it should predict at each step

chunk_length = floor(K1*n1/16); %input size at each step
pred_length = floor(K2*n1/16); %how many samples to predict
N_chunk = floor(length(zx)/chunk_length)-1;
Tspan = (0:chunk_length-1)*dt;
output = zeros(chunk_length, N_chunk, nn(1)); % 5 =nnotes

d = 8; %We chose 8 over 10 to reduce the run-time
a = diag(sqrt(1:d-1),1);
rho_0 = zeros(d,d);
rho_0(4,4) = 1;

%Osc dynamics
for i = 1:N_chunk
   for note = 1:nn(1)
       input = zx(note, 1+(i-1)*chunk_length: i*chunk_length);
       [T,rho] = ode45(@(t,rho) Quantum_osc(t,rho,Tspan,input,a, param), Tspan, reshape(rho_0,[d^2,1]));
       for j = 1:chunk_length
            D = reshape(rho(j,:), [d,d]);
            output(j,i,note) = real(trace(D * (a + a')));
        end
   end
end


%% Training

%Training Y's
Y = zeros(pred_length, N_chunk, nn(1));
for i = 1:N_chunk
    for note = 1:nn(1)
        Y(:,i,note) = zx(note, i*chunk_length + 1: i*chunk_length + pred_length);
    end
end


gamma = 0.01;
output_combined = zeros(nn(1)*chunk_length, N_chunk);
for i=1:N_chunk
    output_combined(:,i) = [output(:,i,1);output(:,i,2);output(:,i,3);output(:,i,4);output(:,i,5)];  
end
%%
A = output_combined' /( output_combined * output_combined' + gamma * eye(nn(1)*chunk_length) );

W = zeros(pred_length, nn(1) * chunk_length, nn(1));
%%
for i = 1:nn(1)
   W(:,:,i) = Y(:,:,i) * A; 
end
%%
e1 = 0;
for i = 1:N_chunk
e1 =  e1 + norm( Y(:,i,1) - W(:,:,1) * output_combined(:,i),1);
end
%%
y2_pred = reshape(W(:,:,2)*output_combined, [pred_length*N_chunk,1]);
y2 = reshape(Y(:,:,2),[pred_length*N_chunk,1]);
%%
figure()
plot(y2);
hold on
plot(y2_pred)
%%
prediction_song = zeros(nn(1),n1);
prediction_song(:,1:chunk_length) = zx(:,1:chunk_length);


n_chunk = floor((n1 - chunk_length)/pred_length)+1;
for i = 1:n_chunk
    output2 = zeros(chunk_length,nn(1));
    for note = 1:nn(1)
        input = prediction_song(note,1+(i-1)*pred_length: (i-1) * pred_length + chunk_length);
        [T,rho] = ode45(@(t,rho) Quantum_osc(t,rho,Tspan,input,a, param), Tspan, reshape(rho_0,[d^2,1]));
        for j = 1:chunk_length
            D = reshape(rho(j,:), [d,d]);
            output2(j,note) = real(trace(D * (a + a')));
        end
    end
    output_combined2 = reshape(output2, [nn(1)*chunk_length,1]);
    for Note = 1:nn(1)
        prediction_song(Note, (i-1)*pred_length+chunk_length+1: i*pred_length+chunk_length) = W(:,:,Note) * output_combined2;
    end
end
%% Plotting the output note

figure()
subplot(5,1,1)
plot(dt*(1:length(prediction_song))/100,prediction_song(1,:),'black')
xlabel('t');
ylabel('Amplitude')
title('C-note')
ylim([-1.1, 1.1])
xlim([0, ds/100])

subplot(5,1,2)
plot(dt*(1:length(prediction_song))/100,prediction_song(2,:),'green')
xlabel('t');
ylabel('Amplitude')
title('D-note')
ylim([-1.1, 1.1])
xlim([0, ds/100])


subplot(5,1,3)
plot(dt*(1:length(prediction_song))/100,prediction_song(3,:),'m')
xlabel('t');
ylabel('Amplitude')
title('E-note')
ylim([-1.1, 1.1])
xlim([0, ds/100])


subplot(5,1,4)
plot(dt*(1:length(prediction_song))/100,prediction_song(4,:),'red')
xlabel('t');
ylabel('Amplitude')
title('F-note')
ylim([-1.1, 1.1])
xlim([0, ds/100])

subplot(5,1,5)
plot(dt*(1:length(prediction_song))/100, prediction_song(5,:),'blue')
xlabel('t');
ylabel('Amplitude')
title('G-note')
ylim([-1.1, 1.1])
xlim([0, ds/100])
