% Full Plot
clc
clear all
close all

load QuantumnessData.mat
load CorrectedI.mat
I = reshape(Quant, [140,1]); Error = reshape(Er, [140,1]);
extraInputs = {'interpreter','latex','fontsize',14};

i4 = []; i10=[];
for i = 1:140
if( mod(i,7) == 1)
i4 = [i4, i];
end
if( mod(i,7) == 0)
i10 = [i10, i];
end
end

subplot(1,2,1)

plot(I,Error, 'o', 'LineWidth', 1.5);
hold on
plot(I(i10), Error(i10), 'ro', 'LineWidth', 1.5);
plot(I(i4), Error(i4), 'go', 'LineWidth', 1.5);
hold off
xlabel('Average Quantumness (I)', extraInputs{:});
ylabel('Test Error', extraInputs{:});


subplot(1,2,2)
load Negativities.mat
plot(neg,Error, 'o', 'LineWidth', 1.5);
hold on
plot(neg(i10), Error(i10), 'ro', 'LineWidth', 1.5);
plot(neg(i4), Error(i4), 'go', 'LineWidth', 1.5);
hold off
xlabel('Wigner Negativity', extraInputs{:});
ylabel('Test Error', extraInputs{:});

legend('d\in\{5,...,9\}','d=10','d=4');

%%

up = median(I);
low = median(I);
q_ind = []; c_ind = [];
for i = 1:length(I)
if(I(i)>=up)
q_ind = [q_ind, i];
end
if(I(i) <= low)
c_ind = [c_ind, i];
end
end

m_q = mean(Error(q_ind)); m_c = mean(Error(c_ind));
s = std(Error(c_ind));
t_stat = (m_c - m_q)/(s/sqrt(length(c_ind)));
%p = normcdf(-t_stat);
[h,p] = ttest2(Error(q_ind), Error(c_ind));
e_max = max(Error);
e_min = min(Error);
xbins=linspace(e_min,e_max,35);
figure;
subplot(2,1,1)
hist(Error(q_ind), xbins);
ylabel('More Quantum', extraInputs{:});

subplot(2,1,2)
hist(2,1,2)
hist(Error(c_ind), xbins);
ylabel('Less Quantum', extraInputs{:});
xlabel('Error', extraInputs{:});
%%
figure();
subplot(1,3,1)
plot(I,Error, 'o', 'LineWidth', 1.5);
hold on
plot(I(i10), Error(i10), 'ro', 'LineWidth', 1.5);
plot(I(i4), Error(i4), 'go', 'LineWidth', 1.5);
hold off
xlabel('Average Quantumness (I)', extraInputs{:});
ylabel('Test Error', extraInputs{:});
legend('d\in\{5,...,9\}','d=10','d=4');
xlim([0,45]);
title('Quantumness vs Error', extraInputs{:})

subplot(1,3,2)
D = 0:139;
D = mod(D,7)+4*ones(1,140);
plot(D, Error, 'o', 'LineWidth', 1.5);
xlim([4,10]);
ylabel('Error', extraInputs{:});
xlabel('Dimension', extraInputs{:});
title('Error vs Dimension', extraInputs{:})


subplot(1,3,3)
plot(D, I, 'ko', 'LineWidth', 1.5);
xlim([4,10]);
ylabel('Quantumness', extraInputs{:});
xlabel('Dimension', extraInputs{:});
title('Quantumness vs Dimension', extraInputs{:})
%%
figure()

subplot(1,2,2)
plot(-neg, I, 'ro', 'LineWidth', 1.5);
xlabel('Initial Wigner Negativity', extraInputs{:});
ylabel('Average Quantumness', extraInputs{:});
title('Wigner Negativity vs Quantumness', extraInputs{:})

subplot(1,2,1)
load Negativities.mat
plot(-neg,Error, 'o', 'LineWidth', 1.5);
hold on
plot(-neg(i10), Error(i10), 'ro', 'LineWidth', 1.5);
plot(-neg(i4), Error(i4), 'go', 'LineWidth', 1.5);
hold off
xlabel('Initial Wigner Negativity', extraInputs{:});
ylabel('Test Error', extraInputs{:});
legend('d\in\{5,...,9\}','d=10','d=4');
title('Wigner Negativity vs Error', extraInputs{:})