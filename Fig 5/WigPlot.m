%6 Wigner Functions of Random State Evolution
extraInputs = {'interpreter','latex','fontsize',14};
figure();
x = zeros(200,200); p = zeros(200,200);
xvals = linspace(-4,4,200);
for i = 1:200
    for j = 1:200
        x(i,j) = xvals(i);
        p(i,j) = xvals(j);
    end
end

subplot(4,2,1)
load('m1.mat');
mesh(x,p,Wig);
zlim([min(min(Wig)), max(max(Wig))]);
formattedText = {'\fontsize{15}\color{black}\bft=0'; '\fontsize{10}\color{gray}\rmMix'}; 
title(formattedText);

subplot(4,2,2)
load('m2.mat');
mesh(x,p,Wig);
zlim([min(min(Wig)), max(max(Wig))]);
formattedText = {'\fontsize{15}\color{black}\bft=2';'\fontsize{10}\color{gray}\rmMix'}; 
title(formattedText);

subplot(4,2,3)
load('c1.mat');
mesh(x,p,Wig);
zlim([min(min(Wig)), max(max(Wig))]);
formattedText = {'\fontsize{10}\color{gray}\rmCat'};
title(formattedText);


subplot(4,2,4)
load('c2.mat');
mesh(x,p,Wig);
zlim([min(min(Wig)), max(max(Wig))]);
title(formattedText);

subplot(4,2,5)
load('coh1.mat');
mesh(x,p,Wig);
zlim([min(min(Wig)), max(max(Wig))]);
formattedText = {'\fontsize{10}\color{gray}\rmCoherent'};
title(formattedText);

subplot(4,2,6)
load('coh2.mat');
mesh(x,p,Wig);
zlim([min(min(Wig)), max(max(Wig))]);
title(formattedText);


subplot(4,2,7)
load('k1.mat');
mesh(x,p,Wig);
xlabel('X', extraInputs{:})
ylabel('P', extraInputs{:})
zlabel('W', extraInputs{:})
zlim([min(min(Wig)), max(max(Wig))]);
formattedText = {'\fontsize{10}\color{gray}\rmKet 6'};
title(formattedText);


subplot(4,2,8)
load('k2.mat');
mesh(x,p,Wig);
zlim([min(min(Wig)), max(max(Wig))]);
title(formattedText);