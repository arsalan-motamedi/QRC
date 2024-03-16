% Animation
extraInputs = {'interpreter','latex','fontsize',14};
myVideo = VideoWriter('GridData'); %open video file
myVideo.FrameRate = 1; 
open(myVideo)

K = ["02", "05", "07", "1", "12"];
kappa = ["02","03","05","1","2","3"];

valK = [0.02,0.05,0.07,0.1,0.12];
valkappa=[0.02,0.03,0.05,0.1,0.2,0.3];

combEr = [];
combQ = [];

%% Plot in a loop and grab frames
for i=1:1:length(K)
    for j=1:length(kappa)
        s = strcat('QuantumnessData-K',K(i),'k',kappa(j),'.mat');
        load(s);
        plot(Quant, Er, 'ob')
        Title = strcat('$K=0.$',K(i),', ','$\kappa=0.$',kappa(j));
        title(Title,extraInputs{:});
        pause(2) %Pause and grab frame
        frame = getframe(gcf); %get frame
        writeVideo(myVideo, frame);
        combEr = [combEr; Er];
        combQ = [combQ;Quant];
    end
end

close(myVideo)

display(corrcoef(combEr,combQ));
save('Data.mat', 'combEr', 'combQ');

%%
figure();
plot(log(combQ),combEr,'ob');
xlabel('$\log$ (Quantumness)', extraInputs{:});
ylabel('Test Error', extraInputs{:});
title('Random State Taining', extraInputs{:});
hold on
for i=1:1:length(K)
    for j=1:length(kappa)
        s = strcat('QuantumnessData-K',K(i),'k',kappa(j),'.mat');
        load(s);
        A = Quant(:);
        B = Er(:);
        %for index=[22]
        %plot(log(A(index)),B(index) ,'*r');
        %end
        %for index=[13]
        %plot(log(A(index)), B(index), '*g');
        %end
    end
end

%% Finding the Wigner functions of the best 10 states

index = findKSmallestIndices(combEr,10);
GoodStates = zeros(10,625);
for c=1:10
    i = mod(index(c,1),7);
    if (i==0)
        i=7;
    end
    j = index(c,2);
    GoodStates(c,:) = mrho(i,j,:);
end
save('GoodStates.mat','GoodStates');

index = findKSmallestIndices(-combEr,10);
NoGoodStates = zeros(10,625);
for c=1:10
    i = mod(index(c,1),7);
    if (i==0)
        i=7;
    end
    j = index(c,2);
    NoGoodStates(c,:) = mrho(i,j,:);
end
save('NoGoodStates.mat','NoGoodStates');
%% (K,kappa) values of best points

x = combQ(:);
y = combEr(:);
N = 15;
a = min(x);
b = max(x);
X = linspace(a,b,N);
meanX = zeros(1,N-1);
for i =1:N-1
    meanX(i) = (X(i)+X(i+1))/2;
end

ParamList = zeros(2,N-1);

for i = 1:N-1
    idx = ListBetween(x, X(i), X(i+1));
    if(isempty(idx))
        ParamList(:,i) = ['empty','empty'];
    else
        [~,ind] = min(y(idx));
        j = idx(ind);
        [j,~] = ind2sub([210,5],j);
        j = (j - mod(j,7))/7;
        j_kappa = mod(j+1,6);
        if(j_kappa==0)
            j_kappa = 6;
        end
        j_K = (j+1-j_kappa)/6+1;
        ParamList(:,i) = [valK(j_K);valkappa(j_kappa)];
    end
end

figure();
plot(log(meanX),ParamList(1,:),'--o');
hold on
plot(log(meanX),ParamList(2,:),'-.o');
hold off
legend('K','kappa')
title('Worst $K$ and $\kappa$ values',extraInputs{:});
xlabel('$\log$(quantumness)',extraInputs{:});
ylabel('$K$ and $\kappa$',extraInputs{:});

function [idx] = ListBetween(z, a, b)
% helper to find values in z that between a,b (a<b)

idx = [];

for i = 1:length(z)
    if(z(i)>=a && z(i)<=b)
        idx = [idx, i];
    end
end
end
