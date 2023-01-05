clear
close all
warning('off')
tic
%% Actual values for stiffness updating variables, each alpha represents
% relative change from nominal stiffness parameter.
% 21 measurement channel
npar=15; % No parameters
Ns=20;
nn=5;
mm=5;
for i=1:npar
thetaE(i)=1+0*i;
end
% thetaE(1)=0.8; thetaE(4)=0.9; thetaE(6)=0.8; thetaE(10)=0.8; 
% thetaE(12)=0.7; thetaE(15)=1.2; 
thetaE(1)=0.8; thetaE(3)=0.9; thetaE(6)=0.8; 
thetaE(8)=0.7; thetaE(10)=0.9; thetaE(11)=0.8; 
thetaE(12)=0.7; thetaE(15)=1.2; 
% alpha_act = [0.05; 0.05; -0.05; -0.10; 0.10; -0.15;
%     0.15; -0.05; -0.10; 0.10; -0.20;
%     -0.30;0.60;-0.30;0.60;]; 


%% Load Sensitivity matrix
for R=1:1
    alpha_act = thetaE';
    LoadStructure;
    n_modes = 5; % Number of measured modes
    modeIndex = 1:n_modes; % Indexes of these measured modes

    %% Assemble structure matrices
    structModel.M0 = M0;
    structModel.K0 = K0;
    structModel.K_j = K_j;

%% Simulate "experimental data"
    [psiExpAll,lambdaExp] = eigs(K_act,M0,n_modes,'sm');
    [lambdaExp,dummyInd] = sort((diag(sqrt(lambdaExp)/2/pi)),'ascend') ;
    lambdaExp = lambdaExp(modeIndex);
    psiExpAll = psiExpAll(:,dummyInd(modeIndex));
    psi_m = psiExpAll(measDOFs,:); % 21 DOFs
% Normalize the mode shape vectors as unity norm
for j = 1:n_modes
    psi_m(:,j) = psi_m(:,j) / norm(psi_m(:,j));
end
freq(:,R)=lambdaExp(1:n_modes);
phi(:,:,R)=psi_m(:,1:n_modes);
end
freqmu=freq;
modelmu=phi;

%% add noise 2%
for i=1:Ns
    freqtrue(:,i)=freqmu.*(1+unifrnd(-0.02,0.02,n_modes,1)); %  frequencies-- 
    modeltrue1(:,:,i)=modelmu.*(1+unifrnd(-0.02,0.02,21,n_modes));%  mode shapes
end
%%

modeltrue=remodel(modeltrue1);
s1=var(freqtrue',1);
for j=1:mm
    for i=1:Ns
        aa(:,i)=modeltrue(:,j,i);
    end
    mu=mean(aa')';
    mu=repmat(mu,1,Ns);
    kk=aa-mu;
    for i=1:Ns
        kkk(i)=(norm(kk(:,i)))^2;
    end
    s2(j)=sum(kkk)/Ns;
end

save('freqtrue_3','freqtrue')
save('modeltrue_3','modeltrue')
save('s1_3','s1') 
save('s2_3','s2')
freqmu=mean(freqtrue,2);
for jj=1:nn
    subplot(4,3,jj)
    plot(1:Ns,freqtrue(jj,:)','-bo',[0,Ns],[freqmu(jj) freqmu(jj)],'r--','linewidth',2)
end  

legend('Measurement','Mean')


toc