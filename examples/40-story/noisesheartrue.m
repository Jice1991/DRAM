
%%
%shear building
clear;clc;
global H E L P B
Ns=20; % No. measurement
a=0.0; % coefficient
nn=8; % No. frequency
mm=8;% No. mode shape
npar=40; % no. story
H=0.6;E0=3.3e10;L=80;P=2500;B=0.4;

thetaE=ones(1,npar);
% mu=0.05;
for i=1:npar
thetaE(i)=1+a*i;
end
thetaE(1)=0.8; thetaE(12)=0.9; thetaE(25)=0.8; thetaE(38)=0.8; 
% for i=1:npar
%    true(:,i)=normrnd(thetaE(i),mu,Ns,1);
% end
true=thetaE;
A=B*H; % area
I=B*H*H*H/12; % moment of inertial
LL=ones(1,npar);  
LL(:)=L/npar; % length of each story  

% stiffness  and mass matrix
K=zeros(npar,npar);
k=zeros(npar+1,1);
M=zeros(npar,npar);
for j=1:1
 for R=1:npar
    E=true(j,R)*E0;
    k(R)=shear_Stiffness(E,I,LL(R));
    K=shear_KAssembly(K,k,R);
    m=shear_Mass(P,A,LL(R));
    M=shear_MAssembly(M,m,R);
end

% find eigenvalue and eigenvector
[V,D]=eig(K,M);
F=diag(sqrt(D))./(2*pi);
% normalization
% V1=[V(:,1)./V(1,1) V(:,2)./V(1,2) V(:,3)./V(1,3) V(:,4)./V(1,4)];
V1=[V(:,1)./V(1,1) V(:,2)./V(1,2) V(:,3)./V(1,3) V(:,4)./V(1,4) V(:,5)./V(1,5) V(:,6)./V(1,6)...
    V(:,7)./V(1,7) V(:,8)./V(1,8) V(:,9)./V(1,9) V(:,10)./V(1,10)];
freqtrue(:,j)=F(1:nn);
modeltrue1(:,:,j)=V1(:,1:mm);
end
freqmu=freqtrue;
modelmu=modeltrue1;
% ana = remodel(modelmu);
% test_mode = modelmu([1:2:39],:);
% save('test_mode','test_mode')
%% add noise 2%
for i=1:Ns
    freqtrue(:,i)=freqmu.*(1+unifrnd(-0.0125,0.0125,8,1)); %  frequencies-- 
    modeltrue1(:,:,i)=modelmu.*(1+unifrnd(-0.0125,0.0125,40,8));%  mode shapes
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
% modeltrue = modeltrue([1:5 16:20 26:30 36:40],:,:); % 20 measured DOF 
% save('freqtrue_exp8','freqtrue')
% save('modeltrue_exp8','modeltrue')
% save('s1_exp8','s1') 
% save('s2_exp8','s2')
freqmu=mean(freqtrue,2);
for jj=1:nn
    subplot(4,3,jj)
    plot(1:Ns,freqtrue(jj,:)','-bo',[0,Ns],[freqmu(jj) freqmu(jj)],'r--','linewidth',2)
end  
legend('Measurement','Mean')

