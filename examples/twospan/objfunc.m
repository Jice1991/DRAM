function  y=objfunc(x,freqtrue,modeltrue,s1,s2)
Ns=20;nn=8;mm=8;npar=30;

load freqtrue_2span
load modeltrue_2span
load s1_2span
load s2_2span

global H E L P B
H=0.25;E0=2.8e7;L=3.8;P=2300;B=0.15; % structural property
A=B*H;
I=B*H*H*H/12;
LL=ones(1,npar);    
LL(:)=L/npar;
K=zeros(2*(npar+1));
M=zeros(2*(npar+1));
% ModifyM=zeros(2*(npar+1)-2);
% ModifyK=zeros(2*(npar+1)-2);
for R=1:npar
    E=x(R)*E0;
    k=Beam1D2Node_Stiffness(E,I,LL(R));
    K=Beam1D2Node_Assembly(K,k,R,R+1);
    m=Beam1D2Node_Mass(P,A,LL(R));
    M=Beam1D2Node_MS_Assembly(M,m,R,R+1);
end
% for i=1:(2*(npar+1)-3)
%     for n=1:(2*(npar+1)-3)
%         ModifyM(i,n)=M(i+1,n+1);
%         ModifyK(i,n)=K(i+1,n+1);
%     end
% end
% ModifyM(2*(npar+1)-2,1:2*(npar+1)-3)=M(2*(npar+1),2:(2*(npar+1)-2));
% ModifyM(1:2*(npar+1)-3,2*(npar+1)-2)=M(2:(2*(npar+1)-2),2*(npar+1));
% ModifyK(1:2*(npar+1)-3,2*(npar+1)-2)=K(2:(2*(npar+1)-2),2*(npar+1));
% ModifyK(2*(npar+1)-2,1:2*(npar+1)-3)=K(2*(npar+1),2:(2*(npar+1)-2));
% ModifyM(2*(npar+1)-2,2*(npar+1)-2)=M(2*(npar+1),2*(npar+1));
% ModifyK(2*(npar+1)-2,2*(npar+1)-2)=K(2*(npar+1),2*(npar+1));
ModifyM=M;
ModifyK=K;
ModifyM([1,npar+1,2*npar+1],:)=[];
ModifyM(:,[1,npar+1,2*npar+1])=[];
ModifyK([1,npar+1,2*npar+1],:)=[];
ModifyK(:,[1,npar+1,2*npar+1])=[];
% 求解简支梁的固有频率和模态振型；
[V,D]=eig(ModifyK,ModifyM);
F=diag(sqrt(D))./(2*pi);
ModifyV=[V(:,1)./V(1,1) V(:,2)./V(1,2) V(:,3)./V(1,3) V(:,4)./V(1,4) V(:,5)./V(1,5) V(:,6)./V(1,6)...
    V(:,7)./V(1,7) V(:,8)./V(1,8) V(:,9)./V(1,9) V(:,10)./V(1,10)];
% ModifyV1=ModifyV(2:2:(size(ModifyV,1)-1),:);
No = size(ModifyV,1);
ModifyV1=ModifyV([2:2:(No-1)/2,((No+1)/2+1):2:(No-1)],:);

freq1=F(1:nn);
model1=ModifyV1(:,1:mm);
freq1=repmat(freq1,1,Ns);
D1=(freqtrue-freq1).^2; % Original Code

model1=remodel(model1);
for i=1:mm
    for j=1:Ns
        a=modeltrue(:,i,j);
        b=model1(:,i);
        D2(i,j)= (a-b)'*(a-b)/s2(i);
    end
end
 
y=-Ns*log(2*pi)-Ns/2*log(sum(s1))-Ns/2*log(sum(s2))-0.5*(sum(sum(D1')./s1) +  sum(sum(D2'))); 


