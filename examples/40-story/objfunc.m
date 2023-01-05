function  y=objfunc(x,freqtrue,modeltrue,s1,s2)
Ns=20;nn=8;mm=8;npar=40;
load freqtrue_exp8
load modeltrue_exp8
load s1_exp8
load s2_exp8
global H E L P B
H=0.6;E0=3.3e10;L=80;P=2500;B=0.4;
A=B*H; % area
I=B*H*H*H/12; % moment of inertial
LL=ones(1,npar);  
LL(:)=L/npar; % length of each story  

K=zeros(npar,npar);
k=zeros(npar+1,1);
M=zeros(npar,npar);

for R=1:npar
    E=x(R)*E0;   
    k(R)=shear_Stiffness(E,I,LL(R));
    K=shear_KAssembly(K,k,R);
    m=shear_Mass(P,A,LL(R));
    M=shear_MAssembly(M,m,R);
end

% find eigenvalue and eigenvactor
[V,D]=eig(K,M);
F=diag(sqrt(D))./(2*pi);
V1=[V(:,1)./V(1,1) V(:,2)./V(1,2) V(:,3)./V(1,3) V(:,4)./V(1,4) V(:,5)./V(1,5) V(:,6)./V(1,6)...
    V(:,7)./V(1,7) V(:,8)./V(1,8) V(:,9)./V(1,9) V(:,10)./V(1,10)];
% V1=[V(:,1)./V(1,1) V(:,2)./V(1,2) V(:,3)./V(1,3) V(:,4)./V(1,4)];
freq1=F(1:nn);
model1=V1(:,1:mm);
freq1=repmat(freq1,1,Ns);
D1=(freqtrue-freq1).^2; % Original Code
% D1=zeros(nn,Ns);
% for i=1:nn
%      for j=1:Ns
%             D1(i,j)=((freqtrue(i,j)-freq1(i,j))/freqtrue(i,j))^2;
%      end
% end

model1=remodel(model1);
for i=1:mm
    for j=1:Ns
        a=modeltrue(:,i,j);
        b=model1(:,i);
      D2(i,j)= (a-b)'*(a-b)./s2(i); % Original Code
%       D2(i,j)= ((a'*b)^2/((a'*a)*(b'*b)))./s2(i);
    end
end
 
 y=-Ns*log(2*pi)-Ns/2*log(sum(s1))-Ns/2*log(sum(s2))-0.5*(sum(sum(D1')./s1) +  sum(sum(D2')));  % Original Code



