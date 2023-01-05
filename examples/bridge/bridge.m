%%
% <html><a href="../index.html">MCMC toolbox</a> » <a href="../examples.html">Examples</a> » Three modes example</html>

%% 
% Metropolis-Hastings MCMC may have troubles if the target
% distribution has multiple modes. Here we test with a target, which is
% 4 dimensional mixed Gaussian with 3 distinct modes. By changing the
% scaling of DRAM, we find all the modes of this example.
tic
clear model options params

nsimu = 15000; % how many simulations
npar = 15;

model.ssfun= 'objfunc'; % model.ssfun: -2*log(likelihood) function (twice 2)

% DRAM, first large, then small
% tip: adjust qcov, drscale, initial value; initial value is important!
% sensitive to initial values
options.nsimu    = nsimu;
options.method   = 'dram'; % 'dram','am','dr', 'ram' or 'mh'
options.qcov     = eye(npar)*0.05;  % proposal covariance
options.adaptint  = 1000; % interval for adaptation, if 'dram' or 'am' used,DEFAULT adaptint = 200
options.drscale  =  0.1; % propotional to qcov 20~0.01;
% options.adascale = (2.4 / sqrt(npar)) * 1;
for i=1:npar, params{i} = {sprintf('x_{%d}',i), 1, 0.5, 1.5}; end % 0.5 is initial value.[0.5 1.5]

[res,chain] = mcmcrun(model,[],params,options);

%%
figure(1); clf; mcmcplot(chain,[],res,'chainpanel')

figure(2); clf; mcmcplot(chain(0.5*nsimu:end,:),[],res,'denspanel')

figure(3); clf; mcmcplot(chain(0.5*nsimu:end,:),[],res,'hist',[],'normal') % nbins,('normal', 'lognor', or 'kernel')

figure(4); clf; mcmcplot(chain,[],res,'acf')

toc
save('chain')

%% plot
figure(5)
npar = 15;
cm = jet(8);
indx=[1 3 6 8 10 11 12 15];
index=1:1:npar;
index(:,indx)=[];
for i=1:npar
    set(gca,'colororder',cm)
    plot(chain(1:end,indx,1),'-','LineWidth',1);  
    hold on
end
% plot([0,nsimu],[0.9,0.9],'r--','markersize',15,'LineWidth',2);
% plot([0,nsimu],[0.8,0.8],'r--','markersize',15,'LineWidth',2);
% legend('\fontsize{15}\bf\theta_1','\fontsize{15}\bf\theta_1_2'...
%         ,'\fontsize{15}\bf\theta_2_5','\fontsize{15}\bf\theta_3_8')
xlabel('No. iteration','fontsize',20,'fontname','Times');
ylabel('Stiffness parameters','fontsize',20,'fontname','Times');
set(gca,'fontsize',20);
set(gca,'ylim',[0.5 1.5],'ytick',[0.5:0.2:1.5])

figure (6)
plot(chain(1:end,index,1),'b-','LineWidth',1);
% plot([0,nsimu],[0.9,0.9],'r--','markersize',15,'LineWidth',2);
% plot([0,nsimu],[0.8,0.8],'r--','markersize',15,'LineWidth',2);
% legend('\fontsize{15}\bf\theta_1','\fontsize{15}\bf\theta_1_2'...
%         ,'\fontsize{15}\bf\theta_2_5','\fontsize{15}\bf\theta_3_8')
xlabel('No. iteration','fontsize',20,'fontname','Times');
ylabel('Stiffness parameters','fontsize',20,'fontname','Times');
set(gca,'fontsize',20);
set(gca,'ylim',[0.5 1.5],'ytick',[0.5:0.2:1.5])



