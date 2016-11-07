%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*
%*%* modify the reaction norm approach given in Lande 2009 to incorporate maternal effects
%*%* also optionally incorporates selection on the costs of plasticity

%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*
%*%* uses parameters from Lande 2009 and Chevin & Lande 2010 where appropriate
%*%* SOME CONDITIONS FOR THE MODEL ASSUMPTIONS
%*%* phi=(Gbb*delta^2)/(Gaa + Gbb*delta^2) NEAR 1
%*%* sig_xi^2/delta^2 AND gamma*E(sig_z^2) SMALL
%*%*  ...where sig_z^2 = Gaa + 2*Gab*eps_tT + Gbb*exp_tT^2 + sig_e2

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%*%* clear memory, close windows
clear;
close;

%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*
%*%* parameter definitions
tau=0.25;       % fraction of a generation between critical period of development and selection 
rho=0.50;		% environmental auto-correlation
  
xi_bar=0;		% mean of stationary autocorrelated noise
sig_xi=0.001;		% standard deviation of stationary autocorrelated noise
    
omega2=40;		% width of fitness function, 3=> strong stabilizing selection
sigma_e2=1;    	% variance in residual component of phenotypic variation
A=0;			% intercept (elevation) of evolution in reference environment
B=2;			% slope of reaction norm in reference environment
Wmax=1;			% maximal fitness
Gaa=0.1;		% variance of elevation
Gbb=0.045;		% variance of plasticity
mbarinit=0.001;    % initial value of m
Gmm=0.045;        % variance of maternal effect
omega_b_2=100;
omega_m_2=100;
freq=0.5;
Ut=10;			% time at which step change of size delta occurs
epsinit=0;      % initial value of environment
delta=10;		% size of step change/amp of sine wave
nrp=1000000-1;		% length of simulation
%nrp=1000000;
%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*
%*%* define things to receive variables
abar=zeros(nrp,1);
bbar=zeros(nrp,1);
mbar=zeros(nrp,1);
zbar=zeros(nrp,1);
zsbr=zeros(nrp,1);
gaz=zeros(nrp,1);
gbz=zeros(nrp,1);
gmzstar=zeros(nrp,1);
sigma_z2=zeros(nrp,1);
fitness=zeros(nrp,1);
fitvar=zeros(nrp,1);
fitexp=zeros(nrp,1);
epst=zeros(nrp,1);
epstT=zeros(nrp,1);
sigz=zeros(nrp,1);
xi_devel=zeros(nrp,1);
xi_select=zeros(nrp,1);

%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*
%*%* initiate
%%Need to double check these
mbar(1)=mbarinit;
abar(1)=0.001;
bbar(1)=0.001;
zbar(1)=0;
sigma_z2(1)=(Gaa+Gbb*epsinit^2.0+Gmm*epsinit^2.0+sigma_e2);
epstT(1)=epsinit;
epst(1)=epsinit;


GG=[Gaa 0 0; 0 Gbb 0; 0 0 Gmm];


%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*
%*%* generate environmental sequence
  xi_devel(1)=sig_xi*normrnd(0, 1);

  for k=1:nrp
     if(rho==0.0)
     xi_select(k)=sig_xi*normrnd(0,1);
     xi_devel(k)=sig_xi*normrnd(0,1);
     else
    tmp=(rho*xi_devel(k)) + (sig_xi*sqrt(1-(rho^2))*normrnd(0, 1));
    xi_select(k)=tmp;
    j=1/tau;
    if(j>=2)
      %*%* not most computationally efficient yet
      %*%* 1/tau == j must be a positive integer, repeat j-1 times to get AC environment
      for i=1:(j-1)
        tmp=(rho*tmp) + (sig_xi*sqrt(1-(rho^2))*normrnd(0,1));
      end  
    xi_devel(k+1)=(rho*tmp) + (sig_xi*sqrt(1-(rho^2))*normrnd(0,1));
    end
     end
  end

%%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*
%*%*loop - the guts of the code
for t=2:nrp
  %*%*%*%*%*%*%*%* obtain environments during development eps_tT and at selection eps_t
%  eps_t=epsinit+(t>Ut)*delta;
  eps_t=sin(freq*t);
  epst(t)=eps_t;
  eps_t=epst(t)+xi_select(t);
  eps_t1=epst(t-1)+xi_select(t-1);
%  eps_tT=epsinit+((t-tau)>Ut)*delta;
  eps_tT=sin(freq*(t-tau));
  epstT(t)=eps_tT;
  eps_tT=epstT(t) + xi_devel(t);
  eps_tT1=epstT(t-1) + xi_devel(t-1);
    
%%%%%%%%%%%% calculate sigma_z2 %%%%%%%%%%%%%%%%%%
sigma_z2(t)=Gaa+Gbb*eps_tT^2.0+Gmm*eps_tT1^2.0+sigma_e2;
  %*%*%*%*%*%*%*%* calculate gamma
  gamma = 1/(omega2 + sigma_z2(t));
  gamma_b=1/(omega_b_2+Gbb);
  gamma_m=1/(omega_m_2+Gmm);
  
  %*%*%*%*%*%*%*%* calculate change (del) in elevation & plasticity
zbar(t) = abar(t) + bbar(t)*eps_tT +mbar(t)*eps_tT1;
theta = A+B*eps_t;
tmp=zbar(t)-theta;
mat1=-1/omega2*[tmp; tmp*eps_tT; tmp*eps_tT1];
mat3=-1/omega2*[0; omega2*bbar(t)/omega_b_2; omega2*mbar(t)/omega_m_2];
beta = mat1+mat3;
GG=[Gaa 0 0; 0 Gbb 0; 0 0 Gmm];
del  = GG * beta;

%*%*%*%*%*%*%*%* calculate mean fitness
fitness(t) = Wmax*sqrt(gamma*gamma_b*gamma_m*omega2*omega_b_2*omega_m_2)*exp(-(gamma/2)*((zbar(t)-theta)^2)- .5 * gamma_m * (mbar(t)^2) - .5 * gamma_b * (bbar(t)^2));
fitvar(t) = sqrt(gamma*gamma_b*gamma_m*omega2*omega_b_2*omega_m_2);
fitexp(t) = exp(-(gamma/2)*((zbar(t)-theta)^2)- .5 * gamma_m * (mbar(t)^2) - .5 * gamma_b * (bbar(t)^2));		 

%*%*%*%*%*%*%*%* update elevation & plasticity
abar(t+1) = abar(t) + del(1);
bbar(t+1) = bbar(t) + del(2);
mbar(t+1) = mbar(t) + del(3);
    
end

%subplot(3,2,1);semilogx(zbar);xlabel('generations','fontsize',14);ylabel('phenotype/elevation','fontsize',14); axis tight;
%subplot(3,2,2);semilogx(abar);xlabel('generations','fontsize',14);ylabel('additive genetic','fontsize',14); axis tight;
%subplot(3,2,3);semilogx(bbar);xlabel('generations','fontsize',14);ylabel('plasticity','fontsize',14); axis tight;
%subplot(3,2,4);semilogx(mbar);xlabel('generations','fontsize',14);ylabel('maternal effect','fontsize',14); axis tight;
%subplot(3,2,5);semilogx(fitness);xlabel('generations','fontsize',14);ylabel('log(mean fitness)','fontsize',14); axis tight;

subplot(3,2,1);plot(zbar(nrp-99:nrp));xlabel('generations','fontsize',14);ylabel('phenotype/elevation','fontsize',14); axis tight;
subplot(3,2,2);plot(abar(nrp-99:nrp));xlabel('generations','fontsize',14);ylabel('additive genetic','fontsize',14); axis tight;
subplot(3,2,3);plot(bbar(nrp-99:nrp));xlabel('generations','fontsize',14);ylabel('plasticity','fontsize',14); axis tight;
subplot(3,2,4);plot(mbar(nrp-99:nrp));xlabel('generations','fontsize',14);ylabel('maternal effect','fontsize',14); axis tight;
subplot(3,2,5);plot(fitness(nrp-99:nrp));xlabel('generations','fontsize',14);ylabel('log(mean fitness)','fontsize',14); axis tight;

%subplot(2,4,1);semilogx(zbar);xlabel('generations','fontsize',14);ylabel('phenotype/elevation','fontsize',14);axis tight;
%subplot(2,4,2);semilogx(abar);xlabel('generations','fontsize',14);ylabel('additive genetic','fontsize',14);axis tight;
%subplot(2,4,3);semilogx(bbar);xlabel('generations','fontsize',14);ylabel('plasticity','fontsize',14);axis tight;
%subplot(2,4,4);semilogx(mbar);xlabel('generations','fontsize',14);ylabel('maternal effect','fontsize',14);axis tight;
%subplot(2,4,5);semilogx(fitness);xlabel('generations','fontsize',14);ylabel('log(mean fitness)','fontsize',14);axis tight;
%subplot(2,4,6);semilogx(fitvar);xlabel('generations','fontsize',14);ylabel('log(variance factor)','fontsize',14);axis tight;
%subplot(2,4,7);semilogx(fitexp);xlabel('generations','fontsize',14);ylabel('log(adaptation factor)','fontsize',14);axis tight;

%subplot(2,4,1);plot(zbar);xlabel('generations','fontsize',14);ylabel('phenotype/elevation','fontsize',14);
%subplot(2,4,2);plot(abar);xlabel('generations','fontsize',14);ylabel('additive genetic','fontsize',14);
%subplot(2,4,3);plot(bbar);xlabel('generations','fontsize',14);ylabel('plasticity','fontsize',14);
%subplot(2,4,4);plot(mbar);xlabel('generations','fontsize',14);ylabel('maternal effect','fontsize',14);
%subplot(2,4,5);plot(B*epst);xlabel('generations','fontsize',14);ylabel('target','fontsize',14);
%subplot(2,4,6);plot(fitness);xlabel('generations','fontsize',14);ylabel('mean fitness','fontsize',14);
%subplot(2,4,7);plot(fitvar);xlabel('generations','fontsize',14);ylabel('variance factor','fontsize',14);
%subplot(2,4,8);plot(fitexp);xlabel('generations','fontsize',14);ylabel('adaptation factor','fontsize',14);


