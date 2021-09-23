%--------------------------------------------------------------------------
% Usage:  [draws] = StepOne(Prior,Data,Hyp,MCMC)
%--------------------------------------------------------------------------
% Purpose: Run the time series model estimate
%--------------------------------------------------------------------------
% Input: 
% Prior: A structure with prior data and factors
% Data:  A structure with the data on the testing period
% Hyp:   A structure with the model hyperparameters
% MCMC:  A structure with the MCMC settings
%--------------------------------------------------------------------------
% draws: A structure with the model draws for states and parameters
%--------------------------------------------------------------------------
% daniele.bianchi@unibocconi.it
%--------------------------------------------------------------------------


function [MCMC_B, MCMC_BKF, MCMC_K, MCMC_L, MCMC_Q, MCMC_R, MCMC_prob] = StepUnique(Prior,Data,Hyp,MCMC)


V_0          = Hyp.B_0;
a_0          = Hyp.a_0*0.5;
b_0          = Hyp.b_0*0.5;
lambda_0     = Hyp.beta_0;
P            = Data.NInstr;

Yprior       = Prior.Y;
Xprior       = Prior.X;
Ys           = Data.RAReturns;
[N]          = size(Ys,2);
Xs           = Data.X;
[T, K]       = size(Xs);                       % K is the number of factors without the intercept
temp         = eye(K)/(Xprior'*Xprior);        %inv(X'X)         
bOLS         = temp*Xprior'*Yprior;
resid        = Yprior - Xprior*bOLS;
M            = MCMC.Ndraws;
s2OLS        = diag((resid'*resid)/(T - K));
VbOLS        = zeros(K,K,N);
for i=1:N
VbOLS(:,:,i) = s2OLS(i)*temp;
end

Qprior2      = [repmat(Hyp.Qvals2(1,1:K+1),N,1)];
Qprior1      = [Hyp.Qvals1(1,1:K),10];

%--------------------------------------------------------------------------
% Room saving
%--------------------------------------------------------------------------
    
    MCMC_R              =  zeros(T,N,M);         % Returns volatilities
    MCMC_lambda         =  zeros(M,T,N,2);       % Component of the chi2 approximation
    MCMC_Q              =  zeros(M,N,K+1);       % State equations volatilities
    MCMC_K              =  zeros(M,T,N,K+1);     % Breaks factors k
    MCMC_B              =  zeros(M,T,N,K);       % Draws of the betas from the posterior
    MCMC_BKF            =  zeros(M,T,N,K);       %
    pki                 =  zeros(M,N,K);
    pKiR                =  zeros(M,N);
    lambda_draw         =  zeros(T,K-P,M);
    
for vr=1:N
    
%--------------------------------------------------------------------------
% Drawing the starting values from the priors
%--------------------------------------------------------------------------

        MCMC_B(1,:,vr,:)       =  ones(T,1)*reshape(bOLS(:,vr)',1,K);   % Educated prior from the OLS for the betas
        MCMC_lambda(1,:,vr,:)  =  3*ones(T,2);                          % 
        MCMC_K(1,:,vr,:)       =  binornd(1,0.05,T,K+1);                % First drawing of the breaks
        MCMC_Q(1,vr,:)         =  Qprior2(vr,:);                        % Initial value for volatilities
        MCMC_R(:,vr,1)         =  ones(T,1)*s2OLS(vr,:);                % Initial value for returns volatilities
    
        pki(1,vr,:)            =  betarnd(Hyp.p0a,Hyp.p0b,1,K);         % Break probabilities for each factor
        pKiR(1,vr)             =  betarnd(Hyp.p0sva,Hyp.p0svb); 
        
end

sig_inv       =  zeros(T-1,M-1);
sig_inv(:,1)  =  ones(T-1,1).*1/b_0*randgamma(a_0,1);
a = 0;
aR = 0;
for i=1:(M-1)
        
%--------------------------------------------------------------------------
% Main loop
%--------------------------------------------------------------------------
    
    for vr=1:N
        clc;
        display(['Percentage of the loop computed: ', num2str([i/M]*100)])
        display(['Asset number: ', num2str(vr) , ' checking breaks %: ', num2str([a./T aR./T])])

   
        Qs        = squeeze(MCMC_Q(i,vr,1:K));                      % Select the state volatility for the ith MCMC draws
        Ks        = squeeze(MCMC_K(i,:,vr,1:K));                    % Select the values for the ith MCMC draws for the breaks of the factors
        bs        = squeeze(MCMC_B(i,:,vr,:));                      % Select the values for the ith MCMC draws for the betas
        Fbsiconst = [zeros(1,1+K);zeros(K,1),eye(K)];             
        Fbsi      = permute(repmat(Fbsiconst,[1,1,T]),[3,1,2]);        
        
        hbsi      = zeros(T,1,K+1);
        Y_stars   = zeros(T,1);
        Gams      = zeros(T,K+1,K+1);
        
        for qs=1:T;
            hbsi(qs,1,:)  = [1,Xs(qs,:)]; 
            Y_stars(qs,:) = Ys(qs,vr)-Xs(qs,:)*reshape(MCMC_B(i,qs,vr,:),1,K)';
            Gams(qs,:,:)  = [reshape(MCMC_R(qs,vr,i),1,1)^0.5,zeros(1,K);zeros(K,1),diag(Qs.^(1/2))];
        end;  

% -------------------------------------------------------------------------        
% Drawing the breaks for the betas given the rest
% -------------------------------------------------------------------------

        Z      = [Y_stars,bs];
        nbs    = 1;                                               
        K_post = breaks_sampler_aff(Ys(:,vr),zeros(T,1),hbsi,Z,0,zeros(T,size(Z,2)),Fbsi,Gams,[0;bs(1,:)'],[s2OLS(vr),zeros(1,K);zeros(K,1),diag(Qs)],Ks,squeeze(pki(i,vr,:)),nbs);
        
        K_post = squeeze(K_post(nbs,:,:)); 
        MCMC_K(i+1,:,vr,1:K)  =  K_post;                             
        
% -------------------------------------------------------------------------                     
% Updating the posterior breaks probabilities
% -------------------------------------------------------------------------
        
        a                = sum(K_post,1);
        pki(i+1,vr,:)    = betarnd(Hyp.p0a+a,Hyp.p0b+T-a);                 
    
% -------------------------------------------------------------------------
% Drawing the betas conditional on the returns volatilities, the breaks, etc..
% -------------------------------------------------------------------------

        Hkf                 =  hbsi(:,:,2:1+K);             
        P0                  =  eye(K);
        b0                  =  mvnrnd(bOLS(:,vr),squeeze(VbOLS(:,:,vr)))';

        [bt, b_post ,~ ,~]  = forward_filtering_backward_sampler(Ys(:,vr), Hkf, [], [], MCMC_R(:,vr,i), zeros(K,1), eye(K), [], K_post,diag(Qs), b0, P0, 1); 
      
        
        b_post                 = squeeze(b_post);
        MCMC_B(i+1,:,vr,:)     = b_post;                               
        MCMC_BKF(i+1,:,vr,:)   = bt(2:T+1,:);                          

% -------------------------------------------------------------------------
% Drawing the volatilities of the betas given the rest
% -------------------------------------------------------------------------
                
        bb                  = [b0';b_post];
        Y_star              = diff(bb,1,1);  
        sQ                  = Y_star'*Y_star;    
        par1                = 0.5*(size(Y_star,1) + Qprior1(1:K));
        par2                = 0.5*(Qprior2(vr,1:K).*Qprior1(1:K) + diag(sQ)'); 
        
        Qif_inv             = zeros(1,K);
        Qif                 = zeros(1,K);
        
        for j=1:K
        Qif_inv(j)          =     1/par2(j) * randgamma(par1(j),1); %random draw from a gamma (a_post, 1/b_post)
        Qif(j)              =     1/Qif_inv(j);
        end
                
        MCMC_Q(i+1,vr,1:K)   =   Qif;

% -------------------------------------------------------------------------        
% Drawing the breaks for the volatility given the rest
% -------------------------------------------------------------------------
        
         Qs                  = MCMC_Q(i,vr,K+1);                % given the volatility for the sigmat state
         Ks                  = reshape(MCMC_K(i,:,vr,K+1),T,1); % given the state for the volatility dynamics
         Sigma2              = reshape(MCMC_R(:,vr,i),T,1);     % given the returns volatilities 
         Fbsi                = ones(T,1,1);
         hbsi                = ones(T,1);
         Y_sts1              = zeros(T,1);
        
         Y_stscal            = (Ys(:,vr) - sum(Xs.*squeeze(MCMC_B(i+1,:,vr,:)),2)).^2;
         Y_sts1(Y_stscal==0) = 0;
         Y_sts1(Y_stscal~=0) = log(Y_stscal);
         Gams                = ones(T,1)*Qs.^(1/2);

        
        Y_sts                = Y_sts1-squeeze(MCMC_lambda(i,:,vr,1))'; 
        Z                    = log(Sigma2);
        nbs                  = 1;
        
        % initial value of the volatility
        
        LnSigma0             = (normrnd(2,sqrt(10),1,1));     
        
        % drawing the breaks
        
        K_post                  = breaks_sampler_aff_R(Y_sts,zeros(T,1),hbsi,Z,MCMC_lambda(i,:,vr,2)'.^0.5,zeros(T,size(Z,2)),Fbsi,Gams,LnSigma0,Qs,Ks,pKiR(i,vr),nbs);%log(s2OLS)
        MCMC_K(i+1,:,vr,K+1)    = K_post';                       % setting the drawn state for the volatility

% -------------------------------------------------------------------------                     
% Updating the posterior breaks probabilities
% -------------------------------------------------------------------------
                
        aR                    = sum(K_post);
        pKiR(i+1,vr)          = betarnd(Hyp.p0sva+aR,Hyp.p0svb+T-aR);
        
% -------------------------------------------------------------------------
% Drawing the conditional volatilities given the rest
% -------------------------------------------------------------------------
               
        Hkf                 = hbsi(:,1);
        P0                  = 1;
        
        [~, LnSigma2_post ,~ ,~]  = forward_filtering_backward_sampler(Y_sts, Hkf, [], [], MCMC_lambda(i,:,vr,2)', 0, 1, [], K_post', Qs, LnSigma0, P0, 1); 
        
        MCMC_R(:,vr,i+1)          = exp(LnSigma2_post');
        
% -------------------------------------------------------------------------
% Selecting the component of the mixture of ten normals for the vol of vol
% -------------------------------------------------------------------------
        
        [res , ~]                 = mixtures(LnSigma2_post',Y_sts1,'O');
        
        MCMC_lambda(i+1,:,vr,1)   = res(:,1);
        MCMC_lambda(i+1,:,vr,2)   = res(:,2);
        
% -------------------------------------------------------------------------
% Given the component drawing the volatility of log(sigma)
% -------------------------------------------------------------------------

        ll                   = [LnSigma0;LnSigma2_post'];
        Y_star               = diff(ll);
        sQ                   = Y_star'*Y_star;
        
        par1R                = 0.5*(size(Y_star,1)+ Qprior1(1,K));
        par2R                = 0.5*(Qprior2(vr,K+1)*Qprior1(1,K+1)+sQ);
        
        QifR_inv  = 1/par2R * randgamma(par1R,1); %random draw from a gamma (a_post, 1/b_post)
        QifR = 1/QifR_inv;
        
        MCMC_Q(i+1,vr,K+1)      = QifR;
        
        
    end
  
  
% -------------------------------------------------------------------------
% Drawing the risk premia
% -------------------------------------------------------------------------

 for t=1:T-1  
        Y       =   Ys(t+1,:)';
        indvar  =   squeeze(MCMC_B(i+1,t,:,2:end-P));  % - P since the last P elements are instrumental variables
        X       =   [ones(size(Y,1),1), indvar];
        XtX     =   X'*X;
        V_inv   =   eye(K-P)/V_0;
        P_post  =   V_inv+sig_inv(t,i)*XtX;
        V_post  =   eye(K-P)/P_post;

        m_post  =   V_post*(V_inv*lambda_0+sig_inv(t,i)*X'*Y);
           
        lambda_draw(t,:,i)  =  (m_post + V_post*mvnrnd(zeros(size(X,2),1),eye(size(X,2)))')';
        beta                =  squeeze(lambda_draw(t,:,i))';    
        
        a_post  =   a_0 + N;
        b_post  =   b_0 + (Y-X*beta)'*(Y-X*beta); 
        
        % sigma draws 
        
        sig_inv(t,i+1) =   1/b_post * randgamma(a_post,1); %random draw from a gamma (a_post, 1/b_post)
         
  end
%     

end

% -------------------------------------------------------------------------
% Storing the results
% -------------------------------------------------------------------------

MCMC_B    = MCMC_B(MCMC.b_seldr,:,:,:);
MCMC_BKF  = MCMC_BKF(MCMC.b_seldr,:,:,:);   
MCMC_R    = MCMC_R(:,:,MCMC.b_seldr);
MCMC_Q    = MCMC_Q(MCMC.b_seldr,:,:);
MCMC_K    = MCMC_K(MCMC.b_seldr,:,:,:);
MCMC_L    = lambda_draw(:,:,MCMC.b_seldr);
MCMC_prob = pki(MCMC.b_seldr,:,:);
   

   
end

