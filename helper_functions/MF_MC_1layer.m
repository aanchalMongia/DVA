function [ X, U, V] = MF_MC_1layer(y,M,sizeX,rankr   )

% Matrix Completion via Factorization
% min nuclear-norm(X) subject to ||y - M(X)||_2<err

% Inputs
% X - matrix to be estimated
% M - masking operator, applied to vectorized form of X
% y - sampled entries


% Copyright (c) Angshul Majumdar 2010
rng(0);
X = randn(sizeX);
%outsweep = 50%500; for 100k
%alpha =1;%1.5 for 100k 


%R=X_tr>0;
%load('Xbase_1m.mat');%X_tr=reshape( M(y,2) ,sizeX); X=baselinePrediction(X_tr,M, 1,1);

%X=Xbase;


[F,S,T]=svd(X);
%load('svd_Xbase10m.mat');
U = F(:,1:rankr);% *S(1:rankr(1),1:rankr(1)); temp = mldivide(U1,X); [F,S,T] = lansvd(temp,rankr(2),'L');U2 = F(:,1:rankr(2))*S(1:rankr(2),1:rankr(2));V = mldivide(U2,temp);
V=S(1:rankr,1:rankr)*T(:,1:rankr)';


%[F,S,T]=lansvd(X,rankr,'L');
%U = F(:,1:rankr)*S(1:rankr,1:rankr);
%V = mldivide(U,X);


%X=U*V;
x=X(:);
MA=0;
min=10;
cost=[];
global alpha 
outsweep=15;

for out = 1:outsweep
out;
        x = x + (1/alpha)*M(y - M(x,1),2); %x(x<0)=1; no use..for single layer no neg vales in x
        X = reshape(x,sizeX);
         U = mrdivide(X, V);  %U(U<0)=0.0001;
        
         X=U*V;
         x = X(:) + (1/alpha)*M(y - M(X(:),1),2); %x(x<0)=1;
          X = reshape(x,sizeX);
          V = mldivide(U,X);   %V(V<0)=0.0001;
    
          X = U*V;
          x = X(:);
        
       
end

end