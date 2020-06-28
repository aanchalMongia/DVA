
function [X] = MF_MC_2layer(y,M,sizeX,rankr, mu1,mu2, Lr, Lc )%M_te,full,Xbase,alpha,outsweep)

% Matrix Completion via Doubly Fraph regularized Deep matrix Factorization

% Inputs
% X - matrix to be estimated
% M - masking operator, applied to vectorized form of X
% y - sampled entries


X_tr=reshape( M(y,2) ,sizeX);
X=X_tr;
%X=rand(sizeX);

[F,S,T]=svd(X);
%U1 = F(:,1:rankr(1));% *S(1:rankr(1),1:rankr(1)); temp = mldivide(U1,X); [F,S,T] = lansvd(temp,rankr(2),'L');U2 = F(:,1:rankr(2))*S(1:rankr(2),1:rankr(2));V = mldivide(U2,temp);
V1=S(1:rankr(1),1:rankr(1))*T(:,1:rankr(1))';

[F2,S2,T2]=svd(V1);
U2=F2(:,1:rankr(2));    %U2(U2<0)=0.0001;
V=S2(1:rankr(2),1:rankr(2))*T2(:,1:rankr(2))';    %V(V<0)=0.0001;



MA=0;
min=10;
%outsweep=20;
%alpha=2;
global alpha
%alpha=2;
outsweep=15;
cost1=[]; cost2=[]; cost3=[];

for out = 1:outsweep
out;
    
        x = X(:);
        x = x + (1/alpha)*M(y - M(x,1),2); %x(x<0)=0;    
        B = reshape(x,sizeX);
           
        try
        U1 =sylvester(mu1*Lr, U2*V*(U2*V)'+eye(rankr(1)), B*(U2*V)'); %U1(U1<0)=0.0001;
        %U1=normc(U1);
       
        %X = U1*U2*V;    x = X(:);   x = x + (1/alpha)*M(y - M(x,1),2); %x(x<0)=0;   
        %X = reshape(x,sizeX);
            
        U2=pinv(U1)*X*pinv(V); %U2(U2<0)=0.0001;
        %U2=normc(U2);
            
        %X = U1*U2*V;    x = X(:);   x = x + (1/alpha)*M(y - M(x,1),2); %x(x<0)=0;  
        %X = reshape(x,sizeX);
          
        V =sylvester((U1*U2)'*U1*U2, mu2*Lc, (U1*U2)'*B );   %V(V<0)=0.0001;
    
        X = U1*U2*V;    X(X<0)=0;
        
           
           catch
               fprintf('!')
           end
           
           cost1=[cost1; norm(y-M(X(:),1),'fro')];
           cost2=[cost2; mu1*trace(U1'*Lr*U1)];
           cost3=[cost3; mu2*trace(V*Lc*V')];
end
%figure;  plot(cost1+cost2+cost3); xlabel('No of interations'); ylabel('Cost function');


end
%save('Xfinal.mat','Xfinal');
