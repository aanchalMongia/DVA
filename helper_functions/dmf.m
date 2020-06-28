
function [Xfinal] = dmf(y,M,sizeX,rankr)

% Matrix Completion via Factorization
% min nuclear-norm(X) subject to ||y - M(X)||_2<err

% Inputs
% X - matrix to be estimated
% M - masking operator, applied to vectorized form of X
% y - sampled entries
rng(0);
X = randn(sizeX);

%x = rand(prod(sizeX),1);
%outsweep = 100;
%alpha =1;


[F,S,T]=svd(X);
%U1 = F(:,1:rankr(1));% *S(1:rankr(1),1:rankr(1)); temp = mldivide(U1,X); [F,S,T] = lansvd(temp,rankr(2),'L');U2 = F(:,1:rankr(2))*S(1:rankr(2),1:rankr(2));V = mldivide(U2,temp);
V1=S(1:rankr(1),1:rankr(1))*T(:,1:rankr(1))';

[F2,S2,T2]=svd(V1);
U2=F2(:,1:rankr(2));
V=S2(1:rankr(2),1:rankr(2))*T2(:,1:rankr(2))';



%X = U1*(U2*V);

%U=randn(sizeX(1),rankr);
%U2=randn()
%V=randn(rankr,sizeX(2));

x=X(:);
%[U1,S1,V1] = lansvd(reshape(x,sizeX),rankr,'L');
%V = S1(1:rankr,1:rankr)*V1(:,1:rankr)';
froErr=0;
min=10000;

global alpha 
outsweep=15;
%apkha=4; outsweep=50;
for out = 1:outsweep
    out;
       x = x + (1/alpha)*M(y - M(x,1),2); %x(x<0)=1;
        X = reshape(x,sizeX);
        U1 = mrdivide(X, U2*V); %U1(U1<0)=0.0001;
        U1=normc(U1);
       
        X = U1*U2*V;
        x = X(:);
        x = x + (1/alpha)*M(y - M(x,1),2); %x(x<0)=1;
        X = reshape(x,sizeX);
        
        %kp= kron(V',U1); %u2=lsqr(kp,X(:),[],30); U2=reshape(u2,size(U2));
        %u2=mldivide(kp,X(:) );
     try
        U2=pinv(U1)*X*pinv(V); %U2(U2<0)=0.0001;
         U2=normc(U2);
        
        
        X = U1*U2*V;
        x = X(:);
        x = x + (1/alpha)*M(y - M(x,1),2); %x(x<0)=1;
        X = reshape(x,sizeX);
          
        V = mldivide(U1*U2,X); %V(V<0)=0.0001;
    
        X = U1*U2*V;
        x = X(:);
        
    %{
    froErr=norm(X-dataXorig, 'fro')/norm(dataXorig,'fro');
      
        if(froErr<min)
            min=froErr;
            Xfinal=X;
            
        end
   %}  
      
        catch
            fprintf('SVD didnt converge!')
        end

end
 Xfinal=X; 