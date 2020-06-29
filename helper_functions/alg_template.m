function Yhat=alg_template(Y,predictionMethod,test_ind,left_out)
%alg_template predicts DTIs based on the prediction method selected in
%start.m or sensitivity_analysis.m
%
% INPUT:
%  Y:                     interaction matrix
%  predictionMethod:      method to use for prediction
%  test_indices:          indices of the test set instances
%  left_out:              in case of S1: left_out = test_indices
%                         in case of S2: left_out = left out drugs
%                         in case of S3: left_out = left out targets
%
% OUTPUT:
%  Yhat:                  prediction scores matrix

    % Parameters
    global Sd Sv % fetches Sd and Sv which were defined in getdata() 
        
    predFn = str2func(['alg_'  predictionMethod]);
    Yhat = predFn(Y,Sd,Sv,test_ind,left_out);

end

%MC
function Yhat=alg_mc(Y,~,~,test_ind,~)

%global rank

  W = ones(size(Y));          % weight matrix W
  W(test_ind) = 0;            % set W=0 for test instances
    
IDX = find(W);
M = opRestriction(numel(W),IDX);
y = M(Y(:),1);
%r=rank;
[Yhat] = IST_MC2(y,M,size(W),0,0); %mask changed and lansvd changed, with NN constraint

end

%MF
function Yhat=alg_mf(matrix,~,~,test_ind,~)

global k


Mo = ones(size(matrix));          % weight matrix W
Mo(test_ind) = 0;            % set W=0 for test instances
IDX = find(Mo);
M = opRestriction(numel(Mo),IDX); 
y=M(matrix(:),1);

Yhat=MF_MC_1layer(y,M,size(matrix),k);
end


%DMF
function Yhat=alg_dmf(matrix,~,~,test_ind,~)

global k1 k2 


Mo = ones(size(matrix));          % weight matrix W
Mo(test_ind) = 0;            % set W=0 for test instances
IDX = find(Mo);
M = opRestriction(numel(Mo),IDX); 
y=M(matrix(:),1);

% add to path dmf function in extras/Baseline_dmf(2LMF)

Yhat=dmf(y,M,size(matrix),[k1 k2]);
%Yhat=dgrdmf_3layer(y,M,size(matrix),[k1 k2 k3], mu1,mu2, Lr, Lc);

end

%GRMC
function Yhat=alg_grmc_admm(matrix,Sd,St,test_ind,~)

global lamda mu1 mu2 pp ;

nu1=0.5; nu2=0.5;  

Z=matrix';
Y=matrix;

    Sd = preprocess_PNN(Sd,pp);
    St = preprocess_PNN(St,pp);

    
% Laplacian Matrices    
Dd = diag(sum(Sd)); Lr = Dd - Sd;  
if(det(Dd)==0)
    Dd=0.1*eye(size(Dd))+Dd;
end
Lr = (Dd^(-0.5))*Lr*(Dd^(-0.5));

Dt = diag(sum(St)); Lc = Dt - St;
if(det(Dt)==0)
    Dt=0.1*eye(size(Dt))+Dt;
end
Lc = (Dt^(-0.5))*Lc*(Dt^(-0.5));
    %}
  W = ones(size(matrix));          % weight matrix W
  W(test_ind) = 0;            % set W=0 for test instances
  IDX = find(W);
   
   A = opRestriction(numel(W),IDX); %ii1=sqrt(nu1)*ones(size(matrix));
   A2=opDiag(numel(W),sqrt(nu1));%opMask(ii1(:));%opDirac(numel(W));%opDiag(size(matrix,1), sqrt(nu1));%ii2=sqrt(nu2)*ones(size(matrix));
   A3=opDiag(numel(W),sqrt(nu2));%opMask(ii2(:));%A3= opDiag(numel(matrix), sqrt(nu2));
   AA=opStack(A,A2,  A3  );

   AA_1BMC=[W;W;W];%[W; sqrt(nu1)*eye([size(W,2) size(W,2)]); sqrt(nu2)*eye([size(W,2) size(W,2)]) ];
   
   ys=[]; ys2=[]; ys3=[];
  
for i=1:20
    %YY=AA(matrix(:),1);  
    
    r2=sqrt(nu1).*Z'; r3=sqrt(nu2).*Y;
    YY=[A(matrix(:),1); r2(:); r3(:) ];
    YY_1BMC=[W.*matrix; W.*matrix;  W.*matrix];%[W.*matrix; sqrt(nu1)*Z';  sqrt(nu2)*Y];
    
     
    [X] = IST_MC2( YY , AA,[ 1*size(matrix,1) size(matrix,2)],0,lamda); 
    
   Z=sylvester(nu1*eye(size(Z,1)),mu1*Lr,nu1*X');
   Y=sylvester(nu2*eye(size(Y,1)),mu2*Lc,nu2*X);
   
   ys=[ys; norm(YY-AA(X(:),1),'fro')];
   ys2=[ys2; lamda*norm(X(:),1)];
   ys3=[ys3; mu1*trace(X'*Lr*X)+mu2*trace(X*Lc*X')];
end

  %figure;  plot(ys+ys2+ys3);
  
Yhat=X; 
   
end

%GRBMC
function Xhat=alg_gr1bmc_ppxa(Y,Sd,St,test_ind,~) %named alg_ppxa_emilie earlier


rng(0);

M = ones(size(Y));          % weight matrix or binary mask
M(test_ind) = 0;
  

theta=5;

noi=20;


global pp mu1 mu2 lamda %lamda=theta here, look at equation
%pp=2;
%nu4=1; nu5=1;

%lamda=0.1;

Sd = preprocess_PNN(Sd,pp);
St = preprocess_PNN(St,pp);   
% Laplacian Matrices    
Dd = diag(sum(Sd)); Lr = Dd - Sd;  
if(det(Dd)==0)
    Dd=0.1*eye(size(Dd))+Dd;
end
Lr = (Dd^(-0.5))*Lr*(Dd^(-0.5));

Dt = diag(sum(St)); Lc = Dt - St;
if(det(Dt)==0)
    Dt=0.1*eye(size(Dt))+Dt;
end
Lc = (Dt^(-0.5))*Lc*(Dt^(-0.5));


X1=randn(size(Y));
X2=randn(size(Y));
X3=randn(size(Y));
X4=randn(size(Y));
X5=randn(size(Y));

Xhat2=randn(size(Y));

X_prev = (1/theta)*(X1+X2+X3+X4+X5);

%binary mask
Rop=opRestriction(numel(Y),find(M(:)==1));
R=opToMatrix(Rop); %(dim of R is #elements2keep X  #totalEelementsInX)

 A = eye(size(R,2)) + theta*(R'*R);


   
   for k=1:noi
        
      
       
        b = theta*R'*Rop(Y(:),1) + X1(:);
        tic
        xhat1 = lsqr(A,b);%lsqr(A,b,1e-100,30);%lsqr(A,b);
        %toc
        Xhat1 = reshape(xhat1, size (M));   
        
        
        tic
        %%%% made Xhat1 upate iterative
        
        [U, S, V] = svd(X2);
        sigma = sign(S).*max(0,abs(S)-(theta*lamda) );  
        Xhat2=U*(sigma)*V';
        
        %toc
        tic
        %T=min(max(X3,0),1);
        %Xhat3= (T>=0.5);
        
        Xhat3=min(max(X3,0),1);

        Xhat4 = X4*pinv(2*theta*mu1*Lc + eye(size(Lc,1))); %assuming that Lc is symmetric
        
        %toc
        tic
        %Xhat4 = sylvester(mu1*eye(size(X4,1)),theta*Lc,(1/2)*X4); % 
        
        Xhat5 = pinv(2*theta*mu2*Lr + eye(size(Lr,1)))*X5;
        
        %toc
        %Xhat5=(sylvester(mu2*eye(size(X5,2)),theta*Lr,(1/2)*X5'))';
        
        Xhat=(1/theta)*(Xhat1+Xhat2+Xhat3+Xhat4+Xhat5);
        
        X1=X1+(2*Xhat-X_prev-Xhat1); %   mid var
        X2=X2+(2*Xhat-X_prev-Xhat2);
        X3=X3+(2*Xhat-X_prev-Xhat3);
        X4=X4+(2*Xhat-X_prev-Xhat4);
        X5=X5+(2*Xhat-X_prev-Xhat5);
        
        X_prev = Xhat ; %Xhat1 0.5087
        
        %keep this for debug
       %crit(k) = 1/2*norm(Rop(Y(:),1)-R*Xhat(:))^2 + ...
       %    lamda.*sum(abs(svd(Xhat(:)))) + mu1*trace(Xhat'*Lr*Xhat)+mu2*trace(Xhat*Lc*Xhat');
    end
   
  Xhat= X_prev;
  %plot(crit)%keep this for debug
end

%GRMF
function y3=alg_grmf(Y,Sd,St,test_ind,~)
%alg_grmf predicts DTIs based on the algorithm described in the following paper: 
% Ali Ezzat, Peilin Zhao, Min Wu, Xiao-Li Li and Chee-Keong Kwoh
% (2016) Drug-target interaction prediction with graph-regularized matrix factorization
%
% INPUT:
%  Y:           interaction matrix
%  Sd:          pairwise row similarities matrix
%  St:          pairwise column similarities matrix
%  test_ind:    indices of test set instances
%
% OUTPUT:
%  y3:          prediction matrix
%

    % parameters
    global num_iter p k lambda_l lambda_d lambda_t

    %[Sd,St]=addNewSimilarities(Y, [], Sd,St);
    % preprocessing Sd & St
    % (Sparsification of matrices via p-nearest-neighbor graphs)
    Sd = preprocess_PNN(Sd,p);
     St = preprocess_PNN(St,p);

    % Laplacian Matrices
    Dd = diag(sum(Sd));
    Ld = Dd - Sd;
     Ld = (Dd^(-0.5))*Ld*(Dd^(-0.5));
    Dt = diag(sum(St));
    Lt = Dt - St;
     Lt = (Dt^(-0.5))*Lt*(Dt^(-0.5));

    % (W)GRMF
    [A,B] = initializer(Y,k);	% initialize A & B
    W = ones(size(Y));          % weight matrix W
    W(test_ind) = 0;            % set W=0 for test instances
    [A,B] = alg_grmf_predict(Y,A,B,Ld,Lt,lambda_l,lambda_d,lambda_t,num_iter,W);    % update A & B

    % compute prediction matrix
    y3 = A*B';
end
function [A,B]=alg_grmf_predict(Y,A,B,Ld,Lt,lambda_l,lambda_d,lambda_t,num_iter,W)
%alg_grmf_predict performs alternating least squares for GRMF
%
% INPUT:
%  Y:           interaction matrix
%  A:           drug latent feature matrix
%  B:           target latent feature matrix
%  Ld:          drug graph Laplacian
%  Lt:          target graph Laplacian
%  lambda_ldt:  regularization parameters
%  num_iter:    number of iterations for alternating least squares
%  W:           weight matrix
%
% OUTPUT:
%  A:           updated drug latent feature matrix
%  B:           updated target latent feature matrix
%
    
    K = size(A,2);
    lambda_d_Ld = lambda_d*Ld;          % to avoid 
    lambda_t_Lt = lambda_t*Lt;          % repeated matrix 
    lambda_l_eye_K = lambda_l*eye(K);   % multiplications

    % if no weight matrix is supplied or W is an all-ones matrix...
    if nargin < 10 || isequal(W,ones(size(W)))

        %%%%%%%%%%%%
        %%% GRMF %%%
        %%%%%%%%%%%%

        for z=1:num_iter
            A = (Y*B  - lambda_d_Ld*A) / (B'*B + lambda_l_eye_K);
            B = (Y'*A - lambda_t_Lt*B) / (A'*A + lambda_l_eye_K);
        end
        
        
    else

        %%%%%%%%%%%%%
        %%% WGRMF %%%
        %%%%%%%%%%%%%

        H = W .* Y;
        for z=1:num_iter
%             % for readability...
%             A_old = A;
%             for i=1:size(A,1)
%                 A(i,:) = (H(i,:)*B - lambda_d*Ld(i,:)*A_old) / (B'*diag(W(i,:))*B + lambda*eye(k));
%             end
%             B_old = B;
%             for j=1:size(B,1)
%                 B(j,:) = (H(:,j)'*A - lambda_t*Lt(j,:)*B_old) / (A'*diag(W(:,j))*A + lambda*eye(k));
%             end

            % equivalent, less readable, faster
            A_old = A;
            HB_minus_alpha_Ld_A_old = H*B - lambda_d_Ld*A_old;
            for a=1:size(A,1)
                A(a,:) = HB_minus_alpha_Ld_A_old(a,:) / (B'*diag(W(a,:))*B + lambda_l_eye_K);
            end
            
            B_old = B;
            HtA_minus_beta_Lt_B_old = H'*A - lambda_t_Lt*B_old;
            for b=1:size(B,1)
                B(b,:) = HtA_minus_beta_Lt_B_old(b,:) / (A'*diag(W(:,b))*A + lambda_l_eye_K);
            end
        end
    end
    
end
function [S,p]=preprocess_PNN(S,p)
%preprocess_PNN sparsifies S by keeping, for each drug/target, the "p"
% nearest neighbors (NNs) and discarding the rest. 

    NN_mat = zeros(size(S));
    for j=1:length(NN_mat)
        [~,indx] = sort(S(j,:),'descend');
        indx = indx(1:p+1);     % keep drug/target j and its "p" NNs
        NN_mat(j,indx) = 1;
    end
    NN_mat = (NN_mat+NN_mat')/2;
    S = NN_mat .* S;

end
function [A,B]=initializer(Y,k)
%initializer initializes the A and B latent feature matrices for either
% of the CMF or GRMF algorithms.
%
% INPUT:
%  Y:   interaction matrix
%  k:   number of latent features
%
% OUTPUT:
%  A:   latent feature matrix for drugs
%  B:   latent feature matrix for targets
%

    [u,s,v] = svds(Y,k);
    A = u*(s^0.5);
    B = v*(s^0.5);

%     % Alternative: Use non-negative matrix factorization
%     k = min(k, min(size(Y)));
%     [A,B] = nnmf(Y,k);
%     B = B';

end












    