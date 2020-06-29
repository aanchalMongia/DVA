function Yhat=alg_template_old(Y,predictionMethod,test_ind,left_out)
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
    global Sd St % fetches Sd and St which were defined in getdata()
    
    
        
    predFn = str2func(['alg_'  predictionMethod]);
    Yhat = predFn(Y,Sd,St,test_ind,left_out);

end

function Xhat=alg_onebmc(Y,Sd,St,test_ind,~)
M = ones(size(Y));          % weight matrix or binary mask
M(test_ind) = 0;
  

global theta1 theta2 theta3
theta=5;
theta1=theta;
theta2=theta;
theta3=theta;
noi=20;


global pp nu4 nu5 lambda %lamda=theta here, look at equation
pp=2;
nu4=1; nu5=1;

lambda=0.1;

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

X_prev=(1/theta)*(X1+X2+X3+X4+X5);

%binary mask
Rop=opRestriction(numel(Y),find(M(:)==1));
R=opToMatrix(Rop); %(dim of R is #elements2keep X  #totalEelementsInX)

%{
I=opRestriction(numel(Y), 1:numel(Y)); % since we need to keep all elements (I*X=X)
I=opToMatrix(I); %(dim of I is #elements2keep=allElements X  #totalEelementsInX, so square Identity matrix of dim numel(X)*numel(X))
%R = ones(size(Y));          % weight matrix W
%R(test_ind) = 0;    

Yvec=Y(:);
X1vec=X1(:);
%}

 A = eye(size(R,2)) + theta*R'*R;


   for k=1:noi
        
       %Xhat1vec=lsqr( [sqrt(theta/2)*R; sqrt(1/2)*I] ,[sqrt(theta/2)*Yvec;sqrt(1/2)*X1vec] ); Xhat1=reshape(Xha1tvec, size(Y))
        %Xhat1=mldivide(sqrt(theta/2)*R; sqrt(1/2)*eye(size(R))] , [sqrt(theta/2)*Y;sqrt(1/2)*X1]);
        
        % you need to solve (I + theta*R’*R)x = theta*R’*y + x_1
       
        b = theta*R'*Rop(Y(:),1) + X1(:);
            %theta*R'*Y(:) + X1(:); %%%% y=Rop(Y(:),1); its same things if u apply to Y or X. Sicne we dont have X we apply to Y
        xhat1 = conjgrad(A,b);%lsqr(A,b,1e-100,30);%lsqr(A,b);
        Xhat1 = reshape(xhat1, size (M));   
        
        %%%% made Xhat1 upate iterative
        
        [U, S, V]=svd(X2);
        sigma= sign(S).*max(0,abs(S)-(theta3/2) ); % ? check
        Xhat2=U*(sigma)*V';
        
        T=min(max(X3,0),1);
        Xhat3= (T>=0.5);
        
        %Xhat3=min(max(X3,0),1);
        
        Xhat4=sylvester(nu4*eye(size(X4,1)),theta1*Lc,(1/2)*X4); % 
        Xhat5=(sylvester(nu5*eye(size(X5,2)),theta2*Lr,(1/2)*X5'))';
        
        Xhat=(1/theta)*(X1+X2+X3+X4+X5);
        
        X1=X1+(2*Xhat-Xhat1-X_prev); %   mid var
        X2=X2+(2*Xhat-Xhat2-X_prev);
        X3=X3+(2*Xhat-Xhat3-X_prev);
        X4=X4+(2*Xhat-Xhat4-X_prev);
        X5=X5+(2*Xhat-Xhat5-X_prev);
        
        X_prev=Xhat5; %Xhat1 0.5087
    end
   
  Xhat=round(X_prev);
end
function Xhat=alg_ppxa_emilie(Y,Sd,St,test_ind,~)
M = ones(size(Y));          % weight matrix or binary mask
M(test_ind) = 0;
  

theta=5;

noi=20;


global pp mu1 mu2 lamda %lamda=theta here, look at equation
%pp=2;
%nu4=1; nu5=1;

%lambda=0.1;

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

%{
I=opRestriction(numel(Y), 1:numel(Y)); % since we need to keep all elements (I*X=X)
I=opToMatrix(I); %(dim of I is #elements2keep=allElements X  #totalEelementsInX, so square Identity matrix of dim numel(X)*numel(X))
%R = ones(size(Y));          % weight matrix W
%R(test_ind) = 0;    

Yvec=Y(:);
X1vec=X1(:);
%}

 A = eye(size(R,2)) + theta*(R'*R);


   
   for k=1:noi
        
       %Xhat1vec=lsqr( [sqrt(theta/2)*R; sqrt(1/2)*I] ,[sqrt(theta/2)*Yvec;sqrt(1/2)*X1vec] ); Xhat1=reshape(Xha1tvec, size(Y))
        %Xhat1=mldivide(sqrt(theta/2)*R; sqrt(1/2)*eye(size(R))] , [sqrt(theta/2)*Y;sqrt(1/2)*X1]);
        
        % you need to solve (I + theta*R’*R)x = theta*R’*y + x_1
       
        b = theta*R'*Rop(Y(:),1) + X1(:);
        tic
        xhat1 = conjgrad(A,b);%lsqr(A,b,1e-100,30);%lsqr(A,b);
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


function Yhat=alg_np(Y,Sd,St,~,~)
%alg_np predicts DTIs based on the Nearest Profile algorithm described in the following paper: 
% Yoshihiro Yamanishi, Michihiro Araki, Alex Gutteridge, Wataru Honda and Minoru Kanehisa,
% (2008) Prediction of drug–target interaction networks from the integration of chemical and genomic spaces

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Nearest Profile (NP) %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Drug
    Sd(logical(eye(length(Sd)))) = 0;   % remove self-similarities
    [maxx, indx] = max(Sd);             % get nearest neighbor for each drug
    for i=1:length(Sd)
        Sd(i, :) = 0;                   % reset all similarities to 0...
        Sd(i, indx(i)) = maxx(i);       % except that of the nearest neighbor
    end
    yd = Sd * Y;

    % Target
    St(logical(eye(length(St)))) = 0;   % remove self-similarities
    [maxx, indx] = max(St);             % get nearest neighbor for each target
    for j=1:length(St)
        St(j, :) = 0;                   % reset all similarities to 0...
        St(j, indx(j)) = maxx(j);       % except that of the nearest neighbor
    end
    yt = (St * Y')';

    % Final Result
    Yhat = (yd + yt) / 2;
end

function Yhat=alg_wp(Y,Sd,St,~,~)
%alg_wp predicts DTIs based on the Weighted Profile algorithm described in the following paper: 
% Yoshihiro Yamanishi, Michihiro Araki, Alex Gutteridge, Wataru Honda and Minoru Kanehisa,
% (2008) Prediction of drug–target interaction networks from the integration of chemical and genomic spaces

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Weighted Profile (WP) %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    yd = bsxfun(@rdivide, Sd * Y, sum(Sd,2));   yd(Y==1) = 1;   % Drug
    yt = bsxfun(@rdivide, Y * St, sum(St));     yt(Y==1) = 1;   % Target
    Yhat = (yd + yt) / 2;
end


function Yhat=alg_rls_wnn(Y,ka,kb,~,~)
%alg_rls_wnn predicts DTIs based on the algorithm described in the following paper: 
% Twan van Laarhoven, Elena Marchiori,
% (2013) Predicting drug–target interactions for new drug compounds using a
%           weighted nearest neighbor profile 
% 
% Code below is adapted from the code available at this website:
% http://cs.ru.nl/~tvanlaarhoven/drugtarget2013/

    %%%%%%%%%%%
    %%% WNN %%%
    %%%%%%%%%%%

    global eta
    %eta = 0.7;     %default
    Y = preprocess_WNN(Y,ka,kb,eta);


    %%%%%%%%%%%
    %%% GIP %%%
    %%%%%%%%%%%

    global alpha
    %alpha = 0.5;   %default
    ka = alpha*ka + (1-alpha)*getGipKernel(Y);
    kb = alpha*kb + (1-alpha)*getGipKernel(Y');


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Regularized Least Squares (RLS-kron) %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    global sigma
	%sigma = 1;     %default
	[va,la] = eig(ka);
	[vb,lb] = eig(kb);
	l = kron(diag(lb)',diag(la));
	l = l ./ (l + sigma);
	m1 = va' * Y * vb;
	m2 = m1 .* l;
	Yhat = va * m2 * vb';
end

function Yhat=alg_nbi(Y,Sd,St,~,~)
%alg_nbi predicts DTIs based on the algorithm described in the following paper: 
% Feixiong Cheng, Chuang Liu, Jing Jiang, Weiqiang Lu, Weihua Li, Guixia Liu, Weixing Zhou, Jin Huang, Yun Tang
% (2012) Prediction of Drug-Target Interactions and Drug Repositioning via Network-Based Inference

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Network-based Inference (NBI) %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % normalize Sd and St:
    Sd = Sd ./ (sum(Sd,2) * sum(Sd));
    St = St ./ (sum(St,2) * sum(St));
    % based on Equation (3) from:
    % Wenhui Wang, Sen Yang, Jing Li
    % (2013) Drug target predictions based on heterogeneous graph inference

    % NBI
    global alpha
    %alpha = 0.5;   %default
    Yhat = Y;
    Yhat = (alpha * Sd * Yhat * St) + ((1 - alpha) * Y);
end


function y3=alg_kbmf2k(Y,Sd,St,~,~)
%alg_kbmf2k predicts DTIs based on the algorithm described in the following paper:
% Mehmet Gönen
% (2012) Predicting drug–target interactions from chemical and genomic kernels using Bayesian matrix factorization

    %--------------------------------------------------------------------

    % parameters
    global rs
    params.R = rs;

    addpath('kbmf2k');

    %--------------------------------------------------------------------
	
	%  where yij = +1 if drug compound di interacts with target protein tj and yij = -1 otherwise
	Y = 2 * (Y > 0.5) - 1;

    % training
    state = alg_kbmf_regression_train(Sd, St, Y, params);

    % predict
    y3 = alg_kbmf_regression_test(Sd, St, state);

end

% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------

% Mehmet Gonen (mehmet.gonen@aalto.fi)
% Helsinki Institute for Information Technology HIIT
% Department of Information and Computer Science
% Aalto University School of Science

function prediction = alg_kbmf_regression_test(Kx, Kz, state)
%     addpath('kbmf2k');
%     addpath('kbmf1mkl1k');
%     addpath('kbmf1k1mkl');
%     addpath('kbmf2mkl');
    prediction = state.parameters.test_function(Kx, Kz, state);
    prediction = prediction.Y.mean;
end

% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------

% Mehmet Gonen (mehmet.gonen@aalto.fi)
% Helsinki Institute for Information Technology HIIT
% Department of Information and Computer Science
% Aalto University School of Science


function state = alg_kbmf_regression_train(Kx, Kz, Y, otherParameters, varargin)
%     addpath('kbmf2k');
%     addpath('kbmf1mkl1k');
%     addpath('kbmf1k1mkl');
%     addpath('kbmf2mkl');


    % external parameters
    if isfield(otherParameters,'R'), parameters.R = otherParameters.R; end
    if isfield(otherParameters,'margin'), parameters.margin = otherParameters.margin; end


    Px = size(Kx, 3);
    Pz = size(Kz, 3);
    %is_supervised = all(~isnan(Y(:)));
    is_supervised = 0;

    parameters.alpha_lambda = 1;
    parameters.beta_lambda = 1;
    if Px > 1 || Pz > 1
        parameters.alpha_eta = 1;
        parameters.beta_eta = 1;
    end
    parameters.iteration = 50;
    parameters.progress = 1;
    parameters.seed = 1606;
    parameters.sigmag = 0.1;
    if Px > 1 || Pz > 1
        parameters.sigmah = 0.1;
    end
    parameters.sigmay = 1.0;

    if is_supervised == 1
        if Px == 1 && Pz == 1
            train_function = @kbmf2k_supervised_regression_variational_train_with_bound;
            test_function = @kbmf2k_supervised_regression_variational_test;
        elseif Px > 1 && Pz == 1
            train_function = @kbmf1mkl1k_supervised_regression_variational_train;
            test_function = @kbmf1mkl1k_supervised_regression_variational_test;
        elseif Px == 1 && Pz > 1
            train_function = @kbmf1k1mkl_supervised_regression_variational_train;
            test_function = @kbmf1k1mkl_supervised_regression_variational_train;
        elseif Px > 1 && Pz > 1
            train_function = @kbmf2mkl_supervised_regression_variational_train;
            test_function = @kbmf2mkl_supervised_regression_variational_test;
        end
    else
        if Px == 1 && Pz == 1
            train_function = @kbmf2k_semisupervised_regression_variational_train;
            test_function = @kbmf2k_semisupervised_regression_variational_test;
        elseif Px > 1 && Pz == 1
            train_function = @kbmf1mkl1k_semisupervised_regression_variational_train;
            test_function = @kbmf1mkl1k_semisupervised_regression_variational_test;
        elseif Px == 1 && Pz > 1
            train_function = @kbmf1k1mkl_semisupervised_regression_variational_train;
            test_function = @kbmf1k1mkl_semisupervised_regression_variational_test;
        elseif Px > 1 && Pz > 1
            train_function = @kbmf2mkl_semisupervised_regression_variational_train;
            test_function = @kbmf2mkl_semisupervised_regression_variational_test;
        end
    end

    for i = 1:2:nargin - 4
        parameters.(varargin{i}) = varargin{i + 1};
    end
    
    parameters.train_function = train_function;
    parameters.test_function = test_function;

    state = train_function(Kx, Kz, Y, parameters);
end


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
[Sd,St]=addNewSimilarities(Y, [], Sd,St);
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

function y3=alg_cmf(Y,Sd,St,test_ind,~)
%alg_cmf predicts DTIs based on the algorithm described in the following paper:
% Xiaodong Zheng, Hao Ding, Hiroshi Mamitsuka and Shanfeng Zhu
% (2013) Collaborative Matrix Factorization with Multiple Similarities for Predicting Drug-Target Interactions
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
    global num_iter k lambda_l lambda_d lambda_t

    [A,B] = initializer(Y,k);	% initialize A & B
    W = ones(size(Y));          % weight matrix W
    W(test_ind) = 0;            % set W=0 for test instances
    [A,B] = alg_cmf_predict(Y,A,B,Sd,St,lambda_l,lambda_d,lambda_t,num_iter,W);     % update A & B

    % compute prediction matrix
    y3 = A*B';

end

function [A,B]=alg_cmf_predict(Y,A,B,Sd,St,lambda_l,lambda_d,lambda_t,num_iter,W)
%alg_cmf_predict performs alternating least squares for CMF
%
% INPUT:
%  Y:           interaction matrix
%  A:           drug latent feature matrix
%  B:           target latent feature matrix
%  Sd:          pairwise row similarities matrix
%  St:          pairwise column similarities matrix
%  lambda_ldt:  regularization parameters
%  num_iter:    number of iterations for alternating least squares
%  W:           weight matrix
%
% OUTPUT:
%  A:           updated drug latent feature matrix
%  B:           updated target latent feature matrix
%
    
    K = size(A,2);
    lambda_d_Sd = lambda_d*Sd;          % to avoid 
    lambda_t_St = lambda_t*St;          % repeated matrix 
    lambda_l_eye_K = lambda_l*eye(K);   % multiplications

    % if no weight matrix is supplied or W is an all-ones matrix...
    if nargin < 10 || isequal(W,ones(size(W)))
        AtA = A'*A;
        BtB = B'*B;
        for z=1:num_iter
            A = (Y*B + lambda_d_Sd*A)  / (BtB + lambda_l_eye_K + lambda_d*(AtA));
            AtA = A'*A;
            B = (Y'*A + lambda_t_St*B) / (AtA + lambda_l_eye_K + lambda_t*(BtB));
            BtB = B'*B;
        end
        
    else
        H = W .* Y;
        for z=1:num_iter
%             % for readability...
%             A_old = A;
%             lambda_d_A_oldt_A_old = lambda_d*(A_old'*A_old);
%             for a=1:size(A,1)
%                 A(a,:) = (H(a,:)*B + lambda_d_Sd(a,:)*A_old) / (B'*B + lambda_l_eye_k + lambda_d_A_oldt_A_old);
%             end
%             B_old = B;
%             lambda_t_B_oldt_B_old = lambda_t*(B_old'*B_old);
%             for b=1:size(B,1)
%                 B(b,:) = (H(:,b)'*A + lambda_t_St(b,:)*B_old) / (A'*A + lambda_l_eye_k + lambda_t_B_oldt_B_old);
%             end

            % equivalent, less readable, faster
            A_old = A;
            HB_plus_lambda_d_Sd_A_old = H*B + lambda_d_Sd*A_old;
            lambda_l_eye_k_plus_lambda_d_A_oldt_A_old = lambda_l_eye_K + lambda_d*(A_old'*A_old);
            for a=1:size(A,1)
                A(a,:) = HB_plus_lambda_d_Sd_A_old(a,:) / (B'*diag(W(a,:))*B + lambda_l_eye_k_plus_lambda_d_A_oldt_A_old);
            end
            B_old = B;
            HtA_plus_lambda_t_St_B_old = H'*A + lambda_t_St*B_old;
            lambda_l_eye_k_plus_lambda_t_B_oldt_B_old = lambda_l_eye_K + lambda_t*(B_old'*B_old);
            for b=1:size(B,1)
                B(b,:) = HtA_plus_lambda_t_St_B_old(b,:) / (A'*diag(W(:,b))*A + lambda_l_eye_k_plus_lambda_t_B_oldt_B_old);
            end

        end
    end
    
end


function Yhat=alg_sitar(Y,Sd,St,test_ind,~)
%alg_nbi predicts DTIs based on the algorithm described in the following paper: 
% Feixiong Cheng, Chuang Liu, Jing Jiang, Weiqiang Lu, Weihua Li, Guixia Liu, Weixing Zhou, Jin Huang, Yun Tang
% (2012) Prediction of Drug-Target Interactions and Drug Repositioning via Network-Based Inference
%
% INPUT:
%  Y:           interaction matrix
%  Sd:          pairwise row similarities matrix
%  St:          pairwise column similarities matrix
%  test_ind:    indices of test set instances
%
% OUTPUT:
%  Yhat:        prediction matrix
%

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Similarity-based Inference of drug-TARgets (SITAR) %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % +ve set
    pos_indices = find(Y == 1);     % +ve instances
    pos_indices = pos_indices(:);   % force into column vector

    % generate feature vectors
    r = 0.5;
    [ds, ts] = ind2sub(size(Y), 1:numel(Y));
    [pds, pts] = ind2sub(size(Y), pos_indices);
    feat_vectors = (Sd(ds, pds).^r) .* (St(ts, pts).^(1-r));
%     feat_vectors = zeros(numel(Y), length(pos_indices));
%     for i=1:numel(Y)
%         d = ds(i);
%         t = ts(i);
%         feat_vectors(i,:) = (Sd(d, pds).^r) .* (St(t, pts).^(1-r));
%     end

    % train
    train_ind = 1:numel(Y);
    train_ind(test_ind) = [];
    train_ind = train_ind(:);
    model = compact(fitcsvm(feat_vectors(train_ind,:), Y(train_ind), 'KernelFunction', 'rbf'));
    %model = compact(fitcsvm(feat_vectors(train_ind,:), Y(train_ind), 'KernelFunction', 'rbf', 'BoxConstraint', 10));
    %model = compact(TreeBagger(500, feat_vectors(train_ind,:), Y(train_ind)));

    % predict
    [~, scores] = predict(model, feat_vectors(test_ind,:));
    Yhat = Y;
    Yhat(test_ind) = scores(:,2);

end

function Yhat=alg_laprls(Y,Sd,St,~,~)
%alg_laprls predicts DTIs based on the algorithm described in the following paper: 
% Zheng Xia, Ling-Yun Wu, Xiaobo Zhou, Stephen TC Wong,
% (2010) Semi-supervised drug-protein interaction prediction from heterogeneous biological spaces
% 
% Code adapted from supplementary material of Laarhoven 2011
%
% INPUT:
%  Y:           interaction matrix
%  Sd:          pairwise row similarities matrix
%  St:          pairwise column similarities matrix
%
% OUTPUT:
%  Yhat:        prediction matrix
%

    % Parameters as per the above paper
    ga1 = 1;
    gb1 = 1;
    ga2 = 0.01;
    gb2 = 0.01;
    ba  = 0.3;
    bb  = 0.3;

    sa = Sd;
    sb = St;
    ka = Y*Y';
    kb = Y'*Y;
    wa = (ga1*sa + ga2*ka) / (ga1+ga2);
    wb = (gb1*sb + gb2*kb) / (gb1+gb2);

    da = diag(sum(wa));
    db = diag(sum(wb));
    la = da^-0.5 * (da - wa) * da^-0.5;
    lb = db^-0.5 * (db - wb) * db^-0.5;

    fa = wa / (wa + ba*la*wa) * Y;
    fb = wb / (wb + bb*lb*wb) * Y';

    Yhat = (fa+fb') / 2;
end

function Yhat=alg_mc(Y,Sd,St,test_ind,~)

global rank

  W = ones(size(Y));          % weight matrix W
  W(test_ind) = 0;            % set W=0 for test instances
    
IDX = find(W);
M = opRestriction(numel(W),IDX);
y = M(Y(:),1);
%r=rank;
[Yhat] = IST_MC2(y,M,size(W),0,0); %mask changed and lansvd changed, with NN constraint

end

function Yhat=alg_grmc2(matrix,Sd,St,test_ind,~)

[Sd,St]=addNewSimilarities(matrix, [], Sd,St);


%[~,Sd,St,~,~]=getdata('nr','cosine','single', test_ind);
%Sd(Sd==NaN)=0;
%St(St==NaN)=0;


global lamda nu1 nu2 mu1 mu2 pp method;
%nu1=0.5; nu2=0.5; mu1=0.5; mu2=0.5; lambda=0.5;

Z=matrix';
Y=matrix;

% Laplacian Matrices
%Lr = (single(full(  laplacian(graph(Sd,'upper','omitselfloops')))));
%Lc = (single(full(  laplacian(graph(St,'upper','omitselfloops')))));
    
    % preprocessing Sd & St
    % (Sparsification of matrices via p-nearest-neighbor graphs)
    
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
    
     if strcmp(method,'NNM')
         [X] = IST_MC2( YY , AA,[ 1*size(matrix,1) size(matrix,2)],0,lamda); 
    
     elseif strcmp(method,'1BMC')
        %temp=ones(size(matrix));t2=AA(temp(:),1);AAnew=AA(t2,2);
        
        [X] = OBMC( matrix ,W, [ 1*size(matrix,1) size(matrix,2)]); 
     end
   Z=sylvester(nu1*eye(size(Z,1)),mu1*Lr,nu1*X');
   Y=sylvester(nu2*eye(size(Y,1)),mu2*Lc,nu2*X);
   
   ys=[ys; norm(YY-AA(X(:),1),'fro')];
   ys2=[ys2; lamda*norm(X(:),1)];
   ys3=[ys3; mu1*trace(X'*Lr*X)+mu2*trace(X*Lc*X')];
end

  %figure;  plot(ys+ys2+ys3);
  
Yhat=X; 
   
end

function Yhat=alg_admm(matrix,Sd,St,test_ind)

%[~,Sd,St,~,~]=getdata('nr','cosine','single', test_ind);
%Sd(Sd==NaN)=0;
%St(St==NaN)=0;


global lamda mu1 mu2 pp ;

nu1=0.5; nu2=0.5;  

Z=matrix';
Y=matrix;

% Laplacian Matrices
%Lr = (single(full(  laplacian(graph(Sd,'upper','omitselfloops')))));
%Lc = (single(full(  laplacian(graph(St,'upper','omitselfloops')))));
    
    % preprocessing Sd & St
    % (Sparsification of matrices via p-nearest-neighbor graphs)
    
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



function Yhat=alg_dgrdmf(matrix,Sd,St,test_ind,left_out)

%[Sd,St]=addNewSimilarities(matrix, left_out, Sd,St);


global k1 k2 k3 mu1 mu2 pp;


Mo = ones(size(matrix));          % weight matrix W
Mo(test_ind) = 0;            % set W=0 for test instances
IDX = find(Mo);
M = opRestriction(numel(Mo),IDX); 
y=M(matrix(:),1);

%pp=2;
% preprocessing Sd & St (Sparsification of matrices via p-nearest-neighbor graphs)
Sd = preprocess_PNN(Sd,pp);
St = preprocess_PNN(St,pp);

% Laplacian Matrices    
Dd = diag(sum(Sd)); Lr = Dd - Sd;  
if(det(Dd)==0) Dd=0.1*eye(size(Dd))+Dd; end
Lr = (Dd^(-0.5))*Lr*(Dd^(-0.5));

Dt = diag(sum(St)); Lc = Dt - St;
if(det(Dt)==0)  Dt=0.1*eye(size(Dt))+Dt; end
Lc = (Dt^(-0.5))*Lc*(Dt^(-0.5));

Yhat=dgrdmf_2layer(y,M,size(matrix),[k1 k2], mu1,mu2, Lr, Lc);
%Yhat=dgrdmf_3layer(y,M,size(matrix),[k1 k2 k3], mu1,mu2, Lr, Lc);

end

% same as dgrdmf just name is different so that diff parameters are fetched
function Yhat=alg_mgrdmf(matrix,Sd,St,test_ind,left_out)

[Sd,St]=addNewSimilarities(matrix, left_out, Sd,St);
      
global k1 k2 k3 mu1 mu2 pp;

Mo = ones(size(matrix));          % weight matrix W
Mo(test_ind) = 0;            % set W=0 for test instances
IDX = find(Mo);
M = opRestriction(numel(Mo),IDX); 
y=M(matrix(:),1);

%pp=2;
% preprocessing Sd & St (Sparsification of matrices via p-nearest-neighbor graphs)
Sd = preprocess_PNN(Sd,pp);
St = preprocess_PNN(St,pp);

% Laplacian Matrices    
Dd = diag(sum(Sd)); Lr = Dd - Sd;  
if(det(Dd)==0) Dd=0.1*eye(size(Dd))+Dd; end
Lr = (Dd^(-0.5))*Lr*(Dd^(-0.5));

Dt = diag(sum(St)); Lc = Dt - St;
if(det(Dt)==0)  Dt=0.1*eye(size(Dt))+Dt; end
Lc = (Dt^(-0.5))*Lc*(Dt^(-0.5));

Yhat=dgrdmf_2layer(y,M,size(matrix),[k1 k2], mu1,mu2, Lr, Lc);
%Yhat=dgrdmf_3layer(y,M,size(matrix),[k1 k2 k3], mu1,mu2, Lr, Lc);

end

function Yhat=alg_mf(matrix,Sd,St,test_ind,~)

global k


Mo = ones(size(matrix));          % weight matrix W
Mo(test_ind) = 0;            % set W=0 for test instances
IDX = find(Mo);
M = opRestriction(numel(Mo),IDX); 
y=M(matrix(:),1);

Yhat=MF_MC_1layer(y,M,size(matrix),k);
end

function Yhat=alg_dmf(matrix,Sd,St,test_ind,~)

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

function Yhat=alg_gr_rbm(Y,Sd,St,test_ind,~)

python_path=py.sys.path;
module_path='D:\AanchalMongia_phdClg\Phd\Python_workspace\RBM\';
if count(python_path,module_path)==0
    insert(python_path, int32(0), module_path);
end

global nh ne lamda1 lamda2 lamda

%returned_obj=py.RBM_DTI.sample();
A=randn(6,8);
returned_obj2=py.RBM_DTI_pv_noReg.train(Y(:)',uint64(size(Y,1)),  uint64(nh),uint64(ne), uint64(lamda1), uint64(lamda2), uint64(lamda), Sd(:)', A(:)');
data = double(py.array.array('d',py.numpy.nditer(returned_obj2))); %d is for double, see link below on types
Yhat = reshape(data,[size(Y,1) size(Y,2)]);
%system(['py -3 ' module_path 'RBM_CF.py'])
 
%{

  W = ones(size(Y));          % weight matrix W
  W(test_ind) = 0;            % set W=0 for test instances
    
IDX = find(W);
M = opRestriction(numel(W),IDX);
y = M(Y(:),1);
r=rank;
[Yhat] = IST_MC(y,M,size(W),r); %mask changed and lansvd changed, with NN constraint
%}
end


function Yhat=alg_ensemble(Y,Sd,St,test_indices,~)
%alg_ensemble predicts DTIs based on the algorithm described in the
%following paper (but without the dimensionality reduction):
% Ali Ezzat, Peilin Zhao, Min Wu, Xiao-Li Li and Chee-Keong Kwoh
% (2017) Drug-Target Interaction Prediction using Ensemble Learning and
%           Dimensionality Reduction 
%
% INPUT:
%  Y:                     interaction matrix
%  test_indices:          indices of the test set instances
%
% OUTPUT:
%  Yhat:                  prediction matrix (including test set scores)
 
    % Parameters
    global numLearners drugFeatureVectors targetFeatureVectors r
 
    % instances to be excluded from training set
    exclude_indices = test_indices;
 
    % generate models
    Yhat = 0;
    for c = 1:numLearners
        % feature subspacing
        numDrugFeatures = size(drugFeatureVectors, 2);
        numTargetFeatures = size(targetFeatureVectors, 2);
        drugFeatures = randperm(numDrugFeatures, floor(numDrugFeatures*r));
        targetFeatures = randperm(numTargetFeatures, floor(numTargetFeatures*r));
        drugFeatVectors = drugFeatureVectors(:, drugFeatures);
        targetFeatVectors = targetFeatureVectors(:, targetFeatures);
 
        % train
        [patterns, labels, ~] = generateTrainingSet(Y, exclude_indices, drugFeatVectors, targetFeatVectors);
        predModel.model = compact(fitctree(patterns, labels, 'Prior', 'uniform'));
        Yhat = Yhat + predictor(predModel, test_indices, 'ensemdt', drugFeatVectors, targetFeatVectors);
    end
    clear patterns labels drugFeatVectors targetFeatVectors
end

function Yhat=predictor(predModel,test_indices,algorithm,drugFeatVectors,targetFeatVectors)
%predictor is a helper function that produces prediction scores for
%feature-based DTI prediction methods

    global batchSize

    % to hold prediction scores
    predScores = zeros(length(test_indices), 1);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % NOTE:                                                            %
    % the entire set of feature vectors of the test set can't fit into %
    % memory, so we do predictions of test set instances in batches    %
    % (default batch size = 10000)                                     %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % collect predictions
    for counter = 0:batchSize:length(test_indices)
        % current batch of test set instances
        start_index = counter + 1;
        end_index = min(counter + batchSize, length(test_indices));
        if start_index > length(test_indices), break; end

        % prepare testing instances
        testingFeatureVectors = generateFeatures(test_indices(start_index:end_index), drugFeatVectors, targetFeatVectors);

        % get prediction scores
        predScores(start_index:end_index) = predictWithModel(algorithm, predModel, testingFeatureVectors);

        clear testingFeatureVectors
    end

    Yhat = zeros(size(drugFeatVectors, 1), size(targetFeatVectors, 1));
    Yhat(test_indices) = predScores;
end

function scores=predictWithModel(algorithm,predModel,X)

    switch algorithm
        % DECISION TREE, RANDOM FOREST, SVM ----------
        case {'dt', 'rf', 'svm'}
            % prediction scores (probabilities)
            [~, scores] = predict(predModel.model, X);
            scores = scores(:,2);

        % OUR METHOD ---------------------------------
        otherwise
            % predicted labels
            scores = predict(predModel.model, X);
            if iscell(scores), scores = str2double(scores); end
    end
end

function Yhat=alg_grmc(Y,Sd,St,~,~)
%Y=Y';

global lambda_l lambda_d lambda_t


% Keep 20% for training, the rest for validation
[y_train, mask_train, y_val, mask_val, y_test, mask_test] = split_observed(Y, [.9, .1, 0]);
params.size_X = size(Y);

%% Normalize data to zero mean and keep the linear transformation details
y_lims_init = [min(y_train), max(y_train)];

mean_train = mean(y_train);

y_train = y_train - mean_train;
y_val = y_val - mean_train;
%y_test = y_test - mean_train;

y_lims_scaled = [min(y_train), max(y_train)];

prob_params.Lr = (single(full(  laplacian(graph(Sd,'upper','omitselfloops')))));
prob_params.Lc = (single(full(  laplacian(graph(St,'upper','omitselfloops')))));

% DATASETS and masks:
prob_params.size_X = params.size_X;
prob_params.mask_val = mask_val;
prob_params.mask_test = mask_test;
prob_params.A_op = @(x) sample_sparse(x, mask_train);
prob_params.At_op = @(x) sample_sparse_t(x, mask_train);
prob_params.AtA_op = @(x) sample_sparse_AtA(x, mask_train);


%% SOLVER PARAMETERS 
solver_params.maxit = 100;
solver_params.verbose = 3;
solver_params.tol_abs = 2e-6;
solver_params.tol_rel = 1e-5;

% need the scaling used for preprocessing to calculate error correctly
solver_params.y_lims_init = y_lims_init;
solver_params.y_lims_scaled = y_lims_scaled;

% MOST IMPORTANT: use verbose = 1 to set rho accordingly (depends on tolerances)
%solver_params.rho_ADMM = .005000;
%solver_params.rho_ADMM = .2 * geomean([max(1e-3,prob_params.gamma_n), geomean([max(1e-3,norm(y_train)), max(1e-3,prob_params.gamma_r), max(1e-3,prob_params.gamma_c)])]);

% for small matrices use false!
solver_params.svds = false;

%%%%% %% disp('With graph regularization')
%% Solve the problem using graphs
prob_params.gamma_n = lambda_l;%.01;
prob_params.gamma_r = lambda_d;%.003;
prob_params.gamma_c = lambda_t;%.003;
solver_params.rho_ADMM = .009;

[X_MC_graphs, stat_MC_graphs] = MC_solve_ADMM(y_train, y_val, y_test, prob_params, solver_params);

Yhat=X_MC_graphs;
end




function Yhat=alg_TMF(Y,Sd,St,test_indices, left_out)%test_indices,~)
%left_out=[5;6;8;13;34;51]; %s2nr
%left_out=[];%s3gpcr
global ds cv_setting

%dataset='NR'%'GPCR'%'NR'%load(['..\' dataset '.mat'])%Scenario = cv_setting%'S1';

if(ds==4 )
    nComp=8;
    if(strcmp( cv_setting,'S1') )
        nComp=4;
    end
elseif( ds==1 )
    nComp=50;
else
    nComp=18;
end
%nCV = n;%10;   
d_s=Sd; t_s=St;

Adj=Y;%DTI;
Sim_row=(d_s+d_s')/2;
Sim_col=(t_s+t_s')/2;

ratio_=nComp;
Debug=true;

fixRandom = Debug;
if fixRandom
    %disp('Debug')
    rand('state',1234567890); % fix random seed for debug
end

if strcmp(cv_setting,'S1')
Parameters.lamda_ar = 1;         Parameters.lamda_ac = 1;
Parameters.lamda_r =.5;           Parameters.lamda_c = .5; 
%PerformS1(Adj, Sim_row,Sim_col, nCV,ratio_,Parameters);
if ratio_ >1
    latent_dim = ratio_; % directly use ratio_ as latent dim
    %disp('directly use ratio_ as the latent dim')
end
[U_r,S_r,~]=TurnSimBySVD(Sim_row);
[U_c,S_c,~]=TurnSimBySVD(Sim_col);
 Observed_Flag = ones(size(Adj) ); %contains the combination of S2, S3, S4 block
 Observed_Flag(test_indices) = 0;
 %Observed_Flag(ID_Tst) = 0;   Observed_Flag(ID_Tst_neg) = 0;
 TempMat =  Adj;
 TempMat(test_indices)=0;
 %TempMat (ID_Tst)  = 0;
 [trn_row_s1,trn_col_s1] = SplitOutPureS1Scenarios(TempMat);
 if ratio_<=1 % using the low-rank
        latent_dim=fix(rank(TempMat)* ratio_);
 end
 clear TempMat;
 LabelMat= Adj; 
 % ONE SHULD NOT BE FEEDING TR DATA TO GetPrior function prior=fraction of
 % 1s
 Prior = GetPrior(Adj); 
LabelMat(test_indices)  = Prior; 
 %LabelMat(ID_Tst)  = Prior;      LabelMat(ID_Tst_neg)  = Prior;        
Predicted_LabelMat =PredictS1_Optim...
        ( LabelMat (trn_row_s1,trn_col_s1),Observed_Flag (trn_row_s1,trn_col_s1) ,...
        U_r,S_r, U_c,S_c, ...
        trn_row_s1, trn_col_s1,trn_row_s1, trn_col_s1, latent_dim,Parameters);
    
 final=zeros(size(Y)); final(trn_row_s1,trn_col_s1)=Predicted_LabelMat;
    
    %{
%test_indices was merged from ID_Tst and ID_Tst_neg in nfold.m file (cv setting 1)
 ID_Tst=test_indices(1:9);
    ID_Tst_neg=test_indices(10:end);
 
    MARK_POS = +9;
    MARK_NEG = -9;

    % 5 evaluation for only pure S1
     MARK_mat = zeros(size(Adj) ); 
    MARK_mat(ID_Tst) = MARK_POS;
    MARK_mat(ID_Tst_neg) = MARK_NEG;

    %idx_pos = find(MARK_mat(trn_row_s1,trn_col_s1) ==MARK_POS);
    %idx_neg = find(MARK_mat(trn_row_s1,trn_col_s1) ==MARK_NEG);
    
    %another way of evaluation: all test indices (not neccesarilty pure
    ones): gives same auc n aupr as it should
    idx_pos = find(MARK_mat ==MARK_POS);
    idx_neg = find(MARK_mat ==MARK_NEG);
    Predicted_LabelMat=final;
    
    TrueScore = Predicted_LabelMat(idx_pos);
    FalseScore = Predicted_LabelMat(idx_neg);
    
    EstimationAUC(TrueScore,FalseScore,2000,0)
        %}
elseif (strcmp(cv_setting, 'S2') )

 Parameters.lamda_ar = 0.05;         Parameters.lamda_ac = 0.05;
 Parameters.lamda_r = 0.5;           Parameters.lamda_c = 0.5;

    OutputScores = zeros(size(Adj));
    if ratio_ >1
    latent_dim = ratio_; % directly use ratio_ as latent dim
    %disp('directly use ratio_ as the latent dim')
    end
    [U_row,S_row,~]=TurnSimBySVD(Sim_row);
    [U_col,S_col,~]=TurnSimBySVD(Sim_col);
    
    ID_Tst_r=left_out;
    b=1:(size(Adj,1)); b(left_out)=[]; ID_Trn_r=b';
    
    TrnLabelMat= Adj(ID_Trn_r,:); % may contain all-zero columns, which cause S4 prediction
    TstLabelMat = Adj(ID_Tst_r,:);%%%%
    %TrnLabelMat= Adj;
    
    % core function
    if ratio_<=1 % using the low-rank
        latent_dim=fix(rank(TrnLabelMat)* ratio_);
    end
    Predicted_LabelMat =PredictS2_Optim(TrnLabelMat,U_row,S_row,ID_Trn_r, ID_Tst_r, U_col,S_col,latent_dim,Parameters);
    %feed in full matrix so that it knows which aretrue test values for it
    %to correctly calculate auc  n aupr in the below function (i.e. coment y2(test_ind) = 0;  line in nfold.m)
    %Measure_S2(Predicted_LabelMat,TrnLabelMat,TstLabelMat)
 
    final=zeros(size(Adj)); final(ID_Tst_r,:)=Predicted_LabelMat;
elseif(strcmp( cv_setting, 'S3') )

        Adj=Y'; Sim_col=(d_s+d_s')/2; Sim_row=(t_s+t_s')/2;
        Parameters.lamda_ar = 0.5;          Parameters.lamda_ac = 0.5;
        Parameters.lamda_r = 0.05;          Parameters.lamda_c = 0.05;
  
    OutputScores = zeros(size(Adj));
    if ratio_ >1
        latent_dim = ratio_; % directly use ratio_ as latent dim
        %disp('directly use ratio_ as the latent dim')
    end
    [U_row,S_row,~]=TurnSimBySVD(Sim_row);
    [U_col,S_col,~]=TurnSimBySVD(Sim_col);
    
     ID_Tst_r=left_out;
     b=1:(size(Adj,1)); b(left_out)=[]; ID_Trn_r=b';
    
  
    TrnLabelMat= Adj(ID_Trn_r,:); % may contain all-zero columns, which cause S4 prediction
    %TstLabelMat = Adj(ID_Tst_r,:);
    %TrnLabelMat= Adj;
    
    if ratio_<=1 % using the low-rank
        latent_dim=fix(rank(TrnLabelMat)* ratio_);
    end
    Predicted_LabelMat =PredictS2_Optim(TrnLabelMat,U_row,S_row,ID_Trn_r, ID_Tst_r, U_col,S_col,latent_dim,Parameters);
    %feed in full matrix so that it knows which aretrue test values for it
    %to correctly calculate auc  n aupr in the below function (i.e. coment y2(test_ind) = 0;  line in nfold.m)
    %Measure_S2(Predicted_LabelMat,TrnLabelMat,TstLabelMat)
     final=zeros(size(Adj)); final(ID_Tst_r,:)=Predicted_LabelMat; final=final';
 
end
        
        Yhat= final;
end

function [AUC_S2,AUPR_S2] = Measure_S2(Predicted_Scores,trn_Adj,tst_Adj)
%% global measure
threshold  =1;
degrees_= sum(trn_Adj,1); % to remove S4 cases, if occured.
%Scores  = Predicted_Scores(:,degrees_>=threshold);
%Label_mat = tst_Adj(:,degrees_>=threshold);
Scores  = Predicted_Scores;
Label_mat = tst_Adj;

TrueScore = Scores( Label_mat(:)>0);
FalseScore= Scores( Label_mat(:)==0);

[AUC_S2, AUPR_S2 ]=EstimationAUC(TrueScore,FalseScore,2000,1);
end

function [U,S,V]=TurnSimBySVD(Sim)
[U,S,V]= svd(Sim);
[U,S,V,~]=CleanSVD(U,S,V);

end



function [trn_row_s1,trn_col_s1] = SplitOutPureS1Scenarios(AdjMat)
Col_flag = (sum(AdjMat,1) > 0);
Row_flag = ( sum(AdjMat,2) > 0);

Col_flag_mat = repmat(Col_flag, size(AdjMat,1),1);
Row_flag_mat = repmat(Row_flag, 1, size(AdjMat,2));

% 0 in Flag is the case of S4, while other entries==1 in row are the
% cases of S2, and other entries==1 in col are the cases of S3. The
% entries ==2 are the cases of S1;
ScenarioFlag = Col_flag_mat + Row_flag_mat;

S1_Flag= 2;

% 1: Pure S1
[r_s1,c_s1]=find(ScenarioFlag==S1_Flag);
trn_row_s1 = unique(r_s1);
trn_col_s1 = unique(c_s1);
end


function Prior = GetPrior(Adj)
DiscreteList = unique(Adj(:));
count_ = zeros(length(DiscreteList),1);
for n=1: length(DiscreteList)
    count_(n) =length( find(Adj(:) == DiscreteList(n) ) ) ;
end
Prior = sum(count_ .* DiscreteList)/ sum(count_);
tol_ = 1e-5;
if Prior< tol_
    Prior = tol_;
end
end
    
    