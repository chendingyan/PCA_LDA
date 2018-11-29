% to be completed
function W = LDA(x,y)
    count = hist(y,unique(y));
    %compute the M matrix, the result should be a 340*340 matrix, with E1,
    %E2...E_class_num in it. Each E should be 5*5 matrix of 1/5
    M = zeros(sum(count));
    class_num = length(count);
    for i = 1:class_num
        index = find(y == i);
        len = length(index);
        temp = ones(len)/len;
        M(index, index) = temp;
    end
    N = length(x(:,1));
    x = x';
    % Data Centering using matrix multiplication
    X = x *(eye(N) - 1/N * ones(N,1) * ones(N,1)');
    % PCA step 1:Compute dot product matrix
    S = X'*X;
    % compute the within-class scatter matrix and perform eigen-analysis 
    I = eye(length(M));
    Sw = (I-M)*S*(I-M);
    Sw = 0.5*(Sw+transpose(Sw));
    [V_Sw,lambda_Sw] = eigs(Sw,N-length(count));
    % data whitening
    U = X*(I-M)*V_Sw*inv(lambda_Sw);
    % compute the between-class scatter matrix 
    Sb=transpose(U)*X*M*transpose(X)*U;
    Sb = 0.5*(Sb+transpose(Sb));
    % compute largest C-1 eigenvectors and eigenvalues of Sb
    [evec_bd,eval_bd] = eigs(Sb,length(count)-1);
    % total transformation
    W=U*evec_bd;

    