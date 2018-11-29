% to be completed
function U = wPCA(x,d)
    %compute number of samples as N
    N = length(x(:,1));
    x = x';
    % Data Centering using matrix multiplication
    X = x *(eye(N) - 1/N * ones(N,1) * ones(N,1)');
    % PCA step 1:Compute dot product matrix
    S = X'*X;
    S = 0.5*(S + S');
    % PCA step 2:Perform eigenanalysis of S 
    [v,lambda] = eigs(S,d);
    % PCA step 3:Compute eigenvectors
    % U is the first d dimension U_d = [u1, u2,...,ud]
    U = X * v * lambda^(-0.5);
    % Do whiten step
    U = U * lambda ^-0.5;