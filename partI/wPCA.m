% to be completed
function U = PCA(x,d)

%compute number of samples as N
  size_array = size(x);
  N = size_array(1);
  % Data Centering using matrix multiplication
  X = x'*(eye(N) - 1/N * ones(N,1) * ones(N,1)');
  % PCA step 1:Compute dot product matrix
  S = X'*X;
  % PCA step 2:Perform eigenanalysis of S 
  [v,l] = eig(S);
  % PCA step 3:Compute eigenvectors
  u = X * v * l^-0.5;
  % U is the first d dimension U_d = [u1, u2,...,ud]
  U = v(:,1:d);
  U = u_d * l^-0.5;