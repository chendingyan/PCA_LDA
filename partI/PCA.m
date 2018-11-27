% to be completed
function U = PCA(x,d)

  x=x';
  %compute number of samples as N
  N = size(x)(2)
  %first we need to normalize 
  X = x*(eye(N) - 1/N * ones(N,1) * ones(N,1)')
  % PCA step 1:Compute dot product matrix
  S = X'*X;
  % PCA step 2:Perform eigenanalysis of S 
  [v,l] = eig(S);
  % PCA step 3:Compute eigenvectors
  u = X * v * l^-0.5;
  % U is the first d dimension U_d = [u1, u2,...,ud]
  U = v(:,1:d);