function F = EstimateFundamentalMatrix(x1, x2)
%% EstimateFundamentalMatrix
% Estimate the fundamental matrix from two image point correspondences 
% Inputs:
%     x1 - size (N x 2) matrix of points in image 1
%     x2 - size (N x 2) matrix of points in image 2, each row corresponding
%       to x1
% Output:
%    F - size (3 x 3) fundamental matrix with rank 2

N=8;
% Construct 8x9 matrix A
A = [x1(1:N,1).*x2(1:N,1) x1(1:N,1).*x2(1:N,2) x1(1:N,1) x1(1:N,2).*x2(1:N,1) ...
    x1(1:N,2).*x2(1:N,2) x1(1:N,2) x2(1:N,1) x2(1:N,2) ones(N,1)];
% Solving via SVD
[~,~,Vt] = svd(A);
V = Vt';
Fr = V(:,9);
F = reshape(Fr,3,3);
% Applying rank constraint for rank(F)=2
[U,D,Vt] = svd(F);
D_tilde = D;
D_tilde(3,3) = 0;
F = U*D_tilde*Vt;
% Normalize the fundamental matrix
F = F/norm(F);


