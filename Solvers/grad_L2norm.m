
function [ g ] = grad_L2norm( x, y, A, At )
%
% Compute the gradient of the L2norm of the data term ||A(X) - Y ||_fro
%
%   INPUT:
%   - x: input vector (or matrix)
%   - A/At: operator and its transpose
%   - y: data vector (or matrix)



x=x(:);
N=size(x,1);

g=A(x)-y;
g=At(g);
g=g(:);



end



