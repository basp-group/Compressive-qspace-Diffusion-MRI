function [ y ] = Phi( x,d,v,A )

% This operator applies matrix A(m,d) to matrix X(d,v) A*X and then unfolds
% the result in a vector
X=reshape(x,d,v);
Y=A*X;

y=Y(:);

end

