function [ Aty ] = Phit( y,m,v,At )

% This operator applies matrix At(d,m) to matrix Y(m,v) At*Y and then unfolds
% the result in a vector

% m: size of measurement vector // d: number of directions // v: number of
% voxels

Y=reshape(y,m,v);

Aty=At*Y;

Aty=Aty(:);



end

