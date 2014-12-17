function x = NormL1_project_pos(x,weights,tau)

if isreal(x)
   x(x<0)=0; 
   x = oneProjector(x,weights,tau);
else
   xa  = abs(x);
   idx = xa < eps;
   xc  = oneProjector(xa,weights,tau);
   xc  = xc ./ xa; xc(idx) = 0;
   x   = x .* xc;
end
