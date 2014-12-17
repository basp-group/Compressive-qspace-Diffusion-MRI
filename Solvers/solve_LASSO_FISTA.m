function [sol, info] = solve_LASSO_FISTA(y, param) 


% [sol, info] = solve_LASSO_FISTA(y, param)
% solve nonconvex, linear measurements' Fordward-Backward

% Optional input arguments
if ~isfield(param, 'verbose'), param.verbose = 1; end
if ~isfield(param, 'rel_obj'), param.rel_obj = 1e-3; end
if ~isfield(param, 'max_iter'), param.max_iter = 200; end
if ~isfield(param, 'gamma'), param.gamma = 1e-2; end


sol=param.solini; u=sol; v=u; told=1;
iter = 0; prev_sol = param.solini;


y=y(:);

% Main loop
while 1
    
    %
    if param.verbose >= 1
        fprintf('Iteration %i:\n', iter);
    end
    
    % gradient descent
    g=grad_L2norm(u,y,param.A,param.At);
    dummy=u-param.gamma*g; e_norm=0.5*norm(param.A(dummy)-y).^2;
    
    % proximal operator
    sol=NormL1_project_pos(dummy, param.weights, param.K);
    
    
    
    %verify stopping criteria
    rel_var = norm(sol - prev_sol)/norm(prev_sol);
    Vrnorm(iter+1)=rel_var;
    
    
    
    if param.verbose >= 1
        fprintf('  0.5*(||Ax-b||_2)^2 = %e, rel_var = %e\n', ...
            e_norm, rel_var);
    end
    
    if (rel_var < param.rel_obj || norm(sol-prev_sol)==0)
        info.crit_ncFB = 1;% 1 == 'REL_NORM';
        break;
    elseif iter >= param.max_iter
        info.crit_ncFB = 0;% 0 == 'MAX_IT';
        break;
    end
    
    % update FISTA
    t = (1+sqrt(1+4*told^2))/2;
    u = v; 

    g=grad_L2norm(sol,y,param.A,param.At);
    dummy=sol-param.gamma*g;
    v=NormL1_project_pos(dummy, param.weights, param.K);


    u = v + (told-1)/t * (v - u);


    % Update number of iteration
    told = t;
    
    % Update variables
    iter = iter + 1;
    prev_sol = sol;
     
end

info.niter=iter;

% Log
if param.verbose>=1
    
    fprintf('\n Solution found:\n');
    
    % Residual
    dummy = param.A(sol); res = norm(y(:)-dummy(:), 2);
    fprintf(' ||y-Ax||_2=%e\n', res);
    
    % Stopping criterion
    fprintf(' %i iterations\n', iter);
    fprintf(' Stopping criterion: %s \n\n', info.crit_ncFB);
    
end

end