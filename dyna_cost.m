function [f,c,fx,fu,fxx,fxu,fuu,cx,cu,cxx,cxu,cuu] = dyna_cost(dyn, cst, x, u, i)
% combine dynamics (continuous) and cost functions  for iLQG
% dyn, cst: function handlers of
dt = 0.02;
op.sdt = 0.02;
disdyn = @(x,u)discrete_dynamics(dyn, x, u, op);
if nargout == 2 
    f = disdyn(x, u);
    c = cst(x, u, i);
    if ~isnan(u)
        c = c*dt;
    end
else
    [dimx, N] = size(x);
    dimu = size(u,1);
    cx = zeros(dimx,N);
    cxx= zeros(dimx,dimx,N);
    cu = zeros(dimu,N);
    cuu= zeros(dimu,dimu,N);
    cxu= zeros(dimx,dimu,N);
    fx = zeros(dimx,dimx,N);
    fu = zeros(dimx,dimu,N);
    for n = 1:N-1
        % linearize dynamics
        %[ff, f_x, f_u] = f(x(:,n), u(:,n));
        xu = [x(:,n);u(:,n)];
        ff = @(xu)disdyn(xu(1:dimx),xu(dimx+1:end));
        J = get_jacobian_fd(ff,xu);
        fx(:,:,n) = J(:,1:dimx);
        fu(:,:,n) = J(:,dimx+1:end);
        

        % quadratize cost, adjust for dt
        [~,l_x,l_xx,l_u,l_uu,l_ux] = cst(x(:,n), u(:,n), i);
        %q0(    n) = dt*l0;
        cx (  :,n) = l_x*dt;
        cxx (:,:,n) = l_xx*dt;
        cu (  :,n) = l_u*dt;
        cuu (:,:,n) = l_uu*dt;
        cxu (:,:,n) = l_ux'*dt;
        
    end
    
    [~,cx(:,N),cxx(:,:,N)] = cst(x(:,N), NaN, i);
    
    
    fx(:,:,N) = NaN(dimx,dimx);
    fu(:,:,N) = NaN(dimx,dimu);
    
    
    [fxx,fxu,fuu] = deal([]);

    
    
    [f,c] = deal([]);
end
end

