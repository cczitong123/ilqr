classdef ILQRController < OptController
    %ILQRCONTROLLER Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        
        syml
        symlf
        symlx
        symlxx
        symlu
        symluu
        symlux
        symlfx
        symlfxx
        
        lders
        lfders
        
        f
        j
    end
    
    methods
        function oc = ILQRController(robot_model, task)
            
            
            %call superclass constructor
            oc@OptController(robot_model, task);
            
            oc.opt_param = [];
            oc.opt_param.umax = oc.umax;
            oc.opt_param.umin = oc.umin;
            %po.lambda_init = 0.01;
            oc.opt_param.lambda_max  = 0.05;
            oc.opt_param.iter_max = 100;
            oc.opt_param.online_plotting = 0;
            oc.opt_param.online_printing = 1;
            oc.opt_param.dcost_converge = 10^-5;
            oc.opt_param.sim_solver = 'rk4';
            % global_optimal by multiple inits
            if( ~isfield(task.task_param, 'multiple_inits'))
                oc.task.task_param.multiple_inits = 1;
            end
            
            % register dynamics, the function should be able to output
            % jacobians
            oc.f = @(x, u) robot_model.dynamics_with_jacobian( x, u );
            %fprintf(1,'robot dynamics loaded.')
            
            % setup task specific cost function
            
            
            %task.cost_param.epsilon = 0 ; %10^-8 ;
            %task.cost_param.dt = oc.cdt;
            %task.cost_param.fd = 0;
            %task.cost_param.w = oc.task_param.w;
            %task.cost_param.model = oc.robot_model;
            %task.cost_param.rege_ratio = oc.task_param.rege_ratio;
            
            %oc.j = @(x,u,t) j_reach_mw_mcc1dof_md( x, u, t, oc.cost_param );
            
            %fprintf(1,'cost function setup done.')
            
%             xv=sym('x',[8 1]);
%             uv=sym('u',[3 1]);
%             oc.syml = symfun(task.cost_run(xv,uv),[xv;uv]);
%             oc.symlx = jacobian(oc.syml, xv);
%             oc.symlxx = hessian(oc.syml, xv);
%             oc.symlu = jacobian(oc.syml,uv);
%             oc.symluu = hessian(oc.syml,uv);
%             oc.symlux = jacobian(gradient(oc.syml,uv),xv);
%             oc.lders.lxtmp = matlabFunction(oc.symlx);
%             oc.lders.lxxtmp = matlabFunction(oc.symlxx);
%             oc.lders.lutmp = matlabFunction(oc.symlu);
%             oc.lders.luutmp = matlabFunction(oc.symluu);
%             oc.lders.luxtmp = matlabFunction(oc.symlux);
%             oc.lders.lx =@(x,u)oc.lders.lxtmp (x(1),x(2),x(3),x(4),x(5),x(6),x(7),x(8),u(1),u(2),u(3));
%             oc.lders.lxx =@(x,u)oc.lders.lxxtmp (x(1),x(2),x(3),x(4),x(5),x(6),x(7),x(8),u(1),u(2),u(3));
%             oc.lders.lu =@(x,u)oc.lders.lutmp (x(1),x(2),x(3),x(4),x(5),x(6),x(7),x(8),u(1),u(2),u(3));
%             oc.lders.luu =@(x,u)oc.lders.luutmp (x(1),x(2),x(3),x(4),x(5),x(6),x(7),x(8),u(1),u(2),u(3));
%             oc.lders.lux =@(x,u)oc.lders.luxtmp (x(1),x(2),x(3),x(4),x(5),x(6),x(7),x(8),u(1),u(2),u(3));
%             
%             
%             oc.symlf = symfun(task.cost_final(xv),xv);
%             oc.symlfx = jacobian(oc.symlf,xv);
%             oc.symlfxx= hessian(oc.symlf,xv);
%             oc.lfders.lfxtmp = matlabFunction(oc.symlfx);
%             oc.lfders.lfxxtmp = matlabFunction(oc.symlfxx);
%             oc.lfders.lx = @(x,u) oc.lfders.lfxtmp(x(1),x(2),x(3),x(4),x(5),x(6),x(7),x(8));
%             oc.lfders.lxx= @(x,u) oc.lfders.lfxxtmp(x(1),x(2),x(3),x(4),x(5),x(6),x(7),x(8));
        end
        
        function [l, l_x, l_xx, l_u, l_uu, l_ux] = l_run(oc,x,u)
            l = oc.task.cost_run(x,u);
            l_x = oc.lders.lx(x,u);
            l_xx=oc.lders.lxx(x,u);
            l_u = oc.lders.lu(x,u);
            l_uu = oc.lders.luu(x,u);
            l_ux = oc.lders.lux(x,u);
        end
        
        function [l, l_x, l_xx ] = l_final(oc,x)
            l = oc.task.cost_final(x);
            l_x = oc.lfders.lx(x);
            l_xx= oc.lfders.lxx(x);
        end
        
        
        % todo:
        function [x, u, L, cost, iter] = ilqr_analytic(oc, dt, N, x0, u0, p)
            % todo: finish this ilqr using analytical derivatives of f,j
            % unpack parameter struct
            if isfield(p,'lambda_init'        ),lambda_init         = p.lambda_init;
            else
                lambda_init        = 1e-2    ;
            end
            if isfield(p,'lambda_factor'      ),lambda_factor       = p.lambda_factor      ;
            else
                lambda_factor      = sqrt(10);
            end
            if isfield(p,'lambda_max'         ),lambda_max          = p.lambda_max         ;
            else
                lambda_max         = 1e-2    ;
            end
            if isfield(p,'dcost_converge'     ),dcost_converge      = p.dcost_converge  ;
            else
                dcost_converge     = 1e-9    ;
            end
            if isfield(p,'iter_max'           ),iter_max            = p.iter_max     ;
            else
                iter_max           = 100     ;
            end
            if isfield(p,'online_printing'    ),online_printing     = p.online_printing;
            else
                online_printing    = 1       ;
            end
            if isfield(p,'online_plotting'    ),online_plotting     = p.online_plotting;
            else
                online_plotting    = 0       ;
            end
            if isfield(p,'umax'               ),umax                = p.umax ;
            else
                umax               = inf     ;
            end
            if isfield(p,'umin'               ),umin                = p.umin     ;
            else
                umin               =-inf     ;
            end
            if isfield(p,'sim_solver'             ),sim_solver              = p.sim_solver  ;
            else
                sim_solver             ='rk4'  ;
            end
            %f = oc.f;
            % initialise everything
            dimX = size(x0, 1);   % state dimensionality
            dimU = size(u0, 1);   % command dimensionality
            ps=oc.sim_param;ps.solver=sim_solver; % simulation parameters
            
            if size(u0,2)==1                                   % initial command sequence
                u = repmat(u0,1,N-1);                               % ...if one initial command given, use as full initial command sequence
            else                                               %
                u    = u0;
            end                                       % ...if full command sequence given, use as initial command sequence
            L    = zeros(dimU, dimX, N-1);                      % initial feedback gain sequence
            x    = simulate_feedforward ( x0, oc.f, u, ps );       % initial state sequence
            cost = oc.cost_total ( x, u ); % initial cost
            t    = (0:N-1)*dt;                                  % time, N timestemps
            
            % initialise other matrices (for speed)
            A  = zeros(dimX,dimX,N-1);
            B  = zeros(dimX,dimU,N-1);
            I  =   eye(dimX         );
            q0 = zeros(        1,N-1);
            q  = zeros(     dimX,N-1);
            Q  = zeros(dimX,dimX,N-1);
            r  = zeros(     dimU,N-1);
            R  = zeros(dimU,dimU,N-1);
            P  = zeros(dimU,dimX,N-1);
            s0 = zeros(        1,N  );
            s  = zeros(     dimX,N  );
            S  = zeros(dimX,dimX,N  );
            l  = zeros(     dimU,N-1);
            
            % initialise lambda
            lambda = lambda_init;
            % initialise update flag
            update = 1;
            
            lambdas = NaN(iter_max);
            costs = NaN(iter_max);
            
            % set up plotting
            if online_plotting
                online_plotting_fig=figure;clf,set(online_plotting_fig,'Name','iLQR: Cost and Lambda Convergence'),set(online_plotting_fig,'NumberTitle','off')
                subplot(1,2,1),hold on,grid on,set(gca,'Yscale','linear'),title('cost'  ),ylabel('cost'  ),xlabel('iteration')
                subplot(1,2,2),hold on,grid on,set(gca,'Yscale','log')   ,title('lambda'),ylabel('lambda'),xlabel('iteration')
            end
            
            % main loop
            for iter = 1:iter_max
                %------ STEP 1: approximate dynamics and cost along new trajectory'
                if update
                    update = 0;
                    % compute LQ approximation
                    for n = 1:N-1
                        % linearize dynamics, adjust for dt
                        %[ff, f_x, f_u] = f(x(:,n), u(:,n));
                        [~, f_x, f_u] = oc.f(x(:,n), u(:,n));
                        A(:,:,n) = I + dt*f_x;
                        B(:,:,n) =     dt*f_u;
                        
                        % quadratize cost, adjust for dt
                        [l0,l_x,l_xx,l_u,l_uu,l_ux] = oc.l_run(x(:,n), u(:,n));
                        q0(    n) = dt*l0;
                        q (  :,n) = dt*l_x;
                        Q (:,:,n) = dt*l_xx;
                        r (  :,n) = dt*l_u;
                        R (:,:,n) = dt*l_uu;
                        P (:,:,n) = dt*l_ux;
                    end
                    
                    % initialise value function approximation with final cost
                    [s0(N),s(:,N),S(:,:,N)] = oc.l_final(x(:,N));
                end
                
                %------ STEP 2: compute optimal control law and cost-to-go
                diverged=0;
                for n = N-1:-1:1
                    % compute shortcuts g,G,H
                    g = r(  :,n) + B(:,:,n)'*s(  :,n+1);
                    G = P(:,:,n) + B(:,:,n)'*S(:,:,n+1)*A(:,:,n);
                    H = R(:,:,n) + B(:,:,n)'*S(:,:,n+1)*B(:,:,n);
                    
                    % check for divergence
                    if any(~isfinite(H)), diverged=1; break; end % divergence: EXIT
                    
                    % find control law
                    [l(:,n), L(:,:,n)] = u_optimal(g,G,H,u(:,n),umin,umax,lambda);
                    
                    % update cost-to-go approximation
                    S (:,:,n) = Q (:,:,n) + A(:,:,n)'*S (:,:,n+1)*A(:,:,n) +    L(:,:,n)'*H*L(:,:,n) + L(:,:,n)'*G + G'*L(:,:,n);
                    s (  :,n) = q (  :,n) + A(:,:,n)'*s (  :,n+1)          +    L(:,:,n)'*H*l(  :,n) + L(:,:,n)'*g + G'*l(:,n);
                    s0(    n) = q0(    n) +           s0(    n+1)          + .5*l(  :,n)'*H*l(  :,n) + l(  :,n)'*g;
                    
                end
                if diverged, fprintf('Optimal control law and cost-to-go computation diverged. '); break; end% divergence: EXIT
                
                %------ STEP 3A: new control sequence, trajectory, cost
                % simulate linearized system to compute new control
                dx = zeros(dimX,1);
                unew = zeros(dimU,N-1);
                for n=1:N-1
                    du = l(:,n) + L(:,:,n)*dx;
                    du = min(max(du+u(:,n),umin),umax) - u(:,n);
                    % todo: constriant on x
                    dx = A(:,:,n)*dx + B(:,:,n)*du ;
                    unew(:,n) = u(:,n) + du;
                end
                
                %------ STEP 3B: simulate system to compute new trajectory and cost
                xnew    = simulate_feedforward ( x0, oc.f, unew, ps );
                costnew = oc.cost_total ( xnew, unew );% initial cost
                
                %------ STEP 4: Levenberg-Marquardt method
                if costnew<cost
                    % decrease lambda (get closer to Newton method)
                    lambda = lambda / lambda_factor;
                    dcost  = cost-costnew; % decrease in cost
                    % update x, u, cost
                    u    = unew;
                    x    = xnew;
                    cost = costnew;
                    % check for convergence
                    if iter>1 && dcost<dcost_converge
                        if online_printing==1
                            fprintf('Cost improvement threshold reached. ');
                        end
                        break; % improvement too small: EXIT
                    end
                    % flag update
                    update = 1;
                else
                    % increase lambda (get closer to gradient descent)
                    lambda = lambda * lambda_factor;
                    if lambda>lambda_max
                        fprintf('Max. lambda reached. '); break; % lambda too large: EXIT
                    end
                end
                
                % plot/print stuff out
                if online_printing==1,fprintf('Iteration = %d; Cost = %.4f; log(Lambda) = %.1f\n',iter,cost,log10(lambda));end
                if online_plotting==1
                    lambdas(iter)=lambda; costs(iter)=cost;
                    set(0,'CurrentFigure',online_plotting_fig),
                    subplot(1,2,1),plot(0:iter-1,costs);
                    subplot(1,2,2),plot(0:iter-1,lambdas);
                end
            end
            
            if iter==iter_max
                fprintf('Max. number of iterations reached. ');
            end
            % print final result if necessary
            if online_printing==2
                fprintf('Iterations = %d;  Cost = %.4f', iter, cost);
            end
            
            fprintf('\n');
            % close figure
            if online_plotting==1, close(online_plotting_fig),end
        end
        
        function cost = cost_total(oc,x,u)
            
            cost_final = oc.l_final(x(:,end));
            cost_run = zeros(1,oc.Nt-1);
            for i=1:oc.Nt - 1
                cost_run(i) = oc.l_run(x(:,i),u(:,i))*oc.cdt ;
            end
            cost = sum(cost_run) + cost_final;
            
        end
        
    end
    
    
    methods (Static)
        % Function for solving (finite horizon) optimal control problems with the Iterative Linear Quadratic Regulator method.
        %
        %    [x, u, L, cost] = ilqr(f, j, dt, N, x0, u0, p )
        %
        % in:
        %     f  - dynamics function (function handle)
        %     j  - cost function (function handle)
        %     dt - time step
        %     N  - number of steps
        %     x0 - start state
        %     u0 - initial command sequence
        %     p  - parameter struct (optionally) containing:
        %      .lambda_init     - initial value of Levenberg-Marquardt lambda
        %      .lambda_factor   - factor for multiplying or dividing lambda
        %      .lambda_max      - threshold on lambda (exit if lambda exceeds this value)
        %      .iter_max        - threshold on number of iterations (exit if exceeded)
        %      .dcost_converge  - threshold on relative improvement in cost (exit if improvement less than this value)
        %      .online_printing - print: {0:never, 1:every iter, 2:final}
        %      .online_plotting - plot:  {0:never, 1:every iter, 2:final}
        %      .umax            - maximum command
        %      .umin            - minimum command
        %      .solver          - simulation solver
        %      (any fields that are unspecified will be set by default values).
        %
        % out:
        %     x      - predicted optimal state sequence
        %     u      - optimal feed-forward command sequence
        %     L      - optimal feedback gain sequence
        %     cost   - predicted cost
        %
        %> \author Matthew Howard (MH) matthew.howard@ed.ac.uk
        %> \date 19/06/11 19:14:08
        %
        function [result] = run_multiple(f, j, dt, N, x0, u0, p)
            %tic
                % set ilqr parameters
                
                Ninit = 8; % number of initial commands
                u0s = zeros(3,Ninit);
                tempu = rand(3,Ninit);
                %u0s(1,:) = umin(1)+(umax(1)-umin(1))*tempu(1,:);
                u00 = u0;
                u01 = [0;0.1;0];
                u0s(1,:) = p.umin(1)+(p.umax(1)-p.umin(1))*tempu(1,:);
                u0s(2,:) = p.umin(2)+(p.umax(2)-p.umin(2))*tempu(2,:);
                u0s(3,:) = p.umin(3)+(p.umax(3)-p.umin(3))*tempu(3,:);
                
                u0s = [u00,u01,u0s];
                Ninit=Ninit + 2; % total inits is 5
                
                xs = cell(Ninit,1);
                us = cell(Ninit,1);
                Ls = cell(Ninit,1);
                costs = zeros(Ninit,1);
                %costs_forward = zeros(Ninit,1);
                
                for i = 1:Ninit
                    u0 = u0s(:,i) ;
                    result_tmp = ILQRController.ilqr(f,j,dt,N,x0,u0,p);
                    
                    % run controller on plant
                    %ppi = []; ppi.xn = xx; ppi.un = uu; ppi.Ln = LL;
                    %Pi = @(x,n)pi_ilqr(x,n,ppi);
                    Ls{i} = result_tmp.L;
                    xs{i}=result_tmp.x;
                    us{i}=result_tmp.u;
                    costs(i)=result_tmp.cost;
                    %[xs{i},us{i}] = simulate_feedback_time_indexed ( x0, f, Pi, ps );
                    % evaluate cost of trajectory on plant
                    %costs(i) = evaluate_trajectory_cost_fh(xs{i},us{i},j,ps);
                    %costs_forward(i) = cost_predicted;
                end
                [cost,index_cost] = min(costs) ;
                %[~,index_cost_forward] = min(costs_forward) ;
                
                %if index_cost ~= index_cost_forward
                %    fprintf('optimal feedback trajectory is not refered to optimal feedforward solution');
                %end
                result.x = xs{index_cost};
                result.u = us{index_cost};
                result.L = Ls{index_cost};
                result.u0 = u0s(:,index_cost); %return the u0 that generate the best result
                result.Ninit = Ninit;
                result.u0s = u0s;
                result.cost = cost;
                result.costs = costs;
                %toc
                fprintf(1,'Cost (evaluated on plant) = %f\n',cost);
        end
        
        % result:
        %   x
        %   u
        %   L
        %   cost
        %   iter
        %   lambdas
        %   costs
        function [result] = ilqr(f, j, dt, N, x0, u0, p )
            
            % unpack parameter struct
            if isfield(p,'lambda_init'        ),lambda_init         = p.lambda_init;
            else
                lambda_init        = 1e-2    ;
            end
            if isfield(p,'lambda_factor'      ),lambda_factor       = p.lambda_factor      ;
            else
                lambda_factor      = 1.6;
            end
            if isfield(p,'lambda_max'         ),lambda_max          = p.lambda_max         ;
            else
                lambda_max         = 1e10    ;
            end
            if isfield(p,'dcost_converge'     ),dcost_converge      = p.dcost_converge  ;
            else
                dcost_converge     = 1e-4    ;
            end
            if isfield(p,'iter_max'           ),iter_max            = p.iter_max     ;
            else
                iter_max           = 100     ;
            end
            if isfield(p,'online_printing'    ),online_printing     = p.online_printing;
            else
                online_printing    = 1       ;
            end
            if isfield(p,'online_plotting'    ),online_plotting     = p.online_plotting;
            else
                online_plotting    = 1       ;
            end
            if isfield(p,'umax'               ),umax                = p.umax ;
            else
                umax               = inf     ;
            end
            if isfield(p,'umin'               ),umin                = p.umin     ;
            else
                umin               =-inf     ;
            end
            if isfield(p,'solver'             ),solver              = p.solver  ;
            else
                solver             ='rk4'  ;
            end
            
            % initialise everything
            dimX = size(x0, 1);   % state dimensionality
            dimU = size(u0, 1);   % command dimensionality
            ps=[];ps.dt=dt;ps.solver=solver; % simulation parameters
            ps.T = p.T;
            if size(u0,2)==1                                   % initial command sequence
                u = repmat(u0,1,N-1);                               % ...if one initial command given, use as full initial command sequence
            else                                               %
                u    = u0;
            end                                       % ...if full command sequence given, use as initial command sequence
            L    = zeros(dimU, dimX, N-1);                      % initial feedback gain sequence
            x    = simulate_feedforward ( x0, f, u, ps );       % initial state sequence
            cost = evaluate_trajectory_cost_fh ( x, u, j, ps ); % initial cost
            t    = (0:N-1)*dt;                                  % time
            
            % initialise other matrices (for speed)
            A  = zeros(dimX,dimX,N-1);
            B  = zeros(dimX,dimU,N-1);
            I  =   eye(dimX         );
            q0 = zeros(        1,N-1);
            q  = zeros(     dimX,N-1);
            Q  = zeros(dimX,dimX,N-1);
            r  = zeros(     dimU,N-1);
            R  = zeros(dimU,dimU,N-1);
            P  = zeros(dimU,dimX,N-1);
            s0 = zeros(        1,N  );
            s  = zeros(     dimX,N  );
            S  = zeros(dimX,dimX,N  );
            l  = zeros(     dimU,N-1);
            
            % initialise lambda
            lambda = lambda_init;
            % initialise update flag
            update = 1;
            
            lambdas = NaN(iter_max,1);
            costs = NaN(iter_max,1);
            
            % set up plotting
            if online_plotting
                online_plotting_fig=figure;clf,set(online_plotting_fig,'Name','iLQR: Cost and Lambda Convergence'),set(online_plotting_fig,'NumberTitle','off')
                subplot(1,2,1),hold on,grid on,set(gca,'Yscale','linear'),title('cost'  ),ylabel('cost'  ),xlabel('iteration')
                subplot(1,2,2),hold on,grid on,set(gca,'Yscale','log')   ,title('lambda'),ylabel('lambda'),xlabel('iteration')
            end
            
            % main loop
            for iter = 1:iter_max
                %------ STEP 1: approximate dynamics and cost along new trajectory'
                if update
                    update = 0;
                    % compute LQ approximation
                    for n = 1:N-1
                        % linearize dynamics, adjust for dt
                        %[ff, f_x, f_u] = f(x(:,n), u(:,n));
                        [~, f_x, f_u] = f(x(:,n), u(:,n));
                        A(:,:,n) = I + dt*f_x;
                        B(:,:,n) =     dt*f_u;
                        
                        % quadratize cost, adjust for dt
                        [l0,l_x,l_xx,l_u,l_uu,l_ux] = j(x(:,n), u(:,n), t(n));
                        q0(    n) = dt*l0;
                        q (  :,n) = dt*l_x;
                        Q (:,:,n) = dt*l_xx;
                        r (  :,n) = dt*l_u;
                        R (:,:,n) = dt*l_uu;
                        P (:,:,n) = dt*l_ux;
                    end
                    
                    % initialise value function approximation with final cost
                    [s0(N),s(:,N),S(:,:,N)] = j(x(:,N), NaN, t(N));
                end
                
                %------ STEP 2: compute optimal control law and cost-to-go
                diverged=0;
                for n = N-1:-1:1
                    % compute shortcuts g,G,H
                    g = r(  :,n) + B(:,:,n)'*s(  :,n+1);
                    G = P(:,:,n) + B(:,:,n)'*S(:,:,n+1)*A(:,:,n);
                    H = R(:,:,n) + B(:,:,n)'*S(:,:,n+1)*B(:,:,n);
                    
                    % check for divergence
                    if any(~isfinite(H)), diverged=1; break; end % divergence: EXIT
                    
                    % find control law
                    [l(:,n), L(:,:,n)] = u_optimal(g,G,H,u(:,n),umin,umax,lambda);
                    
                    % update cost-to-go approximation
                    S (:,:,n) = Q (:,:,n) + A(:,:,n)'*S (:,:,n+1)*A(:,:,n) +    L(:,:,n)'*H*L(:,:,n) + L(:,:,n)'*G + G'*L(:,:,n);
                    s (  :,n) = q (  :,n) + A(:,:,n)'*s (  :,n+1)          +    L(:,:,n)'*H*l(  :,n) + L(:,:,n)'*g + G'*l(:,n);
                    s0(    n) = q0(    n) +           s0(    n+1)          + .5*l(  :,n)'*H*l(  :,n) + l(  :,n)'*g;
                    
                end
                if diverged, fprintf('Optimal control law and cost-to-go computation diverged. '); break; end% divergence: EXIT
                
                %------ STEP 3A: new control sequence, trajectory, cost
                % simulate linearized system to compute new control
                dx = zeros(dimX,1);
                unew = zeros(dimU,N-1);
                for n=1:N-1
                    du = l(:,n) + L(:,:,n)*dx;
                    du = min(max(du+u(:,n),umin),umax) - u(:,n);
                    dx = A(:,:,n)*dx + B(:,:,n)*du ;
                    unew(:,n) = u(:,n) + du;
                end
                
                %------ STEP 3B: simulate system to compute new trajectory and cost
                xnew    = simulate_feedforward ( x0, f, unew, ps );
                costnew = evaluate_trajectory_cost_fh ( xnew, unew, j, ps );% initial cost
                
                %------ STEP 4: Levenberg-Marquardt method
                if costnew<cost
                    % decrease lambda (get closer to Newton method)
                    lambda = lambda / lambda_factor;
                    dcost  = cost-costnew; % decrease in cost
                    % update x, u, cost
                    u    = unew;
                    x    = xnew;
                    cost = costnew;
                    % check for convergence
                    if iter>1 && dcost<dcost_converge
                        if online_printing==1
                            fprintf('Cost improvement threshold reached. ');
                        end
                        break; % improvement too small: EXIT
                    end
                    % flag update
                    update = 1;
                else
                    % increase lambda (get closer to gradient descent)
                    lambda = lambda * lambda_factor;
                    if lambda>lambda_max
                        fprintf('Max. lambda reached. '); break; % lambda too large: EXIT
                    end
                end
                lambdas(iter)=lambda; costs(iter)=cost;
                % plot/print stuff out
                if online_printing==1,fprintf('Iteration = %d; Cost = %.4f; log(Lambda) = %.1f\n',iter,cost,log10(lambda));end
                if online_plotting==1
                    
                    set(0,'CurrentFigure',online_plotting_fig),
                    subplot(1,2,1),plot(0:iter-1,costs);
                    subplot(1,2,2),plot(0:iter-1,lambdas);
                end
            end
            
            if iter==iter_max
                fprintf('Max. number of iterations reached. ');
            end
            % print final result if necessary
            if online_printing==2
                fprintf('Iterations = %d;  Cost = %.4f', iter, cost);
            end
            
            fprintf('\n');
            % close figure
            if online_plotting==1, close(online_plotting_fig),end
            
            result.x = x;
            result.u = u;
            result.L = L;
            result.cost = cost;
            result.iteration = iter;
            result.lambdas = lambdas;
            result.costs = costs;
        end
        
        function [result] = ilqr_sim(f, j, dt, N, x0, u0, p )
            
            % unpack parameter struct
            if isfield(p,'lambda_init'        ),lambda_init         = p.lambda_init;
            else
                lambda_init        = 1e-2    ;
            end
            if isfield(p,'lambda_factor'      ),lambda_factor       = p.lambda_factor      ;
            else
                lambda_factor      = sqrt(10);
            end
            if isfield(p,'lambda_max'         ),lambda_max          = p.lambda_max         ;
            else
                lambda_max         = 1e-2    ;
            end
            if isfield(p,'dcost_converge'     ),dcost_converge      = p.dcost_converge  ;
            else
                dcost_converge     = 1e-9    ;
            end
            if isfield(p,'iter_max'           ),iter_max            = p.iter_max     ;
            else
                iter_max           = 100     ;
            end
            if isfield(p,'online_printing'    ),online_printing     = p.online_printing;
            else
                online_printing    = 1       ;
            end
            if isfield(p,'online_plotting'    ),online_plotting     = p.online_plotting;
            else
                online_plotting    = 1       ;
            end
            if isfield(p,'umax'               ),umax                = p.umax ;
            else
                umax               = inf     ;
            end
            if isfield(p,'umin'               ),umin                = p.umin     ;
            else
                umin               =-inf     ;
            end
            if isfield(p,'solver'             ),solver              = p.solver  ;
            else
                solver             ='rk4'  ;
            end
            
            % initialise everything
            dimX = size(x0, 1);   % state dimensionality
            dimU = size(u0, 1);   % command dimensionality
            ps=[]; ps.dt = 0.001; ps.solver=solver; % simulation parameters
            ps.T = p.T;
            ps0.dt = dt;ps0.solver = solver;ps0.T = p.T;
            tsim = 0:0.001:p.T;
            t    = (0:N-1)*dt;                                  % time
            if size(u0,2)==1                                   % initial command sequence
                u = repmat(u0,1,N-1);                               % ...if one initial command given, use as full initial command sequence
            else                                               %
                u    = u0;
            end                                       % ...if full command sequence given, use as initial command sequence
            L    = zeros(dimU, dimX, N-1);                      % initial feedback gain sequence
            usim = scale_controlSeq(u,t,tsim);
            x = simulate_feedforward ( x0, f, u, ps0 );
            xsim = simulate_feedforward ( x0, f, usim, ps );       % initial state sequence
            cost = evaluate_trajectory_cost_fh ( xsim, usim, j, ps ); % initial cost
            
            
            % initialise other matrices (for speed)
            A  = zeros(dimX,dimX,N-1);
            B  = zeros(dimX,dimU,N-1);
            I  =   eye(dimX         );
            q0 = zeros(        1,N-1);
            q  = zeros(     dimX,N-1);
            Q  = zeros(dimX,dimX,N-1);
            r  = zeros(     dimU,N-1);
            R  = zeros(dimU,dimU,N-1);
            P  = zeros(dimU,dimX,N-1);
            s0 = zeros(        1,N  );
            s  = zeros(     dimX,N  );
            S  = zeros(dimX,dimX,N  );
            l  = zeros(     dimU,N-1);
            
            % initialise lambda
            lambda = lambda_init;
            % initialise update flag
            update = 1;
            
            lambdas = NaN(iter_max,1);
            costs = NaN(iter_max,1);
            
            % set up plotting
            if online_plotting
                online_plotting_fig=figure;clf,set(online_plotting_fig,'Name','iLQR: Cost and Lambda Convergence'),set(online_plotting_fig,'NumberTitle','off')
                subplot(1,2,1),hold on,grid on,set(gca,'Yscale','linear'),title('cost'  ),ylabel('cost'  ),xlabel('iteration')
                subplot(1,2,2),hold on,grid on,set(gca,'Yscale','log')   ,title('lambda'),ylabel('lambda'),xlabel('iteration')
            end
            
            % main loop
            for iter = 1:iter_max
                %------ STEP 1: approximate dynamics and cost along new trajectory'
                if update
                    update = 0;
                    % compute LQ approximation
                    for n = 1:N-1
                        % linearize dynamics, adjust for dt
                        %[ff, f_x, f_u] = f(x(:,n), u(:,n));
                        [~, f_x, f_u] = f(x(:,n), u(:,n));
                        A(:,:,n) = I + dt*f_x;
                        B(:,:,n) =     dt*f_u;
                        
                        % quadratize cost, adjust for dt
                        [l0,l_x,l_xx,l_u,l_uu,l_ux] = j(x(:,n), u(:,n), t(n));
                        q0(    n) = dt*l0;
                        q (  :,n) = dt*l_x;
                        Q (:,:,n) = dt*l_xx;
                        r (  :,n) = dt*l_u;
                        R (:,:,n) = dt*l_uu;
                        P (:,:,n) = dt*l_ux;
                    end
                    
                    % initialise value function approximation with final cost
                    [s0(N),s(:,N),S(:,:,N)] = j(x(:,N), NaN, t(N));
                end
                
                %------ STEP 2: compute optimal control law and cost-to-go
                diverged=0;
                for n = N-1:-1:1
                    % compute shortcuts g,G,H
                    g = r(  :,n) + B(:,:,n)'*s(  :,n+1);
                    G = P(:,:,n) + B(:,:,n)'*S(:,:,n+1)*A(:,:,n);
                    H = R(:,:,n) + B(:,:,n)'*S(:,:,n+1)*B(:,:,n);
                    
                    % check for divergence
                    if any(~isfinite(H)), diverged=1; break; end % divergence: EXIT
                    
                    % find control law
                    [l(:,n), L(:,:,n)] = u_optimal(g,G,H,u(:,n),umin,umax,lambda);
                    
                    % update cost-to-go approximation
                    S (:,:,n) = Q (:,:,n) + A(:,:,n)'*S (:,:,n+1)*A(:,:,n) +    L(:,:,n)'*H*L(:,:,n) + L(:,:,n)'*G + G'*L(:,:,n);
                    s (  :,n) = q (  :,n) + A(:,:,n)'*s (  :,n+1)          +    L(:,:,n)'*H*l(  :,n) + L(:,:,n)'*g + G'*l(:,n);
                    s0(    n) = q0(    n) +           s0(    n+1)          + .5*l(  :,n)'*H*l(  :,n) + l(  :,n)'*g;
                    
                end
                if diverged, fprintf('Optimal control law and cost-to-go computation diverged. '); break; end% divergence: EXIT
                
                %------ STEP 3A: new control sequence, trajectory, cost
                % simulate linearized system to compute new control
                dx = zeros(dimX,1);
                unew = zeros(dimU,N-1);
                for n=1:N-1
                    du = l(:,n) + L(:,:,n)*dx;
                    du = min(max(du+u(:,n),umin),umax) - u(:,n);
                    dx = A(:,:,n)*dx + B(:,:,n)*du ;
                    unew(:,n) = u(:,n) + du;
                end
                
                %------ STEP 3B: simulate system to compute new trajectory and cost
                
                usim = scale_controlSeq(unew,t,tsim);
                xnew    = simulate_feedforward ( x0, f, unew, ps0 );
                
                xsim    = simulate_feedforward ( x0, f, usim, ps );
                
                costnew = evaluate_trajectory_cost_fh ( xsim, usim, j, ps );% initial cost
                
                %------ STEP 4: Levenberg-Marquardt method
                if costnew < cost
                    % decrease lambda (get closer to Newton method)
                    lambda = lambda / lambda_factor;
                    dcost  = cost-costnew; % decrease in cost
                    % update x, u, cost
                    u    = unew;
                    x    = xnew;
                    cost = costnew;
                    % check for convergence
                    if iter>1 && dcost<dcost_converge
                        if online_printing==1
                            fprintf('Cost improvement threshold reached. ');
                        end
                        break; % improvement too small: EXIT
                    end
                    % flag update
                    update = 1;
                else
                    % increase lambda (get closer to gradient descent)
                    lambda = lambda * lambda_factor;
                    if lambda>lambda_max
                        fprintf('Max. lambda reached. '); break; % lambda too large: EXIT
                    end
                end
                lambdas(iter)=lambda; costs(iter)=cost;
                % plot/print stuff out
                if online_printing==1,fprintf('Iteration = %d; Cost = %.4f; log(Lambda) = %.1f\n',iter,cost,log10(lambda));end
                if online_plotting==1
                    
                    set(0,'CurrentFigure',online_plotting_fig),
                    subplot(1,2,1),plot(0:iter-1,costs);
                    subplot(1,2,2),plot(0:iter-1,lambdas);
                end
            end
            
            if iter==iter_max
                fprintf('Max. number of iterations reached. ');
            end
            % print final result if necessary
            if online_printing==2
                fprintf('Iterations = %d;  Cost = %.4f', iter, cost);
            end
            
            fprintf('\n');
            % close figure
            if online_plotting==1, close(online_plotting_fig),end
            
            result.x = x;
            result.u = u;
            result.L = L;
            result.cost = cost;
            result.iteration = iter;
            result.lambdas = lambdas;
            result.costs = costs;
        end
        
        function [result] = ilqr2(f, j_t,j_f, dt, N, x0, u0, p )
            
            % unpack parameter struct
            if isfield(p,'lambda_init'        ),lambda_init         = p.lambda_init;
            else
                lambda_init        = 1e-2    ;
            end
            if isfield(p,'lambda_factor'      ),lambda_factor       = p.lambda_factor      ;
            else
                lambda_factor      = sqrt(10);
            end
            if isfield(p,'lambda_max'         ),lambda_max          = p.lambda_max         ;
            else
                lambda_max         = 1e-2    ;
            end
            if isfield(p,'dcost_converge'     ),dcost_converge      = p.dcost_converge  ;
            else
                dcost_converge     = 1e-9    ;
            end
            if isfield(p,'iter_max'           ),iter_max            = p.iter_max     ;
            else
                iter_max           = 100     ;
            end
            if isfield(p,'online_printing'    ),online_printing     = p.online_printing;
            else
                online_printing    = 1       ;
            end
            if isfield(p,'online_plotting'    ),online_plotting     = p.online_plotting;
            else
                online_plotting    = 1       ;
            end
            if isfield(p,'umax'               ),umax                = p.umax ;
            else
                umax               = inf     ;
            end
            if isfield(p,'umin'               ),umin                = p.umin     ;
            else
                umin               =-inf     ;
            end
            if isfield(p,'solver'             ),solver              = p.solver  ;
            else
                solver             ='rk4'  ;
            end
            
            % initialise everything
            dimX = size(x0, 1);   % state dimensionality
            dimU = size(u0, 1);   % command dimensionality
            ps=[];ps.dt=dt;ps.solver=solver; % simulation parameters
            ps.T = p.T;
            if size(u0,2)==1                                   % initial command sequence
                u = repmat(u0,1,N-1);                               % ...if one initial command given, use as full initial command sequence
            else                                               %
                u    = u0;
            end                                       % ...if full command sequence given, use as initial command sequence
            L    = zeros(dimU, dimX, N-1);                      % initial feedback gain sequence
            x    = simulate_feedforward ( x0, f, u, ps );       % initial state sequence
            cost = ILQRController.traj_cost_fh( x, u, j_t,j_f, ps ); % initial cost
            t    = (0:N-1)*dt;                                  % time
            
            % initialise other matrices (for speed)
            A  = zeros(dimX,dimX,N-1);
            B  = zeros(dimX,dimU,N-1);
            I  =   eye(dimX         );
            q0 = zeros(        1,N-1);
            q  = zeros(     dimX,N-1);
            Q  = zeros(dimX,dimX,N-1);
            r  = zeros(     dimU,N-1);
            R  = zeros(dimU,dimU,N-1);
            P  = zeros(dimU,dimX,N-1);
            s0 = zeros(        1,N  );
            s  = zeros(     dimX,N  );
            S  = zeros(dimX,dimX,N  );
            l  = zeros(     dimU,N-1);
            
            % initialise lambda
            lambda = lambda_init;
            % initialise update flag
            update = 1;
            
            lambdas = NaN(iter_max,1);
            costs = NaN(iter_max,1);
            
            % set up plotting
            if online_plotting
                online_plotting_fig=figure;clf,set(online_plotting_fig,'Name','iLQR: Cost and Lambda Convergence'),set(online_plotting_fig,'NumberTitle','off')
                subplot(1,2,1),hold on,grid on,set(gca,'Yscale','linear'),title('cost'  ),ylabel('cost'  ),xlabel('iteration')
                subplot(1,2,2),hold on,grid on,set(gca,'Yscale','log')   ,title('lambda'),ylabel('lambda'),xlabel('iteration')
            end
            
            % main loop
            for iter = 1:iter_max
                %------ STEP 1: approximate dynamics and cost along new trajectory'
                if update
                    update = 0;
                    % compute LQ approximation
                    for n = 1:N-1
                        % linearize dynamics, adjust for dt
                        %[ff, f_x, f_u] = f(x(:,n), u(:,n));
                        [~, f_x, f_u] = f(x(:,n), u(:,n));
                        A(:,:,n) = I + dt*f_x;
                        B(:,:,n) =     dt*f_u;
                        
                        % quadratize cost, adjust for dt
                        [l0,l_x,l_xx,l_u,l_uu,l_ux] = j_t(x(:,n), u(:,n), t(n));
                        q0(    n) = dt*l0;
                        q (  :,n) = dt*l_x;
                        Q (:,:,n) = dt*l_xx;
                        r (  :,n) = dt*l_u;
                        R (:,:,n) = dt*l_uu;
                        P (:,:,n) = dt*l_ux;
                    end
                    
                    % initialise value function approximation with final cost
                    [s0(N),s(:,N),S(:,:,N)] = j_f(x(:,N), u(:,N), t(N));
                end
                
                %------ STEP 2: compute optimal control law and cost-to-go
                diverged=0;
                for n = N-1:-1:1
                    % compute shortcuts g,G,H
                    g = r(  :,n) + B(:,:,n)'*s(  :,n+1);
                    G = P(:,:,n) + B(:,:,n)'*S(:,:,n+1)*A(:,:,n);
                    H = R(:,:,n) + B(:,:,n)'*S(:,:,n+1)*B(:,:,n);
                    
                    % check for divergence
                    if any(~isfinite(H)), diverged=1; break; end % divergence: EXIT
                    
                    % find control law
                    [l(:,n), L(:,:,n)] = u_optimal(g,G,H,u(:,n),umin,umax,lambda);
                    
                    % update cost-to-go approximation
                    S (:,:,n) = Q (:,:,n) + A(:,:,n)'*S (:,:,n+1)*A(:,:,n) +    L(:,:,n)'*H*L(:,:,n) + L(:,:,n)'*G + G'*L(:,:,n);
                    s (  :,n) = q (  :,n) + A(:,:,n)'*s (  :,n+1)          +    L(:,:,n)'*H*l(  :,n) + L(:,:,n)'*g + G'*l(:,n);
                    s0(    n) = q0(    n) +           s0(    n+1)          + .5*l(  :,n)'*H*l(  :,n) + l(  :,n)'*g;
                    
                end
                if diverged, fprintf('Optimal control law and cost-to-go computation diverged. '); break; end% divergence: EXIT
                
                %------ STEP 3A: new control sequence, trajectory, cost
                % simulate linearized system to compute new control
                dx = zeros(dimX,1);
                unew = zeros(dimU,N-1);
                for n=1:N-1
                    du = l(:,n) + L(:,:,n)*dx;
                    du = min(max(du+u(:,n),umin),umax) - u(:,n);
                    dx = A(:,:,n)*dx + B(:,:,n)*du ;
                    unew(:,n) = u(:,n) + du;
                end
                
                %------ STEP 3B: simulate system to compute new trajectory and cost
                xnew    = simulate_feedforward ( x0, f, unew, ps );
                costnew = ILQRController.traj_cost_fh ( xnew, unew, j_t,j_f, ps );% initial cost
                
                %------ STEP 4: Levenberg-Marquardt method
                if costnew<cost
                    % decrease lambda (get closer to Newton method)
                    lambda = lambda / lambda_factor;
                    dcost  = cost-costnew; % decrease in cost
                    % update x, u, cost
                    u    = unew;
                    x    = xnew;
                    cost = costnew;
                    % check for convergence
                    if iter>1 && dcost<dcost_converge
                        if online_printing==1
                            fprintf('Cost improvement threshold reached. ');
                        end
                        break; % improvement too small: EXIT
                    end
                    % flag update
                    update = 1;
                else
                    % increase lambda (get closer to gradient descent)
                    lambda = lambda * lambda_factor;
                    if lambda>lambda_max
                        fprintf('Max. lambda reached. '); break; % lambda too large: EXIT
                    end
                end
                lambdas(iter)=lambda; costs(iter)=cost;
                % plot/print stuff out
                if online_printing==1,fprintf('Iteration = %d; Cost = %.4f; log(Lambda) = %.1f\n',iter,cost,log10(lambda));end
                if online_plotting==1
                    
                    set(0,'CurrentFigure',online_plotting_fig),
                    subplot(1,2,1),plot(0:iter-1,costs);
                    subplot(1,2,2),plot(0:iter-1,lambdas);
                end
            end
            
            if iter==iter_max
                fprintf('Max. number of iterations reached. ');
            end
            % print final result if necessary
            if online_printing==2
                fprintf('Iterations = %d;  Cost = %.4f', iter, cost);
            end
            
            fprintf('\n');
            % close figure
            if online_plotting==1, close(online_plotting_fig),end
            
            result.x = x;
            result.u = u;
            result.L = L;
            result.cost = cost;
            result.iteration = iter;
            result.lambdas = lambdas;
            result.costs = costs;
        end
        
        function cost = traj_cost_fh( x, u, l_t,l_f, p)
            if nargin<4
                error 'Too few parameters given'
            end
            
            [~ ,Nu ] = size(u);
            [~, Nx] = size(x);
            
            if Nx~=Nu+1
                error 'Bad dimensionality of x and u'
            end
            
            dt = p.dt;
            t  = 0:dt:p.T;
            
            % integrate running cost & add final cost
            %cost = dt*sum(l(x(:,1:end-1),u,t)) + l(x(:,end),u,nan);
            cost_final = l_f(x(:,end),u(:,end),t(:,end));
            cost_run = zeros(1,Nu);
            for i=1:Nu
                cost_run(i) = l_t(x(:,i),u(:,i),t(i))*dt ;
            end
            cost = sum(cost_run) + cost_final;
        end
        
    end
    
    
end

% todo:
function [ l, l_x, l_xx, l_u, l_uu, l_ux ] = j_fd( fl, x,u,t )
%J_REACH_MW_MCC1DOF_MD cost with energy cost measured by mechanical work.
%Fot 1dof MACCEPA-VD with motor dynamics.
%   x : column vector, 8 elements
%   u : command input
%   t : time
%
%   by Fan WU 29 Jun 2016. fan.wu@kcl.ac.uk


% fl:
l = fl(x,u,t);

if nargout>1
    
    % finite difference
    flJ=@(x,u,t)J_cost_fd ( fl, x, u, t );
    [l_x ,l_u      ] = flJ ( x, u, t );
    flH =@(x,u,t)H_cost_fd  ( flJ, x, u, t );
    [l_xx,l_uu,l_ux] = flH  ( x, u, t );
else
    error('missing or wrong value of p.fd')
end


end

function [l,L] = u_optimal(g, G, H, u, umin, umax, lambda)

% eigenvalue decomposition, modify eigenvalues
[V,D] = eig(H);
d = diag(D);
d(d<0) = 0;
d = d + lambda;

% inverse modified Hessian, unconstrained control law
H1 = V*diag(1./d)*V';
l = -H1*g;
L = -H1*G;

% enforce constraints
l = min(max(l+u,umin),umax) - u;
L((l+u==umin)|(l+u==umax),:) = 0;
end

function x = simulate_feedforward ( x0, f, u, p )

N = size(u,2)+1; % N of time index
x = nan(size(x0,1),N); x(:,1)=x0; % initialise x
for n=1:N-1
	x(:,n+1) = simulate_step ( f, x(:,n), u(:,n), p );
end

end

function xn = simulate_step ( f, x, u, p )

switch p.solver
	case 'euler'
		dt=p.dt;
		xn = x + dt*f(x,u); % euler step
	case 'rk4'
		dt = p.dt;
        
		g1 = dt*f(x            ,u);
		g2 = dt*f(x+.5*g1,u);
		g3 = dt*f(x+.5*g2,u);
		g4 = dt*f(x+   g3,u);
		xn = x + (1/6)*(g1 + 2*g2 + 2*g3 + g4);
    case 'ode45'
        tspan = [0 p.dt];
        odefun = @(t,y)f(y,u);
        [~,y] = ode45(odefun,tspan,x);
        xn = y(end,:)';
	otherwise
		error('Unknown solver.')
end
end

function cost = evaluate_trajectory_cost_fh ( x, u, l, p )

if nargin<4
   error 'Too few parameters given'
end

[dimu ,Nu ] = size(u);
[dimx, Nx] = size(x);

if Nx~=Nu+1
   error 'Bad dimensionality of x and u'
end

dt = p.dt;
t  = 0:dt:p.T;

% integrate running cost & add final cost
%cost = dt*sum(l(x(:,1:end-1),u,t)) + l(x(:,end),u,nan);
cost_final = l(x(:,end), NaN, t(end));
cost_run = zeros(1,Nu);
for i=1:Nu
   cost_run(i) = l(x(:,i),u(:,i),t(i))*dt ; 
end
cost = sum(cost_run) + cost_final;

end
