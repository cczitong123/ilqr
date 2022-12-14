% Function for evaluating the cost along a trajectory under a given cost
% function and finite horizen
% 
% in:
%     x     - matrix containing trajectory through state space
%     u     - matrix containing trajectory of commands
%     l     - function handle to cost function
%     p     - simulation struct
%
% out:
%     cost   - cost incurred
%
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

