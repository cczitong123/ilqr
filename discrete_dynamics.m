function xnew = discrete_dynamics(f, x, u, op)
psim.dt = op.sdt;
psim.solver='rk4';
xnew = simulate_step(f, x, u, psim);
end

