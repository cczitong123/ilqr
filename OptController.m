classdef OptController
    %OCreach optimal control for reaching task
    %   A wrapper of all optimal control functionalities
    
    properties
        robot_model
        %:
        cdt
        dt
        Nt
        umax
        umin
        qmax
        qmin
        
        task % task object
        %task_type
        %task_param % target_q,  T, Nt,with_motor_dynamics = 1;
        %global_optimal = 1;
        %opt_solver
        
        sim_param = [];
        cost_param= [];
        % - w
        % - rege_ratio
        
        opt_param = []; % param for optimization algorithm
        
        %%%% task definition
        position0
        x0
        u0
        
        %%%%
        %%%%
        
        
        %%%%
        
        %j % running cost function handler
        
        
        result
    end
    
    methods
        function oc = OptController(robot_model, task)
            % define dynamics function, robot model object, cost functions,
            % optimisation algorithms
            % with motor dynamics, carefully choose false, since
            % robot dynamics now prefer with motor dynamics,
            % without motor dynamics model may be not useable anymore
            oc.robot_model = robot_model;
            oc.task = task;
            
            
            oc.cdt = robot_model.cdt;
            oc.dt = robot_model.cdt;
            oc.umax = robot_model.umax;
            oc.umin = robot_model.umin;
            oc.qmax = robot_model.qmax;
            oc.qmin = robot_model.qmin;
            
            oc.x0 = task.task_param.x0;
            
            
            
            
            if(isfield(task.task_param, 'sim_solver'))
                oc.sim_param.solver = task_param.sim_solver;
            else
                oc.sim_param.solver = 'rk4';%default to rk4
            end
            
            oc.sim_param.dt = oc.robot_model.cdt;
            oc.sim_param.N = oc.task.task_param.Nt;
            oc.Nt = oc.task.task_param.Nt;
        end
    end
    
end

