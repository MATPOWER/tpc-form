classdef math_model_pf_tpc < mp.math_model_pf & wgv.mm_shared_pfcpf_tpc
% wgv.math_model_pf_tpc - Power flow (PF) **math model** for TPC formulation.
%
% Provides formulation-specific and PF-specific subclasses for elements
% and implements formulation-specific node balance constraints.
%
% Overrides the default solve_opts() method.

%   MATPOWER
%   Copyright (c) 2021-2025, Power Systems Engineering Research Center (PSERC)
%   by Wilson Gonzalez Vanegas, Universidad Nacional de Colombia Sede Manizales
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://matpower.org for more info.

     properties
         
     end

    methods
        function obj = math_model_pf_tpc()
            %

            obj@mp.math_model_pf();
            obj.element_classes = { @wgv.mme_bus_pf_tpc, @wgv.mme_gen_pf_tpc, ...
                @wgv.mme_load_pf_tpc, @wgv.mme_branch_pf_tpc, @wgv.mme_shunt_pf_tpc};
        end

        function tag = form_tag(obj)
            %

            if obj.userdata.tpc.quad
                form = 'quad';
            else
                form = 'lin';
            end

            tag = ['tpc-' form];
        end

        function name = form_name(obj)
            %
            
            if obj.userdata.tpc.quad
                form = 'QUAD';
            else
                form = 'LIN';
            end

            name = ['TPC-' form];
        end

        function obj = add_node_balance_constraints(obj, nm, dm, mpopt)
            % Copy auxiliary data for the sake of sintaxis
            ad = obj.aux_data; obj.aux_data = [];
            
            % Index vector for pv and pq buses
            pvq = [ad.pv; ad.pq];

            % Vector of initial voltage angles and magnitudes
            obj.userdata.u0 = obj.init_u(dm, nm);       %% store u0 to userdata for subsequent utilization

            % Find ports of branches and shunts connected to each bus (for computing bus injection)
            ports = obj.ports_by_bus(nm);
            
            %% Define quadratic forms for complex power bus injections at PV-PQ buses
            id_vars = [pvq; ad.pq + nm.node.N];    % the vector of variables is [theta_pv; theta_pq; lnvm_pq]  (npv+2np variables)
            ad.id_vars = id_vars;
            
            [ad.Qbus_pvq, ad.Cbus_pvq, ad.kbus_pvq] = obj.bus_complex_injection(nm, ports(pvq), id_vars);   

            % Active power at PV buses            
            C_Ppv = real(ad.Cbus_pvq(1:ad.npv,:));
            k_Ppv = real(ad.kbus_pvq(1:ad.npv));
            Pbus_pv = real(ad.Sbus(ad.pv));

            % Active power at PQ buses            
            C_Ppq = real(ad.Cbus_pvq(ad.npv+1:ad.npv+ad.npq,:));
            k_Ppq = real(ad.kbus_pvq(ad.npv+1:ad.npv+ad.npq));
            Pbus_pq = real(ad.Sbus(ad.pq));

            % Reactive power at PQ buses            
            C_Qpq = imag(ad.Cbus_pvq(ad.npv+1:ad.npv+ad.npq,:));
            k_Qpq = imag(ad.kbus_pvq(ad.npv+1:ad.npv+ad.npq));
            Qbus_pq = imag(ad.Sbus(ad.pq));
            
            
            if obj.userdata.tpc.lin  % linear tpc-based formulation
                Q_Ppv = [];
                Q_Ppq = [];
                Q_Qpq = []; 
            else                     % quadratic tpc-based formulation
                Q_Ppv = cellfun(@(x)(real(x)), ad.Qbus_pvq(1:ad.npv), 'UniformOutput', false);
                Q_Ppq = cellfun(@(x)(real(x)), ad.Qbus_pvq(ad.npv+1:ad.npv+ad.npq), 'UniformOutput', false);
                Q_Qpq = cellfun(@(x)(imag(x)), ad.Qbus_pvq(ad.npv+1:ad.npv+ad.npq), 'UniformOutput', false);
            end

            %% Power balance constraints
            if obj.userdata.tpc.quad   %% quadratic tpc-based formulation
                
                % Build params for creating set of quadratic constraints
                Q = vertcat(Q_Ppv, Q_Ppq, Q_Qpq);
                C = [C_Ppv; C_Ppq; C_Qpq];
                k = [k_Ppv; k_Ppq; k_Qpq];
                l = [Pbus_pv; Pbus_pq; Qbus_pq] - k;
                
                obj.qcn.add(obj.var,'PQmis_pvq', Q, C, l, l);
            else                        %% linear tpc-based formulation

                % Build params for creating set of quadratic constraints
                C_P = [C_Ppv; C_Ppq];
                l_P = [Pbus_pv - k_Ppv; Pbus_pq - k_Ppq];
                l_Q = Qbus_pq - k_Qpq;                
                obj.lin.add(obj.var, 'Pmis_pvq', C_P, l_P, l_P);
                obj.lin.add(obj.var, 'Qmis_pq', C_Qpq, l_Q, l_Q);
            end

            %% resume auxiliary data
            obj.aux_data = ad;
        end

        function u0 = init_u(obj, dm, nm)
            %
            
            % Initialize vector of voltage angles and log-magnitudes
                
            u0 = spalloc(2 * nm.node.N, 1, 2 * nm.node.N);
            
            % %% 1-phase side
            if ~isempty(dm.source.bus)
                dme_bus = dm.elements.bus;
                u0(dme_bus.on) = deg2rad(dme_bus.tab.va(dme_bus.on));          %% \theta according to 1-phase bus info
                u0(dme_bus.on+dme_bus.n) = log(dme_bus.tab.vm(dme_bus.on));    %% ln(vm) according to 1-phase bus info
            end

            if ~isempty(dm.source.gen)
                dme_gen = dm.elements.gen;
                id_offset = nm.node.idx.N.bus;                   % [ \theta_1p ...]
                [gen_buses, ~] = find(nm.elements.gen.C); 
                u0(id_offset + gen_buses) = ...
                    log(dme_gen.tab.vm_setpoint(dme_gen.on));            %% ln(vm) according to 1-phase gen info
            end            
        end

        function ports = ports_by_bus(obj, nm)
            %
                
            ports = {};            

            %% 1-phase side
            % Find ports connected to each 1-phase bus: branches
            if nm.elements.has_name('branch')
                np = sum(nm.elements.branch.C,2);                         %% number of ports connected to each bus
                [ports, ~] = find(nm.elements.branch.C');
                ports = ports + min(nm.port.idx.i1.branch) - 1;       %% port numbers w.r.t. full network
                ports = mat2cell(ports, np);               
            end

            % Find ports connected to each 1-phase bus: shunts
            if nm.elements.has_name('shunt')
                [ports_shunt1p, id_bus_shunt] = find(nm.elements.shunt.C');
                ports_shunt1p = ports_shunt1p + min(nm.port.idx.i1.shunt) - 1;     %% port numbers w.r.t. full network
                
                for i = 1:length(id_bus_shunt)
                    b = id_bus_shunt(i);                    
                    ports{b} = [ports{b}; ports_shunt1p(i)];   %% Add shunt port to corresponding bus indices
                end           
            end                      
        end

        function [Q, C, k] = bus_complex_injection(obj, nm, ports, id_vars)
            %

            % Quadratic matrices
            if obj.userdata.tpc.quad
                [Q, MQ, kQ] = obj.rearrange_quadratic_terms(nm.Qu, ports, id_vars, 2 * nm.node.N);
            else
                Q = []; MQ = {0}; kQ = {0};
            end

            % Linear components
            [MM, kM] = obj.rearrange_linear_terms(nm.M, ports, id_vars, 2 * nm.node.N);
            C = cell2mat(MM) + cell2mat(MQ);

            % Constant components
            kk = cellfun(@(x)(sum(nm.s(x))), ports, 'UniformOutput', false);
            k = cell2mat(kk) + cell2mat(kM) + cell2mat(kQ);
        end

        function [QQ, MQ, kQ] = rearrange_quadratic_terms(obj, W, id_W, id_vars, nvars)
            % 

            if obj.userdata.tpc.lin
                QQ = [];
                MQ = [];
                kQ = [];
            else
                % vector of initial voltage angles and magnitudes
                u0 = obj.userdata.u0';

                not_id_vars = setdiff((1:nvars)',id_vars); %% indices of known voltage angles ad magnitudes

                if nargin < 5
                    if ~isempty(id_vars)
                        error('wgv.math_model_pf_tpc.sum_quadratic_matrices: arguments ''idx'' and ''nvars'' must be provided \n');
                    end
                end
                
                % sum matrices in the form of [row col val] (vertical stack)
                Qid_W = cellfun(@(x)(vertcat(W{x})), id_W, 'UniformOutput', false);

                if ~isempty(id_vars) %% Sparse matrices
                    
                    % Sum quadratic matrices bus injections via sparse command
                    Qid_W = cellfun(@(x)(sparse(x(:,1), x(:,2), x(:,3), nvars, nvars)), Qid_W, 'UniformOutput', false);

                    % ------ Quadratic matrices ------                    
                    QQ = cellfun(@(x)(x(id_vars,id_vars)), Qid_W, 'UniformOutput', false);    % It leads to twice the linear terms but it calcels out due to (1/2)*x'*Q*x

                    % ------ Linear terms from quadratic matrices------
                    MQ = cellfun(@(x)(u0(not_id_vars)*x(not_id_vars, id_vars)), Qid_W, 'UniformOutput', false);
                    
                    % ------ Constant terms from quadratic matrices ------
                    kQ = cellfun(@(x)((1/2)*u0(not_id_vars)*x(not_id_vars, not_id_vars)*u0(not_id_vars)'), Qid_W, 'UniformOutput', false); 
                else
                    QQ = Qid_W;   %% 3-column arrays of the form [row col val]

                    % MQ and kQ calculation in [row col val]-form to be implemented!!
                end
            end
        end

        function [MM, kM] = rearrange_linear_terms(obj, M, id_M, id_vars, nvars)
            %

            % vector of initial voltage angles and magnitudes
            u0 = obj.userdata.u0;

            not_id_vars = setdiff((1:nvars)',id_vars); %% indices of known voltage angles ad magnitudes

            % ------ Linear terms from matrix M ------
            MM = cellfun(@(x)(sum(M(x, id_vars), 1)), id_M, 'UniformOutput', false);

            % ------ Constant terms from matrix M ------
            kM = cellfun(@(x)(sum(M(x, not_id_vars), 1)*u0(not_id_vars)), id_M, 'UniformOutput', false);
        end
    
        function [f, J] = node_balance_equations(obj, x)
            %
            
            if nargout > 1    % Return the Jacobian               
               [f, ~, ~, J] = obj.qcn.eval(obj.var, x, 'PQmis_pvq');
            else              % Return f               
               f = obj.qcn.eval(obj.var, x, 'PQmis_pvq');
            end
        end


        % -----------------------------------------------------------------
        % ************ OVERRIDDEN METHODS FROM PARENT CLASSES ************
        % -----------------------------------------------------------------

        function nm = network_model_x_soln(obj, nm)
            % THIS METHOD OVERRIDES THE CORRESPONDING PARENT METHOD. 
            % Convert solved state from math model to network model solution.
            % ::
            %
            %   nm = mm.network_model_x_soln(nm)
            %
            % Input:
            %   nm (mp.net_model) : network model object
            %
            % Output:
            %   nm (mp.net_model) : updated network model object
            %
            % Calls convert_x_m2n() which is defined in a subclass of
            % mp.mm_shared_pfcpf allowing it to be shared across math 
            % models for different tasks (PF and CPF).
            
            [nm.soln.u, nm.soln.z, nm.soln.x] = ...    % sonln.u (difference with parent)
                obj.convert_x_m2n(obj.soln.x, nm);     % third argument not passed to update zvars
        end
        
        function opt = solve_opts(obj, nm, dm, mpopt)
            %
            
            if nm.userdata.tpc.lin
                opt = struct( ...
                    'verbose',  mpopt.verbose, ...
                    'leq_opt',  struct('thresh', 1e5) );
            else
                opt = struct( ...
                    'verbose',  mpopt.verbose, ...
                    'max_it', mpopt.pf.nr.max_it);                
            end
        end

        function display(obj)
            % Display the math model object.
            %
            % This method is called automatically when omitting a semicolon
            % on a line that retuns an object of this class.
            %
            % Displays the details of the variables, constraints, costs, and
            % math model elements.
            %
            % See also mp_idx_manager.

            if obj.userdata.tpc.lin 
                ftype = 'lin';
            elseif obj.userdata.tpc.quad
                ftype = 'quad';
            else
                ftype = obj.form_tag();
            end

            fprintf('MATH MODEL CLASS : %s\n', class(obj));
            fprintf('    TASK NAME               : %s\n', obj.task_name());
            fprintf('    TASK TAG                : %s\n', obj.task_tag());
            fprintf('    FORMULATION NAME        : %s\n', obj.form_name());
            fprintf('    FORMULATION TAG         : %s\n', obj.form_tag());
            fprintf('    FORMULATION TYPE        : %s\n', ftype);
            display@opt_model(obj)

            %% elements
            fprintf('\nELEMENTS\n')
            fprintf('========\n')
            fprintf('  name              class\n');
            fprintf(' ----------------  --------------------\n');
            for k = 1:length(obj.elements)
                mme = obj.elements{k};
                fprintf('  %-13s     %s\n', mme.name, class(mme));
            end
        end


    end     %% methods
end         %% classdef