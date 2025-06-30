classdef task_pf < mp.task
% wgv.task_pf - |MATPOWER| task for power flow (PF) with TPC 
% formulation (wgv.form_tpc).
%
% Provides task implementation for the power flow problem.
%
% This includes the handling of iterative runs to enforce generator
% reactive power limits, if requested.
%
% mp.task_pf Properties:
%   * tag - task tag 'PF'
%   * name - task name 'Power Flow'
%   * dc - ``true`` if using DC network model
%   * lnv_quad  ``true`` if using LNV_QUAD network model
%   * iterations - total number of power flow iterations
%   * ref - current ref node indices
%   * ref0 - initial ref node indices
%   * va_ref0 - initial ref node voltage angles
%   * fixed_q_idx - indices of fixed Q gens
%   * fixed_q_qty - Q output of fixed Q gens
%
% mp.task_pf Methods:
%   * run_pre - set :attr:`dc` property
%   * next_dm - optionally iterate to enforce generator reactive limits
%   * enforce_q_lims - implementation of generator reactive limits
%   * network_model_class_default - select default network model constructor
%   * network_model_build_post - initialize properties for reactive limits
%   * network_model_x_soln - correct the voltage angles if necessary
%   * math_model_class_default - select default math model constructor
%
% See also mp.task.

%   MATPOWER
%   Copyright (c) 2021-2025, Power Systems Engineering Research Center (PSERC)
%   by Wilson Gonzalez Vanegas, Universidad Nacional de Colombia Sede Manizales
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://matpower.org for more info.


    properties
        tag = 'PF';             % 
        name = 'Power Flow';    % 
        dc        % ``true`` if using DC network model (from ``mpopt.model``, cached in run_pre())
        lin       % ``true`` if using TPC formulation with linear network model (from ``mpopt.model``, cached in run_pre())
        quad      % ``true`` if using TPC formulation with quadratic network model (from ``mpopt.model``, cached in run_pre())
        iterations              % *(integer)* total number of power flow iterations
        ref                     % *(integer)* current ref node indices
        ref0                    % *(integer)* initial ref node indices
        va_ref0                 % *(double)* initial ref node voltage angles
        fixed_q_idx             % *(integer)* indices of fixed Q gens
        fixed_q_qty             % *(double)* Q output of fixed Q gens
    end

    methods
        %%-----  task methods  -----
        function dm = next_dm(obj, mm, nm, dm, mpopt, mpx)
            % Implement optional iterations to enforce generator reactive
            % limits.
            if ~obj.dc && mpopt.pf.enforce_q_lims
                %% adjust iteration count for previous runs
                obj.iterations = obj.iterations + mm.soln.output.iterations;
                mm.soln.output.iterations = obj.iterations;

                %% enforce Q limits
                [success, dm] = obj.enforce_q_lims(nm, dm, mpopt);
                if ~success                 %% entire task fails if Q lim
                    obj.success = success;  %% enforcement indicates failure
                end
            else        %% don't enforce generator Q limits, once is enough
                dm = [];
            end
        end

        function [success, dm] = enforce_q_lims(obj, nm, dm, mpopt)
            % Used by next_dm() to implement enforcement of generator
            % reactive limits.
            gen_dme = dm.elements.gen;
            [mn, mx, both] = gen_dme.violated_q_lims(dm, mpopt);

            if ~isempty(both)   %% we have some Q limit violations
                if isempty(mn) && isempty(mx)   %% infeasible
                    if mpopt.verbose
                        fprintf('All %d remaining gens exceed their Q limits : INFEASIBLE PROBLEM\n', length(both));
                    end
                    dm = [];
                    success = 0;
                else
                    if mpopt.verbose && ~isempty(mx)
                        fprintf('Gen %d at upper Q limit, converting to PQ bus\n', gen_dme.on(mx));
                    end
                    if mpopt.verbose && ~isempty(mn)
                        fprintf('Gen %d at lower Q limit, converting to PQ bus\n', gen_dme.on(mn));
                    end

                    %% save corresponding limit values
                    obj.fixed_q_qty(mx) = gen_dme.qg_ub(mx);
                    obj.fixed_q_qty(mn) = gen_dme.qg_lb(mn);
                    mx = [mx;mn];

                    %% set qg to binding limit
                    gen_dme.tab.qg(gen_dme.on(mx)) = ...
                        obj.fixed_q_qty(mx) * dm.base_mva;

                    %% convert to PQ bus
                    bus_dme = dm.elements.bus;
                    ref0 = find(bus_dme.type == mp.NODE_TYPE.REF);
                    bidx = bus_dme.i2on(gen_dme.bus(gen_dme.on(mx)));   %% bus of mx
                    if length(ref0) > 1 && any(bus_dme.type(bidx) == mp.NODE_TYPE.REF)
                        error('mp.task_pf.enforce_q_lims: Sorry, MATPOWER cannot enforce Q limits for slack buses in systems with multiple slacks.');
                    end
                    %% set bus type to PQ
                    bus_dme.set_bus_type_pq(dm, bidx);
                    %% potentially pick new slack bus
                    ntv = nm.node_types(nm, dm);        %% node type vector
                    [i1, iN] = nm.get_node_idx('bus');  %% bus node indices
                    btv = ntv(i1:iN);                   %% bus type vector

                    %% indicate if there's been a change in slack bus
                    ref = find(btv == mp.NODE_TYPE.REF);    %% new ref bus indices
                    if mpopt.verbose && ref ~= ref0
                        fprintf('Bus %d is new slack bus\n', ...
                            bus_dme.ID(bus_dme.on(ref)));
                    end

                    %% save indices to list of Q limited gens
                    obj.fixed_q_idx = [obj.fixed_q_idx; mx];

                    %% update dm for next step
                    dm.initialize();
                    dm.update_status();
                    dm.build_params();
                    success = 1;
                end
            else                %% no more Q violations
                dm = [];
                success = 1;
            end
        end

        %%-----  data model methods  -----
        function dm = data_model_build_post(obj, dm, dmc, mpopt)
            % Ensures the existence of all fields for 1-phase and 3-phase
            % data in the .source field of the data model object

            mpc = dm.source;    % Matpower case data

            %% 1-phase model data
            if ~isfield(mpc, 'bus')
                mpc.bus = [];
            end
            if ~isfield(mpc, 'gen')
                mpc.gen = [];
            end
            if ~isfield(mpc, 'branch')
                mpc.branch = [];
            end
            if ~isfield(mpc, 'bus')
                mpc.bus = [];
            end
            if ~isfield(mpc, 'gencost')
                mpc.gencost = [];
            end            

            % update Matpower case
            dm.source = mpc;
        end

        %%-----  network model methods  -----
        function nm_class = network_model_class_default(obj, dm, mpopt)% OK-WGV
            % Implement selector for default network model constructor
            % depending on ``mpopt.model`` and ``mpopt.pf.v_cartesian``.
            switch upper(mpopt.model)
                case 'AC'
                    if mpopt.pf.v_cartesian
                        nm_class = @mp.net_model_acc;
                    else
                        nm_class = @mp.net_model_acp;
                    end
                case 'DC'
                    nm_class = @mp.net_model_dc;
                case 'TPC'
                    nm_class = @wgv.net_model_tpc;
            end
        end       

        function nm = network_model_build_post(obj, nm, dm, mpopt)
            % Initialize mp.task_pf properties, if non-empty AC case with
            % generator reactive limits enforced.
            if ~obj.dc && mpopt.pf.enforce_q_lims ~= 0 && nm.np ~= 0
                if obj.i_nm == 1
                    [ref, ~, ~] = nm.node_types(obj, dm);
                    gen_dme =  dm.elements.gen;
                    obj.iterations = 0;
                    obj.ref0 = ref;             %% initial ref node indices
                    obj.ref = ref;              %% current ref node indices
                    obj.va_ref0 = nm.get_va(ref);%% initial ref node voltage angles
                    obj.fixed_q_idx = [];       %% indices of fixed Q gens
                    obj.fixed_q_qty = zeros(gen_dme.n, 1);  %% Q output of fixed Q gens
                else        %% update index of ref bus
                    [obj.ref, ~, ~] = nm.node_types(obj, dm);
                end
            end
        end

        function nm = network_model_x_soln(obj, mm, nm)
            % Call superclass :meth:`network_model_x_soln() <mp.task.network_model_x_soln>`
            % then correct the voltage angle if the ref node has been changed.
            nm = network_model_x_soln@mp.task(obj, mm, nm);

            if ~obj.dc && obj.i_nm > 1 && obj.ref ~= obj.ref0
                nn = nm.get_idx('node');
                va = nm.soln.u(nn.i1.bus:nn.iN.bus);
                va = va - va(obj.ref0) + obj.va_ref0;
                nm.soln.u(nn.i1.bus:nn.iN.bus) = va;
            end        
        end

        %%-----  mathematical model methods  -----
        function mm_class = math_model_class_default(obj, nm, dm, mpopt)
            % Implement selector for default mathematical model constructor
            % depending on ``mpopt.model``, ``mpopt.pf.v_cartesian``, and
            % ``mpopt.pf.current_balance``.
            switch upper(mpopt.model)
                case 'AC'
                    if mpopt.pf.v_cartesian
                        if mpopt.pf.current_balance
                            mm_class = @mp.math_model_pf_acci;
                        else
                            mm_class = @mp.math_model_pf_accs;
                        end
                    else
                        if mpopt.pf.current_balance
                            mm_class = @mp.math_model_pf_acpi;
                        else
                            mm_class = @mp.math_model_pf_acps;
                        end
                    end
                case 'DC'
                    mm_class = @mp.math_model_pf_dc;

                case 'TPC'
                    mm_class = @wgv.math_model_pf_tpc;
            end
        end
        
        % -----------------------------------------------------------------
        % ************ OVERRIDDEN METHODS FROM PARENT CLASSES ************
        % -----------------------------------------------------------------
        
        function [d, mpopt] = run_pre(obj, d, mpopt)% OK-WGV
            % Set :attr:`dc` property after calling superclass
            % :meth:`run_pre() <mp.task.run_pre>`.
            [d, mpopt] = run_pre@mp.task(obj, d, mpopt);    %% call parent

            %% cache DC model flag
            obj.dc = strcmp(upper(mpopt.model), 'DC');

            %% cache LIN or QUAD flag for TPC model
            if strcmp(upper(mpopt.model), 'TPC')
                if isfield(mpopt.pf, 'tpc')
                    obj.lin = strcmp(mpopt.pf.tpc.form, 'LIN');
                    obj.quad = strcmp(mpopt.pf.tpc.form, 'QUAD');
                else
                    error('\n wgv.task_pf: field ''tpc'' in MPOPTION.PF not found \n')
                end
            else
                obj.lin  = 0;
                obj.quad = 0;
            end           
        end

        function nm = network_model_build_pre(obj, nm, dm, mpopt)
            % THIS METHOD OVERRIDES THE PARENT'S METHOD (See mp.task_pf)
            % 
            % This method adds fields .lin and .quad to the property 
            % 'userdata' of the network model according to this same 
            % properties of the task object

            if ~obj.dc
                nm.userdata.tpc.lin  = obj.lin;
                nm.userdata.tpc.quad = obj.quad;
            end
        end

        function nm = network_model_update(obj, mm, nm)
            % THIS METHOD OVERRIDES THE CORRESPONDING PARENT METHOD.
            % BASICALLY, IT ALLOWS PASSING THE MATH MODEL TO PORT_INJ_SOLN
            % Update network model state, solution values from math model solution.
            % ::
            %
            %   nm = task.network_model_update(mm, nm)
            %
            % Inputs:
            %   mm (mp.math_model) : mathmatical model object
            %   nm (mp.net_model) : network model object
            %
            % Output:
            %   nm (mp.net_model) : updated network model object
            %
            % Called by run() method.

            %% save network state solution (convert from math model state)
            obj.network_model_x_soln(mm, nm);

            %% save port injection solution
            
            if obj.quad || obj.lin
                nm.port_inj_soln(mm);    %% mm is passed to have access to the methods of subclass mp.sm_quad_constraint
            else
                nm.port_inj_soln();
            end
        end
       
    end     %% methods
end         %% classdef
