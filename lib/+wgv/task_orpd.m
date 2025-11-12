classdef task_orpd < mp.task
% wgv.task_orpd - |MATPOWER| task for optimal reactive power dispatch (ORPD).
%
% Provides task implementation for the optimal reactive power dispatch problem.
%
% wgv.task_orpd Properties:
%   * tag - task tag 'ORPD'
%   * name - task name 'Optimal Reactive Power Dispatch'
%   * dc - ``true`` if using DC network model
%   * lnv_quad  ``true`` if using LNV_QUAD network model
%
% wgv.task_orpd Methods:
%   * run_pre - set :attr:`dc` property
%   * print_soln_header - add printout of objective function value
%   * data_model_class_default - select default data model constructor
%   * data_model_build_post - adjust bus voltage limits, if requested
%   * network_model_class_default - select default network model constructor
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
        tag = 'ORPD';
        name = 'Optimal Reactive Power Dispatch';
        dc  % ``true`` if using DC network model (from ``mpopt.model``, cached in run_pre())
        lin       % ``true`` if using TPC formulation with linear network model (from ``mpopt.model``, cached in run_pre())
        quad      % ``true`` if using TPC formulation with quadratic network model (from ``mpopt.model``, cached in run_pre())
    end

    methods
        %%-----  task methods  -----
        function [d, mpopt] = run_pre(obj, d, mpopt)
            % Set :attr:`dc` property after calling superclass
            % :meth:`run_pre() <mp.task.run_pre>`, then check for
            % unsupported AC OPF solver selection.
            [d, mpopt] = run_pre@mp.task(obj, d, mpopt);     %% call parent            

            %% cache LIN or QUAD flag for TPC model
            if strcmp(upper(mpopt.model), 'TPC')
                if isfield(mpopt.orpd, 'tpc')
                    obj.lin = strcmp(mpopt.orpd.tpc.form, 'LIN');
                    obj.quad = strcmp(mpopt.orpd.tpc.form, 'QUAD');
                else
                    error('\n wgv.task_orpd: field ''tpc'' in MPOPTION.PF not found \n')
                end
            else
                error('wgv.task_orpd: current implementation only works with a linear/quadratic TPC model. Other formulations such as AC or DC are not yet supported.')
            end

            %% check for unsupported QTPC ORPD solver selection
            if obj.quad
                alg = upper(mpopt.orpd.tpc.solver);
                switch alg
                    case 'IPOPT'
                        if ~have_feature('ipopt')
                            error('mp.task_opf.run_pre: MPOPT.opf.ac.solver = ''%s'' requires IPOPT (see https://github.com/coin-or/Ipopt)', alg);
                        end
                    case 'FMINCON'
                        if ~have_feature('fmincon')
                            error('mp.task_opf.run_pre: MPOPT.opf.ac.solver = ''%s'' requires FMINCON (Optimization Toolbox 2.x or later)', alg);
                        end
                    case 'KNITRO'
                        if ~have_feature('knitro')
                            error('mp.task_opf.run_pre: MPOPT.opf.ac.solver = ''%s'' requires Artelys Knitro (see https://www.artelys.com/solvers/knitro/)', alg);
                        end
                    case {'MINOPF', 'PDIPM', 'TRALM', 'SDPOPF'}
                        error('mp.task_opf.run_pre: MPOPT.opf.ac.solver = ''%s'' not supported.', alg);
                end
            end
        end

        function print_soln_header(obj, mpopt, fd)
            % Call superclass :meth:`print_soln_header() <mp.task.print_soln_header>`
            % the print out the objective function value.
            if nargin < 3
                fd = 1;     %% print to stdio by default
            end

            print_soln_header@mp.task(obj, mpopt, fd);
            fprintf(fd, ...
                'Objective Function Value = %.2f $/hr\n', obj.mm.soln.f);
        end

        %%-----  data model methods  -----
        function dm_class = data_model_class_default(obj)
            % Implement selector for default data model constructor.

            dm_class = @wgv.data_model_orpd;
        end

        function dm = data_model_build_post(obj, dm, dmc, mpopt)
            % Call superclass :meth:`data_model_build_post() <mp.task.data_model_build_post>`
            % then adjust bus voltage magnitude limits based on generator
            % ``vm_setpoint``, if requested.

            %% call parent
            dm = data_model_build_post@mp.task(obj, dm, dmc, mpopt);

            if ~obj.dc
                %% if requested, adjust bus voltage magnitude
                %% limits based on generator vm_setpoint
                use_vg = mpopt.opf.use_vg;
                if use_vg
                    dm.set_bus_v_lims_via_vg(use_vg);
                end
            end
        end

        %%-----  network model methods  -----
        function nm_class = network_model_class_default(obj, dm, mpopt)
            % Implement selector for default network model constructor
            % depending on a LIN or QUAD form for TPC models.
            if obj.lin || obj.quad
                nm_class = @wgv.net_model_tpc;
            else
                error('task_orpd.network_model_class_default: curret implementation of ORPD tasks only works for TPC-based formulations.')
            end
        end

        %%-----  mathematical model methods  -----
        function mm_class = math_model_class_default(obj, nm, dm, mpopt)
            % Implement selector for default network model constructor
            % depending on a LIN or QUAD form for TPC models.
            if obj.lin || obj.quad
                mm_class = @wgv.math_model_orpd_tpc;
            else
                error('task_orpd.math_model_class_default: curret implementation of ORPD tasks only works for TPC-based formulations.')
            end
        end

        % -----------------------------------------------------------------
        % ************ OVERRIDDEN METHODS FROM PARENT CLASSES ************
        % -----------------------------------------------------------------

        function nm = network_model_build_pre(obj, nm, dm, mpopt)
            % THIS METHOD OVERRIDES THE PARENT'S METHOD (See mp.task_pf)
            %
            % This method adds fields .lin and .quad to the property
            % 'userdata' of the network model according to this same
            % properties of the task object

            if obj.lin || obj.quad
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
