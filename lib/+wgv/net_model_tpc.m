classdef net_model_tpc < mp.net_model & wgv.form_tpc
% wgv.net_model_tpc - Concrete class for |MATPOWER| TPC **network model** objects.
%
% This network model class and all of its network model element classes are
% specific to the TPC formulation and therefore inherit from wgv.form_tpc.
%
% wgv_net_model_tpc Properties:
%   * va - vector of voltage states (voltage angles :math:`\va`)
%   * lnvm - vector of voltage states (log of voltage magnitudes :math:`\ln(\vm)`)
%   * zr - vector of real part of complex non-voltage states, :math:`\zr`
%   * zi - vector of imaginary part of complex non-voltage states, :math:`\zi`
%
% net_model_tpc Methods:
%   * net_model_tpc - constructor, assign default network model element classes
%   * def_set_types - add voltage and non-voltage variable set types for mp_idx_manager
%   * build_params - build incidence matrices, parameters, add ports for each element
%   * port_inj_soln - compute the network port injections at the solution
%
% See also mp.net_model, wgv.form_tpc, mp.form, mp.nm_element.

%   MATPOWER
%   Copyright (c) 2021-2025, Power Systems Engineering Research Center (PSERC)
%   by Wilson Gonzalez Vanegas, Universidad Nacional de Colombia Sede Manizales
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://matpower.org for more info.

    properties
        va = [];
        lnvm = [];
        zr  = [];
        zi  = [];
    end

    methods
        function obj = net_model_tpc()
            % Constructor, assign default network model element classes.
            % ::
            %
            %   nm = net_model_tpc()
            %
            % This network model class and all of its network model element
            % classes are specific to the TPC formulation and therefore
            % inherit from wgv.form_tpc.

            obj@mp.net_model();
            obj.element_classes = { ...
                @wgv.nme_bus_tpc, @wgv.nme_gen_tpc, ...
                @wgv.nme_load_tpc, @wgv.nme_branch_tpc, ...
                @wgv.nme_shunt_tpc };

            %% Due to a bug related to inheritance in constructors in
            %% Octave 5.2 and earlier (https://savannah.gnu.org/bugs/?52614),
            %% INIT_SET_TYPES() cannot be called directly in the
            %% MP_IDX_MANAGER constructor, as desired.
            %%
            %% WORKAROUND:  INIT_SET_TYPES() is called explicitly as needed
            %%              (if obj.node is empty) in BUILD() and DISPLAY(),
            %%              after object construction, but before object use.
        end

        function obj = def_set_types(obj)
            % Add voltage and non-voltage variable set types for mp_idx_manager.
            % ::
            %
            %   nm.def_set_types()
            %
            % Add the following set types:
            %
            %   - ``'vm'``   - VOLTAGE ANG VARS (va)
            %   - ``'lnvm'`` - LOG-VOLTAGE MAG VARS (lnvm)
            %   - ``'zr'``   - NON-VOLTAGE VARS REAL (zr)
            %   - ``'zi'``   - NON-VOLTAGE VARS IMAG (zi)
            %
            % See also mp.net_model.def_set_types, mp_idx_manager.

            def_set_types@mp.net_model(obj);        %% call parent first
            obj.set_types.va   = 'VOLTAGE ANG VARS (va)';
            obj.set_types.lnvm = 'LN-VOLTAGE MAG VARS (lnvm)';
            obj.set_types.zr   = 'NON-VOLTAGE VARS REAL (zr)';
            obj.set_types.zi   = 'NON-VOLTAGE VARS IMAG (zi)';
        end

        function obj = build_params(obj, nm, dm)
            % Build incidence matrices and parameters, and add ports for each element.
            % ::
            %
            %   nm.build_params(nm, dm)
            %
            % Inputs:
            %   nm (mp.net_model) : network model object
            %   dm (mp.data_model) : data model object
            %
            % Call the parent method to do most of the work, then build
            % the aggregate network model parameters.

            %% call parent to build individual element parameters
            build_params@mp.net_model(obj, nm, dm);

            %% aggregate parameters from individual elements
            obj.Qu = obj.stack_matrix_params('Qu', 1);
            obj.Qz = obj.stack_matrix_params('Qz', 1);
            obj.M  = obj.stack_matrix_params('M', 0);
            obj.N  = obj.stack_matrix_params('N', 0);            
            obj.s  = obj.stack_vector_params('s');
        end

        function [G, Gva, Glnvm, Gzr, Gzi] = nodal_complex_power_balance(obj, x, mm)
            % Compute nodal complex power balance constraints.
            % ::
            %
            %   G = nm.nodal_complex_power_balance(x)
            %   [G, Gv1, Gv2, Gzr, Gzi] = nm.nodal_complex_power_balance(x)
            %
            % Compute constraint function and optionally the Jacobian for
            % the complex power balance equality constraints based on
            % outputs of mp.form_ac.port_inj_power and the node incidence
            % matrix.
            %
            % Input:
            %   x (double) : state vector of real vars ( [u; z] )
            %   mm (object): mathematical model object
            %
            % Outputs:
            %   G (complex double) : nodal complex power balance constraint
            %       function, :math:`\G^\mathrm{kcl}(\x)`
            %   Gva (complex double) : Jacobian w.r.t. voltage angles,
            %       :math:`\G^\mathrm{kcl}_\Va` or :math:`\G^\mathrm{kcl}_\Vr`
            %   Gv2 (complex double) : Jacobian w.r.t. log of voltage magnitudes,
            %       :math:`\G^\mathrm{kcl}_\Vm` or :math:`\G^\mathrm{kcl}_\Vi`
            %   Gzr (complex double) : Jacobian w.r.t. real non-voltage variable,
            %       :math:`\G^\mathrm{kcl}_\zr`
            %   Gzi (complex double) : Jacobian w.r.t. imaginary non-voltage variable,
            %       :math:`\G^\mathrm{kcl}_\zi`
            %
            % See also mp.form_ac.port_inj_power,
            % nodal_complex_power_balance_hess.

            %% real-value vector of ports to block diagonal matrix
            [U_blk, Z_blk] = mm.uz2blkdiag(x, obj);

            %% node incidence matrix
            C = obj.C;

            %% get port power injections with derivatives
            if nargout > 1
                [S, Sva, Slnvm, Szr, Szi] = obj.port_inj_power(U_blk, Z_blk, mm.aux_data);
                Gva   = C * Sva;
                Glnvm = C * Slnvm;
                Gzr   = C * Szr;
                Gzi   = C * Szi;
            else
                S = obj.port_inj_power(U_blk, Z_blk, mm.aux_data);
            end

            %% nodal power balance
            G = C * S;
        end

        function d2G = nodal_complex_power_balance_hess(obj, x, lam, mm)
            % Compute nodal complex power balance Hessian.
            % ::
            %
            %   d2G = nm.nodal_complex_power_balance_hess(x_, lam)
            %
            % Compute the Hessian of the nodal complex power balance
            % constraint. Rather than a full, 3-dimensional Hessian, it
            % computes the Jacobian of the vector obtained by muliplying the
            % transpose of the constraint Jacobian by a vector :math:`\lam`.
            % Based on mp.form_ac.port_inj_power_hess.
            %
            % Inputs:
            %   x (double) : state vector of real vars ( [u; z] )
            %   lam (double) : vector :math:`\lam` of multipliers, one for each node
            %   mm (math model): mathematical model object
            %
            % Output:
            %   d2G (complex double) : sparse Hessian matrix,
            %       :math:`\G^\mathrm{kcl}_{\x\x}(\lam)`
            %
            % See also mp.form_ac.port_inj_power_hess,
            % nodal_complex_power_balance.

            %% get port power injection hessians
            d2G = obj.port_inj_power_hess(x, obj.C' * lam, mm.aux_data);
        end


        function obj = port_inj_soln(obj, mm)
            % Compute the network port injections at the solution.
            % ::
            %
            %  nm.port_inj_soln(mm)
            %
            % Inputs:
            %   mm (math model): mathematical model object
            %
            % Takes the solved network state, computes the port complex 
            % power injections, and saves them in ``nm.soln.gs_``.

            
            % Extract solution and create variables for z-vars
            u = obj.soln.u;     %% u_vars [theta; lnvm]
            z = obj.soln.z;     %% z_vars [zr; zi]
            mm.var.add('z_vars', 2*obj.state.N);  %% Create varset for z-vars            

            % Create quadratic constraints for complex power port injection
            half_s = obj.s / 2;

            % (1) constraints for u-vars
            if mm.userdata.tpc.quad        %% quadratic tpc-based formulation
                id_W = mat2cell((1:obj.port.N)', ones(obj.port.N,1));
                id_vars_u = (1:length(u))';
                nvars_u = length(u);
                Qu = mm.rearrange_quadratic_terms(obj.Qu, id_W, id_vars_u, nvars_u);
                mm.qcn.add(mm.var,'Port_inj_u', Qu, obj.M, -half_s, -half_s, {'u_vars'});
            else                             %% linear tpc-based formulation
                mm.lin.add(mm.var,'Port_inj_u', obj.M, -half_s, -half_s, {'u_vars'});
            end

            % (2) constraints for z-vars
            if mm.userdata.tpc.quad        %% quadratic tpc-based formulation
                id_vars_z = (1:length(z))';
                nvars_z = length(z);
                Qz = mm.rearrange_quadratic_terms(obj.Qz, id_W, id_vars_z, nvars_z);
                mm.qcn.add(mm.var, 'Port_inj_z', Qz, obj.N, -half_s, -half_s, {'z_vars'});
            else                            %% linear tpc-based formulation
                mm.lin.add(mm.var,'Port_inj_z', obj.N, -half_s, -half_s, {'z_vars'});
            end

            %% Compute complex power port injections
            % prepare full vector of vars with updated u-vars
            x = mm.var.params();
            i1u = mm.var.idx.i1.u_vars;
            iNu = mm.var.idx.iN.u_vars;
            i1z = mm.var.idx.i1.z_vars;
            iNz = mm.var.idx.iN.z_vars;
            x(i1u:iNu) = u;
            x(i1z:iNz) = z;
            if mm.userdata.tpc.quad        %% quadratic tpc-based formulation
                obj.soln.gs_ = mm.qcn.eval(mm.var, x, 'Port_inj_u')  + ...
                               mm.qcn.eval(mm.var, x, 'Port_inj_z');
            else                            %% linear tpc-based formulation
                obj.soln.gs_ = mm.lin.eval(mm.var, x, 'Port_inj_u') + ...
                               mm.lin.eval(mm.var, x, 'Port_inj_z');
            end
        end

        % -----------------------------------------------------------------
        % ************ OVERRIDDEN METHODS FROM PARENT CLASSES ************
        % -----------------------------------------------------------------

        function P = stack_matrix_params(obj, name, QnotM)
            % THIS FUNCTION OVERRIDES THE CORRESPONDING METHOD OF PARENT
            % CLASS mp.net_model (due to working with real vector x, and 
            % both cell and matrix parameters, i.e. Qu, Qz, M, and N)
            %          
            %
            % Form network parameter by stacking vertically the corresponding 
            % element parameters.
            % 
            %
            %   P = nm.stack_matrix_params(name, WnotM)
            %
            % Inputs:
            %   name (char array) : name of the parameter of interest            
            %   WnotM (boolean) : true if parameter to be stacked is of
            %                     class 'cell' (i.e. Qu or Qz), false otherwise            
            %
            % Outputs:
            %   P (cell/array) : parameter of interest for the full network            
            %
            % A given parameter (e.g. ``Qu``) for the full network is
            % formed by stacking vertically the corresponding parameters 
            % for each element.
            
            if QnotM
                P1p = {};
                P3p = {};
            else
                P1p = [];
                P3p = [];
            end

            for k = 1:length(obj.elements)
                nme = obj.elements{k};
                Pk = nme.(name);
                if isempty(Pk)
                    if QnotM                        
                        Pk = repmat([1 1 0], nme.nk * nme.np, 1);
                        Pk = mat2cell(Pk, ones(nme.nk * nme.np, 1));
                    else
                        if strcmp(name, 'M')
                            Pk = spalloc(nme.nk * nme.np, 2*obj.node.N, 0);
                        else
                            Pk = spalloc(nme.nk * nme.np, 2*obj.state.N, 0);
                        end                        
                    end
                end
                if ismember('3',nme.name) 
                    P3p = vertcat(P3p,Pk);
                else
                    P1p = vertcat(P1p,Pk);
                end                
            end

            if QnotM
                P = [P1p; P3p];
            else
                P = blkdiag(P1p, P3p);
            end
        end

        function display(obj)
            % Display the network model element object.
            %
            % This method is called automatically when omitting a semicolon
            % on a line that retuns an object of this class.
            %
            % Displays the details of the elements, including total number
            % of elements, nodes per element, ports per element, non-voltage
            % state per element, formulation name, tag, and class, and names
            % and dimensions of the model parameters.

            if obj.userdata.tpc.quad
                linorquad = 'quadratic';
            else
                linorquad = 'linear';
            end

            fprintf('NETWORK MODEL ELEMENT NAME  : %s\n', obj.name);
            fprintf('NETWORK MODEL ELEMENT CLASS : %s\n', class(obj));            
            fprintf('    # OF ELEMENTS           : %d\n', obj.nk);
            fprintf('    # OF NODES/ELEM         : %d\n', obj.nn);
            fprintf('    # OF PORTS/ELEM         : %d\n', obj.np);
            fprintf('    # OF NON-V STATES/ELEM  : %d\n', obj.nz);
            if isa(obj, 'mp.form')
                fprintf('    FORMULATION NAME        : %s\n', obj.form_name());
                fprintf('    FORMULATION TAG         : %s\n', obj.form_tag());
                fprintf('    FORMULATION TYPE        : %s\n', linorquad);
                fprintf('    FORMULATION CLASS       : %s\n', obj.find_form_class());                
                fprintf('    MODEL PARAMETERS');
                model_params = obj.model_params();
                for j = 1:length(model_params)
                    pn = model_params{j};   %% parameter name
                    if j == 1
                        fmt = '%7s : ';
                    else
                        fmt = '%27s : ';
                    end
                    if isa(obj.(pn), "cell")
                        nk = size(obj.(pn),1);
                        nempty = sum(cellfun('isempty', obj.(pn)));
                        if strcmp(pn, 'Qz')
                            if sum(sum(cell2mat(obj.Qz))) == 2*nk
                                nempty = nk;
                            end
                        end
                        if nempty == nk
                            fprintf([fmt '-\n'], pn);
                        else
                            s = 'non-empty';
                            fprintf([fmt '%d / %-5d%s\n'], pn, (nk-nempty), nk, s);
                        end
                    else
                        if isempty(obj.(pn))
                            fprintf([fmt '-\n'], pn);
                        else
                            [m, n] = size(obj.(pn));
                            if ~full(any(any(obj.(pn))))
                                s = '(all zeros)';
                            else
                                s = '';
                            end
                            fprintf([fmt '%d x %-7d%s\n'], pn, m, n, s);
                        end
                    end                    
                end               
            end

            %% nodes and states
            obj.display_set('node');
            obj.display_set('port');
            obj.display_set('state');

            %% variables
            vvars = obj.model_vvars();
            zvars = obj.model_zvars();
            for k = 1:length(vvars)
                obj.display_set(vvars{k});
            end
            for k = 1:length(zvars)
                obj.display_set(zvars{k});
            end

            %% elements
            model_params = obj.model_params();
            fprintf('ELEMENTS\n')
            fprintf('========\n')
            fprintf('  name                  N      np    nz    class, param(m,n))\n');
            fprintf(' ----------------   --------  ----  ----  --------------------\n');
            for k = 1:length(obj.elements)
                nme = obj.elements{k};
                fprintf('  %-13s %11d %5d %5d    %s', nme.name, nme.nk, nme.np, nme.nz, class(nme));

                for j = 1:length(model_params)
                    pn = model_params{j};   %% parameter name
                    if ~isempty(nme.(pn))
                        [m, n] = size(nme.(pn));
                        fprintf(', %s(%d,%d)', pn, m, n);
                    end
                end
                fprintf('\n');
            %     nme
            end

            %% user data
            fields = fieldnames(obj.userdata);
            if ~isempty(fields)
                fprintf('\nUSER DATA\n')
                fprintf('=========\n')
                fprintf('  name                         size       class\n');
                fprintf(' ------------------------   -----------  --------------------\n');
                for k = 1:length(fields)
                    f = obj.userdata.(fields{k});
                    [m, n] = size(f);
                    fprintf('  %-24s  %5dx%-5d   %s\n', fields{k}, m, n, class(f));
                end
            end

        end

        function CD = stack_CD_matrix(obj, name, rnotc)
        % Form network port/state matrices by stacking corresponding element
        % C/D matrices along the main matrix block diagonal. 
        % ::
        %
        %   M = nm.stack_CD_matrix(name, rnotc)
        %
        % Inputs:
        %   name (char array) : name of the parameter of interest ('C' or 'D')                    
        %   vnotz (boolean) : true if formulation works with real variables 
        %                     rather than complex variables
        %
        % Outputs:
        %   M (double) : matrix parameter of interest for the full network
        %   id_row (struct) : struct with fields taken from the element 
        %                     classes used to form the matrix parameter
        %                     of interest. Each field has a vector of the
        %                     form [row_start row_end], with the start/end
        %                     row indices of the matrix parameter being 
        %                     stacked along the matrix block diagonal.
        %   id_col (struct) : struct with fields taken from the element 
        %                     classes used to form the matrix parameter
        %                     of interest. Each field has a vector of the
        %                     form [col_start col_end], with the start/end
        %                     column indices of the matrix parameter being 
        %                     stacked along the matrix block diagonal
        %
        % A given port/state matrix parameter (e.g. ``C``) for the full 
        % network is formed by stacking the corresponding matrix parameters 
        % for each element along the matrix block diagonal.

        if rnotc
            ii = {};
            jj = {};
            ss = {};
            last_i = 0;
            last_j = 0;
            for k = 1:length(obj.elements)
                nme = obj.elements{k};
                Mk = nme.(name);
                Mk = blkdiag(Mk,Mk);
                if ~isempty(Mk)
                    [i, j, s] = find(Mk);
                    ii = horzcat(ii, i + last_i);
                    jj = horzcat(jj, j + last_j);
                    ss = horzcat(ss, s);
                    m = 2 * nme.nk * nme.np;        %% total number of ports for class
                    if strcmp(name,'C')
                        n = m;
                    else
                        n = 2 * nme.nk * nme.nz;    %% total number of states for class
                    end
                else
                    m = nme.nk * nme.np;
                    n = m;
                end
                
                last_i = last_i + m;
                last_j = last_j + n;
            end

            CD = sparse(vertcat(ii{:}), vertcat(jj{:}), vertcat(ss{:}), last_i, last_j);
        else
            if strcmp(name,'C')
                vnotz = 1;
            else
                vnotz = 0;
            end
            CD = obj.stack_matrix_params(name, 0, vnotz);
        end
        end         
    end     %% methods
end         %% classdef
