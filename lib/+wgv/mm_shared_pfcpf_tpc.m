classdef (Abstract) mm_shared_pfcpf_tpc < mp.mm_shared_pfcpf
% wgv.mm_shared_pfcpf_tpc - Mixin class for TPC power flow (PF) 
% **math model** objects.
%
% An abstract mixin class inherited by TPC power flow (PF) **math model**
% objects.

%   MATPOWER
%   Copyright (c) 2021-2025, Power Systems Engineering Research Center (PSERC)
%   by Wilson Gonzalez Vanegas, Universidad Nacional de Colombia Sede Manizales
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://matpower.org for more info.

%     properties
%     end

    methods
        function ad = build_aux_data(obj, nm, dm, mpopt)
            %

            %% call parent
            ad = build_aux_data@mp.mm_shared_pfcpf(obj, nm, dm, mpopt);

            %% Calculate Sbus 
            if nm.elements.has_name('gen')    %% (1-phase buses)
                Cgen = nm.elements.gen.C;
                Dgen = nm.elements.gen.D;
                if nm.elements.has_name('load')
                    Cload = nm.elements.load.C;
                    Sload = nm.elements.load.s;
                else
                    Cload = 0;
                    Sload = 0;
                end
                
                Sbus = Cgen*Dgen'*(ad.zr + 1j * ad.zi) - Cload*Sload;
            else
                Sbus = [];
                Sload = [];
            end            

            ad.Sbus = Sbus;
            ad.Sload = Sload;

            %% add fields .lin and .quad to the property 'userdata' of
            %  the math model according to this same properties of 
            %  the task object

            if isfield(nm.userdata, 'tpc')
                obj.userdata.tpc.lin  = nm.userdata.tpc.lin;
                obj.userdata.tpc.quad = nm.userdata.tpc.quad;
            end            
        end

        function obj = add_system_vars_pf(obj, nm, dm, mpopt)
            %

            %% get model variables
            vvars = nm.model_vvars();

            %% voltage angles
            obj.add_system_varset_pf(nm, vvars{1}, 'pv');
            obj.add_system_varset_pf(nm, vvars{1}, 'pq');

            %% voltage magnitudes
            obj.add_system_varset_pf(nm, vvars{2}, 'pq');
        end

        function obj = add_system_varset_pf(obj, nm, vvar, typ)
            %
            ad = obj.aux_data;
            st = nm.(vvar);
            d = st.data;
            mmx_i1 = obj.var.N + 1;
            for k = 1:st.NS
                name = st.order(k).name;
                idx = st.order(k).idx;
                ii = ad.node_type_by_elm(k).(typ);
                nii = length(ii);
                if isempty(idx)
                    obj.add_var([name '_' typ], nii, d.v0.(name)(ii), d.vl.(name)(ii), d.vu.(name)(ii));
                else
                    if all(cell2mat(idx) == 1)
                        dim = size(st.idx.N.(name));
                        if dim(end) == 1, dim(end) = []; end    %% delete trailing 1
                        obj.init_indexed_name('var', [name '_' typ], num2cell(dim));
                    end
                    sc = struct('type', {'{}', '()'}, 'subs', {idx, {ii}});
                    v0 = subsref(d.v0.(name), sc);
                    vl = subsref(d.vl.(name), sc);
                    vu = subsref(d.vu.(name), sc);
                    obj.add_var([name '_' typ], idx, nii, v0, vl, vu);
                end
            end
            mmx_iN = obj.var.N;
            if ad.(['n' typ])
                obj.aux_data.var_map{end+1} = ...
                    {vvar, [], [], ad.(typ), mmx_i1, mmx_iN, []};
            end
        end

        function [u, z, x] = convert_x_m2n(obj, mmx, nm, only_u)
            % convert_x_m2n - Convert math model state to network model state.
            % ::
            %
            %   x = mm.pf_convert(mmx, nm)
            %   [v, z] = mm.pf_convert(mmx, nm)
            %   [v, z, x] = mm.pf_convert(mmx, nm)
            %   ... = mm.pf_convert(mmx, nm, only_v)

            %% update v_, z_ from mmx
            nm_vars = obj.update_nm_vars(mmx, nm);
            u = [nm_vars.va; nm_vars.lnvm];
            z  = [nm_vars.zr; nm_vars.zi];

            %% update z, if requested
            if nargin < 4 || ~only_u
                z = obj.update_z(nm, u, z);
            end

            %% prepare return values
            if nargout < 2
                u = [u; z];
            elseif nargout > 2
                x = [u; z];
            end
        end

        function z = update_z(obj, nm, u, z)
            % update_z - Update/allocate active/reactive injections at slack/PV nodes.
            %
            % Update/allocate slack know active power injections and slack/PV
            % node reactive power injections.

            ad = obj.aux_data;
            if nm.elements.has_name('load')
                Cload = nm.elements.load.C;
            else
                Cload = [];
            end
            Sload = Cload * ad.Sload;
            rpv = [ad.ref; ad.pv];
            id_vars = (1 : 2*nm.node.N)';    %% full vector of u-vars [\theta; lnvm]

            % Find ports of branches and shunts connected to each bus (for computing bus injection)
            ports = obj.ports_by_bus(nm);

            % Define quadratic forms for complex power bus injections at REF-PV buses
            [Qbusrpv, Cbus_rpv, kbus_rpv] = obj.bus_complex_injection(nm, ports(rpv), id_vars);

            % coefficient matrix for power injection states for full network
            CC = nm.C * nm.get_params([], 'N') * blkdiag(nm.D, nm.D)';

            %% Add constraints for Pmis at REF and Qmis at REF/PV buses
            % Linear part goes first:
            
            % Active power at REF nodes
            C_Pref = real(Cbus_rpv(1:ad.nref,:));
            k_Pref = real(kbus_rpv(1:ad.nref));            

            % Reactive power at PV/REF nodes
            C_Qrpv = imag(Cbus_rpv);
            k_Qrpv = imag(kbus_rpv);   

            
            % Quadratic part, if present, goes afterwards:
            if obj.userdata.tpc.quad
                Q_Pref = cellfun(@(x)(real(x)), Qbusrpv(1:ad.nref), 'UniformOutput', false);   % Active power at REF nodes
                Q_Qrpv = cellfun(@(x)(imag(x)), Qbusrpv, 'UniformOutput', false);              % Reactive power at PV/REF nodes
            else
                Q_Pref = [];
                Q_Qrpv = [];
            end
    
            % Add set of variables            
            obj.var.add('u_vars', 2*nm.node.N, u); % u_vars = [theta; lnvm]

            % Add set of constraints
            if obj.userdata.tpc.quad         %% quadratic tpc-based formulation                
                obj.qcn.add(obj.var,'Pmis_ref', Q_Pref, C_Pref, -k_Pref, -k_Pref, {'u_vars'});
                obj.qcn.add(obj.var,'Qmis_rpv', Q_Qrpv, C_Qrpv, -k_Qrpv, -k_Qrpv, {'u_vars'});
            else                              %% linear tpc-based formulation                
                obj.lin.add(obj.var,'Pmis_ref', C_Pref, -k_Pref, -k_Pref, {'u_vars'});
                obj.lin.add(obj.var,'Qmis_rpv', C_Qrpv, -k_Qrpv, -k_Qrpv, {'u_vars'});
            end

            %% ----- Update active power at slack nodes -----
            % coefficient matrix for active power injection states for slack nodes
            CCref = CC(ad.ref, :);
            jr = find(any(CCref(:,1:size(nm.D,1)), 1));   %% indices of corresponding states            

            % allocate active power at slack nodes to 1st direct inj state:
            % find all z (except first one) with direct injection at each
            % slack node
            [i, j] = find(CCref(:,1:size(nm.D,1)));
            if size(i, 2) > 1, i = i'; j = j'; end
            ij = sortrows([i j]);       %% 1st state comes 1st for each node
            [~, k1] = unique(ij(:, 1), 'first');%% index of 1st entry for each node
            % all included states that are not 1st at their node
            jn = unique(ij(~ismember(1:length(i), k1), 2));
            jfirst = setdiff(jr, jn);
            
            if ~isempty(jn) % if we have extra states (more than 1) for any node(s)
                [genbus, ~] = find(nm.elements.gen.C);
                Pref_jn = sparse(genbus(jn), ones(length(jn),1), z(jn), length(genbus), 1);   %% total generation except first gen 
                Pref_jn = full(Pref_jn(Pref_jn ~= 0));
            else
                Pref_jn = 0;
            end

            %% compute port active power injections at REF buses except generators
            % prepare full vector of vars with updated u-vars
            x = obj.var.params();
            i1 = obj.var.idx.i1.u_vars;
            iN = obj.var.idx.iN.u_vars;
            x(i1:iN) = u;
            if obj.userdata.tpc.quad        %% quadratic tpc-based formulation
                Pbus_ref = obj.qcn.eval(obj.var, x, 'Pmis_ref') + real(Sload(ad.ref));
            else                             %% linear tpc-based formulation
                Pbus_ref = obj.lin.eval(obj.var, x, 'Pmis_ref') + real(Sload(ad.ref));
            end
            
            Pgen_ref = CCref(:,jfirst) \ - Pbus_ref;
            
            % allocate active power to the first gen at each REF bus            
            z(jfirst) = Pgen_ref - Pref_jn;

            %% ----- Update reactive power at PV/Slack nodes -----
            % coefficient matrix for reactive power injection states for PV/REF nodes
            CCrpv = CC(rpv, :);
            jrpv = find(any(CCrpv(:,1:size(nm.D,1)), 1));   %% indices of corresponding states            

            % compute port reactive power injections at PV/REF buses except generators
            if obj.userdata.tpc.quad        %% quadratic tpc-based formulation
                Qbusrpv = obj.qcn.eval(obj.var, x, 'Qmis_rpv') + imag(Sload(rpv));
            else                             %% linear tpc-based formulation
                Qbusrpv = obj.lin.eval(obj.var, x, 'Qmis_rpv') + imag(Sload(rpv));
            end

            % allocate reactive power according to the number of gens at each node
            if nm.elements.has_name('gen')
                Cgen = nm.elements.gen.C;
                ngen = nm.elements.gen.nk;
            else
                Cgen = [];
                ngen = 0;
            end           
            Qgen = [];            
            ngens = full(sum(Cgen, 2)); ngens(ad.pq) = [];
            for g = 1:length(ngens)
                Qgen = [Qgen; repmat(Qbusrpv(g), ngens(g), 1)];
            end

            % find indices of z-vars related to reactive injections for gens 
            % according to the rpv indices
            if nm.elements.has_name('gen')
                Dgen = nm.elements.gen.D;
            else
                Dgen = [];
            end
            CD = Cgen * Dgen;
            [Q_id, ~] = find(CD(rpv,:)'); 
            
            % compute weights to distribute reactive injections in gens
            % located at the same bus. In this implementation, reactive
            % power is allocated proportionally to gen active power
            weights = [];
            for b = rpv'
                idg = logical(full(Cgen(b,:)));
                Pbus = sum(abs(z(idg)));
                weights = [weights; abs(z(idg)./Pbus)];
            end
            weights(isnan(weights)) = 1;
            z(ngen + Q_id) = weights.*Qgen;
        end
    end     %% methods
end         %% classdef





















