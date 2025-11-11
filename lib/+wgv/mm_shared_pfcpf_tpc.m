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
                        
            if nm.userdata.ishybrid
                %% Add indices of non-voltage variables for buslinks
                nz_bl = sum(nm.state.idx.N.buslink);
                if nm.elements.has_name('gen')
                    nz1p = nm.state.idx.N.gen;
                else
                    nz1p = 0;
                end
                nz3p = nz_bl;
                if nm.elements.has_name('gen3p')
                    nz3p = nz3p + sum(nm.state.idx.N.gen3p);
                end
                nz = nz1p + nz3p;
                ad.idzbl = (1:2*nz_bl)'+2*nz1p;

                %% Compute [pv pq; pq + nb#p], with # holding single-phase and three-phase indices
                nb1p = nm.node.idx.N.bus;
                nb3p = nm.node.idx.N.bus3p(1);
                nb = nb1p + 3*nb3p;

                pv1p = dm.elements.bus.on(dm.elements.bus.tab.type == 2);
                npv1p = numel(pv1p);
                pq1p = dm.elements.bus.on(dm.elements.bus.tab.type == 1);
                npq1p = numel(pq1p);
                pv3p = dm.elements.bus3p.on(dm.elements.bus3p.tab.type == 2);
                npv3p = numel(pv3p);
                pv3p = repmat(pv3p,1,3) + [zeros(npv3p,1) ones(npv3p,1)*nb3p ones(npv3p,1)*2*nb3p] + 2*nb1p;
                pq3p = dm.elements.bus3p.on(dm.elements.bus3p.tab.type == 1);
                npq3p = numel(pq3p);
                pq3p = repmat(pq3p,1,3) + [zeros(npq3p,1) ones(npq3p,1)*nb3p ones(npq3p,1)*2*nb3p] + 2*nb1p;                

                vars_1p_3p_u = [[pv1p; pq1p; pq1p+nb1p];              % single-phase
                              [pv3p(:); pq3p(:); pq3p(:)+3*nb3p]];  % three-phase

                %% Compute permutation matrix for mapping pvq variables from 1p/3p order to system order
                npv3p = numel(pv3p);
                npq3p = numel(pq3p);
                npv = npv1p+npv3p;
                npq = npq1p+npq3p;
                pm_phase_pvq_to_syst_u = [
                    speye(npv1p)   spalloc(npv1p,2*npq+npv3p,0)
                    spalloc(npv3p,npv1p+2*npq1p,0) speye(npv3p) spalloc(npv3p,2*npq3p,0)
                    spalloc(npq1p,npv1p,0) speye(npq1p) spalloc(npq1p,npq+npv3p+npq3p,0)
                    spalloc(npq3p,npv+2*npq1p,0) speye(npq3p) spalloc(npq3p,npq3p,0)
                    spalloc(npq1p,npv1p+npq1p,0) speye(npq1p) spalloc(npq1p,npv3p+2*npq3p,0)
                    spalloc(npq3p,npv+npq+npq1p,0) speye(npq3p)];

                %% Compute permutation matrix for mapping buslink vars from 1p/3p order to system order
                if nm.elements.has_name('gen')
                    ngen1p = nm.elements.gen.nk;
                else
                    ngen1p = 0;
                end
                if nm.elements.has_name('gen3p')
                    ngen3p = nm.elements.gen3p.nk;
                else
                    ngen3p = 0;
                end
                pm_phase_to_syst_z = [
                    speye(ngen1p) spalloc(ngen1p,ngen1p+2*ngen3p+2*nz_bl,0)
                    spalloc(ngen3p,2*ngen1p,0) speye(ngen3p) spalloc(ngen3p,ngen3p+2*nz_bl,0)
                    spalloc(nz_bl,2*(ngen1p+ngen3p),0) speye(nz_bl) spalloc(nz_bl,nz_bl,0)];
                pm_phase_to_syst_z = repmat(pm_phase_to_syst_z,2,1);

                %% Compute permutation matrix for mapping all voltage angles and magnitudes from system order to 1p/3p
                nb3p = 3*nb3p;
                pm_all_va_lnvm_u = [
                    speye(nb1p) spalloc(nb1p,2*nb3p+nb1p,0)
                    spalloc(nb1p,nb,0) speye(nb1p) spalloc(nb1p,nb3p,0)
                    spalloc(nb3p,nb1p,0) speye(nb3p) spalloc(nb3p,nb,0)
                    spalloc(nb3p,nb+nb1p,0) speye(nb3p)];               
                
                % Compute permutation matrix for mapping all states from real/imag order to 1p/3p order
                pm_all_1p3p_z = [
                    speye(nz1p) spalloc(nz1p,nz+nz3p,0)
                    spalloc(nz3p,2*nz1p,0) speye(nz3p) spalloc(nz3p,nz3p,0)
                    spalloc(nz1p,nz1p,0) speye(nz1p) spalloc(nz1p,2*nz3p,0)                                        
                    spalloc(nz3p,nz+nz1p,0) speye(nz3p)];

                vars_1p_3p_z = (1:2*nz_bl)'+2*nb;                

                ad.pv1p = pv1p;
                ad.pv3p = pv3p;
                ad.vars_1p_3p_pvq = [vars_1p_3p_u; vars_1p_3p_z];
                ad.pm_all_va_lnvm_u = pm_all_va_lnvm_u;
                ad.pm_phase_to_syst_pvq = blkdiag(pm_phase_pvq_to_syst_u, speye(2*nz_bl));
                ad.pm_phase_to_syst_all = blkdiag(pm_all_va_lnvm_u, pm_all_1p3p_z);
                ad.pm_phase_to_syst_z = pm_phase_to_syst_z;
                ad.pm_all_1p3p_z = pm_all_1p3p_z;
            end            

            %% Calculate Sbus (1-phase buses)
            if nm.elements.has_name('bus')
                nb1p = nm.node.idx.N.bus;
            end

            if nm.elements.has_name('gen')    %% (1-phase buses)
                if nm.userdata.ishybrid                    
                    Cgen1p = nm.elements.gen.C(1:nb1p,:);                    
                else
                    Cgen1p = nm.elements.gen.C;
                end
                Dgen1p = nm.elements.gen.D;
            else
                Cgen1p = [];
                Dgen1p = [];
            end

            if nm.elements.has_name('load')                
                if nm.userdata.ishybrid
                    Cload1p = nm.elements.load.C(1:nb1p,:);
                else
                    Cload1p = nm.elements.load.C;
                end
                Sload1p = nm.elements.load.s;
            else
                Cload1p = [];
                Sload1p = [];
            end
            
            if ~isempty(Cgen1p)
                if ~isempty(Cload1p)
                    Sbus1p = Cgen1p*Dgen1p'*(ad.zr + 1j * ad.zi) - Cload1p*Sload1p;
                else
                    Sbus1p = Cgen1p*Dgen1p'*(ad.zr + 1j * ad.zi);
                end                
            else
                if isempty(Cload1p)
                    if nm.elements.has_name('bus')
                        Sbus1p = zeros(nb1p,1);
                    else
                        Sbus1p = [];
                    end                    
                else
                    Sbus1p = - Cload1p*Sload1p;
                end                
            end
            

            %% Calculate Sbus (3-phase buses)
            if nm.elements.has_name('gen3p')    %% (1-phase buses)
                Cgen3p = nm.elements.gen3p.C;                
                Dgen3p = nm.elements.gen3p.D;
            else
                Cgen3p = [];
                Dgen3p = [];
            end

            if nm.elements.has_name('load3p')
                Cload3p = nm.elements.load3p.C;
                Sload3p = nm.elements.load3p.s;
            else
                Cload3p = [];
                Sload3p = [];
            end            
            
            if ~isempty(Cgen3p)
                if ~isempty(Cload3p)
                    Sbus3p = Cgen3p*Dgen3p'*(ad.zr + 1j * ad.zi) - Cload3p*Sload3p;
                else
                    Sbus3p = Cgen3p*Dgen3p'*(ad.zr + 1j * ad.zi);
                end                
            else
                if isempty(Cload3p)
                    if nm.elements.has_name('bus3p')
                        Sbus3p = zeros(nm.node.idx.N.bus3p,1);
                    else
                        Sbus3p = [];
                    end                    
                else
                    Sbus3p = - Cload3p*Sload3p;
                end                
            end          

            %%
            if nm.userdata.ishybrid
                ad.Sbus = [Sbus1p; Sbus3p(nb1p+1:end)];
            else
                ad.Sbus = [Sbus1p; Sbus3p];
            end
            
            ad.Sload = [Sload1p; Sload3p];
            

            %% Add fields .lin and .quad to the property 'userdata' of
            %  the math model according to this same properties of 
            %  the task object
            
            if isfield(nm.userdata, 'tpc')
                obj.userdata.tpc.lin  = nm.userdata.tpc.lin;
                obj.userdata.tpc.quad = nm.userdata.tpc.quad;
            end

            % Use initial value of variables at this point
            ad.var_source = 'u0';

            % Update aux data struct
            obj.aux_data = ad;
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

            %% 
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
                Cload1p = nm.elements.load.C;
            else
                Cload1p = [];
            end
            if nm.elements.has_name('load3p')
                Cload3p = nm.elements.load3p.C;
            else
                Cload3p = [];
            end
            Cload = blkdiag(Cload1p, Cload3p);
            Sload = Cload * ad.Sload;
            rpv = [ad.ref; ad.pv];

            if nm.userdata.ishybrid                
                id_vars = ad.vars_1p_3p_pvq;
            else                
                pvq = [ad.pv; ad.pq];
                id_vars = [pvq; ad.pq + nm.node.N];    % the vector of variables is [theta_pv; theta_pq; lnvm_pq]  (npv+2np variables)
            end            

            % Find ports of branches and shunts connected to each bus (for computing bus injection)
            ports = obj.ports_by_bus(nm);

            % Define quadratic forms for complex power bus injections at REF-PV buses
            ad.var_source = 'solx';
            [Qbusrpv, Cbus_rpv, kbus_rpv] = obj.bus_complex_injection(nm, ports(rpv), id_vars, ad);

            % coefficient matrix for power injection states for full network buses
            if nm.userdata.ishybrid
                CC = nm.C * nm.N * ad.pm_all_1p3p_z';
            else
                CC = nm.C * nm.N;
            end
            

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
            id_vars = [ad.pv; ad.pq; ad.pq + nm.node.N];
            if nm.userdata.ishybrid
                i1 = nm.state.idx.i1.buslink(1);
                iN = nm.state.idx.iN.buslink(end);
                nb = nm.node.N;
                nz = nm.state.N;
                id_vars = [id_vars; (i1:iN)' + 2*nb; (i1:iN)' + 2*nb+nz];
                vars = [u; z];
            else
                vars = u;
            end
            mm_copy = obj.copy();
            mm_copy.var.add('vars', numel(id_vars), vars(id_vars));

            % Add set of constraints
            if obj.userdata.tpc.quad         %% quadratic lnvm-based formulation                
                mm_copy.qcn.add(mm_copy.var,'Pmis_ref', Q_Pref, C_Pref, -k_Pref, -k_Pref, {'vars'});
                mm_copy.qcn.add(mm_copy.var,'Qmis_rpv', Q_Qrpv, C_Qrpv, -k_Qrpv, -k_Qrpv, {'vars'});
            else                              %% linear lnvm-based formulation                
                mm_copy.lin.add(mm_copy.var,'Pmis_ref', C_Pref, -k_Pref, -k_Pref, {'vars'});
                mm_copy.lin.add(mm_copy.var,'Qmis_rpv', C_Qrpv, -k_Qrpv, -k_Qrpv, {'vars'});
            end

            %% Locate and classify generators at PV and Ref nodes
            gen_info = obj.locate_generators(nm,CC);

            %% ----- Update active power at slack nodes -----
            % coefficient matrix for active power injection states for slack nodes
            CCref = CC(ad.ref, :);            
            
            if obj.userdata.tpc.quad        %% quadratic lnvm-based formulation
                Pbus_ref = mm_copy.qcn.eval(mm_copy.var, mm_copy.var.params(), 'Pmis_ref') + real(Sload(ad.ref));
            else                             %% linear lnvm-based formulation
                Pbus_ref = mm_copy.lin.eval(mm_copy.var, mm_copy.var.params(), 'Pmis_ref') + real(Sload(ad.ref));
            end
            
            % allocate active power to generators connected to Ref nodes 
            idz = gen_info.Ref.idz_single;
            if ~isempty(gen_info.Ref.idz_multiple)
                all1st = cell2mat(cellfun(@(x)(x(1)),gen_info.Ref.idz_multiple,'UniformOutput',false));
                idz = [idz; all1st];
            end
            Ptotal_ref = real(CCref(:,idz)) \ -Pbus_ref;  % If there are multiple gens connected at the same bus, the active power of the first one is computed
            
            % multiple generators
            if ~isempty(gen_info.Ref.idz_multiple)
                % distribute active power along all gens acording to Pmax
                [~,id_multiple_ref] = ismember(all1st,idz);
                Pbus = Ptotal_ref(id_multiple_ref);
                for b = 1:numel(Pbus) 
                    idgens = gen_info.Ref.idz_multiple{b};
                    Pmax = gen_info.Ref.Pmax_multiple{b};
                    Pmax(isinf(Pmax)) = 0;
                    if any(Pmax)
                        weights = Pmax/sum(Pmax);
                    else
                        weights = 1/length(Pmax);
                    end
                    z(idgens) = Pbus(b)*weights;
                end
            end

            % single generators            
            if ~isempty(gen_info.Ref.idz_single)
                [~,idgens] = ismember(gen_info.Ref.idz_single,idz);
                z(gen_info.Ref.idz_single) = Ptotal_ref(idgens);
            end

            %% ----- Update reactive power at PV/Slack nodes -----
            nz = nm.state.N;

            % coefficient matrix for reactive power injection states for PV/REF nodes
            CCrpv = CC(rpv, :);

            % extract active power injection of multiple generators connected at PV/Ref nodes
            idz_multiple_rpv = [gen_info.Ref.idz_multiple; gen_info.PV.idz_multiple];
            if ~isempty(idz_multiple_rpv)
                Pg_multiple_rpv = cell(numel(idz_multiple_rpv),1);
                for i = 1:numel(idz_multiple_rpv)
                    Pg_multiple_rpv{i} = z(idz_multiple_rpv{i});
                end
            end

            % compute port reactive power injections at PV/REF buses except generators
            if obj.userdata.tpc.quad        %% quadratic lnvm-based formulation
                Qbus_rpv = mm_copy.qcn.eval(mm_copy.var, mm_copy.var.params(), 'Qmis_rpv') + imag(Sload(rpv));
            else                             %% linear lnvm-based formulation
                Qbus_rpv = mm_copy.lin.eval(mm_copy.var, mm_copy.var.params(), 'Qmis_rpv') + imag(Sload(rpv));
            end

            % allocate reactive power to generators connected to PV/Ref nodes 
            idz = [gen_info.Ref.idz_single; gen_info.PV.idz_single];
            if ~isempty(gen_info.Ref.idz_multiple)
                all1st = cellfun(@(x)(x(1)),gen_info.Ref.idz_multiple,'UniformOutput',false);
            else
                all1st = [];
            end
            if ~isempty(gen_info.PV.idz_multiple)
                all1st_2 = cellfun(@(x)(x(1)),gen_info.PV.idz_multiple,'UniformOutput',false);
                all1st = cell2mat([all1st; all1st_2]);
            end
            idz = [idz; all1st];
            Qtotal_rpv = imag(CCrpv(:,idz+nz)) \ -Qbus_rpv; % If there are multiple gens connected at the same bus, the reactive power of the first one is computed
            
            % multiple generators
            if ~isempty(idz_multiple_rpv)
                % distribute reactive power along all gens acording to injected active power ammounts
                [~,is_bus_multiple_rpv] = ismember(all1st,idz);
                Qbus = Qtotal_rpv(is_bus_multiple_rpv);
                for b = 1:numel(Qbus)
                    idgens = idz_multiple_rpv{b};
                    Pinj = Pg_multiple_rpv{b};
                    weights = Pinj/sum(Pinj);
                    z(idgens+nz) = Qbus(b)*weights;
                end                
            end

            % single generators
            id_single_rpv = [gen_info.Ref.idz_single; gen_info.PV.idz_single];
            if ~isempty(id_single_rpv)
                [~,idgens] = ismember(id_single_rpv,idz);                
                z(id_single_rpv+nz) = Qtotal_rpv(idgens);
            end

            %% Store copy of math model with the additioned variables and constraints
            obj.aux_data.mm_copy = mm_copy;
        end

        function classified_gens = locate_generators(obj,nm,CC)
            %
            
            % Extract info
            ad = obj.aux_data;
            nz = nm.state.N;
            if obj.elements.has_name('gen')
                ng1p = nm.elements.gen.nk;
            else
                ng1p = 0;
            end
            if obj.elements.has_name('gen3p')
                ng3p = nm.elements.gen3p.nk;
            else
                ng3p = 0;
            end
            
            % Initialize returned struct
            es = struct('idz_single', [], ...
                        'idz_multiple', [], ...
                        'is_bus_multiple', []);
            cg = struct('Ref',es,'PV',es);
            cg.Ref.Pmax_multiple = [];
            cg.Ref.Qmax_multiple = [];            

            %% 1) Locate and classify states for generators connected at reference nodes
            CCref = CC(ad.ref,:); %% coefficient matrix

            % Ref nodes with more that one generator connected to them
            cg.Ref.is_bus_multiple = zeros(ad.nref,1);
            id_ref_gt1 = find(sum(real(CCref),2) < -1);
            if ~isempty(id_ref_gt1)
                ngt1 = numel(id_ref_gt1);                
                cg.Ref.idz_multiple = cell(ngt1,1);
                cg.Ref.Pmax_multiple = cell(ngt1,1);
                cg.Ref.Qmax_multiple = cell(ngt1,1);
                for i = 1:ngt1
                    idbus = id_ref_gt1(i);
                    idz = find(CCref(idbus,1:nz));
                    if any(idz <= ng1p)
                        cg.Ref.Pmax_multiple{i} =  nm.zr.data.vu.Pg(idz);
                        cg.Ref.Qmax_multiple{i} =  nm.zi.data.vu.Qg(idz);
                    end
                    if any(idz > ng1p) && ng3p ~= 0
                        Pmax = zeros(size(idz)); Qmax = Pmax;
                        for p = 1:3
                            idlim = (idz >= (p-1)*ng3p+1) & (idz <= p*ng3p); 
                            Pmax(idlim) =  nm.zr.data.vu.Pg3{p}(idlim);
                            Qmax(idlim) =  nm.zi.data.vu.Qg3{p}(idlim);
                        end
                        cg.Ref.Pmax_multiple{i} = Pmax;
                        cg.Ref.Qmax_multiple{i} = Qmax;
                    end                    
                    cg.Ref.idz_multiple{i} = idz;                    
                end
                cg.Ref.is_bus_multiple(id_ref_gt1) = 1;                
            end
            cg.Ref.is_bus_multiple = logical(cg.Ref.is_bus_multiple);

            % Ref nodes with a single generator connected to them
            id_ref_single = setdiff((1:ad.nref)',id_ref_gt1);
            if ~isempty(id_ref_single)
                n1 = numel(id_ref_single);
                cg.Ref.idz_single = zeros(n1,1);
                for i = 1:n1
                    idbus = id_ref_single(i);
                    cg.Ref.idz_single(i) = find(CCref(idbus,1:nz));
                end
            end            

            %% 2) Locate and classify states for generators connected at PV nodes
            CCpv = CC(ad.pv,:); %% coefficient matrix

            % PV nodes with more that one generator connected to them
            cg.PV.is_bus_multiple = zeros(ad.npv,1);
            id_pv_gt1 = find(sum(real(CCpv),2) < -1);
            if ~isempty(id_pv_gt1)
                ngt1 = numel(id_pv_gt1);
                cg.PV.idz_multiple = cell(ngt1,1);
                for i = 1:ngt1
                    idbus = id_pv_gt1(i);
                    cg.PV.idz_multiple{i} = find(CCpv(idbus,1:nz));
                end
                cg.PV.is_bus_multiple(id_pv_gt1) = 1;                
            end
            cg.PV.is_bus_multiple = logical(cg.PV.is_bus_multiple);

            % PV nodes with a single generator connected to them
            id_pv_single = setdiff((1:ad.npv)',id_pv_gt1);
            if ~isempty(id_pv_single)
                n1 = numel(id_pv_single);
                cg.PV.idz_single = zeros(n1,1);
                for i = 1:n1
                    idbus = id_pv_single(i);
                    cg.PV.idz_single(i) = find(CCpv(idbus,1:nz));
                end
            end
          
            %% Create output struct
            classified_gens = cg;
        end
    end     %% methods
end         %% classdef