classdef mme_tdinter_orpd_tpc < mp.mm_element
% wgv.mme_tdinter_orpd_tpc - Math model element for 1-to-3-phase interface for 
%                            TPC-based Optimal Reactive Power Dispacht (ORPD).
%
% Math model element class for 1-to-3-phase interface elements for ORPD task
% under a TPC formulation.
%
% Implements methods for ...

%   MATPOWER
%   Copyright (c) 2021-2025, Power Systems Engineering Research Center (PSERC)
%   by Wilson Gonzalez Vanegas, Universidad Nacional de Colombia Sede Manizales
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://matpower.org for more info.

%     properties
%     end     %% properties

    methods
        function name = name(obj)
            %
            name = 'buslink';
        end        

        function obj = add_constraints(obj, mm, nm, dm, mpopt)
            %

            nme = obj.network_model_element(nm);
            dme = obj.data_model_element(dm);

            %% 1) Voltage constraints

            % variables sets
            vs_va = struct('name', {'Va', 'Va3', 'Va3', 'Va3'}, ...
                            'idx', {{}, {1}, {2}, {3}});
            vs_vm = struct('name', {'LnVm', 'LnVm3', 'LnVm3', 'LnVm3'}, ...
                            'idx', {{}, {1}, {2}, {3}});

            % define linear parameters for voltage constraints
            [A_v, b_v] = obj.voltage_constraints_tdinter(dme, dm.elements.bus.tab.uid, dm.elements.bus3p.tab.uid);

            % add linear constraints for voltage angles and magnitudes
            mm.lin.add(mm.var,'tdinter_va', A_v, b_v, b_v, vs_va);
            mm.lin.add(mm.var,'tdinter_lnvm', A_v, b_v, b_v, vs_vm);

            %% 2) Enchange power factor constraints

            % define linear constraints for per-phase exchange power factor limits
            [A_lead, A_lagg] = obj.exchange_pf_constraints_tdinter(dme, nme);            
            l_lead = -Inf(dme.n,1);
            u_lead = zeros(dme.n,1);
            l_lagg = zeros(dme.n,1);
            u_lagg = Inf(dme.n,1);

            % add quadratic constraints for power factor for leading exchange
            if ~isempty(A_lead)
                for p = 1:(nme.np-1)
                    vs = struct('name',{'Pinter', 'Qinter'}, ...
                                        'idx',{{p}, {p}});

                    mm.lin.add(mm.var,sprintf('min_pf_inter_lead_%d',p), A_lead, l_lead, u_lead, vs);
                end
            end

            if ~isempty(A_lagg)
                for p = 1:(nme.np-1)
                    vs = struct('name',{'Pinter', 'Qinter'}, ...
                                        'idx',{{p}, {p}});

                    mm.lin.add(mm.var,sprintf('min_pf_inter_lagg_%d',p), A_lagg, l_lagg, u_lagg, vs);
                end
            end

            %% 3) Capacity constraints

            % define quadratic forms for three-phase line capacity
            Q_Smax = obj.capacity_constraints_tdinter(dme, nme);
            B_max = [];
            l_max = dme.smax(dme.smax ~= 0);
            
            % add quadratic constraints for apparent power limits of three-phase lines
            if ~isempty(Q_Smax)
                for p = 1:(nme.np-1)
                    vs = struct('name',{'Pinter', 'Qinter'}, ...
                                        'idx',{{p}, {p}});
                    
                    mm.qcn.add(mm.var,sprintf('Smax_tdinter_%d',p),Q_Smax,B_max,-l_max,l_max,vs);                    
                end
            end
        end

        function [A, b] = voltage_constraints_tdinter(obj, dme, id_bus1p, id_bus3p)
            %
            
            nb1p = length(id_bus1p);
            nb3p = length(id_bus3p);

            [~, bus1p] = ismember(dme.tab.bus, id_bus1p);
            [~, bus3p] = ismember(dme.tab.bus3p, id_bus3p);

            jj1p = repmat(bus1p',3,1); 
            jj3p = repmat(bus3p',3,1) + repmat(nb1p + [0 nb3p 2*nb3p]', 1, length(bus3p));

            jj = [jj1p(:); jj3p(:)];
            ii = repmat((1:3*length(bus1p))',2,1);
            
            vals = [ones(3*dme.n,1); -1*ones(3*dme.n,1)];
            
            nvars = nb1p + 3*nb3p;
            A = sparse(ii, jj, vals, 3*dme.n, nvars);

            b = zeros(3*dme.n,1);
        end

        function [A_lead, A_lagg] = exchange_pf_constraints_tdinter(obj, dme, nme)
            %

            % Leading power factor constraints
            ninter = dme.n;
            if any(dme.min_lead_pf ~= 0)
                id_bound = find(dme.min_lead_pf ~= 0);
                pfact = sqrt((1./dme.min_lead_pf(id_bound).^2) - 1);  % Factor of active power
                
                ii_inter = repmat((1:length(id_bound)),2,1);
                ii = ii_inter(:);

                jj = [id_bound; id_bound+ninter];

                vals = [-pfact; ones(length(id_bound),1)];

                A_lead = sparse(ii,jj,vals,ninter,2*ninter);
            else
                A_lead = [];
            end

            % Lagging power factor constraints
            ninter = dme.n;
            if any(dme.min_lagg_pf ~= 0)
                id_bound = find(dme.min_lagg_pf ~= 0);
                pfact = sqrt((1./dme.min_lagg_pf(id_bound).^2) - 1);  % Factor of active power
                
                ii_inter = repmat((1:length(id_bound)),2,1);
                ii = ii_inter(:);

                jj = [id_bound; id_bound+ninter];

                vals = [pfact; ones(length(id_bound),1)];

                A_lagg = sparse(ii,jj,vals,ninter,2*ninter);
            else
                A_lagg = [];
            end
        end

        function Q_Smax = capacity_constraints_tdinter(obj,dme, nme)
            %

            id_line3p_slim = find(dme.smax ~=0);

            if ~isempty(id_line3p_slim)
                nb_slim = length(id_line3p_slim);
                Sid = [id_line3p_slim' ; id_line3p_slim'+nme.nk];
                Sid = Sid(:);
                S = mat2cell([Sid Sid 2*ones(2*nb_slim,1)], 2*ones(nb_slim,1));
                nvars = 2*nme.nk;
                Q_Smax = cellfun(@(x)(sparse(x(:,1), x(:,2), x(:,3), nvars, nvars)), S, 'UniformOutput', false);
            else
                Q_Smax = [];
            end
        end

        function obj = data_model_update_on(obj, mm, nm, dm, mpopt)
            dme = obj.data_model_element(dm);
            nme = obj.network_model_element(nm);
            
            ss = nm.get_idx('state');

            for p = 1:nme.nz
                %% buslink complex power flows
                pbl = nm.soln.z(ss.i1.buslink(p):ss.iN.buslink(p));
                qbl = nm.soln.z((ss.i1.buslink(p):ss.iN.buslink(p))+nm.state.N);
                
                %% update in the data model
                dme.tab.(sprintf('pg%d_start', p))(dme.on) = pbl * dm.base_kva;
                dme.tab.(sprintf('qg%d_start', p))(dme.on) = qbl * dm.base_kva;
            end
        end
    end     %% methods
end         %% classdef
