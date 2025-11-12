classdef mme_gen3p_orpd_tpc < mp.mm_element
% wgv.mme_gen3p_orpd_tpc - Math model element for 3-phase generator for 
%                          TPC-based Optimal Reactive Power Dispacht (ORPD).
%
% Math model element class for 3-phase gen elements for ORPD task
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
            name = 'gen3p';
       end
        
       function obj = add_constraints(obj, mm, nm, dm, mpopt)
            %
            
            dme = obj.data_model_element(dm);

            %% 1) Voltage constraints

            % variables sets            
            vs_vm = struct('name', {'LnVm3', 'LnVm3', 'LnVm3'}, ...
                            'idx', {{1}, {2}, {3}});
            vs_va = struct('name', {'Va3', 'Va3', 'Va3'}, ...
                            'idx', {{1}, {2}, {3}});
            vs_pg = struct('name', {'Pg3', 'Pg3', 'Pg3'}, ...
                            'idx', {{1}, {2}, {3}});
            vs_qg = struct('name', {'Qg3', 'Qg3', 'Qg3'}, ...
                            'idx', {{1}, {2}, {3}});

            % define linear parameters for voltage equality constraints
            [A_v, b_v] = obj.voltage_constraints_gen3p(dme, dm.elements.bus3p.tab.uid);

            % add linear constraints for equality of voltage angles and magnitudes
            if mpopt.orpd.tpc.gen3p_const.vm
                mm.lin.add(mm.var,'bal_lnvm_gen3p', A_v, b_v, b_v, vs_vm);
            end
            if mpopt.orpd.tpc.gen3p_const.va
                mm.lin.add(mm.var,'bal_va_gen3p', A_v, b_v, b_v, vs_va);
            end

            % define linear parameters for power equality constraints
            [A_v, b_v] = obj.power_constraints_gen3p(dme);

            % add linear constraints for equality of active and reactive power          
            if mpopt.orpd.tpc.gen3p_const.pg
                mm.lin.add(mm.var,'bal_pg_gen3p', A_v, b_v, b_v, vs_pg);
            end
            if mpopt.orpd.tpc.gen3p_const.qg
                mm.lin.add(mm.var,'bal_qg_gen3p', A_v, b_v, b_v, vs_qg);
            end
        end

        function [A, b] = voltage_constraints_gen3p(obj, dme, id_bus3p)
            %            
            
            nb3p = length(id_bus3p);
            
            [~, bus3p] = ismember(dme.tab.bus, id_bus3p);
             
            jj3p = repmat(bus3p',4,1) + repmat([0 nb3p nb3p 2*nb3p]', 1, length(bus3p));

            jj = jj3p(:);
            ii3p = repmat((1:2*length(bus3p)),2,1);
            ii = ii3p(:);
            
            vals = repmat([1;-1], 1, 2*dme.n);
            
            nvars = 3*nb3p;
            A = sparse(ii, jj, vals(:), 2*dme.n, nvars);

            b = zeros(2*dme.n,1);
        end

        function [A, b] = power_constraints_gen3p(obj, dme)
            %            
            
            ng3p = dme.n;
            
            idgen3p = (1:ng3p)';
             
            jj3p = repmat(idgen3p',4,1) + repmat([0 ng3p ng3p 2*ng3p]', 1, ng3p);

            jj = jj3p(:);
            ii3p = repmat((1:2*length(idgen3p)),2,1);
            ii = ii3p(:);
            
            vals = repmat([1;-1], 1, 2*dme.n);
            
            nvars = 3*ng3p;
            A = sparse(ii, jj, vals(:), 2*dme.n, nvars);

            b = zeros(2*dme.n,1);
        end

        function obj = data_model_update_on(obj, mm, nm, dm, mpopt)
            %

            dme = obj.data_model_element(dm);
            nme = obj.network_model_element(nm);

            x = nm.soln.x;

            for p = 1:nme.nz
                %% generator active and reactive power
                pg = x(mm.var.idx.i1.Pg3(p):mm.var.idx.iN.Pg3(p));
                qg = x(mm.var.idx.i1.Qg3(p):mm.var.idx.iN.Qg3(p));

                %% update in the data model                
                dme.tab.(sprintf('pg%d', p))(dme.on) = pg * dm.base_kva;
                dme.tab.(sprintf('qg%d', p))(dme.on) = qg * dm.base_kva;
            end
       end
    end     %% methods
end         %% classdef
