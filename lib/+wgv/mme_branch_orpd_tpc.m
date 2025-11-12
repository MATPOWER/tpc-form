classdef mme_branch_orpd_tpc < wgv.mme_branch_pf_tpc
% wgv.mme_branch_orpd_tpc - Math model element for branch for TPC-based
%                           Optimal Reactive Power Dispacht (ORPD).
%
% Math model element class for single-phase branch elements for ORPD task
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
        function obj = add_constraints(obj, mm, nm, dm, mpopt)
            %

            nme = obj.network_model_element(nm);
            dme = obj.data_model_element(dm);

            %% 1) Port injection constraints

            % variables sets
            vs_zr = struct('name',{'Va', 'LnVm', 'Pbranch', 'Pbranch'}, ...
                            'idx',{{}, {}, {1}, {2}});
            vs_zi = struct('name',{'Va', 'LnVm', 'Qbranch', 'Qbranch'}, ...
                            'idx',{{}, {}, {1}, {2}});
            
            % define cuadratic forms for branch port injections
            [Q_P, Q_Q, B_P, B_Q] = obj.port_inj_branch(nme, nm.elements.bus.nk);
            l = nme.s;

            % add quadratic constraints for branch active and reactive power port injections
            mm.qcn.add(mm.var, 'Pbranch_inj', Q_P, B_P, -real(l), -real(l), vs_zr);
            mm.qcn.add(mm.var, 'Qbranch_inj', Q_Q, B_Q, -imag(l), -imag(l), vs_zi);

            %% 2) Capacity constraints

            % variable sets
            vs_from = struct('name',{'Pbranch', 'Qbranch'}, ...
                             'idx',{{1}, {1}});
            vs_to = struct('name',{'Pbranch', 'Qbranch'}, ...
                           'idx',{{2}, {2}});
            
            % define quadratic forms for branch capacity
            Q_Smax = obj.capacity_constraints_branch(dme, nme);
            B_max = [];
            l_max = dme.rate_a(dme.rate_a ~= 0);
            
            % add quadratic constraints for apparent power limits of branches
            if ~isempty(Q_Smax)
                mm.qcn.add(mm.var,'Smax_br_fr',Q_Smax,B_max,-l_max,l_max,vs_from);
                mm.qcn.add(mm.var,'Smax_br_to',Q_Smax,B_max,-l_max,l_max,vs_to);
            end
        end

        function [Qr, Qi, Br, Bi] = port_inj_branch(obj, nme, nb)
            %
            
            % quadratic terms
            nvars = 2*nb + nme.np*nme.nk;                
            Q = cellfun(@(x)(sparse(x(:,1), x(:,2), x(:,3), nvars, nvars)), nme.Qu, 'UniformOutput', false);
            Qr = cellfun(@(x)(real(x)), Q, 'UniformOutput', false);
            Qi = cellfun(@(x)(imag(x)), Q, 'UniformOutput', false);

            % linear terms            
            Br = [real(nme.M)  -1*speye(nme.np*nme.nk)];
            Bi = [imag(nme.M)  -1*speye(nme.np*nme.nk)];
        end

        function Q_Smax = capacity_constraints_branch(obj,dme, nme)
            %

            id_branch_slim = find(dme.rate_a ~=0);

            if ~isempty(id_branch_slim)
                nb_slim = length(id_branch_slim);
                Sid = [id_branch_slim' ; id_branch_slim'+nme.nk];
                Sid = Sid(:);
                S = mat2cell([Sid Sid 2*ones(2*nb_slim,1)], 2*ones(nb_slim,1));
                nvars = 2*nme.nk;
                Q_Smax = cellfun(@(x)(sparse(x(:,1), x(:,2), x(:,3), nvars, nvars)), S, 'UniformOutput', false);
            else
                Q_Smax = [];
            end
        end
    end     %% methods
end         %% classdef