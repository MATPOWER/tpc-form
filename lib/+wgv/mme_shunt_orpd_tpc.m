classdef mme_shunt_orpd_tpc < wgv.mme_shunt_pf_tpc
% wgv.mme_shunt_orpd_tpc - Math model element for shunt for TPC-based
%                          Optimal Reactive Power Dispacht (ORPD).
%
% Math model element class for single-phase shunt elements for ORPD task
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

            % variables sets
            vs_zr = struct('name',{'Va', 'LnVm', 'Pshunt'}, ...
                            'idx',{{}, {}, {}});
            vs_zi = struct('name',{'Va', 'LnVm', 'Qshunt'}, ...
                            'idx',{{}, {}, {}});
            
            % define cuadratic forms for shunt port injections
            [Q_P, Q_Q, B_P, B_Q] = obj.port_inj_shunt(nme, nm.elements.bus.nk);
            l = nme.s;

            % add quadratic constraints for branch active and reactive power port injections
            mm.qcn.add(mm.var, 'Pshunt_inj', Q_P, B_P, -real(l), -real(l), vs_zr);
            mm.qcn.add(mm.var, 'Qshunt_inj', Q_Q, B_Q, -imag(l), -imag(l), vs_zi);
        end

        function [Qr, Qi, Br, Bi] = port_inj_shunt(obj, nme, nb)
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
    end     %% methods
end         %% classdef
