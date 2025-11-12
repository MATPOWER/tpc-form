classdef mme_shunt3p_orpd_tpc < wgv.mme_shunt3p_pf_tpc
% wgv.mme_shunt3p_orpd_tpc - Math model element for 3-phase shunt for TPC-based 
%                            Optimal Reactive Power Dispacht (ORPD).
%
% Math model element class for 3-phase shunt elements for ORPD task
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
            vs_zr = struct('name',{'Va3', 'Va3', 'Va3', ...
                                   'LnVm3', 'LnVm3', 'LnVm3', ...
                                   'Pshunt3p', 'Pshunt3p', 'Pshunt3p', ...
                                   }, ...
                            'idx',{{1}, {2}, {3}, {1}, {2}, {3}, ...
                                   {1} {2}, {3}} ...
                                   );
            vs_zi = struct('name',{'Va3', 'Va3', 'Va3', ...
                                   'LnVm3', 'LnVm3', 'LnVm3', ...
                                   'Qshunt3p', 'Qshunt3p', 'Qshunt3p', ...
                                   }, ...
                            'idx',{{1}, {2}, {3}, {1}, {2}, {3}, ...
                                   {1} {2}, {3}} ...
                                   );
            
            % define cuadratic forms for shunt port injections
            [Q_P, Q_Q, B_P, B_Q] = obj.port_inj_shunt3p(nme, 3*nm.elements.bus3p.nk);
            l = nme.s;

            % add quadratic constraints for branch active and reactive power port injections
            mm.qcn.add(mm.var, 'Pshunt3p_inj', Q_P, B_P, -real(l), -real(l), vs_zr);
            mm.qcn.add(mm.var, 'Qshunt3p_inj', Q_Q, B_Q, -imag(l), -imag(l), vs_zi);
        end

        function [Qr, Qi, Br, Bi] = port_inj_shunt3p(obj, nme, nb3p)
            %
            
            % quadratic terms
            nvars = 2*nb3p + nme.np*nme.nk;                
            Q = cellfun(@(x)(sparse(x(:,1), x(:,2), x(:,3), nvars, nvars)), nme.Qu, 'UniformOutput', false);
            Qr = cellfun(@(x)(real(x)), Q, 'UniformOutput', false);
            Qi = cellfun(@(x)(imag(x)), Q, 'UniformOutput', false);

            % linear terms            
            Br = [real(nme.M)  -1*speye(nme.np*nme.nk)];
            Bi = [imag(nme.M)  -1*speye(nme.np*nme.nk)];
        end
        
    end     %% methods
end         %% classdef
