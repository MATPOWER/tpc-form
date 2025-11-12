classdef mme_line3p_orpd_tpc < wgv.mme_line3p_pf_tpc
% wgv.mme_line3p_orpd_tpc - Math model element for 3-phase line for TPC-based 
%                           Optimal Reactive Power Dispacht (ORPD).
%
% Math model element class for 3-phase line elements for ORPD task
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
            vs_zr = struct('name',{'Va3', 'Va3', 'Va3', ...
                                   'LnVm3', 'LnVm3', 'LnVm3', ...
                                   'Pline3p', 'Pline3p', 'Pline3p', 'Pline3p', 'Pline3p', 'Pline3p', ...
                                   }, ...
                            'idx',{{1}, {2}, {3}, {1}, {2}, {3}, ...
                                   {1} {2}, {3}, {4}, {5}, {6}} ...
                                   );
            vs_zi = struct('name',{'Va3', 'Va3', 'Va3', ...
                                   'LnVm3', 'LnVm3', 'LnVm3', ...                                   
                                   'Qline3p', 'Qline3p', 'Qline3p', 'Qline3p', 'Qline3p', 'Qline3p' ...
                                   }, ...
                            'idx',{{1}, {2}, {3}, {1}, {2}, {3}, ...
                                   {1} {2}, {3}, {4}, {5}, {6}} ...
                                   );
            
            % define cuadratic forms for branch port injections
            [Q_P, Q_Q, B_P, B_Q] = obj.port_inj_line3p(nme, 3*nm.elements.bus3p.nk);
            l = nme.s;

            % add quadratic constraints for branch active and reactive power port injections
            mm.qcn.add(mm.var, 'Pline3p_inj', Q_P, B_P, -real(l), -real(l), vs_zr);
            mm.qcn.add(mm.var, 'Qline3p_inj', Q_Q, B_Q, -imag(l), -imag(l), vs_zi);

            %% 2) Capacity constraints

            % define quadratic forms for three-phase line capacity
            Q_Smax = obj.capacity_constraints_line3p(dme, nme);
            B_max = [];
            l_max = dme.smax(dme.smax ~= 0);
            
            % add quadratic constraints for apparent power limits of three-phase lines
            if ~isempty(Q_Smax)
                for p = 1:nme.np
                    vs = struct('name',{'Pline3p', 'Qline3p'}, ...
                                        'idx',{{p}, {p}});

                    mm.qcn.add(mm.var,sprintf('Smax_line3p_%d',p),Q_Smax,B_max,-l_max,l_max,vs);
                end
            end
        end

        function [Qr, Qi, Br, Bi] = port_inj_line3p(obj, nme, nb3p)
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

        function Q_Smax = capacity_constraints_line3p(obj,dme, nme)
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
    end     %% methods
end         %% classdef
