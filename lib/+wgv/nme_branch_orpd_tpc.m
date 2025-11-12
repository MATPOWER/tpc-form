classdef nme_branch_orpd_tpc < wgv.nme_branch_tpc
% wgv.nme_branch_orpd_tpc - Network model element for branch for TPC-based 
%                           Optimal Reactive Power Dispatch (ORPD)
%
% Inherits from wgv.nme_branch_tpc and implements a method for adding
% non-voltage variables that model the active and reactive power injections
% at each port of the single-phase branch (two ports).

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
        function nz = nz(obj)
            %
            nz = 2;     %% 2 complex non-voltage states per element
        end        

        function obj = add_zvars(obj, nm, dm, idx)
            %
            
            dme = obj.data_model_element(dm);

            smin = -dme.rate_a;
            smax = dme.rate_a;

            id_nolim = (smin==0 & smax==0);
            smin(id_nolim) = -Inf;
            smax(id_nolim) = Inf;

            p = idx{1};

            if p == 1
                nm.init_indexed_name('zr', 'Pbranch', {obj.nz});
                nm.init_indexed_name('zi', 'Qbranch', {obj.nz});
            end
            nm.add_var('zr', 'Pbranch', {p}, obj.nk, 0, smin, smax);
            nm.add_var('zi', 'Qbranch', {p}, obj.nk, 0, smin, smax);
        end
    end
end