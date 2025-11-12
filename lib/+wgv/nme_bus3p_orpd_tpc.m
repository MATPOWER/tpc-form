classdef nme_bus3p_orpd_tpc < mp.nme_bus3p & wgv.form_tpc
% wgv.nme_bus3p_orpd_tpc - Network model element for 3-phase bus for TPC-based
%                          Optimal Reactive Power Dispatch (ORPD) 
%
% Adds voltage variables ``Va3`` and ``LnVm3`` to the network model and inherits
% from wgv.form_tpc.

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
        function obj = add_vvars(obj, nm, dm, idx)
            %
            dme = obj.data_model_element(dm);
            nb = obj.nk;
            p = idx{1};

            %% prepare angle bounds for ref buses
            ref = dme.type == mp.NODE_TYPE.REF;
            va_lb = -pi * ones(nb, 1);
            va_ub =  pi * ones(nb, 1);
            lnvm_start = log(dme.(sprintf('vm%d_start', p)));
            va_start = dme.(sprintf('va%d_start', p));           
            va_lb(ref) = va_start(ref);
            va_ub(ref) = va_start(ref);
            
            %% prepare voltage magnitude bounds
            lnvm_lb = log(dme.vmin);
            lnvm_ub = log(dme.vmax);

            if p == 1
                nm.init_indexed_name('va', 'Va3', {obj.nn});
                nm.init_indexed_name('lnvm', 'LnVm3', {obj.nn});
            end
            nm.add_var('va', 'Va3', {p}, nb, va_start, va_lb, va_ub);
            nm.add_var('lnvm', 'LnVm3', {p}, nb, lnvm_start, lnvm_lb, lnvm_ub);
        end
    end     %% methods
end         %% classdef