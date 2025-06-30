classdef nme_bus_tpc < mp.nme_bus & wgv.form_tpc
% wgv.nme_bus_tpc - Network model element for bus for TPC formulations.
%
% Adds transformed voltage variables ``Va`` and ``LnVm`` to the network model 
% and inherits from wgv.form_tpc.

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
        function obj = add_vvars(obj, nm, dm, idx)% OK-WGV
            %
            dme = obj.data_model_element(dm);
            nb = obj.nk;

            %% prepare angle bounds for ref buses
            va_lb = -Inf(nb, 1);
            va_ub =  Inf(nb, 1);            
            k = find(dme.type == mp.NODE_TYPE.REF);
            va_lb(k) = dme.va_start(k);
            va_ub(k) = dme.va_start(k);

            nm.add_var('va', 'Va', nb, dme.va_start, va_lb, va_ub);
            nm.add_var('lnvm', 'LnVm', nb, log(dme.vm_start), log(dme.vm_lb), log(dme.vm_ub));
        end
    end     %% methods
end         %% classdef
