classdef nme_load_tpc < mp.nme_load & wgv.form_tpc
% wgv.nme_load_tpc - Network model element for load for TPC formulations.
%
% Builds the parameter :math:`\sv` and inherits from wgv.form_tpc.

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
        function obj = build_params(obj, nm, dm)
            %
            build_params@mp.nme_load(obj, nm, dm);  %% call parent

            dme = obj.data_model_element(dm);

            %% constant complex power demand
            obj.s = dme.pd + 1j * dme.qd;            
        end
    end     %% methods
end         %% classdef
