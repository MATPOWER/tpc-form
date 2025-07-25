classdef mme_shunt_pf_tpc < mp.mme_shunt
% wgv.mme_shunt_pf_tpc - Math model element for shunt for TPC power flow.
%
% Math model element class for shunt elements for TPC power flow problems.
%
% Implements method for updating the output data in the corresponding data
% model element for in-service shunts from the math model solution.

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
        function obj = data_model_update_on(obj, mm, nm, dm, mpopt)
            %

            %% shunt complex power consumption
            pp = nm.get_idx('port');
            S = nm.soln.gs_(pp.i1.shunt:pp.iN.shunt);

            %% update in the data model
            dme = obj.data_model_element(dm);
            dme.tab.p(dme.on) = real(S) * dm.base_mva;
            dme.tab.q(dme.on) = imag(S) * dm.base_mva;
        end
    end     %% methods
end         %% classdef
