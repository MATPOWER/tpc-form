classdef mme_gen_pf_tpc < mp.mme_gen
% wgv.mme_gen_pf_tpc - Math model element for generator for TPC power flow.
%
% Math model element class for generator elements for TPC power flow problems.
%
% Implements method for updating the output data in the corresponding data
% model element for in-service generators from the math model solution.

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

            %% generator active power
            ss = nm.get_idx('state');

            %% check for hybrid case
            if nm.userdata.ishybrid
                z = mm.aux_data.pm_all_1p3p_z' * nm.soln.z; 
            else
                z = nm.soln.z;
            end
            pg = z(ss.i1.gen:ss.iN.gen);
            qg = z(ss.iN.gen+1:2*ss.iN.gen);

            %% update in the data model
            dme = obj.data_model_element(dm);
            dme.tab.pg(dme.on) = pg * dm.base_mva;
            dme.tab.qg(dme.on) = qg * dm.base_mva;
        end
    end     %% methods
end         %% classdef
