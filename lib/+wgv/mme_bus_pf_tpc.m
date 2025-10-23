classdef mme_bus_pf_tpc < mp.mme_bus
% wgv.mme_bus_pf_tpc - Math model element for bus for TPC power flow.
%
% Math model element class for bus elements for TPC power flow problems.
%
% Implements method for updating the output data in the corresponding data
% model element for in-service buses from the math model solution.

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
            
            nn = nm.get_idx('node');           

            %% check for hybrid case
            if nm.userdata.ishybrid
                u = mm.aux_data.pm_all_va_lnvm_u*nm.soln.u;
            else
                u = nm.soln.u;                
            end

            %% bus voltages (angles and log of magnitudes)
            va = u(nn.i1.bus:nn.iN.bus);
            lnvm = u(nn.iN.bus+1:2*nn.iN.bus);

            %% update in the data model
            dme = obj.data_model_element(dm);
            dme.tab.va(dme.on) = va * 180/pi;
            dme.tab.vm(dme.on) = exp(lnvm);
        end
    end     %% methods
end         %% classdef
