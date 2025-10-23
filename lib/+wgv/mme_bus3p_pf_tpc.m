classdef mme_bus3p_pf_tpc < mp.mm_element
% wgv.mme_bus3p_pf_tpc - Math model element for 3-phase bus for TPC
%                        power flow.
%
% Math model element class for 3-phase bus elements for power flow task
% under a TPC formulation.
%
% Implements method for updating the output data in the corresponding data
% model element for in-service 3-phase buses from the math model solution.

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
        function name = name(obj)
            %
            name = 'bus3p';
        end

        function obj = data_model_update_on(obj, mm, nm, dm, mpopt)
            %

            dme = obj.data_model_element(dm);
            nme = obj.network_model_element(nm);

            nn = nm.get_idx('node');

            %% check for hybrid case            
            u = nm.soln.u;
            nbus = nm.node.N;
            
            rot = [0 -2*pi/3 2*pi/3];
            for p = 1:nme.nn
                %% complex bus voltages                
                va_p = u(nn.i1.bus3p(p):nn.iN.bus3p(p)) + rot(p);
                lnvm_p = u(nn.i1.bus3p(p)+nbus:nn.iN.bus3p(p)+nbus);

                %% update in the data model
                dme.tab.(sprintf('va%d', p))(dme.on) = va_p * 180/pi;
                dme.tab.(sprintf('vm%d', p))(dme.on) = exp(lnvm_p);
            end
        end
    end     %% methods
end         %% classdef
