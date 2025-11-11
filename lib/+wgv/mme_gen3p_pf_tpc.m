classdef mme_gen3p_pf_tpc < mp.mm_element
% wgv.mme_gen3p_pf_tpc - Math model element for 3-phase generator for 
%                        TPC power flow.
%
% Math model element class for 3-phase generator elements for power flow
% tasks under a TPC formulation.
%
% Implements method for updating the output data in the corresponding data
% model element for in-service 3-phase generators from the math model solution.

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
            name = 'gen3p';
        end

        function obj = data_model_update_on(obj, mm, nm, dm, mpopt)
            %

            dme = obj.data_model_element(dm);
            nme = obj.network_model_element(nm);

            ss = nm.get_idx('state');

            for p = 1:nme.nz
                %% generator active power
                ngen3p = sum(ss.N.gen3p);
                pg = nm.soln.z(ss.i1.gen3p(p):ss.iN.gen3p(p));
                qg = nm.soln.z(ss.i1.gen3p(p)+ngen3p:ss.iN.gen3p(p)+ngen3p);

                %% update in the data model
                dme.tab.(sprintf('pg%d', p))(dme.on) = pg * dm.base_kva;
                dme.tab.(sprintf('qg%d', p))(dme.on) = qg * dm.base_kva;
            end
        end
    end     %% methods
end         %% classdef
