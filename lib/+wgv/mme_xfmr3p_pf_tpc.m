classdef mme_xfmr3p_pf_tpc < mp.mm_element
% wgv.mme_xfmr3p_pf_tpc - Math model element for 3-phase transformer 
%                         for TPC power flow.
%
% Math model element base class for 3-phase transformer elements for power 
% flow task under the TPC formulation.
%
% Implements method for updating the output data in the corresponding data
% model element for in-service 3-phase transformer from the math model solution.

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
            name = 'xfmr3p';
        end

        function obj = data_model_update_on(obj, mm, nm, dm, mpopt)
            %

            dme = obj.data_model_element(dm);
            nme = obj.network_model_element(nm);

            pp = nm.get_idx('port');
            nn = nme.np/2;

            %% branch active power flow
            for p = 1:nn
                s_fr = nm.soln.gs_(pp.i1.xfmr3p(p):pp.iN.xfmr3p(p));
                s_to = nm.soln.gs_(pp.i1.xfmr3p(p+nn):pp.iN.xfmr3p(p+nn));

                %% update in the data model
                dme.tab.(sprintf('pl%d_fr', p))(dme.on) = real(s_fr) * dm.base_kva;
                dme.tab.(sprintf('ql%d_fr', p))(dme.on) = imag(s_fr) * dm.base_kva;
                dme.tab.(sprintf('pl%d_to', p))(dme.on) = real(s_to) * dm.base_kva;
                dme.tab.(sprintf('ql%d_to', p))(dme.on) = imag(s_to) * dm.base_kva;
            end
        end
    end     %% methods
end         %% classdef
