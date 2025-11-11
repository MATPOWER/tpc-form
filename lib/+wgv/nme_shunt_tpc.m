classdef nme_shunt_tpc < mp.nme_shunt & wgv.form_tpc
% wgv.nme_shunt_tpc - Network model element for shunt for TPC formulations.
%
% Implements building of the shunt parameters :math:`\WQu`, math:`\M`, and
% :math:`\s`, and inherits from :class:`mp.wgv_form_tpc`.

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
            build_params@mp.nme_shunt(obj, nm, dm); %% call parent

            dme = obj.data_model_element(dm);
            Ysh = dme.gs + 1j * dme.bs;             %% shunt admittances
            nsh = obj.nk;                           %% number of shunt elements            
            nb = nm.node.idx.N.bus;                         %% total number of buses
            id_sh = dme.bus;                        %% shunt buses             

            %% 1) Qu parameter            
            ii = id_sh + nb;
            jj = ii;
            Qu_vals = 4*conj(Ysh);
            obj.Qu = mat2cell([ii jj Qu_vals], ones(nsh, 1));

            %% 2) M parameter
            ii = (1:nsh)';
            jj = id_sh + nb;
            obj.M = sparse(ii, jj, 2*conj(Ysh), nsh, 2*nb);

            %% 3) s parameter
            obj.s = conj(Ysh);
        end     

    end     %% methods
end         %% classdef