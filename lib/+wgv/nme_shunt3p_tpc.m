classdef nme_shunt3p_tpc < mp.nm_element & wgv.form_tpc
% wgv.nme_shunt3p_tpc - Network model element for three-phase shunt element
%                       with TPC formulation
%

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
            name = 'shunt3p';
        end

        function np = np(obj)
            %
            np = 3;     %% this is a 3 port element
        end

        function obj = build_params(obj, nm, dm)
            %
            build_params@mp.nm_element(obj, nm, dm); %% call parent
            dme = obj.data_model_element(dm);

            gs = [dme.gs1 dme.gs2 dme.gs3];
            bs = [dme.bs1 dme.bs2 dme.bs3];
            Ysh = gs + 1j * bs;                     %% shunt admittances
            nsh = obj.nk;                           %% number of shunt elements            
            nb = (nm.node.N)/3;                     %% total number of buses
            id_sh = dme.bus;                        %% shunt buses             

            %% 1) Qu parameter            
            ii = repmat(id_sh,1,3) + repmat([3*nb 4*nb 5*nb],nsh,1);
            jj = ii;
            Qu_vals = 4*conj(Ysh);
            obj.Qu = mat2cell([ii(:) jj(:) Qu_vals(:)], ones(3*nsh, 1));

            %% 2) M parameter
            jj = ii(:);
            ii = (1:3*nsh)';
            M_vals = 2*conj(Ysh);
            obj.M = sparse(ii, jj, M_vals(:), 3*nsh, 2*nm.node.N);

            %% 3) s parameter
            obj.s = conj(Ysh(:));
        end     

    end     %% methods
end         %% classdef