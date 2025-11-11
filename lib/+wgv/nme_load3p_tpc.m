classdef nme_load3p_tpc < mp.nm_element & wgv.form_tpc
% wgv.nme_load3p_tpc - Network model element for 3-phase load.
%
% Implements the network model element for 3-phase load elements, with
% 3 ports per 3-phase load. 
%
% This class is included for the sake of completeness. Notice that this is 
% an exact copy of mp.nme_load3p, because loads are known quantities that
% do not require any TPC approximation.

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
            name = 'load3p';
        end

        function np = np(obj)
            %
            np = 3;     %% this is a 3 port element
        end

        function obj = build_params(obj, nm, dm)
            %
            build_params@mp.nm_element(obj, nm, dm);    %% call parent
            dme = obj.data_model_element(dm);

            %% constant complex power demand
            pd = [dme.pd1 dme.pd2 dme.pd3];
            qd = pd .* tan(acos( [dme.pf1 dme.pf2 dme.pf3] ));

            obj.s = pd(:) + 1j * qd(:);
        end
    end     %% methods
end         %% classdef
