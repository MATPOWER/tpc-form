classdef dmce_tdinter_orpd_mpc2 < mp.dmc_element
% wgv.dmce_tdinter_orpd_mpc2 - Data model converter element for 1-to-3-phase buslink 
%                              Optimal Reactive Power Dispatch for MATPOWER case v2.

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
            name = 'buslink';
        end

        function df = data_field(obj)
            %
            df = 'buslink';
        end

        function vmap = table_var_map(obj, dme, mpc)
            %
            vmap = table_var_map@mp.dmc_element(obj, dme, mpc);

            %% mapping for each name, default is {'col', []}
            vmap.uid{2}         = 1;
            vmap.name           = {'cell', ''};     %% empty char
            vmap.status{2}      = 4;
            vmap.source_uid     = {'cell', ''};     %% empty char
            vmap.bus{2}         = 2;
            vmap.bus3p{2}       = 3;
            vmap.min_lead_pf{2} = 5;
            vmap.min_lagg_pf{2} = 6;
            vmap.smax{2}        = 7;
        end
    end     %% methods
end         %% classdef
