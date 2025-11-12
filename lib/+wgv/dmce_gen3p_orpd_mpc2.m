classdef dmce_gen3p_orpd_mpc2 < mp.dmc_element
% wgv.dmce_gen3p_orpd_mpc2 - Data model converter element for 3-phase generator 
%                            for Optimal Reactive Power Dispatch (ORPD) tasks
%                            for |MATPOWER| case v2.

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

        function df = data_field(obj)
            %
            df = 'gen3p';
        end

        function vmap = table_var_map(obj, dme, mpc)
            %
            vmap = table_var_map@mp.dmc_element(obj, dme, mpc);

            %% mapping for each name, default is {'col', []}
            vmap.uid{2}             = 1;
            vmap.name               = {'cell', ''};     %% empty char
            vmap.status{2}          = 3;
            vmap.source_uid         = {'cell', ''};     %% empty char
            vmap.bus{2}             = 2;
            vmap.vm1_setpoint{2}    = 4;
            vmap.vm2_setpoint{2}    = 5;
            vmap.vm3_setpoint{2}    = 6;
            vmap.pg1{2}             = 7;
            vmap.pg2{2}             = 8;
            vmap.pg3{2}             = 9;
            vmap.qg1{2}             = 10;
            vmap.qg2{2}             = 11;
            vmap.qg3{2}             = 12;
            vmap.pmin{2}            = 13;
            vmap.pmax{2}            = 14;
            vmap.qmin{2}            = 15;
            vmap.qmax{2}            = 16;
            vmap.a_p{2}             = 17;
            vmap.b_p{2}             = 18;
            vmap.c_p{2}             = 19;
            vmap.a_q{2}             = 20;
            vmap.b_q{2}             = 21;
            vmap.c_q{2}             = 22;
        end
    end     %% methods
end         %% classdef
