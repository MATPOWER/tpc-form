classdef dme_bus3p_orpd < mp.dme_bus3p
% wgv.dme_bus3p - Data model element for 3-phase bus for Optimal Reactive 
%                 Power Dispatch (ORPD) tasks
%
% Implements the data element model for 3-phase bus elements for ORPD.
%
% Adds the following columns in the main data table, found in the
% :attr:`tab` property:
%
%   ===========  =========  ========================================
%   Name         Type       Description
%   ===========  =========  ========================================
%   ``vmin``     *double*   minimum voltage magnitude *(p.u.)*
%   ``vmax``     *double*   miximum voltage magnitude *(p.u.)*
%   ===========  =========  ========================================

%   MATPOWER
%   Copyright (c) 2021-2025, Power Systems Engineering Research Center (PSERC)
%   by Wilson Gonzalez Vanegas, Universidad Nacional de Colombia Sede Manizales
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://matpower.org for more info.

    properties        
        vmin        % minimum per-phase voltage magnitude (p.u.) for buses that are on
        vmax        % miximum per-phase voltage magnitude (p.u.) for buses that are on
    end     %% properties

    methods        

        function names = main_table_var_names(obj)
            %
            names = horzcat( main_table_var_names@mp.dm_element(obj), ...
                {'type', 'base_kv', 'vm1', 'vm2', 'vm3', 'va1', 'va2', 'va3', 'vmin', 'vmax'});
        end

        function obj = build_params(obj, dm)
            %

            %% initialize voltage from bus table
            bus = obj.tab;
            obj.va1_start = bus.va1(obj.on) * pi/180;
            obj.va2_start = bus.va2(obj.on) * pi/180;
            obj.va3_start = bus.va3(obj.on) * pi/180;
            obj.vm1_start = bus.vm1(obj.on);
            obj.vm2_start = bus.vm2(obj.on);
            obj.vm3_start = bus.vm3(obj.on);
            obj.vmin = bus.vmin(obj.on);
            obj.vmax = bus.vmax(obj.on);

            %% set PV buses without online voltage controls to PQ
            i = find(obj.type == mp.NODE_TYPE.PV & ~obj.vm_control);
            obj.type(i) =  mp.NODE_TYPE.PQ; %% direct assignment to type
                                            %% property (as opposed to use of
                                            %% set_bus_type_pv() method)
                                            %% keeps it from propagating to
                                            %% output tables
        end
    end     %% methods
end         %% classdef
