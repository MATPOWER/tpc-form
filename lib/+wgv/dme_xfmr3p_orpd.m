classdef dme_xfmr3p_orpd < mp.dme_xfmr3p
% wgv.dme_xfmr3p_orpd - Data model element for 3-phase transformer for Optimal Reactive 
%                       Power Dispatch (ORPD) tasks
%
% Implements the data element model for 3-phase transformer elements for ORPD.
%
% Adds the following columns in the main data table, found in the
% :attr:`tab` property:
%
%   ===========  =========  =============================================
%   Name         Type       Description
%   ===========  =========  =============================================
%   ``smax``     *double*   maximum per-phase apparent power *(kVA)* for
%                           three-phase transformer
%   ===========  =========  =============================================

%   MATPOWER
%   Copyright (c) 2021-2025, Power Systems Engineering Research Center (PSERC)
%   by Wilson Gonzalez Vanegas, Universidad Nacional de Colombia Sede Manizales
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://matpower.org for more info.

    properties
        smax      % maximum per-phase apparent power (kVA) for three-phase transformers that are on
    end     %% properties

    methods
        function names = main_table_var_names(obj)
            %
            names = horzcat( main_table_var_names@mp.dm_element(obj), ...
                {'bus_fr', 'bus_to', 'r', 'x', 'base_kva', 'base_kv', 'tm', ...
                 'pl1_fr', 'ql1_fr', 'pl2_fr', 'ql2_fr', 'pl3_fr', 'ql3_fr', ...
                 'pl1_to', 'ql1_to', 'pl2_to', 'ql2_to', 'pl3_to', 'ql3_to', ...
                 'smax'});
        end

        function obj = build_params(obj, dm)
            %

            obj.r = obj.tab.r(obj.on);
            obj.x = obj.tab.x(obj.on);
            obj.base_kva = obj.tab.base_kva(obj.on);
            obj.base_kv  = obj.tab.base_kv(obj.on);
            obj.tm = obj.tab.tm(obj.on);
            obj.smax = obj.tab.smax(obj.on) / dm.base_kva;
        end
    end     %% methods
end         %% classdef
