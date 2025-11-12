classdef dme_shunt3p_orpd < mp.dme_shunt3p
% wgv.dme_shunt3p_orpd - Data model element for 3-phase shunt for Optimal 
%                        Reactive Power Dispatch (ORPD) tasks.
%
% Implements the data element model for 3-phase shunt elements for ORPD
%
% Adds the following columns in the main data table, found in the
% :attr:`tab` property:
%
%   ========  =========  =====================================
%   Name      Type       Description
%   ========  =========  =====================================
%   ``qmin``  *double*   minimum per-phase reactive power injection, 
%                        *(kVAr)*, specified at off-nominal voltage
%                        operation
%   ``qmax``  *double*   maximum per-phase reactive power injection, 
%                        *(kVAr)*, specified at off-nominal voltage
%                        operation
%   ========  =========  =====================================
%
% .. [#] *Nominal* means for a voltage of 1 p.u.

%   MATPOWER
%   Copyright (c) 2021-2025, Power Systems Engineering Research Center (PSERC)
%   by Wilson Gonzalez Vanegas, Universidad Nacional de Colombia Sede Manizales
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://matpower.org for more info.

    properties
        qmin    % minimum per-phase reactive power injection for three-phase shunts that are on
        qmax    % maximum per-phase reactive power injection for three-phase shunts that are on
    end     %% properties

    methods
        function names = main_table_var_names(obj)
            %
            names = horzcat( main_table_var_names@mp.dm_element(obj), ...
                {'bus', 'gs1', 'gs2', 'gs3', 'bs1', 'bs2', 'bs3', ...
                 'p1', 'p2', 'p3', 'q1', 'q2', 'q3', 'qmin', 'qmax'});
        end
        
        function obj = build_params(obj, dm)
            %
            obj.gs1 = obj.tab.gs1(obj.on) / dm.base_kva;
            obj.gs2 = obj.tab.gs2(obj.on) / dm.base_kva;
            obj.gs3 = obj.tab.gs3(obj.on) / dm.base_kva;
            obj.bs1 = obj.tab.bs1(obj.on) / dm.base_kva;
            obj.bs2 = obj.tab.bs2(obj.on) / dm.base_kva;
            obj.bs3 = obj.tab.bs3(obj.on) / dm.base_kva;
            obj.qmin = obj.tab.qmin(obj.on) / dm.base_kva;
            obj.qmax = obj.tab.qmax(obj.on) / dm.base_kva;
        end
    end     %% methods
end         %% classdef
