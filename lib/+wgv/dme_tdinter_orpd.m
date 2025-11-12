classdef dme_tdinter_orpd < wgv.dme_tdinter
% wgv.dme_tdinter_orpd - Data model element for 1-to-3-phase interface for
%                        Optimal Reactive Power Dispatch (ORPD)
%
% Implements the data element model for 1-to-3-phase interface elements for ORPD
%
% Adds the following columns in the main data table, found in the
% :attr:`tab` property:
%
%   ===============  =========  ========================================
%   Name             Type       Description
%   ===============  =========  ========================================
%   ``min_lead_pf``  *double*   minimum power factor for leading exchange
%   ``min_lagg_pf``  *double*   minimum power factor for lagging exchange
%   ``smax``         *double*   maximum per-phase apparent power (kVA)
%   ===============  =========  ========================================

%   MATPOWER
%   Copyright (c) 2021-2025, Power Systems Engineering Research Center (PSERC)
%   by Wilson Gonzalez Vanegas, Universidad Nacional de Colombia Sede Manizales
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://matpower.org for more info.

    properties        
        min_lead_pf   %% maximum power factor for leading power interchange at the TD Interface
        min_lagg_pf   %% maximum power factor for lagging power interchange at the TD Interface
        smax          %% maximum per-phase apparent power (kVA) for TD Interfaces
    end     %% properties

    methods
        function names = main_table_var_names(obj)
            %
            names = horzcat( main_table_var_names@mp.dm_element(obj), ...
                {'bus', 'bus3p', 'min_lead_pf', 'min_lagg_pf', 'smax'});
        end

        function obj = build_params(obj, dm)
            %

            %% check for matching base_kv
            base_kv1 = dm.elements.bus.tab.base_kv(obj.bus);
            base_kv3 = dm.elements.bus3p.tab.base_kv(obj.bus3p);
            if any(base_kv1 ~= base_kv3)
                error('mp.dme_tdinter.build_params: tdinter objects can only link buses with identical base_kv');
            end

            %% check for at least one PQ bus per buslink
            type1 = dm.elements.bus.tab.type(obj.bus);
            type3 = dm.elements.bus3p.tab.type(obj.bus3p);
            if any(type1 ~= mp.NODE_TYPE.PQ & type3 ~= mp.NODE_TYPE.PQ)
                error('mp.dme_tdinter.build_params: at least one of the buses linked by a tdinter must be of type PQ');
            end

            %% check for balanced voltage magnitudes for REF and PV buses
            %% on 3-phase side
            refpv = (type3 == mp.NODE_TYPE.REF | type3 == mp.NODE_TYPE.PV);
            vm1_ref_pv = dm.elements.bus3p.tab.vm1(obj.bus3p(refpv));
            vm2_ref_pv = dm.elements.bus3p.tab.vm2(obj.bus3p(refpv));
            vm3_ref_pv = dm.elements.bus3p.tab.vm3(obj.bus3p(refpv));
            if any(vm1_ref_pv ~= vm2_ref_pv | vm2_ref_pv ~= vm3_ref_pv)
                error('mp.dme_tdinter.build_params: tdinter objects can only link to REF or PV buses with balanced voltage magnitudes');
            end

            %%  Import bounds on leading/lagging power factors and capacity at the TD Interface
            obj.min_lead_pf = obj.tab.min_lead_pf(obj.on);
            obj.min_lagg_pf = obj.tab.min_lagg_pf(obj.on);
            obj.smax = obj.tab.smax(obj.on);
        end
    end     %% methods
end         %% classdef
