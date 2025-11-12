classdef dme_gen3p_orpd < mp.dme_gen3p
% wgv.dme_gen3p_opd - Data model element for 3-phase generator for Optimal 
%                     Reactive Power Dispatch (ORPD) tasks
%
% Implements the data element model for 3-phase generator elements for ORPD
%
% Adds the following columns in the main data table, found in the
% :attr:`tab` property:
%
%   ================  =========  =====================================
%   Name              Type       Description
%   ================  =========  =====================================
%   ``pmin``          *double*   minimum per-phase active power *(kW)*
%   ``pmax``          *double*   maximum per-phase active power *(kW)*
%   ``qmin``          *double*   minimum per-phase reactive power *(kVAr)*
%   ``qmax``          *double*   maximum per-phase reactive power *(kVAr)*
%   ``a_p``           *double*   second-order coefficient of active power 
%                                generation cost function *(USD/kW^2)*
%   ``b_p``           *double*   first-order coefficient of active power 
%                                generation cost function *(USD/kW)*
%   ``c_p``           *double*   independent coefficient of active power 
%                                generation cost function *(USD)*
%   ``a_q``           *double*   second-order coefficient of reactive power 
%                                generation cost function *(USD/kVAR^2)*
%   ``b_q``           *double*   first-order coefficient of reactive power 
%                                generation cost function *(USD/kVAr)*
%   ``c_q``           *double*   independent coefficient of reactive power 
%                                generation cost function *(USD)*
%   ================  =========  =====================================

%   MATPOWER
%   Copyright (c) 2021-2025, Power Systems Engineering Research Center (PSERC)
%   by Wilson Gonzalez Vanegas, Universidad Nacional de Colombia Sede Manizales
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://matpower.org for more info.

    properties
        pmin    % minimum per-phase active power (kW) for gens that are on
        pmax    % maximum per-phase active power (kW) for gens that are on
        qmin    % minimum per-phase reactive power (kVAr) for gens that are on
        qmax    % maximum per-phase reactive power (kVAr) for gens that are on
        a_p     % second-order active power coefficient (USD/kW^2) for gens that are on
        b_p     % first-order active power coefficient (USD/kW) for gens that are on
        c_p     % independent active power coefficient (USD) for gens that are on
        a_q     % second-order reactive power coefficient (USD/kVAr^2) for gens that are on
        b_q     % first-order reactive power coefficient (USD/kVAr) for gens that are on
        c_q     % independent reactive power coefficient (USD) for gens that are on
    end     %% properties

    methods
        function names = main_table_var_names(obj)
            %
            names = horzcat( main_table_var_names@mp.dm_element(obj), ...
                {'bus', 'vm1_setpoint', 'vm2_setpoint', 'vm3_setpoint', ...
                'pg1', 'pg2', 'pg3', 'qg1', 'qg2', 'qg3', ...
                'pmin', 'pmax', 'qmin', 'qmax', ...
                'a_p', 'b_p', 'c_p', 'a_q', 'b_q', 'c_q'});
        end

        function obj = build_params(obj, dm)
            %
            base_kva = dm.base_kva;

            gen = obj.tab;

            %% get generator parameters
            obj.pg1_start = gen.pg1(obj.on) / base_kva;
            obj.pg2_start = gen.pg2(obj.on) / base_kva;
            obj.pg3_start = gen.pg3(obj.on) / base_kva;
            obj.qg1_start = gen.qg1(obj.on) / base_kva;
            obj.qg2_start = gen.qg2(obj.on) / base_kva;
            obj.qg3_start = gen.qg3(obj.on) / base_kva;
            obj.vm1_setpoint = gen.vm1_setpoint(obj.on);
            obj.vm2_setpoint = gen.vm2_setpoint(obj.on);
            obj.vm3_setpoint = gen.vm3_setpoint(obj.on);
            obj.pmin = gen.pmin(obj.on) / base_kva; 
            obj.pmax = gen.pmax(obj.on) / base_kva;
            obj.qmin = gen.qmin(obj.on) / base_kva;
            obj.qmax = gen.qmax(obj.on) / base_kva;
            obj.a_p = gen.a_p(obj.on) / base_kva^ 2;
            obj.b_p = gen.b_p(obj.on) / base_kva;
            obj.c_p = gen.c_p(obj.on);
            obj.a_q = gen.a_q(obj.on) / base_kva^2;
            obj.b_q = gen.b_q(obj.on) / base_kva;
            obj.c_q = gen.c_q(obj.on);

            obj.apply_vm_setpoint(dm);
        end
    end     %% methods
end         %% classdef
