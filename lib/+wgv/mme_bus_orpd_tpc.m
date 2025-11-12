classdef mme_bus_orpd_tpc < mp.mme_bus
% wgv.mme_bus_orpd_tpc - Math model element for bus for TPC-based Optmila 
%                        reactive Power Dispatch (ORPD).
%
% Math model element class for bus elements for TPC-based ORPD.
%
% Implements method for updating the output data in the corresponding data
% model element for in-service buses from the math model solution.

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
        function obj = add_costs(obj, mm, nm, dm, mpopt)
            %
                
            dme = obj.data_model_element(dm);
            if nm.userdata.ishybrid && isfield(mpopt.orpd.tpc,'buses_T_Vref')
                [~,id_bus] = ismember(mpopt.orpd.tpc.buses_T_Vref,dme.tab.uid);
                LnVm_ref = log(dme.tab.vm(id_bus));

                %% Voltage deviation "costs"
                if isfield(mpopt.orpd.tpc,'cost_dev_Vref')
                    cost_Vref = mpopt.orpd.tpc.cost_dev_Vref;
                else
                    cost_Vref = 1;
                end

                % define quadratic form for voltage deviation objective function
                [Q_dev, c_dev, k_dev] = obj.voltage_deviation_bus(dme, id_bus, LnVm_ref, cost_Vref);

                % add quadratic costs to objective function
                mm.qdc.add(mm.var, 'Vm_dev_Trans', Q_dev, c_dev, k_dev, {'LnVm'});
            end
        end

        function [Q_dev, c_dev, k_dev] = voltage_deviation_bus(obj, dme, id_bus, LnVm_ref, cost)
            %
            
            Q_dev = zeros(dme.n,1);
            c_dev = zeros(dme.n,1);
            k_dev = zeros(dme.n,1);

            Q_dev(id_bus) = 2*cost*ones(length(id_bus),1);
            c_dev(id_bus) = -2*cost*LnVm_ref;
            k_dev(id_bus) = cost*LnVm_ref.^2;
        end

        function obj = data_model_update_on(obj, mm, nm, dm, mpopt)
            %
            
            x = nm.soln.x;

            %% bus voltages (angles and log of magnitudes)
            va = x(mm.var.idx.i1.Va:mm.var.idx.iN.Va);
            lnvm = x(mm.var.idx.i1.LnVm:mm.var.idx.iN.LnVm);

            %% update in the data model
            dme = obj.data_model_element(dm);
            dme.tab.va(dme.on) = va * 180/pi;
            dme.tab.vm(dme.on) = exp(lnvm);
        end
    end     %% methods
end         %% classdef
