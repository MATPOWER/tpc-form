classdef mme_bus3p_orpd_tpc < wgv.mme_bus3p_pf_tpc
% wgv.mme_bus3p_orpd_tpc - Math model element for 3-phase bus for TPC-based
%                          Optimal Reactive Power Dispatch (ORPD).
%
% Math model element class for 3-phase bus elements for ORPD task
% under a TPC formulation.
%
% Implements methods for adding voltage deviations terms to the objective 
% function, for forming an interior initial point, and for updating
% the output data in the corresponding data model element for in-service
% buses from the math model solution.

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

            if strcmp(mpopt.orpd.tpc.Vref,'DEFAULT')
                LnVm_ref = zeros(dme.n,3);
            else
                Vref = mpopt.orpd.tpc.Vref;
                if size(Vref,1) == dme.n && size(Vref,2) == 3
                    LnVm_ref = log(Vref);
                else
                    error('wgv.mme_bus3p_orpd_tpc.add_costs: mpopt.orpd.tpc.Vref must be a %dx3 matrix of reference voltages', dme.n);
                end                
            end
            
            %% Voltage deviation "costs" 
            
            % variable set
            vs = struct('name',{'LnVm3', 'LnVm3', 'LnVm3'}, ...
                         'idx',{{1}, {2}, {3}});
            
            % define quadratic form for voltage deviation objective function
            [Q_dev, c_dev, k_dev] = obj.voltage_deviation_bus3p(dme, LnVm_ref);

            % add quadratic costs to objective function
            mm.qdc.add(mm.var, 'Vm_dev', Q_dev, c_dev, k_dev, vs);
        end

        function [Q_dev, c_dev, k_dev] = voltage_deviation_bus3p(obj, dme, LnVm_ref)
            %

            Q_dev = repmat(2*ones(dme.n,1),3,1);
            c_dev = -2*LnVm_ref(:);
            k_dev = LnVm_ref(:).^2;
        end
    end     %% methods
end         %% classdef
