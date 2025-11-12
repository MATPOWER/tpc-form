classdef mme_gen_orpd_tpc < mp.mme_gen
% wgv.mme_gen_orpd_tpc - Math model element for generator for TPC-based 
%                        Optimal Reactive Power Dispatch (ORPD).
%
% Math model element class for generator elements for TPC-based ORPD.
%
% Implements method for updating the output data in the corresponding data
% model element for in-service generators from the math model solution.

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
            if nm.userdata.ishybrid                
                Pg_ref = dme.tab.pg;

                %% Active power generation deviation "costs"
                if isfield(mpopt.orpd.tpc,'cost_dev_Pgen')
                    cost_Pgen_ref = mpopt.orpd.tpc.cost_dev_Pgen;
                else
                    cost_Pgen_ref = 1;
                end

                % define quadratic form for voltage deviation objective function
                [Q_dev, c_dev, k_dev] = obj.Pgen_deviation(dme, Pg_ref, cost_Pgen_ref);

                % add quadratic costs to objective function
                mm.qdc.add(mm.var, 'Pg_dev_Trans', Q_dev, c_dev, k_dev, {'Pg'});
            end
        end

        function [Q_dev, c_dev, k_dev] = Pgen_deviation(obj, dme, Pg_ref, cost)
            %
           
            Q_dev = 2*cost*ones(dme.n,1);
            c_dev = -2*cost*Pg_ref;
            k_dev = cost*Pg_ref.^2;
        end
        function obj = data_model_update_on(obj, mm, nm, dm, mpopt)
            %

            x = nm.soln.x;
            
            pg = x(mm.var.idx.i1.Pg:mm.var.idx.iN.Pg);
            qg = x(mm.var.idx.i1.Qg:mm.var.idx.iN.Qg);

            %% update in the data model
            dme = obj.data_model_element(dm);
            dme.tab.pg(dme.on) = pg * dm.base_mva;
            dme.tab.qg(dme.on) = qg * dm.base_mva;
        end
    end     %% methods
end         %% classdef
