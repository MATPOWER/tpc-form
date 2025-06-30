classdef nme_gen_tpc < mp.nme_gen & wgv.form_tpc
% wgv.nme_gen_tpc - Network model element for gen for TPC formulations.
%
% Adds non-voltage state variables ``Pg`` and ``Qg`` to the network model
% and builds the parameter :math:`\NN`.

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
        function obj = add_zvars(obj, nm, dm, idx)
            %
            ng = obj.nk;
            dme = obj.data_model_element(dm);
            nm.add_var('zr', 'Pg', ng, dme.pg_start, dme.pg_lb, dme.pg_ub);
            nm.add_var('zi', 'Qg', ng, dme.qg_start, dme.qg_lb, dme.qg_ub);
        end

        function obj = build_params(obj, nm, dm)% IN PROGRESS - WGV
            %
            build_params@mp.nme_gen(obj, nm, dm);   %% call parent
            
            nzgen = nm.state.idx.N.gen;           % total number of states for 1-phase generators            

            % Indices for N parameter
            ii = 1:nzgen;            
            ii = repmat(ii,2,1); ii = ii(:);         
            jj = (1:nzgen);
            jj = [jj; nzgen + jj]; jj = jj(:);
            Nvals = repmat([-1; -1j], 1, nzgen);
            
            % Build N parameter            
            obj.N = sparse(ii, jj, Nvals(:), nzgen, 2*nzgen);
        end
    end     %% methods
end         %% classdef
