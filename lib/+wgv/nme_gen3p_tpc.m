classdef nme_gen3p_tpc < mp.nme_gen3p & wgv.form_tpc
% wgv.nme_gen3p_tpc - Network model element for 3-phase generator with
%                      TPC formulation
%
% Sets parameter N for three-phase generators and inherits from wgv.form_tpc.

%   MATPOWER
%   Copyright (c) 2021-2025, Power Systems Engineering Research Center (PSERC)
%   by Wilson Gonzalez Vanegas, Universidad Nacional de Colombia Sede Manizales
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://matpower.org for more info.
    
    methods
        function obj = build_params(obj, nm, dm)
            %
            build_params@mp.nme_gen3p(obj, nm, dm);   %% call parent

            nzgen3p = sum(nm.state.idx.N.gen3p);    % total number of states for 3-phase generators            

            % Indices for N parameter
            ii = 1:nzgen3p;            
            ii = repmat(ii,2,1); ii = ii(:);         
            jj = (1:nzgen3p);
            jj = [jj; nzgen3p + jj]; jj = jj(:);
            Nvals = repmat([-1; -1j], 1, nzgen3p);
            
            % Build N parameter
            obj.N = sparse(ii, jj, Nvals(:), nzgen3p, 2*nzgen3p);
        end
    end
end         %% classdef