classdef nme_gen3p_orpd_tpc < wgv.nme_gen3p_tpc
% wgv.nme_gen3p_orpd_tpc - Network model element for 3-phase generator for 
%                          TPC-based optimal Reactive Power Dispatch (ORPD)
%                      
%
% Inherits from wgv.form_tpc and overrides method add_zvars().

%   MATPOWER
%   Copyright (c) 2021-2025, Power Systems Engineering Research Center (PSERC)
%   by Wilson Gonzalez Vanegas, Universidad Nacional de Colombia Sede Manizales
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://matpower.org for more info.
    
    methods
        function obj = add_zvars(obj, nm, dm, idx)
            %
            p = idx{1};
            ng = obj.nk;
            dme = obj.data_model_element(dm);

            pmin = dme.pmin;
            pmax = dme.pmax;
            qmin = dme.qmin;
            qmax = dme.qmax;

            id_no_plim = (pmin==0 & pmax==0);
            pmin(id_no_plim) = 0;
            pmax(id_no_plim) = Inf;

            id_no_qlim = (qmin==0 & qmax==0);
            qmin(id_no_plim) = -Inf;
            qmax(id_no_plim) = Inf;

            pg_start = dme.(sprintf('pg%d_start', p));
            qg_start = dme.(sprintf('qg%d_start', p));

            if p == 1
                nm.init_indexed_name('zr', 'Pg3', {obj.nz});
                nm.init_indexed_name('zi', 'Qg3', {obj.nz});
            end
            nm.add_var('zr', 'Pg3', {p}, ng, pg_start, pmin, pmax);
            nm.add_var('zi', 'Qg3', {p}, ng, qg_start, qmin, qmax);
        end
    end
end         %% classdef