classdef nme_shunt3p_orpd_tpc < wgv.nme_shunt3p_tpc
% wgv.nme_shunt3p_orpd_tpc - Network model element for 3-phase shunt for
%                            TPC-based Optimal Reactive Power Dispatch (ORPD)
%
% Inherits from wgv.nme_shunt3p_tpc and implements a method for adding
% non-voltage variables that model the active and reactive power injections
% at each port of the three-phase shunt (three ports).

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
        function nz = nz(obj)
            %
            nz = 3;     %% 6 complex non-voltage states per element
        end

        function obj = add_zvars(obj, nm, dm, idx)
            %
            
            dme = obj.data_model_element(dm);

            p = idx{1};

            if p == 1
                nm.init_indexed_name('zr', 'Pshunt3p', {obj.nz});
                nm.init_indexed_name('zi', 'Qshunt3p', {obj.nz});
            end
            nm.add_var('zr', 'Pshunt3p', {p}, obj.nk, 0, -Inf, Inf);
            nm.add_var('zi', 'Qshunt3p', {p}, obj.nk, 0, -Inf, Inf);
        end
    end
end