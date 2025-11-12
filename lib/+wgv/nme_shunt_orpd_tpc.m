classdef nme_shunt_orpd_tpc < wgv.nme_shunt_tpc
% wgv.nme_shunt_orpd_tpc - Network model element for shunt for TPC-based 
%                          Optimal Reactive Power Dispatch (ORPD)
%
% Inherits from wgv.nme_shunt_tpc and implements a method for adding
% non-voltage variables that model the active and reactive power injections
% at the single-phase shunt port.

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
            nz = 1;     %% 2 complex non-voltage states per element
        end        

        function obj = add_zvars(obj, nm, dm, idx)
            %

            nm.add_var('zr', 'Pshunt', obj.nk, 0, -Inf, Inf);
            nm.add_var('zi', 'Qshunt', obj.nk, 0, -Inf, Inf);
        end
    end
end