classdef xt_tpc < mp.extension
% wgv.xt_tpc - |MATPOWER| extension to use an approximate formulation in
%              Transformer Polar Coordinates (TPC) for power flow tasks.
%
% For power flow tasks, overrides the mp.task_pf class by the wgv.task_pf
% class included in the +wgv package, which is part of the tpc-form
% repository.
%
% No changes are required for the converter, model, and element classes.
%
% See also mp.extension.

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
        function task_class = task_class(obj, task_class, mpopt)
            % Adds the task class for PF under TPC formulations
            
            task_class = @wgv.task_pf;
        end
    end
end