classdef form_tpc < mp.form
% form_tpc - Base class for |MATPOWER| Transformed Polar Coordinates
%            (TPC)-based  **formulations**.
%
% Used as a mix-in class for all **network model element** classes
% with an TPC network model formulation. That is, each concrete network model
% element class with a LIN or QUAD  TPC-formulation must inherit, at least 
% indirectly, from both mp.nm_element and wgv.form_tpc.
%
% wgv.form_tpc defines the complex port injection as a quadratic form
% of the state variables X, that is, the voltage angles va and the logarithm 
% of the voltage magnitudes \nu, and the non-voltage states Z.
%
% See also mp.form, mp.form_ac, mp.nm_element.

%   MATPOWER
%   Copyright (c) 2021-2025, Power Systems Engineering Research Center (PSERC)
%   by Wilson Gonzalez Vanegas, Universidad Nacional de Colombia Sede Manizales
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://matpower.org for more info.

    properties
        % 
        Qu = {};

        %
        Qz = {};

        % 
        M = [];

        %
        N = [];

        %
        s = [];
    end

    methods
        function name = form_name(obj)
            % Get user-readable name of formulation, i.e. ``'TPC'``.
            %
            % See :meth:`mp.form.form_name`.

            name = 'TPC';
        end

        function tag = form_tag(obj)
            % Get short label of formulation, i.e. ``'tpc'``.
            %
            % See :meth:`mp.form.form_tag`.

            tag = 'tpc';
        end

        function params = model_params(obj)
            % Get cell array of names of model parameters, i.e. ``{'Wu', 'Wz', 'M', 'N', 's'}``.
            %
            % See :meth:`mp.form.model_params`.

           params = {'Qu', 'Qz', 'M', 'N', 's'};
        end

        function vtypes = model_vvars(obj)
            % Get cell array of names of voltage state variables, i.e. ``{'lnvm', 'va'}``.
            %
            % See :meth:`mp.form.model_vvars`.

            vtypes = {'va', 'lnvm'};
        end

        function vtypes = model_zvars(obj)
            % Get cell array of names of non-voltage state variables, i.e. ``{'zr', 'zr'}``.
            %
            % See :meth:`mp.form.model_zvars`.

            vtypes = {'zr', 'zi'};
        end
    end     %% methods
end         %% classdef
