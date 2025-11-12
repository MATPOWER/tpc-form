function varargout = run_orpd(varargin)
% run_orpd  Run an optimal reactive power dispatch (ORPD) using a quadratic
%           network formulation in tranformed polar coordinates (QTPC).
% ::
%
%   run_orpd(d, mpopt, 'mpx', wgv.xt_tpc_3p_orpd)
%   task = run_orpd(...)
%
% This is the main function used to run optimal reactive power dipatch (ORPD)
% problems via the **flexible** |*MATPOWER*| **framework**.
%
% This function is a simple wrapper around run_mp, calling it
% with the first argument set to ``@wgv.task_orpd``.
%
% Inputs:
%   d : data source specification, currently assumed to be a |MATPOWER|
%       case name or case struct (``mpc``)
%   mpopt (struct) : *(optional)* |MATPOWER| options struct
%
%       Additional optional inputs can be provided as *<name>, <val>* pairs,
%       with the following options:
%
%       - ``'print_fname'`` - file name for saving pretty-printed output
%       - ``'soln_fname'`` - file name for saving solved case
%       - ``'mpx'`` - |MATPOWER| extension or cell array of |MATPOWER|
%         extensions to apply
%
% Output:
%   task (wgv.task_orpd) : task object containing the solved run including the
%       data, network, and mathematical model objects.
%
% See also run_mp, wgv.task_orpd.

%   MATPOWER
%   Copyright (c) 2021-2025, Power Systems Engineering Research Center (PSERC)
%   by Wilson Gonzalez Vanegas, Universidad Nacional de Colombia Sede Manizales
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://matpower.org for more info.

[varargout{1:nargout}] = run_mp(@wgv.task_orpd, varargin{:});
