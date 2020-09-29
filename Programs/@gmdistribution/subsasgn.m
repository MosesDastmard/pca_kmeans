function [varargout] = subsasgn(varargin)
%SUBSASGN Subscripted reference for a gmdistribution object.
%   Subscript assignment is not allowed for a gmdistribution object.

%   Copyright 2007 The MathWorks, Inc.
%   $Revision: 1.1.6.2 $  $Date: 2008/02/29 13:12:35 $

error('stats:gmdistribution:subsasgn:NotAllowed',...
      'Subscripted assignments are not allowed.')