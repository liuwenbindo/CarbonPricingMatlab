function y=flipud(x)
%       ******************************************************************
%       *                          ADMAT - 2.0                           *
%       *              Copyright (c) 2008-2009 Cayuga Research           *
%       *                Associates, LLC. All Rights Reserved.           *
%       ******************************************************************


y.val=flipud(x.val);
y.derivH=flipud(x.derivH);

y=class(y,'derivH');
