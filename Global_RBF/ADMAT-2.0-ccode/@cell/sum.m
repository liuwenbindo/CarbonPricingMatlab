function sout=sum(s1, a)

%       ******************************************************************
%       *                          ADMAT - 2.0                           *
%       *              Copyright (c) 2008-2009 Cayuga Research           *
%       *                Associates, LLC. All Rights Reserved.           *
%       ******************************************************************


global ADhess;

if nargin  == 1
    a = 1;
end

[m,n]=size(s1);

if ADhess
    [p,q]=size(s1{1});
    if p==1 || q==1         % each one in the cell is a vector or scalar
        sout = zeros(m,n);
        for i=1:m
            for j=1:n
                sout(i,j)=sum(s1{i,j});
            end
        end
    else
        sout = cell(m,n);
        if a == 1
            for i=1:m
                for j=1:n
                    sout{i,j}=sum(s1{i,j})';
                end
            end
        else    % a == 2
            for i=1:m
                for j=1:n
                    sout{i,j}=sum(s1{i,j});
                end
            end
        end
        
    end
    
else   % Jacobian matrix
    [p,q]=size(s1{1});
    if p==1 || q==1         % each one in the cell is a vector or scalar
        sout = zeros(max(m,n),1);
        for i=1:max(m,n)
                sout(i)=sum(s1{i});
        end
    else
        if a == 1  % sum by rows
            sout = zeros(q, max(m,n));
            for i=1:max(m,n)
                    sout(:,i)=sum(s1{i})';
            end
        else    % a == 2 % sum by columns
            sout = zeros(p, max(m,n));
            for i=1:max(m,n)
                    sout(:,i)=sum(s1{i},2);
            end
        end
        
    end
%     sout = zeros(max(m,n),1);
%     for i=1:max(m,n)
%         sout(:,i)=sum(s1{i,j})';
%     end
    
end
