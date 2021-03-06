function adjmax(i)
%
%  04/2007 -- rearrage the program for readibilty
%  04/2007 -- correct the computation of derivative
%   01/2010 -- add nondifferentiable points detecting.
%
%       ******************************************************************
%       *                          ADMAT - 2.0                           *
%       *              Copyright (c) 2008-2009 Cayuga Research           *
%       *                Associates, LLC. All Rights Reserved.           *
%       ******************************************************************


global tape;
global globp;


if isempty(tape(i).arg2vc)
    y = tape(tape(i).arg1vc).val;
    [m, n] = size(y);
    [temp,I]=max(y);
    if m == 1 || n == 1
        tape(tape(i).arg1vc).W(I,:)=tape(tape(i).arg1vc).W(I,:)+tape(i).W;
    else                 % y is a matrix
        for j = 1 : globp
            tape(tape(i).arg1vc).W(I(j),:,j)=tape(tape(i).arg1vc).W(I(j),:,j)+...
                tape(i).W(I(j),:);
        end
    end

elseif (size(tape(i).val,1)==1  || size(tape(i).val,2)==1)     % output is a vector
    x = tape(tape(i).arg2vc).val;
    y = tape(tape(i).arg1vc).val;

    I1=find(y - x > 0);
    I2=find(y - x== 0);

    I3=find(y - x < 0);
    if length(y)==1               % y is a scalar
        if ~isempty(I1)
            tape(tape(i).arg1vc).W=tape(tape(i).arg1vc).W+sum(tape(i).W(I1,:));
        end
        if ~isempty(I3)
            tape(tape(i).arg2vc).W(I3,:)=tape(tape(i).arg2vc).W(I3,:)+tape(i).W(I3,:);
        end
        if ~isempty(I2)
            if isequal(tape(tape(i).arg1vc).W, sum(tape(i).W(I2,:))./2) && ...
                    isequal(tape(tape(i).arg2vc).W(I2,:),tape(i).W(I2,:)/2)

                tape(tape(i).arg1vc).W=tape(tape(i).arg1vc).W+sum(tape(i).W(I2,:))./2;
                tape(tape(i).arg2vc).W(I2,:)=tape(tape(i).arg2vc).W(I2,:)+tape(i).W(I2,:)/2;
            else
                error('Nondifferentiable points in max() was detected.');
            end

        end
    elseif length(x)==1           % x is a scalar
        if ~isempty(I3)
            tape(tape(i).arg2vc).W=tape(tape(i).arg2vc).W+sum(tape(i).W(I3,:));
        end
        if ~isempty(I1)
            tape(tape(i).arg1vc).W(I1,:)=tape(tape(i).arg1vc).W(I1,:)+tape(i).W(I1,:);
        end
        if ~isempty(I2)
            if isequal(tape(tape(i).arg2vc).W, sum(tape(i).W(I2,:))/2) && ...
                    isequal(tape(tape(i).arg1vc).W(I2,:), tape(i).W(I2,:)/2)

                tape(tape(i).arg2vc).W=tape(tape(i).arg2vc).W+sum(tape(i).W(I2,:))/2;
                tape(tape(i).arg1vc).W(I2,:)=tape(tape(i).arg1vc).W(I2,:)+tape(i).W(I2,:)/2;
            else
                error('Nondifferentiable points in max() was detected.');
            end

        end
    else                         % both x and y are vectors
        if ~isempty(I3)
            tape(tape(i).arg2vc).W(I3,:)=tape(tape(i).arg2vc).W(I3,:)+tape(i).W(I3,:);
        end
        if ~isempty(I1)
            tape(tape(i).arg1vc).W(I1,:)=tape(tape(i).arg1vc).W(I1,:)+tape(i).W(I1,:);
        end
        if ~isempty(I2)
            if isequal(tape(tape(i).arg2vc).W(I2,:),tape(i).W(I2,:)/2) && ...
                    isequal(tape(tape(i).arg1vc).W(I2,:),tape(i).W(I2,:)/2)

                tape(tape(i).arg2vc).W(I2,:)=tape(tape(i).arg2vc).W(I2,:)+tape(i).W(I2,:)/2;
                tape(tape(i).arg1vc).W(I2,:)=tape(tape(i).arg1vc).W(I2,:)+tape(i).W(I2,:)/2;
            else
                error('Nondifferentiable points in max() was detected.');
            end

        end
    end

else                                   % tape(i).val is a matrix
    x = tape(tape(i).arg2vc).val;
    y = tape(tape(i).arg1vc).val;
    [I1,I1j]=find(y - x> 0);
    [I2,I2j]=find(y - x == 0);
    % non differentiable points checking
    if ~isempty(I2) || ~isempty(I2j)
        error('Nondifferentiable points in max() was detected.');
    end
    [I3,I3j]=find(y - x < 0);
    if length(y)==1                    % y is a scalar
        if ~isempty(I1)
            for j=1:length(I1)
                tape(tape(i).arg1vc).W=tape(tape(i).arg1vc).W+...
                    sum(tape(i).W(I1(j),I1j(j),:));
            end
        end
        if ~isempty(I3)
            for j=1:length(I3)
                tape(tape(i).arg2vc).W(I3(j),I3j(j),:) = tape(tape(i).arg2vc).W(I3(j),I3j(j),:)+ ...
                    tape(i).W(I3(j),I3j(j),:);
            end
        end
        if ~isempty(I2)
            for j = 1 : lenth(I2)
                if isequal(tape(tape(i).arg1vc).W ,sum(tape(i).W(I2(j),I2j(j),:))./2) && ...
                        isequal(tape(tape(i).arg2vc).W(I2(j),I2j(j),:), tape(i).W(I2(j),I2j(j),:)/2)

                    tape(tape(i).arg1vc).W = tape(tape(i).arg1vc).W +...
                        sum(tape(i).W(I2(j),I2j(j),:))./2;
                    tape(tape(i).arg2vc).W(I2(j),I2j(j),:)=tape(tape(i).arg2vc).W(I2(j),I2j(j),:) + ...
                        tape(i).W(I2(j),I2j(j),:)/2;
                else
                    error('Nondifferentiable points in max() was detected.');
                end

            end
        end
    elseif length(x)==1                 % x is a scalar
        if ~isempty(I3)
            for j=1:length(I3)
                tape(tape(i).arg2vc).W=tape(tape(i).arg2vc).W+ ...
                    sum(tape(i).W(I3(j),I3j(j),:));
            end
        end
        if ~isempty(I1)
            for j=1:length(I1)
                tape(tape(i).arg1vc).W(I1(j),I1j(j),:)=tape(tape(i).arg1vc).W(I1(j),I1j(j),:)+...
                    tape(i).W(I1(j),I1j(j),:);
            end
        end
        if ~isempty(I2)
            for j=1:length(I2)
                if isequal(tape(tape(i).arg2vc).W, sum(tape(i).W(I2(j),I2j(j),:))/2) && ...
                        isequal(tape(tape(i).arg1vc).W(I2(j),I2j(j),:),tape(i).W(I2(j),I2j(j),:)/2)

                    tape(tape(i).arg2vc).W=tape(tape(i).arg2vc).W+sum(tape(i).W(I2(j),I2j(j),:))/2;
                    tape(tape(i).arg1vc).W(I2(j),I2j(j),:)=tape(tape(i).arg1vc).W(I2(j),I2j(j),:)+...
                        tape(i).W(I2(j),I2j(j),:)/2;
                else
                    error('Nondifferentiable points in max() was detected.');
                end

            end
        end
    else                                    % both x and y are matrices
        if ~isempty(I3)
            for j=1:length(I3)
                tape(tape(i).arg2vc).W(I3(j),I3j(j),:)=tape(tape(i).arg2vc).W(I3(j),I3j(j),:)+...
                    tape(i).W(I3(j),I3j(j),:);
            end
        end
        if ~isempty(I1)
            for j=1:length(I1)
                tape(tape(i).arg1vc).W(I1(j),I1j(j),:)=tape(tape(i).arg1vc).W(I1(j),I1j(j),:)+...
                    tape(i).W(I1(j),I1j(j),:);
            end
        end
        if ~isempty(I2)
            for j = 1 : length(I2)
                if isequal(tape(tape(i).arg2vc).W(I2(j),I2j(j),:), tape(i).W(I2(j),I2j(j),:)/2) && ...
                        isequal(tape(tape(i).arg1vc).W(I2(j),I2j(j),:),tape(i).W(I2(j),I2j(j),:)/2)

                    tape(tape(i).arg2vc).W(I2(j),I2j(j),:) = tape(tape(i).arg2vc).W(I2(j),I2j(j),:) + ...
                        tape(i).W(I2(j),I2j(j),:)/2;
                    tape(tape(i).arg1vc).W(I2(j),I2j(j),:) = tape(tape(i).arg1vc).W(I2(j),I2j(j),:) + ...
                        tape(i).W(I2(j),I2j(j),:)/2;
                else
                    error('Nondifferentiable points in max() was detected.');
                end

            end
        end
    end
end

