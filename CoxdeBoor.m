% ----------------------------------------------------------------------

% B-spline Lab is free software: you can redistribute it and/or modify it under the terms of the 
% GNU Lesser General Public License as published by the Free Software Foundation, either version 3 
% of the License, or (at your option) any later version.
% 
% B-spline Lab is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without 
% even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU 
% Lesser General Public License for more details.
% 
% You should have received a copy of the GNU Lesser General Public License along with B-spline Lab. 
% If not, see http://www.gnu.org/licenses/.

% -----------------------------------------------------------------------

function [Evals, Basis] = CoxdeBoor(Degree,Xi,From,To)

Knots = unique(Xi(From:To));
Sub = 30;

% One overlap (Sub-1) and use the knot-spans
Evals = zeros((Sub-1)*(length(Knots)-1)+1,1);
for j = 1:length(Knots)-1
  if j == 1
    Evals(1:Sub) = linspace(Knots(j),Knots(j+1),Sub);
  else
    Evals((Sub-1)*(j-1)+1:(Sub-1)*j+1) = linspace(Knots(j),Knots(j+1),Sub);
  end  
end

% B-splines are defined on half-open knot intervals, so the last knot is not included.
Evals(end) = Evals(end) - 1e-4;

NumberBasis = length(Xi) - Degree - 1;
k = NumberBasis + Degree;
m = length(Evals);
TempBasis = zeros(k,m);

for B = 1:k
  % All values in Evals
  for n = 1:m
    t = Evals(n);
    if Xi(B) <= t && t < Xi(B+1)
      TempBasis(B,n) = 1;
    end
  end
end

Basis = Recurse(Degree,Xi,Evals,1,TempBasis);

% -----------------------------------------------------------------------

function CurrentTempBasis = Recurse(Degree,Xi,Evals,Current,PrevTempBasis)
k = size(PrevTempBasis,1)-1;
m = size(PrevTempBasis,2);
CurrentTempBasis = zeros(k,m);

for B = 1:k
  for n = 1:m
    t = Evals(n);
    
    if Xi(B+Current) - Xi(B) ~= 0
      CurrentTempBasis(B,n) = ((t-Xi(B))/(Xi(B+Current)-Xi(B))) * PrevTempBasis(B,n);
    end
    
    if Xi(B+Current+1) - Xi(B+1) ~= 0
      CurrentTempBasis(B,n) = CurrentTempBasis(B,n) + ((Xi(B+Current+1)-t)/(Xi(B+Current+1)-Xi(B+1))) * PrevTempBasis(B+1,n);
    end
  end
end

if Current ~= Degree
  CurrentTempBasis = Recurse(Degree,Xi,Evals,Current+1,CurrentTempBasis);
end

% -----------------------------------------------------------------------
