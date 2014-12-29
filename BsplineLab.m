% -----------------------------------------------------------------------
%
% B-spline Lab
% P.J. Barendrecht
%
% -----------------------------------------------------------------------

function BsplineLab(varargin)

if nargin == 0
  clear all
end

close all
clc

fprintf(['--------------------\n', ...
  'B-spline Lab\nPress ''H'' for help.\n', ...
  '--------------------\n\n'])

H = figure('Position', [150 100 1100 800],'NumberTitle','off','Name','B-spline Lab [Vertex]');
set(H, 'MenuBar','none');
set(H,'ToolBar','none');
set(H, 'Renderer','OpenGL');
hold on

set(H,'WindowButtonMotionFcn',@MouseMove);
set(H,'WindowButtonDownFcn',@MouseClick);
set(H,'WindowScrollWheelFcn',@MouseScroll);
set(H,'KeyPressFcn',@KeyPress )

global Degree Xi Points Weights

if nargin == 0
  % Default
  Degree = 3;
  Xi = [0,1,1,1,2,3,3,3,4];
elseif nargin == 2
  Degree = varargin{1};
  Xi = varargin{2};
else
  error('Wrong number of arguments.')
end

Points = [];
Weights = [];
MinX = 1;
MaxX = 9;
MinY = 1;
MaxY = 9;
axis( [ MinX-1, MaxX+1, MinY-1, MaxY+1 ] )
axis off

global SelectMode Grabbed Snap Labels Joints Index

SelectMode = 1;
Grabbed = 0;
Snap = 0;
Labels = 1;
Joints = 0;
Index = [];

SetupCurve()

% -----------------------------------------------------------------------

function SetupCurve()
global Basis Degree Xi Points Weights Index

[~, Basis] = CoxdeBoor(Degree,Xi,Degree+1,length(Xi)-Degree);
NumberPoints = length(Xi) - Degree - 1;

if size(Points,1) < NumberPoints
  AddPoints = NumberPoints - size(Points,1);
  disp(['Please add ', num2str(AddPoints), ' control points.'])
  set(gcf,'Pointer','FullCrossHair')
elseif size(Points,1) > NumberPoints
  DeletePoints = size(Points,1) - NumberPoints;
  Points = Points(1:end-DeletePoints,:);
  Index = [];
  Weights = Weights(1:end-DeletePoints);
  disp([':: Deleted ', num2str(DeletePoints), ' control points.'])
  UpdateCurve()
else
  UpdateCurve()
end

% -----------------------------------------------------------------------

function UpdateCurve()
global Basis Degree Xi Points Weights Index Labels Joints SelectMode

cla

% Control polygon
plot( Points(:,1), Points(:,2), '-', 'Color', .6*[1 1 1], 'LineWidth', 1 );

% Control points
plot( Points(:,1), Points(:,2), 'ko', 'MarkerSize', 12, 'MarkerFaceColor', [0 0 0] );

switch Labels
  case 1
    % Control point labels
    for P = 1:length(Points)
      text( Points(P,1), Points(P,2), ['    P_{', num2str(P), '}'], 'FontSize', 12 )
    end
  case 2
    % Weights
    for P = 1:length(Points)
      text( Points(P,1), Points(P,2), ['    <', num2str(Weights(P)), '>'], 'FontSize', 12 )
    end
  case 3
    % Blossom labels
    for P = 1:length(Points)
      Label = strtrim( sprintf('%.4g  ', Xi(P+1:P+Degree)) );
      text( Points(P,1), Points(P,2), ['    (', Label, ')' ], 'FontSize', 12 )
    end
  case 4
    % Knot spans/intervals
    if mod(Degree,2) == 0
      % Knot spans associated with control points
      for P = 1:length(Points)
        text( Points(P,1), Points(P,2), ['    [', num2str( Xi(round(Degree/2)+1+P)-Xi(round(Degree/2)+P) ), ']'], 'FontSize', 12 )
      end
    else
      % Knot spans associated with edges
      for E = 1:length(Points)-1
        text( (Points(E,1)+Points(E+1,1))/2, (Points(E,2)+Points(E+1,2))/2, ['[', num2str( Xi(round(Degree/2)+1+E)-Xi(round(Degree/2)+E) ), ']'], 'BackGroundColor', .8*[1 1 1], 'FontSize', 12 )
      end
    end
end

Denominator = Weights * Basis;
X = (Weights.*Points(:,1)' * Basis) ./ Denominator;
Y = (Weights.*Points(:,2)' * Basis) ./ Denominator;

plot( X, Y, 'Color', [255 128 0]/255, 'LineWidth', 3 );

if Joints
  Knots = unique(Xi((Degree+1):(end-Degree)));
  Pts = [];
  for Val = Knots
    Joint = ExactEval(Val,0);
    plot( Joint(1), Joint(2), 's', 'Color', [255 128 0]/255, 'MarkerSize', 10, 'MarkerFaceColor', [255 128 0]/255 )
    text( Joint(1), Joint(2)-.3, ['t = ', num2str(Val)], 'FontSize', 12 )
  end
end

if SelectMode == 1
  plot( Points(Index,1), Points(Index,2), 'o', 'MarkerSize', 16, 'Color', .9*[1 1 1], 'LineWidth', 3 );
else
  plot( [Points(Index,1) Points(Index+1,1)], [Points(Index,2) Points(Index+1,2)], 'Color', .9*[1 1 1], 'LineWidth', 3 );
end

% -----------------------------------------------------------------------

function MouseMove(~,~)
global Grabbed MoveAlong Snap Points Index

if Grabbed
  CurrentPoint = get(gca,'CurrentPoint');
  Pt = CurrentPoint(1,1:2);
  
  if Snap
    % Snap to grid
    Steps = 4;
    Pt = round(Steps*Pt) / Steps;
  end
  
  switch MoveAlong
    case 0
      Points(Index,:) = Pt;
    case 1
      Points(Index,1) = Pt(1);
    case 2
      Points(Index,2) = Pt(2);
  end
  
  UpdateCurve()
end

% -----------------------------------------------------------------------

function MouseClick(~,~)
global SelectMode Grabbed Snap Points Weights Xi Degree Index

% Left click
if strcmpi(get(gcf,'SelectionType'), 'Normal')
  
  if Grabbed
    Grabbed = 0;
    Index = [];
    disp(':: Releasing control point')
    UpdateCurve();
  else
    if strcmpi( get(gcf,'Pointer'),'FullCrossHair')
      
      CurrentPoint = get(gca,'CurrentPoint');
      Pt = CurrentPoint(1,1:2);
      
      if Snap
        % Snap to grid
        Steps = 4;
        Pt = round(Steps*Pt) / Steps;
      end
      
      Points(end+1,:) = Pt;
      Weights(end+1) = 1;
      plot( Points(end,1), Points(end,2), 'ko', 'MarkerSize', 12, 'MarkerFaceColor', [0 0 0] );
      disp(':: Added a new control point')
      
      if size(Points,1) == (length(Xi)-Degree-1)
        set(gcf,'Pointer','Cross')
        UpdateCurve()
      end
      
    end
  end
  
  % Right click
elseif strcmpi(get(gcf,'SelectionType'), 'Alt')
  
  if SelectMode == 1
    disp(':: Selecting control point')
    Index = FindClosest();
    UpdateCurve()
    
  elseif SelectMode == 2
    disp(':: Selecting edge')
    Index = FindClosest();
    UpdateCurve()
  end
  
end

% -----------------------------------------------------------------------

function MouseScroll(~,Event)
global SelectMode Weights Index

if SelectMode == 1
  Scale = 10;
  Value = Event.VerticalScrollCount / Scale;
  Weights(Index) = Weights(Index) + Value;
  if Weights(Index) < 0
    warning('Weight value negative!')
  end
  UpdateCurve()
end

% -----------------------------------------------------------------------

function KeyPress(~,Event)
global Xi Degree Labels Joints Points Weights Index SelectMode Grabbed Snap MoveAlong

switch Event.Character
  case {'a', 'A'}
    % Toggle axis visibility
    if strcmpi(get(gca, 'Visible'), 'On')
      axis off;
    else
      axis on;
    end
  case {'b', 'B'}
    % Plot the basis
    PlotBasis();
  case {'c', 'C'}
    % Show Bezier extraction matrices [TODO]
  case {'d', 'D'}
    % Change the degree
    Deg = input('Degree (e.g. quadratic = 2): ');
    if length(Xi) < (2 + 2*Deg)
      warning('Knot vector does not satisfy minimal length for this degree!')
    else
      Degree = Deg;
      SetupCurve()
    end
  case {'e', 'E'}
    % Set select mode to Edges
    SelectMode = 2;
    Index = [];
    set(gcf,'Name','B-spline Lab [Edge]');
    disp(':: Select Mode [Edge]')
    UpdateCurve();
  case {'f', 'F'}
    % Print graphical overview of the support of the basis functions
    PrintSupports()
  case {'g', 'G'}
    % Grab control point
    if (SelectMode == 1 && length(Index))
      disp(':: Grabbing control point')
      Grabbed = 1;
      MoveAlong = 0;
    end
  case {'h', 'H'}
    % Plot hodograph of the curve
    Hodograph()
  case {'i', 'I'}
    % Display information about the curve
    Knts = strtrim(sprintf('%.4g  ', Xi));
    Wts = strtrim(sprintf('%.4g  ', Weights));
    fprintf(['\nKnotvector: [', Knts, ']\n', ...
      'Degree: ', num2str(Degree), '\n', ...
      'Number of control points: ', num2str(size(Points,1)), '\n', ...
      'Weights: ', Wts, '\n\n'])
  case {'j', 'J'}
    % Toggle visibility of the joints (images of the knots)
    Joints = ~Joints;
    UpdateCurve()
  case {'k', 'K'}
    % Change the knotvector
    Xi = input('Knotvector (between [ and ]): ');
    SetupCurve()
  case {'l', 'L'}
    % Cycle through label types (0 no labels, 1 control point labels, 2 weights, 3 blossom labels, 4 knot spans)
    Labels = mod(Labels+1,5);
    switch Labels
      case 1
        disp(':: Displaying control point labels')
      case 2
        disp(':: Displaying weights')
      case 3
        disp(':: Displaying blossom (polar form) labels')
      case 4
        disp(':: Displaying knot spans (intervals)')
    end
    UpdateCurve()
  case {'m', 'M'}
    % Global midpoint knot refinement (Lane-Riesenfeld) [TODO]
  case {'n', 'N'}
    % Placeholder
  case {'o', 'O'}
    % Order elevation [TODO]
  case {'p', 'P'}
    % Set exact position of selected control point
    if SelectMode == 1
      Points(Index,:) = input('Position (between [ and ]): ');
      UpdateCurve()
    end
  case {'q', 'Q'}
    % Quit
    close all
    clear all
    clc
  case {'r', 'R'}
    % Local refinement (knot insertion, Boehm) [TODO]
  case {'s', 'S'}
    % Save the curve in SVG format using plot2svg
    % See http://www.mathworks.nl/matlabcentral/fileexchange/7401-scalable-vector-graphics-svg-export-of-figures
    % Also http://www.mathworks.nl/matlabcentral/fileexchange/23629-exportfig
    plot2svg
  case {'t', 'T'}
    % Evaluate curve at parameter value t (scaffolding approach)
    Scaffold = strcmpi(get(gcf,'CurrentModifier'), 'Shift');
    Val = input('Evaluate curve at parameter value: ');
    if Val < Xi(Degree+1) || Val >= Xi(end-Degree)
      warning('Parameter value out of range!')
    end
    Pt = ExactEval(Val,Scaffold)
    plot(Pt(1), Pt(2), 'rx', 'MarkerSize', 12, 'LineWidth', 2)
    text(Pt(1), Pt(2)-.3, ['t = ', num2str(Val)], 'Color', [1 0 0], 'FontSize', 12 )
  case {'u', 'U'}
    % Reset weights to unity
    Weights = ones(1,length(Weights));
    UpdateCurve()
  case {'v', 'V'}
    % Set select mode to Vertices
    SelectMode = 1;
    Index = [];
    set(gcf,'Name','B-spline Lab [Vertex]');
    disp(':: Select Mode [Vertex]')
    UpdateCurve();
  case {'w', 'W'}
    % Set exact weight of selected control point
    if SelectMode == 1
      Weights(Index) = input('Weight: ');
      UpdateCurve();
    end
  case {'x', 'X'}
    % Toggle constrained movement in X direction
    if Grabbed
      MoveAlong = ~MoveAlong;
    end
  case {'y', 'Y'}
    % Toggle constrained movement in Y direction
    if Grabbed
      MoveAlong = 2*(~MoveAlong);
    end
  case {'z', 'Z'}
    % Show curve in homogeneous space (along with perspective projection) [TODO]
    
  case {'#'}
    % Toggle snap to "grid"
    Snap = ~Snap;
  case {'+'}
    % Increase knot-span associated with currently selected vertex or edge
    if (mod(Degree,2) == 0 && SelectMode == 1 && length(Index)) || (mod(Degree,2) == 1 && SelectMode == 2 && length(Index))
      Xi = [Xi(1:round(Degree/2)+Index) Xi(round(Degree/2)+1+Index:end) + .05];
      SetupCurve()
    end
  case {'-'}
    % Decrease knot-span
    if (mod(Degree,2) == 0 && SelectMode == 1 && length(Index)) || (mod(Degree,2) == 1 && SelectMode == 2 && length(Index))
      if Xi(round(Degree/2)+Index) <= Xi(round(Degree/2)+1+Index) - .05
        Xi = [Xi(1:round(Degree/2)+Index) Xi(round(Degree/2)+1+Index:end) - .05];
        SetupCurve()
      end
    end
  case {'?'}
    % Show help [TODO]
    disp('Help')
end

% -----------------------------------------------------------------------

function Pt = ExactEval(Val, Scaffolding)
global Xi Degree Points

% Using de Boor's algorithm for exact evaluation (generalised de Casteljau)
% This is for B-splines, not for NURBS. Use perspective projection? [TODO]

% Determine relevant segment
Upper = length(Xi) - Degree;

l = 1;
while Val >= Xi(l)
  l = l+1;
  % Check
  if l == Upper
    break
  end
end
l = l-1;

CurrentPoints = Points((l-Degree):l,:);

for k = 1:Degree
  NewPoints = zeros(size(CurrentPoints,1)-1,2);
  
  for n = 1:size(NewPoints,1)
    Alpha = (Val-Xi(l-Degree+k+(n-1)))/(Xi(l+n)-Xi(l-Degree+k+(n-1)));
    NewPoints(n,:) = (1-Alpha)*CurrentPoints(n,:) + Alpha*CurrentPoints(n+1,:);
  end
  
  if Scaffolding
    plot( NewPoints(:,1), NewPoints(:,2), 'x-', 'MarkerSize', 12, 'Color', .4*[1 1 1], 'LineWidth', 1);
  end
  
  CurrentPoints = NewPoints;
  
end

Pt = CurrentPoints;

% -----------------------------------------------------------------------

function Index = FindClosest()
global Points SelectMode;

CurrentPoint = get(gca,'CurrentPoint');
Pt = CurrentPoint(1,1:2);

if SelectMode == 1
  
  Distances = zeros(1,length(Points));
  
  for P = 1:length(Distances)
    Distances(P) = norm( Points(P,:)-Pt, 2);
  end
  
else
  
  Distances = zeros(1,length(Points)-1);
  Px = Pt(1);
  Py = Pt(2);
  
  for E = 1:length(Distances)
    Ax = Points(E,1);
    Ay = Points(E,2);
    Bx = Points(E+1,1);
    By = Points(E+1,2);
    % Difference
    Dx = Bx - Ax;
    Dy = By - Ay;
    Val = ((Px-Ax)*Dx + (Py-Ay)*Dy)/(Dx^2 + Dy^2);
    
    if Val < 0
      Val = 0;
    elseif Val > 1
      Val = 1;
    end
    
    Distances(E) = norm( (1-Val)*Points(E,:)+Val*Points(E+1,:)-Pt, 2);
  end
  
end

[~,Index] = min( Distances );

% -----------------------------------------------------------------------

function Hodograph()
global Points Xi Degree

% Consider weights [TODO]

NewPts = zeros(length(Points)-1,2);

for Q = 1:length(NewPts)
  % Do not use first and last knot
  NewPts(Q,:) = (Degree/(Xi(Q+1+Degree)-Xi(Q+1)))*(Points(Q+1,:)-Points(Q,:));
end

H = figure('NumberTitle','off','Name','Hodograph');
set(H, 'MenuBar','none');
set(H,'ToolBar','none');
set(H, 'Renderer','OpenGL');
hold on
set(H,'KeyPressFcn',@HodoKeyPress )

% Control polygon
plot( NewPts(:,1), NewPts(:,2), '-', 'Color', .6*[1 1 1], 'LineWidth', 1 );

% Control points
plot( NewPts(:,1), NewPts(:,2), 'ko', 'MarkerSize', 12, 'MarkerFaceColor', [0 0 0] );

[~, NewBasis] = CoxdeBoor(Degree-1,Xi(2:end-1),Degree,length(Xi)-Degree-1);

X = NewPts(:,1)' * NewBasis;
Y = NewPts(:,2)' * NewBasis;

plot( X, Y, 'Color', [68 170 0]/255, 'LineWidth', 3 );
plot( 0,0, '+', 'MarkerSize', 14, 'LineWidth', 2, 'Color', [68 170 0]/255 );

for Q = 1:length(NewPts)
  text( NewPts(Q,1), NewPts(Q,2), ['   Q_{', num2str(Q), '}'], 'FontSize', 12 )
end

axis tight
axis off
axis equal

% -----------------------------------------------------------------------

function HodoKeyPress(~,Event)
switch Event.Character
  case {'q', 'Q'}
    % Quit
    close
  case {'s', 'S'}
    % Save the curve in SVG format using plot2svg
    plot2svg
end

% -----------------------------------------------------------------------

function PrintSupports()
global Xi Degree

Spaces = 0;

for B = 1:length(Xi)-Degree-1
  Support = Xi(B:B+Degree+1);
  fprintf('\n')
  fprintf(repmat(' ',1,Spaces))
  fprintf('|')
  PrevKnot = Support(1);
  
  for K = 2:length(Support)
    CurrentKnot = Support(K);
    if CurrentKnot == PrevKnot
      fprintf('|')
    else
      fprintf('-----|')
    end
    PrevKnot = CurrentKnot;
  end
  
  if Xi(B+1) == Xi(B)
    Spaces = Spaces + 1;
  else
    Spaces = Spaces + 6;
  end
end

fprintf('\n\n')

% -----------------------------------------------------------------------

function PlotBasis()
global Xi Degree Points Weights

if strcmpi(get(gcf,'CurrentModifier'), 'Shift')
  From = 1;
  To = length(Xi);
else
  From = Degree+1;
  To = length(Xi)-Degree;
end

[Evals, Basis] = CoxdeBoor(Degree,Xi,From,To);

% Colours in hexadecimal format
Colours = ['AA0088'; 'D40000'; 'FF8000'; 'FFCC00'; '44AA00'; '2CA089'; '3771C8'];

for B = 1:length(Points)
  k = mod(B,length(Colours))+1;
  Col = hex2dec([Colours(k,1:2); Colours(k,3:4); Colours(k,5:6)])'/255;
  plot(Points(B,1), Points(B,2), 'o', 'MarkerSize', 16, 'Color', Col, 'LineWidth', 3 );
end

% Estimation for window dimensions
Width = 1000;
Height = Width / (length(Points) - Degree + 1) + 100;

H = figure('NumberTitle','off','Name','Basis');
set(H, 'MenuBar','none');
set(H,'ToolBar','none');
set(H,'Position',[700 200 Width Height]);
set(H, 'Renderer','OpenGL');
hold on
set(H,'KeyPressFcn',@HodoKeyPress )

Denominator = Weights * Basis;

for B = 1:length(Points)
  k = mod(B,length(Colours))+1;
  Col = hex2dec([Colours(k,1:2); Colours(k,3:4); Colours(k,5:6)])'/255;
  %plot(Evals, Basis(B,:), 'Color', Col, 'LineWidth', 2)
  plot(Evals, Weights(B)*Basis(B,:) ./ Denominator, 'Color', Col, 'LineWidth', 2)
end

axis tight

% -----------------------------------------------------------------------
