function [Rec3dProject] = tom_Rec3dAlignmentEngine(Rec3dControl, Rec3dProject)
%TOM_REC3DALIGNMENTENGINE is a module of TOM_REC3DGUI.
%
%   Rec3dProject = tom_Rec3dAlignmentEngine(Rec3dControl,Rec3dProject)
%
%   It calculates reconstruction parameter for both singleaxis and dualaxis
%   reconstruction projects.
%
%PARAMETERS
%
%  INPUT
%   Rec3dControl     control-structure of Rec3d
%   Rec3dProject     main-structure of Rec3d
%  
%  OUTPUT
%   Rec3dProject     aktualized main-structure of tom_Rec3dGUI
%
%EXAMPLE
%
%REFERENCES
%
%SEE ALSO
%   TOM_REC3DEULERANGLES   TOM_REC3DEULERANGLESRESIDUAL
%   TOM_REC3DNEWORIGIN2
%   TOM_REC3DRIGIDBODYALIGNMENT   TOM_REC3DRIGIDBODYRESIDUAL
%   TOM_REC3DALIGNMENTRESIDUAL
%
%   01/01/07 ME
%
%   Nickell et al., 'TOM software toolbox: acquisition and analysis for electron tomography',
%   Journal of Structural Biology, 149 (2005), 227-234.
%
%   Copyright (c) 2004-2007
%   TOM toolbox for Electron Tomography
%   Max-Planck-Institute of Biochemistry
%   Dept. Molecular Structural Biology
%   82152 Martinsried, Germany
%   http://www.biochem.mpg.de/tom


% -------------------------------------------------------------------------
% Rec3dControl -> AlignmentEngine
% -------------------------------------------------------------------------
ProtocolMode     = Rec3dControl.ProtocolMode;
% -------------------------------------------------------------------------
% Rec3dProject -> AlignmentEngine
% -------------------------------------------------------------------------
% PROJECT
TiltingGeometry  = Rec3dProject.PROJECT.TiltingGeometry;
NumOfProj1       = Rec3dProject.PROJECT.NumOfProj1;
RefProj1         = Rec3dProject.PROJECT.RefProj1;
Tiltangles1      = Rec3dProject.PROJECT.Tiltangles1;
MarkersOnProj1   = Rec3dProject.PROJECT.MarkersOnProj1;
Markerfile1      = Rec3dProject.PROJECT.Markerfile1;
NumOfProj2       = Rec3dProject.PROJECT.NumOfProj2;
RefProj2         = Rec3dProject.PROJECT.RefProj2;
Tiltangles2      = Rec3dProject.PROJECT.Tiltangles2;
MarkersOnProj2   = Rec3dProject.PROJECT.MarkersOnProj2;
Markerfile2      = Rec3dProject.PROJECT.Markerfile2;
NumOfMarkers     = Rec3dProject.PROJECT.NumOfMarkers;
Imdim            = Rec3dProject.PROJECT.Imdim;
% ALIGNMENT
AlignmentMethod  = Rec3dProject.ALIGNMENT.AlignmentMethod;
ReferenceMarker  = Rec3dProject.ALIGNMENT.ReferenceMarker;
% -------------------------------------------------------------------------


% -------------------------------------------------------------------------
% Check
% -------------------------------------------------------------------------
% Check Origin
if strcmp(TiltingGeometry, 'singleaxis')
    if (Markerfile1.Value(2, RefProj1, ReferenceMarker)) == (-1) || ...
       (Markerfile1.Value(3, RefProj1, ReferenceMarker)) == (-1)
         % ERROR
         msgbox('Alignment Origin is not defined!', 'Do Alignment', 'error');
         return;
    end
end   
if strcmp(TiltingGeometry, 'dualaxis')
    if (Markerfile1.Value(2, RefProj1, ReferenceMarker)) == (-1) || ...
       (Markerfile1.Value(3, RefProj1, ReferenceMarker)) == (-1) || ...
       (Markerfile2.Value(2, RefProj2, ReferenceMarker)) == (-1) || ...
       (Markerfile2.Value(3, RefProj2, ReferenceMarker)) == (-1)
         % ERROR
         msgbox('Alignment Origin is not defined!', 'Do Alignment', 'error');
         return;
    end
end    
% -------------------------------------------------------------------------


% ProtocolMode
if strcmp(ProtocolMode, 'on')
    disp(' ');
    disp('Do Alignment');
    disp(' ');
end


% -------------------------------------------------------------------------
% ALIGNMENT
% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
% TiltingGeometry 'singleaxis'
% -------------------------------------------------------------------------
if strcmp(TiltingGeometry, 'singleaxis')

% -------------------------------------------------------------------------
% AlignmentMethod 'rigidbody'
% -------------------------------------------------------------------------
if strcmp(AlignmentMethod, 'rigidbody')

% Define alignment origin of tiltseries1
Origin1 = [Markerfile1.Value(2, RefProj1, ReferenceMarker) ...
           Markerfile1.Value(3, RefProj1, ReferenceMarker) ...
           Imdim/2+1];

% Calculate 'rigidbody' alignment of tiltseries1
[Markerfile1.Value, Rho1, Sigma1, x1, y1, z1] = tom_Rec3dRigidBodyAlignment...
(Markerfile1.Value, ReferenceMarker, RefProj1, Origin1, Imdim, ProtocolMode);

% Combine 3d marker coordinates of tiltseries1
m3d1(1,1:NumOfMarkers) = x1(1:NumOfMarkers);
m3d1(2,1:NumOfMarkers) = y1(1:NumOfMarkers);
m3d1(3,1:NumOfMarkers) = z1(1:NumOfMarkers);

% Move tiltseries1 3d marker coordinates of ReferenceMarker to [0 0 0]
m3d1(1,1:NumOfMarkers) = m3d1(1,1:NumOfMarkers) - Origin1(1) + (Imdim/2+1); 
m3d1(2,1:NumOfMarkers) = m3d1(2,1:NumOfMarkers) - Origin1(2) + (Imdim/2+1); 
m3d1(3,1:NumOfMarkers) = m3d1(3,1:NumOfMarkers) - Origin1(3) + (Imdim/2+1);

% Tiltaxis, translations, isoscale of tiltseries1
for k = 1:NumOfProj1
    Tiltaxis1(k) = (180/pi)*Rho1;
    tx1(k) = Markerfile1.Value(7, k, 1);
    ty1(k) = Markerfile1.Value(8, k, 1);
    isoscale1(k) = 1;
end

% Calculate ResidualMatrix and WarpAlignment of tiltseries1
[ResidualMatrix1, WarpAlignment1] = tom_Rec3dRigidBodyResidual...
(MarkersOnProj1,NumOfProj1,NumOfMarkers,Origin1,Imdim,Tiltaxis1,Tiltangles1,m3d1,tx1,ty1);

end
% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
% AlignmentMethod 'freetilt'
% -------------------------------------------------------------------------
if strcmp(AlignmentMethod, 'freetilt')
    % ERROR
    msgbox('Alignment Method freetilt not implemented!', ...
                   'Do Alignment', 'error');
    return;
end
% -------------------------------------------------------------------------

end
% -------------------------------------------------------------------------


% -------------------------------------------------------------------------
% TiltingGeometry 'dualaxis'
% -------------------------------------------------------------------------
if strcmp(TiltingGeometry, 'dualaxis')

% -------------------------------------------------------------------------
% AlignmentMethod 'rigidbody'
% -------------------------------------------------------------------------
if strcmp(AlignmentMethod, 'rigidbody')

% Define alignment origin of tiltseries1
Origin1 = [Markerfile1.Value(2, RefProj1, ReferenceMarker) ...
           Markerfile1.Value(3, RefProj1, ReferenceMarker) ...
           Imdim/2+1];

% Define alignment origin of tiltseries2
Origin2 = [Markerfile2.Value(2, RefProj2, ReferenceMarker) ...
           Markerfile2.Value(3, RefProj2, ReferenceMarker) ...
           Imdim/2+1];

% Calculate 'rigidbody' alignment of tiltseries1
[Markerfile1.Value, Rho1, Sigma1, x1, y1, z1] = tom_Rec3dRigidBodyAlignment...
(Markerfile1.Value, ReferenceMarker, RefProj1, Origin1, Imdim, ProtocolMode);

% Calculate 'rigidbody' alignment of tiltseries2
[Markerfile2.Value, Rho2, Sigma2, x2, y2, z2] = tom_Rec3dRigidBodyAlignment...
(Markerfile2.Value, ReferenceMarker, RefProj2, Origin2, Imdim, ProtocolMode);

% Combine 3d marker coordinates of tiltseries1
m3d1(1,1:NumOfMarkers) = x1(1:NumOfMarkers);
m3d1(2,1:NumOfMarkers) = y1(1:NumOfMarkers);
m3d1(3,1:NumOfMarkers) = z1(1:NumOfMarkers);

% Combine 3d marker coordinates of tiltseries2
m3d2(1,1:NumOfMarkers) = x2(1:NumOfMarkers);
m3d2(2,1:NumOfMarkers) = y2(1:NumOfMarkers);
m3d2(3,1:NumOfMarkers) = z2(1:NumOfMarkers);

% Move tiltseries1 3d marker coordinates of ReferenceMarker to [0 0 0]
m3d1(1,1:NumOfMarkers) = m3d1(1,1:NumOfMarkers) - Origin1(1) + (Imdim/2+1); 
m3d1(2,1:NumOfMarkers) = m3d1(2,1:NumOfMarkers) - Origin1(2) + (Imdim/2+1); 
m3d1(3,1:NumOfMarkers) = m3d1(3,1:NumOfMarkers) - Origin1(3) + (Imdim/2+1);

% Move tiltseries2 3d marker coordinates of ReferenceMarker to [0 0 0]
m3d2(1,1:NumOfMarkers) = m3d2(1,1:NumOfMarkers) - Origin2(1) + (Imdim/2+1); 
m3d2(2,1:NumOfMarkers) = m3d2(2,1:NumOfMarkers) - Origin2(2) + (Imdim/2+1); 
m3d2(3,1:NumOfMarkers) = m3d2(3,1:NumOfMarkers) - Origin2(3) + (Imdim/2+1);

% Calculate 3d rotation between tiltseries1 and tiltseries2
[Psi, Theta, Phi, RotMatrix] = tom_Rec3dEulerAngles...
(m3d1, m3d2, NumOfMarkers, ProtocolMode);

% Calculate new alignment origin of tiltseries2
Origin2 = tom_Rec3dNewOrigin2...
(Origin1, Origin2, Imdim, RotMatrix, ProtocolMode);

% Calculate 'rigidbody' alignment of tiltseries1
[Markerfile1.Value, Rho1, Sigma1, x1, y1, z1] = tom_Rec3dRigidBodyAlignment...
(Markerfile1.Value, ReferenceMarker, RefProj1, Origin1, Imdim, ProtocolMode);

% Calculate 'rigidbody' alignment of tiltseries2
[Markerfile2.Value, Rho2, Sigma2, x2, y2, z2] = tom_Rec3dRigidBodyAlignment...
(Markerfile2.Value, ReferenceMarker, RefProj2, Origin2, Imdim, ProtocolMode);

% Combine 3d marker coordinates of tiltseries1
m3d1(1,1:NumOfMarkers) = x1(1:NumOfMarkers);
m3d1(2,1:NumOfMarkers) = y1(1:NumOfMarkers);
m3d1(3,1:NumOfMarkers) = z1(1:NumOfMarkers);

% Combine 3d marker coordinates of tiltseries2
m3d2(1,1:NumOfMarkers) = x2(1:NumOfMarkers);
m3d2(2,1:NumOfMarkers) = y2(1:NumOfMarkers);
m3d2(3,1:NumOfMarkers) = z2(1:NumOfMarkers);

% Move tiltseries1 3d marker coordinates of ReferenceMarker to [0 0 0]
m3d1(1,1:NumOfMarkers) = m3d1(1,1:NumOfMarkers) - Origin1(1) + (Imdim/2+1); 
m3d1(2,1:NumOfMarkers) = m3d1(2,1:NumOfMarkers) - Origin1(2) + (Imdim/2+1); 
m3d1(3,1:NumOfMarkers) = m3d1(3,1:NumOfMarkers) - Origin1(3) + (Imdim/2+1);

% Move tiltseries2 3d marker coordinates of ReferenceMarker to [0 0 0]
m3d2(1,1:NumOfMarkers) = m3d2(1,1:NumOfMarkers) - Origin2(1) + (Imdim/2+1); 
m3d2(2,1:NumOfMarkers) = m3d2(2,1:NumOfMarkers) - Origin2(2) + (Imdim/2+1); 
m3d2(3,1:NumOfMarkers) = m3d2(3,1:NumOfMarkers) - Origin2(3) + (Imdim/2+1);

% Tiltaxis, translations, isoscale of tiltseries1
for k = 1:NumOfProj1
    Tiltaxis1(k) = (180/pi)*Rho1;
    tx1(k) = Markerfile1.Value(7, k, 1);
    ty1(k) = Markerfile1.Value(8, k, 1);
    isoscale1(k) = 1;
end

% Tiltaxis, translations, isoscale of tiltseries2
for k = 1:NumOfProj2
    Tiltaxis2(k) = (180/pi)*Rho2;
    tx2(k) = Markerfile2.Value(7, k, 1);
    ty2(k) = Markerfile2.Value(8, k, 1);
    isoscale2(k) = 1;
end

% Calculate ResidualMatrix and WarpAlignment of tiltseries1
[ResidualMatrix1, WarpAlignment1] = tom_Rec3dRigidBodyResidual...
(MarkersOnProj1,NumOfProj1,NumOfMarkers,Origin1,Imdim,Tiltaxis1,Tiltangles1,m3d1,tx1,ty1);

% Calculate ResidualMatrix and WarpAlignment of tiltseries2
[ResidualMatrix2, WarpAlignment2] = tom_Rec3dRigidBodyResidual...
(MarkersOnProj2,NumOfProj2,NumOfMarkers,Origin2,Imdim,Tiltaxis2,Tiltangles2,m3d2,tx2,ty2);

end
% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
% AlignmentMethod 'freetilt'
% -------------------------------------------------------------------------
if strcmp(AlignmentMethod, 'freetilt')
    % ERROR
    msgbox('AlignmentMethod freetilt not implemented!', ...
                   'Do Alignment', 'error');
    return;
end
% -------------------------------------------------------------------------

end
% -------------------------------------------------------------------------





% -------------------------------------------------------------------------
% RECONSTRUCTION / Parameter
% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
% TiltingGeometry 'singleaxis'
% -------------------------------------------------------------------------
if strcmp(TiltingGeometry, 'singleaxis')
    
    % Reconstruction parameter of Tiltseries1
    for k=1:NumOfProj1
    
    % Tiltangles
    Tiltangles(k) = Tiltangles1(k);
    % Tiltaxis
    Tiltaxis(k) = Tiltaxis1(k) + 90;
    % ProjDir
    ProjDir(k) =  Tiltaxis1(k) + 90;
    % tx
    tx(k) = tx1(k);
    % ty
    ty(k) = ty1(k);
    % isoscale
    isoscale(k) = isoscale1(k);
    
    end
     
end
% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
% TiltingGeometry 'dualaxis'
% -------------------------------------------------------------------------
if strcmp(TiltingGeometry, 'dualaxis')
    
    % Reconstruction parameter of Tiltseries1
    for k=1:NumOfProj1
    
    % Tiltangles
    Tiltangles(k) = Tiltangles1(k);
    % Tiltaxis
    Tiltaxis(k) = Tiltaxis1(k) + 90;
    % ProjDir
    ProjDir(k) =  Tiltaxis1(k) + 90;
    % tx
    tx(k) = tx1(k);
    % ty
    ty(k) = ty1(k);
    % isoscale
    isoscale(k) = isoscale1(k);

    end
    
    % Reconstruction parameter of Tiltseries2
    for k=1:NumOfProj2
    
    % Tiltangles
    a = cos((pi/180)*Theta)*cos((pi/180)*Tiltangles2(k)) - ...
        sin((pi/180)*Theta)*sin((pi/180)*Tiltangles2(k))*cos((pi/180)*(-Tiltaxis2(k)+Psi));

    b = sqrt(1-a*a);

    Tiltangles(k+NumOfProj1) = (180/pi)*atan2(b,a);
    
    % Tiltaxis
    a = cos((pi/180)*Theta)*sin((pi/180)*Tiltangles2(k)) + ...
        sin((pi/180)*Theta)*cos((pi/180)*Tiltangles2(k))*cos((pi/180)*(-Tiltaxis2(k)+Psi));

    b = sin((pi/180)*Theta)*sin((pi/180)*(-Tiltaxis2(k)+Psi));

    Tiltaxis(k+NumOfProj1) = (180/pi)*atan2(b,a) + Tiltaxis2(k) + 90;
    
    % ProjDir
    a = sin((pi/180)*Theta)*cos((pi/180)*Tiltangles2(k)) + ...
        cos((pi/180)*Theta)*sin((pi/180)*Tiltangles2(k))*cos((pi/180)*(-Tiltaxis2(k)+Psi));

    b = sin((pi/180)*Tiltangles2(k))*sin((pi/180)*(-Tiltaxis2(k)+Psi));

    ProjDir(k+NumOfProj1) = -(180/pi)*(atan2(b,a) + (pi/180)*Phi) + 90;
    
    % tx
    tx(k+NumOfProj1) = tx2(k);
    % ty
    ty(k+NumOfProj1) = ty2(k);
    % isoscale
    isoscale(k+NumOfProj1) = isoscale2(k);

    end

end
% -------------------------------------------------------------------------


% -------------------------------------------------------------------------
% ALIGNMENT / Residuals
% -------------------------------------------------------------------------
% Calculate AlignmentResidual of tiltseries1
[AveragePerProjection1, AveragePerMarker1] = ...
tom_Rec3dAlignmentResidual(ResidualMatrix1, NumOfMarkers, NumOfProj1);

% Calculate AlignmentResidual of tiltseries2
if strcmp(TiltingGeometry, 'dualaxis')
    [AveragePerProjection2, AveragePerMarker2] = ...
    tom_Rec3dAlignmentResidual(ResidualMatrix2, NumOfMarkers, NumOfProj2);
end

% Calculate EulerAnglesResidual
if strcmp(TiltingGeometry, 'dualaxis')
    [EulerAnglesResidual, ...
     MaximumResidual, ...
     ResidualSpheres, ...
     AverageResidualSphere] = tom_Rec3dEulerAnglesResidual...
    (m3d1, m3d2, NumOfMarkers, RotMatrix, ProtocolMode);
end
% -------------------------------------------------------------------------


% -------------------------------------------------------------------------
% Check
% -------------------------------------------------------------------------
% Check Translations
if strcmp(TiltingGeometry, 'singleaxis')
    for k=1:NumOfProj1
         if tx1(k) > Imdim || ty1(k) > Imdim
         % ERROR
         msgbox('Alignment false! Click more marker or delete projections!', ...
                'Do Alignment', 'error');
         return;
         end
    end
end
if strcmp(TiltingGeometry, 'dualaxis')
    for k=1:NumOfProj1
         if tx1(k) > Imdim || ty1(k) > Imdim
         % ERROR
         msgbox('Alignment false! Click more marker or delete projections!', ...
                'Do Alignment', 'error');
         return;
         end  
    end
    for k=1:NumOfProj2
         if tx2(k) > Imdim || ty2(k) > Imdim
         % ERROR
         msgbox('Alignment false! Click more marker or delete projections!', ...
                'Do Alignment', 'error');
         return;
         end
    end
end    
% -------------------------------------------------------------------------


% -------------------------------------------------------------------------
% AlignmentEngine -> Rec3dProject
% -------------------------------------------------------------------------
% HEAD
Rec3dProject.ProjectStatus = 'aligned';
% ALIGNMENT
Rec3dProject.ALIGNMENT.Origin1 = Origin1;
Rec3dProject.ALIGNMENT.m3d1 = m3d1;
Rec3dProject.ALIGNMENT.Tiltaxis1 = Tiltaxis1;
Rec3dProject.ALIGNMENT.tx1 = tx1;
Rec3dProject.ALIGNMENT.ty1 = ty1;
Rec3dProject.ALIGNMENT.isoscale1 = isoscale1;
Rec3dProject.ALIGNMENT.WarpAlignment1 = WarpAlignment1;
if strcmp(TiltingGeometry, 'dualaxis')
    Rec3dProject.ALIGNMENT.Origin2 = Origin2;
    Rec3dProject.ALIGNMENT.m3d2 = m3d2;
    Rec3dProject.ALIGNMENT.Tiltaxis2 = Tiltaxis2;
    Rec3dProject.ALIGNMENT.tx2 = tx2;
    Rec3dProject.ALIGNMENT.ty2 = ty2;
    Rec3dProject.ALIGNMENT.isoscale2 = isoscale2;
    Rec3dProject.ALIGNMENT.WarpAlignment2 = WarpAlignment2;
    Rec3dProject.ALIGNMENT.RotMatrix = RotMatrix;
    Rec3dProject.ALIGNMENT.Psi = Psi;
    Rec3dProject.ALIGNMENT.Theta = Theta;
    Rec3dProject.ALIGNMENT.Phi = Phi;
end
% ALGRESIDUALS
Rec3dProject.ALGRESIDUALS.ResidualMatrix1 = ResidualMatrix1;
Rec3dProject.ALGRESIDUALS.Sigma1 = Sigma1;
Rec3dProject.ALGRESIDUALS.AveragePerProjection1 = AveragePerProjection1;
Rec3dProject.ALGRESIDUALS.AveragePerMarker1 = AveragePerMarker1;
if strcmp(TiltingGeometry, 'dualaxis')
    Rec3dProject.ALGRESIDUALS.ResidualMatrix2 = ResidualMatrix2;
    Rec3dProject.ALGRESIDUALS.Sigma2 = Sigma2;
    Rec3dProject.ALGRESIDUALS.AveragePerProjection2 = AveragePerProjection2;
    Rec3dProject.ALGRESIDUALS.AveragePerMarker2 = AveragePerMarker2;
    Rec3dProject.ALGRESIDUALS.EulerAnglesResidual = EulerAnglesResidual;
    Rec3dProject.ALGRESIDUALS.MaximumResidual = MaximumResidual;
    Rec3dProject.ALGRESIDUALS.ResidualSpheres = ResidualSpheres;
    Rec3dProject.ALGRESIDUALS.AverageResidualSphere = AverageResidualSphere;
end
% PARAMETER
Rec3dProject.PARAMETER.Tiltangles = Tiltangles;
Rec3dProject.PARAMETER.Tiltaxis = Tiltaxis;
Rec3dProject.PARAMETER.ProjDir = ProjDir;
Rec3dProject.PARAMETER.tx = tx;
Rec3dProject.PARAMETER.ty = ty;
Rec3dProject.PARAMETER.isoscale = isoscale;
% -------------------------------------------------------------------------


% MESSAGE
uiwait(msgbox('Alignment done!', ...
              'Do Alignment', 'warn'));

          