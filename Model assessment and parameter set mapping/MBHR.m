%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function is used for finding the minimum bounding hyper-rectangle
% containing typical trajectories of the scheduling variable of an LPV
% model.
% 
% Input: The reduced scheduling variable data obtained from "PSM.m".
%
% reference : LPV control of nonlinear systems, Bardia Sharif.
%
% Bardia Sharif
% Eindhoven University of Technology
% August 2018
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Delta  =   MBHR(delta)
%% set Yalmip options
ops          = sdpsettings('savesolveroutput', 1, ...
               'savesolverinput' , 1, ...
               'verbose'         , 1, ...
               'solver'          ,'sedumi');

%% find the number of the scheduling varibles
m   =   size(delta,1);
%% Depending on the number of scheduling variables find the minimum bounding hyper rectangle.

if m==1
P1_data     =   delta(1,:);       % data for the fist schedling variable.
P1_low      =   min(P1_data);     % lower bound on the first scheduling variable.
P1_high     =   max(P1_data);     % upper bound on the first scheduling variable. 
p1          =   tvreal('p1',[P1_low P1_high]);  % first scheduling variable        
% the coordinates of the vertices of the convex hull
Xi          =   P1_data;
Xi          =   Xi';
n           =   size(Xi,1);

% set up the LMI for centering

lambda          =   sdpvar(1);          % the radius squared of the circle including the points
C_circle        =   sdpvar(1);          % the center of the circle including the points
 
cons    = lambda>=0;

for  i = 1:n
    vert        =   Xi(i);
   
    lmi_term    =   [lambda             (vert - C_circle)';...
                     (vert - C_circle)   eye(1)];
                 
   cons         =   [cons,lmi_term>=0];
                 
end

sol                 =   optimize(cons,lambda,ops);

X_offset            = value(C_circle(1));
X_recentered        = Xi -X_offset;
P1_new_low          = X_recentered(1);
P1_new_high         = X_recentered(2);
P1_new              = tvreal('P1_new',[P1_new_low P1_new_high]);
Delta.P1            = P1_new;  

else
    if m==2
        P1_data     =   delta(1,:);       % data for the fist schedling variable.
        P1_low      =   min(P1_data);     % lower bound on the first scheduling variable.
        P1_high     =   max(P1_data);     % upper bound on the first scheduling variable. 

        P2_data     =   delta(2,:);       % data for the second schedling variable.
        P2_low      =   min(P2_data);     % lower bound on the second scheduling variable.
        P2_high     =   max(P2_data);     % upper bound on the second scheduling variable.

        p1          =   tvreal('p1',[P1_low P1_high]);  % first scheduling variable
        p2          =   tvreal('p2',[P2_low P2_high]);  % second scheduling variable
% the coordinates of the vertices of the convex hull
        Xi  =   P1_data;%(P_set_tight);
        Xi  =   Xi';
        Yi  =   P2_data;%(P_set_tight);
        Yi  =   Yi';
        n   =   size(Xi,1);

% set up the LMI for centering

      lambda          =   sdpvar(1);          % the radius squared of the circle including the points
      C_circle        =   sdpvar(2,1);        % the center of the circle including the points
 
      cons            =   lambda>=0;

      for  i = 1:n
          vert        =   [Xi(i);Yi(i)];
   
          lmi_term    =   [lambda             (vert - C_circle)';...
                           (vert - C_circle)   eye(2)];
                 
          cons     =   [cons,lmi_term>=0];
                 
      end
        
      sol                 =   optimize(cons,lambda,ops);
      
X_offset            = value(C_circle(1));
Y_offset            = value(C_circle(2));
X_recentered        = Xi -X_offset;
Y_recentered        = Yi- Y_offset;

% angle grid
tt = degtorad(0:1:360);

for i = 1:length(tt)

rot_data = [cos(tt(i)) -sin(tt(i));sin(tt(i)) cos(tt(i))]*[X_recentered';Y_recentered'];

P1_rot_low      =   min(rot_data(1,:));     % lower bound on the first scheduling variable.
P1_rot_high     =   max(rot_data(1,:));     % upper bound on the first scheduling variable. 

P2_rot_low      =   min(rot_data(2,:));     % lower bound on the second scheduling variable.
P2_rot_high     =   max(rot_data(2,:));     % upper bound on the second scheduling variable.

P1_rot_range            =   [P1_rot_low P1_rot_high];        
P2_rot_range            =   [P2_rot_low P2_rot_high];
vertices_extreme        =   combvec(P1_rot_range,P2_rot_range);                      % vertices of the new set
x_rot_new               =   vertices_extreme(1,:);                            
y_rot_new               =   vertices_extreme(2,:);
P_set_rot_new           =   convhull(x_rot_new,y_rot_new);

area(i) = polyarea(x_rot_new(P_set_rot_new),y_rot_new(P_set_rot_new));

end

j = find(area == min(area),1);

SS          = sin(tt(j)); 
CC          = cos(tt(j));                        % rotation matrix giving the best rotation for finding the rectangle of least area
R           = [CC -SS;SS CC];   
delta_star  = R*delta;                           % the rotated scheduling variable.
P1_new_low  =  min(delta_star(1,:));
P1_new_high =  max(delta_star(1,:));

P2_new_low  =   min(delta_star(2,:));
P2_new_high =   max(delta_star(2,:));

P1_new      =   tvreal('P1_new',[P1_new_low P1_new_high]);
P2_new      =   tvreal('P2_new',[P2_new_low P2_new_high]);

P_rot       =   R\[P1_new;P2_new];

Delta.p1    =   P_rot(1);
Delta.p2    =   P_rot(2);


    
    else 
        if m==3 
            P1_data     =   delta(1,:);       % data for the fist schedling variable.
            P1_low      =   min(P1_data);     % lower bound on the first scheduling variable.
            P1_high     =   max(P1_data);     % upper bound on the first scheduling variable. 

            P2_data     =   delta(2,:);       % data for the second schedling variable.
            P2_low      =   min(P2_data);     % lower bound on the second scheduling variable.
            P2_high     =   max(P2_data);     % upper bound on the second scheduling variable.

            P3_data     =   delta(3,:);       % data for the third schedling variable.
            P3_low      =   min(P3_data);     % lower bound on the third scheduling variable.
            P3_high     =   max(P3_data);     % upper bound on the third scheduling variable.

            p1          =   tvreal('p1',[P1_low P1_high]);  % first scheduling variable
            p2          =   tvreal('p2',[P2_low P2_high]);  % second scheduling variable
            p3          =   tvreal('p3',[P3_low P3_high]);  % third scheduling variable
        % the coordinates of the vertices of the convex hull
            Xi  =   P1_data;
            Xi  =   Xi';
            Yi  =   P2_data;
            Yi  =   Yi';
            Zi  =   P3_data;
            Zi  =   Zi'; 
            n   =   size(Xi,1);

        % set up the LMI for centering

        lambda          =   sdpvar(1);          % the radius squared of the sphere including the points
        C_circle        =   sdpvar(3,1);        % the center of the sphere including the points
 
        cons    = lambda>=0;

        for  i = 1:n
        vert        =   [Xi(i);Yi(i);Zi(i)];
   
        lmi_term    =   [lambda             (vert - C_circle)';...
                         (vert - C_circle)   eye(2)];
                 
        cons        =   [cons,lmi_term>=0];
                 
        end

        sol                 =   optimize(cons,lambda,ops);
        
        X_offset            =   value(C_circle(1));
        Y_offset            =   value(C_circle(2));
        Z_offset            =   value(C_circle(3));
        X_recentered        =   Xi -X_offset;
        Y_recentered        =   Yi- Y_offset;
        Z_recentered        =   Zi- Z_offset;
        
        ttx                 =   degtorad(0:1:360);
        tty                 =   degtorad(0:1:360);
        ttz                 =   degtorad(0:1:360);
        
    for i = 1:length(ttx)
        
    Cx              =   cos(ttx(i));
    Sx              =   sin(ttx(i));
    Rx              =   [1 0 0;0 Cx -Sx;0 Sx Cx];
    rot_data        =   Rx*[X_recentered';Y_recentered';Z_recentered'];
        for j = 1:length(tty) 
           Cy              =   cos(tty(j));
           Sy              =   sin(tty(j));
           Ry              =   [Cy 0 Sy;0 1 0;-Sy 0 Cy];
           rot_data        =   Ry*rot_data;
            for k = 1:length(ttz)
               Cz             =   cos(ttz(k));
               Sz             =   sin(ttz(k));
               Rz             =   [Cz -Sz 0;Sz Cz 0;0 0 1];
               rot_data       =   Rz*rot_data;
               P1_rot_low     =   min(rot_data(1,:));     % lower bound on the first scheduling variable.
               P1_rot_high    =   max(rot_data(1,:));     % upper bound on the first scheduling variable. 

               P2_rot_low     =   min(rot_data(2,:));     % lower bound on the second scheduling variable.
               P2_rot_high    =   max(rot_data(2,:));     % upper bound on the second scheduling variable.
    
               P1_rot_range            =   [P1_rot_low P1_rot_high];           
               P2_rot_range            =   [P2_rot_low P2_rot_high];
               vertices_extreme        =   combvec(P1_rot_range,P2_rot_range);                      % vertices of the new set
               x_rot_new               =   vertices_extreme(1,:);                            
               y_rot_new               =   vertices_extreme(2,:);
               P_set_rot_new           =   convhull(x_rot_new,y_rot_new);
              
               area = polyarea(x_rot_new(P_set_rot_new),y_rot_new(P_set_rot_new));
                
               set(k).rotangle  =    [ttx(i);tty(j);ttz(k)];
               set(k).area      =    area;
               
            end
            jj                =    find(set.area== min(set.area),1);
            set2(j).rotangle  =    [set(jj).rotangle(1);set(jj).rotangle(2);set(jj).rotangle(3)];
            set2(j).area      =    set(jj).area;
            
        end
        
        ii                    =   find(set2.area== min(set2.area),1);
        set3(i).rotangle      =  [set2(ii).rotangle(1);set2(ii).rotangle(2);set2(ii).rotangle(3)];
        set3(i).area          =  set2(ii).area;
               

    end

        kk            =  find(set3.area== min(set3.area),1);
        rotangle      =  [set3(kk).rotangle2(1);set3(kk).rotangle(2);set3(kk).rotangle(3)];
        area          =  set3(kk).area;
        
        SSx           =  sin(rotangle(1)); 
        CCx           =  cos(rotangle(1));                                  % rotation matrix giving the best rotation around the X axis
        Rx            =  [1 0 0;0 CCx -SSx;0 SSx CCx];

        SSy           =  sin(rotangle(2)); 
        CCy           =  cos(rotangle(2));                                  % rotation matrix giving the best rotation around the Y axis
        Ry            =  [CCy 0 SSy;0 1 0;-SSy 0 CCy];

        
        SSz           =  sin(rotangle(3)); 
        CCz           =  cos(rotangle(3));                                  % rotation matrix giving the best rotation around the Y axis
        Rz            =  [CCz -SSz 0;SSz CCz 0;0 0 1];
  
        R             =  Rz*Ry*Rx;  
        
        delta_star    =  R*delta;                           % the rotated scheduling variable.
        P1_new_low    =  min(delta_star(1,:));
        P1_new_high   =  max(delta_star(1,:));

        P2_new_low    =   min(delta_star(2,:));
        P2_new_high   =   max(delta_star(2,:));
        
        P3_new_low    =  min(delta_star(3,:));
        P3_new_high   =  max(delta_star(3,:));

        P1_new        =   tvreal('P1_new',[P1_new_low P1_new_high]);
        P2_new        =   tvreal('P2_new',[P2_new_low P2_new_high]);
        P3_new        =   tvreal('P3_new',[P3_new_low P3_new_high]);

        P_rot         =   R\[P1_new;P2_new;P3_new];

        Delta.p1      =   P_rot(1);
        Delta.p2      =   P_rot(2);
        Delta.p3      =   P_rot(3);
        else
            Delta   =   [];
            disp('The dimension of the parameter space is too high for this method!');
        end
    end
end

end