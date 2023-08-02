function varargout = Falling(varargin)
% FALLING MATLAB code for Falling.fig
%      FALLING, by itself, creates a new FALLING or raises the existing
%      singleton*.
%
%      H = FALLING returns the handle to a new FALLING or the handle to
%      the existing singleton*.
%
%      FALLING('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in FALLING.M with the given input arguments.
%
%      FALLING('Property','Value',...) creates a new FALLING or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before Falling_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to Falling_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help Falling

% Last Modified by GUIDE v2.5 24-Mar-2014 15:39:59

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @Falling_OpeningFcn, ...
                   'gui_OutputFcn',  @Falling_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before Falling is made visible.
function Falling_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to Falling (see VARARGIN)

% Choose default command line output for Falling
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes Falling wait for user response (see UIRESUME)
% uiwait(handles.figure1);
grid on
xlim([0, 100])
ylim([0, 100])


% --- Outputs from this function are returned to the command line.
function varargout = Falling_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in compute.
function compute_Callback(hObject, eventdata, handles)
% hObject    handle to compute (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%RETREIVE DATA FROM INPUTS%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

handles.mxrockarray1  = [];
handles.myrockarray1  = [];
handles.mvrockarray1  = [];
handles.mvxrockarray1 = [];
handles.mvyrockarray1 = [];
handles.mvrrockarray1 = [];
handles.mtranslationalKEarray1 = [];
handles.mrotationalKEarray1    = [];
handles.mbouncehtarray1        = [];

iterations1 = handles.iterations1;
xpoint      = handles.xpoint1;
ypoint      = handles.ypoint1;


t_step     = str2num(get(handles.t_step,'String'));
g          = str2num(get(handles.g,'String'));
vfactor    = str2num(get(handles.vscaling,'String'));

radius     = str2num(get(handles.radius,'String'));
mass       = str2num(get(handles.mass,'String'));
num        = str2num(get(handles.num,'String'));

xi         = str2num(get(handles.xi,'String'));
xstd       = str2num(get(handles.xistd,'String'));
yi         = str2num(get(handles.yi,'String'));
ystd       = str2num(get(handles.yistd,'String'));

vxi        = str2num(get(handles.vxi,'String'));
vxstd      = str2num(get(handles.vxistd,'String'));
vyi        = str2num(get(handles.vyi,'String'));
vystd      = str2num(get(handles.vyistd,'String'));
vrot       = str2num(get(handles.vri,'String'));
vrotstd    = str2num(get(handles.vrstd,'String'));

howmany    = str2num(get(handles.numberofcoordinates,'String'));

rt         = str2num(get(handles.rt,'String'));
rt_std     = str2num(get(handles.rtstd,'String'));
rn         = str2num(get(handles.rn,'String'));
rn_std     = str2num(get(handles.rnstd,'String'));
mewr       = str2num(get(handles.rollingfriction,'String'));
properties = zeros(4,(2*(howmany-1)));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%END RETREIVING DATA FROM INPUT%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%CREATING PROPERTIES%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Beggining plotting
num_point=howmany;

%FITTING DATA THROUGH POINTS OF SLOPE
m_slope   = zeros((num_point-1),1);
b_slope   = zeros((num_point-1),1);
div       = xpoint(end)-xpoint(1);
%CALCULATING ANGLES FOR EACH CELL
angle     = zeros((num_point-1),1);
angle_num = zeros((num_point-1),1);

for q = 1:(num_point-1)
    angle(q,1)     = radtodeg(atan((abs(ypoint(q)-ypoint(q+1)))/(abs(xpoint(q+1)-xpoint(q)))));
    angle_num(q,1) = q;
end

for r = 1:(num_point-1)      
    properties(1,((r*2)-1)) = angle(r);
    properties(2,((r*2)-1)) = mewr(r,1);
    properties(3,((r*2)-1)) = rn(r,1);
    properties(4,((r*2)-1)) = rt(r,1); 
    properties(5,((r*2)-1)) = r;
    properties(6,((r*2)-1)) = xpoint(r,1);
    properties(7,((r*2)-1)) = ypoint(r,1);
    properties(8,((r*2)-1)) = ((abs(g))*(cos(degtorad(angle(r)))*(mewr(r,1))));
    
    properties(1,(r*2)) = 0;
    properties(2,(r*2)) = 0;
    properties(3,(r*2)) = rn_std(r,1);
    properties(4,(r*2)) = rt_std(r,1);
    properties(5,(r*2)) = 0;
    properties(6,(r*2)) = xpoint((r+1),1);
    properties(7,(r*2)) = ypoint((r+1),1);
    properties(8,(r*2)) = 0;
end
properties;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%END CREATING PROPERTIES%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%START CALCULATION LOOPS%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for u = 1:num
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%START SETUP OF INITIAL PARAMETERS%%%%%%%%%%%%%%%%
    xi1   = xi+ xstd*(abs(randn(1,1)));    
    yi1   = yi+ ystd*(abs(randn(1,1)));                   
    vxi1  = (vxi+ vxstd*randn(1,1));
    vyi1  = (vyi+ vystd*randn(1,1));
    vrot1 = (vrot+vrotstd*randn(1,1));
    properties1 = properties;      
    
    for i = 1:(num_point-1)    
        properties1(3,((i*2)-1)) = properties1(3,((i*2)-1)) +  (properties1(3,(i*2)) * randn(1,1));
        properties1(4,((i*2)-1)) = properties1(4,((i*2)-1)) +  (properties1(4,(i*2)) * randn(1,1));
        properties1(3,(i*2))     = 0;
        properties1(4,(i*2))     = 0;
    end
    
    x  = xi1;  
    y  = yi1;
    vx = vxi1; 
    vy = vyi1;
    vr = vrot1;
    
    for h = 1:(num_point-1)
         if  x(1,1) >= properties1(6,((h*2)-1)) &  x(1,1) <= properties1(6,(h*2))
         slope_num(1,1)= h;
         end
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %FINISHED SETUP OF INITIAL PARAMETERS. X AND Y FRAME OF REFERENCE%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%START INDIVIDUAL ROCK CALCULATIONS%%%%%%%%%%%%%%%
   
    for o = 1:iterations1-1
       
        %%%%%%%%%%%%%%%%%PROJECTILE TRAJECTORY CODE%%%%%%%%%%%%%%%%%%%
        vystep     = vy(o,1) + (g)*(t_step);                          
        vy(o+1,1)  = vystep;                                          
        vxstep     = vx(o,1) + (0)*(t_step);                              
        vx(o+1,1)  = vxstep;                       
        vr(o+1,1)  = vr(o,1);
        xstep      = x(o,1)+ (vx(o,1)*t_step);                             
        x(o+1,1)   = xstep;        
        ystep      = y(o,1) + ((vy(o+1,1) + vy(o,1))*(t_step)*0.5);
        y(o+1,1)   = ystep;
        %%%%%%%%%%%%END PROJECTILE TRAJECTORY CODE%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %%%%%%%%%%%%%%%%SLOPE LOCATION BELOW ROCK CODE%%%%%%%%%%%%%%%
        for p = 1:(num_point-1)
            
            if  x(o+1,1) >= properties1(6,((p*2)-1)) &  x(o+1,1) <= properties1(6,(p*2))
                if  properties1(7,(p*2)) <  properties1(7,((p*2)-1))
                    y_impact_location = (properties1(6,(p*2))- x(o+1,1)) * (tan(degtorad(properties1(1,((p*2)-1))))) + properties1(7,(p*2));
                    slope_num(o+1,1)  = p;
                else 
                    y_impact_location = (x(o+1,1)-properties1(6,((p*2)-1))) * (tan(degtorad(properties1(1,((p*2)-1))))) + properties1(7,((p*2)-1));
                    slope_num(o+1,1)  = p;
                end 
                break
            end
            
        end 
        %%%%%%%%%%%END SLOPE LOCATION BELOW ROCK CODE%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %%%%%%%%IMPACT LOCATION AND RESTITUTION CODE%%%%%%%%%%%%%%%%
        if y(o+1,1) - y_impact_location <=0.1
            
            if  properties1(7,(p*2)) <  properties1(7,((p*2)-1))  %Downhill 
                if vx(o+1,1) > 0 %Vx is positive
                    vtb     = abs((vy(o+1,1)*(sin((degtorad(properties1(1,((p*2)-1))))))))+(vx(o+1,1)*(cos((degtorad(properties1(1,((p*2)-1)))))));  
                    vnb     = (vy(o+1,1)*cos((degtorad(properties1(1,((p*2)-1))))))+(vx(o+1,1)*sin((degtorad(properties1(1,((p*2)-1))))));                       
                    vta     = vtb*(properties1(4,((p*2)-1)));                                                 
                    vna     = (abs(vnb))*((properties1(3,((p*2)-1)))/(1+(((abs(vnb)))/vfactor)^2));
                    vx(o,1) = (vna*sin((degtorad(properties1(1,((p*2)-1))))))+(vta*cos((degtorad(properties1(1,((p*2)-1))))));                
                    vy(o,1) = (vna*cos((degtorad(properties1(1,((p*2)-1))))))+((-1)*(vta*sin((degtorad(properties1(1,((p*2)-1)))))));
                    vr(o,1) = ((((vx(o,1))^2 + (vy(o,1))^2)^(0.5))/radius);
                else    %Vx is negative        
                    vtb     = abs(vy(o+1,1)*(sin((degtorad(properties1(1,((p*2)-1)))))))+(vx(o+1,1)*(cos((degtorad(properties1(1,((p*2)-1)))))));  
                    vnb     = (vy(o+1,1)*cos((degtorad(properties1(1,((p*2)-1))))))+(vx(o+1,1)*sin((degtorad(properties1(1,((p*2)-1))))));                       
                    vta     = vtb*(properties1(4,((p*2)-1)));                                                 
                    vna     = (abs(vnb))*((properties1(3,((p*2)-1)))/(1+(((abs(vnb)))/vfactor)^2));
                    vx(o,1) = (vna*sin((degtorad(properties1(1,((p*2)-1))))))+(vta*cos((degtorad(properties1(1,((p*2)-1))))));                
                    vy(o,1) = (vna*cos((degtorad(properties1(1,((p*2)-1))))))+((-1)*(vta*sin((degtorad(properties1(1,((p*2)-1)))))));
                    vr(o,1) = ((((vx(o,1))^2 + (vy(o,1))^2)^(0.5))/radius);
                end     
            else %Uphill
                if vx(o+1,1)>0 %Vx is positive
                    vtb     = (vy(o+1,1)*(sin((degtorad(properties1(1,((p*2)-1)))))))+(vx(o+1,1)*(cos((degtorad(properties1(1,((p*2)-1)))))));   
                    vnb     = (vy(o+1,1)*cos((degtorad(properties1(1,((p*2)-1))))))+((-1)*(vx(o+1,1)*sin((degtorad(properties1(1,((p*2)-1))))))); 
                    vta     = vtb*(properties1(4,((p*2)-1)));
                    vna     = (abs(vnb))*((properties1(3,((p*2)-1)))/(1+(((abs(vnb)))/vfactor)^2));
                    vx(o,1) = ((-1)*(vna*sin((degtorad(properties1(1,((p*2)-1)))))))+(vta*cos((degtorad(properties1(1,((p*2)-1))))));                
                    vy(o,1) = (vna*cos((degtorad(properties1(1,((p*2)-1))))))+((vta*sin((degtorad(properties1(1,((p*2)-1)))))));
                    vr(o,1) = ((((vx(o,1))^2 + (vy(o,1))^2)^(0.5))/radius);
                else %Vx is negative
                    vtb     = (vy(o+1,1)*(sin((degtorad(properties1(1,((p*2)-1)))))))+(vx(o+1,1)*(cos((degtorad(properties1(1,((p*2)-1)))))));   
                    vnb     = (vy(o+1,1)*cos((degtorad(properties1(1,((p*2)-1))))))+abs(vx(o+1,1)*sin((degtorad(properties1(1,((p*2)-1)))))); 
                    vta     = (vtb*(properties1(4,((p*2)-1))));
                    vna     = (abs(vnb))*((properties1(3,((p*2)-1)))/(1+(((abs(vnb)))/vfactor)^2));
                    vx(o,1) = ((-1)*(vna*sin((degtorad(properties1(1,((p*2)-1)))))))+(vta*cos((degtorad(properties1(1,((p*2)-1))))));                
                    vy(o,1) = (vna*cos((degtorad(properties1(1,((p*2)-1))))))+((vta*sin((degtorad(properties1(1,((p*2)-1)))))));
                    vr(o,1) = ((((vx(o,1))^2 + (vy(o,1))^2)^(0.5))/radius);            
                end
                
            end                           
            %%%%%%%%END IMPACT LOCATION AND RESTITUTION CODE%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
        
           if vna >0.1
               %%%%%%%%TRAJECTORY AFTER IMPACT CODE%%%%%%%%%%%%%%%%%%%%%
               for a =  o:(iterations1-1)
                   vyf         = vy(a,1) + (g)*(t_step);                               
                   vy(a+1,1)   = vyf;                                                 
                   vxf         = vx(a,1)+(0)*(t_step);                                     
                   vx(a+1,1)   = vxf;
                   vr(a+1,1)   = vr(a,1);
                   xf          = x(a,1)+ (vx(a,1)*t_step);                                 
                   x(a+1,1)    = xf;
           
                   if properties1(7,(p*2)) <  properties1(7,((p*2)-1))
                       yf = y(a,1) + ((vy(a+1,1) + vy(a,1))* (0.5) * (t_step));
                   else
                       yf = y(a,1) + abs((vy(a+1,1) + vy(a,1))* (0.5) * (t_step));
                   end    
                   y(a+1,1) = yf;  
           
                   for p = 1:(num_point-1)
               
                       if  x(a+1,1) >= properties1(6,((p*2)-1)) &  x(a+1,1) <= properties1(6,(p*2))
                            if  properties1(7,(p*2)) <  properties1(7,((p*2)-1))
                                y_impact_location = (properties1(6,(p*2))- x(a+1,1)) * (tan(degtorad(properties1(1,((p*2)-1))))) + properties1(7,(p*2));
                                slope_num(a+1,1)  = p;
                            else 
                                y_impact_location = (x(a+1,1)-properties1(6,((p*2)-1))) * (tan(degtorad(properties1(1,((p*2)-1))))) + properties1(7,((p*2)-1));
                                slope_num(a+1,1)  = p;
                            end 
                        break
                       end
                   end
           
                   if vxf < 0.1 | vyf < 0.1
                       break
                   end
           
               end 
           else  
               %%%%%%%%%%%SLIDING AFTER IMPACT CODE%%%%%%%%%%%%%%%%%%%%%
               vxsmatrix(1,1) = vta; 
               for s = o:iterations1-1  
 
                   c = (1+(s-o));
                   n = (2+(s-o));
  
                   if vxsmatrix(c,1)>0
                       if  properties1(7,(p*2)) <  properties1(7,((p*2)-1))
                       part1 = (0.75)*(vxsmatrix(c,1))^2;
                       part2 = (((abs(g))*(vxsmatrix(c,1))*(t_step)) * ((cos(degtorad(properties(1,((p*2)-1)))) * properties(2,((p*2)-1))) - (sin(degtorad(properties(1,((p*2)-1)))))));
                       vxsmatrix(n,1) = (abs((4*(part1-part2))/3))^0.5;
                           if vxsmatrix(n,1)<0.1
                               vxsmatrix = vxsmatrix((1:c),1);
                               vr = vr((1:(s)),1);
                               x  = x((1:(s)),1);
                               y  = y((1:(s)),1);
                               slope_num = slope_num((1:(s)),1);
                               break  
                           end
                       else
                       part1 = (0.75)*(vxsmatrix(c,1))^2;
                       part2 = (((abs(g))*(vxsmatrix(c,1))*(t_step)) * ((sin(degtorad(properties(1,((p*2)-1)))))  +(cos(degtorad(properties(1,((p*2)-1)))) * properties(2,((p*2)-1)))));
                       vxsmatrix(n,1) = (abs((4*(part1-part2))/3))^0.5;
                           if vxsmatrix(n,1)<0.1
                               vxsmatrix = vxsmatrix((1:c),1);
                               vr = vr((1:(s)),1);
                               x  = x((1:(s)),1);
                               y  = y((1:(s)),1);
                               slope_num = slope_num((1:(s)),1);
                               break  
                           end  
                       end
                   else
                       if  properties1(7,(p*2)) <  properties1(7,((p*2)-1))
                       part1 = (0.75)*(vxsmatrix(c,1))^2;
                       part2 = (((abs(g))*(vxsmatrix(c,1))*(t_step)) * ((sin(degtorad(properties(1,((p*2)-1)))))  +(cos(degtorad(properties(1,((p*2)-1)))) * properties(2,((p*2)-1)))));
                       vxsmatrix(n,1) = (abs((4*(part1-part2))/3))^0.5;
                           if vxsmatrix(n,1)>-0.1
                               vxsmatrix = vxsmatrix((1:c),1);
                               vr = vr((1:(s)),1);
                               x  = x((1:(s)),1);
                               y  = y((1:(s)),1);
                               slope_num = slope_num((1:(s)),1);
                               break
                           end
                       else                       
                       part1 = (0.75)*(vxsmatrix(c,1))^2;
                       part2 = (((abs(g))*(vxsmatrix(c,1))*(t_step)) * ((cos(degtorad(properties(1,((p*2)-1)))) * properties(2,((p*2)-1))) - (sin(degtorad(properties(1,((p*2)-1)))))));
                       vxsmatrix(n,1) = (abs((4*(part1-part2))/3))^0.5;
                           if vxsmatrix(n,1)>-0.1
                               vxsmatrix = vxsmatrix((1:c),1);
                               vr = vr((1:(s)),1);
                               x  = x((1:(s)),1);
                               y  = y((1:(s)),1);
                               slope_num = slope_num((1:(s)),1);
                               break
                           end
                       end
                   end
           
                   distanceslid      = t_step*vxsmatrix(c,1);
                   xdistanceslid     = (distanceslid*(cos((degtorad(properties1(1,((p*2)-1)))))));      
                   x(s+1,1)          = x(s)+xdistanceslid;  
                   ydistanceslid     = (distanceslid*(sin((degtorad(properties1(1,((p*2)-1)))))));
           
                   %CHECKING IF CELL IS ASCENDING OR DESCENDING
                   if  properties1(7,(p*2)) <  properties1(7,((p*2)-1))
                       y(s+1,1) = y(s)-ydistanceslid;
                   else 
                       y(s+1,1) = y(s)+ydistanceslid;
                   end 
                   slope_num(s+1,1) = p;   
                   
                   if  properties1(7,(p*2)) <  properties1(7,((p*2)-1))
                       if vxsmatrix(n,1)>0
                       vx(s+1,1) = (vxsmatrix(n,1)*cos((degtorad(properties1(1,((p*2)-1))))));                
                       vy(s+1,1) = -(vxsmatrix(n,1)*sin((degtorad(properties1(1,((p*2)-1))))));
                       vr(s+1,1) = ((((vx(s,1))^2 + (vy(s,1))^2)^(0.5))/radius);
                       else
                       vx(s+1,1) = (vxsmatrix(n,1)*cos((degtorad(properties1(1,((p*2)-1))))));                
                       vy(s+1,1) = (vxsmatrix(n,1)*sin((degtorad(properties1(1,((p*2)-1))))));
                       vr(s+1,1) = ((((vx(s,1))^2 + (vy(s,1))^2)^(0.5))/radius);
                       end
                   else
                       if vxsmatrix(n,1)>0
                       vx(s+1,1) = (vxsmatrix(n,1)*cos((degtorad(properties1(1,((p*2)-1))))));                
                       vy(s+1,1) = (vxsmatrix(n,1)*sin((degtorad(properties1(1,((p*2)-1))))));
                       vr(s+1,1) = ((((vx(s,1))^2 + (vy(s,1))^2)^(0.5))/radius);
                       else
                       vx(s+1,1) = (vxsmatrix(n,1)*cos((degtorad(properties1(1,((p*2)-1))))));                
                       vy(s+1,1) = -(vxsmatrix(n,1)*sin((degtorad(properties1(1,((p*2)-1))))));
                       vr(s+1,1) = ((((vx(s,1))^2 + (vy(s,1))^2)^(0.5))/radius);  
                       end
                   end      

                   %CHECKING IF ROCK HAS REACHED END OF CELL
                   if vxsmatrix(n,1)>0
                       if x(s+1,1) > properties1(6,(p*2))
                           p = p+1;
                           if p ==0 
                               break
                           end
                           if p>(num_point-1)
                               break
                           end                           
                       end
                   else
                       if x(s+1,1) < properties1(6,((p*2)-1))
                           p = p-1;
                           if p == 0
                               break
                           end
                           if p>(num_point-1)
                               break
                           end 
                       end
                   end  
               end   
               break
           end    
        end
    end
    
vprojectile = ((vx).^2+(vy).^2).^(0.5);

slope_num
%COMPUTING BOUNCE HEIGHT
for j = 1:length(y)
    if properties1(7,((slope_num(j)*2))) < properties1(7,((slope_num(j)*2)-1))
        bounce(j,1) = (((properties1(6,((slope_num(j)*2))) - x(j,1))*tan(degtorad(properties1(1,((slope_num(j)*2)-1))))) + properties(7,((slope_num(j)*2))));
    else
        bounce(j,1) = properties1(7,((slope_num(j)*2))) - ((((properties1(6,((slope_num(j)*2))))) - (x(j,1)))*tan(degtorad(properties(1,((slope_num(j)*2)-1)))));
    end   
end


bounce   = bounce((1:j),1);
bounceht = y-bounce;  
   
%COMPUTING ENERGIES AND BOUNCE HEIGHTS
translationKE = ((0.5).*(mass).*((vprojectile).^2))/1000;
rotationKE    = ((0.5).*((0.5).*(mass).*(radius).^2).*(vr).^2)/1000;

%CREATING ARRAYS OF DATA FOR ANALYSIS IN EXCEL
handles.mxrockarray1((1:length(x)),u)           = x;
handles.myrockarray1((1:length(y)),u)           = y;
handles.mvrockarray1((1:length(vprojectile)),u) = vprojectile;
handles.mvxrockarray1((1:length(vx)),u)         = vx;
handles.mvyrockarray1((1:length(vy)),u)         = vy;
handles.mvrrockarray1((1:length(vr)),u)         = vr;
handles.mtranslationalKEarray1((1:length(translationKE)),u) = translationKE;
handles.mrotationalKEarray1((1:length(rotationKE)),u)       = rotationKE;
handles.mbouncehtarray1((1:length(bounceht)),u)             = bounceht; 

%PLOTTING 
plot(x,y)
xlabel('X (m)')
ylabel('Y (m)')
grid on
hold all
guidata(hObject,handles)
end



function iterations_Callback(hObject, eventdata, handles)
% hObject    handle to iterations (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of iterations as text
%        str2double(get(hObject,'String')) returns contents of iterations as a double


% --- Executes during object creation, after setting all properties.
function iterations_CreateFcn(hObject, eventdata, handles)
% hObject    handle to iterations (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function t_step_Callback(hObject, eventdata, handles)
% hObject    handle to t_step (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of t_step as text
%        str2double(get(hObject,'String')) returns contents of t_step as a double


% --- Executes during object creation, after setting all properties.
function t_step_CreateFcn(hObject, eventdata, handles)
% hObject    handle to t_step (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function g_Callback(hObject, eventdata, handles)
% hObject    handle to g (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of g as text
%        str2double(get(hObject,'String')) returns contents of g as a double


% --- Executes during object creation, after setting all properties.
function g_CreateFcn(hObject, eventdata, handles)
% hObject    handle to g (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function vscaling_Callback(hObject, eventdata, handles)
% hObject    handle to vscaling (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of vscaling as text
%        str2double(get(hObject,'String')) returns contents of vscaling as a double


% --- Executes during object creation, after setting all properties.
function vscaling_CreateFcn(hObject, eventdata, handles)
% hObject    handle to vscaling (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function radius_Callback(hObject, eventdata, handles)
% hObject    handle to radius (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of radius as text
%        str2double(get(hObject,'String')) returns contents of radius as a double


% --- Executes during object creation, after setting all properties.
function radius_CreateFcn(hObject, eventdata, handles)
% hObject    handle to radius (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function mass_Callback(hObject, eventdata, handles)
% hObject    handle to mass (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of mass as text
%        str2double(get(hObject,'String')) returns contents of mass as a double


% --- Executes during object creation, after setting all properties.
function mass_CreateFcn(hObject, eventdata, handles)
% hObject    handle to mass (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function xi_Callback(hObject, eventdata, handles)
% hObject    handle to xi (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of xi as text
%        str2double(get(hObject,'String')) returns contents of xi as a double


% --- Executes during object creation, after setting all properties.
function xi_CreateFcn(hObject, eventdata, handles)
% hObject    handle to xi (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function num_Callback(hObject, eventdata, handles)
% hObject    handle to num (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of num as text
%        str2double(get(hObject,'String')) returns contents of num as a double


% --- Executes during object creation, after setting all properties.
function num_CreateFcn(hObject, eventdata, handles)
% hObject    handle to num (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function yi_Callback(hObject, eventdata, handles)
% hObject    handle to yi (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of yi as text
%        str2double(get(hObject,'String')) returns contents of yi as a double


% --- Executes during object creation, after setting all properties.
function yi_CreateFcn(hObject, eventdata, handles)
% hObject    handle to yi (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function xistd_Callback(hObject, eventdata, handles)
% hObject    handle to xistd (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of xistd as text
%        str2double(get(hObject,'String')) returns contents of xistd as a double


% --- Executes during object creation, after setting all properties.
function xistd_CreateFcn(hObject, eventdata, handles)
% hObject    handle to xistd (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function yistd_Callback(hObject, eventdata, handles)
% hObject    handle to yistd (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of yistd as text
%        str2double(get(hObject,'String')) returns contents of yistd as a double


% --- Executes during object creation, after setting all properties.
function yistd_CreateFcn(hObject, eventdata, handles)
% hObject    handle to yistd (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function vxi_Callback(hObject, eventdata, handles)
% hObject    handle to vxi (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of vxi as text
%        str2double(get(hObject,'String')) returns contents of vxi as a double


% --- Executes during object creation, after setting all properties.
function vxi_CreateFcn(hObject, eventdata, handles)
% hObject    handle to vxi (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function vyi_Callback(hObject, eventdata, handles)
% hObject    handle to vyi (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of vyi as text
%        str2double(get(hObject,'String')) returns contents of vyi as a double


% --- Executes during object creation, after setting all properties.
function vyi_CreateFcn(hObject, eventdata, handles)
% hObject    handle to vyi (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function vxistd_Callback(hObject, eventdata, handles)
% hObject    handle to vxistd (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of vxistd as text
%        str2double(get(hObject,'String')) returns contents of vxistd as a double


% --- Executes during object creation, after setting all properties.
function vxistd_CreateFcn(hObject, eventdata, handles)
% hObject    handle to vxistd (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function vyistd_Callback(hObject, eventdata, handles)
% hObject    handle to vyistd (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of vyistd as text
%        str2double(get(hObject,'String')) returns contents of vyistd as a double


% --- Executes during object creation, after setting all properties.
function vyistd_CreateFcn(hObject, eventdata, handles)
% hObject    handle to vyistd (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function vri_Callback(hObject, eventdata, handles)
% hObject    handle to vri (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of vri as text
%        str2double(get(hObject,'String')) returns contents of vri as a double


% --- Executes during object creation, after setting all properties.
function vri_CreateFcn(hObject, eventdata, handles)
% hObject    handle to vri (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function vrstd_Callback(hObject, eventdata, handles)
% hObject    handle to vrstd (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of vrstd as text
%        str2double(get(hObject,'String')) returns contents of vrstd as a double


% --- Executes during object creation, after setting all properties.
function vrstd_CreateFcn(hObject, eventdata, handles)
% hObject    handle to vrstd (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function rn_Callback(hObject, eventdata, handles)
% hObject    handle to rn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of rn as text
%        str2double(get(hObject,'String')) returns contents of rn as a double


% --- Executes during object creation, after setting all properties.
function rn_CreateFcn(hObject, eventdata, handles)
% hObject    handle to rn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function rnstd_Callback(hObject, eventdata, handles)
% hObject    handle to rnstd (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of rnstd as text
%        str2double(get(hObject,'String')) returns contents of rnstd as a double


% --- Executes during object creation, after setting all properties.
function rnstd_CreateFcn(hObject, eventdata, handles)
% hObject    handle to rnstd (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function rt_Callback(hObject, eventdata, handles)
% hObject    handle to rt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of rt as text
%        str2double(get(hObject,'String')) returns contents of rt as a double


% --- Executes during object creation, after setting all properties.
function rt_CreateFcn(hObject, eventdata, handles)
% hObject    handle to rt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function rtstd_Callback(hObject, eventdata, handles)
% hObject    handle to rtstd (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of rtstd as text
%        str2double(get(hObject,'String')) returns contents of rtstd as a double


% --- Executes during object creation, after setting all properties.
function rtstd_CreateFcn(hObject, eventdata, handles)
% hObject    handle to rtstd (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function rollingfriction_Callback(hObject, eventdata, handles)
% hObject    handle to rollingfriction (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of rollingfriction as text
%        str2double(get(hObject,'String')) returns contents of rollingfriction as a double


% --- Executes during object creation, after setting all properties.
function rollingfriction_CreateFcn(hObject, eventdata, handles)
% hObject    handle to rollingfriction (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in plottype.
function plottype_Callback(hObject, eventdata, handles)
% hObject    handle to plottype (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns plottype contents as cell array
%        contents{get(hObject,'Value')} returns selected item from plottype


% --- Executes during object creation, after setting all properties.
function plottype_CreateFcn(hObject, eventdata, handles)
% hObject    handle to plottype (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function xcoordinateinput_Callback(hObject, eventdata, handles)
% hObject    handle to xcoordinateinput (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of xcoordinateinput as text
%        str2double(get(hObject,'String')) returns contents of xcoordinateinput as a double


% --- Executes during object creation, after setting all properties.
function xcoordinateinput_CreateFcn(hObject, eventdata, handles)
% hObject    handle to xcoordinateinput (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function ycoordinateinput_Callback(hObject, eventdata, handles)
% hObject    handle to ycoordinateinput (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ycoordinateinput as text
%        str2double(get(hObject,'String')) returns contents of ycoordinateinput as a double


% --- Executes during object creation, after setting all properties.
function ycoordinateinput_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ycoordinateinput (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function numberofcoordinates_Callback(hObject, eventdata, handles)
% hObject    handle to numberofcoordinates (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of numberofcoordinates as text
%        str2double(get(hObject,'String')) returns contents of numberofcoordinates as a double




% --- Executes during object creation, after setting all properties.
function numberofcoordinates_CreateFcn(hObject, eventdata, handles)
% hObject    handle to numberofcoordinates (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmenu1.
function popupmenu1_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu1


% --- Executes during object creation, after setting all properties.
function popupmenu1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in data_to_display.
function data_to_display_Callback(hObject, eventdata, handles)
% hObject    handle to data_to_display (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns data_to_display contents as cell array
%        contents{get(hObject,'Value')} returns selected item from data_to_display


% --- Executes during object creation, after setting all properties.
function data_to_display_CreateFcn(hObject, eventdata, handles)
% hObject    handle to data_to_display (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in manualinput.
function manualinput_Callback(hObject, eventdata, handles)
% hObject    handle to manualinput (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

howmany    = str2num(get(handles.numberofcoordinates,'String'));
xpoint1    = str2num(get(handles.xcoordinateinput,'String'));
ypoint1    = str2num(get(handles.ycoordinateinput,'String'));
iterations = str2num(get(handles.iterations,'String'));


%FITTING DATA THROUGH POINTS OF SLOPE
num_point = length(xpoint1);
m_slope   = zeros((num_point-1),1);
b_slope   = zeros((num_point-1),1);
div       = xpoint1(end)-xpoint1(1);
%CALCULATING ANGLES FOR EACH CELL
angle     = zeros((num_point-1),1);
angle_num = zeros((num_point-1),1);

for q = 1:(num_point-1)
    angle(q,1)     = radtodeg(atan((abs(ypoint1(q)-ypoint1(q+1)))/(abs(xpoint1(q+1)-xpoint1(q)))));
    angle_num(q,1) = q;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%FITTING LINES AND CALCULATING SLOPE GEOMETRY %%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Y = Mx + B FORM 
for w = 1:(num_point-1)
    fit          = polyfit((xpoint1(w:(w+1))),(ypoint1(w:(w+1))),1);        %Polyfit to find m and b
    m1_slope     = fit(1);                    
    b1_slope     = fit(2);
    m_slope(1,w) = m1_slope;                                              %Creating m matrix
    b_slope(1,w) = b1_slope;                                              %Creating b matrix
end

%PRE-ALLOCATING FOR PLOT
x_slope = zeros(iterations,1);
y_slope = zeros(iterations,1);

for e = 1:(num_point-1)
    iterat_segment      = round(((xpoint1((e+1),1)-xpoint1(e,1))/div)*iterations);     %Assigning number of iterations per segment
    iterat_segment1(1,e)= iterat_segment;                                            %Storing number of iterations per segment in a matrix
    x_segment           = linspace(xpoint1(e,1),xpoint1((e+1),1),iterat_segment);      %Calculating x for segment
    y_segment           = ((m_slope(1,e).*x_segment) + b_slope(1,e));                %Calculating y for segment
    x_slope(((sum(iterat_segment1(1,(1:e)))-iterat_segment)+1):(sum(iterat_segment1(1,(1:e)))))     = x_segment;     %Concatenating x for segments together
    y_slope(((sum(iterat_segment1(1,(1:e)))-iterat_segment)+1):(sum(iterat_segment1(1,(1:e)))))     = y_segment;     %Concatenating y for segments together
    xyslope_num(((sum(iterat_segment1(1,(1:e)))-iterat_segment)+1):(sum(iterat_segment1(1,(1:e))))) = e;             %Tracking what segment each x and y belongs to
end

%FLIPPING DATA TO COLUMNS
iterat_segment1   = iterat_segment1';
xyslope_num       = xyslope_num';
%REMOVING 0 ELEMENTS
x_slope(x_slope==0) = [];
y_slope(y_slope==0) = [];
%RE-CALCULATING NUMBER OF ITERATIONS AND PLOTTING PROFILE

handles.iterations1  = length(x_slope);
handles.xpoint1      = xpoint1;
handles.ypoint1      = ypoint1;
guidata(hObject,handles)
plot(x_slope,y_slope,'k')
hold on
grid on   
guidata(hObject,handles);


% --- Executes on button press in graphicalinput.
function graphicalinput_Callback(hObject, eventdata, handles)
% hObject    handle to graphicalinput (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

iterations = str2num(get(handles.iterations,'String'));
[xpoint,ypoint] = ginput;

%FITTING DATA THROUGH POINTS OF SLOPE
num_point = length(xpoint);
m_slope   = zeros((num_point-1),1);
b_slope   = zeros((num_point-1),1);
div       = xpoint(end)-xpoint(1);
%CALCULATING ANGLES FOR EACH CELL
angle     = zeros((num_point-1),1);
angle_num = zeros((num_point-1),1);

for q = 1:(num_point-1)
    angle(q,1)     = radtodeg(atan((abs(ypoint(q)-ypoint(q+1)))/(abs(xpoint(q+1)-xpoint(q)))));
    angle_num(q,1) = q;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%FITTING LINES AND CALCULATING SLOPE GEOMETRY %%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Y = Mx + B FORM 
for w = 1:(num_point-1)
    fit          = polyfit((xpoint(w:(w+1))),(ypoint(w:(w+1))),1);        %Polyfit to find m and b
    m1_slope     = fit(1);                    
    b1_slope     = fit(2);
    m_slope(1,w) = m1_slope;                                              %Creating m matrix
    b_slope(1,w) = b1_slope;                                              %Creating b matrix
end

%PRE-ALLOCATING FOR PLOT
x_slope = zeros(iterations,1);
y_slope = zeros(iterations,1);

for e = 1:(num_point-1)
    iterat_segment      = round(((xpoint((e+1),1)-xpoint(e,1))/div)*iterations);     %Assigning number of iterations per segment
    iterat_segment1(1,e)= iterat_segment;                                            %Storing number of iterations per segment in a matrix
    x_segment           = linspace(xpoint(e,1),xpoint((e+1),1),iterat_segment);      %Calculating x for segment
    y_segment           = ((m_slope(1,e).*x_segment) + b_slope(1,e));                %Calculating y for segment
    x_slope(((sum(iterat_segment1(1,(1:e)))-iterat_segment)+1):(sum(iterat_segment1(1,(1:e)))))     = x_segment;     %Concatenating x for segments together
    y_slope(((sum(iterat_segment1(1,(1:e)))-iterat_segment)+1):(sum(iterat_segment1(1,(1:e)))))     = y_segment;     %Concatenating y for segments together
    xyslope_num(((sum(iterat_segment1(1,(1:e)))-iterat_segment)+1):(sum(iterat_segment1(1,(1:e))))) = e;             %Tracking what segment each x and y belongs to
end

%FLIPPING DATA TO COLUMNS
iterat_segment1   = iterat_segment1';
xyslope_num       = xyslope_num';
%REMOVING 0 ELEMENTS
x_slope(x_slope==0) = [];
y_slope(y_slope==0) = [];
%RE-CALCULATING NUMBER OF ITERATIONS AND PLOTTING PROFILE


iterations1         = length(x_slope);
handles.iterations2 = iterations1;
handles.xpoint2     = xpoint;
handles.ypoint2     = ypoint;
guidata(hObject,handles)

plot(x_slope,y_slope,'k')
hold on
grid on   

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%END FITTING LINES AND CALCULATING SLOPE GEOMETRY %%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% --- Executes on button press in pushbutton4.
function pushbutton4_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
cla,clc
xlim([0, 100])
ylim([0, 100])




% --- Executes on button press in computegraphical.
function computegraphical_Callback(hObject, eventdata, handles)
% hObject    handle to computegraphical (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%RETREIVE DATA FROM INPUTS%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

handles.gxrockarray2  = [];
handles.gyrockarray2  = [];
handles.gvrockarray2  = [];
handles.gvxrockarray2 = [];
handles.gvyrockarray2 = [];
handles.gvrrockarray2 = [];
handles.gtranslationalKEarray2  = [];
handles.grotationalKEarray2     = [];
handles.gbouncehtarray2         = [];

iterations1  = handles.iterations2;
xpoint1      = handles.xpoint2;
ypoint1      = handles.ypoint2;

xpoint = xpoint1;
ypoint = ypoint1;

t_step     = str2num(get(handles.t_step,'String'));
g          = str2num(get(handles.g,'String'));
vfactor    = str2num(get(handles.vscaling,'String'));

radius     = str2num(get(handles.radius,'String'));
mass       = str2num(get(handles.mass,'String'));
num        = str2num(get(handles.num,'String'));

xi         = str2num(get(handles.xi,'String'));
xstd       = str2num(get(handles.xistd,'String'));
yi         = str2num(get(handles.yi,'String'));
ystd       = str2num(get(handles.yistd,'String'));

vxi        = str2num(get(handles.vxi,'String'));
vxstd      = str2num(get(handles.vxistd,'String'));
vyi        = str2num(get(handles.vyi,'String'));
vystd      = str2num(get(handles.vyistd,'String'));
vrot       = str2num(get(handles.vri,'String'));
vrotstd    = str2num(get(handles.vrstd,'String'));

howmany    = str2num(get(handles.numberofcoordinates,'String'));

rt         = str2num(get(handles.rt,'String'));
rt_std     = str2num(get(handles.rtstd,'String'));
rn         = str2num(get(handles.rn,'String'));
rn_std     = str2num(get(handles.rnstd,'String'));
mewr       = str2num(get(handles.rollingfriction,'String'));
properties = zeros(4,(2*(howmany-1)));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%END RETREIVING DATA FROM INPUT%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%CREATING PROPERTIES%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Beggining plotting
num_point = howmany;

%FITTING DATA THROUGH POINTS OF SLOPE
m_slope   = zeros((num_point-1),1);
b_slope   = zeros((num_point-1),1);
div       = xpoint(end)-xpoint(1);
%CALCULATING ANGLES FOR EACH CELL
angle     = zeros((num_point-1),1);
angle_num = zeros((num_point-1),1);

for q = 1:(num_point-1)
    angle(q,1)     = radtodeg(atan((abs(ypoint(q)-ypoint(q+1)))/(abs(xpoint(q+1)-xpoint(q)))));
    angle_num(q,1) = q;
end

for r = 1:(num_point-1)      
    properties(1,((r*2)-1)) = angle(r);
    properties(2,((r*2)-1)) = mewr(r,1);
    properties(3,((r*2)-1)) = rn(r,1);
    properties(4,((r*2)-1)) = rt(r,1); 
    properties(5,((r*2)-1)) = r;
    properties(6,((r*2)-1)) = xpoint(r,1);
    properties(7,((r*2)-1)) = ypoint(r,1);
    properties(8,((r*2)-1)) = ((abs(g))*(cos(degtorad(angle(r)))*(mewr(r,1))));
    
    properties(1,(r*2)) = 0;
    properties(2,(r*2)) = 0;
    properties(3,(r*2)) = rn_std(r,1);
    properties(4,(r*2)) = rt_std(r,1);
    properties(5,(r*2)) = 0;
    properties(6,(r*2)) = xpoint((r+1),1);
    properties(7,(r*2)) = ypoint((r+1),1);
    properties(8,(r*2)) = 0;
end
properties;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%END CREATING PROPERTIES%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%START CALCULATION LOOPS%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for u = 1:num
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%START SETUP OF INITIAL PARAMETERS%%%%%%%%%%%%%%%%
    xi1   = xi+ xstd*(abs(randn(1,1)));    
    yi1   = yi+ ystd*(abs(randn(1,1)));                   
    vxi1  = (vxi+ vxstd*randn(1,1));
    vyi1  = (vyi+ vystd*randn(1,1));
    vrot1 = (vrot+vrotstd*randn(1,1));
    properties1 = properties;      
    
    for i = 1:(num_point-1)    
        properties1(3,((i*2)-1)) = properties1(3,((i*2)-1)) +  (properties1(3,(i*2)) * randn(1,1));
        properties1(4,((i*2)-1)) = properties1(4,((i*2)-1)) +  (properties1(4,(i*2)) * randn(1,1));
        properties1(3,(i*2))     = 0;
        properties1(4,(i*2))     = 0;
    end
    
    x  = xi1;  
    y  = yi1;
    vx = vxi1; 
    vy = vyi1;
    vr = vrot1;
    
    for h = 1:(num_point-1)
         if  x(1,1) >= properties1(6,((h*2)-1)) &  x(1,1) <= properties1(6,(h*2))
         slope_num(1,1)= h;
         end
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %FINISHED SETUP OF INITIAL PARAMETERS. X AND Y FRAME OF REFERENCE%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%START INDIVIDUAL ROCK CALCULATIONS%%%%%%%%%%%%%%%
   
    for o = 1:iterations1-1
       
        %%%%%%%%%%%%%%%%%PROJECTILE TRAJECTORY CODE%%%%%%%%%%%%%%%%%%%
        vystep     = vy(o,1) + (g)*(t_step);                          
        vy(o+1,1)  = vystep;                                          
        vxstep     = vx(o,1) + (0)*(t_step);                              
        vx(o+1,1)  = vxstep;                       
        vr(o+1,1)  = vr(o,1);
        xstep      = x(o,1)+ (vx(o,1)*t_step);                             
        x(o+1,1)   = xstep;        
        ystep      = y(o,1) + ((vy(o+1,1) + vy(o,1))*(t_step)*0.5);
        y(o+1,1)   = ystep;
        %%%%%%%%%%%%END PROJECTILE TRAJECTORY CODE%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %%%%%%%%%%%%%%%%SLOPE LOCATION BELOW ROCK CODE%%%%%%%%%%%%%%%
        for p = 1:(num_point-1)
            
            if  x(o+1,1) >= properties1(6,((p*2)-1)) &  x(o+1,1) <= properties1(6,(p*2))
                if  properties1(7,(p*2)) <  properties1(7,((p*2)-1))
                    y_impact_location = (properties1(6,(p*2))- x(o+1,1)) * (tan(degtorad(properties1(1,((p*2)-1))))) + properties1(7,(p*2));
                    slope_num(o+1,1)  = p;
                else 
                    y_impact_location = (x(o+1,1)-properties1(6,((p*2)-1))) * (tan(degtorad(properties1(1,((p*2)-1))))) + properties1(7,((p*2)-1));
                    slope_num(o+1,1)  = p;
                end 
                break
            end
            
        end 
        %%%%%%%%%%%END SLOPE LOCATION BELOW ROCK CODE%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %%%%%%%%IMPACT LOCATION AND RESTITUTION CODE%%%%%%%%%%%%%%%%
        if y(o+1,1) - y_impact_location <=0.1
            
            if  properties1(7,(p*2)) <  properties1(7,((p*2)-1))  %Downhill 
                if vx(o+1,1) > 0 %Vx is positive
                    vtb     = abs((vy(o+1,1)*(sin((degtorad(properties1(1,((p*2)-1))))))))+(vx(o+1,1)*(cos((degtorad(properties1(1,((p*2)-1)))))));  
                    vnb     = (vy(o+1,1)*cos((degtorad(properties1(1,((p*2)-1))))))+(vx(o+1,1)*sin((degtorad(properties1(1,((p*2)-1))))));                       
                    vta     = vtb*(properties1(4,((p*2)-1)));                                                 
                    vna     = (abs(vnb))*((properties1(3,((p*2)-1)))/(1+(((abs(vnb)))/vfactor)^2));
                    vx(o,1) = (vna*sin((degtorad(properties1(1,((p*2)-1))))))+(vta*cos((degtorad(properties1(1,((p*2)-1))))));                
                    vy(o,1) = (vna*cos((degtorad(properties1(1,((p*2)-1))))))+((-1)*(vta*sin((degtorad(properties1(1,((p*2)-1)))))));
                    vr(o,1) = ((((vx(o,1))^2 + (vy(o,1))^2)^(0.5))/radius);
                else    %Vx is negative        
                    vtb     = abs(vy(o+1,1)*(sin((degtorad(properties1(1,((p*2)-1)))))))+(vx(o+1,1)*(cos((degtorad(properties1(1,((p*2)-1)))))));  
                    vnb     = (vy(o+1,1)*cos((degtorad(properties1(1,((p*2)-1))))))+(vx(o+1,1)*sin((degtorad(properties1(1,((p*2)-1))))));                       
                    vta     = vtb*(properties1(4,((p*2)-1)));                                                 
                    vna     = (abs(vnb))*((properties1(3,((p*2)-1)))/(1+(((abs(vnb)))/vfactor)^2));
                    vx(o,1) = (vna*sin((degtorad(properties1(1,((p*2)-1))))))+(vta*cos((degtorad(properties1(1,((p*2)-1))))));                
                    vy(o,1) = (vna*cos((degtorad(properties1(1,((p*2)-1))))))+((-1)*(vta*sin((degtorad(properties1(1,((p*2)-1)))))));
                    vr(o,1) = ((((vx(o,1))^2 + (vy(o,1))^2)^(0.5))/radius);
                end     
            else %Uphill
                if vx(o+1,1)>0 %Vx is positive
                    vtb     = (vy(o+1,1)*(sin((degtorad(properties1(1,((p*2)-1)))))))+(vx(o+1,1)*(cos((degtorad(properties1(1,((p*2)-1)))))));   
                    vnb     = (vy(o+1,1)*cos((degtorad(properties1(1,((p*2)-1))))))+((-1)*(vx(o+1,1)*sin((degtorad(properties1(1,((p*2)-1))))))); 
                    vta     = vtb*(properties1(4,((p*2)-1)));
                    vna     = (abs(vnb))*((properties1(3,((p*2)-1)))/(1+(((abs(vnb)))/vfactor)^2));
                    vx(o,1) = ((-1)*(vna*sin((degtorad(properties1(1,((p*2)-1)))))))+(vta*cos((degtorad(properties1(1,((p*2)-1))))));                
                    vy(o,1) = (vna*cos((degtorad(properties1(1,((p*2)-1))))))+((vta*sin((degtorad(properties1(1,((p*2)-1)))))));
                    vr(o,1) = ((((vx(o,1))^2 + (vy(o,1))^2)^(0.5))/radius);
                else %Vx is negative
                    vtb     = (vy(o+1,1)*(sin((degtorad(properties1(1,((p*2)-1)))))))+(vx(o+1,1)*(cos((degtorad(properties1(1,((p*2)-1)))))));   
                    vnb     = (vy(o+1,1)*cos((degtorad(properties1(1,((p*2)-1))))))+abs(vx(o+1,1)*sin((degtorad(properties1(1,((p*2)-1)))))); 
                    vta     = (vtb*(properties1(4,((p*2)-1))));
                    vna     = (abs(vnb))*((properties1(3,((p*2)-1)))/(1+(((abs(vnb)))/vfactor)^2));
                    vx(o,1) = ((-1)*(vna*sin((degtorad(properties1(1,((p*2)-1)))))))+(vta*cos((degtorad(properties1(1,((p*2)-1))))));                
                    vy(o,1) = (vna*cos((degtorad(properties1(1,((p*2)-1))))))+((vta*sin((degtorad(properties1(1,((p*2)-1)))))));
                    vr(o,1) = ((((vx(o,1))^2 + (vy(o,1))^2)^(0.5))/radius);            
                end
                
            end                           
            %%%%%%%%END IMPACT LOCATION AND RESTITUTION CODE%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
        
           if vna >0.1
               %%%%%%%%TRAJECTORY AFTER IMPACT CODE%%%%%%%%%%%%%%%%%%%%%
               for a =  o:(iterations1-1)
                   vyf         = vy(a,1) + (g)*(t_step);                               
                   vy(a+1,1)   = vyf;                                                 
                   vxf         = vx(a,1)+(0)*(t_step);                                     
                   vx(a+1,1)   = vxf;
                   vr(a+1,1)   = vr(a,1);
                   xf          = x(a,1)+ (vx(a,1)*t_step);                                 
                   x(a+1,1)    = xf;
           
                   if properties1(7,(p*2)) <  properties1(7,((p*2)-1))
                       yf = y(a,1) + ((vy(a+1,1) + vy(a,1))* (0.5) * (t_step));
                   else
                       yf = y(a,1) + abs((vy(a+1,1) + vy(a,1))* (0.5) * (t_step));
                   end    
                   y(a+1,1) = yf;  
           
                   for p = 1:(num_point-1)
               
                       if  x(a+1,1) >= properties1(6,((p*2)-1)) &  x(a+1,1) <= properties1(6,(p*2))
                            if  properties1(7,(p*2)) <  properties1(7,((p*2)-1))
                                y_impact_location = (properties1(6,(p*2))- x(a+1,1)) * (tan(degtorad(properties1(1,((p*2)-1))))) + properties1(7,(p*2));
                                slope_num(a+1,1)  = p;
                            else 
                                y_impact_location = (x(a+1,1)-properties1(6,((p*2)-1))) * (tan(degtorad(properties1(1,((p*2)-1))))) + properties1(7,((p*2)-1));
                                slope_num(a+1,1)  = p;
                            end 
                        break
                       end
                   end
           
                   if vxf < 0.1 | vyf < 0.1
                       break
                   end
           
               end 
           else  
               %%%%%%%%%%%SLIDING AFTER IMPACT CODE%%%%%%%%%%%%%%%%%%%%%
               vxsmatrix(1,1) = vta; 
               for s = o:iterations1-1  
 
                   c = (1+(s-o));
                   n = (2+(s-o));
  
                   if vxsmatrix(c,1)>0
                       if  properties1(7,(p*2)) <  properties1(7,((p*2)-1))
                       part1 = (0.75)*(vxsmatrix(c,1))^2;
                       part2 = (((abs(g))*(vxsmatrix(c,1))*(t_step)) * ((cos(degtorad(properties(1,((p*2)-1)))) * properties(2,((p*2)-1))) - (sin(degtorad(properties(1,((p*2)-1)))))));
                       vxsmatrix(n,1) = (abs((4*(part1-part2))/3))^0.5;
                           if vxsmatrix(n,1)<0.1
                               vxsmatrix = vxsmatrix((1:c),1);
                               vr = vr((1:(s)),1);
                               x  = x((1:(s)),1);
                               y  = y((1:(s)),1);
                               slope_num = slope_num((1:(s)),1);
                               break  
                           end
                       else
                       part1 = (0.75)*(vxsmatrix(c,1))^2;
                       part2 = (((abs(g))*(vxsmatrix(c,1))*(t_step)) * ((sin(degtorad(properties(1,((p*2)-1)))))  +(cos(degtorad(properties(1,((p*2)-1)))) * properties(2,((p*2)-1)))));
                       vxsmatrix(n,1) = (abs((4*(part1-part2))/3))^0.5;
                           if vxsmatrix(n,1)<0.1
                               vxsmatrix = vxsmatrix((1:c),1);
                               vr = vr((1:(s)),1);
                               x  = x((1:(s)),1);
                               y  = y((1:(s)),1);
                               slope_num = slope_num((1:(s)),1);
                               break  
                           end  
                       end
                   else
                       if  properties1(7,(p*2)) <  properties1(7,((p*2)-1))
                       part1 = (0.75)*(vxsmatrix(c,1))^2;
                       part2 = (((abs(g))*(vxsmatrix(c,1))*(t_step)) * ((sin(degtorad(properties(1,((p*2)-1)))))  +(cos(degtorad(properties(1,((p*2)-1)))) * properties(2,((p*2)-1)))));
                       vxsmatrix(n,1) = (abs((4*(part1-part2))/3))^0.5;
                           if vxsmatrix(n,1)>-0.1
                               vxsmatrix = vxsmatrix((1:c),1);
                               vr = vr((1:(s)),1);
                               x  = x((1:(s)),1);
                               y  = y((1:(s)),1);
                               slope_num = slope_num((1:(s)),1);
                               break
                           end
                       else                       
                       part1 = (0.75)*(vxsmatrix(c,1))^2;
                       part2 = (((abs(g))*(vxsmatrix(c,1))*(t_step)) * ((cos(degtorad(properties(1,((p*2)-1)))) * properties(2,((p*2)-1))) - (sin(degtorad(properties(1,((p*2)-1)))))));
                       vxsmatrix(n,1) = (abs((4*(part1-part2))/3))^0.5;
                           if vxsmatrix(n,1)>-0.1
                               vxsmatrix = vxsmatrix((1:c),1);
                               vr = vr((1:(s)),1);
                               x  = x((1:(s)),1);
                               y  = y((1:(s)),1);
                               slope_num = slope_num((1:(s)),1);
                               break
                           end
                       end
                   end

           
                   distanceslid      = t_step*vxsmatrix(c,1);
                   xdistanceslid     = (distanceslid*(cos((degtorad(properties1(1,((p*2)-1)))))));      
                   x(s+1,1)          = x(s)+xdistanceslid;  
                   ydistanceslid     = (distanceslid*(sin((degtorad(properties1(1,((p*2)-1)))))));
           
                   %CHECKING IF CELL IS ASCENDING OR DESCENDING
                   if  properties1(7,(p*2)) <  properties1(7,((p*2)-1))
                       y(s+1,1) = y(s)-ydistanceslid;
                   else 
                       y(s+1,1) = y(s)+ydistanceslid;
                   end 
                   slope_num(s+1,1) = p;    
                   
                   if  properties1(7,(p*2)) <  properties1(7,((p*2)-1))
                       if vxsmatrix(n,1)>0
                       vx(s+1,1) = (vxsmatrix(n,1)*cos((degtorad(properties1(1,((p*2)-1))))));                
                       vy(s+1,1) = -(vxsmatrix(n,1)*sin((degtorad(properties1(1,((p*2)-1))))));
                       vr(s+1,1) = ((((vx(s,1))^2 + (vy(s,1))^2)^(0.5))/radius);
                       else
                       vx(s+1,1) = (vxsmatrix(n,1)*cos((degtorad(properties1(1,((p*2)-1))))));                
                       vy(s+1,1) = (vxsmatrix(n,1)*sin((degtorad(properties1(1,((p*2)-1))))));
                       vr(s+1,1) = ((((vx(s,1))^2 + (vy(s,1))^2)^(0.5))/radius);
                       end
                   else
                       if vxsmatrix(n,1)>0
                       vx(s+1,1) = (vxsmatrix(n,1)*cos((degtorad(properties1(1,((p*2)-1))))));                
                       vy(s+1,1) = (vxsmatrix(n,1)*sin((degtorad(properties1(1,((p*2)-1))))));
                       vr(s+1,1) = ((((vx(s,1))^2 + (vy(s,1))^2)^(0.5))/radius);
                       else
                       vx(s+1,1) = (vxsmatrix(n,1)*cos((degtorad(properties1(1,((p*2)-1))))));                
                       vy(s+1,1) = -(vxsmatrix(n,1)*sin((degtorad(properties1(1,((p*2)-1))))));
                       vr(s+1,1) = ((((vx(s,1))^2 + (vy(s,1))^2)^(0.5))/radius);  
                       end
                   end
                   %CHECKING IF ROCK HAS REACHED END OF CELL
                   if vxsmatrix(n,1)>0
                       if x(s+1,1) > properties1(6,(p*2))
                           p = p+1;
                           if p == 0 
                               break
                           end
                           if p>(num_point-1)
                               break
                           end
                       end
                   else
                       if x(s+1,1) < properties1(6,((p*2)-1))
                           p = p-1;          
                           if p == 0 
                               break
                           end
                           if p>(num_point-1)
                               break
                           end
                       end
                   end
               end   
               break
           end    
        end
    end
vprojectile = ((vx).^2+(vy).^2).^(0.5);
%COMPUTING BOUNCE HEIGHT
for j = 1:length(y)
    if properties1(7,((slope_num(j)*2))) < properties1(7,((slope_num(j)*2)-1))
        bounce(j,1) = (((properties1(6,((slope_num(j)*2))) - x(j,1))*tan(degtorad(properties1(1,((slope_num(j)*2)-1))))) + properties(7,((slope_num(j)*2))));
    else
        bounce(j,1) = properties1(7,((slope_num(j)*2))) - ((((properties1(6,((slope_num(j)*2))))) - (x(j,1)))*tan(degtorad(properties(1,((slope_num(j)*2)-1)))));
    end   
end

bounce   = bounce((1:j),1);
bounceht = y-bounce;  
   
%COMPUTING ENERGIES AND BOUNCE HEIGHTS
translationKE = ((0.5).*(mass).*((vprojectile).^2))/1000;
rotationKE    = ((0.5).*((0.5).*(mass).*(radius).^2).*(vr).^2)/1000;

%CREATING ARRAYS OF DATA FOR ANALYSIS IN EXCEL
handles.gxrockarray2((1:length(x)),u)           = x;
handles.gyrockarray2((1:length(y)),u)           = y;
handles.gvrockarray2((1:length(vprojectile)),u) = vprojectile;
handles.gvxrockarray2((1:length(vx)),u)         = vx;
handles.gvyrockarray2((1:length(vy)),u)         = vy;
handles.gvrrockarray2((1:length(vr)),u)         = vr;
handles.gtranslationalKEarray2((1:length(translationKE)),u) = translationKE;
handles.grotationalKEarray2((1:length(rotationKE)),u)       = rotationKE;
handles.gbouncehtarray2((1:length(bounceht)),u)             = bounceht; 

%PLOTTING 
plot(x,y)
xlabel('X (m)')
ylabel('Y (m)')
grid on
hold all

guidata(hObject,handles)
end



% --- Executes on button press in plotdropdownselection.
function plotdropdownselection_Callback(hObject, eventdata, handles)
% hObject    handle to plotdropdownselection (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

num     = str2num(get(handles.num,'String'));
string1 = get(handles.plottype,'Value');
string2 = get(handles.plottype,'String');
string  = string2(string1);

s1 = 'Rock Velocity - Graphical';
s2 = 'Rock X Velocity - Graphical';
s0 = 'Rock Y Velocity - Graphical';
s3 = 'Rock Rotational Velocity - Graphical';
s5 = 'Rock Kinetic Energy - Graphical';
s7 = 'Rock Rotational Energy - Graphical';
s8 = 'Rock Total Energy - Graphical'
s9 = 'Rock Bounce Height - Graphical';


gxrockarray   = handles.gxrockarray2;
gvrockarray   = handles.gvrockarray2;
gvxrockarray  = handles.gvxrockarray2;
gvyrockarray  = handles.gvyrockarray2;
gvrrockarray  = handles.gvrrockarray2;
gtranslationalKEarray = handles.gtranslationalKEarray2;
grotationalKEarray    = handles.grotationalKEarray2;
gtotalenergy          = gtranslationalKEarray + grotationalKEarray;
gbouncehtarray        = handles.gbouncehtarray2;

iterations  = size(gxrockarray);
iterations1 = iterations(1);


if strcmp(string,s1)==1
    for yy = 1:num
        for zz = 1:iterations1
            if gxrockarray(zz,yy)==0
            break
            end
        end
    figure(1)
    plot(gxrockarray(1:(zz-1),yy),gvrockarray(1:(zz-1),yy))
    xlabel('X (m)')
    ylabel('Translational Velocity (m/s)')
    grid on
    hold all
    end
elseif strcmp(string,s2)==1
    for yy = 1:num
        for zz = 1:iterations1
            if gxrockarray(zz,yy)==0
            break
            end
        end
    figure(2)
    plot(gxrockarray(1:(zz-1),yy),gvxrockarray(1:(zz-1),yy))
    xlabel('X (m)')
    ylabel('Translational X Velocity (m/s)')
    grid on
    hold all
    end    
elseif strcmp(string,s0)==1
    for yy = 1:num
        for zz = 1:iterations1
            if gxrockarray(zz,yy)==0
            break
            end
        end
    figure(3)
    plot(gxrockarray(1:(zz-1),yy),gvyrockarray(1:(zz-1),yy))
    xlabel('X (m)')
    ylabel('Translational Y Velocity (m/s)')
    grid on
    hold all
    end    
elseif strcmp(string,s3)==1
    for yy = 1:num  
        for zz = 1:iterations1
            if gxrockarray(zz,yy)==0
            break
            end
        end
    figure(4)
    plot(gxrockarray(1:(zz-1),yy),gvrrockarray(1:(zz-1),yy))
    xlabel('X (m)')
    ylabel('Rotational velocty (rad/s)')
    grid on
    hold all
    end
elseif strcmp(string,s5)==1
    for yy = 1:num  
        for zz = 1:iterations1
            if gxrockarray(zz,yy)==0
            break
            end
        end
    figure(5)
    plot(gxrockarray(1:(zz-1),yy),gtranslationalKEarray(1:(zz-1),yy))
    xlabel('X (m)')
    ylabel('Kinetic Energy (kJ)')
    grid on
    hold all
    end
elseif strcmp(string,s7)==1
    for yy = 1:num   
        for zz = 1:iterations1
            if gxrockarray(zz,yy)==0
            break
            end
        end
    figure(6)
    plot(gxrockarray(1:(zz-1),yy),grotationalKEarray(1:(zz-1),yy)) 
    xlabel('X (m)')
    ylabel('Rotational Energy (kJ)')
    grid on
    hold all
    end
elseif strcmp(string,s8)==1
    for yy = 1:num   
        for zz = 1:iterations1
            if gxrockarray(zz,yy)==0
            break
            end
        end
    figure(7)
    plot(gxrockarray(1:(zz-1),yy),gtotalenergy(1:(zz-1),yy)) 
    xlabel('X (m)')
    ylabel('Total Energy (kJ)')
    grid on
    hold all
    end
elseif strcmp(string,s9)==1
    for yy = 1:num
        for zz = 1:iterations1
            if gxrockarray(zz,yy)==0
            break
            end
        end
    figure(8)
    plot(gxrockarray(1:(zz-1),yy),gbouncehtarray(1:(zz-1),yy))
    xlabel('X (m)')
    ylabel('Bounce Height Above Slope (m)')
    grid on
    hold all
    end
end


% --- Executes on button press in plotmanualdropdownselection.
function plotmanualdropdownselection_Callback(hObject, eventdata, handles)
% hObject    handle to plotmanualdropdownselection (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

num      = str2num(get(handles.num,'String'));
string1  = get(handles.plottype,'Value');
string2  = get(handles.plottype,'String');
string   = string2(string1);


s2  = 'Rock Velocity - Manual';
s1  = 'Rock X Velocity - Manual';
s3  = 'Rock Y Velocity - Manual';
s4  = 'Rock Rotational Velocity - Manual';
s6  = 'Rock Kinetic Energy - Manual';
s8  = 'Rock Rotational Energy - Manual';
s9  = 'Rock Total Energy - Manual'
s10 = 'Rock Bounce Height - Manual';

mxrockarray   = handles.mxrockarray1;
mvrockarray   = handles.mvrockarray1;
mvxrockarray  = handles.mvxrockarray1;
mvyrockarray  = handles.mvyrockarray1;
mvrrockarray  = handles.mvrrockarray1;
mtranslationalKEarray  = handles.mtranslationalKEarray1;
mrotationalKEarray     = handles.mrotationalKEarray1;
mtotalenergy           = mtranslationalKEarray + mrotationalKEarray;
mbouncehtarray         = handles.mbouncehtarray1;

iterations  = size(mxrockarray);
iterations1 = iterations(1);


if strcmp(string,s2)==1
    for yy = 1:num   
        for zz = 1:iterations1
            if mxrockarray(zz,yy)==0
            break
            end
        end
    figure(1)
    plot(mxrockarray(1:(zz-1),yy),mvrockarray(1:(zz-1),yy))
    xlabel('X (m)')
    ylabel('Translational Velocity (m/s)')
    grid on
    hold all
    end
elseif strcmp(string,s1)==1
    for yy = 1:num   
        for zz = 1:iterations1
            if mxrockarray(zz,yy)==0
            break
            end
        end
    figure(2)
    plot(mxrockarray(1:(zz-1),yy),mvxrockarray(1:(zz-1),yy))
    xlabel('X (m)')
    ylabel('Translational X Velocity (m/s)')
    grid on
    hold all
    end        
elseif strcmp(string,s3)==1
    for yy = 1:num   
        for zz = 1:iterations1
            if mxrockarray(zz,yy)==0
            break
            end
        end
    figure(3)
    plot(mxrockarray(1:(zz-1),yy),mvyrockarray(1:(zz-1),yy))
    xlabel('X (m)')
    ylabel('Translational Y Velocity (m/s)')
    grid on
    hold all
    end 
elseif strcmp(string,s4)==1
    for yy = 1:num
        for zz = 1:iterations1
            if mxrockarray(zz,yy)==0
            break
            end
        end
    figure(4)
    plot(mxrockarray(1:(zz-1),yy),mvrrockarray(1:(zz-1),yy))
    xlabel('X (m)')
    ylabel('Rotational velocty (rad/s)')
    grid on
    hold all
    end
elseif strcmp(string,s6)==1
    for yy = 1:num
        for zz = 1:iterations1
            if mxrockarray(zz,yy)==0
            break
            end
        end
    figure(5)
    plot(mxrockarray(1:(zz-1),yy),mtranslationalKEarray(1:(zz-1),yy)) 
    xlabel('X (m)')
    ylabel('Kinetic Energy (kJ)')
    grid on
    hold all
    end
elseif strcmp(string,s8)==1
    for yy = 1:num
        for zz = 1:iterations1
            if mxrockarray(zz,yy)==0
            break
            end
        end
    figure(6)
    plot(mxrockarray(1:(zz-1),yy),mrotationalKEarray(1:(zz-1),yy))
    xlabel('X (m)')
    ylabel('Rotational Energy (kJ)')
    grid on
    hold all
    end
elseif strcmp(string,s9)==1
    for yy = 1:num
        for zz = 1:iterations1
            if mxrockarray(zz,yy)==0
            break
            end
        end
    figure(7)
    plot(mxrockarray(1:(zz-1),yy),mtotalenergy(1:(zz-1),yy))
    xlabel('X (m)')
    ylabel('Total Energy (kJ)')
    grid on
    hold all
    end
elseif strcmp(string,s10)==1
    for yy = 1:num 
        for zz = 1:iterations1
            if mxrockarray(zz,yy)==0
            break
            end
        end
    figure(8)
    plot(mxrockarray(1:(zz-1),yy),mbouncehtarray(1:(zz-1),yy))
    xlabel('X (m)')
    ylabel('Bounce Height Above Slope (m)')
    grid on
    hold all
    end
end


% --- Executes on button press in displaymanualchartdata.
function displaymanualchartdata_Callback(hObject, eventdata, handles)
% hObject    handle to displaymanualchartdata (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

mxrockarray  = handles.mxrockarray1;
myrockarray  = handles.myrockarray1;
mvrockarray  = handles.mvrockarray1;
mvxrockarray = handles.mvxrockarray1;
mvyrockarray = handles.mvyrockarray1;
mvrrockarray = handles.mvrrockarray1;
mtranslationalKEarray = handles.mtranslationalKEarray1;
mrotationalKEarray    = handles.mrotationalKEarray1;
mtotalenergy          = mtranslationalKEarray + mrotationalKEarray;
mbouncehtarray        = handles.mbouncehtarray1;

string1 = get(handles.data_to_display,'Value');
string2 = get(handles.data_to_display,'String');
string  = string2(string1);


s1 = 'Rock X Location - Manual';
s2 = 'Rock Y Location - Manual';
s3 = 'Rock Velocity - Manual';
s4 = 'Rock X Velocity - Manual';
s5 = 'Rock Y Velocity - Manual';
s6 = 'Rock Rotational Velocity - Manual';
s7 = 'Rock Kinetic Energy - Manual';
s8 = 'Rock Rotational Energy - Manual';
s0 = 'Rock Total Energy - Manual'
s9 = 'Rock Bounce Height - Manual';

if strcmp(string,s1)==1
    set(handles.uitable2, 'Data', mxrockarray);
    csvwrite('Manual Rock X Locations',mxrockarray)
elseif strcmp(string,s2)==1
    set(handles.uitable2, 'Data', myrockarray);
    csvwrite('Manual Rock Y Locations',myrockarray)
elseif strcmp(string,s3)==1 
    set(handles.uitable2, 'Data', mvrockarray);
    csvwrite('Manual Rock Velocities',mvrockarray)
elseif strcmp(string,s4)==1
    set(handles.uitable2, 'Data', mvxrockarray);
    csvwrite('Manual Rock X Velocities',mvxrockarray)
elseif strcmp(string,s5)==1
    set(handles.uitable2, 'Data', mvyrockarray);
    csvwrite('Manual Rock Y Velocities',mvyrockarray)
elseif strcmp(string,s6)==1
    set(handles.uitable2, 'Data', mvrrockarray);
    csvwrite('Manual Rock Rotational Velocities',mvrrockarray)
elseif strcmp(string,s7)==1
    set(handles.uitable2, 'Data', mtranslationalKEarray);
    csvwrite('Manual Rock Translational Energies',mtranslationalKEarray)
elseif strcmp(string,s8)==1
    set(handles.uitable2, 'Data', mrotationalKEarray);
    csvwrite('Manual Rock Rotational Energies',mrotationalKEarray)
elseif strcmp(string,s0)==1
    set(handles.uitable2, 'Data', mtotalenergy);
    csvwrite('Manual Rock Total Energies',mtotalenergy)
elseif strcmp(string,s9)==1
    set(handles.uitable2, 'Data', mbouncehtarray);
    csvwrite('Manual Rock Bounce Heights',mbouncehtarray)
end


% --- Executes on button press in displaygraphicalchartdata.
function displaygraphicalchartdata_Callback(hObject, eventdata, handles)
% hObject    handle to displaygraphicalchartdata (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


gxrockarray  = handles.gxrockarray2;
gyrockarray  = handles.gyrockarray2;
gvrockarray  = handles.gvrockarray2;
gvxrockarray = handles.gvxrockarray2;
gvyrockarray = handles.gvyrockarray2;
gvrrockarray = handles.gvrrockarray2;
gtranslationalKEarray  = handles.gtranslationalKEarray2;
grotationalKEarray     = handles.grotationalKEarray2;
gtotalenergy           = gtranslationalKEarray + grotationalKEarray;
gbouncehtarray         = handles.gbouncehtarray2;

string1 = get(handles.data_to_display,'Value');
string2 = get(handles.data_to_display,'String');
string  = string2(string1);

s1 = 'Rock X Location - Graphical';
s2 = 'Rock Y Location - Graphical';
s3 = 'Rock Velocity - Graphical';
s4 = 'Rock X Velocity - Graphical';
s5 = 'Rock Y Velocity - Graphical';
s6 = 'Rock Rotational Velocity - Graphical';
s7 = 'Rock Kinetic Energy - Graphical';
s8 = 'Rock Rotational Energy - Graphical';
s0 = 'Rock Total Energy - Graphical'
s9 = 'Rock Bounce Height - Graphical';

if strcmp(string,s1)==1
    set(handles.uitable2, 'Data', gxrockarray);
    csvwrite('Graphical Rock X Locations',gxrockarray)
elseif strcmp(string,s2)==1
    set(handles.uitable2, 'Data', gyrockarray);
    csvwrite('Graphical Rock Y Locations',gyrockarray)
elseif strcmp(string,s3)==1 
    set(handles.uitable2, 'Data', gvrockarray);
    csvwrite('Graphical Rock Velocities',gvrockarray)
elseif strcmp(string,s4)==1
    set(handles.uitable2, 'Data', gvxrockarray);
    csvwrite('Graphical Rock X Velocities',gvxrockarray)
elseif strcmp(string,s5)==1
    set(handles.uitable2, 'Data', gvyrockarray);
    csvwrite('Graphical Rock Y Velocities',gvyrockarray)
elseif strcmp(string,s6)==1
    set(handles.uitable2, 'Data', gvrrockarray);
    csvwrite('Graphical Rock Rotational Velocities',gvrrockarray)
elseif strcmp(string,s7)==1
    set(handles.uitable2, 'Data', gtranslationalKEarray);
    csvwrite('Graphical Rock Translational Energies',gtranslationalKEarray)
elseif strcmp(string,s8)==1
    set(handles.uitable2, 'Data', grotationalKEarray);
    csvwrite('Graphical Rock Rotational Energies',grotationalKEarray)
elseif strcmp(string,s0)==1
    set(handles.uitable2, 'Data', gtotalenergy);
    csvwrite('Graphical Rock Total Energies',gtotalenergy)
elseif strcmp(string,s9)==1
    set(handles.uitable2, 'Data', gbouncehtarray);
    csvwrite('Graphical Rock Bounce Heights',gbouncehtarray)
end


% --- Executes on button press in gridoff.
function gridoff_Callback(hObject, eventdata, handles)
% hObject    handle to gridoff (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
grid off

% --- Executes on button press in gridon.
function gridon_Callback(hObject, eventdata, handles)
% hObject    handle to gridon (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
grid on
