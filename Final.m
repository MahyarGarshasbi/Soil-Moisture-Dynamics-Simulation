clc
clear
close all


Zr         = 0.5 * 1000 ;                          %root depth ( mm )
Kc.initial = 0.6  ;                                %onion
Kc.mid     = 1    ;                                %onion
Kc.dev     = Kc.initial:(Kc.mid-Kc.initial)/45:Kc.mid ;
Kc.end     = 1    ;                                %onion
Ks         = 1000 ;                                %mm/day loam sandy
phi        = 0.42 ;                                %loam sandy
beta       = 12.7 ;                                %loam sandy  ? 
s.h        = 0.08 ;                                %loam sandy
s.w        = 0.11 ;                                %loam sandy
s.star     = 0.31 ;                                %loam sandy
s.fc       = 0.52 ;                                %loam sandy
s.target   = 0.50;
ET.W       = 3    ;
ET.o       = xlsread('data96-97.xlsx','N6:N248') ; %penman-monteith
h          = xlsread('data96-97.xlsx','l6:l248') ; %rainfall
delta      = 0.5;                                  %interception of grass vegetation
h          = max ( 0 , h - delta ) ;
Ky         = 1.1;                                  %yield response factor
area       = 4046.8564224 ;                        %1acre ,(m^2)
yield.p    = 17 ;                                  %yield per acre ( ton )
%potential evapotranspiration
for t=1:243
    if t <= 20
        ET.P(t) = Kc.initial * ET.o(t);
    elseif t > 20 && t <= 65
        ET.P(t) = Kc.dev (t-20) * ET.o(t);
    elseif t > 65 && t <= 193
        ET.P(t) = Kc.mid     * ET.o(t);
    else
        ET.P(t) = Kc.end     * ET.o(t);
    end
end


nw         = ET.W / (phi * Zr);                   % ?w
n          = ET.P / (phi * Zr);                   % ?
m          = Ks  / (phi * Zr *(exp(beta*(1-s.fc))-1)) ;

%% 1 : micro irrigation 

for t=1:243 
if t == 1 
    S1(1)=s.star;
end
S1 (t) = min ( 1 ,(( h( t ) /(phi * Zr)) + S1( t )) );
if S1 ( t ) < s.star 
    I1 ( t ) =( s.star - S1(t))* (phi * Zr) ;
    S1 ( t ) = s.star ;
end    
flag  = 0;
k     = 0;
while flag == 0 
    k   = k+1 ;
    if k == 1
    s.k1 ( k , t+1 ) = S1(t);
    end
    s.b =( S1(t) +  s.k1( k , t+1 )) /2 ;
    if s.b > 0      && s.b <= s.h 
         p                =  0 ;
         s.k1 ( k+1 , t+1 )= S1(t) - p ;
         ET.a1(t)          = p * phi * Zr ;
    end
    if s.b > s.h    && s.b <= s.w 
        p                 = nw * ((s.b - s.h) / ( s.w - s.h )) ;
        s.k1 ( k+1 , t+1 ) = S1(t) - p ;
        ET.a1(t)           = p * phi * Zr ;
    end
    if s.b > s.w    && s.b <= s.star 
         p                = nw + ( n(t) - nw ) * (( s.b - s.w )/( s.star - s.w )) ;
         s.k1 ( k+1 , t+1 )= S1(t) - p ;
         ET.a1(t)          = p * phi * Zr ;
    end 
    if s.b > s.star && s.b <= s.fc 
         p                = n(t);
         s.k1 ( k+1 , t+1 )= S1(t) - p ;
         ET.a1(t)          = p * phi * Zr ;
    end     
    if s.b  > s.fc  &&  s.b <= 1 
         p                = n(t) + m * ( exp( beta * ( s.b - s.fc ) ) -1) ;
         s.k1 ( k+1 , t+1 )= S1(t) - p ;
         ET.a1(t)          = n(t) * phi * Zr ;
    end
   

    if abs ( s.k1( k+1 ,t+1 ) - s.k1( k , t+1 ) ) < 0.000001 
    flag   = 1;
    end 
end
    S1(t+1) = s.k1( k+1 , t+1 );
end



%% 2 : micro Deficit irrigation 

for t=1:243 
if t == 1 
    S2(1)=s.star;
end
S2 (t) = min ( 1 ,(( h( t ) /(phi * Zr)) + S2( t )) );
if S2 ( t ) < 0.7 * s.star 
    I2 ( t ) =( 0.7 * s.star - S2(t))* (phi * Zr) ;
    S2 ( t ) = 0.7 * s.star ;
end    
flag  = 0;
k     = 0;
while flag == 0 
    k   = k+1 ;
    if k == 1
    s.k2 ( k , t+1 ) = S2(t);
    end
    s.b =( S2(t) +  s.k2( k , t+1 )) /2 ;
    if s.b > 0      && s.b <= s.h 
         p                =  0 ;
         s.k2 ( k+1 , t+1 )= S2(t) - p ;
         ET.a2(t)          = p * phi * Zr ;
    end
    if s.b > s.h    && s.b <= s.w 
        p                 = nw * ((s.b - s.h) / ( s.w - s.h )) ;
        s.k2 ( k+1 , t+1 ) = S2(t) - p ;
        ET.a2(t)           = p * phi * Zr ;
    end
    if s.b > s.w    && s.b <= s.star 
         p                = nw + ( n(t) - nw ) * (( s.b - s.w )/( s.star - s.w )) ;
         s.k2 ( k+1 , t+1 )= S2(t) - p ;
         ET.a2(t)          = p * phi * Zr ;
    end 
    if s.b > s.star && s.b <= s.fc 
         p                = n(t);
         s.k2 ( k+1 , t+1 )= S2(t) - p ;
         ET.a2(t)          = p * phi * Zr ;
    end     
    if s.b  > s.fc  &&  s.b <= 1 
         p                = n(t) + m * ( exp( beta * ( s.b - s.fc ) ) -1) ;
         s.k2 ( k+1 , t+1 )= S2(t) - p ;
         ET.a2(t)          = n(t) * phi * Zr ;
    end
   

    if abs ( s.k2( k+1 ,t+1 ) - s.k2( k , t+1 ) ) < 0.000001 
    flag   = 1;
    end 
end
    S2(t+1) = s.k2( k+1 , t+1 );
end

%% 3 : traditional  irrigation 

for t=1:243 
if t == 1 
    S3(1)=s.star;                              
end
S3 (t) = min ( 1 ,(( h( t ) /(phi * Zr)) + S3( t )) );
if S3 ( t ) < s.star
   
    I3 ( t ) =( s.target - S3(t))* (phi * Zr) ;
    S3 ( t ) = s.target ;
end    
flag  = 0;
k     = 0;
while flag == 0 
    k   = k+1 ;
    if k == 1
    s.k3 ( k , t+1 ) = S3(t);
    end
    s.b =( S3(t) +  s.k3( k , t+1 )) /2 ;
    if s.b > 0      && s.b <= s.h 
         p                =  0 ;
         s.k3 ( k+1 , t+1 )= S3(t) - p ;
         ET.a3(t)          = p * phi * Zr ;
    end
    if s.b > s.h    && s.b <= s.w 
        p                 = nw * ((s.b - s.h) / ( s.w - s.h )) ;
        s.k3 ( k+1 , t+1 ) = S3(t) - p ;
        ET.a3(t)           = p * phi * Zr ;
    end
    if s.b > s.w    && s.b <= s.star 
         p                = nw + ( n(t) - nw ) * (( s.b - s.w )/( s.star - s.w )) ;
         s.k3 ( k+1 , t+1 )= S3(t) - p ;
         ET.a3(t)          = p * phi * Zr ;
    end 
    if s.b > s.star && s.b <= s.fc 
         p                = n(t);
         s.k3 ( k+1 , t+1 )= S3(t) - p ;
         ET.a3(t)          = p * phi * Zr ;
    end     
    if s.b  > s.fc  &&  s.b <= 1 
         p                = n(t) + m * ( exp( beta * ( s.b - s.fc ) ) -1) ;
         s.k3 ( k+1 , t+1 )= S3(t) - p ;
         ET.a3(t)          = n(t) * phi * Zr ;
    end
   

    if abs ( s.k3( k+1 ,t+1 ) - s.k3( k , t+1 ) ) < 0.000001 
    flag   = 1;
    end 
end
    S3(t+1) = s.k3( k+1 , t+1 );
end
%% 4 : traditional Deficit irrigation 

for t=1:243 
if t == 1 
    S4(1)=s.star;
end
S4 (t) = min ( 1 ,(( h( t ) /(phi * Zr)) + S4( t )) );
if S4 ( t ) < 0.7 * s.star 
    I4 ( t ) =( s.star - S4(t))* (phi * Zr) ;
    S4 ( t ) = s.star ;
end    
flag  = 0;
k     = 0;
while flag == 0 
    k   = k+1 ;
    if k == 1
    s.k4 ( k , t+1 ) = S4(t);
    end
    s.b =( S4(t) +  s.k4( k , t+1 )) /2 ;
    if s.b > 0      && s.b <= s.h 
         p                =  0 ;
         s.k4 ( k+1 , t+1 )= S4(t) - p ;
         ET.a4(t)          = p * phi * Zr ;
    end
    if s.b > s.h    && s.b <= s.w 
        p                 = nw * ((s.b - s.h) / ( s.w - s.h )) ;
        s.k4 ( k+1 , t+1 ) = S4(t) - p ;
        ET.a4(t)           = p * phi * Zr ;
    end
    if s.b > s.w    && s.b <= s.star 
         p                = nw + ( n(t) - nw ) * (( s.b - s.w )/( s.star - s.w )) ;
         s.k4 ( k+1 , t+1 )= S4(t) - p ;
         ET.a4(t)          = p * phi * Zr ;
    end 
    if s.b > s.star && s.b <= s.fc 
         p                = n(t);
         s.k4 ( k+1 , t+1 )= S4(t) - p ;
         ET.a4(t)          = p * phi * Zr ;
    end     
    if s.b  > s.fc  &&  s.b <= 1 
         p                = n(t) + m * ( exp( beta * ( s.b - s.fc ) ) -1) ;
         s.k4 ( k+1 , t+1 )= S4(t) - p ;
         ET.a4(t)          = n(t) * phi * Zr ;
    end
   

    if abs ( s.k4( k+1 ,t+1 ) - s.k4( k , t+1 ) ) < 0.000001 
    flag   = 1;
    end 
end
    S4(t+1) = s.k4( k+1 , t+1 );
end

%% NO irrigation
for t=1:243 
if t == 1 
    S5(1)=s.star;
end
S5 (t) = min ( 1 ,(( h( t ) /(phi * Zr)) + S5( t )) );

    
flag  = 0;
k     = 0;
while flag == 0 
    k   = k+1 ;
    if k == 1
    s.k5 ( k , t+1 ) = S5(t);
    end
    s.b =( S5(t) +  s.k5( k , t+1 )) /2 ;
    if s.b > 0      && s.b <= s.h 
         p                =  0 ;
         s.k5 ( k+1 , t+1 )= S5(t) - p ;
         ET.a5(t)          = p * phi * Zr ;
    end
    if s.b > s.h    && s.b <= s.w 
        p                 = nw * ((s.b - s.h) / ( s.w - s.h )) ;
        s.k5 ( k+1 , t+1 ) = S5(t) - p ;
        ET.a5(t)           = p * phi * Zr ;
    end
    if s.b > s.w    && s.b <= s.star 
         p                = nw + ( n(t) - nw ) * (( s.b - s.w )/( s.star - s.w )) ;
         s.k5 ( k+1 , t+1 )= S5(t) - p ;
         ET.a5(t)          = p * phi * Zr ;
    end 
    if s.b > s.star && s.b <= s.fc 
         p                = n(t);
         s.k5 ( k+1 , t+1 )= S5(t) - p ;
         ET.a5(t)          = p * phi * Zr ;
    end     
    if s.b  > s.fc  &&  s.b <= 1 
         p                = n(t) + m * ( exp( beta * ( s.b - s.fc ) ) -1) ;
         s.k5 ( k+1 , t+1 )= S5(t) - p ;
         ET.a5(t)          = n(t) * phi * Zr ;
    end
   

    if abs ( s.k5( k+1 ,t+1 ) - s.k5( k , t+1 ) ) < 0.000001 
    flag   = 1;
    end 
end
    S5(t+1) = s.k5( k+1 , t+1 );
end

%% 

alfa1 = 1 - Ky*(1-(sum(ET.a1)/sum(ET.P))) ;
alfa2 = 1 - Ky*(1-(sum(ET.a2)/sum(ET.P))) ;
alfa3 = 1 - Ky*(1-(sum(ET.a3)/sum(ET.P))) ;
alfa4 = 1 - Ky*(1-(sum(ET.a4)/sum(ET.P))) ;
alfa5 = 1 - Ky*(1-(sum(ET.a5)/sum(ET.P))) ;

% yield ( ton ) 

yield.i1 = alfa1 * yield.p ;
yield.i2 = alfa2 * yield.p ; 
yield.i3 = alfa3 * yield.p ; 
yield.i4 = alfa4 * yield.p ; 
yield.i5 = alfa5 * yield.p ;

% total irrigation ( 1000*m^3 ) 

I.i1 = sum( I1 ) * area * 10^-6 ;
I.i2 = sum( I2 ) * area * 10^-6 ;
I.i3 = sum( I3 ) * area * 10^-6 ;
I.i4 = sum( I4 ) * area * 10^-6 ;

t = 1 : 244 ;

figure

plot ( t , S1 , t , S3 , t , S5  )

figure

plot ( t , S1 , t , S2 , t , S5  )

figure

plot ( t , S3 , t , S4 , t , S5  )



clear alfa1 alfa2 alfa3 alfa4 alfa5 beta delta area    
clear flag k t Ks Kc Ky m n nw p phi s Zr I1 I2 I3 I4



