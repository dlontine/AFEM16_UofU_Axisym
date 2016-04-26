#----------------------------------------------------------------------------#
# ---------------------- Analytical Slolutions ------------------------------#
#----------------------------------------------------------------------------#
def A_Plate_Point_Fellipa_R(r,z,E,v,P,R,h,Ri):
    D = E*h**3/(12*(1-v**2))
    a = Ri
    b = R
    u_r = P * a**2 * (1+v) * (b**2 + r**2 * (1 - 2*v)) / (E * (b**2 - a**2) * r)
    return u_r
    
def A_Plate_Point_Fellipa_Z(r,z,E,v,P,R,h):
    D = E*h**3/(12*(1-v**2))
    u_z = -P/(16*math.pi*D) * ((3+v)/(1+v)*(R**2-r**2) + 2*r**2*math.log(r/R))
    return u_z
    

def PointLoadCenterDiscAnalyticUr(r,z,E,v,P,R,h,Ri):
    D = E*h**3/(12*(1-v**2))
    a = Ri
    b = R
    u_r = P * a**2 * (1+v) * (b**2 + r**2 * (1 - 2*v)) / (E * (b**2 - a**2) * r)
    return u_r
    
def A_Thick_Infinite_Cyl(E,v,P,OD,h,X=None,Y=None,inD=None,**kwargs):
    #Function robustness items:    
    if Xcoord is None:
        Xcoord=0.0
    if Ycoord is None:
        Ycoord=0.0 
    u_z=0.0 #All z displacement is fixed due to boundary conditions (plane strain)
    #From Fellipa "Verification Problems" eq 7.2:
    a=OD/2.0
    b=inD/2.0
    num=a**2*(1+v)*(b**2+r**2*(1-2*v))
    den=E*(b**2-a**2)*r
    u_r=P*num/den
    return u_r,u_z

#Roymech solutions:

def A_Plate_Point_Clamped(E,v,P,RO,h,z,r,RI):
    D = E*h**3/(12*(1-v**2))
    u_z = P*r**2/(16*math.pi*D)
    return u_z

def A_Plate_Point_Pinned(E,v,P,RO,h,z,r,RI):
    D = E*h**3/(12*(1-v**2))
    u_z = (5+v)*P*RO**4 / (64*(1+v)*D)
    return u_z

def A_Plate_Pressure_Clamped(E,v,P,RO,h,z,r,RI):
    D = E*h**3/(12*(1-v**2))
    u_z = P*RO**4/(64*D)
    return u_z

def A_Plate_Pressure_Pinned(E,v,P,RO,h,z,r,RI):
    D = E*h**3/(12*(1-v**2))
    u_z = (5+v)*P*r**4/(64*(1+v)*D)
    return u_z

def A_Washer_Point_Clamped(E,v,P,RO,h,z,r,RI):
    a = RO
    b = RI
    c = a/b
    t = h
    k = -.0016*c**6 + .0233*c**5 + -.1285*c**4 + .3072*c**3 - .2544*c**2 + .051
    u_z = k * P * a**2 / E * t**3
    return u_z

def A_Washer_Point_Pinned(E,v,P,RO,h,z,r,RI):
    a = RO
    b = RI
    c = a/b
    t = h
    k = 0.0111*c**6 - 0.1724*c**5 + 1.0195*c**4 - 2.7879*c**3 + 3.1547*c**2 -1.1484
    u_z = k * P * a**2 / E * t**3
    return u_z

def A_Washer_Pressure_Clamped(E,v,P,RO,h,z,r,RI):
    a = RO
    b = RI
    c = a/b
    t = h
    k = -0.0015*c**6 + 0.0230*c**5 + -0.1289*c**4 + .3166*c**3 + -0.2812*c**2 + 0.0733
    u_z = k * P * a**4 / (E * t**3)
    return u_z

def A_Washer_Pressure_Pinned(E,v,P,RO,h,z,r,RI):
    a = RO
    b = RI
    c = a/b
    t = h
    k = 0.01*c**6 + -.1585*c**5 + .9563*c**4 + -2.6988*c**3 + 3.2063*c**2 + -1.4443
    u_z = k * P * a**4 / (E * t**3)
    return u_z