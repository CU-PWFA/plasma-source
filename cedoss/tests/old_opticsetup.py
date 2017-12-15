
if choice==3:#Most simple, but optic is close to gas jet
    f = 0.030 #m
    #  30mm is good, but expensive
    #  0.00625 0.0125 0.025 0.050 0.100
    #  0.015 0.020 0.030 0.040 0.075 0.150
    q_y = GB.Prop_Init_q(wavelength, w0, -.5, 1)
    q_z = GB.Prop_Init_q(wavelength, w0, -.5, 1)
    
    GB.Prop_Cylindrical_FreeSpace(q_y,q_z,.500,l_step)
    GB.Prop_CylindricalLens(q_z,q_y,1.000)
    GB.Prop_CylindricalLens(q_y,q_z,1.000)
    
    GB.Prop_Cylindrical_FreeSpace(q_y,q_z,.5-f,l_step) 
    GB.Prop_CylindricalLens(q_y,q_z,-2*f)
    
    GB.Prop_Cylindrical_FreeSpace(q_y,q_z,.4,l_step)
    q_x=q_y; q_y=q_z

if choice==2: #A bunch, but it works well
    q_y = GB.Prop_Init_q(wavelength, w0, -.5, 1)
    q_z = GB.Prop_Init_q(wavelength, w0, -.5, 1)
    
    GB.Prop_Cylindrical_FreeSpace(q_y,q_z,.500,l_step)
    GB.Prop_CylindricalLens(q_z,q_y,1.000)
    GB.Prop_CylindricalLens(q_y,q_z,1.000)
    
    GB.Prop_Cylindrical_FreeSpace(q_y,q_z,.150,l_step) 
    GB.Prop_CylindricalLens(q_y,q_z,.200)
    
    GB.Prop_Cylindrical_FreeSpace(q_y,q_z,.030,l_step) 
    GB.Prop_CylindricalLens(q_y,q_z,-.2)
    
    GB.Prop_Cylindrical_FreeSpace(q_y,q_z,.020,l_step) 
    GB.Prop_CylindricalLens(q_y,q_z,-.2)
    
    GB.Prop_Cylindrical_FreeSpace(q_y,q_z,.21,l_step) 
    GB.Prop_CylindricalLens(q_y,q_z,-.1)
    
    GB.Prop_Cylindrical_FreeSpace(q_y,q_z,.4,l_step)
    q_x=q_y; q_y=q_z

if choice==1:
    q_y = GB.Prop_Init_q(wavelength, w0, -.5, 1)
    q_z = GB.Prop_Init_q(wavelength, w0, -.5, 1)
    
    GB.Prop_Cylindrical_FreeSpace(q_y,q_z,.5,l_step)
    GB.Prop_CylindricalLens(q_z,q_y,1)
    GB.Prop_Cylindrical_FreeSpace(q_y,q_z,.1,l_step)
    
    GB.Prop_CylindricalLens(q_y,q_z,.3)
    GB.Prop_Cylindrical_FreeSpace(q_y,q_z,.159,l_step)

    GB.Prop_CylindricalLens(q_y,q_z,.018)
    GB.Prop_Cylindrical_FreeSpace(q_y,q_z,.4,l_step)
    q_x=q_y; q_y=q_z

if choice==0: #1 meter, 3 lens
    #Focuses wz over a distance of 1m with 2 cyl. lens to .35mm
    q_y = GB.Prop_Init_q(wavelength, w0, -.5, 1)
    q_z = GB.Prop_Init_q(wavelength, w0, -.5, 1)
    
    GB.Prop_CylindricalLens(q_y, q_z, .8)
    GB.Prop_Cylindrical_FreeSpace(q_y,q_z,.35,l_step)
    GB.Prop_CylindricalLens(q_y, q_z, -.105)
    GB.Prop_Cylindrical_FreeSpace(q_y,q_z,.15,l_step)
    
    GB.Prop_CylindricalLens(q_z, q_y, 1)
    GB.Prop_Cylindrical_FreeSpace(q_y,q_z,1,l_step)
    q_x=q_y; q_y=q_z