import numpy as np
import IPython


""" -------------------- Diffractometer Calculations Master Functions -------------------- """

def hklFromAngles_2plus2(E, delta, gamma, omega, alpha, U, B, mode=1):
    """
    Calculate the hkl vector for a given set of angles.
    mode 1,3,5: horizontal, fixed incident angle
    mode 2,4,6: vertical, fixed incident angle
    """
    
    wavelength = 12.3984/E
    K = 2*np.pi/wavelength
    
    delta = np.deg2rad(delta)
    gamma = np.deg2rad(gamma)
    omega = np.deg2rad(omega)
    alpha = np.deg2rad(alpha)
    
    """rotation matrices"""
    if mode==1 or mode==3 or mode==5: # sample surface horizontal
        gamma = -gamma # sign convention for modes 1,3 and 5
        Delta = np.array([[1, 0, 0],
                        [0, np.cos(delta), -np.sin(delta)],
                        [0, np.sin(delta), np.cos(delta)]])

        Gamma = np.array([[np.cos(gamma), -np.sin(gamma), 0],
                        [np.sin(gamma), np.cos(gamma), 0],
                        [0, 0, 1]])
        
        Omega = np.array([[np.cos(omega), -np.sin(omega), 0],
                        [np.sin(omega), np.cos(omega), 0],
                        [0, 0, 1]])
        
        Alpha = np.array([[1, 0, 0],
                        [0, np.cos(alpha), -np.sin(alpha)],
                        [0, np.sin(alpha), np.cos(alpha)]])


    elif mode==2 or mode==4 or mode==6: # sample surface vertical
        Delta = np.array([[np.cos(delta), np.sin(delta), 0],
                        [-np.sin(delta), np.cos(delta), 0],
                        [0, 0, 1]])

        Gamma = np.array([[1, 0, 0],
                        [0, np.cos(gamma), -np.sin(gamma)],
                        [0, np.sin(gamma), np.cos(gamma)]])

        Omega = np.array([[np.cos(omega), np.sin(omega), 0],
                        [-np.sin(omega), np.cos(omega), 0],
                        [0, 0, 1]])

        Alpha = np.array([[1, 0, 0],
                        [0, np.cos(alpha), -np.sin(alpha)],
                        [0, np.sin(alpha), np.cos(alpha)]])

    
    """ calculate H """
    UBH = np.dot(Gamma, Delta) - np.identity(3)
    UBH = np.dot(UBH, np.array([0,K,0]))
    UBH = np.dot(np.linalg.inv(Alpha), UBH)
    UBH = np.dot(np.linalg.inv(Omega), UBH)

    Uinv = np.linalg.inv(U)
    Binv = np.linalg.inv(B)
    
    H = np.dot(Uinv,UBH)
    H = np.dot(Binv, H)
    
    Q = np.linalg.norm(UBH)
    tt = np.rad2deg(np.arccos( np.cos(-gamma) * np.cos(delta) ))
    
    return H, Q, tt




def hklToAngles_2plus2(hkl, E, alpha, U, B, mode=1):
    """
    Calculates the diffractometer angles for a given crystal (UB mat) and hkl
    mode 1: horizontal, fixed incident angle
    mode 2: vertical, fixed incident angle
    """
    
    wavelength = 12.3984/E
    K = 2*np.pi/wavelength
    
    UBH = U.dot(np.dot(B,hkl))
    Q = np.linalg.norm(UBH)

    """Test for the direct beam """
    Uinv = np.linalg.inv(U)
    Binv = np.linalg.inv(B)
    Y = np.dot(Uinv, np.array([0,1,0]))
    Y = np.dot(Binv, Y)
    
    print("The reciprocal lattice vector parallel to the x-ray beam at omega = 0 is: [%.4f %.4f %.4f]" % (Y[0], Y[1], Y[2]))
    
    if mode == 1:
        delta, gamma, omega, angout = diffAngles_mode1(UBH,K,alpha)
    elif mode == 2:
        delta, gamma, omega, angout = diffAngles_mode2(UBH,K,alpha)
        
    tt = np.arccos( np.cos(gamma) * np.cos(delta) )

    delta=np.rad2deg(delta)
    gamma=np.rad2deg(gamma)
    omega=np.rad2deg(omega)
    tt = np.rad2deg(tt)
    angout=np.rad2deg(angout)
        
    return delta, gamma, omega, tt, angout, Q




def hklToAngles_4C(hkl, N, E, U, B, mode, *args):
    wavelength = 12.3984/E

    if mode==1:
        alp = np.deg2rad(args[0])
        # print('\n/!\\ /!\\ /!\\ The input incident angle is assumed to be in degrees. /!\\ /!\\ /!\\  \n')

    UBH = U.dot(np.dot(B,hkl))
    UBN = U.dot(np.dot(B,N))
    Q = np.linalg.norm(UBH)
    
    """ 2theta """
    dspacing = 2*np.pi/np.linalg.norm(B.dot(hkl))
    ttheta = np.arcsin(wavelength / (2*dspacing))*2

    if mode==1:
        """ fixed incident angle """
        """ (I) find azimuthale angle for fixed alpha from Mochrie """ 
        """ (i) find N component parallel and perpendicular to hkl """
        Npar = np.dot( UBN/np.linalg.norm(UBN),UBH/np.linalg.norm(UBH) )
        Nperp = np.linalg.norm( UBN/np.linalg.norm(UBN) - Npar*UBH/np.linalg.norm(UBH) )

        """ (ii) calculate psi """
        temp = (-np.sin(alp)+Npar*np.sin(ttheta/2)) / (Nperp*np.cos(ttheta/2))
        """ catch unobtainable alphas """
        if temp>1:
            temp = 1
        elif temp<-1:
            temp = -1

        psi = np.arccos(temp)
        print('psi = '+str(np.rad2deg(psi)))


        """ (II) Calulation of the diffractometer angle for the given azimuth psi """
        """ (i) construct the matrix T """
        """ t1 parallel to hkl """
        t1 = -UBH/np.linalg.norm(UBH) # negative because of rotation convention
        """ t2 in the plane of hkl and N """
        t2 = np.cross(UBH, UBN)
        t2 = np.cross(t2,t1)
        t2 = -t2/np.linalg.norm(t2) # negative for same reason
        """ t3 vector normal to t1 and t2 """
        t3 = np.cross(t1,t2)
        T = np.transpose(np.array([t1,t2,t3]))
        Tinv = np.transpose(T) # T is orthogonal (from paper R_0=Tinv)

        """ (ii) calculate angles """
        Psi = rotation_matrix([1,0,0],psi)
        R = Psi@Tinv

        chi = np.arctan2( np.sqrt(R[2,0]**2+R[2,1]**2), R[2,2] )
        phi = np.arctan2( -R[2,1], R[2,0] )
        omega = np.arctan2( -R[1,2], -R[0,2] )
        theta = ttheta/2 + omega
        
        return ttheta, theta, chi, phi, omega, Q










""" -------------------- Tools and low level functions -------------------- """

def UBmat(a, aa, N):
    """Find the length and angles of the real and reciprocal vectors"""
    a0,a1,a2 = vectorFromLengthAndAngles(a,aa)
    a,aa,b,ba = lengthAndAngles(a0,a1,a2)
    
    """Orthonormalize the reciprocal lattice vectors (B matrix)"""
    B,N = ortho_recip(N,a,aa,b,ba)
    
    """Rotation angle and normal with respect to z axis (U matrix)"""
    z = np.array([0,0,1])
    angle = np.arccos( np.dot(z,N) /np.linalg.norm(z)/np.linalg.norm(N))  #radian
    rotN = np.cross(N,z)/np.linalg.norm(np.cross(N,z)) #axis of rotation
    U = rotation_matrix(rotN, angle)
    
    return U, B
    
    
    
    
def vectorFromLengthAndAngles(a,aa):
    """
    Takes a set of lattice vectors lengths and angles and calculates the cartesian 
    vectors.
    a   : Length of lattice vectors.
    aa  : Angles between lattice vectors [Between a3 and a2, a3 and a1, a1
        and a2].
    b   : Length of reciprocal vectors.
    ba  : Angles between reciprocal vectors.
    a1,a2,a3    : Lattice vectors given in cubic coordinates as [x y z].
    """
    aa = np.deg2rad(aa)

    a0 = a[0]*  np.array([1,0,0])
    a1 = a[1]* np.array([np.cos(aa[2]), np.sin(aa[2]), 0 ])
    v = np.sqrt(1-np.cos(aa[0])**2 - np.cos(aa[1])**2 - np.cos(aa[2])**2
        + 2*np.cos(aa[0])*np.cos(aa[1])*np.cos(aa[2]) )
    a2 = a[2]* np.array([np.cos(aa[1]), (np.cos(aa[0])-np.cos(aa[1])*np.cos(aa[2]))/np.sin(aa[2]),
        v/np.sin(aa[2])] )

    return a0, a1, a2

    


def lengthAndAngles(a0,a1,a2):
    """
    Takes a set of lattice vectors, and calculates the
    length of and angles between both the vectors and their corresponding
    reciprocal vectors.
    
    a   : Length of lattice vectors.
    aa  : Angles between lattice vectors [Between a3 and a3, a3 and a1, a1
        and a2].
    b   : Length of reciprocal vectors.
    ba  : Angles between reciprocal vectors.
    a1,a2,a3    : Lattice vectors given in cubic coordinates as [x y z].
    """

    a = np.array([np.linalg.norm(a0), np.linalg.norm(a1), np.linalg.norm(a2) ])
    aa = np.array([np.arccos(np.dot(a1,a2)/a[1]/a[2]), 
                   np.arccos(np.dot(a2,a0)/a[2]/a[0]), 
                   np.arccos(np.dot(a0,a1)/a[0]/a[1])])
    
    nor = np.dot(a0,np.cross(a1,a2))
    b0 = 2*np.pi*np.cross(a1,a2)/nor
    b1 = 2*np.pi*np.cross(a2,a0)/nor
    b2 = 2*np.pi*np.cross(a0,a1)/nor
    
    b = np.array([np.linalg.norm(b0), np.linalg.norm(b1), np.linalg.norm(b2) ])
    ba = np.array([np.arccos(np.dot(b1,b2)/b[1]/b[2]),
                   np.arccos(np.dot(b2,b0)/b[2]/b[0]),
                   np.arccos(np.dot(b0,b1)/b[0]/b[1])])
    
    aa = np.rad2deg(aa)
    ba = np.rad2deg(ba)
                   
    return a, aa, b, ba
    
    


def rotation_matrix(axis, theta):
    """
    Return the rotation matrix associated with counterclockwise rotation about
    the given axis by theta radians.
    """
    axis = np.asarray(axis)
    theta = np.asarray(theta)
    axis = axis/np.sqrt(np.dot(axis, axis))
    a = np.cos(theta/2.0)
    b, c, d = -axis*np.sin(theta/2.0)
    aa, bb, cc, dd = a*a, b*b, c*c, d*d
    bc, ad, ac, ab, bd, cd = b*c, a*d, a*c, a*b, b*d, c*d
    return np.array([
        [aa+bb-cc-dd, 2*(bc+ad), 2*(bd-ac)],
        [2*(bc-ad), aa+cc-bb-dd, 2*(cd+ab)],
        [2*(bd+ac), 2*(cd-ab), aa+dd-bb-cc]
        ])


                     

def ortho_recip(H,a,aa,b,ba):
    """
    Rewrite a reciprocal lattice vector in terms of a cartesian
    basic.
 
    [Hw,B] = cartesian(H,a,aa,b,ba);
    
    Hw       : [h;k;l], reciprocal letter in the orthogonal coordinate
            system.
    B       : The orthonomalization matrix.
    H      : [hw;kw;lw] reciprocal lattice vector in original coordiante
            system.
    a,aa,b,ba: Output of lengthAndAngles.m.
    """

    ba = np.deg2rad(ba)
    aa= np.deg2rad(aa)
    
    B = np.array([ [b[0], b[1]*np.cos(ba[2]), b[2]*np.cos(ba[1])],
                   [0, b[1]*np.sin(ba[2]), -b[2]*np.sin(ba[1])*np.cos(aa[0])],
                   [0, 0, 2*np.pi/a[2]] ])
    Hw = np.dot(B,H)
    
    return B, Hw




def diffAngles_mode1(H,K,alpha):
    """
    mode 1: horizontal geometry, fixed incident angle
    inputs:
        H: orthonormalized hkl
        K: wavevector
        alpha: incident angle
        
    outputs:
        delta, gamma, omega: diffractometer angles
        angOut: outgoing angle
        
        
    based on Simon's code
    """
    
    alpha=np.deg2rad(alpha)

    h = H[0]/K
    k = H[1]/K
    l = H[2]/K
    
    Y = -.5*(h**2+k**2+l**2)
    Z = (l+np.sin(alpha)*Y)/np.cos(alpha)
    X = np.sqrt( h**2+k**2 - (np.cos(alpha)*Y+np.sin(alpha)*Z)**2 ) # Sign is +/-, but + gives positive delta.
    W = np.cos(alpha)*Y+np.sin(alpha)*Z
    
    gamma = -np.arctan2(-X,Y+1)
    delta = np.arctan2(Z*np.sin(gamma),X)
    omega = np.arctan2(h*W-k*X,h*X+k*W)    
    angOut = np.arcsin(l-np.sin(alpha))
    
    return delta, gamma, omega, angOut



def diffAngles_mode2(H,K,alpha):
    """
    mode 2: vertical geometry, fixed incident angle
    inputs:
        H: orthonormalized hkl
        K: wavevector
        alpha: incident angle
        
    outputs:
        delta, gamma, omega: diffractometer angles
        angOut: outgoing angle
        
        
    based on Simon's code
    """
    
    alpha=np.deg2rad(alpha)

    h = H[0]/K
    k = H[1]/K
    l = H[2]/K

    Y = -.5*(h**2+k**2+l**2)
    Z = (l+np.sin(alpha)*Y)/np.cos(alpha)
    X = np.sqrt( h**2+k**2 - (np.cos(alpha)*Y+np.sin(alpha)*Z)**2 ) # Sign is +/-, but + gives positive delta.
    W = np.cos(alpha)*Y+np.sin(alpha)*Z

    omega = np.arctan2(k*X-h*W, h*X+k*W)
    delta = np.arcsin(X)
    gamma = np.arctan2(Z, Y+1)
    angOut = np.arcsin(l-np.sin(alpha))

    return delta, gamma, omega, angOut










"""------------ interface functions for the GUI ------------"""
def main_4C(hkl, a, aa, N, E, *args):

    alp = args[0]

    wavelength = 12.3984/E
    K = 2*np.pi/wavelength

    a0,a1,a2 = vectorFromLengthAndAngles(a,aa)
    a,aa,b,ba = lengthAndAngles(a0,a1,a2)
    U,B = UBmat(a, aa, N)
    return hklToAngles_4C(hkl, N, E, U, B, 1, alp)


def main_22C(hkl, a, aa, N, E, mode, *args):

    alp = args[0]

    wavelength = 12.3984/E
    K = 2*np.pi/wavelength

    a0,a1,a2 = vectorFromLengthAndAngles(a,aa)
    a,aa,b,ba = lengthAndAngles(a0,a1,a2)
    U,B = UBmat(a, aa, N)

    return hklToAngles_2plus2(hkl, E, alp, U, B, mode)


def main_22C_hkl(E, delta, gamma, omega, alpha, a, aa, N, mode, *args):
    wavelength = 12.3984/E
    K = 2*np.pi/wavelength

    a0,a1,a2 = vectorFromLengthAndAngles(a,aa)
    a,aa,b,ba = lengthAndAngles(a0,a1,a2)
    U,B = UBmat(a, aa, N)
    return hklFromAngles_2plus2(E, delta, gamma, omega, alpha, U, B, mode)