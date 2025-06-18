import matplotlib.pyplot as plt, numpy as np, arcpy as ac, shutil as sh
import os
from scipy.constants import golden_ratio as gold

ac.env.overwriteOutput = True
    # Globe parameters:
m = 300000000                                # scale of inscribed circle of output model 
Du = Dv = np.radians(10)                    # paralel and meridian intervals (10°)  
du = dv = np.radians(1)                     # graticule sampling resolution (1°)
radius = 6378137                            # reference sphere radius
u_gamma, u_sigma = 52.6226, 10.8123         # |latitudes| of face vertices 
o = 50                                      # graticule offset (relative) to face projection center
u_alpha = 26.5651      

def grat_lim(K: float, offset: float, plus: bool = True):
    '''Compute graticule limits [rad].'''
    if plus:                                                            # set 1st and 3rd operator to +, otherwise to -
        return np.radians(K + offset - np.mod(K + offset, 10))                        
    return np.radians(K - offset - np.mod(K - offset, 10))

def globe_face(ub: list[float], vb: list[float], R: float, K: np.array, 
              umin: float, umax: float, vmin: np.float64, vmax: np.float64, 
              Du: float, Dv: float, du: float, dv: float, 
              proj: callable, cw_deg: float, face_num: int, p: ac.mp.ArcGISProject,
              M: int, ref_lyt_cent: np.array):
    '''Create a face, plot it to .SVG with graticule and continents, import the face and graticule into an ArcGIS map.
    - cw_deg = CW rotation [deg] of a face'''          
    print(f'\nFace {face_num}:') 
    K, ub, vb = np.radians(K), np.radians(ub), np.radians(vb)                                   # Convert projection centers and face verticesto radians
    Xp, Yp, Xm, Ym = graticule(R, K[0], K[1], umin, umax, vmin, vmax, Du, Dv, du, dv, proj)     # Project graticule to gnomonic projection in oblique aspect     

    conts = []                                                                                  # create continent data
    for name in [r"DataSmall\Australia.txt", r"DataSmall\Antarctica.txt", r"DataSmall\America.txt", r"DataSmall\EurasiaAfrica.txt"]:            
        polygons = load(name)                                                                   # load a continent [rad]       
        xc, yc = continent(polygons, R, K[0], K[1], proj)                                       # Project to gnomonic projection in oblique aspect
        conts.append((xc, yc))
    print('Continents plotted onto face.')

    XC = np.concatenate([np.append(xc, np.nan) for xc, _ in conts])                             # join continents into per coordinate variables
    YC = np.concatenate([np.append(yc, np.nan) for _, yc in conts])                                          
    XB, YB = boundary(ub, vb, R, K[0], K[1], proj)                                              # Project face boundary
    positioned = rotate_export_plot(Xm, Ym, Xp, Yp, XC, YC, XB, YB, cw_deg, face_num)           # rotate the face, plot into .SVG    
    to_arcpy(positioned, face_num, R, K, p, cw_deg, M, ref_lyt_cent)                            # import data into ArcGIS Pro geodatabase and face´s predefined map  

def graticule(R, uk, vk, umin, umax, vmin, vmax, Du, Dv, du, dv, proj):
    '''Project Graticule to gnomonic projection in oblique aspect'''
    counter = 0                                         
    for u in np.arange(umin, umax + Du, Du):            # Create paralels
        vp = np.arange(vmin, vmax + dv, dv)             # Longitude of paralels´ points       
        up = [u] * len(vp)       
        
        sp, dp = uv_to_sd(up, vp, uk, vk)               # Convert to oblique aspect                
        xp, yp = proj(R, sp, dp)                        # Project to gnomonic projection 

        if counter == 0:                                # Add parallel to the list of parallels
            XP, YP = np.array([xp]), np.array([yp])
        else:
            XP, YP = np.append(XP, np.array([xp]), axis = 0), np.append(YP, np.array([yp]), axis = 0)
        counter += 1

    counter = 0
    for v in np.arange(vmin, vmax + Dv, Dv):            # Create meridians
        um = np.arange(umin, umax + du, du)             # Longitude of meridians´ points
        vm = [v] * len(um)  
        
        sm, dm = uv_to_sd(um, vm, uk, vk)               # Convert to oblique aspect
        xm, ym = proj(R, sm, dm)                        # Project to gnomonic projection 

        if counter == 0:                                # Add parallel to the list of parallels
            XM, YM = np.array([xm]), np.array([ym])
        else:
            XM, YM = np.append(XM, np.array([xm]), axis = 0), np.append(YM, np.array([ym]), axis = 0)
        counter += 1
    
    print('Graticule projected to gnomonic projection.')  
    return XP, YP, XM, YM

def uv_to_sd(ur, vr, ukr, vkr):      
    '''Convert WGS-84 coordinates to carthographic coordinates\n- input [rad], output [rad]'''
    dvr = vkr - vr
    sr = np.arcsin(np.sin(ur) * np.sin(ukr) + (np.cos(ur) * np.cos(ukr)) * np.cos(dvr))         # Transform latitude
    
    nom = np.sin(dvr) * np.cos(ur)                                                              # Transform longitude (quadrant adjusment)
    denom = ((np.cos(ur) * np.sin(ukr)) * np.cos(dvr)) - np.sin(ur) * np.cos(ukr)
    dr = np.arctan2(nom, denom)

    return sr, dr

def gnom(R, sr, dr):     
    '''Project to gnomonic projection\n - input [rad], output [rad]'''
    x = R * np.tan((np.pi / 2) - sr) * np.cos(dr)      
    y = -R * np.tan((np.pi / 2) - sr) * np.sin(dr)      
    return x, y

def load(filename):
    '''load polygon data into np.array [rad]'''
    with open(filename, encoding = 'utf-8') as f:
        polygons = np.array([line_to_pnt(f.readline())])
        l = f.readline()

        while l: 
            polygons = np.append(polygons, np.array([line_to_pnt(l)]), axis = 0)
            l = f.readline()
    
    return np.radians(polygons) 

def line_to_pnt(line):
    '''Convert <lat><lon> row to np.array'''
    s = line.split()
    return np.array([float(s[0]), float(s[1])])

def continent(C, R, uk, vk, proj):
    '''Project a continent to gnomonic projection in oblique aspect'''
    u, v = C[:, 0], C[:, 1]                     # Extract coordinates from input file as "<longitude> <latitude>"
    s, d = uv_to_sd(u, v, uk, vk)               # Transform to oblique aspect
    
        # Find points: s < s_min:
    s_min = np.radians(20);                     # minimal latitude for drawing [rad]
    indices = np.where(s >= s_min)              # Search indices with s < s_min
    s, d = s[indices], d[indices]               # only keep points with values >= than s_min

    XC, YC = proj(R, s, d)                      # Project to gnomonic projection
    return XC, YC

def boundary(u, v, R, uk, vk, proj):
    '''project face boundary to gnomonic projection in oblique aspect'''
    s, d = uv_to_sd(u, v, uk, vk)               # Transform to oblique aspect
    XB, YB = proj(R, s, d)                      # Transform to gnomonic projection
    return XB, YB

def rotate_export_plot(Xm, Ym, Xp, Yp, XC, YC, XB, YB, cw_deg, face_num):
    '''Rotate the face, plot into .SVG.
    - Return properly rotated faces.''' 
    for pair in ((Xm, Ym), (Xp, Yp)):           # rotate graticule
        shp = pair[0].shape    
        for i in np.arange(shp[0]):             
            for j in np.arange(shp[1]):         # rotate each point (set to 90° when calling globe_face)
                pair[0][i, j], pair[1][i, j] = rotate(np.array([pair[0][i, j], pair[1][i, j]]).transpose())        

    for pair in ((XC, YC), (XB, YB)):           # rotate continents and boundary
        shp = pair[0].shape    
        for i in np.arange(shp[0]): 
            pair[0][i], pair[1][i] = rotate(np.array([pair[0][i], pair[1][i]]).transpose())

    plt.plot(Xm.transpose(), Ym.transpose(), 'k')           # plot meridians                
    plt.plot(Xp.transpose(), Yp.transpose(), 'k')           # plot parallels 
    plt.plot(XC.transpose(), YC.transpose(), 'b')           # plot continents
    plt.plot(XB.transpose(), YB.transpose(), 'g')           # plot face boundary

    range = max(max(abs(XB)), max(abs(YB))) + 0.1           # plot range 0.1 rad bigger than face range
    plt.xlim([-range, range])                               # set the plot range
    plt.ylim([-range, range])

    plt.xticks()                                            # enable plot ticks
    plt.yticks()

    ax = plt.gca()
    ax.set_aspect('equal', adjustable = 'box')              # set equal axis scale

    plt.savefig(fr'faces\face{face_num}.svg')               # export face to .SVG
    plt.clf()                                               # reset plot figure
    print('Face exported to .SVG.')
    return Xm, Ym, Xp, Yp, XC, YC, XB, YB

def rotate(xy_column: np.array) -> np.array:
    '''Rotate a point by 270 degrees (CCW).'''
    phi = np.radians(270)                                               # CCW rotation angle
    matrix = np.round(np.array([[np.cos(phi), -np.sin(phi)],            # rotation matrix            
                                [np.sin(phi), np.cos(phi)]]),              
                      10)
    return matrix.dot(xy_column)

def to_arcpy(positioned, face_num, R, K, p, cw_deg, M, ref_lyt_cent):          
    '''import data into ArcGIS Pro geodatabase and face´s predefined map'''       
        # create face´s spatial reference:
    K = np.degrees(K)
    if face_num in (1, 2, 3, 4, 5, 11):                                                                # gnomonic polar north in oblique aspect, parametrized radius and cartographic pole
        sr = ac.SpatialReference(text = fr'PROJCS["face{face_num}_N", GEOGCS["GCS_WGS_1984", DATUM["D_WGS_1984", SPHEROID["WGS_1984",{R},298.257223563]], PRIMEM["Greenwich",0.0], UNIT["Degree",0.0174532925199433]],\
            PROJECTION["Gnomonic"], PARAMETER["False_Easting",0.0], PARAMETER["False_Northing",0.0], PARAMETER["Longitude_Of_Center",{K[1]}], PARAMETER["Latitude_Of_Center",{K[0]}], UNIT["Meter",1.0]]')        
    else:                                                                                            # gnomonic polar south in oblique aspect     
        sr = ac.SpatialReference(text = fr'PROJCS["face{face_num}_S", GEOGCS["GCS_WGS_1984", DATUM["D_WGS_1984", SPHEROID["WGS_1984",{R},298.257223563]], PRIMEM["Greenwich",0.0], UNIT["Degree",0.0174532925199433]],\
            PROJECTION["Gnomonic"], PARAMETER["False_Easting",0.0], PARAMETER["False_Northing",0.0], PARAMETER["Longitude_Of_Center",{K[1]}],PARAMETER["Latitude_Of_Center",{K[0]}], UNIT["Meter",1.0]]')        
        
    group_n = f'face{face_num}'                                                                         # group name: name of face´s map
    m = p.listMaps(group_n)[0]                                                                          # face´s map
    m.spatialReference = sr                                                                             # set map coordinate system
    clipper = rf'outputProj\globeFaces.gdb\boundary_{face_num}'                                         # face boundary geodatabase location
        # create pentagon boundary: 
    bnd_pts = ac.Array()                                                                                
    for i in range(len(positioned[6])):                                                                 # positioned[6] = XB, [7] = YB
        bnd_pts.add(ac.Point(positioned[6][i], positioned[7][i]))
    bnd_pts.add(bnd_pts[0])
    ac.management.CopyFeatures(ac.Polygon(bnd_pts, sr), clipper)                                        # copy geometry into geodatabase feature class
        # create arcpy.Geometry -> feature class from graticule:
    names = ('d','s')                                                                                   # graticule symbols          
    pairs = ((positioned[0], positioned[1]),                                                            # (Xm, Ym)
             (positioned[2], positioned[3]))                                                            # (Xp, Yp)
    for i in range(len(pairs)):                                                                                          
        x_arr, y_arr = pairs[i]
        lines = []

        for j in range(x_arr.shape[0]):
            arr = ac.Array()
            for k in range(x_arr.shape[1]):
                arr.add(ac.Point(x_arr[j, k], y_arr[j, k]))
            lines.append(ac.Polyline(arr, sr))

        ac.Clip_analysis(lines, clipper, rf'outputProj\globeFaces.gdb\{names[i]}_{face_num}')           # clip graticule to the pentagon boundary

    names = ('boundary', 'd', 's')                                                                      # all created layers
    sym_names = ('boundary', 'graticule', 'graticule')                                                  # symbology names          
        # apply symbology to the face:
    for i in range(3):                                                                                  
        layer_n = f'{names[i]}_{face_num}'                                                              # layer name
        layer_loc = os.path.join(os.getcwd(), rf'outputProj\globeFaces.gdb\{layer_n}')                  # layer location
        m.addDataFromPath(layer_loc)                                                                    # add graticule and boundary to current face                       
        sym_m = p.listMaps('symbology')[0]                                                              # find symbology map
        ac.management.ApplySymbologyFromLayer(m.listLayers(layer_n)[0], sym_m.listLayers(sym_names[i])[0])          # apply the symbology
        #set proper layer order, clip map with face boundary:
    bound = m.listLayers('boundary*')[0]                                                                # face boundary layer                                                                 
    m.moveLayer(m.listLayers()[0], bound)                                                               # make it top layer
    print('Data imported and visualized.')
    m.clipLayers(bound)                                                                                 # clip map to face boundary

    l = p.listLayouts()[0]                                                                              # map layout
    mf = l.listElements('mapframe_element', group_n)[0]                                                 # face´s map in layout
    bound_ext = ac.Describe(bound).extent
        # calculate globe model parameters:
    rho = R / M                                                                                         # radius of inscribed circle on the face
    print(f'rho: {round(rho * 1000, 2)} mm')
    a = (rho / (np.sqrt(25 + 10 * np.sqrt(5)) / 10)) * 1000                                             # face edge length
    print(f'face edge length: {round(a, 2)} mm')     
    face_scale = M / gold                                                                               # face_scale honoring globe scale where radius of its inscribed sphere r = R / M
    print(f'face scale = 1 : {int(face_scale)}')
        # set map frame properties: 
    mf.camera.setExtent(bound_ext)                                                                      # set map frame extent to face boundary extent                                                    
    mf.camera.heading = -cw_deg                                                                         # rotate map content (visible)          
    mf.elementRotation = -cw_deg                                                                        # rotate layout (invisible)
    mf.camera.scale = face_scale
    mf.elementWidth  = bound_ext.width * 1000 / mf.camera.scale                                         # fit map frame dimensions to face´s scale and extent
    mf.elementHeight = bound_ext.height * 1000 / mf.camera.scale
            # compute map frame center based on M = 100 000 000 with X, Y minimum face vertex coordinates both equal to 10 mm:
    ref_a = 92.67975559889058                                                                           # reference face edge length for M = 100 000 000
    lyt_cent = ref_lyt_cent * (a / ref_a)                                                               # adjust map frame center 
    mf.elementPositionX, mf.elementPositionY = lyt_cent[0], lyt_cent[1]

    print('Face scaled and placed into dodecahedral net.')

def fit_to_layout(p):
    '''Shift dodecahedral net to keep smallest x, y vertices 10 mm from layout edge. Choose smallest ISO portait (A5-A0) format to keep the net 10 or more mm from the edge.'''
    print()
    l = p.listLayouts()[0]                                                                  # map layout
    mfs = [l.listElements('mapframe_element', f'face{i}')[0] for i in range(1, 13)]         # face map frames in numerical order
    centers = np.array([[mf.elementPositionX, mf.elementPositionY] for mf in mfs])
        # set map frame anchors to min/max coordinates:
    mfs[3].setAnchor('BOTTOM_MID_POINT'); mfs[4].setAnchor('BOTTOM_MID_POINT'); mfs[6].setAnchor('TOP_MID_POINT'); mfs[7].setAnchor('TOP_MID_POINT')            
    x_min = mfs[4].elementPositionX; y_min = mfs[7].elementPositionY; x_max = mfs[6].elementPositionX; y_max = mfs[3].elementPositionY          # get face vertices min/max coordinates

    if x_min != 10 or y_min != 10:  	                                                    # are face vertices min X, Y both 10 mm?
        x_delta, y_delta = 10 - x_min, 10 - y_min
        for i in range(len(mfs)):                                              
            mf = mfs[i]
            mf.setAnchor('CENTER_POINT')                                                    # set anchor to center point again

            center = centers[i] + np.array([x_delta, y_delta])                              # compute shifted center
            mf.elementPositionX = center[0]; mf.elementPositionY = center[1]                # set adjusted center

        x_max += x_delta                                                                    # adjust face vertices max X, Y
        y_max += y_delta
        print('Globe model min X, Y set to 10 mm.')
   
    iso_portraits = {'A5': (148, 210), 'A4': (210, 297), 'A3': (297, 420), 'A2': (420, 594),  'A1': (594, 841), 'A0': (841, 1189)}
    width, height = 0, 0

    for dim in iso_portraits.values():                                                  # find smallest ISO portrait format to fit into with 10 mm distance
        if x_max <= (dim[0] - 10) and y_max <= (dim[1] - 10):
            width, height = dim[0], dim[1]                                              # ISO portrait found
            break
    if width == 0:                                                                      # ISO portrait not found
        width, height = x_max + 10, y_max + 10                                          # create custom layout

    l.changePageSize(width, height, False)                                              # change layout size
    print('Proper layout size set.')    

def export_map(p):
    '''export the map layout'''
    print('Exporting layout, may take long for high scale models.') 
    p.listLayouts()[0].exportToPDF('outputMap.pdf', 600, image_compression = 'JPEG2000', image_quality = 'NORMAL', embed_fonts = False) 
    p.save()
    print('Map exported, ArcGIS project saved.')
# ========================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================
if not os.path.isdir('faces'):                                  # check for 'faces' folder (for output .SVGs), create it if missing
    os.mkdir('faces')

try:                                                            # create output ArcGIS Pro project:
    sh.rmtree(r'outputProj')                                    # delete existing output project
except FileNotFoundError:                                       # ignore if non existent
    pass
except BaseException as e:
    raise SystemExit(f'Try exiting the folder or deleting it manually:\n {e}')

try:
    sh.copytree('referenceProj', 'outputProj')                  # make a copy of a reference project if present in CWD
except FileNotFoundError:
    raise SystemExit('Input ArcGIS project must be present in /referenceProj.')

p = ac.mp.ArcGISProject(r"outputProj\globeFaces.aprx")          # open reference project copy
print('Output ArcGIS Project initialized.')

    # Coordinates of face vertices [latitude; longitude]: 
A = np.array([-u_gamma, 0]);      B = np.array([-u_gamma, 72]);     C = np.array([-u_gamma, 144]);        D = np.array([-u_gamma, -144])    
E = np.array([-u_gamma, -72]);    F = np.array([-u_sigma, 0]);      G = np.array([u_sigma, 36]);          H = np.array([-u_sigma, 72])
I = np.array([u_sigma, 108]);     J = np.array([-u_sigma, 144]);    K = np.array([u_sigma, 180]);         L = np.array([-u_sigma, -144])
M = np.array([u_sigma, -108]);    N = np.array([-u_sigma, -72]);    O = np.array([u_sigma, -36]);         P = np.array([u_gamma, 36])
Q = np.array([u_gamma, 108]);     R = np.array([u_gamma, 180]);     S = np.array([u_gamma, -108]);        T = np.array([u_gamma, -36])      
    # per face parameters:                                     
        # centers of projection [latitude; longitude]:
centers = [np.array(i) for i in [[u_alpha, 0], [u_alpha, 72], [u_alpha, 144], [u_alpha, -144], [u_alpha, -72], [-u_alpha, 36], [-u_alpha, 108], [-u_alpha, 180], [-u_alpha, -108], [-u_alpha, -36], [90.0, 0], [-90.0, 0]]]       
        # latitudes of boundary pentagons:
ub = [[G[0], P[0], T[0], O[0], F[0], G[0]], [P[0], Q[0], I[0], H[0], G[0], P[0]], [K[0], J[0], I[0], Q[0], R[0], K[0]], [K[0], R[0], S[0], M[0], L[0], K[0]], [T[0], S[0], M[0], N[0], O[0], T[0]], [H[0], G[0], F[0], A[0], B[0], H[0]],
      [J[0], C[0], B[0], H[0], I[0], J[0]], [J[0], C[0], D[0], L[0], K[0], J[0]], [L[0], D[0], E[0], N[0], M[0], L[0]], [A[0], E[0], N[0], O[0], F[0], A[0]], [R[0], Q[0], P[0], T[0], S[0], R[0]], [C[0], B[0], A[0], E[0], D[0], C[0]]]
        # longitudes of boundary pentagons:
vb = [[G[1], P[1], T[1], O[1], F[1], G[1]], [P[1], Q[1], I[1], H[1], G[1], P[1]], [K[1], J[1], I[1], Q[1], R[1], K[1]], [K[1], R[1], S[1], M[1], L[1], K[1]], [T[1], S[1], M[1], N[1], O[1], T[1]], [H[1], G[1], F[1], A[1], B[1], H[1]],
      [J[1], C[1], B[1], H[1], I[1], J[1]], [J[1], C[1], D[1], L[1], K[1], J[1]], [L[1], D[1], E[1], N[1], M[1], L[1]], [A[1], E[1], N[1], O[1], F[1], A[1]], [R[1], Q[1], P[1], T[1], S[1], R[1]], [C[1], B[1], A[1], E[1], D[1], C[1]]] 
        # graticule limits [grad]:
umin = [grat_lim(i[0], o, False) for i in centers]                                                                      # latitude minimums                                                                                                                                                        
umax = [grat_lim(i[0], o) for i in centers]                                                                             # lat. max. 
vmin = np.concatenate(([grat_lim(i[1], o, False) for i in centers[:10]], [np.radians(-180)] * 2), axis = None)          # longitude min. | first 10 elements of centers, last 2 elements are -pi [rad]
vmax = np.concatenate(([grat_lim(i[1], o) for i in centers[:10]], [np.radians(180)] * 2), axis = None)                  # longitude max.
ccw_degs = [0, -72, -144, 144, 72, 0, 72, 144, -144, -72, 0, -36]                                                       # face rotation to fit into globe shell
    # layout face centers obtained by manually assembling a dodecahedral net for M = 100 000 000, X, Y minimum face vertex coordinates are both 10 mm:
ref_lyt_cents = np.array(((206.2879, 400.2061), (334.7686, 493.5598 ), (285.6747 , 644.5757), (126.8776 , 644.5886), (77.8196, 493.5536), (281.2674 , 312.0626 ), 
                        (409.6946 , 218.7081 ), (360.6628 , 67.6864), (201.8509, 67.6908 ), (152.7895 , 218.725), (206.2878, 542.8256), (276.833 , 183.068)))

for i in range(12):         # Create a face, plot it to .SVG, import into ArcGIS PRo         
    globe_face(ub[i], vb[i], radius, centers[i], umin[i], umax[i], vmin[i], vmax[i], Du, Dv, du, dv, gnom, ccw_degs[i], i + 1, p, m, ref_lyt_cents[i])           

fit_to_layout(p)
export_map(p)
del p