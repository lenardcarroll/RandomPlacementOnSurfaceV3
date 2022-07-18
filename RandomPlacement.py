import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import argparse
from shapely.geometry import Point
from shapely.geometry.polygon import Polygon
import random
from pathlib import Path
import math
import warnings
warnings.filterwarnings("ignore", category=RuntimeWarning) 


#Arguments Cheat/Help Sheet
parser = argparse.ArgumentParser()
parser.add_argument("-in", "--input", dest = "input", default = "input.xyz", help="Name of the input structure file, which contains the coordinates.")
parser.add_argument("-mdist", "--mindist", dest = "mindist", default = "2.15", help="Specify the minimum allowed distance between surface atoms and adsorbed atoms")
parser.add_argument("-xdist", "--maxdist", dest = "maxdist", default = "2.50", help="Specify the maximum allowed distance between surface atoms and adsorbed atoms")
parser.add_argument("-adist", "--adist", dest = "adist", default = "4.25", help="Specify the minimum allowed distance between adsorbed atoms")
parser.add_argument("-out", "--output", dest = "output", default = "output.xyz", help="Choose the file you want the final coordinates to be printed out to.")
parser.add_argument("-plot", "--plot", dest = "plot", default = "N", help="If Y is chosen, a simple plot is made of the surface and the random points. Otherwise (N) none is made.")
args = parser.parse_args()

VecAndSurface = open("Vectors.txt","r")
content = VecAndSurface.readlines()

CellVectors = []
for i in range(0,3):
    CellVectors.append(content[i].split())
for i in range(0, len(CellVectors)):
    CellVectors[i][0] = float(CellVectors[i][0])
    CellVectors[i][1] = float(CellVectors[i][1])
    CellVectors[i][2] = float(CellVectors[i][2])

SurfaceAtoms = content[3]
Adsorbates = content[4]

commapos = []
for i in range(len(SurfaceAtoms)):
    if SurfaceAtoms[i]==',':
        commapos.append(i)
SRange = []
if len(commapos)==0:
    x = SurfaceAtoms[:len(SurfaceAtoms)]
    if '-' in x:
        xindex = x.index('-')
        firstval = int(x[:xindex])
        secondval = int(x[xindex+1:])
        for i in range(firstval-1,secondval):
            SRange.append(i)
else:
    for i in range(len(commapos)):
        if i==0:
            x = SurfaceAtoms[:commapos[i]]
            if '-' in x:
                xindex = x.index('-')
                firstval = int(x[:xindex])
                secondval = int(x[xindex+1:])
                for i in range(firstval-1,secondval):
                    SRange.append(i)
            else:
                SRange.append(int(x)-1)
        else:
            x = SurfaceAtoms[commapos[i-1]+1:commapos[i]]
            if '-' in x:
                xindex = x.index('-')
                firstval = int(x[:xindex])
                secondval = int(x[xindex+1:])
                for i in range(firstval-1,secondval):
                    SRange.append(i)
            else:
                SRange.append(int(x)-1)
if len(commapos)>0:
    x =  SurfaceAtoms[commapos[len(commapos)-1]+1:]
    if '-' in x:
        xindex = x.index('-')
        firstval = int(x[:xindex])
        secondval = int(x[xindex+1:])
        for i in range(firstval-1,secondval):
            ORange.append(i)
    else:
        ORange.append(int(x)-1)

commapos = []
for i in range(len(Adsorbates)):
    if Adsorbates[i]==',':
        commapos.append(i)
ORange = []
if len(commapos)==0:
    x = Adsorbates[:len(Adsorbates)]
    if '-' in x:
        xindex = x.index('-')
        firstval = int(x[:xindex])
        secondval = int(x[xindex+1:])
        for i in range(firstval-1,secondval):
            ORange.append(i)
    else:
        if x!='None':
             ORange.append(int(x)-1)

else:
    for i in range(len(commapos)):
        if i==0:
            x = Adsorbates[:commapos[i]]
            if '-' in x:
                xindex = x.index('-')
                firstval = int(x[:xindex])
                secondval = int(x[xindex+1:])
                for i in range(firstval-1,secondval):
                    ORange.append(i)
            else:
                ORange.append(int(x)-1)
        else:
            x = Adsorbates[commapos[i-1]+1:commapos[i]]
            if '-' in x:
                xindex = x.index('-')
                firstval = int(x[:xindex])
                secondval = int(x[xindex+1:])
                for i in range(firstval-1,secondval):
                    ORange.append(i)
            else:
                ORange.append(int(x)-1)
if len(commapos)>0:
    x =  Adsorbates[commapos[len(commapos)-1]+1:]
    if '-' in x:
        xindex = x.index('-')
        firstval = int(x[:xindex])
        secondval = int(x[xindex+1:])
        for i in range(firstval-1,secondval):
            ORange.append(i)
    else:
        ORange.append(int(x)-1)


altCornerVecs = []
CornerVecs = [[0,0,0],[0,1,0],[1,1,0],[1,0,0],[0,0,1],[1,0,1],[0,1,1],[1,1,1]]

df1 = pd.read_csv(args.input,skiprows=2,names=["Atoms","X","Y","Z"],sep="\s+", engine='python')
SurfaceAtoms = []
OtherAtoms = []
for i in range(len(df1)):
    if i in SRange:
        SurfaceAtoms.append({"Atoms": df1["Atoms"].iloc[i], "X": df1["X"].iloc[i], "Y": df1["Y"].iloc[i], "Z": df1["Z"].iloc[i]})
    elif i in ORange:
        OtherAtoms.append({"Atoms": df1["Atoms"].iloc[i], "X": df1["X"].iloc[i], "Y": df1["Y"].iloc[i], "Z": df1["Z"].iloc[i]})

dfSAtoms = pd.DataFrame(SurfaceAtoms)
dfOAtoms = pd.DataFrame(OtherAtoms)

def getCoordinates(Vector,CellVectors):
    return np.matmul(Vector,CellVectors)
def newCoordinates():
    for i in CornerVecs:
        altCornerVecs.append(getCoordinates(i,CellVectors))
newCoordinates()
distVals = [[0,0,0]]
distCombos = [[0,1],[0,2],[0,3],[1,0],[1,3],[2,0],[2,1],[3,1]]
for i in distCombos:
    distVals.append([altCornerVecs[i[1]][0] - altCornerVecs[i[0]][0],altCornerVecs[i[1]][1] - altCornerVecs[i[0]][1],altCornerVecs[i[1]][2] - altCornerVecs[i[0]][2]])

how_many_mol = []
num_lines = []
for p in Path('.').glob('*.atom'):
    with open(p.name) as f:
        lines = f.read()
        first = lines.split('\n', 1)[0]
    how_many_mol.append(int(first))
    num_lines.append(len(list(open(p.name))))
    if len(how_many_mol) == 1:
        df = pd.read_csv(p.name, skiprows=1, names=['Atoms', 'X', 'Y', 'Z'], sep="\s+", engine='python')
    else:
        df2 = pd.read_csv(p.name, skiprows=1, names=['Atoms', 'X', 'Y', 'Z'], sep="\s+", engine='python')
        df = pd.concat([df, df2], axis=0)

#Split all the loaded in molecules (to be added onto the surface) into separate parts.
sum_lines = 0
range_info = []
for i in range(len(num_lines)):
    range_info.append([sum_lines,sum_lines+num_lines[i]-2])
    sum_lines += num_lines[i] - 1
frames = [ df.iloc[i[0]:i[1]+1].copy() for i in range_info ]

#Determines how many times a specific molecule should be placed on the surface
How_Many_Systems = []
for i in range(len(how_many_mol)):
    for j in range(how_many_mol[i]):
        for k in range(len(frames[i])):
            How_Many_Systems.append(i)

dfSAtoms = dfSAtoms.sort_values(by ='X', ascending = True)

coordinates = []
for i in dfSAtoms.index:
    coordinates.append([str(dfSAtoms['Atoms'].iloc[i]), dfSAtoms['X'].iloc[i], dfSAtoms['Y'].iloc[i], dfSAtoms['Z'].iloc[i]])
O_atoms = []
for i in dfOAtoms.index:
    O_atoms.append([str(dfOAtoms['Atoms'].iloc[i]), dfOAtoms['X'].iloc[i], dfOAtoms['Y'].iloc[i], dfOAtoms['Z'].iloc[i]])
O_atomslen = len(O_atoms)
#Added sorted X and Y coordinates to list points
points = []
for i in range(len(dfSAtoms)):
    points.append([dfSAtoms['X'].iloc[i],dfSAtoms['Y'].iloc[i]])

#Add first point into list convex_points
convex_points = [[points[0][0],points[0][1]]]

#Give the first point an angle 90
angles = [90]

#Algorithm starts, check https://github.com/lenardcarroll/myConvexHull.py/blob/main/README.md for an explanation
j = 0
while j>-1:
    altered_angles = []
    k_points = []
    unaltered_angles = []
    i = convex_points[len(convex_points)-1]
    for k in range(len(points)):
        if i!=points[k]:
            altx1y1 = [i[0] - i[0],i[1] - i[1]]
            altx2y2 = [points[k][0] - i[0],points[k][1] - i[1]]
            altx2y2 = altx2y2/np.linalg.norm(altx2y2)
            ang = np.arctan2(altx2y2[1],altx2y2[0]) - np.arctan2(altx1y1[1],altx1y1[0])
            altang = ang - angles[j] + 2*np.pi
            if altang > 2*np.pi:
                altang = altang - 2*np.pi
            altered_angles.append(altang)
            unaltered_angles.append(ang)
            k_points.append(k)
    max_altered_angle = max(altered_angles)
    convex_points_max = points[k_points[altered_angles.index(max_altered_angle)]]
    convex_points.append(convex_points_max)
    convex_points_unaltered_angle = unaltered_angles[altered_angles.index(max_altered_angle)]
    angles.append(convex_points_unaltered_angle)

    if convex_points_max == convex_points[0]:
        j=-1
    if j > 0:
        points.remove(convex_points_max)
        j+=1
    elif j == 0:
        j+=1

#Construct a polygon out of the convex points
polygonsurface = []
for i in convex_points:
    polygonsurface.append((i[0],i[1]))
polygon = Polygon(polygonsurface)

#Find the midpoint of the polygon
sumMidPointx = 0
sumMidPointy = 0
for i in convex_points:
    sumMidPointx += i[0]
    sumMidPointy += i[1]   
MidPoint = [sumMidPointx/len(convex_points),sumMidPointy/len(convex_points)]

#Split the polygon into a bunch of triangles using the polygon vertices
Triangles = []
for i in range(1,len(convex_points)):
    Triangles.append([convex_points[i-1],convex_points[i],MidPoint])

OO_distance = float(args.adist)

#function to calculate the distance between multiple vectors
def dist2(a,b,c):
    if np.abs(np.sqrt((a[0] - b[0])**2 + (a[1] - b[1])**2)-np.sqrt((a[0] - c[0])**2 + (a[1] - c[1])**2)-np.sqrt((c[0] - b[0])**2 + (c[1] - b[1])**2))<0.001:
        return True

multiSurface = []
for i in coordinates:
    for j in distVals:
        multiSurface.append([i[1] + j[0], i[2] + j[1], i[3] + j[2]])

i=0
while i<len(How_Many_Systems):
    TriangleChoice = random.choice(range(len(Triangles)))
    j = Triangles[TriangleChoice]
    x1y1 = [j[0][0]-j[2][0],j[0][1]-j[2][1]]
    x2y2 = [j[1][0]-j[2][0],j[1][1]-j[2][1]]
    MidPoint = [0,0]
    r = 0
    while r>-1:
        a1 = random.uniform(0,1)
        a2 = random.uniform(0,1)
        if (a1+a2)>1:
            r+=1
        else:
            r=-1
    newxy = [a1*x1y1[0]+a2*x2y2[0],a1*x1y1[1]+a2*x2y2[1]]
    newxy = [newxy[0]+j[2][0],newxy[1]+j[2][1]]
    if polygon.contains(Point(newxy[0],newxy[1])) == True:
        z_preliminary = random.uniform(float(args.mindist), float(args.maxdist))
        atom_distances = []
        Error_list = []
        for k in range(len(coordinates)):
            atom_distances.append(np.sqrt((coordinates[k][1]-newxy[0])**2+(coordinates[k][2]-newxy[1])**2))
        closest_member = min(atom_distances)
        if closest_member > z_preliminary:
            i+=0
            continue
        atom_number = 0
        for k in range(len(coordinates)):
            if closest_member == atom_distances[k]:
                atom_number = k
        z_value = np.sqrt((z_preliminary)**2 - (coordinates[atom_number][1]-newxy[0])**2 - (coordinates[atom_number][2]-newxy[1])**2) + coordinates[atom_number][3]
        if math.isnan(z_value) == True:
            i+=0
            continue
        newO_atoms = []
        tempO_atoms = []
        for j in range(len(frames[How_Many_Systems[i]])):
            tempO_atoms.append([frames[How_Many_Systems[i]].iloc[j]['Atoms'],frames[How_Many_Systems[i]].iloc[j]['X']+newxy[0],frames[How_Many_Systems[i]].iloc[j]['Y']+newxy[1],frames[How_Many_Systems[i]].iloc[j]['Z']+z_value])
        for k in O_atoms:
            for j in distVals:
                newO_atoms.append([k[1] + j[0], k[2] + j[1], k[3] + j[2]])    
        if len(newO_atoms)>0:
            for k in range(len(newO_atoms)):
                for l in range(len(tempO_atoms)):
                    value = (np.sqrt((newO_atoms[k][0]-tempO_atoms[l][1])**2+(newO_atoms[k][1]-tempO_atoms[l][2])**2+(newO_atoms[k][2]-tempO_atoms[l][3])**2))
                    if value<OO_distance:
                        Error_list.append("Error")
        if len(Error_list)>0:
            i+=0
            continue
        tempO_atoms = []
        for k in range(len(frames[How_Many_Systems[i]])):
            tempO_atoms.append([frames[How_Many_Systems[i]].iloc[k]['Atoms'],frames[How_Many_Systems[i]].iloc[k]['X']+newxy[0],frames[How_Many_Systems[i]].iloc[k]['Y']+newxy[1],frames[How_Many_Systems[i]].iloc[k]['Z']+z_value])
        for k in range(len(multiSurface)):
            for l in range(len(tempO_atoms)):
                value = np.sqrt((multiSurface[k][0]-tempO_atoms[l][1])**2+(multiSurface[k][1]-tempO_atoms[l][2])**2+(multiSurface[k][2]-tempO_atoms[l][2])**2)
                if value<float(args.mindist):
                    Error_list.append("Error")
        if len(Error_list)>0:
            i+=0
            continue
        elif len(Error_list)==0:
            for j in tempO_atoms:
                O_atoms.append(j)
                i+=1
    else:
        convex_close = []
        for j in convex_points:
            dist = np.abs((newxy[0]-j[0])**2+(newxy[1]-j[1])**2)
            if dist<0.001:
                convex_close.append(dist)
        if len(convex_close) == 0:
            dotp = []
            for k in range(len(convex_points)):
                for l in range(len(convex_points)):
                    if k!=l:
                        if dist2(convex_points[k],convex_points[l],[newxy[0],newxy[1]]) == True:
                            dotp.append(1)
            if len(dotp)==0:
                i+=0
                continue
            else:
                z_preliminary = random.uniform(float(args.mindist), float(args.maxdist))
                atom_distances = []
                Error_list = []
                for k in range(len(coordinates)):
                    atom_distances.append(np.sqrt((coordinates[k][1]-newxy[0])**2+(coordinates[k][2]-newxy[1])**2))
                closest_member = min(atom_distances)
                if closest_member > z_preliminary:
                    i+=0
                    continue
                atom_number = 0
                for k in range(len(coordinates)):
                    if closest_member == atom_distances[k]:
                        atom_number = k
                z_value = np.sqrt((z_preliminary)**2 - (coordinates[atom_number][1]-newxy[0])**2 - (coordinates[atom_number][2]-newxy[1])**2) + coordinates[atom_number][3]
                if math.isnan(z_value) == True:
                    i+=0
                    continue
                newO_atoms = []
                tempO_atoms = []
                for j in range(len(frames[How_Many_Systems[i]])):
                    tempO_atoms.append([frames[How_Many_Systems[i]].iloc[j]['Atoms'],frames[How_Many_Systems[i]].iloc[j]['X']+newxy[0],frames[How_Many_Systems[i]].iloc[j]['Y']+newxy[1],frames[How_Many_Systems[i]].iloc[j]['Z']+z_value])
                for k in O_atoms:
                    for j in distVals:
                        newO_atoms.append([k[1] + j[0], k[2] + j[1], k[3] + j[2]])    
                if len(newO_atoms)>0:
                    for k in range(len(newO_atoms)):
                        for l in range(len(tempO_atoms)):
                            value = (np.sqrt((newO_atoms[k][0]-tempO_atoms[l][1])**2+(newO_atoms[k][1]-tempO_atoms[l][2])**2+(newO_atoms[k][2]-tempO_atoms[l][3])**2))
                            if value<OO_distance:
                                Error_list.append("Error")
                if len(Error_list)>0:
                    i+=0
                    continue
                tempO_atoms = []
                for k in range(len(frames[How_Many_Systems[i]])):
                    tempO_atoms.append([frames[How_Many_Systems[i]].iloc[k]['Atoms'],frames[How_Many_Systems[i]].iloc[k]['X']+newxy[0],frames[How_Many_Systems[i]].iloc[k]['Y']+newxy[1],frames[How_Many_Systems[i]].iloc[k]['Z']+z_value])
                for k in range(len(multiSurface)):
                    for l in range(len(tempO_atoms)):
                        value = np.sqrt((multiSurface[k][0]-tempO_atoms[l][1])**2+(multiSurface[k][1]-tempO_atoms[l][2])**2+(multiSurface[k][2]-tempO_atoms[l][2])**2)
                        if value<float(args.mindist):
                            Error_list.append("Error")
                if len(Error_list)>0:
                    i+=0
                    continue
                elif len(Error_list)==0:
                    for j in tempO_atoms:
                        O_atoms.append(j)
                        i+=1
        else:
            z_preliminary = random.uniform(float(args.mindist), float(args.maxdist))
            atom_distances = []
            Error_list = []
            for k in range(len(coordinates)):
                atom_distances.append(np.sqrt((coordinates[k][1]-newxy[0])**2+(coordinates[k][2]-newxy[1])**2))
            closest_member = min(atom_distances)
            if closest_member > z_preliminary:
                i+=0
                continue
            atom_number = 0
            for k in range(len(coordinates)):
                if closest_member == atom_distances[k]:
                    atom_number = k
            z_value = np.sqrt((z_preliminary)**2 - (coordinates[atom_number][1]-newxy[0])**2 - (coordinates[atom_number][2]-newxy[1])**2) + coordinates[atom_number][3]
            if math.isnan(z_value) == True:
                i+=0
                continue
            newO_atoms = []
            tempO_atoms = []
            for j in range(len(frames[How_Many_Systems[i]])):
                tempO_atoms.append([frames[How_Many_Systems[i]].iloc[j]['Atoms'],frames[How_Many_Systems[i]].iloc[j]['X']+newxy[0],frames[How_Many_Systems[i]].iloc[j]['Y']+newxy[1],frames[How_Many_Systems[i]].iloc[j]['Z']+z_value])
            for k in O_atoms:
                for j in distVals:
                    newO_atoms.append([k[1] + j[0], k[2] + j[1], k[3] + j[2]])    
            if len(newO_atoms)>0:
                for k in range(len(newO_atoms)):
                    for l in range(len(tempO_atoms)):
                        value = (np.sqrt((newO_atoms[k][0]-tempO_atoms[l][1])**2+(newO_atoms[k][1]-tempO_atoms[l][2])**2+(newO_atoms[k][2]-tempO_atoms[l][3])**2))
                        if value<OO_distance:
                            Error_list.append("Error")
            if len(Error_list)>0:
                i+=0
                continue
            tempO_atoms = []
            for k in range(len(frames[How_Many_Systems[i]])):
                tempO_atoms.append([frames[How_Many_Systems[i]].iloc[k]['Atoms'],frames[How_Many_Systems[i]].iloc[k]['X']+newxy[0],frames[How_Many_Systems[i]].iloc[k]['Y']+newxy[1],frames[How_Many_Systems[i]].iloc[k]['Z']+z_value])
            for k in range(len(multiSurface)):
                for l in range(len(tempO_atoms)):
                    value = np.sqrt((multiSurface[k][0]-tempO_atoms[l][1])**2+(multiSurface[k][1]-tempO_atoms[l][2])**2+(multiSurface[k][2]-tempO_atoms[l][2])**2)
                    if value<float(args.mindist):
                        Error_list.append("Error")
            if len(Error_list)>0:
                i+=0
                continue
            elif len(Error_list)==0:
                for j in tempO_atoms:
                    O_atoms.append(j)
                    i+=1
new_O_atoms = []

'''
Place the molecules from the .txt files in the positions randomly generated on the surface. 
'''

for i in range(O_atomslen,len(O_atoms)):
    X = O_atoms[i][1]
    Y = O_atoms[i][2]
    Z = O_atoms[i][3]
    new_O_atoms.append([O_atoms[i][0],X,Y,Z])

f = open(args.output, 'w')
print(len(df1)+len(new_O_atoms),file=f)
print("Surface Generated",file=f)
for k in range(len(df1)):
    print(df1.iloc[k]['Atoms'],df1.iloc[k]['X'],df1.iloc[k]['Y'],df1.iloc[k]['Z'],file=f)
for i in new_O_atoms:
    print(i[0], i[1], i[2], i[3],file=f)
f.close()

if args.plot == 'Y':
    minx = polygonsurface[0][0]
    miny = polygonsurface[0][1]
    maxx = polygonsurface[0][0]
    maxy = polygonsurface[0][1]
    for i in polygonsurface:
        if i[0]<minx:
            minx = i[0]
        if i[0]>maxx:
            maxx = i[0]
        if i[1]<miny:
            miny = i[0]
        if i[1]>maxy:
            maxy = i[1]
    minmaxX = maxx - minx
    minmaxY = maxy - miny
    wh = minmaxX/minmaxY
    X1 = []
    Y1= []
    for i in convex_points:
        X1.append(i[0])
        Y1.append(i[1])
    plt.scatter(X1,Y1,color='blue',lw=3)
    plt.plot(X1,Y1,color='blue',lw=3)
    for i in O_atoms:
        plt.scatter(i[1],i[2],color='red',edgecolors='black',s=250,lw=2)
    plt.gcf().set_size_inches(14.4*wh, 14.4)
    plt.grid(False)
    plt.xlabel("x-coordinates",fontsize=20)
    plt.ylabel("y-coordinates",fontsize=20)
    plt.yticks(fontsize=20)
    plt.xticks(fontsize=20)
    plt.show()
