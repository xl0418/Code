from mpl_toolkits.mplot3d.axes3d import Axes3D
from mpl_toolkits.mplot3d import proj3d
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
from collections import OrderedDict
import pandas as pd
import seaborn as sns
#Make sure these are floating point values:
scale_x = 1.0
scale_y = 2.5
scale_z = 1.0

#Axes are scaled down to fit in scene
max_scale=max(scale_x, scale_y, scale_z)

scale_x=scale_x/max_scale
scale_y=scale_y/max_scale
scale_z=scale_z/max_scale

#Create scaling matrix
scale = np.array([[scale_x,0,0,0],
                  [0,scale_y,0,0],
                  [0,0,scale_z,0],
                  [0,0,0,1]])
print(scale)

def get_proj_scale(self):
    """
    Create the projection matrix from the current viewing position.

    elev stores the elevation angle in the z plane
    azim stores the azimuth angle in the x,y plane

    dist is the distance of the eye viewing point from the object
    point.

    """
    relev, razim = np.pi * self.elev/180, np.pi * self.azim/180

    xmin, xmax = self.get_xlim3d()
    ymin, ymax = self.get_ylim3d()
    zmin, zmax = self.get_zlim3d()

    # transform to uniform world coordinates 0-1.0,0-1.0,0-1.0
    worldM = proj3d.world_transformation(
        xmin, xmax,
        ymin, ymax,
        zmin, zmax)

    # look into the middle of the new coordinates
    R = np.array([0.5, 0.5, 0.5])

    xp = R[0] + np.cos(razim) * np.cos(relev) * self.dist
    yp = R[1] + np.sin(razim) * np.cos(relev) * self.dist
    zp = R[2] + np.sin(relev) * self.dist
    E = np.array((xp, yp, zp))

    self.eye = E
    self.vvec = R - E
    self.vvec = self.vvec / proj3d.mod(self.vvec)

    if abs(relev) > np.pi/2:
    # upside down
      V = np.array((0, 0, -1))
    else:
      V = np.array((0, 0, 1))
    zfront, zback = -self.dist, self.dist

    viewM = proj3d.view_transformation(E, R, V)
    perspM = proj3d.persp_transformation(zfront, zback)
    M0 = np.dot(viewM, worldM)
    M = np.dot(perspM, M0)

    return np.dot(M, scale);

Axes3D.get_proj=get_proj_scale

"""
You need to include all the code above.
From here on you should be able to plot as usual.
"""




dir = 'C:/Liang/Googlebox/Research/Project3/batchsim_results/1e+07/'
scenario = np.array(['LR','MR','HR'])
sce_short = np.array(['L','M','H'])
spatiallevel = np.array(['low','intermediate','high'])
jclabel = np.array(['0','0.5','1'])
plabel = np.array(['0','1e2','1e4','1e6','1e8','Inf'])
cmaps = OrderedDict()

i_sce = 2
csv_name = scenario[i_sce]
folder_name = sce_short[i_sce]

# Compare JC effect when fixing phylogenetic strength
color_group = np.tile(['r','g','b'],7)
for j in range(1,7):
    group = 1
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    nbins = 12
    width = .9
    for i in range(1,4):
        abundance_list = np.zeros(12)
        for rep in range(61,101):
            filename=dir+'spatialpara1e+07%s%i%i/%s%i%irep%i.csv' % (folder_name,i,j,csv_name,i,j,rep)
            df = np.log(pd.read_csv(filename, header=None).values[0]).tolist()
            df_hist, df_bins = np.histogram(df, bins=[i for i in range(13)])
            abundance_list = np.vstack((abundance_list,df_hist))
        abundance_list = abundance_list[1:,:]
        quantile_col = np.zeros((5,1))
        for col in range(12):
            quantile_col = np.hstack((quantile_col,np.quantile(abundance_list[:,col],q=[0,.25,.50,.75,1]).reshape(5,1)))
        quantile_col = quantile_col[:,1:]

        xs = (df_bins[:-1] + df_bins[1:]) / 2
        ax.bar(xs, quantile_col[2,:], zs=group*10,width=width, zdir='y', color=color_group[group], alpha=0.6)
        for no_bars in range(12):
            ax.plot([xs[no_bars], xs[no_bars]], [group*10, group*10], [quantile_col[1,no_bars],
                   quantile_col[3,no_bars]], marker="_",color = color_group[group])

        group += 1

    ax.set_yticks([10,20,30])
    ax.set_yticklabels([0,0.5,1])
    ax.set_xlabel('Abundance (log)')
    ax.set_ylabel('JC strength')
    ax.set_zlabel('Frequencey')
    ax.set_zlim(0,40)
    title = 'Abundance distribution under $\sigma_{\phi}=%s$\nwhen spatial level is %s' % (plabel[j-1],spatiallevel[i_sce])
    fig.suptitle(title, fontsize=16)

    outputname = 'C:/Liang/Googlebox/Research/Project3/batchsim_results/AbundanceDisSp%sP%s.png' % (folder_name,plabel[j-1])
    fig.savefig(outputname)
    plt.close(fig)

# Compare phylogenetic effect when fixing JC strength
color_group_p = np.array([sns.color_palette("Paired")[1],sns.color_palette("Paired")[5]])

for i in range(1,4):
    group = 1
    fig = plt.figure(figsize=(10, 6))
    fig.subplots_adjust(left=0.1, right=1.1,bottom=0,top=1.2)
    ax = fig.add_subplot(111, projection='3d')
    nbins = 12
    width = .9
    for j in range(1,7):
        abundance_list = np.zeros(12)
        for rep in range(61,101):
            filename=dir+'spatialpara1e+07%s%i%i/%s%i%irep%i.csv' % (folder_name,i,j,csv_name,i,j,rep)
            df = np.log(pd.read_csv(filename, header=None).values[0]).tolist()
            df_hist, df_bins = np.histogram(df, bins=[i for i in range(13)])
            abundance_list = np.vstack((abundance_list,df_hist))
        abundance_list = abundance_list[1:,:]
        quantile_col = np.zeros((5,1))
        for col in range(12):
            quantile_col = np.hstack((quantile_col,np.quantile(abundance_list[:,col],q=[0,.25,.50,.75,1]).reshape(5,1)))
        quantile_col = quantile_col[:,1:]

        xs = (df_bins[:-1] + df_bins[1:]) / 2
        ax.bar(xs, quantile_col[2,:], zs=group*10,width=width, zdir='y', color=color_group_p[group%2], alpha=0.8)
        for no_bars in range(12):
            ax.plot([xs[no_bars], xs[no_bars]], [group*10, group*10], [quantile_col[1,no_bars],
                   quantile_col[3,no_bars]], marker="_",color = color_group_p[group%2])

        group += 1

    ax.set_yticks([10,20,30,40,50,60])
    ax.set_yticklabels(plabel)
    ax.set_xlabel('Abundance (log)',fontsize=15)
    ax.set_ylabel('Phylogenetic strength' ,fontsize=15)
    ax.set_zlabel('Frequencey',fontsize=15)
    ax.xaxis.labelpad = 15
    ax.yaxis.labelpad = 25
    ax.zaxis.labelpad = 15

    ax.set_zlim(0,40)
    title = 'Abundance distribution under $\psi=%s$\nwhen spatial level is %s' % (jclabel[i-1],spatiallevel[i_sce])
    fig.suptitle(title, fontsize=16)

    outputname = 'C:/Liang/Googlebox/Research/Project3/batchsim_results/AbundanceDisSp%sJC%s.png' % (folder_name,jclabel[i-1])
    fig.savefig(outputname)
    plt.close(fig)

