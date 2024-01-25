import torch
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.path import Path
from sklearn.decomposition import PCA
import seaborn as sns
import logomaker



### load data
# Mason
pt = torch.load('mason_unique_results.pt', map_location='cpu').numpy()
labels = torch.load('mason_unique_label.pt', map_location='cpu').numpy()
sequ = torch.load('mason_unique_sequence.pt', map_location='cpu').numpy()

# Brij
#pt = torch.load('brij_short_results.pt', map_location='cpu').numpy()
#labels = torch.load('brij_short_label.pt', map_location='cpu').numpy()
#sequ = torch.load('brij_short_sequence.pt', map_location='cpu').numpy()

### compute PCA
pca = PCA(n_components=2)

# load data
# (34146, 2)
pca_transform = pca.fit_transform(pt.reshape(-1, 20*10))

### plot PCA
x_axes = f'PCA1 explained variance: {pca.explained_variance_ratio_[0]*100:.2f}%'
y_axes = f'PCA2 explained variance: {pca.explained_variance_ratio_[1]*100:.2f}%'

X = pd.DataFrame(pca_transform, columns = [x_axes, y_axes])
X['Labels'] = pd.DataFrame(labels, columns=['Labels'])['Labels'].apply(lambda x: 'Binder' if x > 0.5 else 'Non Binder')

#x_b = X.query('Labels == "Binder"')
#x_nb = X.query('Labels == "Non Binder"')
ax = sns.scatterplot(x=x_axes, y=y_axes, data=X, linewidth=0, hue='Labels', s=5)
#ax = sns.scatterplot(x=x_axes, y=y_axes, data=x_nb, linewidth=0, s=5)

# ax.set_axis_off()
ax.set_frame_on(False)
# plt.title('PCA on Integraded Gradients for sequences with implanted signals')
plt.show()



### how you could extract specific sequences

#-0.0370, 0.0512#rechts unten
#-0.1263, 0.1838#links unten
#-0.0051, 0.1493#rechts oben
#-0.0779, 0.2759#links oben

# Mason negative extract
vertices = np.asarray([(-0.1263, 0.1838),#links unten
                       (-0.0370, 0.0512),#rechts unten
                       (-0.0051, 0.1493),#rechts oben
                       (-0.0779, 0.2759)])#links oben

# Mason spike right down extract
vertices = np.asarray([(-0.0735, -0.096),#links unten
                       (0.1508, -0.425),#rechts unten
                       (0.2074, -0.312),#rechts oben
                       (0.0211, -0.04)])#links oben


# Mason spike right down up
vertices = np.asarray([(0.0484, -0.008),#links unten
                       (0.3635, 0.15),#rechts unten
                       (0.3703, 0.295),#rechts oben
                       (0.0513, 0.121)])#links oben

# Brij spike down
vertices = np.asarray([(0.0578, -0.045),#links unten
                       (0.5968, 0.029),#rechts unten
                       (0.5942, 0.098),#rechts oben
                       (0.0628, 0.009)])#links oben


# Brij spike vertical 1
vertices = np.asarray([(-0.0263, 0.089),#links unten
                       (-0.0056, 0.093),#rechts unten
                       (-0.0784, 0.604),#rechts oben
                       (-0.0997, 0.588)])#links oben

# Brij spike diagonal
vertices = np.asarray([(0.1356, 0.127),#links unten
                       (0.2887, 0.264),#rechts unten
                       (0.2648, 0.298),#rechts oben
                       (0.1256, 0.147)])#links oben


### ploting selection

path = Path(vertices)
mask = path.contains_points(pca_transform)
X['mask'] = mask
ax = sns.scatterplot(x=x_axes, y=y_axes, data=X, linewidth=0, hue='mask', s=5)
ax.set_frame_on(False)
plt.show()

### extracting actual sequences
extracted = pca_transform[mask]
ex_sequ_neg = sequ[mask, :]
ex_sequ_pos = sequ[~mask, :]

ex_sequ_pos = sequ
# prepare encoding
alphabet = 'ACDEFGHIKLMNPQRSTVWY'
encoding = {letter:i for i, letter in enumerate(alphabet)}
reverse_encoding = {i:letter for i, letter in enumerate(alphabet)}

sequ_lst_neg = list()
for i in ex_sequ_neg:
    sequ_lst_neg.append([reverse_encoding[aa] for aa in i.tolist()])

sequ_lst_pos = list()
for i in ex_sequ_pos:
    sequ_lst_pos.append(''.join([reverse_encoding[aa] for aa in i.tolist()]))








df_out = pd.DataFrame(np.zeros((10, 20)), columns = [j for j in 'ACDEFGHIKLMNPQRSTVWY'])

df_tmp = pd.DataFrame(sequ_lst_neg)
for i in range(10):
    t = df_tmp[i].value_counts() / df_tmp[i].value_counts().sum()
    t = t.reset_index()
    for _, row in t.iterrows():
        df_out.loc[i, row['index']] = row[i]




#cdrh3_df = pd.DataFrame(gi[i, :, :], columns = [j for j in 'ACDEFGHIKLMNPQRSTVWY'])

### LOGO plot of sequence distribution

plt.clf()

# create Logo object
crp_logo = logomaker.Logo(df_out,
                          shade_below=.5,
                          fade_below=.5,
                          font_name='Arial Rounded MT Bold',
                          color_scheme='skylign_protein')# skylign_protein, dmslogo_funcgroup

# style using Logo methods
crp_logo.style_spines(visible=False)
crp_logo.style_spines(spines=['left', 'bottom'], visible=True)
crp_logo.style_xticks(rotation=0, fmt='%d', anchor=0)

# style using Axes methods
crp_logo.ax.set_ylabel("Frequency", labelpad=-1)
crp_logo.ax.set_xlabel("Position", labelpad=-1)
crp_logo.ax.set_xticklabels('%d'%(x+1) for x in range(10))
crp_logo.ax.xaxis.set_ticks_position('none')
crp_logo.ax.xaxis.set_tick_params(pad=-1)
# crp_logo.ax.set_yticks([-xn, -xn/2, 0, xn/2, xn])

# crp_logo.fig.savefig('books_read.png')
#seq = "".join([reverse_encoding[repr(data[i][j])] for j in range(11)])
#crp_logo.fig.suptitle(f'{seq}    {"Non-Binding" if label[i] < 0.5 else "Binding"}    Prediction: {error[i]:.2f}')

crp_logo.fig.show()
