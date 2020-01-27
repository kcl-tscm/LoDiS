from CNA import *
from matplotlib import pyplot as plt

""" Armand

This is an example file to show how to use the CNA python descriptor to
obtain CNA histograms. This currently works, with any xyz file, and saves the
histograms to a png.
Drawbacks: The r_cut is provided at the start, without changing per frame which
is a no no. Can lead to manyyyy issues, and honestly the Master function should
ALWAYS be sorted before returning, which I will probably try to change in
the github. However, this perfectly works for now, and the code can be stolen
and put anywhere you want to get CNA graphs.

"""

filename = "test_file.xyz"
r_cut = 3.0

MasterKey=Master(filename, r_cut)
MasterKey=sorted(MasterKey)

""" Armand

MasterKey is NOT ordered when produced by CNA, while CNA_Data IS.
This sorted line is to ordered the MasterKey to fit the CNA_Data, to plot the
right data.

"""

End: int = 3
Start: int = 0
Step: int = 1

CNA_Data = CNA_Sig_Frame(filename, MasterKey, R_Cut = r_cut, Frames = End,
                        Skip = Step)

for i in range(End-Start):
    plt.figure(i)
    plt.title('CNAs {} frame:{}'.format(filename[:-4],i))
    plt.xlabel('CNA Classes')
    plt.xticks(np.arange(len(MasterKey)),labels=MasterKey, rotation=-45,
                ha="left", rotation_mode="anchor")
    plt.ylabel('Percentage')
    plt.bar(np.arange(len(MasterKey)), CNA_Data[i], color='r')
    plt.savefig('CNAs_{}_frame_{}.png'.format(filename[:-4],i),
                bbox_inches='tight')
