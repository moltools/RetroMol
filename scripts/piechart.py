import argparse 
import matplotlib.pyplot as plt

sizes = {
    "cyclic peptides": 1733,
    "linear peptides": 783,
    "open-chain polyketides": 765,
    "cyclic peptides & lipopeptides": 467,
    "macrolide lactones": 279,
    "polyene macrolides": 262,
    "linear peptides & lipopeptides": 252,
    "cyclic peptides & microcystins": 245,
    "cyanopeptolins & cyclic peptides": 131,
    "bafilomycins": 93,
    "erythromycins": 91,
    "tylosins": 90,
    "antimycins": 66,
    "cyanopeptolins, cyclic peptides & lipopeptides": 51,
    "lipopeptides": 49,
    "epothilones": 48,
    "lactam bearing macrolide lactones": 34,
    "erythromycins & tylosins": 12,
    "cyclic peptides, linear peptides & lipopeptides": 9,
    "mycrocysints": 6,
    "cyclic peptides & open-chain polyketides": 1,
    "macrolide lactones & open-chain polyketides": 1,
    "macrolide lactones & polyene macrolides": 1,
}

parser = argparse.ArgumentParser()
parser.add_argument("-o", type=str, required=True, help="Output file name")
args = parser.parse_args()

# make pie chart
fig, ax = plt.subplots()
ax.pie(sizes.values(), labels=sizes.keys(), autopct='%1.1f%%')
ax.axis('equal')  # Equal aspect ratio ensures that pie is drawn as a circle.   
# add outer line
plt.gca().add_artist(plt.Circle((0,0),0.70,fc='white'))
# black line around the circle
plt.gca().add_artist(plt.Circle((0,0),0.70,fc='none',ec='black'))
plt.gca().add_artist(plt.Circle((0,0),1.0,fc='none',ec='black'))
plt.savefig(args.o)