import random
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.pyplot import figure
import matplotlib.ticker as ticker


# The participants will be divided into 2-3 groups. Each group will be given one of the following mutation functions:
# Group 1 mutator:
def mutate_sequence1(seq, index_to_mutate = [8, 19, 23, 39]): 

  new_seq = ""

  # randomly choose index
  idx = random.choice(index_to_mutate)

  # ignore the first 3 and last 3 base
  for i in range(3, len(seq)-3):  

    base = seq[i]

    if base == "C": 
      weights = [0.2, 0.2, 0.5, 0.1]
    elif base == "G": 
      weights = [0.2, 0.2, 0.1, 0.5]
    elif base == "A": 
      weights = [0.1, 0.5, 0.2, 0.2]
    else: 
      weights = [0.5, 0.1, 0.2, 0.2]

    new_base = random.choices(["A", "T", "G", "C"], weights, k=1)[0]

    if i == idx: 
      new_seq += new_base
    # add some noise that there is a 0.05 chance of randomly mutating any other base
    else: 
      if random.random() < 0.05: 
        new_seq += new_base
      else: 
        new_seq += base

  return(seq[:3] + new_seq + seq[len(seq)-3:])

# Group 2 mutator:
def mutate_sequence2(seq, index_to_mutate = [23, 32, 46]): 

  new_seq = ""

  # randomly choose index
  idx = random.choice(index_to_mutate)

  # ignore the first 3 and last 3 base
  for i in range(3, len(seq)-3):  

    base = seq[i]

    if base == "C": 
      weights = [0.2, 0.5, 0.2, 0.1]
    elif base == "G": 
      weights = [0.5, 0.2, 0.1, 0.2]
    elif base == "A": 
      weights = [0.1, 0.2, 0.5, 0.2]
    else: 
      weights = [0.2, 0.1, 0.2, 0.5]

    new_base = random.choices(["A", "T", "G", "C"], weights, k=1)[0]

    if i == idx: 
      new_seq += new_base
    # add some noise that there is a 0.05 chance of randomly mutating any other base
    else: 
      if random.random() < 0.05: 
        new_seq += new_base
      else: 
        new_seq += base

  return(seq[:3] + new_seq + seq[len(seq)-3:])


# Next, the participants will explore the extent of the mutations across
# simulations by writing a function to calculate hamming distance. They can
# then implement the following function to visualize the distribution of
# hamming distances across simulations:
def plot_hamming_distribution(hamming_distances, title):
  unique_dists = sorted(list(set(hamming_distances)))
  dist_counts = {}

  for dist in unique_dists:
    dist_counts[dist] = hamming_distances.count(dist)

  figure(figsize=(10, 6), dpi=80)
  plt.bar(dist_counts.keys(), dist_counts.values(), width=0.9)
  plt.xlabel("Hamming distance")
  plt.title(title)
  x = np.arange(0, max(unique_dists)+1, 1)
  plt.xticks(x)
  plt.show()


# Next, participants will try to figure out what the mutator did (using our provided plotting functions)
# Visualize mutations at the DNA level:
def plot_nt_frequencies(orig_seq, sequences): 

  values = np.array([np.array(list(seq)) for seq in sequences])
  bases = ["A", "T", "G", "C"]

  counts = np.zeros((len(bases), len(sequences[1])))
  original = np.zeros((2, len(sequences[1])))
  original[0, ] = range(len(orig_seq))
  original_bases = []

  for i in range(counts.shape[1]): 
    for j, base in enumerate(bases): 
      if base == orig_seq[i]:
        counts[j,i] = 0
        original_bases.append(bases.index(base))
      else: 
        counts[j,i] = list(values[:,i]).count(base) / len(sequences)
  original[1, ] = original_bases

  fig, ax = plt.subplots(1,1)
  fig.set_figheight(3)
  fig.set_figwidth(30)

  ax.set_xticks(list(range(len(sequences))))
  ax.set_xticklabels(list(orig_seq))
  ax.set_yticks(list(range(len(bases))))
  ax.set_yticklabels(bases)

  # ax2 = ax.twiny()
  # ax2.spines["bottom"].set_position(("axes", -0.10))
  # ax2.tick_params('both', length=0, width=0, which='minor')
  # ax2.tick_params('both', direction='in', which='major')
  # ax2.xaxis.set_ticks_position("bottom")
  # ax2.xaxis.set_label_position("bottom")
  
  # ax2.set_xticks(list(range(0,51, 3)))
  # ax2.xaxis.set_major_formatter(ticker.NullFormatter())
  # ax2.xaxis.set_minor_locator(ticker.FixedLocator([0.3, 0.8]))
  # ax2.xaxis.set_minor_formatter(ticker.FixedFormatter(['mammal', 'reptiles']))


  color_map = plt.imshow(counts, aspect='auto')
  color_map.set_cmap("Reds")
  plt.colorbar(fraction=0.01, pad=0.01)

  for i in range(3,len(orig_seq),3):
    plt.axvline(x = i-0.5, color="gray")
  plt.scatter(original[0,], original[1,], color='black')
  plt.show()

# Visualize mutations at the amino acid level:
def reformat_protein_sims(orig_protein, protein_sequences):
  # if not full length: 
  reformatted_seqs = []
  for protein in protein_sequences:

    if len(protein) != len(orig_protein):

      # replace protein with 
      protein += "".join(["x"]*(len(orig_protein)-len(protein)))

    reformatted_seqs.append(protein)
  return(reformatted_seqs)

def plot_AA_frequencies(orig_protein, protein_sequences):
  #reformat sims for plotting
  formatted_protein_sequences = reformat_protein_sims(orig_protein, protein_sequences)

  AAs = ['A','C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y', '_']
  values = np.array([list(protein) for protein in formatted_protein_sequences])
  counts = np.zeros((len(AAs), len(formatted_protein_sequences[1])))

  original = np.zeros((2, len(formatted_protein_sequences[1])))
  original[0, ] = range(len(orig_protein))
  original_AAs = []

  for i in range(counts.shape[1]):
    for j, aa in enumerate(AAs):
      counts[j,i] = list(values[:,i]).count(aa)
      if orig_protein[i] == aa:
        counts[j,i] = 0
        original_AAs.append(AAs.index(aa))
      else:
        counts[j,i] = list(values[:,i]).count(aa) / len(formatted_protein_sequences)
  original[1, ] = original_AAs

  fig, ax = plt.subplots(1,1)
  fig.set_figheight(10)
  fig.set_figwidth(35)

  ax.set_xticks(list(range(len(formatted_protein_sequences))))
  ax.set_xticklabels(list(orig_protein))
  ax.set_yticks(list(range(len(AAs))))
  ax.set_yticklabels(AAs)

  color_map = plt.imshow(counts, aspect='auto')
  color_map.set_cmap("Reds")
  plt.colorbar(pad=0.01)
  for i in range(1,len(orig_protein)):
    plt.axvline(x = i-0.5, color="gray")
  plt.scatter(original[0,], original[1,], color='black')
  plt.show()
