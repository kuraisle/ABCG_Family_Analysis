# Functions for importing the alignment and computing conservation patterns

# read_fasta -- my own little function for reading fasta files output by the MAFFT server
# a_a_entropy -- calculates the entropy for a list of single letter amino acid codes
# find_most_common -- returns the most common value of a list. If there's more than one value meeting that criterion, it will return them both, but that won't be an issue
# group_conservation_column -- for a column, calculates which groups are conserved
# group_conservation_alignment -- basically a wrapper for group_conservation_column
# remove_alignment_gaps -- maps sequence positions to columns in the sequence alignment
# conserved_search -- lets you search through a list of conservation patterns for a group or grouping

def read_fasta(path):
  aligned_sequences = []
  with open(path, 'r') as f:
    aligned_sequences = f.read().split('>')

  aligned_sequences = [line.split('\n') for line in aligned_sequences]
  seq_for_convert = []

  for i in range(1,len(aligned_sequences)):
    tmp_str = ''.join(aligned_sequences[i][1:])
    tmp = (aligned_sequences[i][0], tmp_str)
    seq_for_convert.append(tmp)
  return seq_for_convert

def a_a_entropy(group, chars = {'A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y', '-'}):
  counts = []
  for char in chars:
    count = len([x for x in group if x.upper() == char])
    counts.append((char, count))
  frac = [(x, y/len(group)) for (x,y) in counts if y != 0]
  entropy = -sum([y*math.log2(y) for (x,y) in frac])
  return entropy

def find_most_common(some_list):
  uniques = set(some_list)
  counts = []
  for val in uniques:
    count = len([x for x in some_list if x == val])
    counts.append((val, count))
  max_count = max([y for (x,y) in counts])
  most_common = [x for (x,y) in counts if y == max_count]
  return most_common[0]

def group_conservation_column(column_groups, gap_cutoff = 0.1, group_gap_cutoff = 0.3, entropy_cutoff = 2/3):
  full_groups = []
  for group in column_groups:
    full_groups = full_groups + group[1]
  total_gap_frac = len([x for x in full_groups if x == '-'])/len(full_groups)
  if total_gap_frac >= gap_cutoff:
    return 'Gap'
  else:
    group_gap_fracs = []
    for group in column_groups:
      gap_frac = len([x for x in group if x == '-'])/len(group)
      group_gap_fracs.append(gap_frac)
    if max(group_gap_fracs) >= group_gap_cutoff:
      return 'Gap'
    else:
      position_entropy = a_a_entropy(full_groups)
      if position_entropy < 2/3:
        return 'Totally Conserved'
      else:
        group_conservation = [(group[0], a_a_entropy(group[1])) for group in column_groups]
        if min([y for (x,y) in group_conservation]) > entropy_cutoff:
          return 'Not conserved'
        else:
          conserved_groups = set([group[0] for group in group_conservation if group[1] < entropy_cutoff])
          conserved_a_a = [(group[0], find_most_common(group[1])) for group in column_groups if group[0] in conserved_groups]
          conserved_set = set([group[1] for group in conserved_a_a])
          conservation_pattern = [(a_a, [group[0] for group in conserved_a_a if group[1] == a_a])for a_a in conserved_set]
          return conservation_pattern

def group_conservation_alignment(group_list, gap_cutoff = 0.1, group_gap_cutoff = 0.3, entropy_cutoff = 2/3):
  descriptor_output = []
  align_length = len(group_list[0][1][1][1])
  for i in range(align_length):
    column_list = [(group[0],[seq[1][i] for seq in group[1]]) for group in group_list]
    descriptor = group_conservation_column(column_list, gap_cutoff, group_gap_cutoff, entropy_cutoff)
    descriptor_output.append((i, descriptor))
  return descriptor_output

def remove_alignment_gaps(aligned_sequence, initial_value = 1):
  real_aligned_pairing = []
  for i in range(len(aligned_sequence)):
    if aligned_sequence[i] != '-':
      pairing = (len(real_aligned_pairing) + initial_value,i)
      real_aligned_pairing.append(pairing)
  return real_aligned_pairing

def conserved_search(conserved_list, search_term):
  out = []
  for column in conserved_list:
    pattern = column[1]
    if type(search_term) == list:
      for part in pattern:
        if part[1] == search_term:
          out.append(column)
    elif type(search_term) == str:
      for part in pattern:
        if search_term in part[1]:
          out.append(column)
  return out