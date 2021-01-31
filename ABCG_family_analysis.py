import streamlit as st
import math
import pandas as pd
import numpy as np
import logomaker as lm
import matplotlib.pyplot as plt

st.title('Functional Divergence in Mammalian ABCGs')

@st.cache
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

@st.cache
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

MAFFT_alignment = read_fasta('20200806 MAFFT ABCG small low gap.fasta')
ABCG1_sequences = MAFFT_alignment[:35]
ABCG4_sequences = MAFFT_alignment[35:71]
ABCG2_sequences = MAFFT_alignment[71:103]
ABCG5_sequences = MAFFT_alignment[103:-35]
ABCG8_sequences = MAFFT_alignment[-35:]

alignment_protein_conservation = group_conservation_alignment([('ABCG1', ABCG1_sequences), ('ABCG2', ABCG2_sequences), ('ABCG4', ABCG4_sequences), ('ABCG5', ABCG5_sequences), ('ABCG8', ABCG8_sequences)])
functionally_divergent = [x for x in alignment_protein_conservation if x[1] not in ['Totally Conserved', 'Gap', 'Not conserved']]

def protein_logo(positions):
  """Draws a sequence logo of the positions requested showing differences between ABCG family members

  Arguments:
  positions -- a list of positions within the sequence alignment
  """
  ABCG1 = []
  ABCG2 = []
  ABCG4 = []
  ABCG5 = []
  ABCG8 = []

  for seq in ABCG1_sequences:
    tmp = ''
    for i in positions:
      tmp = tmp + seq[1][i]
    ABCG1.append(tmp)

  for seq in ABCG2_sequences:
    tmp = ''
    for i in positions:
      tmp = tmp + seq[1][i]
    ABCG2.append(tmp)

  for seq in ABCG4_sequences:
    tmp = ''
    for i in positions:
      tmp = tmp + seq[1][i]
    ABCG4.append(tmp)

  for seq in ABCG5_sequences:
    tmp = ''
    for i in positions:
      tmp = tmp + seq[1][i]
    ABCG5.append(tmp)

  for seq in ABCG8_sequences:
    tmp = ''
    for i in positions:
      tmp = tmp + seq[1][i]
    ABCG8.append(tmp)

  fig = plt.figure(figsize=[0.5*len(ABCG1[0]), 5])
  
  ax = plt.subplot2grid((5, 1), (0,0))
  ABCG1_logo = lm.Logo(lm.alignment_to_matrix(ABCG1), ax = ax, color_scheme='black')
  ax.set_xticks(range(len(positions)))
  ax.set_xticklabels(positions)
  ax.xaxis.tick_top()
  ax1 = plt.subplot2grid((5, 1), (1,0))
  ABCG2_logo = lm.Logo(lm.alignment_to_matrix(ABCG2), ax = ax1, color_scheme='black')
  ax1.set_xticks([])
  ax2 = plt.subplot2grid((5, 1), (2,0))
  ABCG4_logo = lm.Logo(lm.alignment_to_matrix(ABCG4), ax = ax2, color_scheme='black')
  ax2.set_xticks([])
  ax3 = plt.subplot2grid((5, 1), (3,0))
  ABCG5_logo = lm.Logo(lm.alignment_to_matrix(ABCG5), ax = ax3, color_scheme='black')
  ax3.set_xticks([])
  ax4 = plt.subplot2grid((5, 1), (4,0))
  ABCG8_logo = lm.Logo(lm.alignment_to_matrix(ABCG8), ax = ax4, color_scheme='black')
  ax4.set_xticks(range(len(positions)))
  
  plt.xticks(rotation = 45, ha = 'right')
  this_conservation_pattern = []
  for i in positions:
    this_conservation_pattern.append(conservation_pattern[i])
  ax4.set_xticklabels(this_conservation_pattern)
  ax4.tick_params(labelsize = 8)

  ax.set_yticks([])
  ax1.set_yticks([])
  ax2.set_yticks([])
  ax3.set_yticks([])
  ax4.set_yticks([])

  ax.set_ylabel('ABCG1', rotation = 0, ha = 'right', fontsize = 20)
  ax1.set_ylabel('ABCG2', rotation = 0, ha = 'right', fontsize = 20)
  ax2.set_ylabel('ABCG4', rotation = 0, ha = 'right', fontsize = 20)
  ax3.set_ylabel('ABCG5', rotation = 0, ha = 'right', fontsize = 20)
  ax4.set_ylabel('ABCG8', rotation = 0, ha = 'right', fontsize = 20)

  conservation_colours = conserved_colours(positions)

  for pos in range(len(conservation_colours[0])):
    ABCG1_logo.highlight_position(p = pos, color = conservation_colours[0][pos])
  
  for pos in range(len(conservation_colours[0])):
    ABCG2_logo.highlight_position(p = pos, color = conservation_colours[1][pos])
  
  for pos in range(len(conservation_colours[0])):
    ABCG4_logo.highlight_position(p = pos, color = conservation_colours[2][pos])
  
  for pos in range(len(conservation_colours[0])):
    ABCG5_logo.highlight_position(p = pos, color = conservation_colours[3][pos])

  for pos in range(len(conservation_colours[0])):
    ABCG8_logo.highlight_position(p = pos, color = conservation_colours[4][pos])
  
  return fig

def conserved_colours(positions):
  """Generate colour palletes for positions in sequence logos according to their conservation

  Arguments:
  positions -- a list of positions within the sequence alignment 
  """
  patterns = []

  for i in positions:
    patterns.append(alignment_protein_conservation[i][1])
    
  colour_scheme = []

  for protein in ['ABCG1', 'ABCG2', 'ABCG4', 'ABCG5', 'ABCG8']:
    prot_colours = []
    for pos in patterns:
      if pos == 'Gap':
        prot_colours.append('grey')
      elif pos == 'Totally Conserved':
        prot_colours.append('darkslategray')
      else:
        pat_list = [x[1] for x in pos]
        representative = [member for sublist in pat_list for member in sublist]
        if protein in representative:
          if len(pat_list) == 1:
            prot_colours.append('green')
          elif len(representative) == 5:
            prot_colours.append('aqua')
          else:
            prot_colours.append('red')
        else:
          prot_colours.append('white')
    colour_scheme.append(prot_colours)
  return colour_scheme

conservation_pattern = []

for pos in alignment_protein_conservation:
  if pos[1] == 'Gap':
    conservation_pattern.append('Gap')
  elif pos[1] == 'Totally Conserved':
    conservation_pattern.append('Totally Conserved')
  elif type(pos[1]) == list:
    pattern = [set(x[1]) for x in pos[1]]
    pattern_str = str(pattern)
    conservation_pattern.append(pattern_str.replace("'", '').replace('[', '').replace(']', ''))
  else:
    conservation_pattern.append(pos)

human_seq_map = {
    'Human ABCG1': remove_alignment_gaps(MAFFT_alignment[19][1]),
    'Human ABCG4': remove_alignment_gaps(MAFFT_alignment[35][1]),
    'Human ABCG2': remove_alignment_gaps(MAFFT_alignment[90][1]),
    'Human ABCG5': remove_alignment_gaps(MAFFT_alignment[103][1]),
    'Human ABCG8': remove_alignment_gaps(MAFFT_alignment[139][1])
}

@st.cache
def position_extract(pos_string):
    out = []
    for pos in pos_string.split(','):
        if '-' in pos:
            rangevals = pos.split('-')
            out = out + list(range(int(rangevals[0]), int(rangevals[1])+1))
        else:
            out.append(int(pos))
    return out

seq_map_df = pd.DataFrame(
    human_seq_map['Human ABCG1'], columns = ['ABCG1', 'Alignment']
).merge(
    pd.DataFrame(
        human_seq_map['Human ABCG4'], columns = ['ABCG4', 'Alignment']
    ),
    on = 'Alignment'
).merge(
    pd.DataFrame(
        human_seq_map['Human ABCG2'], columns = ['ABCG2', 'Alignment']
    ),
    on = 'Alignment'
).merge(
    pd.DataFrame(
        human_seq_map['Human ABCG5'], columns = ['ABCG5', 'Alignment']
    ),
    on = 'Alignment'
).merge(
    pd.DataFrame(
        human_seq_map['Human ABCG8'], columns = ['ABCG8', 'Alignment']
    ),
    on = 'Alignment'
)

# Display components

## Sidebar
search_or_pattern = st.sidebar.selectbox(
    'Do you want to search by position or explore conservation patterns?',
    ['Search by position', 'Explore conservation patterns']
)

if search_or_pattern == 'Search by position':
    relative_position = st.sidebar.selectbox(
        'What would you like your positions relative to?',
        ['Alignment', 'Human ABCG1', 'Human ABCG4', 'Human ABCG2', 'Human ABCG5', 'Human ABCG8']
    )
else:
    ABCG1_group = st.sidebar.selectbox(
        'Which group for ABCG1?',
        ['Group 1', 'Group 2', 'Group 3', 'Group 4', 'Group 5']
    )
    ABCG4_group = st.sidebar.selectbox(
        'Which group for ABCG4?',
        ['Group 1', 'Group 2', 'Group 3', 'Group 4', 'Group 5']
    )
    ABCG2_group = st.sidebar.selectbox(
        'Which group for ABCG2?',
        ['Group 1', 'Group 2', 'Group 3', 'Group 4', 'Group 5']
    )
    ABCG5_group = st.sidebar.selectbox(
        'Which group for ABCG5?',
        ['Group 1', 'Group 2', 'Group 3', 'Group 4', 'Group 5']
    )
    ABCG8_group = st.sidebar.selectbox(
        'Which group for ABCG8?',
        ['Group 1', 'Group 2', 'Group 3', 'Group 4', 'Group 5']
    )



## Main body
if search_or_pattern == 'Search by position':
    if relative_position == 'Alignment':
        logo_pos_string = st.text_input(
            'Which columns do you want a logo of? Separate values with commas. For ranges, use a hyphen',
            '890, 891-905'
        )
        position_list = position_extract(logo_pos_string)
    else:
        logo_pos_str = st.text_input(
            f'Which {relative_position} positions do you want a logo of? Separate values with commas. For ranges, use a hyphen',
            '400-410'
        )
        rel_pos_list = position_extract(logo_pos_str)
        position_list = [x[1] for x in human_seq_map[relative_position] if x[0] in rel_pos_list]



    st.write(protein_logo(position_list))
else:
    cons_pat_str = 'Conservation pattern: '
    grouping_in = [ABCG1_group, ABCG2_group, ABCG4_group, ABCG5_group, ABCG8_group]
    grouping_set = set(grouping_in)
    grouping_tup = [
        ('ABCG1', ABCG1_group),
        ('ABCG2', ABCG2_group),
        ('ABCG4', ABCG4_group),
        ('ABCG5', ABCG5_group),
        ('ABCG8', ABCG8_group)
    ]
    query_list = []
    for group in grouping_set:
        cons_pat_str = cons_pat_str + '('
        group_members = [x[0] for x in grouping_tup if x[1] == group]
        for member in group_members:
            cons_pat_str = cons_pat_str + member + ', '
        cons_pat_str = cons_pat_str + ')'
        query_list.append(group_members)
    
    if len(grouping_set) == 1:
        query_pos = [x[0] for x in alignment_protein_conservation if x[1] == 'Totally Conserved']
    else:
        query_results = functionally_divergent[:]
        for query in query_list:
            query_results = conserved_search(query_results, query)
        query_pos = [x[0] for x in query_results]
    query_seq_df = seq_map_df.loc[seq_map_df.Alignment.isin(query_pos)]

    st.write(cons_pat_str)
    st.write(query_seq_df.set_index('Alignment'))