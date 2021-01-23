# The functions for generating protein sequence logos summarising the sequence alignment
# protein_logo draws the sequence logo
# conserved_colours generates colours for the sequence logo
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

  fig = plt.figure(figsize=[0.5*len(ABCG1[0]), 9])
  
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
  
  plt.xticks(rotation = 70, ha = 'right')
  this_conservation_pattern = []
  for i in positions:
    this_conservation_pattern.append(conservation_pattern[i])
  ax4.set_xticklabels(this_conservation_pattern)

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
  
  fig.tight_layout()

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