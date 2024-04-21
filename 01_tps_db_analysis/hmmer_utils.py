import pandas as pd
import ast

class DomainHit:
  def __init__(self, input_list):
    (domain_name, domain_accession, domain_len, sequence_name, sequence_accession, sequence_len, sequence_evalue, sequence_score, sequence_bias, num, total_num, domain_cevalue, domain_ievalue, domain_score, 
domain_bias, domain_start, domain_end, alignment_start, alignment_end, envelope_start, envelope_end, accuracy, description) = input_list
    self.domain_name = domain_name
    self.domain_accession = domain_accession
    self.domain_len = int(domain_len)
    self.sequence_name = sequence_name
    self.sequence_accession = sequence_accession
    self.sequence_len = int(sequence_len)
    self.sequence_evalue = float(sequence_evalue)
    self.sequence_score = float(sequence_score)
    self.sequence_bias = float(sequence_bias)
    self.num = int(num)
    self.total_num = int(total_num)
    self.domain_cevalue = float(domain_cevalue)
    self.domain_ievalue = float(domain_ievalue)
    self.domain_score = float(domain_score)
    self.domain_bias = float(domain_bias)
    self.domain_start = int(domain_start)
    self.domain_end = int(domain_end)
    self.domain_span = self.domain_end - self.domain_start + 1
    self.domain_coverage = self.domain_span / self.domain_len
    self.alignment_start = int(alignment_start)
    self.alignment_end = int(alignment_end)
    self.envelope_start = int(envelope_start)
    self.envelope_end = int(envelope_end)
    self.accuracy = float(accuracy)
    self.description = description

  def __str__(self):
    return f'{self.domain_name} (start: {self.domain_start}, end: {self.domain_end}, cevalue: {self.domain_cevalue}, coverage: {self.domain_coverage})'
  
def reduce_hits(hits, criterion='default'):
    """
    Reduces the list of hits to non-overlaping list of hits.
    To select from overlaping hits, there are 3 criterions:
    cevalue: select hit with higher cevalue

    coverage: select hit with higher coverage

    default: if there is not big difference in coverage (<30%),
    select hit with higher cevalue, otherwise select hit
    with higher coverage
    """
    reduced_hits = []
    current_hit = None

    for hit in hits: 
        if current_hit is None:
            current_hit = hit
        else:
            if hit.envelope_start > current_hit.envelope_end:
                # If the current hit doesn't overlap with the new hit, add it to the reduced list
                reduced_hits.append(current_hit)
                current_hit = hit
            else:
                # If there is overlap, select the one with the lowest cevalue
                if criterion == 'cevalue':
                    if hit.domain_cevalue < current_hit.domain_cevalue:
                        current_hit = hit

                # If there is overlap, select the one with the highest coverage
                elif criterion == 'coverage':
                    if hit.domain_coverage > current_hit.domain_coverage:
                        current_hit = hit
                        
                # If there is overlap and not big difference in coverage, select the one with highest cevalue
                # Otherwise select the one with highest coverage
                elif criterion == 'default':
                    coverage_diff = hit.domain_coverage - current_hit.domain_coverage
                    if coverage_diff < 0.3:
                        if hit.domain_cevalue < current_hit.domain_cevalue:
                            current_hit = hit
                    else:
                        current_hit = hit

    # Add the last selected hit to the reduced list
    if current_hit is not None:
        reduced_hits.append(current_hit)

    return reduced_hits

def print_architecture(architecture):
    for hit in architecture:
        print(hit.domain_name, hit.envelope_start, hit.envelope_end,hit.domain_coverage, hit.domain_cevalue, hit.domain_ievalue)

def print_comparison(seq_name, results_architectures):
    cevalue_arch, coverage_arch, default_arch = results_architectures[seq_name]
    print('cevalue architecture:')
    print_architecture(cevalue_arch)
    print()
    print('coverage architecture')
    print_architecture(coverage_arch)
    print()
    print('default architecture')
    print_architecture(default_arch)


def get_domain_hit_list(file_name):
    """
    Load the hmmscan results and parse all the hits to a list of DomainHit
    """
    with open(file_name, 'r') as in_file:
        domain_hits = []
        for line in in_file.readlines():
            if line.startswith('#'):
                continue
            line_list = line.split()
            description = " ".join(line_list[22:])
            domain_hit = DomainHit(line_list[:22] + [description])
            domain_hits.append(domain_hit)

    return domain_hits

def get_reduced_architectures(domain_hits):
    """
    Returns dictionary where for every sequence, there is a list of 3 architecture lists:
    0: cevalue architecture
    1: coverage architecture
    2: default architecture
    """

    # Initialize dictionary with list of hits for every sequence
    sequences = [domain_hit.sequence_name for domain_hit in domain_hits]
    sequences = set(sequences)
    init_lists = [[] for i in range(len(sequences))]
    results = dict(zip(sequences, init_lists))

    # Put every hit to corresponding sequence
    for domain_hit in domain_hits:
        results[domain_hit.sequence_name].append(domain_hit)

    # Get the architectures
    results_architectures = dict(zip(sequences, init_lists))

    for seq, hits in results.items():
        
        # Sort the domain hits based on start in the sequence
        sorted_hits = sorted(hits, key=lambda x: x.envelope_start, reverse=False)

        # Get the architectures using different criterions
        cevalue_arch = reduce_hits(sorted_hits, criterion='cevalue')
        coverage_arch = reduce_hits(sorted_hits, criterion='coverage')
        default_arch = reduce_hits(sorted_hits, criterion='default')
        results_architectures[seq] = [cevalue_arch, coverage_arch, default_arch]

    return results_architectures

def get_architectures_df(results_architectures, criterion="default", supfam=False):
    # Simplify the architectures to list of the domain names in the order of occurance
    sequences = results_architectures.keys()
    sequences = set(sequences)
    init_lists = [[] for i in range(len(sequences))]

    short_cevalue_archs = dict(zip(sequences, init_lists))
    short_partial_cevalue_archs = dict(zip(sequences, init_lists)) # When domain is covered by less than 50%, put suffix 'partial'

    short_coverage_archs = dict(zip(sequences, init_lists))
    short_partial_coverage_archs = dict(zip(sequences, init_lists)) # When domain is covered by less than 50%, put suffix 'partial'

    short_default_archs = dict(zip(sequences, init_lists))
    short_partial_default_archs = dict(zip(sequences, init_lists)) # When domain is covered by less than 50%, put suffix 'partial'

    if not supfam:
        for seq, (cevalue_arch, coverage_arch, default_arch) in results_architectures.items():

            short_cevalue_archs[seq] = [hit.domain_accession for hit in cevalue_arch]
            short_partial_cevalue_archs[seq] = [hit.domain_accession if hit.domain_coverage > 0.5 else hit.domain_accession+'_partial' for hit in cevalue_arch]

            short_coverage_archs[seq] = [hit.domain_accession for hit in coverage_arch]
            short_partial_coverage_archs[seq] = [hit.domain_accession if hit.domain_coverage > 0.5 else hit.domain_accession+'_partial' for hit in coverage_arch]

            short_default_archs[seq] = [hit.domain_accession for hit in default_arch]
            short_partial_default_archs[seq] = [hit.domain_accession if hit.domain_coverage > 0.5 else hit.domain_accession+'_partial' for hit in default_arch]
    else:
        for seq, (cevalue_arch, coverage_arch, default_arch) in results_architectures.items():

            short_cevalue_archs[seq] = [hit.domain_name for hit in cevalue_arch]
            short_partial_cevalue_archs[seq] = [hit.domain_name if hit.domain_coverage > 0.5 else hit.domain_name+'_partial' for hit in cevalue_arch]

            short_coverage_archs[seq] = [hit.domain_name for hit in coverage_arch]
            short_partial_coverage_archs[seq] = [hit.domain_name if hit.domain_coverage > 0.5 else hit.domain_name+'_partial' for hit in coverage_arch]

            short_default_archs[seq] = [hit.domain_name for hit in default_arch]
            short_partial_default_archs[seq] = [hit.domain_name if hit.domain_coverage > 0.5 else hit.domain_name+'_partial' for hit in default_arch]
         
    if criterion == 'cevalue':
        # Create a df with the architectures as strings
        for seq, arch_list in short_partial_cevalue_archs.items():
            short_partial_cevalue_archs[seq] = str(arch_list)
        pfam_df = pd.DataFrame.from_dict(short_partial_cevalue_archs, orient='index').reset_index()
        pfam_df.columns = ['id', 'architecture']

    elif criterion == 'coverage':
        # Create a df with the architectures as strings
        for seq, arch_list in short_partial_coverage_archs.items():
            short_partial_coverage_archs[seq] = str(arch_list)
        pfam_df = pd.DataFrame.from_dict(short_partial_coverage_archs, orient='index').reset_index()
        pfam_df.columns = ['id', 'architecture']

    elif criterion == 'default':
        # Create a df with the architectures as strings
        for seq, arch_list in short_partial_default_archs.items():
            short_partial_default_archs[seq] = str(arch_list)
        pfam_df = pd.DataFrame.from_dict(short_partial_default_archs, orient='index').reset_index()
        pfam_df.columns = ['id', 'architecture']

    pfam_df['architecture_l'] = pfam_df['architecture'].apply(lambda x: ast.literal_eval(x))
    return pfam_df

def join_pfam_supfam_dfs(pfam_df, supfam_df):
    df = pd.merge(pfam_df, supfam_df, how='outer',left_on='id', right_on='id',suffixes=('_pfam','_supfam'))

    # Fill mising architectures with ''
    df['architecture_pfam'] = df['architecture_pfam'].fillna('')
    df['architecture_supfam'] = df['architecture_supfam'].fillna('')

    # Fill missing architectures with []
    df['architecture_l_pfam'] = df['architecture_l_pfam'].apply(lambda x: x if isinstance(x, list) else [])
    df['architecture_l_supfam'] = df['architecture_l_supfam'].apply(lambda x: x if isinstance(x, list) else [])

    # Create a column with # domains
    df['n_doms_pfam'] = df['architecture_l_pfam'].apply(lambda x: len(x))
    df['n_doms_supfam'] = df['architecture_l_supfam'].apply(lambda x: len(x))

    # Create a column with information whether it contains at least one partial domain
    df['contains_pfam_partial'] = df['architecture_pfam'].apply(lambda x: '_partial' in x)
    df['contains_supfam_partial'] = df['architecture_supfam'].apply(lambda x: '_partial' in x)

    return df

def filter_partial_df(df):
    pfam_len_1 = df['n_doms_pfam'] == 1
    pfam_has_partial = df['contains_pfam_partial'] == True
    pfam_no_arch = df['architecture_pfam'] == ''
    supfam_len_1 = df['n_doms_supfam'] == 1
    supfam_has_partial = df['contains_supfam_partial'] == True
    supfam_no_arch = df['architecture_supfam'] == ''

    # one partial Pfam and one partial Supfam
    pfam_partial_supfam_partial = pfam_len_1 & pfam_has_partial & supfam_len_1 & supfam_has_partial

    # one partial Pfam and no Supfam
    pfam_partial_no_supfam = pfam_len_1 & pfam_has_partial & supfam_no_arch

    # one partial supfam and no pfam
    supfam_partial_no_pfam = supfam_len_1 & supfam_has_partial & pfam_no_arch

    return df[~(pfam_partial_supfam_partial | pfam_partial_no_supfam | supfam_partial_no_pfam)]