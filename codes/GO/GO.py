# --- Load Modules ---
import numpy as np
import pandas as pd
import matplotlib as mpl
mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['ps.fonttype'] = 42
mpl.rcParams['font.sans-serif'] = 'Arial'
import matplotlib.pyplot as plt
import seaborn as sns
from pprint import pprint


# --- Load Genes ---
nanog_256_genes = pd.read_csv('../../Supplemental_Information/Table_S5_Nanog_target_genes/Nanog_pre_MBT_target.csv')
nanog_1k_new_genes = pd.read_csv('../../Supplemental_Information/Table_S5_Nanog_target_genes/Nanog_MBT_target.csv')
nanog_dome_new_genes = pd.read_csv('../../Supplemental_Information/Table_S5_Nanog_target_genes/Nanog_post_MBT_target.csv')
nanog_256_genes_in_10k = set(nanog_256_genes['Gene'].values)
nanog_1k_new_genes_in_10k = set(nanog_1k_new_genes['Gene'].values)
nanog_dome_new_genes_in_10k = set(nanog_dome_new_genes['Gene'].values)
nanog_256_genes_in_10kb_cluster, nanog_256_genes_in_10kb_remove_cluster = set(), set()
for index, row in nanog_256_genes.iterrows():
    if row['Nanog binding cluster (NBC) target']:
        nanog_256_genes_in_10kb_cluster.add(row['Gene'])
    else:
        nanog_256_genes_in_10kb_remove_cluster.add(row['Gene'])

# --- perform g:GOSt ---
from gprofiler import GProfiler
query={
        'pre-MBT': list(nanog_256_genes_in_10k),
        'MBT': list(nanog_1k_new_genes_in_10k),
        'post-MBT': list(nanog_dome_new_genes_in_10k),
        'NBC target pre-MBT': list(nanog_256_genes_in_10kb_cluster),
        'other target pre-MBT': list(nanog_256_genes_in_10kb_remove_cluster),
    }

gp = GProfiler(return_dataframe=True)
GO_output = gp.profile(
    organism='drerio',
    sources=['GO:CC','GO:MF','GO:BP'],
    query=query,
    all_results=True,
    no_evidences=False,
)

GO_output['enrichment'] = GO_output['intersection_size'] / GO_output['query_size'] / (GO_output['term_size'] / GO_output['effective_domain_size'])
GO_output['-log10 P-value'] = -np.log10(GO_output['p_value'])


# --- merge ---
def linked(gene_list1, gene_list2, threshold=0.7):
    return True if len(set(gene_list1) & set(gene_list2)) / len(
        set(gene_list1) | set(gene_list2)) > threshold else False


def term_merge(GO_results,
               query_name,
               index,
               selected_term_index,
               threshold={
                   'merge': 0.7,
                   'term_size': (15, 500),
                   'p_value': 0.001
               }):
    current_selected_term_index = set([index])
    current_index = [index]
    next_index = []
    _round = 0
    while current_index:
        _round += 1
        #print(f'round: {_round}, term size: {len(current_selected_term_index)}, next index: {current_index}')
        for i in GO_results.index:
            if i in selected_term_index or i in current_selected_term_index or GO_results.loc[
                    i, 'query'] != query_name or GO_results.loc[
                        i, 'term_size'] < threshold['term_size'][
                            0] or GO_results.loc[i, 'term_size'] > threshold[
                                'term_size'][1]:
                continue
            if linked(GO_results.loc[index, 'intersections'],
                      GO_results.loc[i, 'intersections'],
                      threshold=threshold['merge']):
                current_selected_term_index.add(i)
                next_index.append(i)
        current_index, next_index = next_index[:], []
    return current_selected_term_index


def single_query_cluster_gen(GO_results,
                             query_name,
                             threshold={
                                 'merge': 0.7,
                                 'term_size': (15, 500),
                                 'p_value': 0.001
                             }):
    '''
    return: gen(DataFrame)
    '''
    filteres = [
        True if x < threshold['p_value'] and y == query_name
        and threshold['term_size'][0] <= z <= threshold['term_size'][1] else
        False for x, y, z in zip(
            GO_results['p_value'],
            GO_results['query'],
            GO_results['term_size'],
        )
    ]
    filtered_terms = GO_results.loc[filteres, :].sort_values(by='enrichment',
                                                             ascending=False)
    cluster = pd.DataFrame(columns=GO_results.columns)
    selected_term_index = set()
    current_selected_term_index = set()
    for index in filtered_terms.index:
        if index in selected_term_index:
            continue
        #print(index)
        current_selected_term_index = term_merge(GO_results, query_name, index,
                                                 selected_term_index,
                                                 threshold)
        #print(current_selected_term_index)
        yield GO_results.iloc[list(current_selected_term_index
                                   ), :].sort_values(by='enrichment',
                                                     ascending=False)
        selected_term_index = selected_term_index | current_selected_term_index
        current_selected_term_index = set()


def query_cluster(GO_results, query, natives):
    for query_name in query.keys():
        try:
            index, row = next(GO_results.loc[[
                True if x in natives and y == query_name else False
                for x, y in zip(GO_results['native'], GO_results['query'])
            ], ['-log10 P-value', 'enrichment']].sort_values(
                by='enrichment', ascending=False).head(1).iterrows())
            yield query_name, row['-log10 P-value'], row['enrichment']
        except StopIteration:
            yield query_name, np.nan, np.nan


def gen_top_K_cluster(GO_results,
                      query,
                      K=3,
                      threshold={
                          'merge': 0.7,
                          'term_size': (15, 500),
                          'p_value': 0.001
                      }):
    import csv
    import numpy as np
    from collections import defaultdict
    query_name_list = list(query.keys())
    GO_cluster = {
        query_name: {
            'term_id': [],
            '-log10 P-value': [],
            'enrichment score': [],
            'name': [],
        }
        for query_name in query_name_list
    }
    for index, query_name in enumerate(query_name_list):
        with open(f'GO_{query_name}_top{K}_cluster.csv', 'w') as fhd:
            f_csv = csv.writer(fhd)
            f_csv.writerow([
                'source', 'native', 'name', '-log10 p-value',
                'enrichment score', 'intersection genes', 'term size',
                'quer size', 'intersection size', 'effective domain size'
            ])
            go_results_cluster_gen = single_query_cluster_gen(
                GO_output, query_name)
            for i in range(K):
                #print(query_name,i)
                try:
                    cluster = next(go_results_cluster_gen)
                except StopIteration:
                    break
                if cluster.shape[0] == 0:
                    break
                f_csv.writerow(['#', '-----', f'cluster {i+1}'] +
                               ['-----'] * 7)
                go_term_name = ''
                count = 0
                for _, row in cluster.iterrows():
                    count += 1
                    f_csv.writerow([
                        row['source'], row['native'], row['name'],
                        -np.log10(row['p_value']), row['enrichment'],
                        ','.join(row['intersections']), row['term_size'],
                        row['query_size'], row['intersection_size'],
                        row['effective_domain_size']
                    ])
                    if count == 1:
                        go_term_name = row['name']
                term_ids = list(cluster['native'])
                for row in query_cluster(GO_output, query, term_ids):
                    GO_cluster[row[0]]['term_id'].append(set(term_ids))
                    GO_cluster[row[0]]['-log10 P-value'].append(row[1])
                    GO_cluster[row[0]]['enrichment score'].append(row[2])
                    GO_cluster[row[0]]['name'].append(go_term_name)
    #print(GO_cluster)
    merge_cluster = []
    indexes_to_be_merged = set(
        range(len(GO_cluster[query_name_list[0]]['term_id'])))
    ###BFS
    # get edges
    edges = defaultdict(list)
    for i in range(len(indexes_to_be_merged)):
        for j in range(i + 1, len(indexes_to_be_merged)):
            if len(GO_cluster[query_name_list[0]]['term_id'][i]
                   & GO_cluster[query_name_list[0]]['term_id'][j]) > 0:
                edges[i].append(j)
                edges[j].append(i)
    #print(edges)
    # merge
    merged_indexex = set()
    while len(indexes_to_be_merged - merged_indexex) > 0:
        merge_cluster.append([])
        current_index = [list(indexes_to_be_merged - merged_indexex)[0]]
        merged_indexex.add(current_index[0])
        while current_index:
            merge_cluster[-1].extend(current_index)
            next_index = []
            for index in current_index:
                for next_node in edges[index]:
                    if next_node in merged_indexex:
                        continue
                    next_index.append(next_node)
                    merged_indexex.add(next_node)
            current_index = next_index[:]
    #print(merge_cluster)
    for mc in merge_cluster:
        for query_name in query_name_list:
            target_index = mc[np.array([
                GO_cluster[query_name]['enrichment score'][i] for i in mc
            ]).argmax()]
            yield GO_cluster[query_name]['name'][mc[0]], GO_cluster[
                query_name]['-log10 P-value'][target_index], GO_cluster[
                    query_name]['enrichment score'][target_index], query_name


GO_cluster = pd.DataFrame(
    gen_top_K_cluster(GO_output,
                      query,
                      K=3,
                      threshold={
                          'merge': 0.7,
                          'term_size': (15, 500),
                          'p_value': 0.001
                      }),
    columns=['name', '-log10 P-value', 'enrichment score', 'sample'])


GO_cluster.sort_values(by=['enrichment score'],ascending=False,inplace=True)
GO_cluster['enrichment score'] = [x if x <=50 else 50 for x in GO_cluster['enrichment score']]


# --- plot ---
with sns.axes_style('white', rc={
        'xtick.bottom': True,
        'ytick.left': True
}), sns.plotting_context('paper',
                         rc={
                             'axes.titlesize': 14,
                             'axes.labelsize': 12,
                             'xtick.labelsize': 10,
                             'ytick.labelsize': 10,
                             'legend.fontsize': 10
                         }):

    fig, ax = plt.subplots(figsize=(6.4 * 1.5, 4.8))
    sns.scatterplot(
        x='sample',
        y='name',
        size='enrichment score',
        hue='-log10 P-value',
        data=GO_cluster,
        palette='YlGnBu',
        sizes=(0, 700),
        edgecolor='k',
        ax=ax,
    )
    ax.legend(bbox_to_anchor=(1, 0), loc='lower left')
    #ax.set_xticks(range(5))
    #ax.set_xticklabels(
    #    ['pre-MBT', 'MBT', 'post-MBT', 'NBC target pre-MBT', 'othr target pre-MBT'])
    ax.set_xlabel('Nanog target gene class')
    ax.set_ylabel('GO term name')
    fig.tight_layout()
    fig.savefig('Fig3B_Scatterplot_GO_results.pdf', transparent=True)


gradient = np.linspace(0, 1, 256)
gradient = np.vstack((gradient, gradient))
with sns.axes_style('white', rc={
        'xtick.bottom': True,
        'ytick.left': True
}), sns.plotting_context('paper',
                         rc={
                             'axes.titlesize': 10,
                             'axes.labelsize': 8,
                             'xtick.labelsize': 6,
                             'ytick.labelsize': 6,
                             'legend.fontsize': 6
                         }):
    fig, ax = plt.subplots(figsize=(6.4 * .5,4.8 * .3))
    ax.imshow(gradient, aspect='auto', cmap='YlGnBu')
    ax.set_xticks(np.linspace(0,256,6))
    ax.set_xticklabels(np.arange(0,18,3))
    ax.set_yticks([])
    fig.savefig('Fig3B_colorbar.pdf',transparent=True)
    
