import numpy as np
import pandas as pd
import sys

def shifted_geometric_mean(vals, s):
    n = len(vals)
    product = max(vals[0] + 1.0, s)
    for i in range(1, n):
        product *= max(vals[i] + 1.0, s)
    mean = (product)**(1.0 / n) - 1.0
    return mean

def table2(path, out):
    df = pd.read_table(path, delim_whitespace=True)


    def get_optimal(tab, sol):
        if sol['Formulation'] == '-TDPLP':
            return sol['lp-relaxation']
        else:
            return tab[(tab['Instance'] == sol['Instance']) & (tab['Formulation'] == '-TDPLP')]['lp-relaxation'].values[0]
        
    df['optimal'] = df.apply(lambda x : get_optimal(df, x), axis=1)
    df['gap'] = (df['optimal'] - df['lp-relaxation']) / df['lp-relaxation']
    agg_table = df.groupby(['Formulation', 'Graph_type', 'n', 'r']).agg({
    'runtime' : (lambda x : shifted_geometric_mean(x.to_list(), 1)),
    'gap' : (lambda x : shifted_geometric_mean(x.to_list(), 1)),
    'num_cuts' : (lambda x : shifted_geometric_mean(x.to_list(), 1)),
    'separation_time' : (lambda x : shifted_geometric_mean(x.to_list(), 1))
    }).round(2)
    agg_table.to_csv(out, sep=',', na_rep='-')


def table3(path, out):

    df = pd.read_table(path, delim_whitespace=True)
    flow_ip_table = df[df['Formulation'] == '-FLOW'].groupby(['Graph_type', 'n', 'm', 'r']).sum()
    flow_table = df[df['Formulation'].str.strip() == '-FLOWLP'].groupby(['Graph_type', 'n', 'm', 'r']).sum()#.reindex(['Graph_type', 'n', 'm', 'r'], axis=1)
    flow_dp_table = df[df['Formulation'].str.strip() == '-FLOWLP+DP'].groupby(['Graph_type', 'n', 'm', 'r']).sum()
    #print(flow_table)
    
    table3 = pd.DataFrame(index = df.groupby(['Graph_type', 'n', 'm', 'r']).sum().index) 
    table3['Optimal Solution'] = flow_ip_table['lp-relaxation']
    table3['Flow LP Relaxation'] = flow_table['lp-relaxation']
    table3['Flow LP Runtime'] = flow_table['runtime']
    table3['Flow+DP LP Relaxation'] = flow_dp_table['lp-relaxation']
    table3['Flow+DP LP Runtime'] = flow_dp_table['runtime']
    table3['Flow+DP Cuts'] = flow_dp_table['num_cuts']
    table3['Gap Improvement'] = ( table3['Flow+DP LP Relaxation'] - table3['Flow LP Relaxation']  ) * 100.0 / ( table3['Flow+DP LP Relaxation'] )

    table3.to_csv(out, sep=',', na_rep='-')
    
def table4(path, out):

    df = pd.read_table(path, delim_whitespace=True)

    table4 = pd.DataFrame(index=df.groupby(['n', 'm', 'r']).sum().index)
    tcf_df = df[df['Formulation'].str.strip() == '-TCFIP'].groupby(['n', 'm', 'r']).sum()
    block_df = df[df['Formulation'].str.strip() == '-BLOCKIP'].groupby(['n', 'm', 'r']).sum()

    table4['Optimal ObjVal'] = block_df['best_obj']
    table4['TCF Runtime'] = tcf_df['runtime']
    table4['TCF Cuts'] = tcf_df['num_lazy_cuts']
    table4['Block Runtime'] = block_df['runtime']
    table4['Block Cuts'] = block_df['num_lazy_cuts']

    table4.to_csv(out, sep=',', na_rep='-')


def table5(path, out):
    
    df = pd.read_table(path, delim_whitespace=True)
    flow_df = df[df['prob_cut'] == 0].groupby(['n']).agg({'runtime':(lambda x : shifted_geometric_mean(x.to_list(), 1))})

    prim_df = df[(df['cut_type'].str.strip() == '-prim') & (df['r_pct_trees_cut'] == 1)].groupby(['prob_cut', 'n']).agg({'runtime':(lambda x : shifted_geometric_mean(x.to_list(), 1))})
    knap_125_df = df[(df['cut_type'].str.strip() == '-knap') & (df['r_pct_trees_cut'] == 1.25)].groupby(['prob_cut', 'n']).agg({'runtime':(lambda x : shifted_geometric_mean(x.to_list(), 1))})
    knap_150_df = df[(df['cut_type'].str.strip() == '-knap') & (df['r_pct_trees_cut'] == 1.5)].groupby(['prob_cut', 'n']).agg({'runtime':(lambda x : shifted_geometric_mean(x.to_list(), 1))})
    print(flow_df)
    
    table5 = pd.DataFrame(index=df[df['prob_cut'] > 0].groupby(['prob_cut', 'n']).sum().index)
    ns = table5.index.get_level_values('n')
    table5['FLOW unscaled'] = flow_df.loc[ns].values
    table5['r unscaled'] = prim_df['runtime']
    table5['r speed-up'] =  table5['FLOW unscaled'] / prim_df['runtime'] 
    table5['1.25r unscaled'] = knap_125_df['runtime']
    table5['1.25r speed-up'] = table5['FLOW unscaled'] / knap_125_df['runtime']
    table5['1.5r unscaled'] = knap_150_df['runtime'] 
    table5['1.5r speed-up'] = table5['FLOW unscaled'] /  knap_150_df['runtime']
    
    table5.to_csv(out, sep=',', na_rep='-')
    

def table6(path, out):

    df = pd.read_table(path, delim_whitespace=True)

    prim_df = df[(df['cut_type'].str.strip() == '-prim') & (df['r_pct_trees_cut'] == 1)].groupby(['prob_cut', 'n', 'm']).mean()
    knap_125_df = df[(df['cut_type'].str.strip() == '-knap') & (df['r_pct_trees_cut'] == 1.25)].groupby(['prob_cut', 'n', 'm']).mean()
    knap_150_df = df[(df['cut_type'].str.strip() == '-knap') & (df['r_pct_trees_cut'] == 1.5)].groupby(['prob_cut', 'n', 'm']).mean()
    table6 = pd.DataFrame(index=df[df['prob_cut'] > 0].groupby(['prob_cut', 'n', 'm']).sum().index)
    table6['r cuts'] = prim_df['num_cuts']
    table6['r sep. time'] = prim_df['separation_time']
    table6['r nodes'] = prim_df['b&b_nodes']
    table6['1.25r cuts'] = knap_125_df['num_cuts']
    table6['1.25r sep. time'] = knap_125_df['separation_time']
    table6['1.25r nodes'] = knap_125_df['b&b_nodes']
    table6['1.5r cuts'] = knap_150_df['num_cuts']
    table6['1.5r sep. time'] = knap_150_df['separation_time']
    table6['1.5r nodes'] = knap_150_df['b&b_nodes']
    table6.to_csv(out, sep=',', na_rep='-')

def table7(path, out):
    
    df = pd.read_table(path, delim_whitespace=True)

    tri_df = df[df['Formulation'].str.strip()=='-TRI'].groupby(['Graph_type', 'n']).agg({'runtime':(lambda x : shifted_geometric_mean(x.to_list(), 1))})
    flow_df = tri_df = df[df['Formulation'].str.strip()=='-FLOW'].groupby(['Graph_type', 'n']).agg({'runtime':(lambda x : shifted_geometric_mean(x.to_list(), 1))})
    path_df = df[df['Formulation'].str.strip()=='-PATH'].groupby(['Graph_type', 'n']).agg({'runtime':(lambda x : shifted_geometric_mean(x.to_list(), 1))})
    flowcut_df = df[df['Formulation'].str.strip()=='-FLOW+'].groupby(['Graph_type', 'n']).agg({'runtime':(lambda x : shifted_geometric_mean(x.to_list(), 1))})
    table7 = pd.DataFrame(index=df.groupby(['Graph_type', 'n']).sum().index)
    table7['TRI Runtime'] = tri_df['runtime']
    table7['FLOW Runtime'] = flow_df['runtime']
    table7['PATH Runtime'] = path_df['runtime']
    table7['FLOW+ Runtime'] = flowcut_df['runtime']
    table7['TRI Scaled'] = tri_df['runtime'] / table7[['TRI Runtime', 'FLOW Runtime', 'PATH Runtime', 'FLOW+ Runtime']].min(axis=1)
    table7['FLOW Scaled'] = flow_df['runtime'] / table7[['TRI Runtime', 'FLOW Runtime', 'PATH Runtime', 'FLOW+ Runtime']].min(axis=1)
    table7['PATH Scaled'] = path_df['runtime'] / table7[['TRI Runtime', 'FLOW Runtime', 'PATH Runtime', 'FLOW+ Runtime']].min(axis=1)
    table7['FLOW+ Scaled'] = flowcut_df['runtime'] / table7[['TRI Runtime', 'FLOW Runtime', 'PATH Runtime', 'FLOW+ Runtime']].min(axis=1)
    table7[['TRI Runtime', 'TRI Scaled', 'FLOW Runtime', 'FLOW Scaled', 'PATH Runtime', 'PATH Scaled', 'FLOW+ Runtime', 'FLOW+ Scaled' ]].to_csv(out, sep=',', na_rep='-')

def table8(path, out):
    

    df = pd.read_table(path, delim_whitespace=True)
    df['mn-ratio'] = df['m']/df['n']
    df['mn-ratio'] = df['mn-ratio'].round()

    tri_df = df[df['Formulation'].str.strip()=='-TRI'].groupby(['Graph_type', 'mn-ratio']).agg({'runtime':(lambda x : shifted_geometric_mean(x.to_list(), 1))})
    flow_df = tri_df = df[df['Formulation'].str.strip()=='-FLOW'].groupby(['Graph_type', 'mn-ratio']).agg({'runtime':(lambda x : shifted_geometric_mean(x.to_list(), 1))})
    path_df = df[df['Formulation'].str.strip()=='-PATH'].groupby(['Graph_type', 'mn-ratio']).agg({'runtime':(lambda x : shifted_geometric_mean(x.to_list(), 1))})
    flowcut_df = df[df['Formulation'].str.strip()=='-FLOW+'].groupby(['Graph_type', 'mn-ratio']).agg({'runtime':(lambda x : shifted_geometric_mean(x.to_list(), 1))})
    
    table8 = pd.DataFrame(index=df.groupby(['Graph_type', 'mn-ratio']).sum().index)

    table8['TRI Runtime'] = tri_df['runtime']
    table8['FLOW Runtime'] = flow_df['runtime']
    table8['PATH Runtime'] = path_df['runtime']
    table8['FLOW+ Runtime'] = flowcut_df['runtime']
    table8['TRI Scaled'] = tri_df['runtime'] / table8[['TRI Runtime', 'FLOW Runtime', 'PATH Runtime', 'FLOW+ Runtime']].min(axis=1)
    table8['FLOW Scaled'] = flow_df['runtime'] / table8[['TRI Runtime', 'FLOW Runtime', 'PATH Runtime', 'FLOW+ Runtime']].min(axis=1)
    table8['PATH Scaled'] = path_df['runtime'] / table8[['TRI Runtime', 'FLOW Runtime', 'PATH Runtime', 'FLOW+ Runtime']].min(axis=1)
    table8['FLOW+ Scaled'] = flowcut_df['runtime'] / table8[['TRI Runtime', 'FLOW Runtime', 'PATH Runtime', 'FLOW+ Runtime']].min(axis=1)
    table8[['TRI Runtime', 'TRI Scaled', 'FLOW Runtime', 'FLOW Scaled', 'PATH Runtime', 'PATH Scaled', 'FLOW+ Runtime', 'FLOW+ Scaled' ]].to_csv(out, sep=',', na_rep='-')

def table9(path, out):
    
    df = pd.read_table(path, delim_whitespace=True)
    
    table9 = pd.DataFrame(index=df.groupby(['Graph_type', 'n', 'm']).sum().index)
    tri_df = df[df['Formulation'].str.strip()=='-TRI'].groupby(['Graph_type','n', 'm']).agg({'runtime':(lambda x : shifted_geometric_mean(x.to_list(), 1))})
    flow_df = tri_df = df[df['Formulation'].str.strip()=='-FLOW'].groupby(['Graph_type','n', 'm']).agg({'runtime':(lambda x : shifted_geometric_mean(x.to_list(), 1))})
    path_df = df[df['Formulation'].str.strip()=='-PATH'].groupby(['Graph_type', 'n', 'm']).agg({'runtime':(lambda x : shifted_geometric_mean(x.to_list(), 1))})
    flowcut_df = df[df['Formulation'].str.strip()=='-FLOW+'].groupby(['Graph_type', 'n', 'm']).agg({'runtime':(lambda x : shifted_geometric_mean(x.to_list(), 1))})
    

    table9['TRI Runtime'] = tri_df['runtime']
    table9['FLOW Runtime'] = flow_df['runtime']
    table9['PATH Runtime'] = path_df['runtime']
    table9['FLOW+ Runtime'] = flowcut_df['runtime']
    table9['TRI Scaled'] = tri_df['runtime'] / table9[['TRI Runtime', 'FLOW Runtime', 'PATH Runtime', 'FLOW+ Runtime']].min(axis=1)
    table9['FLOW Scaled'] = flow_df['runtime'] / table9[['TRI Runtime', 'FLOW Runtime', 'PATH Runtime', 'FLOW+ Runtime']].min(axis=1)
    table9['PATH Scaled'] = path_df['runtime'] / table9[['TRI Runtime', 'FLOW Runtime', 'PATH Runtime', 'FLOW+ Runtime']].min(axis=1)
    table9['FLOW+ Scaled'] = flowcut_df['runtime'] / table9[['TRI Runtime', 'FLOW Runtime', 'PATH Runtime', 'FLOW+ Runtime']].min(axis=1)

    table9[['TRI Runtime', 'TRI Scaled', 'FLOW Runtime', 'FLOW Scaled', 'PATH Runtime', 'PATH Scaled', 'FLOW+ Runtime', 'FLOW+ Scaled' ]].to_csv(out, sep=',', na_rep='-')



def table10(path, out):
    df = pd.read_table(path, delim_whitespace=True)
    
    table10 = pd.DataFrame(index=df.groupby(['Graph_type', 'n', 'm', 'r']).sum().index)
    tri_df = df[df['Formulation'].str.strip()=='-TRI'].groupby(['Graph_type','n', 'm', 'r']).mean()
    flow_df = df[df['Formulation'].str.strip()=='-FLOW'].groupby(['Graph_type','n', 'm', 'r']).mean()
    path_df = df[df['Formulation'].str.strip()=='-PATH'].groupby(['Graph_type', 'n', 'm', 'r']).mean()
    flowcut_df = df[df['Formulation'].str.strip()=='-FLOW+'].groupby(['Graph_type', 'n', 'm', 'r']).mean()

    table10['TRI Bound'] = tri_df['best_obj']
    table10['TRI Runtime'] = tri_df['runtime']
    table10['TRI Gap'] = tri_df['mip_gap']
    table10['FLOW Bound'] = flow_df['best_obj']
    table10['FLOW Runtime'] = flow_df['runtime']
    table10['FLOW Gap'] = flow_df['mip_gap']
    table10['PATH Bound'] = path_df['best_obj']
    table10['PATH Runtime'] = path_df['runtime']
    table10['PATH Gap'] = path_df['mip_gap']
    table10['FLOW+ Bound'] = flowcut_df['best_obj']
    table10['FLOW+ Runtime'] = flowcut_df['runtime']
    table10['FLOW+ Gap'] = flowcut_df['mip_gap']

    table10.to_csv(out, sep=',', na_rep='-')



if __name__ == '__main__':
    table_num = sys.argv[1]
    table_path = sys.argv[2]
    out_path = sys.argv[3]
    if int(table_num) == 2:
        table2(table_path, out_path)

    elif int(table_num) == 3:
        table3(table_path, out_path)

    elif int(table_num) == 4:
        table4(table_path, out_path)

    elif int(table_num) == 5:
        table5(table_path, out_path)

    elif int(table_num) == 6:
        table6(table_path, out_path)

    elif int(table_num) == 7:
        table7(table_path, out_path)

    elif int(table_num) == 8:
        table8(table_path, out_path)
    
    elif int(table_num) == 9:
        table9(table_path, out_path)

    elif int(table_num) == 10:
        table10(table_path, out_path)

    else:
        exit