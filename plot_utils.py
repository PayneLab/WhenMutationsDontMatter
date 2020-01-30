import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

'''
@Param df1:
     Dataframe of p-values for the first cancer (x-axis). 
     The first column must include the trans gene name. 
@Param df1_name:
    String. Name of the x-axis. Include the first cancer name. 

@Param df2:
    Dataframe of p-values for the second cancer (y-axis).
    The first column must include the trans gene name. 
    
@Param df2_name:
    String. Name of the y-axis. Include the second cancer name.

@Param cat_df:
    Dataframe of Booleans. First column must be gene names. 
    Can include up to 5 additional columns. Each column represents 
    a certain pathway or group of genes specified with boolean 
    variables (True is in the group). Each group will become another 
    plot graphing the p-values of the genes in that group.
    
@Param save_file_name:
    String. Optional to save the figure. The name of the file to 
    save the figure to.

This function plots the p-values for two cancers where 
x = cancer1_pvals and y = cancer2_pvals.  

The wrap_ttest function will help with creating df1 and df2. 
Both dataframes must have the gene included in the first column 
(omics can be included) and the p-value. Wrap_ttest creates this 
needed dataframe when its parameter return_all = True. 

The cat_df dataframe alows for multiple plots of specific gene groups
to be included. The first plot will always be all genes followed by one 
plot for each column included in cat_df. It is optional to include cat_df.
The arrangement of the subplots follows in rows and columns: (1,2) for
two plots, (2,2) for four plots(2,3) for six plots. The first plot will
allways be the comprehensive p-values for all genes done in the t-test
comparisons for the two cancers.
'''

def binary_pval_plot(df1, df1_name, df2, df2_name, cat_df=None, save_file_name=None):
    # Step 1: Combine Dataframes
    combined = df1.merge(df2, on=df1.columns[0]) # merge 2 pval df
    combined = combined.replace(regex=True,to_replace='_proteomics', value='') # only gene names
    combined_df = combined.rename(columns={combined.columns[1]: df1_name+" p-values", 
                                              combined.columns[2]: df2_name+" p-values"}) # Rename for use with x and y-axis
    if cat_df is not None:
        combined_df = combined.merge(cat_df, left_on= combined.columns[0], right_on= cat_df.columns[0]) # merge pathways
    
    else: # Plots one plot if no cat_df provided
        plt.figure(figsize=(6, 6)) 
        all_pvals = sns.scatterplot(x=combined_df.columns[1], y=combined_df.columns[2], data=combined_df)
        all_pvals.set_title("Comprehensive "+df1_name+" and "+df2_name+ " P-Values")
        all_pvals.set_xscale('log')
        all_pvals.set_yscale('log')
        plt.xlim(1e-5, 1e0) # 0.00005 to 1
        plt.ylim(1e-5, 1e0)
        
        if save_file_name is not None:
            plt.savefig(save_file_name+'.png')
        plt.show()
        plt.clf()        
        plt.close()
        return 0

    # Step 2: Find number of plots needed (dimensions of array for subplots)
    paths = list(cat_df.columns[1:])
    total = len(paths) + 1
    two_plots = False
    if total == 2:
        m_row = 1
        m_col = 2
        two_plots = True
    elif total <= 4:
        m_row = 2
        m_col = 2 
    elif total <= 6:
        m_row = 2
        m_col = 3
    else:
        print("Two many columns in cat_df. Can only plot 5 pathways/groups. (6 total plots)")
        return 1
        
    
    #sns.set(font_scale = 1.1)
    fig, axes = plt.subplots(m_row, m_col, sharex=True, sharey=True) # share x -axis title
    
    #Step 3: Create Multiple Plots
    if two_plots == True:  ### needed because subplots dimensions are a 1D array, not 2D
        # First plot with all p-values
        plt.rcParams['figure.figsize']=(5, 5) #size of plot
        ax = sns.scatterplot(x=combined_df.columns[1], y=combined_df.columns[2], data=combined_df,  
            ax=axes[0])
        ax.set_title("Comprehensive "+df1_name+" and "+df2_name+ " P-Values")
        ax.set_xscale('log')
        ax.set_yscale('log')
        plt.xlim(1e-5, 1e0) # 0.00005 to 1
        plt.ylim(1e-5, 1e0)
        # 1 Catagory plot
        group_name = str(combined_df.columns[4])
        only_p = combined_df.loc[combined_df[group_name] == True]
        ax = sns.scatterplot(x=combined_df.columns[1], y=combined_df.columns[2], data=only_p, color='orange',
            ax=axes[1]).set_title(group_name) #axes only has one number (1D array)

    else:
        # First plot with all p-values
        plt.rcParams['figure.figsize']=(15, 10) #size of plot
        ax = sns.scatterplot(x=combined_df.columns[1], y=combined_df.columns[2], data=combined_df,  
            ax=axes[0,0])
        ax.set_title("Comprehensive "+df1_name+" and "+df2_name+ " P-Values")
        ax.set_xscale('log')
        ax.set_yscale('log')
        plt.xlim(1e-5, 1e0) # 0.00005 to 1
        plt.ylim(1e-5, 1e0)
        
        # Catagory plots
        i = 0
        j = 1
        for e in paths:
            only_p = combined_df.loc[combined_df[e] == True]
            ax = sns.scatterplot(x=combined_df.columns[1], y=combined_df.columns[2], data=only_p, color='orange',
                ax=axes[i, j]).set_title(e)
            # i and j used to set next plot in axes
            if i <= (m_row - 1) and j < (m_col - 1):
                j += 1
            elif j == (m_col - 1):
                i+=1
                j=0
    
    if save_file_name is not None:
        fig.savefig(save_file_name+'.png')
    plt.show()
    plt.clf()        
    plt.close()
    return 0
