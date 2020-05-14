import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import math as math
import scipy.stats
import re
import sys
import statsmodels.stats.multitest
from bokeh.palettes import RdBu
from bokeh.models import LinearColorMapper, ColumnDataSource, ColorBar
from bokeh.models.ranges import FactorRange
from bokeh.plotting import figure, show
from bokeh.io import output_notebook
from bokeh.io import export_png
from bokeh.io import export_svgs


'''
@Param df: Dataframe. Contains column with x-axis categorical variables, y-axis categorical variables,
and columns for circle size and color gradient. 
@Param circle_var. String. Name of column for numeric data to base circle size off of. Designed for p_values. 
@Param color_var. String. Name of column of numeric data to base color gradient off of. Designed for correlation value or difference in medians
@Param x_axis String. Name of column for x-axis categorical labels
@Param y_axis String. Name of column for y-axis categorical labels
@Param x_axis_lab. String. Default is no label. 
@Param y_axis_lab. String. Default is no label. 
@Param plot_width. Default is 10000
@Param plot_heigh. Default is 650.
@Param save_png. Default is to not save a png file. To save file, set save_png to the name (string) you would like to save the file as. *NOTE* inorder to use this function you must download phantomjs onto computer by using conda install phantomjs


This function creates a bokeh map that is heat map with extra variable of size of the circles. 

'''
def plotCircleHeatMap ( df, circle_var, color_var, x_axis, y_axis,plot_width= 1000, plot_height = 650, x_axis_lab = "no_label", y_axis_lab = "no_label", show_plot = True, save_png = "plot.png"):
  
    # circle_var designed for pvalues. Normalized by taking log 10 of values and multiplying by 5 
    #added a new column to make the plot size
    
    df["size2"] = df[circle_var].apply(lambda x: -1*(np.log(x)))
    df['size'] = (df["size2"])*1.5
    #find values to set color bar min/ max as 
    maxval = df[color_var].max()
    minval = df[color_var].min()
    if maxval > abs(minval):
        minval = maxval * -1 
    if maxval < abs(minval):
        maxval = minval * -1
    colors = list((RdBu[9]))
    exp_cmap = LinearColorMapper(palette=colors, low = minval, high = maxval)
    p = figure(x_range = FactorRange(), y_range = FactorRange(), plot_width= plot_width, 
               plot_height=plot_height, 
               toolbar_location=None, tools="hover")

    p.scatter(x_axis,y_axis,source=df, fill_alpha=1,  line_width=0, size="size", 
              fill_color={"field":color_var, "transform":exp_cmap})

    p.x_range.factors = sorted(df[x_axis].unique().tolist())
    p.y_range.factors = sorted(df[y_axis].unique().tolist(), reverse = True)
    p.xaxis.major_label_orientation = math.pi/2
    
    if (x_axis_lab != "no_label" ):
        p.xaxis.axis_label = x_axis_lab
    if (x_axis_lab != "no_label" ):   
        p.yaxis.axis_label = y_axis_lab

    bar = ColorBar(color_mapper=exp_cmap, location=(0,0))
    p.add_layout(bar, "right")
    if show_plot:
        output_notebook()
        show(p)
      
    if save_png != "plot.png":
        export_png(p, filename= save_png)

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
        combined_df = combined_df.merge(cat_df, left_on= combined_df.columns[0], right_on= cat_df.columns[0]) # merge pathways

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

    fig, axes = plt.subplots(m_row, m_col, sharex=True, sharey=True) # share x -axis title
    plt.rcParams['figure.figsize']=(15, 10) #size of plot
    sns.set(font_scale = 1.1)

    #Step 3: Create Multiple Plots
    if two_plots == True:  ### needed because subplots dimensions are a 1D array, not 2D
        # First plot with all p-values
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


import pandas as pd
import scipy.stats
import statsmodels.stats.multitest

'''
@Param df: Dataframe. Contains numeric values (such as proteomics) for linear regression
@Param label_column: String. Name of column that will be your x axis and will be compared to all values in df unless otherwise specified.
@Param alpha: significant level
@Param comparison_columns: columns that will be looped through and used as y axis for linear regression.
All other columns beside label column unless specified here.
@Param correction_method: String. Specifies method of adjustment for multiple testing. See -
https://www.statsmodels.org/stable/generated/statsmodels.stats.multitest.multipletests.html
    - for documentation and available methods.

This function will return a data frame will all significant linear regressions. The data frame includes the comparison, slope, R-squared, and P-value.
'''
import pandas as pd
import numpy as np
import scipy.stats
import statsmodels.stats.multitest
import re
import sys
'''
@Param df: Dataframe.Each column is a different gene/ comparison. Rows contains numeric values (such as proteomics) for correlation test
@Param label_column: String. Name of column that will be your x axis and will be compared to all values in df unless otherwise specified.
@Param alpha: significant level
@Param comparison_columns: columns that will be looped through and used as y axis for correlation test.
All other columns beside label column unless specified here.
@Param correction_method: String. Specifies method of adjustment for multiple testing. See -
https://www.statsmodels.org/stable/generated/statsmodels.stats.multitest.multipletests.html
    - for documentation and available methods.

This function will return a data frame with the columns comparison, the correlation coefficient, and the p value.
'''
def wrap_pearson_corr(df,label_column, alpha=.05,comparison_columns=None,correction_method='bonferroni',return_all = True):


    #df = df.dropna(axis=1, how="all")

    '''If no comparison columns specified, use all columns except the specified labed column'''
    if not comparison_columns:
        comparison_columns = list(df.columns)
        comparison_columns.remove(label_column)
    '''Store comparisons,p-values, correlation in their own array'''
    comparisons = []
    pvals = []
    correlation=[]


    '''Format results in a pandas dataframe'''
    newdf = pd.DataFrame(columns=['Comparison','Correlation','P_value'])
    for gene in comparison_columns:
        #create subset df with interacting gene/ gene (otherwise drop NaN drops everything)
        df_subset = df[[label_column,gene]]
        #do a linear regression to see if it's a meaningful association
        #dropna will remove rows with nan
        df_subset = df_subset.dropna(axis=0, how="any")
        count_row = df_subset.shape[0]
        if count_row > 20:
            x1 = df_subset[[label_column]].values
            y1 = df_subset[[gene]].values
            x1 = x1[:,0]
            y1 = y1[:,0]
            corr, pval = scipy.stats.pearsonr(x1,y1)

            comparisons.append(gene)
            pvals.append(pval)
            correlation.append(corr)


    '''Correct for multiple testing to determine if each comparison meets the new cutoff'''
    results = statsmodels.stats.multitest.multipletests(pvals=pvals, alpha=alpha, method=correction_method)
    reject = results[0]

    if return_all:
        for i in range(0,len(comparisons)):
            newdf = newdf.append({'Comparison': comparisons[i],"Correlation": correlation[i],'P_value': pvals[i]}, ignore_index=True)

    '''Else only add significant comparisons'''
    if (return_all == False):
            for i in range(0, len(reject)):
                if reject[i]:
                    newdf = newdf.append({'Comparison': comparisons[i],"Correlation": correlation[i],'P_value': pvals[i]}, ignore_index=True)

    '''Sort dataframe by ascending p-value'''
    newdf = newdf.sort_values(by='P_value', ascending=True)
    '''If results df is not empty, return it, else return None'''
    return newdf



'''
@Param df1: Dataframe. Contains numeric values (such as proteomics) for linear regression
@Param x_axis: String. Used as the label for the x-axis as well as the column name for the x-axis values.
@Param y_axis:String. Used as the label for the y-axis as well as the column name for the y-axis values.
@Param title: String. Used as title of figure
@Param ra_stats: Boolean. Defalt is False. If true it will print out the linear regression stats.
@Param show_plot: Boolean. Defalt is True. If true it will show the linear regression plot.
@Param save_file_name: String.File will not be saved unless a file name is specified. Plot is saved in current directory as png.

This fuction takes a dataframe with numeric values (such as proteomics) and performs a linear regression analysis between two user specified columns within the dataframe. The function will then create the linear regression graph and can print the graph to the screen and save the figure depending on user input.
'''


def plot_lin_regression(df1,x_axis, y_axis, title, ra_stats = False, show_plot = True, save_file_name = "file_name" ):
    #subset df to have just the x and y axis specified
    df1_subset = df1[[x_axis,y_axis]]
    #drop na values. Can't have in linear regression
    df1_subset = df1_subset.dropna(axis=0, how="any")

    x1 = df1_subset[[x_axis]].values
    y1 = df_gbm_subset[[y_axis]].values
    x1 = x1[:,0]
    y1 = y1[:,0]

    slope, intercept, r_value, p_value, std_err = scipy.stats.linregress(x1,y1)
    if ra_stats:
        print ('Slope of regression: %s\nR-squared: %s\nP-value: %s'%(slope, r_value**2, p_value))

    sns.set(style="darkgrid")
    plot = sns.regplot(x=x1, y=y1, data=df1)
    plot.set(xlabel=x_axis, ylabel=y_axis, title=title)
    if show_plot:
        plt.show()
        plt.clf()
        plt.close()

    if save_file_name != "file_name":
        plt.savefig(save_file_name+'.png')


'''
@Param dflist: List. A list of cancer dataframes that contain the mutation type frequnecy (dataframes that are returned from the get_genotype_all_var function)
@Param names_of_df: List. Names of the cancers (names MUST correlate to dflist)
@Param title: String. The title of the graph.
@Param save_to_path: String. The absoloute path to save the figure to. Defualts to saveing in current directory as step_1.png

This function takes a list of dataframes that contain the mutation type frequncy for cancers and plots them.
'''
def figure1_plot_mutations(dflist = None, names_of_df=None, title=None, save_to_path=None):
    number_of_df = len(dflist)

    allLabels = []
    for df in dflist:
        #get the labels for each and make a combined label that they'll all use
        mutation = df["Mutation"]
        labels = list(set(mutation))

        allLabels.append(labels)

    flat_list = [item for sublist in allLabels for item in sublist]
    all_labels = list(set(flat_list))
    all_labels.sort()
    allLabels = all_labels

#     For each df, add na to their labels if it doesn't exist in all_labels
    labels_for_each_df = []
    frequencies_for_each_df = []
    for df in dflist:
        mutation = df["Mutation"].tolist()
        mutationlist = list(set(mutation))
        mutationlist.sort()
        ordered_mut_list = []
        match = True
        mutPosition = 0

        for position in range(len(all_labels)):
            try:

                if mutationlist[mutPosition] == all_labels[position]:
                    ordered_mut_list.append(mutationlist[mutPosition])
                    mutPosition += 1

                else:
                    ordered_mut_list.append("na")

            except IndexError:
                ordered_mut_list.append("na")


        labels_for_each_df.append(ordered_mut_list)

        #get the freq of each mutation type
        freq = []
        for mutation_type in ordered_mut_list:
            freq.append(mutation.count(mutation_type))

        PercentFreq = [x*100 / sum(freq) for x in freq]
        frequencies_for_each_df.append(PercentFreq)


    #Now plot it using arrays
    width = 0.1
    x = np.arange(len(allLabels))
    a4_dims = (25, 13) #dimensions for bigger plot
    fig, ax = plt.subplots(figsize=a4_dims)
    for position in range(0, number_of_df):
        r = ax.bar(x+(width*position), frequencies_for_each_df[position], width,label=names_of_df[position], alpha=.5, linewidth=0)



    ax.set_ylabel('Percent Sample')
    ax.set_title(title)
    ax.set_xticks(x)
    ax.set_xticklabels(allLabels)
    ax.legend()



    plt.setp(ax.get_xticklabels(),rotation='vertical')
    if save_to_path == None:
        plt.savefig("step_1.png")
    else:
        plt.savefig(save_to_path)

    plt.show()
