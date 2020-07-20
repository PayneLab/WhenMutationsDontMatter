import pandas as pd
import seaborn as sns
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
from bokeh.io import output_notebook, export_png, export_svgs
from bokeh.layouts import row


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
@Param font_size. Int. Default is 12pt.
@Param save_png. Default is to not save a png file. To save file, set save_png to the name (string) you would like to save the file as. *NOTE* inorder to use this function you must download phantomjs onto computer by using conda install phantomjs


This function creates a bokeh map that is heat map with extra variable of size of the circles. 

'''
def plotCircleHeatMap ( df, circle_var, color_var, x_axis, y_axis, plot_width= 1000, plot_height = 650, font_size = 12, x_axis_lab = "no_label", y_axis_lab = "no_label", show_plot = True, save_png = "plot.png", legend_min = 1e-6, legend_max = 0.01, show_legend = True):
  
    # circle_var designed for pvalues. Normalized by taking log 10 of values and multiplying by 5 
    #added a new column to make the plot size
    
    df["size2"] = df[circle_var].apply(lambda x: -1*(np.log(x)))
    df['size'] = (df["size2"])*3
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

    #p.x_range.factors = sorted(df[x_axis].unique().tolist())
    p.x_range.factors = df[x_axis].unique().tolist()
    p.y_range.factors = sorted(df[y_axis].unique().tolist(), reverse = True)
    p.xaxis.major_label_orientation = math.pi/2
    
    # font size
    p.axis.major_label_text_font_size = str(font_size)+"pt"
    
    if (x_axis_lab != "no_label" ):
        p.xaxis.axis_label = x_axis_lab
    if (y_axis_lab != "no_label" ):   
        p.yaxis.axis_label = y_axis_lab

    bar = ColorBar(color_mapper=exp_cmap, location=(0,0))
    p.add_layout(bar, "right")
    
    # Create Circle Legend
    circle_legend = create_circle_legend(df, circle_var, color_var, legend_min, legend_max)
    
    if show_plot:  
        if show_legend:
            # Create Circle Legend
            circle_legend = create_circle_legend(df, circle_var, color_var, legend_min, legend_max)
            output_notebook()
            show(row(p, circle_legend))
        else:
            output_notebook()
            show(p)
      
    if save_png != "plot.png":
        export_png(row(p, circle_legend), filename= save_png)
             
        
'''
@Param df: Dataframe. Same as df passed to plotCircleHeatMap.
@Param legend_min: Float. Lowest p-value to include in the legend.
@Param legend_max: Float. Highest p-value to include in the legend.

Returns: df to be used in creating the circle legend. 
'''

def create_circle_legend_df(color_var, legend_min, legend_max):
    
    # Find middle pvals
    # Find difference between exponents of the min and max
    exp_legend_min = np.log10(legend_min)
    exp_legend_max = np.log10(legend_max)
    delta = exp_legend_min - exp_legend_max 
    # Split difference into quarters 
    num = 4
    val = delta / num
    # Find middle exponents
    exp2 = round(exp_legend_min - val)
    pval2 = 1*10**exp2
 
    
    # Foramat scientific notation pvals as strings for y_axis labels  
    max_str = "{:.1e}".format(legend_max, '.2f')

    pval_str_2 = "{:.1e}".format(pval2, '.2f')
    min_str = "{:.1e}".format(legend_min, '.2f')
    
    # max to min
    data = {'P_Value':  [legend_max, pval2, legend_min],
            'y_axis': [max_str, pval_str_2, min_str],
            'x_axis': ['', '', ''],
            color_var: [1.5, 1.5, 1.5]}

    fake_df = pd.DataFrame (data, columns = ['x_axis', 'y_axis', 'P_Value', color_var])
    
    fake_df["size2"] = fake_df['P_Value'].apply(lambda x: -1*(np.log(x)))
    fake_df['size'] = (fake_df["size2"])*3
    
    return fake_df

'''
@Param df: Dataframe. Same as df passed to plotCircleHeatMap.
@Param circle_var: Column Label. Same as passed to plotCircleHeatMap.
@Param color_var: Column Label. Same as passed to plotCircleHeatMap.
@Param x_axis: Column Label. Used on the x-axis.
@Param y_axis: Column Label. Used on the y-axis.
@Param legend_min: Float. Lowest p-value to include in the legend.
@Param legend_max: Float. Highest p-value to include in the legend.

Returns: df to be used in creating the circle legend. 
'''

def create_circle_legend(df, circle_var, color_var, legend_min, legend_max,
                         x_axis = 'x_axis', y_axis = 'y_axis', 
                         plot_height = 200, plot_width = 120):
    # Use the smallest pval
    if df[circle_var].min() < legend_min:
        legend_min = df[circle_var].min()
    circle_df = create_circle_legend_df(color_var, legend_min, legend_max)
    
    maxval = circle_df[color_var].max()
    minval = circle_df[color_var].min()
    if maxval > abs(minval):
        minval = maxval * -1 
    if maxval < abs(minval):
        maxval = minval * -1
    colors = list((RdBu[9]))
    exp_cmap = LinearColorMapper(palette=colors, low = minval, high = maxval)
    
    circle = figure(x_range = FactorRange(), y_range = FactorRange(), plot_width= plot_width, 
               plot_height=plot_height, toolbar_location=None, tools="hover")

    circle.scatter(x_axis, y_axis, source = circle_df, fill_alpha=1,  line_width=0, size="size", 
              fill_color={"field":color_var, "transform":exp_cmap})
    
    circle.x_range.factors = sorted(circle_df[x_axis].unique().tolist())
    circle.y_range.factors = circle_df[y_axis].unique().tolist() # plots in reverse order of circle_df (max to min)
    circle.xaxis.major_label_orientation = math.pi/2
    
    circle.xaxis.axis_label = 'FDR_P-Values'
    
    
    return circle
      


'''
@Param plot_df:
     Dataframe with at least 2 columns of either p-value or correlation 
     values for two cancers. Can include up to 3 additional pathway columns. 
     Each pathway column represents a certain pathway or group of genes 
     specified with boolean variables (True = is in the group). Each pathway  
     column will become another scatterplot. 
       NOTE: Column names will be used as titles for the pathway plots.
       
@Param val:
    String. "corr": Correlation values 
            "pval": p-values
            "median": difference in medians
            This parameter specifies the scale of the plot
    
@Param x_axis_col:
    String. Column label to plot on the x-axis.

@Param y_axis_col:
    String. Column label to plot on the y-axis.
    
@Param plot_pathways:
    List. Contains the column names containing pathway boolean values.
    
@Param hue_col:
    String. Column label to base hue color.

@Param hue_dict:
    Dictionary. A dictionary with values from the hue_col as keys 
    and chosen colors for the plot as values. 
    
@Param narrow:
    Float. Limit of the graph scale. Can be a positive or negative float. 

@Param save_file_name:
    String. Optional to save the figure. The name of the file to
    save the figure to.

Returns:
    Scatterplots.

The wrap_ttest function will help with creating plot_df.
Both dataframes must have the gene included in the first column
(omics can be included) and the p-value. Wrap_ttest creates this
needed dataframe when its parameter return_all = True.

The path_df dataframe alows for multiple plots of specific gene groups
to be included. The first plot will always be all genes followed by one
plot for each column included in path_df. 
'''


def binary_val_plot(plot_df, val, x_axis_col, y_axis_col, title, pathway_cols = None, hue_col=None, 
                    color_dict = None, narrow = 0, save_file_name=None):
    
    f = plt.figure(figsize=(12, 12))
    gs = f.add_gridspec(2, 2)
    sns.set(font_scale = 1.2)
        
    # Create Main plot
    all_pvals = f.add_subplot(gs[0, 0])
    if hue_col is not None:
        all_pvals = sns.scatterplot(x=x_axis_col, y=y_axis_col, data=plot_df,
                                   hue=hue_col, palette=color_dict) 
    else:
        all_pvals = sns.scatterplot(x=x_axis_col, y=y_axis_col, data=plot_df)
    all_pvals.set_title(title)    
    plt.xlabel(x_axis_col)
    plt.ylabel(y_axis_col)
    
    # Set scale of x and y axis
    if narrow != 0:
        n = np.abs(narrow)
        axes = {'x_low': -1*n, 'x_high': n, 'y_low': -1*n, 'y_high': n}
    else:
        if val == 'median':
            axes = {'x_low': -1.5, 'x_high': 1.5, 'y_low': -1.5, 'y_high': 1.5}
        elif val == 'corr':
            axes = {'x_low': -1, 'x_high': 1, 'y_low': -1, 'y_high': 1}
        elif val == 'pval':
            all_pvals.set_xscale('log')
            all_pvals.set_yscale('log')
            axes = {'x_low': 1e-5, 'x_high': 1e0, 'y_low': 1e-5, 'y_high': 1e0}

    plt.xlim(axes['x_low'], axes['x_high']) 
    plt.ylim(axes['y_low'], axes['y_high'])
    
    
    # Create subplots for pathways/groups in path_df
    if pathway_cols is not None:
        i = 0
        j = 1
        for col in pathway_cols:
            bool_serries = plot_df[col]
            group_df = plot_df[bool_serries] # keep genes in pathway
            path = f.add_subplot(gs[i, j])
            path = sns.scatterplot(x=x_axis_col, y=y_axis_col, data=group_df, color = 'orange')
            if val == 'pval':
                path.set_xscale('log')
                path.set_yscale('log')
            plt.xlim(axes['x_low'], axes['x_high']) 
            plt.ylim(axes['y_low'], axes['y_high'])
            path.set_title(col) 
            plt.tight_layout()
            
            if i == j or i > j:
                j += 1
            elif i < j:
                i += 1 
                j -= 1
        
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
def wrap_pearson_corr(df,label_column, alpha=.05,comparison_columns=None,correction_method='bonferroni',return_all = True, return_corrected_pvals = False):


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
    
    newdf = pd.DataFrame()
    for gene in comparison_columns:
        #create subset df with interacting gene/ gene (otherwise drop NaN drops everything)
        df_subset = df[[label_column,gene]]
        #do a linear regression to see if it's a meaningful association
        #dropna will remove rows with nan
        df_subset = df_subset.dropna(axis=0, how="any")
        count_row = df_subset.shape[0]
        if count_row > 30:
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
    corrected_pval = results[1]

    if return_all:
        if return_corrected_pvals:
            for i in range(0,len(comparisons)):
                newdf = newdf.append({'Comparison': comparisons[i],"Correlation": correlation[i], "P_value": corrected_pval[i]}, ignore_index=True)
                  
        if return_corrected_pvals == False:
            for i in range(0,len(comparisons)):
                newdf = newdf.append({'Comparison': comparisons[i],"Correlation": correlation[i],'P_value': pvals[i]}, ignore_index=True)

    '''Else only add significant comparisons'''
    if (return_all == False):
        
        if return_corrected_pvals:
            for i in range(0, len(reject)):
                if reject[i]:
                    newdf = newdf.append({'Comparison': comparisons[i],"Correlation": correlation[i],"P_value": corrected_pval[i]}, ignore_index=True)
        if return_corrected_pvals == False:
                for i in range(0, len(reject)):
                    if reject[i]:
                        newdf = newdf.append({'Comparison': comparisons[i],"Correlation": correlation[i],"P_value": pvals[i]}, ignore_index=True)
    
    
    newdf = newdf.sort_values(by= 'P_value', ascending=True)
    '''If results df is not empty, return it, else return None'''
    return newdf



'''
@Param df1: Dataframe. Contains numeric values (such as proteomics) for linear regression
@Param x_axis: String. Used as the label for the x-axis as well as the column name for the x-axis values.
@Param y_axis:String. Used as the label for the y-axis as well as the column name for the y-axis values.
@Param title: String. Used as title of figure
@Param ra_stats: Boolean. Default is False. If true it will print out the pearson correlation and p-value.
@Param x_coor: Float.  Default is 1.0. Is the x coordinate where the pearson correlation and p-value will print. 
@Param y_coor: Float.  Default is 1.0. Is the y coordinate where the pearson correlation and p-value will print. 
@Param show_plot: Boolean. Default is True. If true it will show the linear regression plot.
@Param save_file_name: String.File will not be saved unless a file name is specified. Plot is saved in current directory as png.

This fuction takes a dataframe with numeric values (such as proteomics) and performs a pearson coorelation analysis between two user specified columns within the dataframe. The function will then create the perason coorelation graph and can print the graph to the screen and save the figure depending on user input.
'''
def plot_pearson(df1,x_axis, y_axis, hue = "none", title = "", ra_stats = False, x_coor= 1.0 , y_coor = 1.0, show_plot = True, save_file_name = "file_name"):
    #format dfs 
    if hue != "none":   
      
        df1_subset = df1[[x_axis,y_axis,hue]]
        df1_subset = df1_subset.dropna(axis=0, how="any")
        count_row = df1_subset.shape[0]
        if count_row > 30:
            x1 = df1_subset[[x_axis]].values
            y1 = df1_subset[[y_axis]].values
            x1 = x1[:,0]
            y1 = y1[:,0]
            corr, pval = scipy.stats.pearsonr(x1,y1)

            sns.set(style="darkgrid")
            plt.rcParams["figure.figsize"] = (30,22)
            graph = sns.lmplot(x= x_axis, y= y_axis, data=df1_subset, hue= hue, fit_reg=False)
            sns.regplot(x=x1, y=y1, data=df1_subset,scatter = False)
            graph.set(title = title)
        else:
            return 0

    if hue == "none":
        
        df1_subset = df1[[x_axis,y_axis]]
        df1_subset = df1_subset.dropna(axis=0, how="any")
        count_row = df1_subset.shape[0]
        if count_row > 30:
            x1 = df1_subset[[x_axis]].values
            y1 = df1_subset[[y_axis]].values
            x1 = x1[:,0]
            y1 = y1[:,0]
            corr, pval = scipy.stats.pearsonr(x1,y1)

            sns.set(style="darkgrid")
            graph = sns.lmplot(x= x_axis, y= y_axis, data=df1_subset, fit_reg=False)
            sns.regplot(x=x1, y=y1, data=df1_subset,scatter = False)
            plt.title(label = title, fontsize = 30)
            plt.xlabel(x_axis, fontsize=20)
            plt.ylabel( y_axis, fontsize=20)
            plt.xticks(fontsize = 17)
            plt.yticks(fontsize = 17)
        else:
            return 0

    if ra_stats:
        pval = myTrunc(pval,5)
        corr = myTrunc(corr,3)
        plt.text(x_coor,y_coor, "Correlation: %s\nPvalue: %s"%(corr,pval), fontsize = 17)
        
    if save_file_name != "file_name":
        plt.tight_layout()
        plt.subplots_adjust(top=0.88)
        plt.savefig(save_file_name+'.png', dpi = 300)
    if show_plot:
        plt.show()
        plt.clf()
        plt.close()


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
    a4_dims = (40, 25) #dimensions for bigger plot
    fig, ax = plt.subplots(figsize=a4_dims)
    for position in range(0, number_of_df):
        r = ax.bar(x+(width*position), frequencies_for_each_df[position], width,label=names_of_df[position], alpha=.5, linewidth=0)



    ax.set_ylabel('Percent Sample', size = 60)
    ax.set_title(title)
    ax.set_xticks(x)
    ax.set_xticklabels(allLabels)
    ax.legend()
    ax.set_xticklabels(allLabels, size = 60)
    ax.set_ylim(0,100)
    ax.legend(prop={'size': 50})
    ax.yaxis.set_tick_params(labelsize=60)


    plt.setp(ax.get_xticklabels(),rotation= 45)
    plt.tight_layout()
    if save_to_path == None:
        plt.savefig("step_1.png",dpi = 300)
    else:
        plt.savefig(save_to_path)

    plt.show()
