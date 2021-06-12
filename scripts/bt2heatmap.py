# -*- coding: utf-8 -*-
"""
Created on Fri Sep 21 15:37:26 2018

@author: qwzhou

This is a heatmap visualization script for meth on gene/repeat/histone regions.
"""

#import sys
#reload(sys)
#sys.setdefaultencoding('utf8')

import numpy as np
import argparse
from collections import OrderedDict
import numpy as np
import matplotlib
matplotlib.use('Agg')
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['svg.fonttype'] = 'none'
matplotlib.rcParams['font.size'] = 14.0
import matplotlib.pyplot as plt
from matplotlib.font_manager import FontProperties
import matplotlib.gridspec as gridspec
from matplotlib import ticker

import sys
sys.setrecursionlimit(1000000) # 设置递归深度，否则数据量较大并且cluster时会出问题
import os
from matplotlib.backends.backend_pdf import PdfPages
import argparse
from pandas import Series,DataFrame  #DataFrame通常来装二维的表格
import pandas as pd      #pandas是流行的做数据分析的包
import seaborn as sns

plt.switch_backend('agg')

y = []
z = []
k = []

dict_data={}
##col
num_samples=0
#row
num_groups=0

group_boundaries=[0]
sample_boundaries=[0] #0,600
group_labels=[]
sample_labels=[]
upstream=[]
downstream=[]
bin_size=[]
body=[] #0, 0
centermode = False
onlybody = False
totalline=0

#Create a dictionary. The keys and values are read from the file. The bond is Chr_ POS, the value is C / C + t
def readmatrix(filename, nsample, label, processed, plotidx):
    if plotidx == 0:
        global num_samples
        num_samples+=1
        if centermode == True:
            body.append(0)
            bin_size.append(10)
        else:
            body.append(2000)
            bin_size.append(10)  
        if onlybody == False:
            upstream.append(2000)
            downstream.append(2000)
        else:
            upstream.append(0)
            downstream.append(0)

    with open(filename,'r')as df:
        for line in df:
            if line.count('\n') == len(line):
                continue
            data = line.strip().split()
            for i in range(1, len(data)):
                if data[i] == 'nan':
                    dict_data.setdefault(data[0],[]).append(float(0))
                else:
                    dict_data.setdefault(data[0],[]).append(float(data[i]))
    
    linelen=len(data)-1
    if plotidx == 0:
        sample_boundaries.append((processed+1)*linelen )
        sample_labels.append(filename)
        #sample_labels.append(label[processed])

    print(processed+1, nsample)
    
    if(processed+1<nsample):
        return
    #print（dict_data）
    if(nsample>1):
        for key, val in dict_data.items():
            if(len(val)/linelen<nsample):
                del dict_data[key] 

    if(len(dict_data)==0):
        print("The merged matrix data is zero")
        exit()
    frame = pd.DataFrame.from_dict(dict_data,orient = 'index') #,columns=label
    global totalline
    totalline += len(dict_data)
    group_boundaries.append(totalline)
    global num_groups
    num_groups+=1
    group_labels.append(filename)

    print(sample_labels, group_boundaries)
    return frame

#Create a dictionary. The keys and values are read from the file. The bond is Chr_ POS, the value is C / C + t
def readmr(filename, nsample, label, processed, plotidx):
    if plotidx == 0:
        global num_samples
        num_samples+=1
        body.append(10)
        bin_size.append(10)
        upstream.append(0)
        downstream.append(0)
    
    with open(filename,'r')as df:
        for line in df:
            #If this line is a newline, skip it. Here we use the length of '\n' to find the empty line
            if line.count('\n') == len(line):
                continue
            #Clear the leading and trailing spaces (if any) for each line, and then separate them with "\t"
            #for kv in [line.strip().split('\t')]:
            kv=line.strip().split('\t')
            #print(kv)
            #按照键，把值写进去
            if(len(kv)>=7):
                dict_data.setdefault(kv[6],[]).append(float(float(kv[4])/int(kv[5])))
            else:
                dict_data.setdefault(kv[0]+"_"+kv[1]+"_"+kv[2],[]).append(float(float(kv[4])/int(kv[5])))

    linelen=1
    if plotidx == 0:
        sample_boundaries.append((processed+1)*linelen )
        sample_labels.append(filename)
        #sample_labels.append(label[processed])

    if(processed+1<nsample):
        return
    #print（dict_data）
    if(nsample>1):
        for key, val in dict_data.items():
            if(len(val)<nsample):
                del dict_data[key] 
    #这是把键读出来成为一个列表
    #columnsname=list(dict_data.keys())
    #print( columnsname, dict_data,len(dict_data), len(columnsname))
    #Create a dataframe, the column name is the key name, that is, Chr_ POS

    if(len(dict_data)==0):
        print("The merged matrix data is zero")
        exit()
    frame = pd.DataFrame.from_dict(dict_data,orient = 'index') #,columns=label
    #frame = pd.DataFrame(dict_data,index='',columns=columnsname).unstack(0)
    #print(frame) #, frame.values, label
    #把DataFrame输出到一个表，要行名字和列名字
    #frame.to_csv('file_out0.txt',index=True,header=True)
    global totalline
    totalline += len(dict_data)
    group_boundaries.append(totalline)
    global num_groups
    num_groups+=1
    group_labels.append(filename)

    print(sample_labels, group_boundaries)
    return frame

def readfile(filename, y, z, k):
    nline=0
    with open(filename, 'r') as fig:
        for line in fig:
            data = line.split()
            if nline == 0:
                y.append(list(map(float,data[1:])))
            elif nline == 1:
                z.append(list(map(float,data[1:])))
            elif nline == 2:
                k.append(list(map(float,data[1:])))
            nline=nline+1
 
######################################################
def find_martrix_max_value(data_matrix):  
    new_data=[]  
    for i in range(len(data_matrix)):  
        new_data.append(max(data_matrix[i]))
    return max(new_data)

#######################################################
def plotline(plotx, ploty, title, label, nsample, legend, filename, scale_ls, index_ls, ylabel, plotn):
    prosamp = 0
    #fig, ax = plt.subplots()
    ax = fig.add_subplot(1,3,plotn)
    while prosamp < nsample:
        ploty[prosamp] = [i*percentage for i in ploty[prosamp]]
        #for i,item in enumerate(y[prosamp]):
        #    if item >6:
        #        y[prosamp][i]=6
        ax.plot(plotx, ploty[prosamp], label=label[prosamp]) #,color="dodgerblue"
        prosamp = prosamp +1
    #dashes = [10, 5, 100, 5]
    #line1.set_dashes(dashes)   #   dash line
    # Remove the plot frame lines. They are unnecessary here.
    #ax.spines['top'].set_visible(False)
    #ax.spines['bottom'].set_visible(False)
    #ax.spines['right'].set_visible(False)
    #ax.spines['left'].set_visible(False)
    ax.xaxis.set_major_formatter(plt.FuncFormatter('{:.0f}'.format))
    #ax.yaxis.set_major_formatter(plt.FuncFormatter('{:.1f}%'.format))
    #plt.grid(True, 'major', 'y', ls='--', lw=.5, c='k', alpha=.3)
    plt.tick_params(axis='both', which='both', bottom=False, top=False,
                labelbottom=True, left=True, right=False, labelleft=True)
    #ax.axes.get_xaxis().set_visible(False)
    if legend == 1:
        plt.legend(loc='best', prop={'size': 12})   #   legend , loc is the legend location
    #plt.axvline(x=xlen/2-1, ls="--", color='black')
    
    #plt.axhline(y=0, xmin=0.05, xmax=0.5, linewidth=8, color='gray')
    #plt.axhline(y=0, xmin=0.5, xmax=0.505, linewidth=8, color='k' )
    #plt.axhline(y=0, xmin=0.505, xmax=0.95, linewidth=8, color='gray') 
    
    print(scale_ls, index_ls, len(plotx))
    #index_ls = ['-200bp','Start', "+200bp"]
    plt.xticks(scale_ls2,index_ls,color='k', size=12)
    ax.set_title(title,size=15)
    ax.set_ylabel(ylabel,size=15)
    #ax.set_ylabel('Methylation Level',size=15)
    maxy = 1
    maxy = find_martrix_max_value(ploty) * 1.1
    
    ax.set_ylim(0.0, maxy)
    plt.setp(ax.xaxis.get_ticklabels(), rotation='45')
    #plt.show()
    
    #filename2=filename + ".png"
    #plt.savefig(filename2, bbox_inches='tight')


def writableFile(string):
    """
    function that tests if a given path is writable
    """
    try:
        open(string, 'w').close()
        os.remove(string)
    except:
        msg = "{} file can't be opened for writing".format(string)
        raise argparse.ArgumentTypeError(msg)
    return string

def check_list_of_comma_values(value):
    if value is None:
        return None
    for foo in value:
        foo = value.split(",")
        if len(foo) < 2:
            raise argparse.ArgumentTypeError("%s is an invalid element of a list of comma separated values. "
                                             "Only argument elements of the following form are accepted: 'foo,bar'" % foo)
    return value

def check_float_0_1(value):
    v = float(value)
    if v < 0.0 or v > 1.0:
        raise argparse.ArgumentTypeError("%s is an invalid floating point value. It must be between 0.0 and 1.0" % value)
    return v

color_options = "', '".join([x for x in plt.colormaps() if not x.endswith('_r')])

def getArgs(args=None):
    parser = argparse.ArgumentParser(add_help=False)
    parser.add_argument("-f", "--mrfile", 
                        default=None, 
                        help="input methylevel files, wildtype.body.c*.txt", 
                        nargs='+')
    parser.add_argument("-m", "--matrixfile", 
                        default=None, 
                        help="input methylevel matrix files, wildtype.GENE.cg.txt", 
                        nargs='+')
    parser.add_argument("-l", "--samplesLabel", 
                        help="the label of the samples", 
                        nargs='+')
    parser.add_argument('-z', '--groupLabels', 
                        help='Labels for the regions plotted in the '
                        'heatmap. If more than one region is being '
                        'plotted, a list of labels separated by spaces is required. '
                        'If a label itself contains a space, then quotes are '
                        'needed. For example, --groupLabels label_1, "label 2". ',
                        nargs='+')
    parser.add_argument("-sl", "--startlabel", 
                        default="start", 
                        help="the start label of the samples", 
                        )
    parser.add_argument("-el", "--endlabel", 
                        default="end", 
                        help="the end label of the samples", 
                        )
    parser.add_argument("-pl", "--centerlabel", 
                        default=None, 
                        help="the center label of the samples", 
                        )
    parser.add_argument("--plotmatrix",
                        default='1x1',
                        help='1x1, default, row x col, order by columun, for exsample, 2x3 :'
                        'file1 file2 file3'
                        'file4 file5 file6',
                        )
    parser.add_argument('--outFileName', '-o',
                        help='Output file name.',
                        metavar='FILENAME',
                        default='myMeth.pdf',
                        type=writableFile,
                        required=True)

    parser.add_argument(
                        "-c", "--colorMap",
                        help='Color map to use for the heatmap. If more than one heatmap is being plotted the color '
                            'of each heatmap can be enter individually (e.g. `--colorMap Reds Blues`).'
                            'The available options are: \'' + color_options + '\'',
                        default=['vlag'],
                        nargs='+')
    parser.add_argument(
                        '--alpha',
                        default=1.0,
                        type=check_float_0_1,
                        help='The alpha channel (transparency) to use for the heatmaps. The default is 1.0 and values '
                            'must be between 0 and 1.')
    parser.add_argument(
                        '--colorList',
                        help='List of colors to use to create a colormap. For example, if `--colorList black,yellow,blue` '
                            'is set (colors separated by comas) then a color map that starts with black, continues to '
                            'yellow and finishes in blue is created. If this option is selected, it overrides the --colorMap '
                            'chosen. The list of valid color names can be seen here: '
                            'http://matplotlib.org/examples/color/named_colors.html  '
                            'The number of transitions is defined by the --colorNumber option.',
                        type=check_list_of_comma_values,
                        nargs='+')
    parser.add_argument(
                        '--colorNumber',
                        help='--colorList is required for an effect. This controls the '
                        'number of transitions from one color to the other. If --colorNumber is '
                        'the number of colors in --colorList then there will be no transitions '
                        'between the colors.',
                        type=int,
                        default=256)
    parser.add_argument(
                        '--missingDataColor',
                        default='white',
                        help='If --missingDataAsZero was not set, such cases '
                        'will be colored in white by default.' )
    parser.add_argument('--sortRegions',
                        help='Whether the heatmap should present the regions sorted. The default is '
                        'to sort in descending order based on the mean value per region.',
                        choices=["descend", "ascend", "no"],
                        default='descend')
    parser.add_argument('--sortUsing',
                        help='Indicate which method should be used for sorting. For each row the method is computed. ',
                        choices=["mean", "median", "max", "min", "sum"],
                        default='mean')
    parser.add_argument('--sortUsingSamples',
                        help='List of sample numbers (order as in matrix), which are used by --sortUsing for sorting. '
                        'If no value is set, it uses all samples. Example: --sortUsingSamples 1 3',
                        type=int, nargs='+')
    parser.add_argument('--linesAtTickMarks',
                        help='Draw dashed lines from all tick marks through the heatmap. '
                        'This is then similar to the dashed line draw at region bounds '
                        'when using a reference point and --sortUsing region_length',
                        action='store_true')
    parser.add_argument('--clusterUsingSamples',
                        help='List of sample numbers (order as in matrix), that are used for clustering by '
                        '--kmeans or --hclust if not given, all samples are taken into account for clustering. '
                        'Example: --ClusterUsingSamples 1 3',
                        type=int, nargs='+')
    parser.add_argument(
                        '--kmeans',
                        help='Number of clusters to compute. When this option is set, the matrix is split into clusters '
                        'using the k-means algorithm. Only works for data that is not grouped, otherwise only the first group will '
                        'be clustered. ',
                        type=int)
    parser.add_argument(
                        '--hclust',
                        help='Number of clusters to compute. When this option is set, then the matrix is split into clusters '
                        'using the hierarchical clustering algorithm, using "ward linkage". Only works for data that is not grouped, otherwise only the first '
                        'group will be clustered. --hclust could be very slow if you have >1000 regions. In those cases, you might prefer --kmeans or if more '
                        'clustering methods are required you can save the underlying matrix and run '
                        'the clustering using  other software. The plotting of the clustering may '
                        'fail with an error if a cluster has very few members compared to the '
                        'total number of regions.',
                        type=int)
    parser.add_argument("-s", "--scale",
                        default=[1.0, 1.0, 1.0],
                        nargs='+',
                        help='Maximum value for the Y-axis. Multiple values, separated by '
                            'spaces can be set for each profile. If the number of yMin values is smaller than'
                            'the number of plots, the values are recycled.',
                        type=float)
#    parser.add_argument("-xl", "--xlabel",
#                        default=False,
#                        nargs='+',
#                        help='[auto, bool(True/False), list-like, int, optional]. If True, plot the row names of the dataframe. If False, dont plot the row names.'
#                        'If list-like, plot these alternate labels as the xticklabels. If an integer, use the row names but'
#                        'plot only every n label. If "auto", try to densely plot non-overlapping labels.')
#    parser.add_argument("-yl", "--ylabel",
#                        default=False,
#                        nargs='+',
#                        help='[auto, bool(True/False), list-like, int, optional]. If True, plot the column names of the dataframe. If False, dont plot the column names.'
#                        'If list-like, plot these alternate labels as the xticklabels. If an integer, use the column names but'
#                        'plot only every n label. If "auto", try to densely plot non-overlapping labels.')
    parser.add_argument( '-t', '--title',
                        help='Title of the plot, to be printed on top of '
                        'the generated image. Leave blank for no title.',
                        nargs='+',
                        default='')
#    parser.add_argument( '-cl', '--cluster',
#        help='Plot a matrix using hierachical clustering to arrange the rows and columns.',
#        default=False)
    parser.add_argument('--zMin',
                        default=[0],
                        help='Values to anchor the colormap',
                        nargs='+',
                        type=float)
    parser.add_argument('--zMax',
                        default=[0.4],
                        help='Values to anchor the colormap, Maximum value for the heatmap.',
                        nargs='+',
                        type=float)
#    parser.add_argument('--zMid',
#                        default=None,
#                        help='The value at which to center the colormap when plotting divergant data. '
#                        'Using this parameter will change the default cmap if none is specified.',
#                        type=float)
    parser.add_argument("-ft", "--image_format",
                        default="pdf",
                        help="The file format, e.g. 'png', 'pdf', 'svg', ... "
                        "The behavior when this is unset is documented under fname.")
#    parser.add_argument("--clustermethod",
#                        default="single",
#                        choices=["single", "complete", "average","weighted",
#                        "centroid", "median", "ward"],
#                        help="Linkage method to use for calculating clusters, 'single', 'complete' "
#                        " 'average', 'weighted', 'centroid', 'median', 'ward' "
#                        " See scipy.cluster.hierarchy.linkage documentation for more information: "
#                        "https://docs.scipy.org/doc/scipy/reference/generated/scipy.cluster.hierarchy.linkage.html"
#                        )
#    parser.add_argument("--zscore",
#                        default=None,
#                        help="Either 0 (rows) or 1 (columns). Whether or not to calculate z-scores for "
#                        "the rows or the columns. Z scores are: z = (x - mean)/std, so values in each row (column) will "
#                        "get the mean of the row (column) subtracted, then divided by the standard deviation of the row (column)."
#                        " This ensures that each row (column) has mean of 0 and variance of 1."
#                        )
    parser.add_argument('--perGroup',
                        help='The default is to plot all groups of regions by '
                        'sample. Using this option instead plots all samples by '
                        'group of regions. Note that this is only useful if you '
                        'have multiple groups of regions. by sample rather than '
                        'group.',
                        action='store_true')
    parser.add_argument("--dpi",
                        default=100,
                        help="Set the DPI to save the figure. default: 100",
                        type=int)
    parser.add_argument("--figsize",
                        default="1.5x11",
                        help="Set the figure size to save the figure. [with]x[height], default: 1.5x11")
    parser.add_argument(
                        '--boxAroundHeatmaps',
                        help='By default black boxes are plot around heatmaps. This can be turned off '
                            'by setting --boxAroundHeatmaps no',
                        default='yes')
    parser.add_argument("--help", "-h", action="help",
                        help="show this help message and exit")

    return parser

def getProfileTicks(referencePointLabel, startLabel, endLabel, idx):
    """
    
    returns the position and labelling of the xticks that
    correspond to the heatmap

    As of deepTools 3, the various parameters can be lists, in which case we then need to index things (the idx parameter)

    As of matplotlib 3 the ticks in the heatmap need to have 0.5 added to them.

    As of matplotlib 3.1 there is no longer padding added to all ticks. Reference point ticks will be adjusted by width/2
    or width for spacing and the last half of scaled ticks will be shifed by 1 bin so the ticks are at the beginning of bins.
    """
    w = bin_size
    b = upstream
    a = downstream
    if idx is not None:
        w = w[idx]
        b = b[idx]
        a = a[idx]

    d = 0
    c=0
    m = body
    if idx is not None:
        m = m[idx]

    if b < 1e5:
        quotient = 1000
        symbol = 'Kb'
    else:
        quotient = 1e6
        symbol = 'Mb'

    if m == 0:
        xticks = [(k / w) for k in [0, b - 0.5 * w, b + a - w]]
        xtickslabel = ['{0:.1f}'.format(-(float(b) / quotient)),
                       referencePointLabel,
                       '{0:.1f}{1}'.format(float(a) / quotient, symbol)]
    else:
        xticks_values = [0]
        xtickslabel = []

        # only if upstream region is set, add a x tick
        if b > 0:
            xticks_values.append(b)
            xtickslabel.append('{0:.1f}'.format(-(float(b) / quotient)))

        xtickslabel.append(startLabel)

        # set the x tick for the body parameter, regardless if
        # upstream is 0 (not set)
        if c > 0:
            xticks_values.append(b + c)
            xtickslabel.append("")

        if d > 0:
            xticks_values.append(b + c + m)
            xtickslabel.append("")

        # We need to subtract the bin size from the last 2 point so they're placed at the beginning of the bin
        xticks_values.append(b + c + m + d - w)
        xtickslabel.append(endLabel)

        if a > 0:
            xticks_values.append(b + c + m + d + a - w)
            xtickslabel.append('{0:.1f}{1}'.format(float(a) / quotient, symbol))

        xticks = [(k / w) for k in xticks_values]
        xticks = [max(x, 0) for x in xticks]

    return xticks, xtickslabel

def get_matrix(mymatrix, group, sample):
    """
    Returns a sub matrix from the large
    matrix. Group and sample are ids,
    thus, row = 0, col=0 get the first group
    of the first sample.

    Returns
    -------
    dictionary containing the matrix,
    the group label and the sample label
    """
    group_start = group_boundaries[group]
    group_end = group_boundaries[group + 1]
    sample_start = sample_boundaries[sample]
    sample_end = sample_boundaries[sample + 1]

    #print('gs', group_start, group_end, group, sample_start, sample_end,sample)
    return {'matrix': np.ma.masked_invalid(mymatrix[group_start:group_end, :][:, sample_start:sample_end]),
            'group': group_labels[group],
            'sample': sample_labels[sample]}

def prepare_layout(mymatrix, heatmapsize, showSummaryPlot, showColorbar, perGroup, colorbar_position):
    """
    prepare the plot layout
    as a grid having as many rows
    as samples (+1 for colobar)
    and as many rows as groups (or clusters) (+1 for profile plot)
    """
    heatmapwidth, heatmapheight = heatmapsize

    numcols = num_samples
    numrows = num_groups
    if perGroup:
        numcols, numrows = numrows, numcols

    # the rows have different size depending
    # on the number of regions contained in the
    if perGroup:
        # heatmap
        height_ratio = np.array([np.amax(np.diff(group_boundaries))] * numrows)
        # scale ratio to sum = heatmapheight
        height_ratio = heatmapheight * (height_ratio.astype(float) / height_ratio.sum())
    else:
        # heatmap
        height_ratio = np.diff(group_boundaries)
        # scale ratio to sum = heatmapheight
        height_ratio = heatmapheight * (height_ratio.astype(float) / height_ratio.sum())

    # convert the height_ratio from numpy array back to list
    height_ratio = height_ratio.tolist()
    # the width ratio is equal for all heatmaps
    width_ratio = [heatmapwidth] * numcols

    if showColorbar:
        if colorbar_position == 'below':
            numrows += 2  # a spacer needs to be added to avoid overlaps
            height_ratio += [4 / 2.54]  # spacer
            height_ratio += [1 / 2.54]
        else:
            numcols += 1
            width_ratio += [1 / 2.54]

    if showSummaryPlot:
        numrows += 2  # plus 2 because a spacer is added
        # make height of summary plot
        # proportional to the width of heatmap
        sumplot_height = heatmapwidth
        spacer_height = heatmapwidth / 8
        # scale height_ratios to convert from row
        # numbers to heatmapheigt fractions
        height_ratio = np.concatenate([[sumplot_height, spacer_height], height_ratio])

    #print(numrows, numcols, height_ratio, width_ratio)
    grids = gridspec.GridSpec(numrows, numcols, height_ratios=height_ratio, width_ratios=width_ratio)

    return grids

def getTicks(reference_point_label, idx, startLabel, endLabel):
    """
    This is essentially a wrapper around getProfileTicks to accomdate the fact that each column has its own ticks.
    """
    if reference_point_label is not None:
        xticks, xtickslabel = getProfileTicks(reference_point_label[idx], startLabel, endLabel, idx)
    else:
        xticks, xtickslabel = getProfileTicks(reference_point_label, startLabel, endLabel, idx)
    return xticks, xtickslabel

def plotheatmaper(mymatrix, outFileName,
               colorMapDict={'colorMap': ['binary'], 'missingDataColor': 'white', 'alpha': 1.0},
               plotTitle='',
               xAxisLabel='', yAxisLabel='', regionsLabel='',
               zMin=None, zMax=None,
               yMin=None, yMax=None,
               averageType='median',
               reference_point_label=None,
               startLabel='TSS', endLabel="TES",
               heatmapHeight=25,
               heatmapWidth=7.5,
               perGroup=False, whatToShow='plot, heatmap and colorbar',
               plot_type='lines',
               linesAtTickMarks=False,
               image_format=None,
               legend_location='upper-left',
               box_around_heatmaps=True,
               label_rotation=0.0,
               dpi=200,
               interpolation_method='auto',
               num_samples=1):

    colorbar_pergroup = True

    if reference_point_label is not None:
        reference_point_label = [reference_point_label] * num_samples
    startLabel = startLabel
    endLabel = endLabel

    matrix_flatten = None
    if zMin is None:
        matrix_flatten = mymatrix.flatten()
        # try to avoid outliers by using np.percentile
        zMin = np.percentile(matrix_flatten, 1.0)
        if np.isnan(zMin):
            zMin = [None]
        else:
            zMin = [zMin]  # convert to list to support multiple entries
    elif 'auto' in zMin:
        matrix_flatten = mymatrix.flatten()
        auto_min = np.percentile(matrix_flatten, 1.0)
        if np.isnan(auto_min):
            auto_min = None
        new_mins = [float(x) if x != 'auto' else auto_min for x in zMin]
        zMin = new_mins
    else:
        new_mins = [float(x) for x in zMin]
        zMin = new_mins

    if zMax is None:
        if matrix_flatten is None:
            matrix_flatten = mymatrix.flatten()
        # try to avoid outliers by using np.percentile
        zMax = np.percentile(matrix_flatten, 98.0)
        if np.isnan(zMax) or zMax <= zMin[0]:
            zMax = [None]
        else:
            zMax = [zMax]
    elif 'auto' in zMax:
        matrix_flatten = mymatrix.flatten()
        auto_max = np.percentile(matrix_flatten, 98.0)
        if np.isnan(auto_max):
            auto_max = None
        new_maxs = [float(x) if x != 'auto' else auto_max for x in zMax]
        zMax = new_maxs
    else:
        new_maxs = [float(x) for x in zMax]
        zMax = new_maxs
    if (len(zMin) > 1) & (len(zMax) > 1):
        for index, value in enumerate(zMax):
            if value <= zMin[index]:
                sys.stderr.write("Warnirng: In bigwig {}, the given zmin ({}) is larger than "
                                 "or equal to the given zmax ({}). Thus, it has been set "
                                 "to None. \n".format(index + 1, zMin[index], value))
                zMin[index] = None

    if yMin is None:
        yMin = [None]
    if yMax is None:
        yMax = [None]
    if not isinstance(yMin, list):
        yMin = [yMin]
    if not isinstance(yMax, list):
        yMax = [yMax]

    plt.rcParams['font.size'] = 14.0
    fontP = FontProperties()

    showSummaryPlot = False
    showColorbar = False

    if whatToShow == 'plot and heatmap':
        showSummaryPlot = True
    elif whatToShow == 'heatmap and colorbar':
        showColorbar = True
    elif whatToShow == 'plot, heatmap and colorbar':
        showSummaryPlot = True
        showColorbar = True

    # colormap for the heatmap
    if colorMapDict['colorMap']:
        cmap = []
        for color_map in colorMapDict['colorMap']:
            cmap.append(plt.get_cmap(color_map))
            cmap[-1].set_bad(colorMapDict['missingDataColor'])  # nans are printed using this color

    if colorMapDict['colorList'] and len(colorMapDict['colorList']) > 0:
        # make a cmap for each color list given
        cmap = []
        for color_list in colorMapDict['colorList']:
            cmap.append(matplotlib.colors.LinearSegmentedColormap.from_list(
                'my_cmap', color_list.replace(' ', '').split(","), N=colorMapDict['colorNumber']))
            cmap[-1].set_bad(colorMapDict['missingDataColor'])  # nans are printed using this color

    if not colorbar_pergroup and (len(cmap) > 1 or len(zMin) > 1 or len(zMax) > 1):
        # position color bar below heatmap when more than one
        # heatmap color is given
        colorbar_position = 'below'
    else:
        colorbar_position = 'side'

    grids = prepare_layout(mymatrix, (heatmapWidth, heatmapHeight),
                           showSummaryPlot, showColorbar, perGroup, colorbar_position)

    # figsize: w,h tuple in inches
    figwidth = heatmapWidth / 2.54
    figheight = heatmapHeight / 2.54
    if showSummaryPlot:
        # the summary plot ocupies a height
        # equal to the fig width
        figheight += figwidth

    numsamples = num_samples
    if perGroup:
        num_cols = num_groups
    else:
        num_cols = numsamples
    total_figwidth = figwidth * num_cols
    if showColorbar:
        if colorbar_position == 'below':
            figheight += 1 / 2.54
        else:
            total_figwidth += 1 / 2.54

    fig = plt.figure(figsize=(total_figwidth, figheight))
    fig.suptitle(plotTitle, y=1 - (0.06 / figheight))

    # color map for the summary plot (profile) on top of the heatmap
    cmap_plot = plt.get_cmap('jet')
    numgroups = num_groups
    if perGroup:
        color_list = cmap_plot(np.arange(num_samples) / num_samples)
    else:
        color_list = cmap_plot(np.arange(numgroups) / numgroups)
    alpha = colorMapDict['alpha']

    first_group = 0  # helper variable to place the title per sample/group
    for sample in range(num_samples):
        #print("sample ", sample, range(num_samples))
        sample_idx = sample
        for group in range(numgroups):
            group_idx = group
            # add the respective profile to the
            # summary plot
            #print(group, sample)
            sub_matrix = get_matrix(mymatrix, group, sample)
            #print('submatrix', sub_matrix)
            if showSummaryPlot:
                if perGroup:
                    sample_idx = sample + 2  # plot + spacer
                else:
                    group += 2  # plot + spacer
                first_group = 1

            if perGroup:
                ax = fig.add_subplot(grids[sample_idx, group])
                # the remainder (%) is used to iterate
                # over the available color maps (cmap).
                # if the user only provided, lets say two
                # and there are 10 groups, colormaps they are reused every
                # two groups.
                cmap_idx = group_idx % len(cmap)
                zmin_idx = group_idx % len(zMin)
                zmax_idx = group_idx % len(zMax)
            else:
                ax = fig.add_subplot(grids[group, sample])
                # see above for the use of '%'
                
                if colorbar_pergroup:
                    cmap_idx = group % len(cmap)
                    zmin_idx = group % len(zMin)
                    zmax_idx = group % len(zMax)
                else:
                    cmap_idx = sample % len(cmap)
                    zmin_idx = sample % len(zMin)
                    zmax_idx = sample % len(zMax)

            if group == first_group and not showSummaryPlot and not perGroup:
                title = sample_labels[sample]
                ax.set_title(title)

            if box_around_heatmaps is False:
                # Turn off the boxes around the individual heatmaps
                ax.spines['top'].set_visible(False)
                ax.spines['right'].set_visible(False)
                ax.spines['bottom'].set_visible(False)
                ax.spines['left'].set_visible(False)
            rows, cols = sub_matrix['matrix'].shape
            # if the number of rows is too large, then the 'nearest' method simply
            # drops rows. A better solution is to relate the threshold to the DPI of the image
            if interpolation_method == 'auto':
                if rows >= 1000:
                    interpolation_method = 'bilinear'
                else:
                    interpolation_method = 'nearest'

            # if np.clip is not used, then values of the matrix that exceed the zmax limit are
            # highlighted. Usually, a significant amount of pixels are equal or above the zmax and
            # the default behaviour produces images full of large highlighted dots.
            # If interpolation='nearest' is used, this has no effect
            sub_matrix['matrix'] = np.clip(sub_matrix['matrix'], zMin[zmin_idx], zMax[zmax_idx])
            img = ax.imshow(sub_matrix['matrix'],
                            aspect='auto',
                            interpolation=interpolation_method,
                            origin='upper',
                            vmin=zMin[zmin_idx],
                            vmax=zMax[zmax_idx],
                            cmap=cmap[cmap_idx],
                            alpha=alpha,
                            extent=[0, cols, rows, 0])
            img.set_rasterized(True)
            # plot border at the end of the regions

            if perGroup:
                ax.axes.set_xlabel(sub_matrix['group'])
                if sample < num_samples - 1:
                    ax.axes.get_xaxis().set_visible(False)
            else:
                ax.axes.get_xaxis().set_visible(False)
                ax.axes.set_xlabel(xAxisLabel)
            ax.axes.set_yticks([])
            if perGroup and group == 0:
                ax.axes.set_ylabel(sub_matrix['sample'])
            elif not perGroup and sample == 0:
                ax.axes.set_ylabel(sub_matrix['group'])

            # Plot vertical lines at tick marks if desired
            if linesAtTickMarks:
                xticks_heat, xtickslabel_heat = getTicks(reference_point_label, sample , startLabel, endLabel)
                xticks_heat = [x + 0.5 for x in xticks_heat]  # There's an offset of 0.5 compared to the profile plot
                if np.ceil(max(xticks_heat)) != float(sub_matrix['matrix'].shape[1]):
                    tickscale = float(sub_matrix['matrix'].shape[1]) / max(xticks_heat)
                    xticks_heat_use = [x * tickscale for x in xticks_heat]
                else:
                    xticks_heat_use = xticks_heat
                for x in xticks_heat_use:
                    ax.axvline(x=x, color='black', linewidth=0.5, dashes=(3, 2))

            # add labels to last block in a column
            if (perGroup and sample == numsamples - 1) or \
               (not perGroup and group_idx == numgroups - 1):

                # add xticks to the bottom heatmap (last group)
                ax.axes.get_xaxis().set_visible(True)
                xticks_heat, xtickslabel_heat = getTicks(reference_point_label, sample , startLabel, endLabel)
                xticks_heat = [x + 0.5 for x in xticks_heat]  # There's an offset of 0.5 compared to the profile plot
                if np.ceil(max(xticks_heat)) != float(sub_matrix['matrix'].shape[1]):
                    tickscale = float(sub_matrix['matrix'].shape[1]) / max(xticks_heat)
                    xticks_heat_use = [x * tickscale for x in xticks_heat]
                    ax.axes.set_xticks(xticks_heat_use)
                else:
                    ax.axes.set_xticks(xticks_heat)
                ax.axes.set_xticklabels(xtickslabel_heat, size=8)

                # align the first and last label
                # such that they don't fall off
                # the heatmap sides
                ticks = ax.xaxis.get_major_ticks()
                ticks[0].label1.set_horizontalalignment('left')
                ticks[-1].label1.set_horizontalalignment('right')

                ax.get_xaxis().set_tick_params(
                    which='both',
                    top=False,
                    direction='out')

                if not colorbar_pergroup and showColorbar and colorbar_position == 'below':
                    # draw a colormap per each heatmap below the last block
                    if perGroup:
                        col = group_idx
                    else:
                        col = sample
                    ax = fig.add_subplot(grids[-1, col])
                    tick_locator = ticker.MaxNLocator(nbins=3)
                    cbar = fig.colorbar(img, cax=ax, alpha=alpha, orientation='horizontal', ticks=tick_locator)
                    labels = cbar.ax.get_xticklabels()
                    ticks = cbar.ax.get_xticks()
                    if ticks[0] == 0:
                        # if the label is at the start of the colobar
                        # move it a bit inside to avoid overlapping
                        # with other labels
                        labels[0].set_horizontalalignment('left')
                    if ticks[-1] == 1:
                        # if the label is at the end of the colobar
                        # move it a bit inside to avoid overlapping
                        # with other labels
                        labels[-1].set_horizontalalignment('right')
                    # cbar.ax.set_xticklabels(labels, rotation=90)
            
            if colorbar_pergroup and sample_idx+1 == numsamples:
                ax = fig.add_subplot(grids[group_idx, -1])
                tick_locator = ticker.MaxNLocator(nbins=3)
                fig.colorbar(img, cax=ax, alpha=alpha, ticks=tick_locator)

    if not colorbar_pergroup and showColorbar and colorbar_position != 'below':
        if showSummaryPlot:
            # we don't want to colorbar to extend
            # over the profiles and spacer top rows
            grid_start = 2
        else:
            grid_start = 0

        ax = fig.add_subplot(grids[grid_start:, -1])
        fig.colorbar(img, cax=ax, alpha=alpha)

    if box_around_heatmaps:
        plt.subplots_adjust(wspace=0.10, hspace=0.025, top=0.85, bottom=0, left=0.04, right=0.96)
    else:
        #  When no box is plotted the space between heatmaps is reduced
        plt.subplots_adjust(wspace=0.05, hspace=0.01, top=0.85, bottom=0, left=0.04, right=0.96)

    plt.savefig(outFileName, bbox_inches='tight', pdd_inches=0, dpi=dpi, format=image_format)
    plt.close()

def sort_groups(mymatrix, sort_using='mean', sort_method='no', sample_list=None):
    """
    Sorts and rearranges the submatrices according to the
    sorting method given.
    """
    if sort_method == 'no':
        return

    if (sample_list is not None) and (len(sample_list) > 0):
        # get the ids that correspond to the selected sample list
        idx_to_keep = []
        for sample_idx in sample_list:
            idx_to_keep += range(sample_boundaries[sample_idx], sample_boundaries[sample_idx + 1])

        matrix = mymatrix[:, idx_to_keep]

    else:
        matrix = mymatrix

    # compute the row average:
    if sort_using == 'mean':
        matrix_avgs = np.nanmean(matrix, axis=1)
    elif sort_using == 'mean':
        matrix_avgs = np.nanmean(matrix, axis=1)
    elif sort_using == 'median':
        matrix_avgs = np.nanmedian(matrix, axis=1)
    elif sort_using == 'max':
        matrix_avgs = np.nanmax(matrix, axis=1)
    elif sort_using == 'min':
        matrix_avgs = np.nanmin(matrix, axis=1)
    elif sort_using == 'sum':
        matrix_avgs = np.nansum(matrix, axis=1)
    else:
        sys.exit("{} is an unsupported sorting method".format(sort_using))

    # order per group
    _sorted_regions = []
    _sorted_matrix = []
    for idx in range(len(group_labels)):
        start = group_boundaries[idx]
        end = group_boundaries[idx + 1]
        order = matrix_avgs[start:end].argsort()
        if sort_method == 'descend':
            order = order[::-1]
        _sorted_matrix.append(mymatrix[start:end, :][order, :])

    return np.vstack(_sorted_matrix)
    #print(np.vstack(_sorted_matrix))
    #self.set_sorting_method(sort_method, sort_using)

def hmcluster(mymatrix, k, evaluate_silhouette=True, method='kmeans', clustering_samples=None):
    matrix = np.asarray(mymatrix)
    matrix_to_cluster = matrix
    if clustering_samples is not None:
        assert all(i > 0 for i in clustering_samples),\
            "all indices should be bigger than or equal to 1."
        assert all(i <= len(sample_labels) for i in
                    clustering_samples),\
            "each index should be smaller than or equal to {}(total "\
            "number of samples.)".format(len(sample_labels))

        clustering_samples = np.asarray(clustering_samples) - 1

        samples_cols = []
        for idx in clustering_samples:
            samples_cols += range(sample_boundaries[idx],
                                   sample_boundaries[idx + 1])
        print(clustering_samples, sample_boundaries[idx], sample_boundaries[idx+1])
        matrix_to_cluster = matrix_to_cluster[:, samples_cols]
    if np.any(np.isnan(matrix_to_cluster)):
        # replace nans for 0 otherwise kmeans produces a weird behaviour
        sys.stderr.write("*Warning* For clustering nan values have to be replaced by zeros \n")
        matrix_to_cluster[np.isnan(matrix_to_cluster)] = 0

    if method == 'kmeans':
        from scipy.cluster.vq import vq, kmeans

        centroids, _ = kmeans(matrix_to_cluster, k)
        # order the centroids in an attempt to
        # get the same cluster order
        cluster_labels, _ = vq(matrix_to_cluster, centroids)

    if method == 'hierarchical':
        # normally too slow for large data sets
        from scipy.cluster.hierarchy import fcluster, linkage
        Z = linkage(matrix_to_cluster, method='ward', metric='euclidean')
        cluster_labels = fcluster(Z, k, criterion='maxclust')
        # hierarchical clustering labels from 1 .. k
        # while k-means labels 0 .. k -1
        # Thus, for consistency, we subtract 1
        cluster_labels -= 1

    # sort clusters
    _clustered_mean = []
    _cluster_ids_list = []
    for cluster in range(k):
        cluster_ids = np.flatnonzero(cluster_labels == cluster)
        _cluster_ids_list.append(cluster_ids)
        _clustered_mean.append(matrix_to_cluster[cluster_ids, :].mean())

    # reorder clusters based on mean
    cluster_order = np.argsort(_clustered_mean)[::-1]
    # create groups using the clustering
    group_labels = []
    group_boundaries = [0]
    _clustered_matrix = []
    cluster_number = 1
    for cluster in cluster_order:
        group_labels.append("cluster_{}".format(cluster_number))
        cluster_number += 1
        cluster_ids = _cluster_ids_list[cluster]
        group_boundaries.append(group_boundaries[-1] +
                                        len(cluster_ids))
        _clustered_matrix.append(matrix[cluster_ids, :])

    mymatrix = np.vstack(_clustered_matrix)

    return mymatrix

def set_group_labels(new_labels):
    """ sets new labels for groups
    """
    global group_labels
    if len(new_labels) != len(group_labels):
        raise ValueError("length new labels != length original labels")
    group_labels = new_labels

def set_sample_labels(new_labels):
    """ sets new labels for groups
    """
    global sample_labels
    if len(new_labels) != len(sample_labels):
        raise ValueError("length new labels != length original labels")
    sample_labels = new_labels

if __name__ == '__main__':
    #args = vars(ap.parse_args())
    args = getArgs().parse_args()
    #methpoint(sys.argv[1],sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5])
    nsample=0
    if(args.mrfile!=None):
        nsample = len(args.mrfile)
    msample=0
    if(args.matrixfile!=None):
        msample = len(args.matrixfile)
    
    Nsample = nsample + msample

    label=[] #args.label
    filename=args.outFileName
    zmin=args.zMin
    zmax=args.zMax
    
    #if args.zMid == None:
    #    zmid = (zmin+zmax)/2
    #else:
    #    zmid=args.zMid
    #print(zmin, zmid, zmax)
    plotmatrix = args.plotmatrix.split("x")
    if nsample > 1:
        if plotmatrix[1] == '1':
            plotmatrix[1] = nsample
    if msample > 1:
        if plotmatrix[1] == '1':
            plotmatrix[1] = msample
    print(plotmatrix)

    if args.centerlabel != None:
        centermode = True

    dict_data={}
    matrix_row = int(plotmatrix[0])
    matrix_col = int(plotmatrix[1])
    if(nsample>0):
        if len(plotmatrix) > 0:
            for idx in range(0, matrix_row): # row         
                dict_data={}
                for x in range(idx*matrix_col, (idx+1)*matrix_col ): # col
                    newframe = readmr(args.mrfile[x], matrix_col, label, x-idx*matrix_col, idx)
                if idx == 0:
                    frame = newframe
                else:
                    frame = pd.concat([frame, newframe], axis=0)
        else:
            for x in range(0, nsample):
                frame = readmr(args.mrfile[x], nsample, label, x, 0)
    if(msample>0):
        if len(plotmatrix) > 0:
            for idx in range(0, matrix_row): # row         
                dict_data={}
                for x in range(idx*matrix_col, (idx+1)*matrix_col ): # col
                    newframe = readmatrix(args.matrixfile[x], matrix_col, label, x-idx*matrix_col, idx)
                if idx == 0:
                    frame = newframe
                else:
                    frame = pd.concat([frame, newframe], axis=0)
        else:
            for x in range(0, msample):
                frame = readmatrix(args.matrixfile[x], msample, label, x, 0)

    print("NxM", sample_boundaries, num_samples, num_groups)

    mydpi=args.dpi
    figsize=args.figsize.split('x')
    if(len(figsize)<2):
        print("Unexpected figure size, should be [with]x[height], eg. 3x10")
    fwith=float(figsize[0])
    fheight=float(figsize[1])
    plt.figure(dpi=mydpi, figsize=[fwith,fheight])
    print(frame.head())
    #color=args.color
    xticklabels = False #args.xlabel
    yticklabels = False #args.ylabel
    zscore = False #args.zscore
    frame=frame.values # for matplotlib imshow
    
    cols=frame.shape[0]
    rows=frame.shape[1]
    print(frame, rows, cols )
    image_format=args.image_format
    colorNumber=256
    missingDataColor=args.missingDataColor
    colorList=args.colorList #black,yellow,blue
    colorMap=args.colorMap
    colormap_dict = {'colorMap': colorMap,
                     'colorList': colorList,
                     'colorNumber': colorNumber,
                     'missingDataColor': missingDataColor,
                     'alpha': args.alpha}
    sortRegions=args.sortRegions
    mykmeans=args.kmeans ## need para
    myhclust=args.hclust
    clusterUsingSamples=args.clusterUsingSamples
    if sortRegions == 'keep':
        sortRegions = 'no'  # These are the same thing
    if mykmeans is not None:
        frame=hmcluster(frame, mykmeans, method='kmeans', clustering_samples=clusterUsingSamples)
    elif myhclust is not None:
        print("Performing hierarchical clustering."
              "Please note that it might be very slow for large datasets.\n")
        hmcluster(frame, myhclust, method='hierarchical', clustering_samples=clusterUsingSamples)
    #print(frame)

    if args.groupLabels:
        set_group_labels(args.groupLabels)

    if args.samplesLabel and len(args.samplesLabel):
        print(args.samplesLabel, "hhh", sample_labels)
        set_sample_labels(args.samplesLabel)

    sortUsing=args.sortUsing
    if sortRegions != 'no':
        sortUsingSamples = []
        if args.sortUsingSamples is not None:
            for i in args.sortUsingSamples:
                if (i > 0 and i <= num_samples):
                    sortUsingSamples.append(i - 1)
                else:
                    exit("The value {0} for --sortSamples is not valid. Only values from 1 to {1} are allowed.".format(args.sortUsingSamples, num_samples))
            print('Samples used for ordering within each group: ', sortUsingSamples)
        print(sortUsing, sortRegions, sortUsingSamples)
        frame = sort_groups(frame, sort_using=sortUsing,
                    sort_method=sortRegions,
                    sample_list=sortUsingSamples)
    print(frame)
    print("num_samples ", num_samples)
    plotheatmaper(frame, filename,
               colorMapDict=colormap_dict,
               plotTitle='',
               xAxisLabel='', yAxisLabel='', regionsLabel='',
               zMin=None, zMax=zmax,
               yMin=None, yMax=None,
               averageType='median',
               reference_point_label=args.centerlabel,
               startLabel=args.startlabel, endLabel=args.endlabel,
               heatmapHeight=25,
               heatmapWidth=7.5,
               perGroup=args.perGroup, whatToShow='heatmap and colorbar',
               plot_type='lines',
               linesAtTickMarks=args.linesAtTickMarks,
               image_format=image_format,
               legend_location='upper-left',
               box_around_heatmaps=args.boxAroundHeatmaps,
               label_rotation=0.0,
               dpi=200,
               interpolation_method='gaussian',
               num_samples=num_samples)
    #['auto', 'nearest', 'bilinear', 'bicubic', 'gaussian']

    mysns=False
    if mysns:
        fig, ax = plt.subplots() #fig.add_subplot(1,1,1)
        #img = ax.imshow(frame,aspect='auto',interpolation="bilinear",
        #                origin='upper')
        #img.set_rasterized(True)
        if args.cluster == True:
            clustermethod = args.clustermethod
            print("cluster:", clustermethod)
            sns.clustermap(frame,figsize=[fwith,fheight], z_score=zscore,
            method=clustermethod, row_cluster=True, col_cluster=False, cmap=colorMap, center=zmid, vmin=zmin, vmax=zmax, xticklabels=xticklabels, yticklabels=yticklabels) #, pivot_kws=None, method='average', metric='euclidean', z_score=None, 
        else:
            sns.heatmap(data=frame,shading="gouraud",
            cmap=colorMap,center=zmid,vmin=zmin,vmax=zmax,annot=False,fmt="d", xticklabels=xticklabels, yticklabels=yticklabels)
        #standard_scale=None, figsize=None, cbar_kws=None, row_cluster=True, col_cluster=True, 
        #row_linkage=None, col_linkage=None, row_colors=None, col_colors=None, mask=None
        #annot 显示注释 fmt='.2f'
        #cbar_pos=(0.02, 0.8, 0.05, 0.18), #(left, bottom, width, height) 图例位置
        plt.tight_layout()
        plt.savefig(filename, dpi=mydpi, format=image_format)
        plt.close()
