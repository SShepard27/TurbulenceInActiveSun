import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt


#############################################################################
############################################################################# 
def plotsetup(font_size = 22):

  '''
  Graham Kerr
  NASA/GSFC & CUA
  Feb 2022

  NAME: plotsetup

  PURPOSE: Creates a dictionary with my preferred properties
           for MPL plots. 

  INPUTS: font_size -- flt
                       default = 22

  OUTPUTS: figprops -- a dictionary that is used with rcParams.update()

  NOTES:  

  '''

  font = {#'family': 'Avenir LT Std',
        'color':  'black',
        'weight': 'medium',
        'size': font_size,
        }

  plot_params = {'ytick.direction': 'in', 
               'xtick.direction': 'in', 
               'xtick.minor.visible': True,
               'ytick.minor.visible': True,
               'xtick.major.size': 10, 'xtick.minor.size': 5,
               'ytick.major.size': 10, 'ytick.minor.size': 5,
               'ytick.right': False,
               'xtick.top': False,
               'ytick.left':True,
               'xtick.bottom':False,
               'ytick.major.width': 1.5,
               'xtick.major.width': 1.5,
               'ytick.minor.width': 1.5,
               'xtick.minor.width': 1.5,
               'axes.linewidth': 1.5,
               'axes.spines.top': False,
               'axes.spines.bottom': True,
               'axes.spines.left': True,
               'axes.spines.right': False,
               'axes.titlepad' : 18 }

  plot_lg_params = {'legend.frameon': False}

  plot_dict = {'font.size':font['size'], 
                 # 'font.family':font['family'], 
                 'font.weight':font['weight'],
                 'ytick.direction': plot_params['ytick.direction'],
                 'xtick.direction': plot_params['xtick.direction'],
                 'xtick.minor.visible': plot_params['xtick.minor.visible'],
                 'ytick.minor.visible': plot_params['ytick.minor.visible'],
                 'ytick.major.size':  plot_params['ytick.major.size'], 
                 'ytick.minor.size':  plot_params['ytick.minor.size'],
                 'xtick.major.size':  plot_params['xtick.major.size'],                                
                 'xtick.minor.size':  plot_params['xtick.minor.size'],
                 'ytick.right': plot_params['ytick.right'],
                 'xtick.top': plot_params['xtick.top'],
                 'ytick.major.width': plot_params['ytick.major.width'],
                 'xtick.major.width': plot_params['xtick.major.width'],
                 'ytick.minor.width': plot_params['ytick.minor.width'],
                 'xtick.minor.width': plot_params['xtick.minor.width'],                    
                 'axes.linewidth': plot_params['axes.linewidth'],
                 'axes.spines.top' : plot_params['axes.spines.top'],
                 'axes.spines.bottom' : plot_params['axes.spines.bottom'],
                 'axes.spines.left' : plot_params['axes.spines.left'],
                 'axes.spines.right' : plot_params['axes.spines.right'],
                 'axes.titlepad' : plot_params['axes.titlepad'],
                 'legend.frameon': plot_lg_params['legend.frameon']}

  return plot_dict

####################################################################
####################################################################
def plotsetup_image(font_size = 22):

  '''
  Graham Kerr
  NASA/GSFC & CUA
  Feb 2022

  NAME: plotsetup_image

  PURPOSE: Creates a dictionary with my preferred properties
           for MPL plots, this version for images (so all spines are included)

  INPUTS: font_size -- flt
                       default = 22

  OUTPUTS: plot_dict -- a dictionary that is used with rcParams.update()

  NOTES:  
    
    
    '''

  font = {#'family': 'Avenir LT Std',
        'color':  'black',
        'weight': 'medium',
        'size': font_size,
        }

  plot_params = {'ytick.direction': 'inout', 
               'xtick.direction': 'inout', 
               'xtick.minor.visible': True,
               'ytick.minor.visible': True,
               'xtick.major.size': 10, 'xtick.minor.size': 5,
               'ytick.major.size': 10, 'ytick.minor.size': 5,
               'ytick.right': True,
               'xtick.top': True,
               'ytick.left':True,
               'xtick.bottom':True,
               'ytick.major.width': 1.5,
               'xtick.major.width': 1.5,
               'ytick.minor.width': 1.5,
               'xtick.minor.width': 1.5,
               'axes.linewidth': 1.5,
               'axes.spines.top': True,
               'axes.spines.bottom': True,
               'axes.spines.left': True,
               'axes.spines.right': True,
               'axes.titlepad' : 18 }

  plot_lg_params = {'legend.frameon': False}

  plot_dict = {'font.size':font['size'], 
                 # 'font.family':font['family'], 
                 'font.weight':font['weight'],
                 'ytick.direction': plot_params['ytick.direction'],
                 'xtick.direction': plot_params['xtick.direction'],
                 'xtick.minor.visible': plot_params['xtick.minor.visible'],
                 'ytick.minor.visible': plot_params['ytick.minor.visible'],
                 'ytick.major.size':  plot_params['ytick.major.size'], 
                 'ytick.minor.size':  plot_params['ytick.minor.size'],
                 'xtick.major.size':  plot_params['xtick.major.size'],                                
                 'xtick.minor.size':  plot_params['xtick.minor.size'],
                 'ytick.right': plot_params['ytick.right'],
                 'xtick.top': plot_params['xtick.top'],
                 'ytick.major.width': plot_params['ytick.major.width'],
                 'xtick.major.width': plot_params['xtick.major.width'],
                 'ytick.minor.width': plot_params['ytick.minor.width'],
                 'xtick.minor.width': plot_params['xtick.minor.width'],                    
                 'axes.linewidth': plot_params['axes.linewidth'],
                 'axes.spines.top' : plot_params['axes.spines.top'],
                 'axes.spines.bottom' : plot_params['axes.spines.bottom'],
                 'axes.spines.left' : plot_params['axes.spines.left'],
                 'axes.spines.right' : plot_params['axes.spines.right'],
                 'axes.titlepad' : plot_params['axes.titlepad'],
                 'legend.frameon': plot_lg_params['legend.frameon']}

  return plot_dict



