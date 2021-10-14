"""
Code written by Jeremy Cushman and updated/documented by Jorge Torres.

Todo:
    * Document what each function does.
    * Update sensitivities, etc.
    * Make sure to add word "Preliminary if using CUORE's unpublished results"

"""

from numpy import cos, sin, exp, pi, sqrt, arcsin
import numpy as np
import matplotlib.pylab as plt
from random import uniform, betavariate
import cmath


def nuM_range(m,hiearchy,includeSterile, p12, p13, p14, nSamples = 1000):

    """
    This function calculates the effective Majorana mass,
    bete decay kinetic mass, and sum mass of all neutrino flavors
    given the lightest neutrino mass.
    The calculation depends on the hiearchy specified (IH, or NH)
    and can possibly include the simple 3+1 sterile neutrino model if desired.

    In the calculation, CP-violating phase, two majorana phases (3 if sterile)
    are randomly drawn from (0, 2pi).
    Mixing parameters are drawn randomly from (min, max) using a beta function
    that strongly favors the end points.


    Suffix notation:
    _0: best fit values
    _3sigma: 3sigma related
    _U: upper boundary
    _L: lower boundary

    Parameters
    ----------
    m : array_like
        1D array of masses to be sampled [units: meV].

    hierarchy : string_like
        "NH": uses normal ordering scenario
        "IH": uses inverted hierarchy scenario

    includeSterile : boolean
        Whether to include paramenters P14 for extra neutrino

    p12, p13, p14 : array_like
        Mixing parameters. They must follow the following format:

        pij = [[angle_bestfit, angle_min, angle_max], [deltaMassSquared_bestfit, deltaMassSqaured_min, deltaMassSqaured_max]]

    nSamples : int
          nSamples specifies how many random sampling points are used. nSamples = 100K takes ~1 hour

    Returns
    -------
    array of arrays
        Arrays containing majorana effective mass,...
    """


    ## Define parameters
    T12_0 = p12[0][0]
    C12_0= cos(T12_0)
    S12_0= sin(T12_0)

    T14_0 = p14[0][0]
    S14_0 = sin(T14_0)

    T13_0 = p13[0][0]
    C13_0 = cos(T13_0)
    S13_0 = sin(T13_0)

    ## Squared mass difference
    dM2_12_0 = p12[1][0]
    dM2_13_0 = p13[1][0]
    dM2_14_0 = p14[1][0]

    mbb_U_0= 1E-3
    mbb_L_0= 1E6
    mb_0 =0
    sumM_0 =0

    mbb_U_3sigma= 1E-3
    mbb_L_3sigma= 1E6
    mb_U_3sigma =1E-3
    mb_L_3sigma =1E6
    sumM_U_3sigma =1E-3
    sumM_L_3sigma =1E6

    ##Define masses from square mass differences depending on the hierarchy scenario
    if hiearchy == 'NH':
        m1_0 = m
        m2_0 = sqrt(m**2 + dM2_12_0)
        m3_0 = sqrt(m**2 + dM2_13_0)
        m4_0 = sqrt(m**2 + dM2_14_0)

    else:
        m1_0 = sqrt(m**2+ dM2_13_0)
        m2_0 = sqrt(m**2 + dM2_13_0 + dM2_12_0)
        m3_0 = m
        m4_0 = sqrt(m**2 + dM2_13_0 + dM2_14_0)
    if not includeSterile:
        m4_0 = 0

    for ii in range(nSamples):
        alpha = uniform(0, 2*pi)
        beta  = uniform(0, 2*pi)
        gamma = uniform(0, 2*pi)

        T12 = p12[0][1] + betavariate(0.1,0.1) * (p12[0][2]-p12[0][1])
        T13 = p13[0][1] + betavariate(0.1,0.1) * (p13[0][2]-p13[0][1])
        T14 = p14[0][1] + betavariate(0.1,0.1) * (p14[0][2]-p14[0][1])

        dM2_12 = p12[1][1] + betavariate(0.1,0.1) * (p12[1][2]-p12[1][1])
        dM2_13 = p13[1][1] + betavariate(0.1,0.1) * (p13[1][2]-p13[1][1])
        dM2_14 = p14[1][1] + betavariate(0.1,0.1) * (p14[1][2]-p14[1][1])

        C12 = cos(T12); S12 = sin(T12)
        C13 = cos(T13); S13 = sin(T13)
        C14 = cos(T14); S14 = sin(T14)

        if hiearchy == 'NH':
            m1 = m
            m2 = sqrt(m**2 + dM2_12)
            m3 = sqrt(m**2 + dM2_13)
            m4 = sqrt(m**2 + dM2_14)

        else:
            m1 = sqrt(m**2+ dM2_13)
            m2 = sqrt(m**2 + dM2_13 + dM2_12)
            m3 = m
            m4 = sqrt(m**2 + dM2_13 + dM2_14)
        if  not includeSterile:
            m4 = 0

        m_bb_0 = abs(\
                m1_0 * C12_0**2 * C13_0**2\
                + m2_0 * S12_0**2 * C13_0**2 * cmath.exp(complex(0,alpha)) \
                + m3_0 * S13_0**2 * cmath.exp(complex(0,beta))\
                + m4_0 * S14_0**2 * cmath.exp(complex(0,gamma))\
                )

        if m_bb_0 > mbb_U_0: mbb_U_0 = m_bb_0
        if m_bb_0 < mbb_L_0: mbb_L_0 = m_bb_0

        m_bb_3sigma = abs(\
                m1 * C12**2 * C13**2 \
                + m2 * S12**2 * C13 **2 * cmath.exp(complex(0,alpha)) \
                + m3 * S13**2 * cmath.exp(complex(0,beta)) \
                + m4 * S14**2 * cmath.exp(complex(0,gamma)) \
                )

        if m_bb_3sigma > mbb_U_3sigma: mbb_U_3sigma = m_bb_3sigma
        if m_bb_3sigma < mbb_L_3sigma: mbb_L_3sigma = m_bb_3sigma

        #### kinetic neutrino mass from beta Decay
        mb_0 = sqrt(\
                    m1_0**2 * C12_0**2 * C13_0**2 \
                    + m2_0**2 * S12_0**2 * C13_0**2 \
                    + m3_0**2 * S13_0**2 \
                    + m4_0**2 * S14_0**2 \
                    )

        mb_3sigma = sqrt(\
                         m1**2 * C12**2 * C13**2 \
                         + m2**2 * S12**2 * C13**2 \
                         + m3**2 * S13**2 \
                         + m4**2 * S14**2 \
                         )

        if mb_3sigma > mb_U_3sigma: mb_U_3sigma = mb_3sigma
        if mb_3sigma < mb_L_3sigma: mb_L_3sigma = mb_3sigma

        #### Sum of neutrino mass ############
        sumM_0 = m1_0 + m2_0 + m3_0 + m4_0
        sumM_3sigma = m1 + m2 + m3 + m4

        if sumM_3sigma > sumM_U_3sigma: sumM_U_3sigma = sumM_3sigma
        if sumM_3sigma < sumM_L_3sigma: sumM_L_3sigma = sumM_3sigma

    return [mbb_L_0, mbb_U_0, mbb_L_3sigma, mbb_U_3sigma], \
           [mb_0, mb_L_3sigma, mb_U_3sigma], \
           [sumM_0, sumM_L_3sigma, sumM_U_3sigma]


def AddExperimentalLimits(IH, NH, xMin, xMax, isotopes=None, yMin=-1, yMax=-1):
    """
    This function adds specified limits and sensitivities to Lobster plot

    Parameters
    ----------
    ax : axes object
        matplotlib object where plot will be drawn.

    xArray : array_like
        This is the x-axis (energy in meV) array.

    y4DArray : array_like
        This is a 4D array containing the central value spread (min, max, respectively) in the first two internal arrays y4DArray[:,1], y4DArray[:,2], and the 3sigma spread values (min, max, respectively) in y4DArray[:,2], y4DArray[:,3].

    col : string
        choose a color based on matplotlib color palettes: https://matplotlib.org/stable/gallery/color/named_colors.html

    xlab : string
        x-axis label

    ylab : string
        y-axis label

    Returns
    -------
    plot object
        It returns a plot object.
    """

    arrowXscale = 1.2
    axSpan = xMax-xMin
    arrColor = '#2B4970'
    if(isotopes==None):
        isotopes = ["Xe", "Te", "Ge"]
        
    ##Te limit from CUORE
    
    if("Te" in isotopes):
        mbb_min_Te = 90
        mbb_max_Te = 305
        A_Te = 130
#         IH.axhspan(mbb_min_Te, mbb_max_Te, xmin = (A_Te)/axSpan, xmax = (A_Te)/axSpan, lw=0,fc=arrColor,ec=arrColor, alpha=0.3)
        Teline = IH.hlines(mbb_min_Te, A_Te-2, A_Te+2, color=arrColor, label='$^{130}$Te limit (CUORE [Prelim.])')
        IH.errorbar(A_Te, mbb_min_Te,yerr = mbb_max_Te-mbb_min_Te-20, lolims=True,  color=arrColor)
        IH.text(A_Te-5, mbb_min_Te, '$^{130}$Te', color=arrColor,fontsize='medium', ha="right", fontweight="book")

    ## Ge limit from GERDA (2020): https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.125.252502
    if("Ge" in isotopes):
        mbb_min_Ge = 79
        mbb_max_Ge = 180
        A_Ge = 76
#         IH.axhspan(mbb_min_Ge, mbb_max_Ge, xmin = (A_Ge)/axSpan, xmax = (A_Ge)/axSpan, lw=0, ec=arrColor,fill=None, alpha=0.3, hatch='\\\\\\')
        Geline = IH.hlines(mbb_min_Ge, A_Ge-2, A_Ge+2, color=arrColor, label='$^{76}$Ge limit (GERDA)', linestyle='-')
        IH.errorbar(A_Ge, mbb_min_Ge, yerr=mbb_max_Ge-mbb_min_Ge-20, lolims=True,  color=arrColor)
        IH.text(A_Ge+3, mbb_min_Ge, '$^{76}$Ge', color=arrColor,fontsize='medium', ha="left", fontweight="book")

    ##Xe limit from KamLand-Zen (2016): https://journals.aps.org/prl/pdf/10.1103/PhysRevLett.117.082503
    if("Xe" in isotopes):
        mbb_min_Xe = 61
        mbb_max_Xe = 165
        A_Xe = 136
        
#         IH.axhspan(mbb_min_Xe, mbb_max_Xe, xmin = (A_Xe)/axSpan, xmax = (A_Xe+10)/axSpan, lw=0, ec='#AA7F39',fill=None, alpha=0.3, hatch='///')
        Xeline = IH.hlines(mbb_min_Xe, A_Xe-2, A_Xe+2, color=arrColor, label='$^{136}$Xe limit (KamLAND-Zen)', linestyle='-')
        IH.errorbar(A_Xe, mbb_min_Xe, yerr=mbb_max_Xe-mbb_min_Xe-20, lolims=True,  color=arrColor)
        IH.text(A_Xe+3, mbb_min_Xe, '$^{136}$Xe', color=arrColor,fontsize='medium', ha="left", fontweight="book")

    ## Mo limit from CUPID-Mo (2021): https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.126.1
    if("Mo" in isotopes):
        mbb_min_Mo = 300
        mbb_max_Mo = 500
        A_Mo = 100
        
#         IH.axhspan(mbb_min_Mo, mbb_max_Mo,lw=0, xmin = 100, xmax = 150, ec='#AA7F39',fill=None, alpha=0.3, hatch='///')
        Moline = IH.hlines(mbb_min_Mo, A_Mo-2, A_Mo+2, color=arrColor, label='$^{100}$Mo limit (CUPID-Mo)', linestyle='-')
        IH.errorbar(A_Mo, mbb_min_Mo, yerr=mbb_max_Mo-mbb_min_Mo-20, lolims=True,  color=arrColor)
        IH.text(A_Mo+3, mbb_min_Mo, '$^{100}$Mo', color=arrColor,fontsize='medium', ha="left", fontweight="book")


    ##
    ## Sensitivities ##
    ##
    
#     IH.axhspan(50, 130,lw=0, fc='#2B4970',ec='#2B4970', alpha=0.3)
#     IH.axhline(50, color='#2B4970', label='CUORE Sensitivity')
#     IH.errorbar(xMin*arrowXscale, 50, yerr=50, lolims=True,  color='#2B4970')
#     IH.text(xMin*pow(arrowXscale,2), 50*1.2, 'CUORE Sensitivity', color='#2B4970',fontsize='small')

#     IH.legend(handles=[Teline, Geline, Xeline], loc=3,prop={'size':9})

    # ON the right panel
#     NH.axhspan(270, 650, lw=0,fc='#2B4970',ec='#AA9B39', alpha=0.3, label='Te-130 Limit')
#     #NH.axhspan(270, 650, lw=0,ec='#AA9B39',fill=None, alpha=0.3, hatch='///')
#     NH.axhline(270, color='#2B4970', label='Te-130 Limit')

#     NH.axhspan(200, 400, lw=0,ec='#AA7F39',fill=None, alpha=0.3, hatch='\\\\\\')
#     NH.axhline(200, color='#AA7F39', label='Ge-76 Limit', linestyle='--')

#     NH.axhspan(120, 250,lw=0, ec='#AA7F39',fill=None, alpha=0.3, hatch='///')
#     NH.axhline(120, color='#AA7F39', label='Xe-136 Limit', linestyle='-')

#     IH.axhspan(50, 130, lw=0,fc='#2B4970',ec='#452F74', alpha=0.3, label='CUORE Sensitivity')
#     IH.axhspan(50, 130,lw=0, ec='#452F74',fill=None, alpha=0.3, hatch='\\\\\\')
#     IH.axhline(50, color='#2B4970', label='CUORE Sensitivity')

    # IH.set_title('3 Flavors, Inverted Hierarchy')
    NH.set_xlabel('Atomic number')

#     IH.set_xlim(xMin, xMax)
    NH.set_xlim(xMin, xMax)

    if yMin>0 and yMax >0:
        IH.set_ylim(yMin, yMax)
    else:
        IH.set_ylim(xMin, xMax)


def massContour(ax, xArray, y4DArray, col, xlab, ylab=''):
    """
    This function makes the "Lobster plot"

    Parameters
    ----------
    ax : axes object
        matplotlib object where plot will be drawn.

    xArray : array_like
        This is the x-axis (energy in meV) array.

    y4DArray : array_like
        This is a 4D array containing the central value spread (min, max, respectively) in the first two internal arrays y4DArray[:,1], y4DArray[:,2], and the 3sigma spread values (min, max, respectively) in y4DArray[:,2], y4DArray[:,3].

    col : string
        choose a color based on matplotlib color palettes: https://matplotlib.org/stable/gallery/color/named_colors.html

    xlab : string
        x-axis label

    ylab : string
        y-axis label

    Returns
    -------
    plot object
        It returns a plot object.
    """

    if y4DArray.shape[1]==3:
        ax.plot(xArray, y4DArray[:,0], color=col)
        ax.fill_between(xArray,y4DArray[:,1], y4DArray[:,2], color=col, alpha=0.2)
    elif y4DArray.shape[1]==4:
        ax.fill_between(xArray,y4DArray[:,0], y4DArray[:,1], color=col, alpha=0.6)
        ax.fill_between(xArray,y4DArray[:,2], y4DArray[:,3], color=col, alpha=0.2)
    else:
        return

    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlabel(xlab)
    ax.set_ylabel(ylab)

def addSensitivity(ax, sens_mbb_min, sens_mbb_max, xMin, experiment_name = "Experiment", color = "C4", hatch = "///"):
    """
    This function adds a sensitivity band with limits [sens_mbb_min, sens_mbb_max] for a given experiment.

    Parameters
    ----------
    ax : axes object
        matplotlib object where plot will be drawn.
        
    sens_mbb_min, sens_mbb_max : float_like
        The sensitivity limits

    experiment_name : string
        Name of the experiment.

    col : string
        choose a color based on matplotlib color palettes: https://matplotlib.org/stable/gallery/color/named_colors.html

    Returns
    -------
    plot object
        It returns a plot object with sensitivities.
    """
        
    arrowXscale = 1.2
    ax.axhspan(sens_mbb_min, sens_mbb_max,lw=0, color = color,fc=color,ec=color, alpha=0.1, hatch = "\\")
    ax.annotate(s='', xy=(1,sens_mbb_min), xytext=(1,sens_mbb_max), arrowprops=dict(arrowstyle='<-', color=color))
    ax.text(xMin*(pow(arrowXscale,2)), sens_mbb_min*1.2, experiment_name, color=color,fontsize='large')
    