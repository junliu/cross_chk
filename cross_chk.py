#!/usr/bin/env python
# -*- coding: UTF-8 -*-
#

import wx
import os
import re
import sys
import time
from wx.lib.pubsub import pub
from threading import Thread

import matplotlib
matplotlib.use('WXAgg')
import matplotlib.gridspec as gridspec
from matplotlib.backends.backend_wxagg import FigureCanvasWxAgg \
as FigureCanvas, NavigationToolbar2WxAgg as NavigationToolbar

sys.path.append('shared')
import numpy as np
from astropy.io import fits as pyfits
import ObsData
import MathLib as ML


matplotlib.rcParams['axes.facecolor'] = '#545352'
matplotlib.rcParams['axes.edgecolor'] = 'white'
matplotlib.rcParams['figure.facecolor'] = '#424140'
matplotlib.rcParams['figure.edgecolor'] = 'white'
matplotlib.rcParams['grid.color'] = 'white'
matplotlib.rcParams['text.color'] = 'white'
matplotlib.rcParams['axes.labelcolor'] = 'white'
matplotlib.rcParams['xtick.color'] = 'white'
matplotlib.rcParams['ytick.color'] = 'white'
matplotlib.rcParams['text.usetex'] = False
matplotlib.rcParams['text.latex.unicode'] = False


class OptDialog(wx.Dialog):

  def __init__(self, parent, *args, **kwds):

    kwds['style'] = wx.DEFAULT_DIALOG_STYLE
    wx.Dialog.__init__(self, parent, *args, **kwds)
    self.p = parent
    self.ST_TELE = wx.StaticText(self, -1, 'Telescope')
    self.ST_DIA = wx.StaticText(self, -1, 'Diameter')
    self.ST_LON = wx.StaticText(self, -1, 'Longtitude')
    self.ST_LAT = wx.StaticText(self, -1, 'Latitude')
    self.C_TELE = wx.Choice(self, -1, choices=['NSRT', 'Effelsberg', 'GBT', 'TM65', 'Custom'])
    self.TC_DIA = wx.TextCtrl(self, -1, '26.0')
    self.TC_LON = wx.TextCtrl(self, -1, '87.178108')
    self.TC_LAT = wx.TextCtrl(self, -1, '43.4709389')
    self.ST_BLANK = wx.StaticText(self, -1, '')
    self.ST_M = wx.StaticText(self, -1, 'Meter')
    self.ST_DEG0 = wx.StaticText(self, -1, 'Deg.')
    self.ST_DEG1 = wx.StaticText(self, -1, 'Deg.')
    self.SZ_TELE_staticbox = wx.StaticBox(self, -1, '')
    self.ST_FREQ = wx.StaticText(self, -1, 'Freq.')
    self.ST_TCAL = wx.StaticText(self, -1, 'Tcal')
    self.TC_FREQ = wx.TextCtrl(self, -1, '4800.0')
    self.TC_TCAL = wx.TextCtrl(self, -1, '1.0')
    self.ST_MHZ = wx.StaticText(self, -1, 'MHz')
    self.ST_K = wx.StaticText(self, -1, 'K')
    self.SZ_BE_staticbox = wx.StaticBox(self, -1, '')
    self.BT_FLAG = wx.ToggleButton(self, -1, 'Flag Options ON')
    self.BT_BEAM = wx.Button(self, -1, 'Beam Width', style=wx.BU_EXACTFIT)
    self.ST_DOFF = wx.StaticText(self, -1, 'D_Offsets')
    self.ST_DWIDTH = wx.StaticText(self, -1, 'D_Width')
    self.ST_ERR = wx.StaticText(self, -1, 'Rela_Err')
    self.ST_SYMM = wx.StaticText(self, -1, 'D_Symmetry')
    self.ST_AVG = wx.StaticText(self, -1, 'D_Average')
    self.TC_BEAM = wx.TextCtrl(self, -1, '550.0')
    self.TC_DOFF = wx.TextCtrl(self, -1, '80.0')
    self.TC_DWIDTH = wx.TextCtrl(self, -1, '10')
    self.TC_ERR = wx.TextCtrl(self, -1, '10')
    self.TC_DSYMM = wx.TextCtrl(self, -1, '10')
    self.TC_DAVG = wx.TextCtrl(self, -1, '15')
    self.ST_AC0 = wx.StaticText(self, -1, '\'')
    self.ST_AC1 = wx.StaticText(self, -1, '\'')
    self.ST_PER0 = wx.StaticText(self, -1, '%')
    self.ST_PER1 = wx.StaticText(self, -1, '%')
    self.ST_PRE2 = wx.StaticText(self, -1, '%')
    self.ST_PER3 = wx.StaticText(self, -1, '%')
    self.SZ_FLAG_staticbox = wx.StaticBox(self, -1, '')
    self.BT_HELP = wx.Button(self, wx.ID_HELP, '', style=wx.BU_EXACTFIT)
    self.BT_CANCEL = wx.Button(self, wx.ID_CANCEL, '', style=wx.BU_EXACTFIT)
    self.BT_OK = wx.Button(self, wx.ID_OK, '', style=wx.BU_EXACTFIT)
    self.SZ_BUTTON_staticbox = wx.StaticBox(self, -1, '')

    self.__set_properties()
    self.__do_layout()

    self.Bind(wx.EVT_CHOICE, self.ON_CTELE, self.C_TELE)
    self.Bind(wx.EVT_TOGGLEBUTTON, self.ON_FLAG, self.BT_FLAG)
    self.Bind(wx.EVT_BUTTON, self.ON_BEAM, self.BT_BEAM)
    self.Bind(wx.EVT_BUTTON, self.ON_HELP, self.BT_HELP)
    self.Bind(wx.EVT_BUTTON, self.ON_CANCEL, self.BT_CANCEL)
    self.Bind(wx.EVT_BUTTON, self.ON_OK, self.BT_OK)

  def __set_properties(self):

    self.SetTitle('Option')

    self.C_TELE.SetSelection(self.C_TELE.FindString(self.p.tele_name))
    self.TC_DIA.SetValue(str(self.p.diam))
    self.TC_LON.SetValue(str(self.p.obs_lon))
    self.TC_LAT.SetValue(str(self.p.obs_lat))

    self.BT_FLAG.SetValue(self.p.flag)
    if self.p.flag:
      self.BT_FLAG.SetBackgroundColour(wx.Colour(0, 255, 0))
      self.BT_FLAG.SetLabel('Flag Options ON')
    else:
      self.BT_FLAG.SetBackgroundColour(wx.Colour(255, 0, 0))
      self.BT_FLAG.SetLabel('Flag Options OFF')

    self.TC_FREQ.SetValue('%0.1f' %self.p.freq)
    self.TC_TCAL.SetValue('%0.1f' %self.p.tcal)

    self.TC_BEAM.SetValue('%0.1f' %self.p.beam)
    self.TC_DOFF.SetValue('%0.1f' %self.p.doff)
    self.TC_DWIDTH.SetValue('%d' %(self.p.dwidth*100))
    self.TC_ERR.SetValue('%d' %(self.p.rerr*100))
    self.TC_DSYMM.SetValue('%d' %(self.p.dsymm*100))
    self.TC_DAVG.SetValue('%d' %(self.p.davg*100))

  def __do_layout(self):

    MainSizer = wx.BoxSizer(wx.VERTICAL)
    self.SZ_BUTTON_staticbox.Lower()
    SZ_BUTTON = wx.StaticBoxSizer(self.SZ_BUTTON_staticbox, wx.HORIZONTAL)
    self.SZ_FLAG_staticbox.Lower()
    SZ_FLAG = wx.StaticBoxSizer(self.SZ_FLAG_staticbox, wx.HORIZONTAL)
    SZ_FLAG_2 = wx.BoxSizer(wx.VERTICAL)
    SZ_FLAG_1 = wx.BoxSizer(wx.VERTICAL)
    SZ_FLAG_0 = wx.BoxSizer(wx.VERTICAL)
    SZ_BTFLAG = wx.BoxSizer(wx.HORIZONTAL)
    self.SZ_BE_staticbox.Lower()
    SZ_BE = wx.StaticBoxSizer(self.SZ_BE_staticbox, wx.HORIZONTAL)
    SZ_BE_2 = wx.BoxSizer(wx.VERTICAL)
    SZ_BE_1 = wx.BoxSizer(wx.VERTICAL)
    SZ_BE_0 = wx.BoxSizer(wx.VERTICAL)
    self.SZ_TELE_staticbox.Lower()
    SZ_TELE = wx.StaticBoxSizer(self.SZ_TELE_staticbox, wx.HORIZONTAL)
    SZ_TELE_2 = wx.BoxSizer(wx.VERTICAL)
    SZ_TELE_1 = wx.BoxSizer(wx.VERTICAL)
    SZ_TELE_0 = wx.BoxSizer(wx.VERTICAL)
    SZ_TELE_0.Add(self.ST_TELE, 1, wx_lpt('l|r|ar|acv'), 5)
    SZ_TELE_0.Add(self.ST_DIA, 1, wx_lpt('l|ar|t|ar|acv'), 5)
    SZ_TELE_0.Add(self.ST_LON, 1, wx_lpt('l|r|t|ar|acv'), 5)
    SZ_TELE_0.Add(self.ST_LAT, 1, wx_lpt('l|r|t|ar|acv'), 5)
    SZ_TELE.Add(SZ_TELE_0, 0, wx_lpt('l|t|e'), 5)
    SZ_TELE_1.Add(self.C_TELE, 0, wx_lpt('b|e|ach|acv'), 5)
    SZ_TELE_1.Add(self.TC_DIA, 1, wx_lpt('b|e|ach|acv'), 5)
    SZ_TELE_1.Add(self.TC_LON, 1, wx_lpt('b|e|ach|acv'), 5)
    SZ_TELE_1.Add(self.TC_LAT, 1, wx_lpt('b|e|ach|acv'), 5)
    SZ_TELE.Add(SZ_TELE_1, 0, wx_lpt('l|r|e'), 5)
    SZ_TELE_2.Add(self.ST_BLANK, 1, wx_lpt('l|r|acv'), 5)
    SZ_TELE_2.Add(self.ST_M, 1, wx_lpt('l|r|t|acv'), 5)
    SZ_TELE_2.Add(self.ST_DEG0, 1, wx_lpt('l|r|t|acv'), 5)
    SZ_TELE_2.Add(self.ST_DEG1, 1, wx_lpt('l|r|t|acv'), 5)
    SZ_TELE.Add(SZ_TELE_2, 0, wx_lpt('r|t|e'), 5)
    MainSizer.Add(SZ_TELE, 0, wx_lpt('l|r|e'), 2)
    SZ_BE_0.Add(self.ST_FREQ, 1, wx_lpt('l|r|ar|acv'), 5)
    SZ_BE_0.Add(self.ST_TCAL, 1, wx_lpt('l|r|t|ar|acv'), 5)
    SZ_BE.Add(SZ_BE_0, 0, wx_lpt('l|t|e'), 5)
    SZ_BE_1.Add(self.TC_FREQ, 1, wx_lpt('b|e|ach|acv'), 5)
    SZ_BE_1.Add(self.TC_TCAL, 1, wx_lpt('b|e|ach|acv'), 5)
    SZ_BE.Add(SZ_BE_1, 1, wx_lpt('l|r|e'), 5)
    SZ_BE_2.Add(self.ST_MHZ, 1, wx_lpt('l|r|acv'), 5)
    SZ_BE_2.Add(self.ST_K, 1, wx_lpt('l|r|t|acv'), 5)
    SZ_BE.Add(SZ_BE_2, 0, wx_lpt('r|t|e'), 5)
    MainSizer.Add(SZ_BE, 0, wx_lpt('l|r|e'), 2)
    SZ_BTFLAG.Add(self.BT_FLAG, 1, wx_lpt('t|ach|acv'), 3)
    MainSizer.Add(SZ_BTFLAG, 0, wx_lpt('l|r|t|e'), 2)
    SZ_FLAG_0.Add(self.BT_BEAM, 1, wx_lpt('l|r|ar|acv'), 5)
    SZ_FLAG_0.Add(self.ST_DOFF, 1, wx_lpt('l|r|t|ar|acv'), 5)
    SZ_FLAG_0.Add(self.ST_DWIDTH, 1, wx_lpt('l|r|t|ar|acv'), 5)
    SZ_FLAG_0.Add(self.ST_ERR, 1, wx_lpt('l|r|t|ar|acv'), 5)
    SZ_FLAG_0.Add(self.ST_SYMM, 1, wx_lpt('l|r|t|ar|acv'), 5)
    SZ_FLAG_0.Add(self.ST_AVG, 1, wx_lpt('l|r|t|ar|acv'), 5)
    SZ_FLAG.Add(SZ_FLAG_0, 0, wx_lpt('l|e'), 5)
    SZ_FLAG_1.Add(self.TC_BEAM, 1, wx_lpt('b|e|ach|acv'), 5)
    SZ_FLAG_1.Add(self.TC_DOFF, 1, wx_lpt('b|e|ach|acv'), 5)
    SZ_FLAG_1.Add(self.TC_DWIDTH, 1, wx_lpt('b|e|ach|acv'), 5)
    SZ_FLAG_1.Add(self.TC_ERR, 1, wx_lpt('b|e|ach|acv'), 5)
    SZ_FLAG_1.Add(self.TC_DSYMM, 1, wx_lpt('b|e|ach|acv'), 5)
    SZ_FLAG_1.Add(self.TC_DAVG, 1, wx_lpt('b|e|ach|acv'), 5)
    SZ_FLAG.Add(SZ_FLAG_1, 0, wx_lpt('l|r|e'), 5)
    SZ_FLAG_2.Add(self.ST_AC0, 1, wx_lpt('l|r|acv'), 5)
    SZ_FLAG_2.Add(self.ST_AC1, 1, wx_lpt('l|r|t|acv'), 5)
    SZ_FLAG_2.Add(self.ST_PER0, 1, wx_lpt('l|r|t|acv'), 5)
    SZ_FLAG_2.Add(self.ST_PER1, 1, wx_lpt('l|r|t|acv'), 5)
    SZ_FLAG_2.Add(self.ST_PRE2, 1, wx_lpt('l|r|t|acv'), 5)
    SZ_FLAG_2.Add(self.ST_PER3, 1, wx_lpt('l|r|t|acv'), 5)
    SZ_FLAG.Add(SZ_FLAG_2, 0, wx_lpt('r|t|e'), 5)
    MainSizer.Add(SZ_FLAG, 0, wx_lpt('l|r|e'), 2)
    SZ_BUTTON.Add(self.BT_HELP, 1, wx_lpt('l|r|b|ach|acv'), 5)
    SZ_BUTTON.Add(self.BT_CANCEL, 1, wx_lpt('l|r|b|ach|acv'), 5)
    SZ_BUTTON.Add(self.BT_OK, 1, wx_lpt('l|r|b|ach|acv'), 5)
    MainSizer.Add(SZ_BUTTON, 0, wx_lpt('l|r|e'), 2)
    self.SetSizer(MainSizer)
    MainSizer.Fit(self)
    self.Layout()

  def ON_CTELE(self, event):
    tele = self.C_TELE.GetString(self.C_TELE.GetCurrentSelection())
    if tele == 'NSRT':
      self.TC_DIA.SetValue('26.0')
      self.TC_LON.SetValue('87.178108')
      self.TC_LAT.SetValue('43.4709389')

    elif tele == 'Effelsberg':
      self.TC_DIA.SetValue('100.0')
      self.TC_LON.SetValue('27.178108')
      self.TC_LAT.SetValue('48.4709389')

    else:
      pass

    event.Skip()

  def ON_FLAG(self, event):
    self.p.flag = self.BT_FLAG.GetValue()
    if self.p.flag:
      self.BT_FLAG.SetBackgroundColour(wx.Colour(0, 255, 0))
      self.BT_FLAG.SetLabel('Flag Options ON')
    else:
      self.BT_FLAG.SetBackgroundColour(wx.Colour(255, 0, 0))
      self.BT_FLAG.SetLabel('Flag Options OFF')

    event.Skip()

  def ON_BEAM(self, event):
    freq = float(self.TC_FREQ.GetValue())
    dia = float(self.TC_DIA.GetValue())
    val = 1.12*3600*300*180/np.pi/freq/dia
    self.TC_BEAM.SetValue('%0.1f' %val)
    self.TC_DOFF.SetValue('%0.1f' %(0.15*val))

  def ON_HELP(self, event):
    print('No help available currenetly')
    event.Skip()

  def ON_CANCEL(self, event):
    event.Skip()

  def ON_OK(self, event):

    self.p.tele_name = self.C_TELE.GetString(self.C_TELE.GetCurrentSelection())
    self.p.diam = float(self.TC_DIA.GetValue())
    self.p.obs_lon = float(self.TC_LON.GetValue())
    self.p.obs_lat = float(self.TC_LAT.GetValue())

    self.p.freq = float(self.TC_FREQ.GetValue())
    self.p.tcal = float(self.TC_TCAL.GetValue())

    self.p.beam = float(self.TC_BEAM.GetValue())
    self.p.doff = float(self.TC_DOFF.GetValue())
    self.p.dwidth = float(self.TC_DWIDTH.GetValue())/100
    self.p.rerr = float(self.TC_ERR.GetValue())/100
    self.p.dsymm = float(self.TC_DSYMM.GetValue())/100
    self.p.davg = float(self.TC_DAVG.GetValue())/100

    event.Skip()


class MPL_Base(wx.Panel):

  def __init__(self, parent):
    wx.Panel.__init__(self,parent=parent, id=-1)
    self.p = parent
    self.Figure = \
    matplotlib.figure.Figure(subplotpars= \
    matplotlib.figure.SubplotParams(0.06, 0.05, 0.98, 0.93, None, None))

    self.FigureCanvas = FigureCanvas(self, -1, self.Figure)
    self.ToolBar = NavigationToolbar(self.FigureCanvas)
    TopBoxSizer = wx.BoxSizer(wx.VERTICAL)
    TopBoxSizer.Add(self.FigureCanvas, proportion =1, flag = wx.EXPAND)
    TopBoxSizer.Add(self.ToolBar, proportion =0, flag = wx.EXPAND)
    self.SetSizer(TopBoxSizer)

  def UpdatePlot(self):
    self.FigureCanvas.draw()


  def gsplot(self):

    ffit = open(os.path.join(self.p.fit_path, self.p.scannum+'.fit'))
    trfit = ffit.read()[:-1].split('\n')
    ffit.close()
    trfit = [e for e in trfit if not e.startswith('#')]
    try:
      srcname = trfit[0].split()[0]
    except IndexError:
      print('No data found in FITS file !!!')
      return

    fdat = open(os.path.join(self.p.fit_path, self.p.scannum+'.dat'))
    trdat = fdat.readlines()
    fdat.close()
    trdat = [e[:-1] for e in trdat if (e[0]!='#' or e[:2]=='#!' or (e[0]=='#' and
      len(e.split())==6))]
    trdat = trdat[1:]
    tr = [e for e in trdat if re.match('#! SUBSCAN', e)]
    keys = [trdat.index(e) for e in tr]
    keys.append(len(trdat))

    subs = len(trfit)

    amp, eamp = np.zeros(subs), np.zeros(subs)
    scandir = np.array([])

    n1 = int((subs+1)**0.5)
    n2 = n1*1
    while n1*n2 < subs+1:
      n2 = n2 +1

    self.Figure.clf()
    gsp = gridspec.GridSpec(n1, n2)

    sub_flag = np.array([True]*subs)
    for i in range(subs):
      ax = self.Figure.add_subplot(gsp[i])

      trf = trfit[i].split()
      scandir = np.append(scandir, trf[2])
      mjd = float(trf[4])
      tsys = float(trf[14])
      amp[i], eamp[i] = float(trf[7]), float(trf[8])
      off = float(trf[9])
      hpbw = float(trf[11])
      az, el, sunang = float(trf[5]), float(trf[6]), float(trf[16])
      del trf

      trd = trdat[keys[i]+2:keys[i+1]]
      trd = [e.split()[0:3][::2] for e in trd]
      trd = np.array(trd).transpose()

      dataX = trd[0].astype(np.float)
      dataY = trd[1].astype(np.float)
      del trd

      if amp[i] != 0:
        fitY = ML.Gauss(1, 1).fGauss(np.array([amp[i],off,hpbw,0,0]), dataX)
      else:
        fitY = dataY*1.0


      ampcolor, offcolor, hpbwcolor = 'white', 'white', 'white'
      titlecolor, titlebackcolor = 'black', '#00FF00'

      xmin, xmax = min(dataX), max(dataX)
      dy = max(dataY) - min(dataY)
      if self.p.flag == False:
        ymin = max(min(dataY), -0.1*dy)
        ymax = min(max(dataY)+0.1*dy, 1.5*amp[i])

      else:
        if (amp[i] <= 0) or (eamp[i]/amp[i] > self.p.rerr):
          ampcolor = 'red'
        elif eamp[i]/amp[i] > 0.5*self.p.rerr:
          ampcolor = 'yellow'
        if np.fabs(off) > self.p.doff:
          offcolor = 'red'
        elif np.fabs(off) > 0.5*self.p.doff:
          offcolor = 'yellow'
        if np.fabs(hpbw/self.p.beam-1) > self.p.dwidth:
          hpbwcolor = 'red'
        elif np.fabs(hpbw/self.p.beam-1) > 0.5*self.p.dwidth:
          hpbwcolor = 'yellow'

      if 'red' in [ampcolor, offcolor, hpbwcolor]:
        sub_flag[i] = False

        ymin = min(dataY)
        ymax = max(dataY) +0.1*dy

      else:
        ymin = max(min(dataY), -0.1*dy)
        ymax = min(max(dataY)+0.1*dy, 1.5*amp[i])

      ax.text(0.05*xmax+0.95*xmin, 0.9*ymax+0.1*ymin, \
          'Amp=%0.4f' %amp[i], \
          fontsize=14, \
          horizontalalignment='left',\
          verticalalignment='top',\
          linespacing=1.5,\
          color=ampcolor)

      ax.text(0.05*xmax+0.95*xmin, 0.815*ymax+0.185*ymin, \
          'Off=%0.2f' %off, \
          fontsize=14, \
          horizontalalignment='left',\
          verticalalignment='top',\
          linespacing=1.5,\
          color=offcolor)

      ax.text(0.05*xmax+0.95*xmin, 0.73*ymax+0.27*ymin, \
          'HPBW=%0.2f' %hpbw, \
          fontsize=14, \
          horizontalalignment='left',\
          verticalalignment='top',\
          linespacing=1.5,\
          color=hpbwcolor)

      ax.text(0.95*xmax+0.05*xmin, 0.9*ymax+0.1*ymin, \
          'Az=%0.1f\nEl=%0.1f\nSun=%0.1f' %(az, el, sunang), \
          fontsize=14, \
          horizontalalignment='right',\
          verticalalignment='top',\
          linespacing=1.5,\
          color='white')

      ax.set_xlim(xmin, xmax)
      ax.set_ylim(ymin, ymax)
      ax.set_xlabel('%s-OFF (arcsec.)' %scandir[i])
      ax.set_ylabel('COUNTS')
      ax.plot(dataX, dataY, color='white', linewidth=1)
      if amp[i] != 0:
        ax.plot(dataX, fitY, color='#00FF00', linewidth=2)
      ax.grid(True)


    ax = self.Figure.add_subplot(gsp[subs])

    _amp0 = amp+eamp
    _amp1 = amp-eamp
    dy = max(_amp0) - min(_amp1)
    ymin = 1.1*min(_amp1) - 0.1*max(_amp0)
    ymax = 1.1*max(_amp0) - 0.1*min(_amp1)
    ax.set_xlim(0, subs+1)
    ax.set_ylim(ymin, ymax)
    del _amp0, _amp1

    ax.set_xlabel('Comparison of Subscans')

    if self.p.flag == True:
      if True not in sub_flag[scandir=='ALON'] or \
          True not in sub_flag[scandir=='ALAT']:
          titlecolor, titlebackcolor = 'white', 'red'

      inv_flag = np.array([not e for e in sub_flag])

      ax.errorbar(np.arange(1,subs+1)[sub_flag], \
          amp[sub_flag], yerr=eamp[sub_flag], fmt='d', \
          ecolor='white', mfc='white', mec='white', \
          ms=8, lw=2)

      ax.errorbar(np.arange(1,subs+1)[inv_flag], \
          amp[inv_flag], yerr=eamp[inv_flag], fmt='o', \
          ecolor='red', mfc='red', mec='red', \
          ms=10, lw=2)
    else:
      ax.errorbar(np.arange(subs)+1, amp, yerr=eamp, fmt='d',\
        ecolor='white', mfc='white', mec='white', \
        ms=8, lw=2)

    ax.plot([0, subs+1], [np.mean(amp)]*2, color='white', linewidth=2)

    self.Figure.suptitle('\nScan: %4s     %s     @%0.1fMHz' \
        %(self.p.scannum, srcname, self.p.freq), \
        weight='bold', ha='center', \
        linespacing=0.3, fontsize=20, \
        backgroundcolor=titlebackcolor,\
        color=titlecolor)

    self.UpdatePlot()
    del ax


class MyFrame(wx.Frame):
  def __init__(self, *args, **kwds):

    kwds['style'] = wx.DEFAULT_FRAME_STYLE
    wx.Frame.__init__(self, *args, **kwds)
    self.SysColor = wx.SystemSettings.GetColour(wx.SYS_COLOUR_BACKGROUND)

    self.CB_SCAN = wx.CheckBox(self, wx.ID_ANY, 'Scan')
    self.CB_DATE = wx.CheckBox(self, wx.ID_ANY, 'Date')
    self.TC_SCAN0 = wx.TextCtrl(self, -1, '0')
    self.TC_SCAN1 = wx.TextCtrl(self, -1, '9999')
    self.TC_Chan = wx.TextCtrl(self, -1, '0.5*(R+L)')
    self.TC_Spike = wx.TextCtrl(self, -1, '0')
    self.DPC0 = wx.DatePickerCtrl(self, -1)
    self.DPC1 = wx.DatePickerCtrl(self, -1)

    self.LB_SCANS = wx.ListBox(self, -1, choices=[], \
        style = wx.LB_SINGLE | wx.LB_HSCROLL)
    self.TC_FITS = wx.TextCtrl(self, -1, os.environ['HOME'])
    self.TC_FIT = wx.TextCtrl(self, -1, os.environ['HOME']+'/fit_tmp/')

    self.BT_OBSMODE = wx_bt(self, 'Observation Mode OFF', \
        self.MODE_ONOFF, style=wx.BU_EXACTFIT)
    self.BT_AFIT = wx_bt(self, 'Auto Fitting OFF', \
        self.FIT_ONOFF, style=wx.BU_EXACTFIT)
    self.BT_EXIT = wx_bt(self, 'EXIT', self.ON_EXIT, \
        style=wx.BU_EXACTFIT)

    self.GSPLOT = MPL_Base(self)

    self.__set_properties()
    self.__do_layout()

    self.Bind(wx.EVT_CHECKBOX, self.CHK_SCAN_DATE)
    self.Bind(wx.EVT_CLOSE, self.OnClose)
    self.Bind(wx.EVT_LISTBOX_DCLICK, self.DCLICK_SCAN, self.LB_SCANS)

    # load scan to the list
    ssu = Static_Scan_Update(self)
    ssu.setDaemon(True)
    ssu.start()

    pub.subscribe(self.GUI_Update, 'updating')


  def __set_properties(self):

    self.SetTitle('Gaussian Fitting for Cross-Scan Observation')
    self.SetSize((1280, 720))
    #self.SetSize((1920, 1080))
    self.CB_SCAN.SetValue(False)
    self.CB_DATE.SetValue(False)
    self.LB_SCANS.SetMinSize((340,100))
    self.FS = wx.SystemSettings_GetFont(wx.SYS_OEM_FIXED_FONT).GetPointSize()
    self.Font = wx.Font(self.FS, wx.MODERN, wx.NORMAL, wx.NORMAL, False)
    self.LB_SCANS.SetFont(self.Font)
    self.BT_EXIT.SetBackgroundColour(wx.Colour(255, 0, 0))

    self.fits_path = self.TC_FITS.GetValue()
    self.fit_path = self.TC_FIT.GetValue()

    self.tele_name = 'NSRT'
    self.diam = 26.0
    self.obs_lon = 87.178108
    self.obs_lat = 43.4709389

    self.freq = 4800.0
    self.tcal = 1.0

    self.flag = True
    self.beam = 1.12*3600*300*180/np.pi/self.freq/self.diam
    self.doff = self.beam*0.15
    self.dwidth = 0.1
    self.rerr = 0.1
    self.dsymm = 0.1
    self.davg = 0.15


  def __do_layout(self):

    SZ0 = wx.BoxSizer(wx.VERTICAL)
    SZ1 = wx.BoxSizer(wx.HORIZONTAL)
    SZ_LEFT = wx.BoxSizer(wx.VERTICAL)
    SZ_PLOT = wx.BoxSizer(wx.VERTICAL)

    SZ0.Add(SZ1, 1, wx_lpt('l|r|b|e'), 5)
    SZ1.Add(SZ_LEFT, 0, wx_lpt('l|r|b|e'), 1)
    SZ1.Add(SZ_PLOT, 4, wx_lpt('e'), 0)

    SZ_LISTSCAN = wx_hsbs(self)
    SZ_FILTER = wx_hsbs(self)
    SZ_LOCATION = wx_hsbs(self)
    SZ_SET = wx_vsbs(self)
    SZ_BUTTON = wx_hsbs(self)

    ########
    SZ_LEFT.Add(SZ_LISTSCAN, 1, wx_lpt('l|r|e'), 2)
    SZ_LEFT.Add(SZ_FILTER, 0, wx_lpt('l|r|e'), 2)
    SZ_LEFT.Add(SZ_LOCATION, 0, wx_lpt('l|r|e'), 2)
    SZ_LEFT.Add(SZ_SET, 0, wx_lpt('l|r|e'), 2)
    SZ_LEFT.Add(SZ_BUTTON, 0, wx_lpt('l|r|e'), 2)

    SZ_PLOT.Add(self.GSPLOT, 1, wx.EXPAND | wx.ALL, 5)

    ########
    SZ_LISTSCAN.Add(self.LB_SCANS, 1, wx_lpt('l|r|b|e|ach|acv'), 5)

    ########
    SZ_FILTER_0 = wx.BoxSizer(wx.VERTICAL)
    SZ_FILTER_1 = wx.BoxSizer(wx.VERTICAL)
    SZ_FILTER_2 = wx.BoxSizer(wx.VERTICAL)
    #SZ_FILTER_3 = wx.BoxSizer(wx.VERTICAL)

    SZ_FILTER_0.Add(self.CB_SCAN, 1, wx_lpt('l|r|ar|acv|e'), 5)
    SZ_FILTER_0.Add(self.CB_DATE, 1, wx_lpt('l|r|t|ar|acv|e'), 5)
    SZ_FILTER_1.Add(self.TC_SCAN0, 1, wx_lpt('b|e|ach|acv'), 5)
    SZ_FILTER_1.Add(self.DPC0, 0, wx_lpt('b|e|ach|acv'), 5)
    SZ_FILTER_2.Add(self.TC_SCAN1, 1, wx_lpt('r|b|e|ach|acv'), 5)
    SZ_FILTER_2.Add(self.DPC1, 0, wx_lpt('r|b|e|ach|acv'), 5)
    #SZ_FILTER_3.Add(self.BT_SCAN, 0, wx_lpt('l|r|b|e|ach|acv'), 5)
    #SZ_FILTER_3.Add(self.BT_DATE, 0, wx_lpt('l|r|b|e|ach|acv'), 5)
    SZ_FILTER.Add(SZ_FILTER_0, 1, wx_lpt('b|e|acv'), 5)
    SZ_FILTER.Add(SZ_FILTER_1, 2, wx_lpt('l|r|e'), 5)
    SZ_FILTER.Add(SZ_FILTER_2, 2, wx_lpt('l|r|e'), 5)
    #SZ_FILTER.Add(SZ_FILTER_3, 1, wx.EXPAND, 0)

    ########
    SZ_LOC_0 = wx.BoxSizer(wx.VERTICAL)
    SZ_LOC_1 = wx.BoxSizer(wx.VERTICAL)
    SZ_LOC_2 = wx.BoxSizer(wx.VERTICAL)

    SZ_LOC_0.Add(wx_st(self, 'FITS'), 1, wx_lpt('l|r|ar|acv'), 5)
    SZ_LOC_0.Add(wx_st(self, 'fit'), 1, wx_lpt('l|r|t|ar|acv'), 5)
    SZ_LOC_1.Add(self.TC_FITS, 1, wx_lpt('b|e|ach|acv'), 5)
    SZ_LOC_1.Add(self.TC_FIT, 1, wx_lpt('b|e|ach|acv'), 5)
    SZ_LOC_2.Add(wx_bt(self, 'Browse', self.BROWSE_FITS, \
        style=wx.BU_EXACTFIT), 0, wx_lpt('l|r|b|e|ach|acv'), 5)
    SZ_LOC_2.Add(wx_bt(self, 'Browse', self.BROWSE_FIT, \
        style=wx.BU_EXACTFIT), 0, wx_lpt('l|r|b|e|ach|acv'), 5)

    SZ_LOCATION.Add(SZ_LOC_0, 0, wx_lpt('t|e'), 5)
    SZ_LOCATION.Add(SZ_LOC_1, 1, wx_lpt('l|r|e'), 5)
    SZ_LOCATION.Add(SZ_LOC_2, 0, wx.EXPAND, 0)

    ########
    SZ_SET_0 = wx.BoxSizer(wx.HORIZONTAL)
    SZ_SET_1 = wx.BoxSizer(wx.HORIZONTAL)

    SZ_SET_0.Add(wx_st(self, 'DeNoise'), 0, wx_lpt('l|r|ar|acv'), 5)
    SZ_SET_0.Add(self.TC_Spike, 0, wx_lpt('b|e|ar|acv'), 5)
    SZ_SET_0.Add(wx_st(self, 'Channel'), 0, wx_lpt('l|r|ar|acv'), 5)
    SZ_SET_0.Add(self.TC_Chan, 1, wx_lpt('b|r|e|ar|acv'), 5)

    SZ_SET_1.Add(self.BT_OBSMODE, 1, wx_lpt('l|r|b|ach|acv'), 5)
    SZ_SET_1.Add(self.BT_AFIT, 1, wx_lpt('l|r|b|ach|acv'), 5)

    SZ_SET.Add(SZ_SET_0, 1, wx.EXPAND, 0)
    SZ_SET.Add(SZ_SET_1, 0, wx.EXPAND, 0)

    ########
    SZ_BUTTON.Add(wx_bt(self, 'Option', self.OPEN_OPTION, \
        style=wx.BU_EXACTFIT), 1, wx_lpt('l|r|b|ach|acv'), 5)
    SZ_BUTTON.Add(wx_bt(self, 'Next', self.FIT_NEXT, \
        style=wx.BU_EXACTFIT), 1, wx_lpt('l|r|b|ach|acv'), 5)
    SZ_BUTTON.Add(self.BT_EXIT, 1, wx_lpt('l|r|b|ach|acv'), 5)

    self.SetSizer(SZ0)
    self.Layout()


  def GUI_Update(self, msg):

    if msg[0] == 'Static_Scan_Update':
      idx, content = msg[1:3]
      self.LB_SCANS.InsertItems([content], idx)

    elif msg[0] == 'Dynamic_Scan_Update':
      idx, content = msg[1:3]
      self.LB_SCANS.InsertItems([content], idx)

    elif msg[0] == 'Auto_Gaussian_Fitting':
      idx = msg[1]
      self.LB_SCANS.Select(idx)
      self.scannum = self.scans[idx]

      self.fit_data()
      self.GSPLOT.gsplot()


  def OnClose(self, event):
    event.Skip()


  def CHK_SCAN_DATE(self, event):

    self.LB_SCANS.SetItems([])
    ssu = Static_Scan_Update(self)
    ssu.setDaemon(True)
    ssu.start()

    event.Skip()

  def MODE_ONOFF(self, event):

    if self.BT_OBSMODE.BackgroundColour == self.SysColor:
      self.BT_OBSMODE.SetBackgroundColour(wx.Colour(0, 255, 0))
      self.BT_OBSMODE.SetLabel('Observation Mode ON')

      fits_path = os.path.join(os.environ['HOME'], 'Raw')
      self.TC_FITS.SetValue(fits_path)
      self.fits_path = fits_path

      self.LB_SCANS.SetItems([])
      self.scans = []
      #ssu = Static_Scan_Update(self)
      #ssu.setDaemon(True)
      #ssu.start()

      dsu = Dynamic_Scan_Update(self)
      dsu.setDaemon(True)
      dsu.start()

    else:
      self.BT_OBSMODE.SetBackgroundColour(self.SysColor)
      self.BT_OBSMODE.SetLabel('Observation Mode OFF')

    event.Skip()


  def FIT_ONOFF(self, event):

    if self.BT_AFIT.BackgroundColour == self.SysColor:
      self.BT_AFIT.SetBackgroundColour(wx.Colour(0, 255, 0))
      self.BT_AFIT.SetLabel('Auto Fitting ON')

      agf = Auto_Gaussian_Fitting(self)
      agf.setDaemon(True)
      agf.start()

    else:
      self.BT_AFIT.SetBackgroundColour(self.SysColor)
      self.BT_AFIT.SetLabel('Auto Fitting OFF')


    event.Skip()


  def BROWSE_FITS(self, event):

    if self.BT_OBSMODE.BackgroundColour == self.SysColor:

      dialog = wx.DirDialog(self, 'Set Raw FITS Location', \
          defaultPath = os.environ['HOME'], \
          style = wx.DD_DEFAULT_STYLE)
      if dialog.ShowModal() == wx.ID_OK:
        dialog.Destroy()

      self.fits_path = dialog.GetPath()
      self.TC_FITS.SetValue(self.fits_path)

      self.LB_SCANS.SetItems([])
      ssu = Static_Scan_Update(self)
      ssu.setDaemon(True)
      ssu.start()

    event.Skip()


  def BROWSE_FIT(self, event):

    dialog = wx.DirDialog(self, 'Set fit file Location', \
        defaultPath = os.environ['HOME'], \
        style = wx.DD_DEFAULT_STYLE)
    if dialog.ShowModal() == wx.ID_OK:
      self.fit_path = dialog.GetPath()
      self.TC_FIT.SetValue(self.fit_path)
      dialog.Destroy()

    event.Skip()


  def OPEN_OPTION(self, event):

    dlg = OptDialog(self, -1)
    val = dlg.ShowModal()

    dlg.CenterOnScreen()

    dlg.Destroy()


  def DCLICK_SCAN(self, event):
    idx = self.LB_SCANS.Selection
    self.scannum = self.scans[idx]
    self.fit_data()
    self.GSPLOT.gsplot()

    event.Skip()


  def ON_EXIT(self, event):

    self.Destroy()

    event.Skip()


  def FIT_NEXT(self, event):
    idx = self.LB_SCANS.Selection
    self.LB_SCANS.Select(idx+1)
    self.scannum = self.scans[idx+1]
    self.fit_data()
    self.GSPLOT.gsplot()

    event.Skip()


  def fit_data(self):

    fits_name = self.scannum+'.FITS'

    kwds = {'chan': self.TC_Chan.Value,
        'despike': float(self.TC_Spike.Value)}
    try:
      rd = ObsData.RawData(self, fits_name, **kwds)
      rd.fit_data()
      del fits_name, kwds, rd
      return 1

    except:
      return 0


class Static_Scan_Update(Thread):

  def __init__(self, parent):
    Thread.__init__(self)
    self.__name__ = 'Static_Scan_Update'
    self.setName(self.__name__)
    self.p = parent

  def run(self):

    scans = []
    self.p.scans = load_scans(self.p.fits_path)
    N = len(self.p.scans)
    if N == 0:
      return

    if self.p.CB_SCAN.Value == True:
      tmp = np.array(self.p.scans)
      scan0 = '%04d' %int(self.p.TC_SCAN0.Value)
      scan1 = '%04d' %int(self.p.TC_SCAN1.Value)

      flt = (tmp >= scan0)*(tmp <= scan1)
      tmp = tmp[flt]
      self.p.scans = list(tmp)
      del scan0, scan1, flt, tmp

    i = 0
    if self.p.CB_DATE.Value == True:
      date0 = wxDateFormatter(self.p.DPC0.Value)
      date1 = wxDateFormatter(self.p.DPC1.Value)

      for scan in self.p.scans:
        fname = os.path.join(self.p.fits_path, '%0s.FITS' %scan)
        content = load_fits_property(fname)
        scan_date = content.split()[3]
        if scan_date >= date0 and scan_date <= date1 and content != None:
          scans.append(content.split()[0])
          wx.CallAfter(pub.sendMessage, 'updating', \
            msg = (self.__name__, i, content))
          i+=1

    else:
      for scan in self.p.scans:
        fname = os.path.join(self.p.fits_path, '%4s.FITS' %scan)
        content = load_fits_property(fname)
        if content != None:
          scans.append(content.split()[0])
          wx.CallAfter(pub.sendMessage, 'updating', \
              msg = (self.__name__, i, content))
          i+=1

    self.p.scans = scans
    #if len(scans) > 0:
    #    self.p.scannum = self.p.scans[0]


class Dynamic_Scan_Update(Thread):

  def __init__(self, parent):
    Thread.__init__(self)
    self.__name__  = 'Dynamic_Scan_Update'
    self.setName(self.__name__)
    self.p = parent
    #self.start()

  def run(self):

    while (self.p.BT_OBSMODE.BackgroundColour==wx.Colour(0, 255, 0)):
      scans = load_scans(self.p.fits_path)

      if self.p.CB_SCAN.Value == True:
        tmp = np.array(scans)
        scan0 = '%04d' %int(self.p.TC_SCAN0.Value)
        scan1 = '%04d' %int(self.p.TC_SCAN1.Value)

        flt = (tmp >= scan0)*(tmp <= scan1)
        tmp = tmp[flt]
        scans = list(tmp)
        del scan0, scan1, flt, tmp

      i = len(self.p.scans)
      if self.p.CB_DATE.Value == True:
        date0 = wxDateFormatter(self.p.DPC0.Value)
        date1 = wxDateFormatter(self.p.DPC1.Value)

      else:
        date0 = '0000-00-00'
        date1 = '9999-99-99'

      for scan in scans:
        if scan not in self.p.scans:
          fname = os.path.join(self.p.fits_path, '%0s.FITS' %scan)
          content = load_fits_property(fname)
          scan_date = content.split()[3]
          if scan_date >= date0 and scan_date <= date1 \
            and content != None:
            self.p.scans.append(scan)
            wx.CallAfter(pub.sendMessage, 'updating', \
              msg = (self.__name__, i, content))
            i+=1

      time.sleep(5)


class Auto_Gaussian_Fitting(Thread):

  def __init__(self, parent):
    Thread.__init__(self)
    self.__name__ = 'Auto_Gaussian_Fitting'
    self.setName(self.__name__)
    self.p = parent
    #self.start()

  def run(self):

    try: idx = self.p.LB_SCANS.Selection
    except: idx = 0
    if idx == -1: idx = 0

    while self.p.BT_AFIT.BackgroundColour==wx.Colour(0, 255, 0):
      if len(self.p.scans) == idx + 1:
        next_fsize = os.path.getsize(os.path.join(self.p.fits_path,
          self.p.scans[idx]+'.FITS'))

        if next_fsize <= 23040.0:
          wx.CallAfter(pub.sendMessage, 'updating', \
            msg = (self.__name__, idx-1))
          time.sleep(15)

        else:
          wx.CallAfter(pub.sendMessage, 'updating', \
            msg = (self.__name__, idx))
          idx += 1
          time.sleep(3)

      elif len(self.p.scans) == idx:
        wx.CallAfter(pub.sendMessage, 'updating', \
          msg = (self.__name__, idx-1))
        time.sleep(10)

      else:
        wx.CallAfter(pub.sendMessage, 'updating', \
          msg = (self.__name__, idx))
        idx += 1
        time.sleep(3)


def load_scans(pth):

  scans = os.listdir(pth)
  scans = [e[0:4] for e in scans if e.endswith('.FITS')]
  scans.sort()

  return scans


def load_fits_property(fname):

  p = {'Jan':1,
      'Feb':2,
      'Mar':3,
      'Apr':4,
      'May':5,
      'Jun':6,
      'Jul':7,
      'Aug':8,
      'Sep':9,
      'Oct':10,
      'Nov':11,
      'Dec':12}

  try:
    fits = pyfits.open(fname)
    #if len(fits) == 3: return None
    #else:
    fits = pyfits.open(fname)[0]
    src = fits.header['object']
    date = fits.header['date_obs'].split()
    time = date[3]
    date = '%04d-%02d-%02d' %(int(date[4]), p[date[1]], int(date[2]))
    scan = os.path.split(fname)[1][0:4]

    out_text = '%4s %-10s %-9s %s' %(scan, src, time, date)

    return out_text

  except IOError:

    return None


def wxDateFormatter(wxdate):
  return '%04d-%02d-%02d' %(wxdate.Year, wxdate.Month+1, wxdate.Day)

def wx_st(parent, label):
  return wx.StaticText(parent, -1, label)


def wx_bt(parent, label, handler, **kwds):

  button = wx.Button(parent, -1, label, **kwds)
  parent.Bind(wx.EVT_BUTTON, handler, button)

  return button


def wx_hsbs(parent, **kwds):

  return wx.StaticBoxSizer(wx.StaticBox(parent, -1, ''), wx.HORIZONTAL)


def wx_vsbs(parent, **kwds):

  return wx.StaticBoxSizer(wx.StaticBox(parent, -1, ''), wx.VERTICAL)


def wx_lpt(line):
  # layout properties
  p = {'a' : '240',
    'at': '0',
    'v' : '8',
    'h' : '4',
    'b': '128',
    't' : '64',
    'e': '8192',
    'l' : '16',
    'r' : '32',
    'al' : '0',
    'ar' : '512',
    'ach' : '256',
    'acv' : '2048'}
  line = line.split('|')
  line = [p[e] for e in line]
  line = ('|').join(line)

  return eval(line)


def main():

  class App(wx.App):

    def OnInit(self):
      #self.frame = MyFrame(None, -1)
      #self.frame.Show()
      #self.SetTopWindow(self.frame)
      MyFrame(None, -1).Show()
      return True

  app = App()
  app.MainLoop()


if __name__ == '__main__':

  main()
