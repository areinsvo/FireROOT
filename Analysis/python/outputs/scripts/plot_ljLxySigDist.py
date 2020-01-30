#!/usr/bin/env python
from __future__ import print_function
import os
from FireROOT.Analysis.Utils import *
from rootpy.plotting import Hist, HistStack, Legend, Canvas
from rootpy.io import root_open

inname = os.path.join(os.getenv('CMSSW_BASE'), 'src/FireROOT/Analysis/python/outputs/rootfiles/ljLxySigDist.root')
outdir = os.path.join(os.getenv('CMSSW_BASE'), 'src/FireROOT/Analysis/python/outputs/plots/ljLxySigDist')
if not os.path.isdir(outdir): os.makedirs(outdir)

from rootpy.plotting.style import set_style
set_style(MyStyle())

c = Canvas()

f = root_open(inname)

samples = 'mXX-150_mA-0p25_lxy-300|mXX-500_mA-1p2_lxy-300|mXX-800_mA-5_lxy-300'.split('|')
samples.extend( 'mXX-100_mA-5_lxy-0p3|mXX-1000_mA-0p25_lxy-0p3'.split('|') )

histCollection = [
    {
        'name': 'lxysig',
        'binning': (20, 0, 20),
        'title': 'muon-type lepton-jet lxy significance;L_{xy}/#sigma_{L_{lxy}};counts/1',
    },
    {
        'name': 'lxysig1',
        'binning': (20, 0, 20),
        'title': 'muon-type lepton-jet lxy significance(min d0Sig>1);L_{xy}/#sigma_{L_{lxy}};counts/1',
    },
    {
        'name': 'lxysig2',
        'binning': (20, 0, 20),
        'title': 'muon-type lepton-jet lxy significance(min d0Sig>1, minIso<0.03);L_{xy}/#sigma_{L_{lxy}};counts/1',
    },
    {
        'name': 'lxysig3',
        'binning': (20, 0, 20),
        'title': 'muon-type lepton-jet lxy significance(min d0Sig>1, minIso<0.1, #Delta#phi>#pi/2);L_{xy}/#sigma_{L_{lxy}};counts/1',
    },
]



for hinfo in histCollection:
    for chan in ['2mu2e', '4mu']:
        hs = []

        histName = 'bkgs__{}__{}'.format(chan, hinfo['name'])
        hstack = getattr(f, histName)
        for h in hstack:
            h.fillstyle='solid'
            h.linewidth=0
            h.legendstyle='F'
            h.fillcolor = bkgCOLORS[h.title]
            h.xaxis.SetRange(1, h.nbins()+1)
        stackError = ErrorBandFromHistStack(hstack)

        hs.append(hstack)
        hs.append(stackError)

        legend = Legend(len(samples)+hstack.GetNhists()+1, pad=c, margin=0.1, entryheight=0.02, textsize=12)

        for h in hstack: legend.AddEntry(h)
        legend.AddEntry(stackError)

        for i, s in enumerate(samples):
            histName = '{}__{}__{}'.format(s, chan, hinfo['name'])
            h = getattr(f, histName).Clone()
            h.xaxis.SetRange(1, h.nbins()+1)
            h.title = s#+' (norm.)'
            # h.Scale(1./h.Integral())
            h.drawstyle = 'hist pmc plc'
            h.legendstyle = 'L'
            h.linewidth = 2
            hs.append(h)
            legend.AddEntry(h)

        xmin, xmax, ymin, ymax = get_limits(hs)
        draw(hs, pad=c, xlimits=(xmin, hs[-1].xaxis.GetBinUpEdge(hs[-1].nbins()+1)), logy=True)
        legend.Draw()
        title = TitleAsLatex('[{}] {}'.format(chan.replace('mu', '#mu'), hinfo['title'].split(';')[0]))
        title.Draw()
        draw_labels('59.74 fb^{-1} (13 TeV)', cms_position='left', extra_text='work-in-progress')

        c.SaveAs('{}/{}_{}.pdf'.format(outdir, hinfo['name'], chan))
        c.Clear()

f.Close()