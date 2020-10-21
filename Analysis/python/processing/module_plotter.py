#!/usr/bin/env python
from __future__ import print_function
import argparse
import os, sys, math
from functools import reduce
import ROOT
from FireROOT.Analysis.Utils import *
from FireROOT.Analysis.samples.signalnumbers import genxsec
from rootpy.plotting import Hist, Legend, Canvas, HistStack, Pad
from rootpy.plotting.shapes import Line
from rootpy.io import root_open

## parser
parser = argparse.ArgumentParser(description="module plotter.")
parser.add_argument("--inname", "-i", type=str, default=None, help='input ROOT file name')
parser.add_argument("--normsig", "-r", type=float, default=-1, help='Normalize signal distributions to a fix value')
parser.add_argument("--normsigxsec", "-x", type=float, default=59.74, help='Normalize signal distributions to a fix xsec [fb]')
parser.add_argument("--dataset", "-d", type=str, default='mc', choices=['data', 'mc', 'all'], help='dataset to plot')
parser.add_argument("--mbase", "-b", type=str, default='modules', choices=['modules', 'proxy',], help='module base name')
parser.add_argument("--histname", "-t", type=str, default=None, help='only plot this histogram')
parser.add_argument("--logx", action='store_true')
parser.add_argument("--logy", action='store_true')
parser.add_argument("--overflow", type=bool, default=True)
parser.add_argument("--underflow", type=bool, default=True)
parser.add_argument("--xmin", type=float, default=None, help='min X')
parser.add_argument("--xmax", type=float, default=None, help='max X')
parser.add_argument("--ymin", type=float, default=None, help='min Y')
parser.add_argument("--ymax", type=float, default=None, help='max Y')
parser.add_argument("--xdiv", type=int, default=None, help='X ndivisions')
parser.add_argument("--ydiv", type=int, default=None, help='Y ndivisions')
parser.add_argument("--vline", "-lx", type=float, default=None, help="draw a vline.")

args = parser.parse_args()
if args.normsig>0 and args.normsigxsec>0:
    sys.exit('You cannot normalize signal to a fixed value and a fixed xsec at the same time.')

inputdir = os.path.join(os.getenv('CMSSW_BASE'), 'src/FireROOT/Analysis/python/outputs/rootfiles/{}/'.format(args.mbase))
assert(os.path.isdir(inputdir))
allfiles = [f for f in os.listdir(inputdir) if f.endswith('.root')]
if args.inname+'.root' not in allfiles:
    sys.exit('Available ROOT files: {}'.format(str(allfiles)))

def get_unique_histnames(fname):
    bkgH, sigH, dataH = [], [], []
    f = root_open(fname)
    for dirpath, dirnames, filenames in f.walk():
        if not filenames: continue
        if dirpath.endswith('data') or dirpath.endswith('bkg'):
            dirobj = reduce(lambda a,b: getattr(a, b), [f]+dirpath.split('/'))
            for fn in filenames:
                h_ = getattr(dirobj, fn)
                if isinstance(h_, HistStack): h_ = h_[0]
                if h_.GetDimension()>1: continue
                if dirpath.endswith('data'): dataH.append(fn)
                if dirpath.endswith('bkg'): bkgH.append(fn)
        if 'sig' in dirpath and not dirnames:
            dirobj = reduce(lambda a,b: getattr(a, b), [f]+dirpath.split('/'))
            for fn in filenames:
                if getattr(dirobj, fn).GetDimension()>1: continue
                sigH.append(fn)
    f.close()
    if args.dataset=='data': return list(set(sigH)&set(dataH))
    elif args.dataset=='mc': return list(set(sigH)&set(bkgH)) if sigH and bkgH else list(set(sigH+bkgH))
    elif args.dataset=='all': return list(set(sigH)&set(bkgH)&set(dataH))



if __name__ ==  '__main__':

    infile = os.path.join(inputdir, args.inname+'.root')
    outdir = os.path.join(os.getenv('CMSSW_BASE'), 'src/FireROOT/Analysis/python/outputs/plots/{}/'.format(args.mbase))
    outdir = os.path.join(outdir, args.inname)
    if not os.path.isdir(outdir): os.makedirs(outdir)

    from rootpy.plotting.style import set_style
    set_style(MyStyle())

    hist_toplot = []
    if args.histname: hist_toplot = [args.histname]
    else: hist_toplot = get_unique_histnames(infile)

    c = Canvas()
    f = root_open(infile)
    for hname in hist_toplot:
        for chan in ['4mu', '2mu2e']:
            if not hasattr(f, 'ch'+chan): continue
            channelDir = getattr(f, 'ch'+chan)

            hs, legItems = [], []
            htitle = None

            if args.dataset=='mc' or args.dataset=='all':
                ## bkg
                if hasattr(channelDir, 'bkg') and hasattr(channelDir.bkg, hname):
                    hstack = getattr(channelDir.bkg, hname)
                    if not htitle: htitle = hstack.title

                    for h in hstack:
                        h.fillstyle='solid'
                        h.linewidth=0
                        h.legendstyle='F'
                        h.fillcolor = bkgCOLORS[h.title]
                        if args.overflow:
                            h = h.merge_bins([(-2, -1),])
                            # h.SetBinContent(h.nbins(), h.GetBinContent(h.nbins())+h.overflow())
                            # h.SetBinError(h.nbins(), math.sqrt(h[h.nbins()].value))
                        if args.underflow: h = h.merge_bins([(0, 1),])

                    hs.append(hstack)


                    if hstack.Integral():
                        stackError = ErrorBandFromHistStack(hstack)
                        hs.append(stackError)
                        for h in hstack:
                            if h.integral(overflow=True)==0: continue
                            legItems.append(h)
                        legItems.append(stackError)

            if args.dataset=='data' or args.dataset=='all':
                ## data
                if hasattr(channelDir, 'data') and hasattr(channelDir.data, hname):
                    hData = getattr(channelDir.data, hname)
                    if not htitle: htitle = hData.title
                    hData.title = 'data'
                    # hData.SetBinContent(hData.nbins(), hData.GetBinContent(hData.nbins())+hData.overflow())
                    # hData.SetBinError(hData.nbins(), math.sqrt(hData[hData.nbins()].value))
                    if args.overflow:  hData = hData.merge_bins([(-2, -1),])
                    if args.underflow: hData = hData.merge_bins([(0, 1),])
                    hData.legendstyle = 'LEP'
                    hs.append(hData)
                    legItems.append(hData)

            ## sig
            if hasattr(channelDir, 'sig'):
                sampleSig = 'mXX-150_mA-0p25_lxy-300|mXX-500_mA-1p2_lxy-300|mXX-800_mA-5_lxy-300'.split('|')
                sampleSig.extend( 'mXX-100_mA-5_lxy-0p3|mXX-1000_mA-0p25_lxy-0p3'.split('|') )

                for i, ds in enumerate(sampleSig):
                    if not hasattr(channelDir.sig, ds): continue
                    dsdir = getattr(channelDir.sig, ds)
                    h=None
                    for k in dsdir.keys():
                        if k.name!=hname: continue
                        h = getattr(dsdir, k.name).Clone()
                    if h is None: continue

                    if not htitle: htitle = h.title
                    
                    h.title = ds
                    if args.normsig>0:
                        h.title += ' (norm.)'
                        if h.integral(overflow=True)>0:
                            h.scale( 1.*args.normsig/h.integral(overflow=True) )
                    if args.normsigxsec>0 and args.normsigxsec!=59.74:
                        h.title += ' (norm. {:g}fb)'.format(args.normsigxsec)
                        mboundstate = int(ds.split('_')[0].replace('mXX-',''))
                        if h.integral(overflow=True)>0:
                            h.scale( args.normsigxsec/genxsec[mboundstate] )
                    if args.overflow:
#                        h = h.merge_bins([(-2, -1),])
                         h.SetBinContent(h.nbins(), h.GetBinContent(h.nbins())+h.overflow())
                         h.SetBinError(h.nbins(), math.sqrt(h[h.nbins()].value))
                    if args.underflow: 
                        h.SetBinContent(1, h.GetBinContent(1)+h.GetBinContent(0))
#                        h = h.merge_bins([(0, 1),])

                    h.drawstyle = 'hist'
                    h.color = sigCOLORS[i]
                    h.legendstyle = 'L'
                    h.linewidth = 2
                    hs.append(h)
                    legItems.append(h)


            if args.dataset=='all':
                # divide canvas to draw ratio
                mainPad = Pad(0, 0.25, 1, 1.)
                mainPad.SetBottomMargin(0.0)
                mainPad.Draw()
                subPad = Pad(0, 0.05, 1, 0.24)
                subPad.SetTopMargin(0.02)
                subPad.SetBottomMargin(0.25)
                subPad.Draw()
            else: mainPad = ROOT.gPad

            mainPad.cd()
            legend = Legend(legItems, pad=mainPad, margin=0.25, leftmargin=0.45, topmargin=0.02, entrysep=0.01, entryheight=0.02, textsize=10)

            xmin_, xmax_, ymin_, ymax_ = get_limits(hs, logx=args.logx, logy=args.logy,)
            if args.xmax is not None: xmax_ = args.xmax
            if args.xmin is not None: xmin_ = args.xmin
            if args.ymax is not None: ymax_ = args.ymax
            if args.ymin is not None: ymin_ = args.ymin


            axes, limits = draw(
                hs, pad=mainPad, logx=args.logx, logy=args.logy,
                xlimits=(xmin_, xmax_), ylimits=(ymin_, ymax_),
                xdivisions=args.xdiv, ydivisions=args.ydiv,
            )
            if htitle and 'phi' in htitle:
                decorate_axis_pi(axes[0])

            if args.logx: axes[0].SetMoreLogLabels()
            mainPad.SetGrid()
            # ROOT.gPad.Update()
            # ROOT.gPad.RedrawAxis('G')
            legend.Draw()
            title = TitleAsLatex('[{}] {}'.format(chan.replace('mu', '#mu'), htitle.split(';')[0]))
            title.Draw()
            draw_labels('59.74 fb^{-1} (13 TeV)', cms_position='left', extra_text='work-in-progress')
            if args.vline is not None:
                vline = Line(args.vline, mainPad.GetUymin(), args.vline, axes[1].GetXmax())
                vline.color='black'
                vline.linewidth=2
                vline.linestyle='dashed'
                vline.Draw()

            if args.dataset=='all':
                # draw ratio on subpad
                subPad.cd()

                _ratio = hData.clone()
                _ratio.Divide(sumHistStack(hstack))
                _ratio.SetMarkerSize(0.8)
                _ratio.yaxis.SetTitle('Data/MC')
                _ratio.yaxis.SetTitleOffset(0.35)
                _ratio.yaxis.CenterTitle()
                _ratio.yaxis.SetTitleSize(0.1)
                _ratio.yaxis.SetLabelSize(0.1)
                _ratio.xaxis.SetTitleSize(0.1)
                _ratio.xaxis.SetLabelSize(0.1)

                _mcerr = sumHistStack(hstack)
                _mcerr.Divide(sumHistStack(hstack))
                _mcerr.SetMarkerSize(0)
                _mcerr.SetFillColor(ROOT.kGray+2)
                _mcerr.SetFillStyle(3254)

                axis, limits = draw([_ratio,], pad=subPad, xdivisions=args.xdiv, ylimits=(0,2), ydivisions=5, )
                _mcerr.Draw('same e2')
                if htitle and 'phi' in htitle: decorate_axis_pi(_ratio.xaxis)
                ratio_line = Line(_ratio.xaxis.GetXmin(), 1, _ratio.xaxis.GetXmax(), 1)
                ratio_line.color = 'red'
                ratio_line.linewidth=2
                ratio_line.Draw()

            c.SaveAs('{}/{}__{}_{}.pdf'.format(outdir, args.dataset, chan, hname))
            c.Clear()

    f.close()
