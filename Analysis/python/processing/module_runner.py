#!/usr/bin/env python
from __future__ import print_function
import argparse
import math, numbers
import os, sys, fnmatch, json
from collections import OrderedDict
from multiprocessing import Pool

from tqdm import tqdm
from FireROOT.Analysis.Utils import *
from FireROOT.Analysis.DatasetMapLoader import (
    DatasetMapLoader,
    SigDatasetMapLoader,
    CentralSignalMapLoader,
)

from rootpy.io import root_open
from rootpy.logger import log
log = log[__name__]


## parser
parser = argparse.ArgumentParser(description="module runner.")
parser.add_argument("--dataset", "-d",  type=str, nargs='*', default='mc', help='dataset type')
parser.add_argument("--sigparam", "-p",  type=str, nargs='*',default=None, help='signal parameters')
parser.add_argument("--module", "-m", type=str, help='module path')
parser.add_argument("--outname", "-o", type=str, default=None, help='output ROOT file name')
parser.add_argument("--maxevents", "-n", type=int, default=-1, help='max number of events to run')
parser.add_argument("--channel", "-ch", nargs='*', default=['2mu2e', '4mu'], choices=['2mu2e', '4mu'], help='channels to run')
parser.add_argument("--create", "-c", type=str, default='recreate', choices=['recreate', 'update'], help='update output by')
parser.add_argument("--mbase", "-b", type=str, default='modules', choices=['modules', 'proxy', 'debug', 'centralSig'], help='module base name')
parser.add_argument("--update-signalsample", type=str, default=None, help='path to updated signal sample json file, will be update to signal DataMap')
parser.add_argument("--private", action='store_true', help='run private signal sample')
args = parser.parse_args()

def args_sanity(args):
    ## dataset
    bkg, sig, data = False, False, False
    if isinstance(args.dataset, str) or 'all' in args.dataset:
        bkg, sig, data = True, True, True
    else:
        if 'mc' in args.dataset: bkg, sig = True, True
        if 'bkg' in args.dataset: bkg = True
        if 'sig' in args.dataset: sig = True
        if 'data' in args.dataset: data = True

    ## module
    moduleBase = os.path.join(os.getenv('CMSSW_BASE'), 'src/FireROOT/Analysis/python/processing/{}'.format(args.mbase))
    assert(os.path.isdir(moduleBase))

    return bkg, sig, data

runbkg, runsig, rundata = args_sanity(args)
_modulebase = 'FireROOT.Analysis.processing.{}'.format(args.mbase)
imp = __import__('{}.{}'.format(_modulebase, args.module), fromlist=['MyEvents', 'histCollection'])



if __name__ == '__main__':

    outdir = os.path.join(os.getenv('CMSSW_BASE'), 'src/FireROOT/Analysis/python/outputs/rootfiles/{}/'.format(args.mbase))
    if not os.path.isdir(outdir): os.makedirs(outdir)

    if args.outname: outname = os.path.join(outdir, '{}.root'.format(args.outname))
    else:            outname = os.path.join(outdir, '{}.root'.format(args.module))
    if args.create == 'update' and not os.path.isfile(outname):
        sys.exit('UPDATE was used not file not recreated yet.')

    if args.mbase!='proxy' and (runbkg or rundata): dml = DatasetMapLoader()

    if runsig:
        if args.private: sdml = SigDatasetMapLoader()
        else:            sdml = CentralSignalMapLoader()


        sampleSig = 'mXX-150_mA-0p25_lxy-300|mXX-500_mA-1p2_lxy-300|mXX-800_mA-5_lxy-300'.split('|')
        sampleSig.extend( 'mXX-100_mA-5_lxy-0p3|mXX-1000_mA-0p25_lxy-0p3'.split('|') )
        if args.sigparam:
            print("Checking matching samples now", args.sigparam)
            sampleSig = []
            for s in args.sigparam:
                if '*' in s or '?' in s:
                    sampleSig.extend( fnmatch.filter(sdml.get_datasets('4mu').keys(), s) )
                else: sampleSig.append(s)
            sampleSig = list(set(sampleSig))
            print("Using these samples:",sampleSig)

        def dofill(pack):
            ds, files, scale, maxevents, channel = pack
#            events_ = imp.MyEvents(files=files, type='MC', maxevents=maxevents, channel=channel)
            events_ = imp.MyEvents(files=files, type='MC', dtag=ds, maxevents=maxevents, channel=channel)
            events_.setScale(scale)

            for chan in events_.channel:
                for hinfo in imp.histCollection:
                    if 'class' in hinfo:
                        events_.bookHisto('{}/{}'.format(chan, hinfo['name']),
                                hinfo['class'](*hinfo['binning'], name='{}__{}__{}'.format(ds, chan, hinfo['name']),
                                    title=hinfo['title'], drawstyle='hist e1', legendstyle='LEP')
                                )
                    else:
                        _dim = 0
                        for x in hinfo['binning']:
                            if isinstance(x, numbers.Number): _dim+=1/3.
                            elif isinstance(x, list): _dim+=1
                        if _dim==1:
                            events_.bookHisto('{}/{}'.format(chan, hinfo['name']),
                                    ROOT.Hist(*hinfo['binning'], name='{}__{}__{}'.format(ds, chan, hinfo['name']),
                                        title=hinfo['title'], drawstyle='hist', legendstyle='L')
                                    )
                        else:
                            events_.bookHisto('{}/{}'.format(chan, hinfo['name']),
                                    ROOT.Hist2D(*hinfo['binning'], name='{}__{}__{}'.format(ds, chan, hinfo['name']),
                                        title=hinfo['title'], drawstyle='colz', legendstyle='F')
                                    )
                events_.Histos['{}/cutflow'.format(chan)].name='{}__{}__cutflow'.format(ds, chan)
                events_.Histos['{}/cutflow'.format(chan)].legendstyle='L'
            events_.process()
            events_.postProcess()
            return ds, events_.histos

        ### signal 4mu
        if '4mu' in args.channel:
            sigDS_4mu_inc, sigSCALE_4mu_inc = sdml.fetch('4mu')
            if args.update_signalsample: sigDS_4mu_inc.update( json.load(open(args.update_signalsample)) )

            sigDS_4mu, sigSCALE_4mu = {}, {}
            for t in sampleSig:
                for k in sigDS_4mu_inc:
                    if not k.startswith(t): continue
                    sigDS_4mu[t] = sigDS_4mu_inc[k]
                for k in sigSCALE_4mu_inc:
                    if not k.startswith(t): continue
                    sigSCALE_4mu[t] = sigSCALE_4mu_inc[k]

            SigHists4mu = OrderedDict()

            for i, ds in enumerate(sampleSig, start=1):
                if ds not in sigDS_4mu or not sigDS_4mu[ds]: continue
                packages = []
                historesult = []
                pool = Pool(processes=min(len(sigDS_4mu[ds]), 2))
                for f in sigDS_4mu[ds]:
                    packages.append((ds, [f], sigSCALE_4mu[ds], args.maxevents, ['4mu',]))
                for res in tqdm(pool.imap_unordered(dofill, packages), total=len(packages)):
                    historesult.append(res[1])
                pool.close()
                pool.join()
                hists = historesult.pop()
                for k in hists:
                    for res in historesult:
                        hists[k].Add(res[k])
                SigHists4mu[ds] = hists
            log.info('channel 4mu filling done')


        ### signal 2mu2e
        if '2mu2e' in args.channel:
            sigDS_2mu2e_inc, sigSCALE_2mu2e_inc = sdml.fetch('2mu2e')
            if args.update_signalsample: sigDS_2mu2e_inc.update( json.load(open(args.update_signalsample)) )

            sigDS_2mu2e, sigSCALE_2mu2e = {}, {}
            for t in sampleSig:
                for k in sigDS_2mu2e_inc:
                    if not k.startswith(t): continue
                    sigDS_2mu2e[t] = sigDS_2mu2e_inc[k]
                for k in sigSCALE_2mu2e_inc:
                    if not k.startswith(t): continue
                    sigSCALE_2mu2e[t] = sigSCALE_2mu2e_inc[k]

            SigHists2mu2e = OrderedDict()

            for i, ds in enumerate(sampleSig, start=1):
                if ds not in sigDS_2mu2e or not sigDS_2mu2e[ds]: continue
                packages = []
                historesult = []
                pool = Pool(processes=min(len(sigDS_2mu2e[ds]), 2))
                for f in sigDS_2mu2e[ds]:
                    packages.append((ds, [f], sigSCALE_2mu2e[ds], args.maxevents, ['2mu2e',]))
                for res in tqdm(pool.imap_unordered(dofill, packages), total=len(packages)):
                    historesult.append(res[1])
                pool.close()
                pool.join()
                hists = historesult.pop()
                for k in hists:
                    for res in historesult:
                        hists[k].Add(res[k])
                SigHists2mu2e[ds] = hists
            log.info('channel 2mu2e filling done')


    if runbkg:
        if args.mbase=='proxy':
            from FireROOT.Analysis.DatasetMapLoader import ProxyEventsBkgDatasetMapLoader
            dml = ProxyEventsBkgDatasetMapLoader()
            bkgDS, bkgMAP, bkgSCALE = dml.fetch()
        else:
            bkgDS, bkgMAP, bkgSCALE = dml.fetch('bkg')

        BkgHists = {}
        for ds, files in tqdm(bkgDS.items()):
            events_ = imp.MyEvents(files=files, outname=outname, type='MC', maxevents=args.maxevents, dtag=ds, channel=args.channel)
            events_.setScale(bkgSCALE[ds])
            for chan in events_.channel:
                for hinfo in imp.histCollection:
                    if 'class' in hinfo: continue # do not support types except TH1, TH2
                    _dim = 0
                    for x in hinfo['binning']:
                        if isinstance(x, numbers.Number): _dim+=1/3.
                        elif isinstance(x, list): _dim+=1
                    if _dim==1:
                        events_.bookHisto('{}/{}'.format(chan, hinfo['name']),
                                ROOT.Hist(*hinfo['binning'], title=hinfo['title'],
                                    drawstyle='hist', fillstyle='solid', linewidth=0, legendstyle='F')
                                )
                    else:
                        events_.bookHisto('{}/{}'.format(chan, hinfo['name']),
                                ROOT.Hist2D(*hinfo['binning'], title=hinfo['title'],
                                    drawstyle='colz', legendstyle='F')
                                )
                events_.Histos['{}/cutflow'.format(chan)].decorate(fillstyle='solid', linewidth=0, legendstyle='F')
            events_.process()
            events_.postProcess()
            BkgHists[ds] = events_.histos
        log.info('background MC done')

    if rundata:
        if args.mbase=='proxy':
            from FireROOT.Analysis.DatasetMapLoader import ProxyEventsDataDatasetMapLoader
            ddml = ProxyEventsDataDatasetMapLoader()
            dataDS, dataMAP = ddml.fetch()
        else:
            dataDS, dataMAP = dml.fetch('data')

        DataHists = {}
        _files = []
        for ds in dataDS:
            if isinstance(dataDS[ds], str): _files.append(dataDS[ds])
            else: _files.extend(dataDS[ds])

        def dofill(pack):
            files, maxevents, channel = pack
            events_ = imp.MyEvents(files=files, outname=outname, type='DATA', maxevents=maxevents, channel=channel)
            for chan in events_.channel:
                for hinfo in imp.histCollection:
                    if 'class' in hinfo:
                        events_.bookHisto('{}/{}'.format(chan, hinfo['name']),
                                hinfo['class'](*hinfo['binning'], name='data__{}__{}'.format(chan, hinfo['name']),
                                    title=hinfo['title'], drawstyle='hist e1', legendstyle='LEP')
                                )
                    else:
                        _dim = 0
                        for x in hinfo['binning']:
                            if isinstance(x, numbers.Number): _dim+=1/3.
                            elif isinstance(x, list): _dim+=1
                        if _dim==1:
                            events_.bookHisto('{}/{}'.format(chan, hinfo['name']),
                                    ROOT.Hist(*hinfo['binning'], name='data__{}__{}'.format(chan, hinfo['name']),
                                        title=hinfo['title'], drawstyle='hist e1', legendstyle='LEP')
                                    )
                        else:
                            events_.bookHisto('{}/{}'.format(chan, hinfo['name']),
                                    ROOT.Hist2D(*hinfo['binning'], name='data__{}__{}'.format(chan, hinfo['name']),
                                        title=hinfo['title'], drawstyle='colz', legendstyle='F')
                                    )
                events_.Histos['{}/cutflow'.format(chan)].name = 'data__{}__cutflow'.format(chan)
                events_.Histos['{}/cutflow'.format(chan)].decorate(drawstyle='hist e1', legendstyle='LEP')
            events_.process()
            events_.postProcess()
            return events_.histos

        packages = []
        historesult = []
        pool = Pool(processes=min(len(_files), 2))

        for f in _files:
            packages.append(([f], args.maxevents, args.channel))
        for res in tqdm(pool.imap_unordered(dofill, packages), total=len(packages)):
            historesult.append(res)
        pool.close()
        pool.join()
        DataHists = historesult.pop()
        for k in DataHists:
            for res in historesult:
                DataHists[k].Add(res[k])
        log.info('data done')

    if not imp.histCollection:
        sys.exit(0)

    log.info('saving to {}'.format(outname))

    f = root_open(outname, args.create)
    if runsig:
        if '4mu' in args.channel:
            try: f.mkdir('ch4mu/sig', recurse=True)
            except: pass
            f.ch4mu.sig.cd()
            for ds, hs in SigHists4mu.items():
                try: f.ch4mu.sig.mkdir(ds)
                except: pass
                getattr(f.ch4mu.sig, ds).cd()
                for h in hs.values():
                    h.SetName( h.GetName().replace('{}__4mu__'.format(ds), '') )
                    if args.create=='update': h.Write('', ROOT.TObject.kOverwrite)
                    else: h.Write()
        if '2mu2e' in args.channel:
            try: f.mkdir('ch2mu2e/sig', recurse=True)
            except: pass
            f.ch2mu2e.sig.cd()
            for ds, hs in SigHists2mu2e.items():
                try: f.ch2mu2e.sig.mkdir(ds)
                except: pass
                getattr(f.ch2mu2e.sig, ds).cd()
                for h in hs.values():
                    h.SetName( h.GetName().replace('{}__2mu2e__'.format(ds), '') )
                    if args.create=='update': h.Write('', ROOT.TObject.kOverwrite)
                    else: h.Write()

    if rundata:
        try:
            for ch in args.channel:
                f.mkdir('ch{}/data'.format(ch), recurse=True)
        except:
            pass
        for h in DataHists.values():
            _, chan, name = h.GetName().split('__')
            getattr(f, 'ch{}'.format(chan)).data.cd()
            h.SetName(name)
            if args.create=='update': h.Write('', ROOT.TObject.kOverwrite)
            else: h.Write()

    if runbkg:
        try:
            for ch in args.channel:
                f.mkdir('ch{}/bkg'.format(ch), recurse=True)
        except:
            pass
        for chan in args.channel:
            getattr(f, 'ch{}'.format(chan)).bkg.cd()

            nameTitles = [] # get unique (name,title)
            for ds, hs in BkgHists.items():
                for k, h in hs.items():
                    if not k.startswith(chan): continue
                    nameTitles.append(
                        ( k.split('/')[-1],
                          ';'.join([h.title, h.axis(0).GetTitle(), h.axis(1).GetTitle()]) )
                    )
            nameTitles = list(set(nameTitles))

            for n, t in nameTitles:
                histName = '{}/{}'.format(chan, n)
                CatHists = mergeHistsFromMapping(extractHistByName(BkgHists, histName), bkgMAP, bkgCOLORS)
                hstack = ROOT.HistStack(list(CatHists.values()), name=n, title=t, drawstyle='HIST')
                if args.create=='update': hstack.Write('', ROOT.TObject.kOverwrite)
                else: hstack.Write()
    f.Close()
