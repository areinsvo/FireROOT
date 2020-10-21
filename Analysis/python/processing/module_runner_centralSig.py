#!/usr/bin/env python
from __future__ import print_function
import argparse
import math, numbers
import os, sys, fnmatch
from multiprocessing import Pool

import tqdm
from FireROOT.Analysis.Utils import *
from FireROOT.Analysis.DatasetMapLoader import CentralSignalMapLoader

from rootpy.io import root_open
from rootpy.logger import log
log = log[__name__]

## parser
parser = argparse.ArgumentParser(description="module runner for central signal samples.")
parser.add_argument("--module", "-m", type=str, help='module name')
parser.add_argument("--outname", "-o", type=str, default=None, help='output ROOT file name')
parser.add_argument("--maxevents", "-n", type=int, default=-1, help='max number of events to run')
parser.add_argument("--create", "-c", type=str, default='recreate', choices=['recreate', 'update'], help='update output by')
parser.add_argument("--sigparam", "-p",  type=str, nargs='*',default=None, help='signal parameters')
parser.add_argument("--mbase", "-b", type=str, default='modules', choices=['modules', 'centralSig'], help='module base name')
parser.add_argument("--channel", "-ch", nargs='*', default=['2mu2e', '4mu'], choices=['2mu2e', '4mu'], help='channels to run')
args = parser.parse_args()

def args_sanity(args):
    # module
    moduleBase = os.path.join(os.getenv('CMSSW_BASE'), 'src/FireROOT/Analysis/python/processing/%s'%args.mbase)
    if not os.path.isdir(moduleBase): sys.exit('Module base not exist!')
    allmodules = [
        fn.split('.')[0] for fn in os.listdir(moduleBase) \
        if os.path.isfile(os.path.join(moduleBase, fn)) \
        and not fn.startswith('_') \
        and fn.endswith('.py')
    ]
    if args.module not in allmodules:
        sys.exit('Available modules: {}'.format(str(allmodules)))

args_sanity(args)
_modulebase = 'FireROOT.Analysis.processing.{}'.format(args.mbase)
try:
    imp = __import__('{}.{}'.format(_modulebase, args.module), fromlist=['MyEvents', 'histCollection'])
except Exception as e:
    print(e)
    print('Fail to import; run `scram b`.')



if __name__ == '__main__':

    outdir = os.path.join(os.getenv('CMSSW_BASE'), 'src/FireROOT/Analysis/python/outputs/rootfiles/centralSig')
    if not os.path.isdir(outdir): os.makedirs(outdir)
    if args.outname:
        outname = os.path.join(outdir, '{}.root'.format(args.outname))
    else:
        outname = os.path.join(outdir, '{}.root'.format(args.module))
    if args.create == 'update' and not os.path.isfile(outname):
        sys.exit('UPDATE was used not file not recreated yet.')

    sdml = CentralSignalMapLoader()
    sigDS_2mu2e, sigSCALE_2mu2e = sdml.fetch('2mu2e')
    sigDS_4mu, sigSCALE_4mu = sdml.fetch('4mu')

    sampleSig = list(set(sigDS_2mu2e.keys()) | set(sigDS_4mu.keys()))
    if args.sigparam:
        _selected = []
        for s in args.sigparam:
            if '*' in s or '?' in s:
                _selected.extend( fnmatch.filter(sampleSig, s) )
            else: _selected.append(s)
        sampleSig = list(set(_selected))
        print(sampleSig)

    def dofill(pack):
        ds, files, scale, maxevents, channel = pack
        events_ = imp.MyEvents(files=files, type='MC', dtag=ds, maxevents=maxevents, channel=channel)
        events_.setScale(scale)

        for chan in events_.channel:
            for hinfo in imp.histCollection:
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

    # tqdm: https://github.com/tqdm/tqdm/issues/484#issuecomment-351001534

    ### signal 4mu
    if '4mu' in args.channel:
        packages = []
        SigHists4mu = []
        pool = Pool(processes=2)

        for i, ds in enumerate(sampleSig, start=1):
            if ds not in sigDS_4mu or not sigDS_4mu[ds]: continue
            packages.append((ds, sigDS_4mu[ds], sigSCALE_4mu[ds], args.maxevents, ['4mu',]))
        for res in tqdm.tqdm(pool.imap_unordered(dofill, packages), total=len(packages)):
            SigHists4mu.append(res)
        pool.close()
        pool.join()
        SigHists4mu = dict(SigHists4mu)
        log.info('channel 4mu filling done')


    ### signal 2mu2e
    if '2mu2e' in args.channel:
        packages = []
        SigHists2mu2e = []
        pool = Pool(processes=2)

        for i, ds in enumerate(sampleSig, start=1):
            if ds not in sigDS_2mu2e or not sigDS_2mu2e[ds]: continue
            packages.append((ds, sigDS_2mu2e[ds], sigSCALE_2mu2e[ds], args.maxevents, ['2mu2e',]))

        for res in tqdm.tqdm(pool.imap_unordered(dofill, packages), total=len(packages)):
            SigHists2mu2e.append(res)
        pool.close()
        pool.join()
        SigHists2mu2e = dict(SigHists2mu2e)
        log.info('channel 2mu2e filling done')


    log.info('saving to {}'.format(outname))
    f = root_open(outname, args.create)
    if '4mu' in args.channel:
        try: f.mkdir('ch4mu')
        except: pass
        chanDir = f.ch4mu
        chanDir.cd()
        for ds, hs in SigHists4mu.items():
            mxxma = ds.rsplit('_', 2)[0]
            lxyctau = ds.split('_', 2)[-1]
            try: chanDir.mkdir(mxxma+'/'+lxyctau, recurse=True)
            except: pass
            getattr(getattr(chanDir, mxxma), lxyctau).cd()
            for h in hs.values():
                h.SetName( h.GetName().replace('{}__4mu__'.format(ds), '') )
                if args.create=='update': h.Write('', ROOT.TObject.kOverwrite)
                else: h.Write()

    if '2mu2e' in args.channel:
        try: f.mkdir('ch2mu2e')
        except: pass
        chanDir = f.ch2mu2e
        chanDir.cd()
        for ds, hs in SigHists2mu2e.items():
            mxxma = ds.rsplit('_', 2)[0]
            lxyctau = ds.split('_', 2)[-1]
            try: chanDir.mkdir(mxxma+'/'+lxyctau, recurse=True)
            except: pass
            getattr(getattr(chanDir, mxxma), lxyctau).cd()
            for h in hs.values():
                h.SetName( h.GetName().replace('{}__2mu2e__'.format(ds), '') )
                if args.create=='update': h.Write('', ROOT.TObject.kOverwrite)
                else: h.Write()
    f.close()
