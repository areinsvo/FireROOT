#!/usr/bin/env python

from FireROOT.Analysis.Events import *
from FireROOT.Analysis.Utils import *


#ARH: need to check ntuples to see where daudr and some of the other daughter variables come from

class MyEvents(SignalEvents):
    def __init__(self, files=None, type='MC', maxevents=-1, channel=['2mu2e', '4mu'], **kwargs):
        super(MyEvents, self).__init__(files=files, type=type, maxevents=maxevents, channel=channel, **kwargs)

    def processEvent(self, event, aux):
        if aux['channel'] not in self.Channel: return
        chan = aux['channel']

        dp_toMu = [p for p in aux['dp'] if p.daupid==13]
        dp_toEl = [p for p in aux['dp'] if p.daupid==11]

        DR_THRESHOLD = 0.4 # maximum dR limit between (gen darkphoton, lepton-jets) matching

        # dp->mu, mu
#ARH: maybe check how many dark photons there are, so we know we aren't double counting any? Probably easier than trying to spot a mistake in how Weinan has filled the dp collection
        for dp in dp_toMu:
            if abs(dp.p4.eta())>2.4: continue

#ARH: interesting, so lxy here is the actual gen distance that the dp travels before decaying
            lxy = (dp.dauvtx - dp.vtx).Rho()
            lz  = (dp.dauvtx - dp.vtx).Z()
            # if abs(lz)>800: continue
            if dp.p4.pt()>30:
                self.Histos['%s/lxyDpToMu__total' % chan].Fill(lxy)
                self.Histos['%s/lepDrDpToMu__total' % chan].Fill(dp.daudr)
                self.Histos['%s/etaDpToMu__total' % chan].Fill(dp.p4.eta())

            #Don't apply pt requirement when plotting eff wrt pt
            self.Histos['%s/ptDpToMu__total' % chan].Fill(dp.p4.pt())             
            self.Histos['%s/ptDp' % chan].Fill(dp.p4.pt())

            mindr, matched = 999., None
            for lj in event.leptonjets:
                if not lj.isMuonType(): continue
#ARH: passSelection includes pt, eta, minTwoTrackDist, max(dz)> 40 for muons, and muon charge cut
                if not lj.passSelection(event): continue
                distance = DeltaR(dp.p4, lj.p4)
                if distance > DR_THRESHOLD: continue
                if distance<mindr:
                    mindr = distance
                    matched = lj

#ARH: efficiency vs pt might not really work because passSelection already includes a pt cut on the LJ
            if matched:
                if dp.p4.pt()>30:
                    self.Histos['%s/lxyDpToMu__match' % chan].Fill(lxy)
                    self.Histos['%s/lepDrDpToMu__match' % chan].Fill(dp.daudr)
                    self.Histos['%s/etaDpToMu__match' % chan].Fill(dp.p4.eta())

                self.Histos['%s/ptDpToMu__match' % chan].Fill(dp.p4.pt())


        # dp->el, el
        for dp in dp_toEl:
            if abs(dp.p4.eta())>2.4: continue

            #Don't apply pt cut when looking at eff with respect to pt
            self.Histos['%s/ptDpToEl__total' % chan].Fill(dp.p4.pt())
            self.Histos['%s/ptDp' % chan].Fill(dp.p4.pt())

            lxy = (dp.dauvtx - dp.vtx).Rho()
            lz  = (dp.dauvtx - dp.vtx).Z()

            if dp.p4.pt()>30: 
                self.Histos['%s/lxyDpToEl__total' % chan].Fill(lxy)
                self.Histos['%s/lepDrDpToEl__total' % chan].Fill(dp.daudr)
                self.Histos['%s/etaDpToEl__total' % chan].Fill(dp.p4.eta())

            mindr, matched = 999., None
            for lj in event.leptonjets:
                if not lj.isEgmType(): continue
                if not lj.passSelection(event): continue
                distance = DeltaR(dp.p4, lj.p4)
                if distance > DR_THRESHOLD: continue
                if distance<mindr:
                    mindr = distance
                    matched = lj

#ARH: efficiency vs pt might not really work because passSelection already includes a pt cut on the LJ 
            if matched is not None:
                if dp.p4.pt()>30:
                    self.Histos['%s/lxyDpToEl__match' % chan].Fill(lxy)
                    self.Histos['%s/lepDrDpToEl__match' % chan].Fill(dp.daudr)
                    self.Histos['%s/etaDpToEl__match' % chan].Fill(dp.p4.eta())
                self.Histos['%s/ptDpToEl__match' % chan].Fill(dp.p4.pt())


            else:
                # matching with PFElectrons
                _mindr, _matched = 999., None
                for e in event.electrons:
                    if e.p4.pt()<10 or abs(e.p4.eta())>2.4: continue
                    distance = DeltaR(dp.p4, e.p4)
                    if distance > DR_THRESHOLD: continue
                    else:
                        _mindr = distance
                        _matched = e
                if _matched:
                    self.Histos['%s/lxyDpToEl__matchEle' % chan].Fill(lxy)

                # matching with PFPhotons
                _mindr, _matched = 999., None
                for g in event.photons:
                    if g.p4.pt()<10 or abs(g.p4.eta())>2.4: continue
                    distance = DeltaR(dp.p4, g.p4)
                    if distance > DR_THRESHOLD: continue
                    else:
                        _mindr = distance
                        _matched = g
                if _matched:
                    self.Histos['%s/lxyDpToEl__matchPho' % chan].Fill(lxy)



histCollection = [
    {
        'name'   : 'lxyDpToMu__total',
        'binning': (100, 0, 500),
        'title'  : 'Z_{d} lxy;lxy [cm];counts/5cm',
    },
    {
        'name'   : 'lxyDpToMu__match',
        'binning': (100, 0, 500),
        'title'  : 'Z_{d} lxy;lxy [cm];counts/5cm',
    },
    {
        'name'   : 'lxyDpToEl__total',
        'binning': (100, 0, 250),
        'title'  : 'Z_{d} lxy;lxy [cm];counts/2.5cm',
    },
    {
        'name'   : 'lxyDpToEl__match',
        'binning': (100, 0, 250),
        'title'  : 'Z_{d} lxy;lxy [cm];counts/2.5cm',
    },
    {
        'name'   : 'ptDpToMu__total',
        'binning': (100, 0, 500),
        'title'  : 'Z_{d} pt;pt [GeV];counts/5 GeV',
    },
    {
        'name'   : 'ptDp',
        'binning': (140, 0, 700),
        'title'  : 'Z_{d} pt;pt [GeV];counts/5 GeV',
    },
    {
        'name'   : 'ptDpToMu__match',
        'binning': (100, 0, 500),
        'title'  : 'Z_{d} pt;pt [GeV];counts/5 GeV',
    },
    {
        'name'   : 'ptDpToEl__total',
        'binning': (100, 0, 500),
        'title'  : 'Z_{d} pt;pt [GeV];counts/5 GeV',
    },
    {
        'name'   : 'ptDpToEl__match',
        'binning': (100, 0, 500),
        'title'  : 'Z_{d} pt;pt [GeV];counts/5 GeV',
    },
    ## matching with PFElectrons
    {
        'name'   : 'lxyDpToEl__matchEle',
        'binning': (100, 0, 250),
        'title'  : 'Z_{d} lxy;lxy [cm];counts/2.5cm',
    },
    {
        'name'   : 'lxyDpToEl__matchPho',
        'binning': (100, 0, 250),
        'title'  : 'Z_{d} lxy;lxy [cm];counts/2.5cm',
    },


    ## opening angle
    {
        'name'   : 'lepDrDpToMu__total',
        'binning': (50, 0, 0.5),
        'title'  : 'Z_{d} #DeltaR(#mu^{+}#mu^{-});#DeltaR(#mu^{+}#mu^{-});counts/0.01',
    },
    {
        'name'   : 'lepDrDpToMu__match',
        'binning': (50, 0, 0.5),
        'title'  : 'Z_{d} #DeltaR(#mu^{+}#mu^{-});#DeltaR(#mu^{+}#mu^{-});counts/0.01',
    },
    {
        'name'   : 'lepDrDpToEl__total',
        'binning': (50, 0, 0.5),
        'title'  : 'Z_{d} #DeltaR(e^{+}e^{-});#DeltaR(e^{+}e^{-});counts/0.01',
    },
    {
        'name'   : 'lepDrDpToEl__match',
        'binning': (50, 0, 0.5),
        'title'  : 'Z_{d} #DeltaR(e^{+}e^{-});#DeltaR(e^{+}e^{-});counts/0.01',
    },
    ## eta
    {
        'name'   : 'etaDpToMu__total',
        'binning': (50, -2.5, 2.5),
        'title'  : 'Z_{d} eta;eta ;counts',
    },
    {
        'name'   : 'etaDpToMu__match',
        'binning': (50, -2.5, 2.5),
        'title'  : 'Z_{d} eta;eta ;counts',
    },
    {
        'name'   : 'etaDpToEl__total',
        'binning': (50, -2.5, 2.5),
        'title'  : 'Z_{d} eta;eta ;counts',
    },
    {
        'name'   : 'etaDpToEl__match',
        'binning': (50, -2.5, 2.5),
        'title'  : 'Z_{d} eta;eta ;counts',
    }
]
