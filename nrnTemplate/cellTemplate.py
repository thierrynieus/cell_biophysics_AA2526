from neuron import h

class HHcell:
    def __init__(self, synPar={}, rndPar={}, RecAll=2, coord={},varDt=False):    
        '''
            Defines a HH cell model
            synPar:
            rndPar:    specifies parameter for random input to the cell ... (needs at least 1 exc syn !!)
            RecAll:      0 save nothing, 1 saves voltage, 2 saves also synaptic current
            coord:
            
        '''
        self.RecAll=RecAll
        self.soma = h.Section(cell=self)
        self.soma.nseg  = 1 
        self.soma.diam  = 10
        self.soma.L  = 10
        self.soma.cm = 1        
        self.soma.Ra = 100	
        #self.soma.celsius=30

        self.soma.insert('hh')
        
        self.soma.ena = 87.39
        self.soma.ek = -84.69
        #self.soma.eca = 129.33
        
        h.celsius = 30
        
        self.I=h.IClamp(self.soma(0.5))

        self.pos = {}
        self.pos['x'] = coord['x']
        self.pos['y'] = coord['y']        

        self.ExcSynList = []
        self.InhSynList = []
        self.connList = []

        # determine neural type
        if  synPar['type'].rfind('EXC')>=0:
            self.whatami='EXC'
        elif synPar['type'].rfind('INH')>=0:
            self.whatami='INH'
        else:
            self.whatami='UNDEF'        
            print("Neither excitatory nor inhibitory neuron. The neural type if undefined and no synapse will be formed!")

        self.weight_e = synPar['weight_e']
        self.weight_i = synPar['weight_i']

        # record spike time stamps
        self.record = {}
        self.record['Iexc'] = []    
        self.record['Iinh'] = []
        self.nc_spike = h.NetCon(self.soma(0.5)._ref_v, None, -20, 0, 1, sec=self.soma)
        #self.nc_spike.delay=0
        self.record['spk'] = h.Vector()
        self.nc_spike.record(self.record['spk'])
        
        # inject synaptic noise     
        self.RandomStim = h.SpikeGenerator(0.5, sec=self.soma) 
        self.synNoise = h.Exp2Syn(self.soma(0.5))
        self.RandomInput = h.NetCon(self.RandomStim,self.synNoise, -20, 0, 1, sec=self.soma) # create connection
        #print(rndPar['weight'])
        self.RandomInput.weight[0] = rndPar['weight']

        # record synaptic noise
        self.nc_spike_rnd = h.NetCon(self.RandomStim,None,-20,0,1,sec=self.soma)
        self.record['spk_rnd'] = h.Vector()
        self.nc_spike_rnd.record(self.record['spk_rnd'])
              
        if varDt:
            # setup the "LOCAL TIME STEP" integrator
            self.Hines = h.CVode()
            self.Hines.active(1)
            self.Hines.use_local_dt(1)

        if RecAll>0:
            self.record['t'] = h.Vector()
            self.record['vm'] = h.Vector()              
            if varDt:        
                self.Hines.record(self.soma(0.5)._ref_v, self.record['vm'], self.record['t'], sec=self.soma)  
            else:           
                self.record['t'].record(h._ref_t)                 
                self.record['vm'].record(self.soma(0.5)._ref_v) 
        
    def createSyn(self,InputWhatAmi):
        '''
        '''
        if InputWhatAmi.rfind('EXC')>=0:
            self.ExcSynList.append(h.Exp2Syn(self.soma(0.5)))
            self.ExcSynList[-1].e=0
            
        elif InputWhatAmi.rfind('INH')>=0:     
            self.InhSynList.append(h.Exp2Syn(self.soma(0.5)))  
            self.InhSynList[-1].e=-65
               
        else:
            print("Since ",InputWhatAmi," is neither an excitatory nor an inhibitory neuron not connection is made!")
        
    def connect_to(self, dest):
        '''
            connect last appended synapse from self to dest (this completes previous call createSyn from dest cell)
        '''
        if self.whatami.rfind('EXC')>=0:   
            #connTMP = h.NetCon(self.soma(0.5)._ref_v,dest.ExcSynList[-1],sec=self.soma)        
            connTMP=h.NetCon(self.soma(0.5)._ref_v, dest.ExcSynList[-1], sec=self.soma)
            connTMP.weight[0] = self.weight_e
            connTMP.delay = 0
            self.connList.append(connTMP)  
            #
            if self.RecAll==2:
                dest.record['Iexc'].append(h.Vector())   
                dest.record['Iexc'][-1].record(dest.ExcSynList[-1]._ref_i, sec=self.soma)                           

        elif self.whatami.rfind('INH')>=0:       
            connTMP=h.NetCon(self.soma(0.5)._ref_v,dest.InhSynList[-1], sec=self.soma)
            connTMP.weight[0] = self.weight_i
            connTMP.delay = 0
            self.connList.append(connTMP)  
            if self.RecAll==2:
                dest.record['Iinh'].append(h.Vector()) 
                dest.record['Iinh'][-1].record(dest.InhSynList[-1]._ref_i,sec=self.soma)
        else:
            print("Since ",self.whatami," is neither an excitatory nor an inhibitory neuron not connection is made!")


    def SetSynNoise(self,rndPar):
        '''
        '''
        # RANDOM STIM
        if 'start' in rndPar: self.RandomStim.start=rndPar['start']  
        if 'end' in rndPar: self.RandomStim.end=rndPar['end']
        if 'fast_invl' in rndPar: self.RandomStim.fast_invl=rndPar['fast_invl']
        if 'slow_invl' in rndPar: self.RandomStim.slow_invl=rndPar['slow_invl']
        if 'burst_len' in rndPar: self.RandomStim.burst_len=rndPar['burst_len']
        if 'delay' in rndPar: self.RandomStim.delay=rndPar['delay']
        if 'noise' in rndPar: self.RandomStim.noise=rndPar['noise']
        if 'seed' in rndPar: self.RandomStim.seed(rndPar['seed'])


    def destroy(self):
        del self.connList
        del self.ExcSynList
        del self.InhSynList
        del self.RandomStim 
        del self.RandomInput
        del self.synNoise
        del self.nc_spike
        del self.record
        del nc_spike
        if RecAll>0: 
            del self.record
        
        
