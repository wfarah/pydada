import numpy as np
import os
from collections import OrderedDict

DEFAULT_HDR_SIZE = 4096
DEFAULT_NFFT = 4096
DEFAULT_ORDER = "TF"

class DadaHeader(OrderedDict):
    def __init__(self, fname):
        super(DadaHeader, self).__init__()
        self.fname = fname
        hdr_str = self.fname.read(DEFAULT_HDR_SIZE)
        self.parse_header(hdr_str)
        if 'ORDER' not in self:
            self['ORDER'] = DEFAULT_ORDER
    
    def parse_header(self, hdr_str):
        for token in hdr_str.split("\n"):
            if token.startswith("HDR_SIZE"):
                hdr_size = int(token.strip().split(" ")[-1])
                if hdr_size == DEFAULT_HDR_SIZE:
                    pass
                if hdr_size < DEFAULT_HDR_SIZE:
                    hdr_str = hdr_str[:hdr_size]
                if hdr_size > DEFAULT_HDR_SIZE:
                    self.fname.seek(0)
                    hdr_str = self.fname.read(hdr_size)
                break
        
        for token in hdr_str.split("\n"):
            token_cleaned = token.replace("\x00","")
            if not token_cleaned:
                pass
            token_split = token_cleaned.strip().replace("\t"," ").split(" ")
            key,val = token_split[0], token_split[-1]
            if not key and not val:
                continue
            
            # Check if value is boolean
            if val.lower() == "false":
                self[key] = False
                continue
            elif val.lower() == "true":
                self[key] = True
                continue
            # Check if value is int
            try:
                self[key] = int(val)
                continue
            except ValueError as e:
                pass
            
            # Check if value is float
            try:
                self[key] = float(val)
                continue
            except ValueError as e:
                pass
            
            # Value is a string
            self[key] = val
            


class DadaFile(object):
    def __init__(self, fname):
        self.filename = fname
        self.file = open(self.filename,"r")
        self.header = self.read_header()
        if 'NANT' not in self.header:
            self.header['NANT'] = 1
        self.nsamples = (os.path.getsize(self.filename) - self.header['HDR_SIZE'])/\
            self.header['NCHAN']/self.header['NDIM']/self.header['NBIT']/self.header['NANT']*8
        self.get_freq_order()
        
    def read_header(self):
        return DadaHeader(self.file)
    
    def get_freq_order(self):
        BW = float(self.header['BW'])
        FREQ = self.header['FREQ']
        NCHAN = self.header['NCHAN']
        if BW > 1:
            tmp = np.linspace(FREQ - (BW - BW/NCHAN)/2, 
                              FREQ + (BW - BW/NCHAN)/2,
                              NCHAN)
        else:
            tmp = np.linspace(FREQ + (BW - BW/NCHAN)/2, 
                              FREQ - (BW - BW/NCHAN)/2,
                              NCHAN)            
        self.freqs = tmp        
        
    
class DadaMux(object):
    def __init__(self, filenames):
        self.dadafiles = [DadaFile(fname) for fname in filenames]
        self.NBIT = self.uniq_hdr("NBIT")
        self.NDIM = self.uniq_hdr("NDIM")
        self.NCHAN = sum([dadafile.header['NCHAN'] for dadafile in self.dadafiles])
        if len(set([dada_file.nsamples for dada_file in self.dadafiles])) != 1:
            raise RuntimeError("Files do not have the same number of samples")
        self.NSAMPLES = self.dadafiles[0].nsamples
        self.ORDER = self.uniq_hdr("ORDER")
        self.NPOL = self.uniq_hdr("NPOL")
        self.TSAMP = self.uniq_hdr("TSAMP") * 1e-6
        
        if self.dadafiles[0].header['BW'] > 0:
            _reverse = False
        else:
            _reverse = True
            
        self.dadafiles.sort(key=lambda dada_file:dada_file.header['FREQ'], reverse=_reverse)
        
        self.freqs = np.array([freq for dadafile in self.dadafiles for freq in list(dadafile.freqs)])
        self.BW = sum([dadafile.header['BW'] for dadafile in self.dadafiles])
        self.fch1 = self.freqs[0]
        
    
    def allhdr_key(self, key):
        return [dadafile.header[key] for dadafile in self.dadafiles]
    
    def uniq_hdr(self, key):
        s = self.allhdr_key(key)
        if len(set(s)) != 1:
            raise ValueError('Expected the same header value for {} for all files. Got these values {}'.\
                             format(key, s))
        return set(s).pop()
        
    def read(self, start=0, nsamples=None):
        if self.NDIM == 2 and self.NBIT == 32:
            dtype = np.complex64
        else:
            raise NotImplementedError("Data type not yet implemented, please it to source code")
        
        if not nsamples:
            nsamples = self.NSAMPLES
        total_samples = nsamples - start
        assert total_samples <= self.NSAMPLES,\
            'Invalid read request. nsamp={} samp_start={} difference={}'.format(nsamples, start, nsamples-start)
        
        # Giant buffer for data
        
        #for dadafile 
        shape = tuple()
        if self.ORDER == "T":
            shape += (self.NCHAN,total_samples,)
        else:
            for key in self.ORDER:
                if key == "F":
                    shape += (self.NCHAN,)
                elif key == "T":
                    shape += (total_samples,)
        shape += (self.NPOL,)
        print shape
        
        #buffer = np.empty(shape=shape,dtype=dtype)
        buffer = BasebandData(shape=shape, dtype=dtype, NCHAN=self.NCHAN, NBIT=self.NBIT, 
                              freqs=self.freqs, fch1=self.fch1, TSAMP=self.TSAMP,BW=self.BW)
        chan_track = 0
        for dadafile in self.dadafiles:
            nchan = dadafile.header['NCHAN']
            if self.ORDER == "T":
                buffer[chan_track:chan_track+nchan] = np.fromfile(dadafile.file,
                        dtype=dtype, count=total_samples*self.NPOL).reshape((nchan,total_samples,self.NPOL))

            elif self.ORDER.index("F") == 0:
                buffer[chan_track:chan_track+nchan] = np.fromfile(dadafile.file, 
                    dtype=dtype, count=nchan*total_samples*self.NPOL).reshape((nchan,total_samples,self.NPOL))
                
            elif self.ORDER.index("F") == 1:
                buffer[:,chan_track:chan_track+nchan] = np.fromfile(dadafile.file, 
                    dtype=dtype, count=nchan*total_samples*self.NPOL).reshape((total_samples,nchan,self.NPOL))
            chan_track += nchan
                
        return np.squeeze(buffer).T
    


class BasebandData(np.ndarray):
    def __new__(cls, *args, **kwargs):
        ndarray_kw = ['shape', 'dtype',  'buffer', 'offset', 'strides', 'order']
        to_ndarray = {}
        to_myclass = {}
        for k,v in kwargs.items():
            if k not in ndarray_kw:
                to_myclass[k] = v
            else:
                to_ndarray[k] = v
        new = np.ndarray.__new__(cls, *args, **to_ndarray)
        for k,v in to_myclass.items():
            setattr(new, k, v)
        return new
    
    def __array_finalize__(self, obj):
        if obj is None: return
        attrs = ['NCHAN','NBIT','freqs','TSAMP','fch1','dedispersed','DM','BW']
        for attr in attrs:
            if hasattr(obj, attr): setattr(self, attr, getattr(obj, attr))
        setattr(obj, 'DM', 0)
    
    def detect(self):
        return np.abs(self)
    
    def get_DM_delays(self, DM):
        delays = DM * 4.148808e3 * ((self.fch1**-2) - (self.freqs**-2))
        return (delays/self.TSAMP).round().astype("int32")
        
    def dedisperse(self,DM,coherent=False, kernel_size=DEFAULT_NFFT):
        setattr(self, "dedispersed", True)
        setattr(self, 'DM', DM)
        if not coherent:
            dmdelays = self.get_DM_delays(DM)
            for ichan in range(self.shape[0]):
                self[ichan] = np.roll(self[ichan], dmdelays[ichan])
                
        if coherent:
            self.dedisperse(DM, coherent=False)
            tdm = 8.3*10**6*100*(self.BW/self.NCHAN)*self.freqs**(-3) * 10**3
            nfft_skip = np.round(max(tdm*self.BW/self.NCHAN))
            nfft_skip = next_power_of_2(int(nfft_skip))
            foff = self.BW/self.NCHAN
            chan_padded = np.zeros((self[0].size + nfft_skip), dtype=self.dtype)
            
            if self[0].size%kernel_size != 0:
                raise RuntimeError("Kernel size (%i) should divide the total "
                                   "number of samples (%i)" %(kernel_size, self[0].size))
            if kernel_size <= nfft_skip:
                raise RuntimeError("Kernel size (%i) should be less than the "
                                  "Wrap-around region (%i)" %(kernel_size, nfft_skip))
                
            kernel_size_updated = kernel_size + nfft_skip
    
            for ichan,freq0 in enumerate(self.freqs):
                chan_padded[-nfft_skip/2:nfft_skip/2] = 0
                chan_padded[nfft_skip/2:-nfft_skip/2] = self[ichan].copy()
                kk = [chan_padded[i:i+kernel_size_updated] 
                      for i in range(0, len(chan_padded), kernel_size_updated-nfft_skip)]
                segments = np.array(kk[:-1])
                
                phasors = get_dispersion_phasors(DM, freq0, self.BW/self.NCHAN, kernel_size_updated)
                transform = np.fft.fftshift(np.fft.fft(segments, n=kernel_size_updated, axis=1))
                
                transform *= np.exp(-1j*phasors)
                segments[:] = np.fft.ifft(np.fft.ifftshift(transform), axis=1)
                
                self[ichan] = segments[:,nfft_skip/2:-nfft_skip/2].flatten()
                
                
                

    
    def time_scrunch(self):
        if self.dtype not in [np.float128, np.float64, np.float32, np.float16, np.float,
                             np.int64, np.int32, np.int16, np.int8]:
            raise RuntimeError("dtype (%s) is not suitable for time scrunching" %self.dtype)
        ts = TimeSeries(shape=(self.shape[1],), dtype=self.dtype, TSAMP=self.TSAMP)
        ts[:] = self.sum(axis=0)
        return ts
    
    
    def spec_fft(self, nfft):
        finefilbank = self.copy()
        finefilbank.NCHAN = self.NCHAN * nfft
        finefilbank.TSAMP = self.TSAMP * nfft
        finefilbank.freqs = np.array([hires_chan for chan in self.freqs 
                              for hires_chan in 
                              np.linspace(-self.BW/self.NCHAN/2. + chan, self.BW/self.NCHAN/2. + chan, nfft)])
        finefilbank.fch1 = finefilbank.freqs[0]
        finefilbank.shape = (finefilbank.NCHAN, self.shape[1]/nfft)
        
        dfft = np.zeros_like(self[0]).reshape(-1, nfft)

        for i,chan in enumerate(self):
            dfft[:] = chan.copy().reshape(-1, nfft)
            finefilbank[i*nfft:(i+1)*nfft] = np.fft.fftshift(np.fft.fft(dfft, axis=1), axes=1).T
        
        return finefilbank
    
    def spec_ifft(self, nfft):
        coarsefilbank = self.copy()
        
        coarsefilbank.NCHAN = self.NCHAN / nfft
        coarsefilbank.TSAMP = self.TSAMP / nfft
        coarsefilbank.freqs = np.mean(self.freqs.reshape(-1,nfft), axis=1)
        coarsefilbank.fch1 = coarsefilbank.freqs[0]
        
        coarsefilbank.shape = (coarsefilbank.NCHAN, coarsefilbank.shape[1]*nfft)
        
        difft = np.zeros((len(self[0]),nfft), dtype=coarsefilbank.dtype)
        print difft.shape
        
        for i,coarsechan in enumerate(coarsefilbank):
            difft[:] = self[i*nfft:(i+1)*nfft].T.copy()
            coarsechan[:] = np.fft.ifft(np.fft.ifftshift(difft, axes=1), axis=1).flatten()
        
        return coarsefilbank
    
    
    

class TimeSeries(np.ndarray):
    def __new__(cls, *args, **kwargs):
        ndarray_kw = ['shape', 'dtype',  'buffer' 'offset', 'strides', 'order']
        to_ndarray = {}
        to_myclass = {}
        for k,v in kwargs.items():
            if k not in ndarray_kw:
                to_myclass[k] = v
            else:
                to_ndarray[k] = v
        new = np.ndarray.__new__(cls, *args, **to_ndarray)
        for k,v in to_myclass.items():
            setattr(new, k, v)
        return new
    
    def __array_finalize__(self, obj):
        if obj is None: return
        attrs = ['TSAMP']
        for attr in attrs:
            if hasattr(obj, attr): setattr(self, attr, getattr(obj, attr))
    
    
    def normalise(self):
        self/=(1.4826*self._get_mad())
    
    def _get_mad(self):
        return np.median(abs(self - np.median(self)))
    
    def remove_baseline(self):
        self -= np.median(self)
        
def get_dispersion_phasors(DM, freq0, bw, nfft): 
    subfreqs = np.linspace(-bw/2.,+bw/2.,nfft)
    phasors = 2*np.pi*DM/(2.410331*10**(-10))*subfreqs*subfreqs/(freq0*freq0*(subfreqs+freq0))
    return phasors

def next_power_of_2(x):  
    return 1 if x == 0 else 2**(x - 1).bit_length()
